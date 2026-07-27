from __future__ import annotations

import io
import signal
import unittest

from cladistica.progress import ERASE_CURRENT_LINE, PipelineProgress


class TTYBuffer(io.StringIO):
    def __init__(self, columns: int = 80) -> None:
        super().__init__()
        self.columns = columns

    def isatty(self) -> bool:
        return True

    def terminal_width(self) -> int:
        return self.columns

    def resize(self, columns: int) -> None:
        self.columns = columns


class ResizeDuringWriteBuffer(TTYBuffer):
    def __init__(self, columns: int = 80) -> None:
        super().__init__(columns)
        self.progress: PipelineProgress | None = None
        self.triggered = False

    def write(self, text: str) -> int:
        written = super().write(text)
        if (
            not self.triggered
            and self.progress is not None
            and text.startswith(ERASE_CURRENT_LINE)
            and "Bootstrap" in text
        ):
            self.triggered = True
            self.resize(40)
            self.progress._handle_sigwinch(getattr(signal, "SIGWINCH", 0), None)
        return written


def rendered_lines(output: str) -> list[str]:
    lines: list[str] = []
    for chunk in output.split(ERASE_CURRENT_LINE)[1:]:
        line = chunk.split("\r", 1)[0].split("\n", 1)[0]
        if line:
            lines.append(line)
    return lines


def progress_for(
    stream: TTYBuffer,
    stages: list[tuple[str, str]] | None = None,
) -> PipelineProgress:
    return PipelineProgress(
        stages or [("bootstrap", "Bootstrap replicates")],
        stream=stream,
        terminal_width_provider=stream.terminal_width,
    )


class ProgressTests(unittest.TestCase):
    def test_output_never_reaches_the_terminal_last_column(self) -> None:
        for columns in (120, 80, 60, 40, 25):
            with self.subTest(columns=columns):
                stream = TTYBuffer(columns)
                progress = progress_for(stream)
                progress.set_progress(
                    "bootstrap",
                    25,
                    "250/1000 bootstrap replicates",
                )
                rendered = rendered_lines(stream.getvalue())[-1]
                self.assertLessEqual(len(rendered), columns - 1)

    def test_resize_from_80_to_40_uses_the_new_width(self) -> None:
        stream = TTYBuffer(80)
        progress = progress_for(stream)
        progress.set_progress("bootstrap", 25, "250/1000 replicates")
        before = len(stream.getvalue())
        stream.resize(40)
        progress.set_progress("bootstrap", 26, "260/1000 replicates")
        delta = stream.getvalue()[before:]
        self.assertTrue(delta.startswith(ERASE_CURRENT_LINE + "\n" + ERASE_CURRENT_LINE))
        self.assertLessEqual(len(rendered_lines(delta)[-1]), 39)

    def test_resize_from_40_to_100_terminates_the_old_logical_line(self) -> None:
        stream = TTYBuffer(40)
        progress = progress_for(stream)
        progress.set_progress("bootstrap", 25, "old-detail")
        before = len(stream.getvalue())
        stream.resize(100)
        progress.set_progress("bootstrap", 26, "new-detail")
        delta = stream.getvalue()[before:]
        self.assertTrue(delta.startswith(ERASE_CURRENT_LINE + "\n" + ERASE_CURRENT_LINE))
        self.assertNotIn("old-detail", delta)
        self.assertIn("new-detail", rendered_lines(delta)[-1])

    def test_every_interactive_draw_starts_with_erase_line(self) -> None:
        stream = TTYBuffer(80)
        progress = progress_for(stream)
        progress.start("bootstrap")
        progress.set_progress("bootstrap", 25, "250/1000")
        progress.set_progress("bootstrap", 50, "500/1000")
        progress.succeed("bootstrap", "1000/1000")
        self.assertEqual(stream.getvalue().count(ERASE_CURRENT_LINE), 4)
        for chunk in stream.getvalue().split(ERASE_CURRENT_LINE)[1:]:
            self.assertFalse(chunk.startswith("\x1b"))

    def test_running_updates_do_not_write_newlines_until_success(self) -> None:
        stream = TTYBuffer(80)
        progress = progress_for(stream)
        progress.start("bootstrap")
        progress.set_progress("bootstrap", 25, "250/1000")
        progress.set_progress("bootstrap", 50, "500/1000")
        self.assertNotIn("\n", stream.getvalue())
        progress.succeed("bootstrap", "1000/1000")
        self.assertEqual(stream.getvalue().count("\n"), 1)
        self.assertTrue(stream.getvalue().endswith("\n"))

    def test_shorter_detail_cannot_leave_an_old_suffix(self) -> None:
        stream = TTYBuffer(100)
        progress = progress_for(stream)
        progress.start("bootstrap", "a very long obsolete detail suffix")
        before = len(stream.getvalue())
        progress.update("bootstrap", "short")
        delta = stream.getvalue()[before:]
        self.assertTrue(delta.startswith(ERASE_CURRENT_LINE))
        self.assertEqual(len(rendered_lines(delta)), 1)
        self.assertIn("short", rendered_lines(delta)[0])
        self.assertNotIn("obsolete", rendered_lines(delta)[0])

    def test_non_tty_never_emits_ansi_escape_sequences(self) -> None:
        stream = io.StringIO()
        progress = PipelineProgress(
            [("align", "MUSCLE marker alignments")],
            stream=stream,
        )
        progress.start("align", "rbcL: 4 sequences")
        progress.set_progress("align", 50, "1/2 markers")
        progress.succeed("align", "2 marker alignments")
        self.assertNotIn("\x1b", stream.getvalue())
        self.assertNotIn("\r", stream.getvalue())

    def test_non_tty_plain_output_is_unchanged(self) -> None:
        stream = io.StringIO()
        with PipelineProgress(
            [("align", "MUSCLE marker alignments"), ("ml", "ML tree")],
            stream=stream,
        ) as progress:
            progress.start("align", "rbcL: 4 sequences")
            progress.succeed("align", "1 marker alignment")
            progress.skip("ml", "Skipped by option")

        self.assertEqual(
            stream.getvalue(),
            "[>>] MUSCLE marker alignments: rbcL: 4 sequences\n"
            "[OK] MUSCLE marker alignments: 1 marker alignment\n"
            "[--] ML tree: Skipped by option\n",
        )

    def test_keyboard_interrupt_replaces_active_line_with_one_failed_line(self) -> None:
        stream = TTYBuffer(80)
        progress = progress_for(stream, [("align", "MUSCLE marker alignments")])
        progress.start("align", "rbcL")
        progress.__exit__(KeyboardInterrupt, KeyboardInterrupt(), None)
        output = stream.getvalue()
        self.assertEqual(progress.stages["align"].state, "failed")
        self.assertEqual(output.count("\n"), 1)
        self.assertTrue(output.endswith("\n"))
        failed_line = rendered_lines(output)[-1]
        self.assertIn("!0% Alignment", failed_line)
        self.assertIn("KeyboardInterrupt", failed_line)

    def test_sigwinch_redraws_the_active_stage_at_the_new_width(self) -> None:
        stream = TTYBuffer(80)
        progress = progress_for(stream)
        progress.start("bootstrap", "250/1000 replicates")
        before = len(stream.getvalue())
        stream.resize(35)
        progress._handle_sigwinch(getattr(signal, "SIGWINCH", 0), None)
        delta = stream.getvalue()[before:]
        self.assertTrue(delta.startswith(ERASE_CURRENT_LINE + "\n" + ERASE_CURRENT_LINE))
        self.assertLessEqual(len(rendered_lines(delta)[-1]), 34)

    def test_width_provider_is_consulted_for_every_draw(self) -> None:
        stream = TTYBuffer(80)
        calls: list[int] = []

        def provider() -> int:
            calls.append(stream.columns)
            return stream.columns

        progress = PipelineProgress(
            [("bootstrap", "Bootstrap replicates")],
            stream=stream,
            terminal_width_provider=provider,
        )
        progress.start("bootstrap")
        progress.set_progress("bootstrap", 25)
        progress.succeed("bootstrap")
        self.assertEqual(calls, [80, 80, 80])

    def test_sigwinch_during_a_stream_write_is_deferred_safely(self) -> None:
        stream = ResizeDuringWriteBuffer(80)
        progress = progress_for(stream)
        stream.progress = progress
        progress.start("bootstrap", "250/1000 replicates")
        lines = rendered_lines(stream.getvalue())
        self.assertTrue(stream.triggered)
        self.assertGreaterEqual(len(lines), 2)
        self.assertLessEqual(len(lines[-1]), 39)

    def test_wide_medium_and_narrow_rendering_priorities(self) -> None:
        cases = {
            120: (
                "Bootstrap replicates",
                "250/1000 replicates",
                "[",
            ),
            70: ("Bootstrap", "250/1000", "["),
            40: ("Bootstrap", "250/1000", "25%"),
            24: ("Bootstrap", "", "25%"),
        }
        for columns, (label, detail, prefix) in cases.items():
            with self.subTest(columns=columns):
                stream = TTYBuffer(columns)
                progress = progress_for(stream)
                progress.set_progress(
                    "bootstrap",
                    25,
                    "250/1000 replicates",
                )
                line = rendered_lines(stream.getvalue())[-1]
                self.assertIn(label, line)
                if detail:
                    self.assertIn(detail, line)
                self.assertTrue(line.startswith(prefix))

    def test_repeated_updates_at_same_integer_percent_are_suppressed(self) -> None:
        stream = TTYBuffer(80)
        progress = progress_for(stream)
        progress.start("bootstrap")
        initial_length = len(stream.getvalue())
        for replicate in range(1, 10):
            progress.set_progress(
                "bootstrap",
                replicate / 100,
                f"{replicate}/1000 replicates",
            )
        self.assertEqual(len(stream.getvalue()), initial_length)
        progress.set_progress("bootstrap", 1, "10/1000 replicates")
        self.assertGreater(len(stream.getvalue()), initial_length)

    def test_disabled_progress_writes_nothing(self) -> None:
        stream = TTYBuffer()
        with PipelineProgress(
            [("model", "ModelFinder")],
            enabled=False,
            stream=stream,
            terminal_width_provider=stream.terminal_width,
        ) as progress:
            progress.start("model")
            progress.set_progress("model", 75)
            progress.succeed("model")
        self.assertEqual(stream.getvalue(), "")


if __name__ == "__main__":
    unittest.main()
