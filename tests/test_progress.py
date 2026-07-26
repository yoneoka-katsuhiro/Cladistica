from __future__ import annotations

import io
import unittest

from cladistica.progress import PipelineProgress


class TTYBuffer(io.StringIO):
    def isatty(self) -> bool:
        return True


class ProgressTests(unittest.TestCase):
    def test_progress_bar_shows_percentage_and_short_stage_name(self) -> None:
        stream = TTYBuffer()
        progress = PipelineProgress([("model", "ModelFinder")], stream=stream)
        progress.set_progress("model", 75)
        self.assertIn(
            "[========================--------]  75% ModelFinder",
            stream.getvalue(),
        )

    def test_approximate_progress_has_tilde(self) -> None:
        stream = TTYBuffer()
        progress = PipelineProgress([("ml", "ML tree search")], stream=stream)
        progress.set_progress("ml", 25, "candidate trees", approximate=True)
        self.assertIn(
            "[========------------------------] ~25% ML tree search - candidate trees",
            stream.getvalue(),
        )

    def test_interactive_output_has_no_animation_or_ansi_cursor_control(self) -> None:
        stream = TTYBuffer()
        with PipelineProgress([("align", "MUSCLE alignment")], stream=stream) as progress:
            progress.start("align", "rbcL")
            progress.set_progress("align", 50, "1/2 markers")
            progress.succeed("align", "2 marker alignments")

        output = stream.getvalue()
        self.assertNotIn("\x1b", output)
        self.assertNotIn("Sequence stream", output)
        self.assertNotIn("+----------------------------------------+", output)
        self.assertFalse(hasattr(progress, "_thread"))
        self.assertTrue(output.endswith("\n"))

    def test_repeated_updates_at_same_integer_percent_are_suppressed(self) -> None:
        stream = TTYBuffer()
        progress = PipelineProgress([("bi", "MrBayes MCMC")], stream=stream)
        progress.start("bi")
        initial_length = len(stream.getvalue())
        for generation in range(1, 10):
            progress.set_progress(
                "bi",
                generation / 100,
                f"{generation}/10,000 generations",
            )
        self.assertEqual(len(stream.getvalue()), initial_length)
        progress.set_progress("bi", 1, "100/10,000 generations")
        self.assertGreater(len(stream.getvalue()), initial_length)

    def test_keyboard_interrupt_marks_active_stage_failed_without_animation(self) -> None:
        stream = TTYBuffer()
        progress = PipelineProgress([("align", "Alignment")], stream=stream)
        progress.start("align")
        progress.__exit__(KeyboardInterrupt, KeyboardInterrupt(), None)
        self.assertEqual(progress.stages["align"].state, "failed")
        self.assertIn("!0% Alignment", stream.getvalue())
        self.assertNotIn("\x1b", stream.getvalue())

    def test_non_tty_progress_uses_plain_status_lines(self) -> None:
        stream = io.StringIO()
        with PipelineProgress(
            [("align", "MUSCLE marker alignments"), ("ml", "ML tree")],
            stream=stream,
        ) as progress:
            progress.start("align", "rbcL: 4 sequences")
            progress.succeed("align", "1 marker alignment")
            progress.skip("ml", "Skipped by option")

        output = stream.getvalue()
        self.assertIn("[>>] MUSCLE marker alignments: rbcL: 4 sequences", output)
        self.assertIn("[OK] MUSCLE marker alignments: 1 marker alignment", output)
        self.assertIn("[--] ML tree: Skipped by option", output)

    def test_disabled_progress_writes_nothing(self) -> None:
        stream = TTYBuffer()
        with PipelineProgress([("model", "ModelFinder")], enabled=False, stream=stream) as progress:
            progress.start("model")
            progress.set_progress("model", 75)
            progress.succeed("model")
        self.assertEqual(stream.getvalue(), "")


if __name__ == "__main__":
    unittest.main()
