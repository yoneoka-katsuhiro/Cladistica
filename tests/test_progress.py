from __future__ import annotations

import io
import unittest

from cladistica.cli import build_parser
from cladistica.progress import (
    ANIMATION_INTERVAL_SECONDS,
    DEFAULT_STAGE_ITEMS,
    STREAM_HEIGHT,
    STREAM_WIDTH,
    PipelineProgress,
    StreamLane,
    normalize_stream_item,
    sequence_windows,
)


class TTYBuffer(io.StringIO):
    def isatty(self) -> bool:
        return True


class ProgressTests(unittest.TestCase):
    def test_demo_command_can_preview_detected_failure(self) -> None:
        args = build_parser().parse_args(["demo", "--fail", "--duration", "3"])
        self.assertTrue(args.fail)
        self.assertEqual(args.duration, 3)

    def test_stream_uses_a_fixed_five_by_forty_canvas(self) -> None:
        progress = PipelineProgress([("align", "Alignment")], stream=io.StringIO())
        progress.feed("Hymenasplenium hondoense")
        progress._initialize_lanes()
        canvas = progress._stream_canvas()
        self.assertEqual(len(canvas), STREAM_HEIGHT)
        self.assertTrue(all(len(row) == STREAM_WIDTH for row in canvas))

    def test_stage_change_replaces_lanes_with_stage_specific_vocabulary(self) -> None:
        progress = PipelineProgress(
            [("accessions", "Accessions"), ("model", "ModelFinder")],
            stream=TTYBuffer(),
        )
        progress.feed("Hymenasplenium hondoense")
        progress._initialize_lanes()
        progress.start("model", "testing models")
        self.assertEqual(progress._active_key, "model")
        self.assertEqual(progress._lanes[0].text, "ModelFinder")
        self.assertIn("GTR", DEFAULT_STAGE_ITEMS["model"])
        self.assertIn("DDBJ", DEFAULT_STAGE_ITEMS["accessions"])

    def test_progress_bar_shows_percentage_and_short_stage_name(self) -> None:
        stream = TTYBuffer()
        progress = PipelineProgress([("model", "ModelFinder")], stream=stream)
        progress._initialize_lanes()
        progress.set_progress("model", 75)
        self.assertIn(
            "[========================--------]  75% ModelFinder",
            stream.getvalue(),
        )

    def test_text_moves_from_left_to_right(self) -> None:
        progress = PipelineProgress([("align", "Alignment")], stream=io.StringIO())
        progress._lanes = [
            StreamLane("Hymenasplenium hondoense", row - 8)
            for row in range(STREAM_HEIGHT)
        ]
        positions_before = [lane.x for lane in progress._lanes]
        progress._advance_stream()
        self.assertEqual(
            [lane.x for lane in progress._lanes],
            [position + 1 for position in positions_before],
        )

    def test_sequence_windows_show_real_bases_without_loading_whole_alignment(self) -> None:
        sequence = "A" * 40 + "C" * 40 + "G" * 40 + "T" * 40
        windows = sequence_windows(sequence)
        self.assertEqual(len(windows), 3)
        self.assertTrue(all(len(window) == STREAM_WIDTH for window in windows))
        self.assertEqual(sequence_windows("ACGTN"), ["ACGTN"])

    def test_failure_breaks_visible_text_into_falling_particles(self) -> None:
        progress = PipelineProgress(
            [("align", "Alignment")],
            stream=TTYBuffer(),
        )
        progress._lanes = [
            StreamLane(f"Taxon_{row}", row * 3)
            for row in range(STREAM_HEIGHT)
        ]
        visible_characters = sum(
            character != " "
            for row in progress._stream_canvas()
            for character in row
        )
        progress.fail("align", "process stopped")
        self.assertTrue(progress._falling)
        for _ in range(20):
            progress._frame += 1
            progress._advance_fall()
        fallen_characters = sum(
            character != " "
            for row in progress._particles or []
            for character in row
        )
        bottom_characters = sum(
            character != " " for character in (progress._particles or [[]])[-1]
        )
        self.assertEqual(fallen_characters, visible_characters)
        self.assertGreater(bottom_characters, 0)

    def test_keyboard_interrupt_also_starts_character_collapse(self) -> None:
        progress = PipelineProgress(
            [("align", "Alignment")],
            stream=TTYBuffer(),
        )
        progress._initialize_lanes()
        progress.start("align")
        progress.__exit__(KeyboardInterrupt, KeyboardInterrupt(), None)
        self.assertTrue(progress._falling)
        self.assertEqual(progress.stages["align"].state, "failed")

    def test_animation_rate_is_lightweight(self) -> None:
        self.assertGreaterEqual(ANIMATION_INTERVAL_SECONDS, 0.08)
        self.assertLessEqual(ANIMATION_INTERVAL_SECONDS, 0.12)

    def test_stream_item_removes_terminal_control_characters(self) -> None:
        self.assertEqual(
            normalize_stream_item("\x1b Hymenasplenium   hondoense "),
            "Hymenasplenium hondoense",
        )

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


if __name__ == "__main__":
    unittest.main()
