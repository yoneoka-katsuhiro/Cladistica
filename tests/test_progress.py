from __future__ import annotations

import io
import unittest

from cladistica.progress import PipelineProgress


class ProgressTests(unittest.TestCase):
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
