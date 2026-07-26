from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from cladistica.concat import concatenate_alignments
from cladistica.io import read_fasta


class ConcatenationTests(unittest.TestCase):
    def test_concatenates_and_fills_missing_partitions(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            input_dir = root / "input"
            input_dir.mkdir()
            (input_dir / "rbcL.fasta").write_text(">sampleA|marker=rbcL\nAAAA\n>sampleB|marker=rbcL\nCCCC\n", encoding="utf-8")
            (input_dir / "matK.fasta").write_text(">sampleA|marker=matK\nGG\n", encoding="utf-8")
            output_dir = root / "output"
            result = concatenate_alignments(input_dir=input_dir, output_dir=output_dir, markers=["rbcL", "matK"])
            records = read_fasta(output_dir / "concatenated_cpDNA.fasta")
            self.assertEqual(result.records_written, 2)
            self.assertEqual(records["sampleA"], "AAAAGG")
            self.assertEqual(records["sampleB"], "CCCC??")


if __name__ == "__main__":
    unittest.main()
