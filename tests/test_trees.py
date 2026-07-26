from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from cladistica.trees import run_tree_analyses


class TreeRunnerTests(unittest.TestCase):
    def test_dry_run_writes_reports_and_safe_partition_labels(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            input_dir = root / "input"
            input_dir.mkdir()
            (input_dir / "concatenated_cpDNA.fasta").write_text(">sampleA\nAACCGG\n>sampleB\nAATTGG\n", encoding="utf-8")
            (input_dir / "partition_coordinates.tsv").write_text(
                "marker\tstart\tend\tlength\ntrnL-F\t1\t2\t2\nrps4-trnS\t3\t6\t4\n",
                encoding="utf-8",
            )
            output_dir = root / "trees"
            result = run_tree_analyses(input_dir=input_dir, output_dir=output_dir, dry_run=True)
            nexus = (output_dir / "00_inputs" / "mrbayes_analysis.nex").read_text(encoding="utf-8")
            self.assertEqual(result.records_written, 0)
            self.assertIn("charset trnL_F", nexus)
            self.assertIn("charset rps4_trnS", nexus)
            self.assertTrue((output_dir / "05_reports" / "alignment_statistics.tsv").exists())


if __name__ == "__main__":
    unittest.main()
