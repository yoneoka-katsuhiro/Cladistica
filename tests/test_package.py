from __future__ import annotations

import tempfile
import unittest
import zipfile
from pathlib import Path

from cladistica.package import package_run_outputs


class PackageTests(unittest.TestCase):
    def test_packages_only_required_flat_files(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            run_dir = root / "run"

            accessions = run_dir / "accessions"
            (accessions / "fasta_by_marker").mkdir(parents=True)
            (accessions / "accession_all.csv").write_text("taxon,rbcL\nA,LC1\n", encoding="utf-8")
            (accessions / "accession_selected.csv").write_text("taxon,rbcL\nA,LC1\n", encoding="utf-8")
            (accessions / "fasta_by_marker" / "rbcL.fasta").write_text(">A\nACGT\n", encoding="utf-8")
            (accessions / "summly.txt").write_text("accessions\n", encoding="utf-8")
            (accessions / "run.log").write_text("accession log\n", encoding="utf-8")

            (run_dir / "alignments" / "intermediate" / "combined_by_marker").mkdir(parents=True)
            (run_dir / "alignments" / "intermediate" / "combined_by_marker" / "rbcL.fasta").write_text(
                ">A\nACGT\n>B\nACGA\n",
                encoding="utf-8",
            )
            (run_dir / "alignments" / "alignment_report.txt").write_text("alignment\n", encoding="utf-8")
            (run_dir / "concatenated").mkdir(parents=True)
            (run_dir / "concatenated" / "concatenated_cpDNA.fasta").write_text(">A\nACGT\n", encoding="utf-8")

            trees = run_dir / "trees"
            (trees / "00_inputs").mkdir(parents=True)
            (trees / "00_inputs" / "iqtree_partitions.txt").write_text("DNA, rbcL = 1-4\n", encoding="utf-8")
            (trees / "00_inputs" / "mrbayes_analysis.nex").write_text("#NEXUS\n", encoding="utf-8")
            (trees / "02_iqtree_ml").mkdir(parents=True)
            (trees / "02_iqtree_ml" / "cladistica_ml.treefile").write_text("(A);\n", encoding="utf-8")
            (trees / "03_mrbayes").mkdir(parents=True)
            (trees / "03_mrbayes" / "mrbayes_analysis.nex.con.tre").write_text("(A);\n", encoding="utf-8")
            (trees / "03_mrbayes" / "mrbayes_analysis.nex.run1.p").write_text("run1\n", encoding="utf-8")
            (trees / "03_mrbayes" / "mrbayes_analysis.nex.run2.p").write_text("run2\n", encoding="utf-8")
            (trees / "05_reports").mkdir(parents=True)
            (trees / "05_reports" / "model_selection_summary.txt").write_text("models\n", encoding="utf-8")

            output_dir = root / "deliverable"
            archive_path = root / "test.zip"
            result = package_run_outputs(run_dir=run_dir, output_dir=output_dir, archive_path=archive_path)

            expected = {
                "accession_all.csv",
                "accession_selected.csv",
                "rbcL.fasta",
                "concatenated.fasta",
                "partitions.txt",
                "BI.nex",
                "ML.tre",
                "BI.tre",
                "run1.p",
                "run2.p",
                "summly.txt",
                "run.log",
            }
            self.assertEqual(result.records_written, len(expected))
            self.assertEqual({path.name for path in output_dir.iterdir()}, expected)
            self.assertEqual(result.warnings, [])
            self.assertIn(">B", (output_dir / "rbcL.fasta").read_text(encoding="utf-8"))
            with zipfile.ZipFile(archive_path) as archive:
                self.assertEqual(set(archive.namelist()), expected)


if __name__ == "__main__":
    unittest.main()
