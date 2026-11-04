from __future__ import annotations

import tempfile
import unittest
import zipfile
from pathlib import Path
from unittest.mock import patch

from cladistica.package import package_run_outputs, write_zip


class PackageTests(unittest.TestCase):
    def test_existing_archive_cannot_destroy_original_input_before_packaging(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            source = root / "source"
            source.mkdir()
            original = root / "original.fasta"
            original.write_bytes(b">irreplaceable\nACGT\n")
            with self.assertRaises(FileExistsError):
                package_run_outputs(run_dir=source, output_dir=root / "output", archive_path=original)
            self.assertEqual(original.read_bytes(), b">irreplaceable\nACGT\n")
            self.assertFalse((root / "output").exists())

    def test_archive_inside_source_or_output_is_rejected_before_copying(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            source, output = root / "source", root / "output"
            source.mkdir()
            for archive in (source / "results.zip", output / "results.zip"):
                with self.subTest(archive=archive):
                    with self.assertRaisesRegex(ValueError, "outside"):
                        package_run_outputs(run_dir=source, output_dir=output, archive_path=archive)
                    self.assertFalse(output.exists())
                    self.assertFalse(archive.exists())

    def test_direct_zip_call_does_not_overwrite_existing_archive(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            source = root / "source"
            source.mkdir()
            archive = root / "original.zip"
            archive.write_bytes(b"previous archive")
            with self.assertRaises(FileExistsError):
                write_zip(source, archive)
            self.assertEqual(archive.read_bytes(), b"previous archive")

    def test_archive_created_during_compression_is_not_overwritten(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            source = root / "source"
            source.mkdir()
            (source / "data").write_bytes(b"data to compress")
            destination = root / "result.zip"
            original_write = zipfile.ZipFile.write

            def concurrent_writer(archive, *args, **kwargs):
                destination.write_bytes(b"a concurrent writer's result")
                return original_write(archive, *args, **kwargs)

            with patch("cladistica.package.zipfile.ZipFile.write", concurrent_writer):
                with self.assertRaises(FileExistsError):
                    write_zip(source, destination)
            self.assertEqual(destination.read_bytes(), b"a concurrent writer's result")
            self.assertEqual(set(root.iterdir()), {source, destination})

    def test_failed_compression_never_publishes_a_partial_archive(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            source = root / "source"
            source.mkdir()
            (source / "data").write_bytes(b"original data")
            destination = root / "result.zip"
            with patch("cladistica.package.zipfile.ZipFile.write", side_effect=OSError("disk full")):
                with self.assertRaisesRegex(OSError, "disk full"):
                    write_zip(source, destination)
            self.assertEqual(list(root.iterdir()), [source])
            self.assertEqual((source / "data").read_bytes(), b"original data")

    def test_packages_flat_results_and_complete_bi_bundle(self) -> None:
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
            self.assertEqual({path.name for path in output_dir.iterdir()}, expected | {"BI_resume"})
            self.assertEqual(result.warnings, [])
            self.assertIn(">B", (output_dir / "rbcL.fasta").read_text(encoding="utf-8"))
            with zipfile.ZipFile(archive_path) as archive:
                self.assertTrue(expected.issubset(archive.namelist()))
                self.assertIn("BI_resume/mrbayes_analysis.nex.run1.p", archive.namelist())


if __name__ == "__main__":
    unittest.main()
