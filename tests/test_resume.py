from __future__ import annotations

import shutil
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from cladistica.download import validate_selected_accessions
from cladistica.io import read_fasta
from cladistica.models import PipelineResult
from cladistica.workflow import (
    parse_partition_file,
    run_accession_survey,
    run_from_fasta,
    run_from_selected_accessions,
)


def fake_tree_analyses(*, input_dir: Path, output_dir: Path, skip_ml: bool, skip_bi: bool, **kwargs: object) -> PipelineResult:
    fasta = next(input_dir.glob("*.fasta"))
    records = read_fasta(fasta)
    length = len(next(iter(records.values())))
    inputs = output_dir / "00_inputs"
    inputs.mkdir(parents=True)
    (inputs / "iqtree_partitions.txt").write_text(f"DNA, cpDNA = 1-{length}\n", encoding="utf-8")
    (inputs / "mrbayes_analysis.nex").write_text("#NEXUS\n", encoding="utf-8")
    reports = output_dir / "05_reports"
    reports.mkdir(parents=True)
    (reports / "analysis_summary.txt").write_text("tree summary\n", encoding="utf-8")
    (reports / "model_selection_summary.txt").write_text("Best-fit model according to BIC: TEST\n", encoding="utf-8")
    if not skip_ml:
        ml = output_dir / "02_iqtree_ml"
        ml.mkdir(parents=True)
        (ml / "cladistica_ml.treefile").write_text("(A,B);\n", encoding="utf-8")
    if not skip_bi:
        bi = output_dir / "03_mrbayes"
        bi.mkdir(parents=True)
        (bi / "mrbayes_analysis.nex").write_text("#NEXUS\n", encoding="utf-8")
        (bi / "mrbayes_analysis.nex.con.tre").write_text("(A,B);\n", encoding="utf-8")
        (bi / "mrbayes_analysis.nex.run1.p").write_text("run1\n", encoding="utf-8")
        (bi / "mrbayes_analysis.nex.run2.p").write_text("run2\n", encoding="utf-8")
    return PipelineResult(str(output_dir), records_written=int(not skip_ml) + int(not skip_bi))


def copy_alignment(input_fasta: Path, output_fasta: Path, command: str) -> None:
    output_fasta.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(input_fasta, output_fasta)


class ResumeWorkflowTests(unittest.TestCase):
    def test_resume_from_concatenated_fasta_without_accession_outputs(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            alignment = root / "matrix.fasta"
            alignment.write_text(">A\nAACCGG\n>B\nAATTGG\n", encoding="utf-8")
            partitions = root / "partitions.txt"
            partitions.write_text("DNA, rbcL = 1-3\nDNA, matK = 4-6\n", encoding="utf-8")
            output = root / "output"

            with patch("cladistica.workflow.run_tree_analyses", side_effect=fake_tree_analyses):
                result = run_from_fasta(
                    project_dir=root,
                    fasta_dir=None,
                    additional_fasta_dir=None,
                    concatenated_fasta=alignment,
                    partition_file=partitions,
                    output_dir=output,
                    markers=None,
                )

            self.assertEqual(result.records_written, 2)
            self.assertTrue((output / "concatenated.fasta").exists())
            self.assertTrue((output / "ML.tre").exists())
            self.assertTrue((output / "run1.p").exists())
            self.assertTrue((output / "run2.p").exists())
            self.assertFalse((output / "accession_selected.csv").exists())
            self.assertEqual(result.warnings, [])

    def test_resume_from_marker_fasta_merges_additional_samples(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            primary = root / "primary"
            additional = root / "additional"
            primary.mkdir()
            additional.mkdir()
            (primary / "rbcL.fasta").write_text(">A|taxon=A\nAACCGG\n", encoding="utf-8")
            (additional / "rbcL.fa").write_text(">B|taxon=B\nAATTGG\n", encoding="utf-8")
            output = root / "output"

            with patch("cladistica.align.run_muscle", side_effect=copy_alignment), patch(
                "cladistica.workflow.run_tree_analyses",
                side_effect=fake_tree_analyses,
            ):
                result = run_from_fasta(
                    project_dir=root,
                    fasta_dir=primary,
                    additional_fasta_dir=additional,
                    concatenated_fasta=None,
                    partition_file=None,
                    output_dir=output,
                    markers=["rbcL"],
                    muscle_command="true",
                )

            merged = read_fasta(output / "rbcL.fasta")
            self.assertEqual(result.records_written, 2)
            self.assertEqual({header.split("|", 1)[0] for header in merged}, {"A", "B"})
            self.assertTrue((output / "concatenated.fasta").exists())

    def test_resume_from_hand_selected_accessions(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            selected = root / "accession_selected.csv"
            selected.write_text("taxon,sample_id,rbcL\nTaxon A,A,AA1\nTaxon B,B,BB1\n", encoding="utf-8")
            all_accessions = root / "accession_all.csv"
            all_accessions.write_text("taxon,sample_id,rbcL\nTaxon A,A,AA1\nTaxon B,B,BB1\n", encoding="utf-8")
            output = root / "output"

            def fake_download(*, output_dir: Path, **kwargs: object) -> PipelineResult:
                output_dir.mkdir(parents=True)
                (output_dir / "rbcL.fasta").write_text(">A|taxon=Taxon_A\nAACCGG\n>B|taxon=Taxon_B\nAATTGG\n", encoding="utf-8")
                (output_dir / "summly.txt").write_text("download\n", encoding="utf-8")
                (output_dir / "run.log").write_text("download log\n", encoding="utf-8")
                return PipelineResult(str(output_dir), records_written=2)

            with patch(
                "cladistica.workflow.download_fasta_from_accession_table",
                side_effect=fake_download,
            ), patch("cladistica.align.run_muscle", side_effect=copy_alignment), patch(
                "cladistica.workflow.run_tree_analyses",
                side_effect=fake_tree_analyses,
            ):
                result = run_from_selected_accessions(
                    project_dir=root,
                    accession_selected=selected,
                    accession_all=all_accessions,
                    additional_fasta_dir=None,
                    output_dir=output,
                    email="researcher@example.org",
                    markers=["rbcL"],
                    muscle_command="true",
                )

            self.assertEqual(result.records_written, 2)
            self.assertTrue((output / "accession_all.csv").exists())
            self.assertTrue((output / "accession_selected.csv").exists())
            self.assertTrue((output / "ML.tre").exists())
            self.assertTrue((output / "run2.p").exists())

    def test_accession_survey_stops_after_all_table(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            output = root / "output"

            def fake_accessions(*, output_dir: Path, **kwargs: object) -> PipelineResult:
                output_dir.mkdir(parents=True)
                (output_dir / "accession_all.csv").write_text(
                    "taxon,sample_id,rbcL,selection_status\nTaxon A,A,AA1,selected\n",
                    encoding="utf-8",
                )
                (output_dir / "summly.txt").write_text("selection\n", encoding="utf-8")
                (output_dir / "run.log").write_text("selection log\n", encoding="utf-8")
                (output_dir / "genbank_query_log.tsv").write_text("query\tcount\nA\t1\n", encoding="utf-8")
                return PipelineResult(str(output_dir), records_written=1)

            with patch("cladistica.workflow.run_accession_pipeline", side_effect=fake_accessions):
                result = run_accession_survey(
                    project_dir=root,
                    genus="Testgenus",
                    outgroups=[],
                    output_dir=output,
                    email="researcher@example.org",
                    markers=["rbcL"],
                )

            self.assertEqual(result.records_written, 1)
            self.assertEqual(
                {path.name for path in output.iterdir()},
                {"accession_all.csv", "summly.txt", "run.log"},
            )

    def test_partition_parser_rejects_out_of_range_coordinates(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "partitions.txt"
            path.write_text("DNA, rbcL = 1-100\n", encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "Invalid partition coordinates"):
                parse_partition_file(path, alignment_length=10)

    def test_selected_table_rejects_multiple_accessions_per_marker(self) -> None:
        rows = [{"taxon": "Taxon A", "sample_id": "A", "rbcL": "AA1;AA2"}]
        with self.assertRaisesRegex(ValueError, "multiple accessions"):
            validate_selected_accessions(rows, ["rbcL"])


if __name__ == "__main__":
    unittest.main()
