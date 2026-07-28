from __future__ import annotations

import io
import sys
import tempfile
import time
import unittest
from pathlib import Path
from subprocess import CompletedProcess, TimeoutExpired
from unittest.mock import patch

from cladistica.progress import PipelineProgress
from cladistica.trees import (
    IQTreeProgressParser,
    MrBayesProgressParser,
    iqtree_cli_flags,
    iqtree_major_version,
    run_checked,
    run_tree_analyses,
    write_mrbayes_nexus,
)


class TreeRunnerTests(unittest.TestCase):
    def test_iqtree_version_detection_and_compatibility_flags(self) -> None:
        version_1 = CompletedProcess(
            ["iqtree", "--version"],
            0,
            stdout="IQ-TREE multicore version 1.6.12",
            stderr="",
        )
        with patch("cladistica.trees.subprocess.run", return_value=version_1):
            self.assertEqual(iqtree_major_version("iqtree"), 1)
            self.assertEqual(iqtree_cli_flags("iqtree"), ("-spp", "-nt", 1))

        version_3 = CompletedProcess(
            ["iqtree3", "--version"],
            0,
            stdout="IQ-TREE multicore version 3.0.1",
            stderr="",
        )
        with patch("cladistica.trees.subprocess.run", return_value=version_3):
            self.assertEqual(iqtree_cli_flags("iqtree3"), ("-p", "-T", 3))

    def test_iqtree_version_detection_failure_uses_current_flags(self) -> None:
        with patch(
            "cladistica.trees.subprocess.run",
            side_effect=TimeoutExpired(["iqtree", "--version"], 20),
        ):
            self.assertEqual(iqtree_major_version("iqtree"), 0)
            self.assertEqual(iqtree_cli_flags("iqtree"), ("-p", "-T", 0))

    def test_keyboard_interrupt_terminates_external_child_process(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            started = time.monotonic()

            def interrupt_after_first_line(line: str) -> None:
                raise KeyboardInterrupt

            with self.assertRaises(KeyboardInterrupt):
                run_checked(
                    [
                        sys.executable,
                        "-c",
                        "import time; print('ready', flush=True); time.sleep(30)",
                    ],
                    root,
                    root / "child.log",
                    on_line=interrupt_after_first_line,
                )
            self.assertLess(time.monotonic() - started, 6)

    def test_iqtree_log_parser_splits_model_bootstrap_and_ml_progress(self) -> None:
        progress = PipelineProgress(
            [
                ("model", "ModelFinder"),
                ("bootstrap", "Bootstrap replicates"),
                ("ml", "ML tree search"),
            ],
            enabled=True,
            stream=io.StringIO(),
        )
        parser = IQTreeProgressParser(progress, bootstrap=2)
        for line in (
            "Selecting individual models for 2 charsets using BIC...",
            "  1 GTR+F+I+G4 100.0 0.1 rbcL",
            "CPU time for ModelFinder: 1.0 seconds",
            "===> START BOOTSTRAP REPLICATE NUMBER 1",
            "TREE SEARCH COMPLETED AFTER 20 ITERATIONS",
            "===> START BOOTSTRAP REPLICATE NUMBER 2",
            "TREE SEARCH COMPLETED AFTER 22 ITERATIONS",
            "Consensus tree written to test.contree",
            "INITIALIZING CANDIDATE TREE SET",
            "Iteration 10 / LogL: -100.0",
            "TREE SEARCH COMPLETED AFTER 30 ITERATIONS",
            "FINALIZING TREE SEARCH",
            "Analysis results written to:",
        ):
            parser(line)
        parser.finish()
        self.assertEqual(progress.stages["model"].percent, 100)
        self.assertEqual(progress.stages["bootstrap"].percent, 100)
        self.assertEqual(progress.stages["ml"].percent, 100)

    def test_mrbayes_log_parser_uses_exact_generation_percentage(self) -> None:
        progress = PipelineProgress(
            [("bi", "MrBayes MCMC"), ("bi_summary", "BI summary")],
            enabled=True,
            stream=io.StringIO(),
        )
        parser = MrBayesProgressParser(progress, ngen=10_000)
        parser("      5000 -- (-100.0) * (-101.0)")
        self.assertEqual(progress.stages["bi"].percent, 50)
        parser("      Analysis completed in 7 seconds")
        parser("      Summarizing parameters in files run1.p and run2.p")
        parser('   Summarizing trees in files "run1.t" and "run2.t"')
        parser("   Calculating tree probabilities...")
        parser.finish()
        self.assertEqual(progress.stages["bi"].percent, 100)
        self.assertEqual(progress.stages["bi_summary"].percent, 100)

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
            self.assertIn("set autoclose=yes nowarn=yes;", nexus)
            self.assertIn("prset applyto=(all) ratepr=variable;", nexus)
            self.assertTrue((output_dir / "05_reports" / "alignment_statistics.tsv").exists())

    def test_single_partition_mrbayes_input_does_not_add_variable_rate_prior(
        self,
    ) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            fasta = root / "single.fasta"
            fasta.write_text(
                ">sampleA\nAACCGG\n>sampleB\nAATTGG\n",
                encoding="utf-8",
            )
            output = root / "single.nex"
            write_mrbayes_nexus(
                output,
                fasta,
                [{"marker": "rbcL", "start": 1, "end": 6, "length": 6}],
                ngen=10_000,
                nruns=2,
                nchains=4,
                burninfrac=0.25,
            )
            nexus = output.read_text(encoding="utf-8")
            self.assertIn("set autoclose=yes nowarn=yes;", nexus)
            self.assertNotIn("ratepr=variable", nexus)


if __name__ == "__main__":
    unittest.main()
