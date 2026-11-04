from __future__ import annotations

import json
import os
from pathlib import Path
import shutil
import tempfile
import unittest

from cladistica.bayes import (
    ANALYSIS, MANIFEST, checkpoint_generation, resume_bayes, sample_paths,
    save_manifest, sha256, start_bayes, verify_bundle, last_parameter_generation,
)
from cladistica.package import package_run_outputs
from cladistica.trees import run_checked, write_mrbayes_nexus
from cladistica.workflow import durable_analysis_workspace


def make_bundle(directory: Path, nruns: int = 2) -> dict:
    directory.mkdir()
    (directory / ANALYSIS).write_text("#NEXUS\n  mcmcp ngen=1000;\n  mcmc;\n")
    (directory / (ANALYSIS + ".ckp")).write_text("#NEXUS\n[ID: 123]\n[generation: 1000]\n")
    for path in sample_paths(directory, nruns):
        path.write_text("[ID: 123]\nGen LnL\n0 -40\n1000 -30\n" if path.suffix == ".p" else "#NEXUS\n[ID: 123]\nbegin trees;\nend;\n")
    state = dict(schema=1, status="complete", ngen=1000, nruns=nruns,
                 mrbayes_version="3.2.7a", analysis_sha256=sha256(directory / ANALYSIS), files={}, history=[])
    save_manifest(directory, state)
    return state


class BayesResumeTests(unittest.TestCase):
    def test_checkpoint_and_all_runs_are_required(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary) / "bundle"
            state = make_bundle(root, nruns=3)
            self.assertEqual(verify_bundle(root, state), 1000)
            (root / (ANALYSIS + ".run3.t")).unlink()
            with self.assertRaisesRegex(ValueError, "samples are missing"):
                verify_bundle(root, state)

    def test_single_run_uses_mrbayes_single_run_basenames(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary) / "bundle"
            state = make_bundle(root, nruns=1)
            self.assertEqual(verify_bundle(root, state), 1000)

    def test_changed_input_is_rejected_before_launch(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary) / "bundle"
            state = make_bundle(root)
            (root / ANALYSIS).write_text("different alignment/model")
            with self.assertRaisesRegex(ValueError, "input/model/settings have changed"):
                verify_bundle(root, state)

    def test_changed_sample_hash_is_rejected(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary) / "bundle"
            state = make_bundle(root)
            sample = root / (ANALYSIS + ".run1.p")
            state["files"] = {sample.name: sha256(sample)}
            sample.write_text(sample.read_text() + "2000 -20\n")
            with self.assertRaisesRegex(ValueError, "file is missing or changed"):
                verify_bundle(root, state)

    def test_sample_from_another_run_is_rejected(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary) / "bundle"
            state = make_bundle(root)
            path = root / (ANALYSIS + ".run2.t")
            path.write_text(path.read_text().replace("123", "456"))
            with self.assertRaisesRegex(ValueError, "IDs do not match"):
                verify_bundle(root, state)

    def test_missing_checkpoint_cannot_resume_from_parameter_files(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary) / "bundle"
            state = make_bundle(root)
            (root / (ANALYSIS + ".ckp")).unlink()
            with self.assertRaisesRegex(ValueError, "Checkpoint .ckp is missing"):
                verify_bundle(root, state)

    def test_target_must_exceed_checkpoint(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            make_bundle(root / "bundle")
            with self.assertRaisesRegex(ValueError, "total target"):
                resume_bayes(bundle=root / "bundle", output_dir=root / "result", ngen=1000)
            self.assertFalse((root / "result").exists())

    def test_interruption_keeps_workspace(self):
        with tempfile.TemporaryDirectory() as temporary:
            output = Path(temporary) / "result"
            with self.assertRaises(KeyboardInterrupt):
                with durable_analysis_workspace(output) as work:
                    (work / "checkpoint").write_text("keep me")
                    raise KeyboardInterrupt
            self.assertEqual((output / "_work/checkpoint").read_text(), "keep me")
            with self.assertRaisesRegex(ValueError, "not empty"):
                with durable_analysis_workspace(output):
                    pass

    def test_packaging_cannot_erase_source(self):
        with tempfile.TemporaryDirectory() as temporary:
            source = Path(temporary) / "source"
            source.mkdir()
            valuable = source / "checkpoint.ckp"
            valuable.write_text("keep me")
            with self.assertRaisesRegex(ValueError, "must not replace"):
                package_run_outputs(run_dir=source, output_dir=source)
            self.assertTrue(valuable.exists())

    def test_packaging_preserves_entire_bi_directory(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            trees = root / "source/trees"
            trees.mkdir(parents=True)
            make_bundle(trees / "03_mrbayes")
            (trees / "03_mrbayes" / (ANALYSIS + ".mcmc")).write_text("diagnostics")
            package_run_outputs(run_dir=root / "source", output_dir=root / "result")
            for path in (trees / "03_mrbayes").iterdir():
                self.assertEqual(path.read_bytes(), (root / "result/BI_resume" / path.name).read_bytes())


@unittest.skipUnless(os.environ.get("CLADISTICA_TEST_MRBAYES") == "1" and shutil.which("mb"),
                     "Set CLADISTICA_TEST_MRBAYES=1 with MrBayes installed to run real integration tests")
class RealMrBayesTests(unittest.TestCase):
    def make_analysis(self, root: Path, target: int, nruns: int = 2) -> Path:
        matrix = Path(__file__).resolve().parents[1] / "examples/tiny.fasta"
        analysis = root / "input.nex"
        write_mrbayes_nexus(analysis, matrix, [{"marker": "demo", "start": 1, "end": 32}],
                           ngen=target, nruns=nruns, nchains=2, burninfrac=0.25)
        return analysis

    def assert_continued(self, source: Path, output: Path, target: int):
        before = {p.name: sha256(p) for p in source.iterdir() if p.is_file()}
        resume_bayes(bundle=source, output_dir=output, ngen=target)
        self.assertEqual(before, {p.name: sha256(p) for p in source.iterdir() if p.is_file()})
        for path in sample_paths(output / "BI_resume", 2):
            if path.suffix == ".p":
                self.assertEqual(last_parameter_generation(path), target)
                generations = [int(line.split()[0]) for line in path.read_text().splitlines()
                               if line.split() and line.split()[0].isdigit()]
                self.assertEqual(generations, list(range(0, target + 1, 1000)))
        self.assertIn("Using samples up to generation", (output / "BI_resume/run_mrbayes.log").read_text())

    def test_completed_run_extends_without_losing_old_samples(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            start_bayes(root / "initial", self.make_analysis(root, 2000), "mb", ngen=2000, nruns=2, runner=run_checked)
            self.assert_continued(root / "initial", root / "extended", 4000)

    def test_interrupted_run_resumes_from_checkpoint(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            def interrupt(line):
                import re
                # Stop after a generation report, when the preceding checkpoint
                # has been written, rather than reading a file during replacement.
                match = re.match(r"\s*(\d+)\s+--", line)
                if match and int(match.group(1)) >= 10000:
                    raise KeyboardInterrupt
            with self.assertRaises(KeyboardInterrupt):
                start_bayes(root / "interrupted", self.make_analysis(root, 1000000), "mb",
                            ngen=1000000, nruns=2, runner=run_checked, on_line=interrupt)
            state = json.loads((root / "interrupted" / MANIFEST).read_text())
            self.assertEqual(state["status"], "interrupted")
            target = checkpoint_generation(root / "interrupted") + 3000
            self.assert_continued(root / "interrupted", root / "recovered", target)


if __name__ == "__main__":
    unittest.main()
