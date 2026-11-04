"""Render actual output from the bundled synthetic example (requires matplotlib)."""
from __future__ import annotations

import json
import io
import re
import os
from pathlib import Path
import sys
import tempfile

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
os.environ.setdefault("MPLCONFIGDIR", str(ROOT / ".matplotlib-cache"))
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from Bio import Phylo
from cladistica.bayes import ANALYSIS, resume_bayes, sha256
from cladistica.workflow import run_from_fasta


def parameter_trace(path: Path):
    rows = []
    for line in path.read_text().splitlines():
        fields = line.split()
        if fields and fields[0].isdigit():
            rows.append((int(fields[0]), float(fields[1])))
    return rows


def main():
    with tempfile.TemporaryDirectory(prefix="cladistica_demo_") as temporary:
        root = Path(temporary)
        initial = root / "initial"
        extended = root / "extended"
        run_from_fasta(project_dir=root, fasta_dir=None, additional_fasta_dir=None,
                       concatenated_fasta=ROOT / "examples/tiny.fasta", partition_file=None,
                       output_dir=initial, ngen=2000, skip_ml=True)
        resume_bayes(bundle=initial, output_dir=extended, ngen=4000)
        bundle = extended / "BI_resume"
        nexus = (bundle / (ANALYSIS + ".con.tre")).read_text()
        # The plot uses topology/branch lengths only; remove MrBayes extended
        # annotation comments that Bio.Nexus cannot parse. Raw output is unchanged.
        tree = Phylo.read(io.StringIO(re.sub(r"\[[^\]]*\]", "", nexus)), "nexus")
        tree.ladderize()
        plt.rcParams.update({"font.size": 11, "font.family": "DejaVu Sans"})
        figure, axes = plt.subplots(1, 2, figsize=(13, 4.5), gridspec_kw={"width_ratios": [1, 1.3]})
        figure.patch.set_facecolor("#f5f7f8")
        for axis in axes:
            axis.set_facecolor("#f5f7f8")
        Phylo.draw(tree, axes=axes[0], do_show=False, show_confidence=False)
        axes[0].set_title("", loc="center")
        axes[0].set_title("Consensus tree · synthetic data", loc="left", weight="bold", pad=20)
        axes[0].set_ylabel("")
        axes[0].set_yticks([])
        axes[0].set_xlabel("Substitutions per site")
        for run, color in ((1, "#176b8c"), (2, "#c36939")):
            values = parameter_trace(bundle / f"{ANALYSIS}.run{run}.p")
            axes[1].plot(*zip(*values), marker="o", color=color, label=f"Run {run}", linewidth=2)
        axes[1].axvline(2000, color="#60666d", linestyle="--", linewidth=1)
        axes[1].set_title("Samples retained across continuation", loc="left", weight="bold", pad=20)
        axes[1].set_xlabel("Generation · 2,000 → 4,000 total")
        axes[1].set_ylabel("Log likelihood")
        axes[1].legend(frameon=False)
        axes[1].grid(axis="y", alpha=.2)
        for axis in axes:
            axis.spines[["top", "right"]].set_visible(False)
        figure.text(.03, .02, "4 synthetic taxa · 32 sites · short software demonstration, not a biological result or convergence assessment", color="#505860", fontsize=10)
        figure.tight_layout(rect=(0, .07, 1, 1))
        images = ROOT / "docs/images"
        images.mkdir(parents=True, exist_ok=True)
        figure.savefig(images / "demo-output.png", dpi=160, facecolor=figure.get_facecolor())
        state = json.loads((bundle / "resume_manifest.json").read_text())
        (images / "demo-provenance.json").write_text(json.dumps({
            "dataset": "examples/tiny.fasta", "synthetic": True,
            "input_sha256": sha256(ROOT / "examples/tiny.fasta"),
            "cladistica_version": state["cladistica_version"], "mrbayes_version": state["mrbayes_version"],
            "initial_generations": 2000, "total_generations": 4000,
            "sample_generations": [0,1000,2000,3000,4000],
            "created_at": state["updated_at"],
            "purpose": "Software demonstration only; no biological inference or convergence claim.",
        }, indent=2) + "\n")


if __name__ == "__main__":
    main()
