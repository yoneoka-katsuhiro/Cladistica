from __future__ import annotations

import shutil
import subprocess
from collections import Counter
from pathlib import Path
import re

from .io import read_fasta, read_tsv, write_tsv
from .models import PipelineResult
from .progress import PipelineProgress


def detect_command(candidates: list[str], explicit: str = "") -> str | None:
    if explicit:
        return explicit if shutil.which(explicit) else None
    for command in candidates:
        found = shutil.which(command)
        if found:
            return found
    return None


def alignment_statistics(fasta_path: Path) -> dict[str, object]:
    records = read_fasta(fasta_path)
    if not records:
        raise ValueError(f"No FASTA records found: {fasta_path}")
    lengths = {len(sequence) for sequence in records.values()}
    if len(lengths) != 1:
        raise ValueError(f"Concatenated FASTA has inconsistent sequence lengths: {fasta_path}")
    length = next(iter(lengths))
    variable = 0
    parsimony = 0
    for index in range(length):
        states = [sequence[index].upper() for sequence in records.values() if sequence[index].upper() in {"A", "C", "G", "T"}]
        counts = Counter(states)
        if len(counts) > 1:
            variable += 1
        if sum(1 for count in counts.values() if count >= 2) >= 2:
            parsimony += 1
    return {
        "samples": len(records),
        "alignment_length": length,
        "variable_sites": variable,
        "parsimony_informative_sites": parsimony,
    }


def read_partitions(input_dir: Path, fasta_path: Path) -> list[dict[str, object]]:
    partition_path = input_dir / "partition_coordinates.tsv"
    if partition_path.exists():
        rows, _ = read_tsv(partition_path)
        return rows
    stats = alignment_statistics(fasta_path)
    return [{"marker": "cpDNA", "start": 1, "end": stats["alignment_length"], "length": stats["alignment_length"]}]


def write_iqtree_partition_file(path: Path, partitions: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = [f"DNA, {partition_label(str(row['marker']))} = {row['start']}-{row['end']}" for row in partitions]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def partition_label(value: str) -> str:
    label = re.sub(r"[^A-Za-z0-9_]+", "_", value).strip("_")
    return label or "partition"


def nexus_taxon_label(value: str) -> str:
    if re.fullmatch(r"[A-Za-z0-9_.-]+", value):
        return value
    return "'" + value.replace("'", "''") + "'"


def write_mrbayes_nexus(
    path: Path,
    fasta_path: Path,
    partitions: list[dict[str, object]],
    ngen: int,
    nruns: int,
    nchains: int,
    burninfrac: float,
) -> None:
    records = read_fasta(fasta_path)
    alignment_length = len(next(iter(records.values())))
    lines = [
        "#NEXUS",
        "",
        "begin data;",
        f"  dimensions ntax={len(records)} nchar={alignment_length};",
        "  format datatype=dna missing=? gap=-;",
        "  matrix",
    ]
    for sample, sequence in records.items():
        lines.append(f"  {nexus_taxon_label(sample)} {sequence}")
    lines.extend(["  ;", "end;", "", "begin mrbayes;"])
    for row in partitions:
        lines.append(f"  charset {partition_label(str(row['marker']))} = {row['start']}-{row['end']};")
    if len(partitions) > 1:
        names = ", ".join(partition_label(str(row["marker"])) for row in partitions)
        lines.append(f"  partition cpDNA = {len(partitions)}: {names};")
        lines.append("  set partition=cpDNA;")
        lines.append("  unlink statefreq=(all) revmat=(all) shape=(all);")
    lines.extend(
        [
            "  lset nst=6 rates=gamma;",
            "  prset statefreqpr=dirichlet(1,1,1,1);",
            f"  mcmcp ngen={ngen} nruns={nruns} nchains={nchains} samplefreq=1000 printfreq=1000 diagnfreq=5000 burninfrac={burninfrac};",
            "  mcmc;",
            "  sump;",
            "  sumt;",
            "end;",
        ]
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def run_checked(command: list[str], cwd: Path, log_path: Path) -> None:
    log_path.parent.mkdir(parents=True, exist_ok=True)
    with log_path.open("w", encoding="utf-8") as log:
        result = subprocess.run(command, cwd=cwd, text=True, stdout=log, stderr=subprocess.STDOUT, check=False)
    if result.returncode != 0:
        raise RuntimeError(f"Command failed ({result.returncode}): {' '.join(command)}. See {log_path}")


def parse_iqtree_selected_models(iqtree_report: Path) -> list[dict[str, object]]:
    if not iqtree_report.exists():
        return []
    rows: list[dict[str, object]] = []
    in_table = False
    for raw in iqtree_report.read_text(encoding="utf-8", errors="replace").splitlines():
        line = raw.strip()
        if line.startswith("List of best-fit models per partition"):
            in_table = True
            continue
        if not in_table:
            continue
        if not line or line.startswith("ID ") or line.startswith("AIC"):
            continue
        if line.startswith("SUBSTITUTION PROCESS"):
            break
        parts = line.split()
        if len(parts) < 10 or not parts[0].isdigit():
            continue
        rows.append(
            {
                "partition_id": parts[0],
                "model": parts[1],
                "log_likelihood": parts[2],
                "aic": parts[3],
                "aicc": parts[6],
                "bic": parts[9],
            }
        )
    return rows


def write_model_selection_summary(iqtree_report: Path, output_dir: Path) -> None:
    rows = parse_iqtree_selected_models(iqtree_report)
    report_dir = output_dir / "05_reports"
    write_tsv(
        report_dir / "selected_models.tsv",
        rows,
        ["partition_id", "model", "log_likelihood", "aic", "aicc", "bic"],
    )
    best_line = ""
    if iqtree_report.exists():
        for raw in iqtree_report.read_text(encoding="utf-8", errors="replace").splitlines():
            if raw.startswith("Best-fit model according to BIC:"):
                best_line = raw.strip()
                break
    lines = [
        "Cladistica model selection summary",
        "",
        best_line or "Best-fit model according to BIC: unavailable",
        "",
        "Per-partition models:",
    ]
    for row in rows:
        lines.append(f"- partition {row['partition_id']}: {row['model']} (BIC {row['bic']})")
    (report_dir / "model_selection_summary.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")


def copy_tree_outputs(work_dir: Path, output_dir: Path) -> None:
    tree_dir = output_dir / "04_trees"
    tree_dir.mkdir(parents=True, exist_ok=True)
    for path in work_dir.rglob("*"):
        if tree_dir in path.parents:
            continue
        if path.suffix in {".treefile", ".contree", ".tre", ".nex"} or path.name.endswith(".con.tre"):
            destination = tree_dir / path.name
            if path.is_file():
                destination.write_bytes(path.read_bytes())


def run_tree_analyses(
    *,
    input_dir: Path,
    output_dir: Path,
    threads: str = "AUTO",
    bootstrap: int = 1000,
    skip_ml: bool = False,
    skip_bi: bool = False,
    iqtree_command: str = "",
    mrbayes_command: str = "",
    ngen: int = 1_000_000,
    nruns: int = 2,
    nchains: int = 4,
    burninfrac: float = 0.25,
    dry_run: bool = False,
    progress: PipelineProgress | None = None,
) -> PipelineResult:
    fasta_candidates = sorted(input_dir.glob("*.fasta")) + sorted(input_dir.glob("*.fa")) + sorted(input_dir.glob("*.fas"))
    if not fasta_candidates:
        raise FileNotFoundError(f"No concatenated FASTA file was found in {input_dir}")
    fasta_path = fasta_candidates[0]
    stats = alignment_statistics(fasta_path)
    partitions = read_partitions(input_dir, fasta_path)
    write_tsv(output_dir / "05_reports" / "alignment_statistics.tsv", [stats], ["samples", "alignment_length", "variable_sites", "parsimony_informative_sites"])
    write_tsv(output_dir / "05_reports" / "partition_statistics.tsv", partitions, ["marker", "start", "end", "length"])
    write_tsv(
        output_dir / "05_reports" / "parameters.tsv",
        [
            {"parameter": "threads", "value": threads},
            {"parameter": "bootstrap", "value": bootstrap},
            {"parameter": "skip_ml", "value": skip_ml},
            {"parameter": "skip_bi", "value": skip_bi},
            {"parameter": "ngen", "value": ngen},
            {"parameter": "nruns", "value": nruns},
            {"parameter": "nchains", "value": nchains},
            {"parameter": "burninfrac", "value": burninfrac},
            {"parameter": "dry_run", "value": dry_run},
        ],
        ["parameter", "value"],
    )
    iqtree = detect_command(["iqtree3", "iqtree2", "iqtree"], iqtree_command)
    mrbayes = detect_command(["mb"], mrbayes_command)
    if not skip_ml and not iqtree and not dry_run:
        raise RuntimeError("IQ-TREE was not found. Install iqtree3/iqtree2/iqtree or use --skip-ml.")
    if not skip_bi and not mrbayes and not dry_run:
        raise RuntimeError("MrBayes command `mb` was not found. Install MrBayes or use --skip-bi.")

    input_copy_dir = output_dir / "00_inputs"
    input_copy_dir.mkdir(parents=True, exist_ok=True)
    copied_fasta = input_copy_dir / fasta_path.name
    copied_fasta.write_text(fasta_path.read_text(encoding="utf-8"), encoding="utf-8")
    partition_file = output_dir / "00_inputs" / "iqtree_partitions.txt"
    write_iqtree_partition_file(partition_file, partitions)
    mrbayes_nexus = output_dir / "00_inputs" / "mrbayes_analysis.nex"
    write_mrbayes_nexus(mrbayes_nexus, copied_fasta, partitions, ngen=ngen, nruns=nruns, nchains=nchains, burninfrac=burninfrac)

    completed = 0
    if not dry_run and not skip_ml:
        if progress:
            progress.start("ml", f"ModelFinder + {bootstrap} bootstrap replicates")
        ml_dir = output_dir / "02_iqtree_ml"
        ml_dir.mkdir(parents=True, exist_ok=True)
        command = [
            iqtree or "iqtree",
            "-s",
            str(copied_fasta),
            "-p",
            str(partition_file),
            "-m",
            "MFP+MERGE",
            "-b",
            str(bootstrap),
            "-T",
            threads,
            "-pre",
            str(ml_dir / "cladistica_ml"),
        ]
        run_checked(command, output_dir, output_dir / "run_iqtree_ml.log")
        write_model_selection_summary(ml_dir / "cladistica_ml.iqtree", output_dir)
        completed += 1
        if progress:
            progress.succeed("ml", "ML.tre created")
    elif progress:
        progress.skip("ml", "Skipped by option")
    if not dry_run and not skip_bi:
        if progress:
            progress.start("bi", f"2 runs x {ngen} generations")
        bi_dir = output_dir / "03_mrbayes"
        bi_dir.mkdir(parents=True, exist_ok=True)
        analysis_copy = bi_dir / mrbayes_nexus.name
        analysis_copy.write_text(mrbayes_nexus.read_text(encoding="utf-8"), encoding="utf-8")
        run_checked([mrbayes or "mb", analysis_copy.name], bi_dir, output_dir / "run_mrbayes.log")
        completed += 1
        if progress:
            progress.succeed("bi", "BI.tre, run1.p, and run2.p created")
    elif progress:
        progress.skip("bi", "Skipped by option")
    copy_tree_outputs(output_dir, output_dir)
    summary = [
        "Cladistica iqtree_mrbayes_runner",
        "",
        f"Input FASTA: {fasta_path}",
        f"Samples: {stats['samples']}",
        f"Alignment length: {stats['alignment_length']}",
        f"Variable sites: {stats['variable_sites']}",
        f"Parsimony-informative sites: {stats['parsimony_informative_sites']}",
        f"IQ-TREE command: {iqtree or ''}",
        f"MrBayes command: {mrbayes or ''}",
        f"Dry run: {dry_run}",
        "",
        "MrBayes convergence screening:",
        "The generated run uses MrBayes diagnostics, but publication-level BI still requires manual trace/ESS inspection in Tracer or an equivalent tool.",
    ]
    (output_dir / "05_reports" / "analysis_summary.txt").write_text("\n".join(summary) + "\n", encoding="utf-8")
    return PipelineResult(str(output_dir), records_written=completed)
