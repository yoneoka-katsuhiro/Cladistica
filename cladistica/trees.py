from __future__ import annotations

import shutil
import subprocess
import sys
import time
from collections import Counter
from collections.abc import Callable
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


def iqtree_major_version(command: str) -> int:
    """Return the IQ-TREE major version, or 0 if it cannot be determined.

    IQ-TREE 2/3 use ``-T`` for threads while IQ-TREE 1 uses ``-nt``. Detecting
    the version lets the runner pick the correct flag and warn when only the
    legacy v1 binary is available.
    """
    try:
        result = subprocess.run(
            [command, "--version"],
            text=True,
            capture_output=True,
            timeout=20,
            check=False,
        )
    except (OSError, subprocess.SubprocessError):
        return 0
    match = re.search(r"version\s+(\d+)", f"{result.stdout}\n{result.stderr}", flags=re.I)
    return int(match.group(1)) if match else 0


def iqtree_cli_flags(command: str) -> tuple[str, str, int]:
    """Return partition flag, thread flag, and detected IQ-TREE major version."""
    major_version = iqtree_major_version(command)
    if major_version == 1:
        return "-spp", "-nt", major_version
    return "-p", "-T", major_version


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
    lines.extend(["  ;", "end;", "", "begin mrbayes;", "  set autoclose=yes nowarn=yes;"])
    for row in partitions:
        lines.append(f"  charset {partition_label(str(row['marker']))} = {row['start']}-{row['end']};")
    if len(partitions) > 1:
        names = ", ".join(partition_label(str(row["marker"])) for row in partitions)
        lines.append(f"  partition cpDNA = {len(partitions)}: {names};")
        lines.append("  set partition=cpDNA;")
        lines.append("  unlink statefreq=(all) revmat=(all) shape=(all);")
        lines.append("  prset applyto=(all) ratepr=variable;")
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


def stop_process(process: subprocess.Popen[str]) -> None:
    if process.poll() is not None:
        return
    process.terminate()
    try:
        process.wait(timeout=5)
    except subprocess.TimeoutExpired:
        process.kill()
        process.wait()


def run_checked(
    command: list[str],
    cwd: Path,
    log_path: Path,
    *,
    on_line: Callable[[str], None] | None = None,
) -> None:
    log_path.parent.mkdir(parents=True, exist_ok=True)
    with log_path.open("w", encoding="utf-8") as log:
        process = subprocess.Popen(
            command,
            cwd=cwd,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            bufsize=1,
        )
        try:
            if process.stdout is None:
                raise RuntimeError("External command output stream was not created.")
            for raw_line in process.stdout:
                log.write(raw_line)
                log.flush()
                if on_line:
                    on_line(raw_line.rstrip())
            returncode = process.wait()
        except BaseException:
            stop_process(process)
            raise
        finally:
            if process.stdout:
                process.stdout.close()
    if returncode != 0:
        raise RuntimeError(f"Command failed ({returncode}): {' '.join(command)}. See {log_path}")


class IQTreeProgressParser:
    def __init__(self, progress: PipelineProgress, bootstrap: int) -> None:
        self.progress = progress
        self.bootstrap = max(1, bootstrap)
        self.phase = "model"
        self.candidate_total = 0
        self.bootstrap_replicate = 0
        self.ml_started = 0.0

    def __call__(self, line: str) -> None:
        if "Selecting individual models for" in line:
            self.progress.set_progress("model", 8, "testing partition models", approximate=True)
        total_match = re.search(r"about\s+(\d+)\s+total partition schemes", line)
        if total_match:
            self.candidate_total = int(total_match.group(1))
        if self.phase == "model":
            model_match = re.match(r"\s*(\d+)\s+([A-Za-z0-9+_.-]+)\s+", line)
            if model_match:
                tested = int(model_match.group(1))
                denominator = max(self.candidate_total, tested + 10)
                percent = min(92, 10 + 82 * tested / denominator)
                self.progress.set_progress(
                    "model",
                    percent,
                    f"{tested} candidates",
                    approximate=True,
                )
        if "CPU time for ModelFinder:" in line:
            self.progress.succeed("model")
            self.phase = "bootstrap"
            self.progress.start("bootstrap")
        replicate_match = re.search(r"START BOOTSTRAP REPLICATE NUMBER\s+(\d+)", line)
        if replicate_match:
            if self.phase == "model":
                self.progress.succeed("model")
                self.progress.start("bootstrap")
            self.phase = "bootstrap"
            self.bootstrap_replicate = int(replicate_match.group(1))
            completed = max(0, self.bootstrap_replicate - 1)
            self.progress.set_progress(
                "bootstrap",
                completed * 100 / self.bootstrap,
                f"{completed}/{self.bootstrap} replicates",
            )
        if "TREE SEARCH COMPLETED" in line and self.phase == "bootstrap":
            completed = min(self.bootstrap_replicate, self.bootstrap)
            self.progress.set_progress(
                "bootstrap",
                completed * 100 / self.bootstrap,
                f"{completed}/{self.bootstrap} replicates",
            )
        if "Consensus tree written to" in line:
            self.progress.succeed("bootstrap", f"{self.bootstrap}/{self.bootstrap} replicates")
            self.phase = "ml"
            self.ml_started = time.monotonic()
            self.progress.start("ml")
        if self.phase != "ml":
            return
        if "INITIALIZING CANDIDATE TREE SET" in line:
            self.progress.set_progress("ml", 25, "candidate trees", approximate=True)
        iteration_match = re.search(r"Iteration\s+(\d+)\s+/", line)
        if iteration_match:
            iteration = int(iteration_match.group(1))
            remaining_match = re.search(
                r"\((\d+)h:(\d+)m:(\d+)s left\)",
                line,
            )
            if remaining_match and self.ml_started:
                hours, minutes, seconds = map(int, remaining_match.groups())
                remaining = hours * 3600 + minutes * 60 + seconds
                elapsed = max(0.1, time.monotonic() - self.ml_started)
                percent = 100 * elapsed / (elapsed + remaining) if remaining else 88
            else:
                percent = 25 + min(60, iteration / 4)
            self.progress.set_progress(
                "ml",
                min(88, percent),
                f"iteration {iteration}",
                approximate=True,
            )
        if "TREE SEARCH COMPLETED" in line:
            self.progress.set_progress("ml", 92, "best tree found", approximate=True)
        if "FINALIZING TREE SEARCH" in line:
            self.progress.set_progress("ml", 97, "finalizing", approximate=True)
        if "Analysis results written to:" in line:
            self.progress.succeed("ml")

    def finish(self) -> None:
        for key in ("model", "bootstrap", "ml"):
            stage = self.progress.stages.get(key)
            if stage and stage.state not in {"success", "skipped"}:
                self.progress.succeed(key)


class MrBayesProgressParser:
    def __init__(self, progress: PipelineProgress, ngen: int) -> None:
        self.progress = progress
        self.ngen = max(1, ngen)
        self.summary_started = False

    def __call__(self, line: str) -> None:
        generation_match = re.match(r"\s*(\d+)\s+--", line)
        if generation_match and not self.summary_started:
            generation = min(int(generation_match.group(1)), self.ngen)
            self.progress.set_progress(
                "bi",
                generation * 100 / self.ngen,
                f"{generation:,}/{self.ngen:,} generations",
            )
        if "Analysis completed in" in line and not self.summary_started:
            self.progress.succeed("bi", f"{self.ngen:,}/{self.ngen:,} generations")
            self.summary_started = True
            self.progress.start("bi_summary")
            self.progress.set_progress(
                "bi_summary",
                10,
                "reading traces",
                approximate=True,
            )
        summary_milestones = (
            ("Summarizing parameters in files", 35, "parameter summary"),
            ("Model parameter summaries over the runs", 55, "convergence statistics"),
            ("Summarizing trees in files", 70, "tree summary"),
            ("Summary statistics for informative taxon bipartitions", 85, "bipartitions"),
            ("Calculating tree probabilities", 95, "consensus tree"),
        )
        for marker, percent, detail in summary_milestones:
            if marker in line:
                if not self.summary_started:
                    self.summary_started = True
                    self.progress.start("bi_summary")
                self.progress.set_progress(
                    "bi_summary",
                    percent,
                    detail,
                    approximate=True,
                )

    def finish(self) -> None:
        bi_stage = self.progress.stages.get("bi")
        if bi_stage and bi_stage.state not in {"success", "skipped"}:
            self.progress.succeed("bi")
        if self.progress.stages.get("bi_summary"):
            self.progress.succeed("bi_summary")


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
            progress.start("model")
        ml_dir = output_dir / "02_iqtree_ml"
        ml_dir.mkdir(parents=True, exist_ok=True)
        partition_flag, threads_flag, major_version = iqtree_cli_flags(
            iqtree or "iqtree"
        )
        if major_version == 1:
            print(
                "WARNING: IQ-TREE version 1 was detected. Falling back to '-spp' "
                "for partitions and '-nt' for threads ('-p' and '-T' are IQ-TREE "
                "2+ syntax). IQ-TREE 2 or 3 is strongly recommended for "
                "partitioned '-m MFP+MERGE' ML analyses.",
                file=sys.stderr,
                flush=True,
            )
        command = [
            iqtree or "iqtree",
            "-s",
            str(copied_fasta),
            partition_flag,
            str(partition_file),
            "-m",
            "MFP+MERGE",
            "-b",
            str(bootstrap),
            threads_flag,
            threads,
            "-pre",
            str(ml_dir / "cladistica_ml"),
        ]
        iqtree_progress = IQTreeProgressParser(progress, bootstrap) if progress else None
        run_checked(
            command,
            output_dir,
            output_dir / "run_iqtree_ml.log",
            on_line=iqtree_progress,
        )
        if iqtree_progress:
            iqtree_progress.finish()
        write_model_selection_summary(ml_dir / "cladistica_ml.iqtree", output_dir)
        completed += 1
    elif progress:
        reason = "Dry run" if dry_run else "Skipped by option"
        for key in ("model", "bootstrap", "ml"):
            progress.skip(key, reason)
    if not dry_run and not skip_bi:
        if progress:
            progress.start("bi")
        bi_dir = output_dir / "03_mrbayes"
        bi_dir.mkdir(parents=True, exist_ok=True)
        analysis_copy = bi_dir / mrbayes_nexus.name
        analysis_copy.write_text(mrbayes_nexus.read_text(encoding="utf-8"), encoding="utf-8")
        mrbayes_progress = MrBayesProgressParser(progress, ngen) if progress else None
        run_checked(
            [mrbayes or "mb", analysis_copy.name],
            bi_dir,
            output_dir / "run_mrbayes.log",
            on_line=mrbayes_progress,
        )
        if mrbayes_progress:
            mrbayes_progress.finish()
        completed += 1
    elif progress:
        reason = "Dry run" if dry_run else "Skipped by option"
        for key in ("bi", "bi_summary"):
            progress.skip(key, reason)
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
