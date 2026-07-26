from __future__ import annotations

import shutil
import subprocess
from collections import OrderedDict
from pathlib import Path

from .config import marker_map
from .io import read_fasta, sample_id_from_header, write_fasta, write_tsv
from .models import PipelineResult
from .progress import PipelineProgress


def detect_muscle(command: str = "") -> str | None:
    if command:
        return command if shutil.which(command) else None
    for name in ("muscle", "muscle5"):
        found = shutil.which(name)
        if found:
            return found
    return None


def marker_fasta_path(base_dir: Path, marker: str) -> Path:
    for suffix in (".fasta", ".fa", ".fas", ".fna"):
        path = base_dir / f"{marker}{suffix}"
        if path.exists():
            return path
    return base_dir / f"{marker}.fasta"


def collect_marker_records(marker: str, genbank_dir: Path, user_dir: Path) -> tuple[OrderedDict[str, str], list[dict[str, object]]]:
    combined: OrderedDict[str, str] = OrderedDict()
    duplicates: list[dict[str, object]] = []
    for source, base_dir in (("genbank", genbank_dir), ("user", user_dir)):
        path = marker_fasta_path(base_dir, marker)
        for header, sequence in read_fasta(path).items():
            sample_id = sample_id_from_header(header)
            existing_header = next((key for key in combined if sample_id_from_header(key) == sample_id), "")
            if existing_header:
                duplicates.append(
                    {
                        "marker": marker,
                        "sample_id": sample_id,
                        "kept_source": source if source == "user" else "previous",
                        "discarded_header": existing_header if source == "user" else header,
                    }
                )
                if source == "user":
                    del combined[existing_header]
                    combined[header] = sequence
                continue
            combined[header] = sequence
    return combined, duplicates


def run_muscle(input_fasta: Path, output_fasta: Path, command: str) -> None:
    output_fasta.parent.mkdir(parents=True, exist_ok=True)
    attempts = [
        [command, "-align", str(input_fasta), "-output", str(output_fasta)],
        [command, "-in", str(input_fasta), "-out", str(output_fasta)],
    ]
    errors: list[str] = []
    for attempt in attempts:
        result = subprocess.run(attempt, text=True, capture_output=True, check=False)
        if result.returncode == 0 and output_fasta.exists():
            return
        errors.append(result.stderr.strip() or result.stdout.strip() or f"exit {result.returncode}")
    raise RuntimeError("MUSCLE failed: " + " | ".join(error for error in errors if error))


def convert_terminal_gaps_to_missing(sequence: str) -> str:
    chars = list(sequence)
    first = next((index for index, char in enumerate(chars) if char not in {"-", "?"}), None)
    if first is None:
        return "?" * len(chars)
    last = len(chars) - 1 - next(index for index, char in enumerate(reversed(chars)) if char not in {"-", "?"})
    for index in range(0, first):
        if chars[index] == "-":
            chars[index] = "?"
    for index in range(last + 1, len(chars)):
        if chars[index] == "-":
            chars[index] = "?"
    return "".join(chars)


def terminal_process_alignment(input_fasta: Path, output_fasta: Path) -> dict[str, object]:
    records = read_fasta(input_fasta)
    processed = {header: convert_terminal_gaps_to_missing(sequence) for header, sequence in records.items()}
    write_fasta(output_fasta, processed)
    lengths = {len(sequence) for sequence in processed.values()}
    return {
        "records": len(processed),
        "alignment_length": next(iter(lengths)) if len(lengths) == 1 else "",
        "length_status": "ok" if len(lengths) <= 1 else "inconsistent",
    }


def align_by_marker(
    *,
    genbank_dir: Path,
    user_dir: Path,
    output_dir: Path,
    markers: list[str] | None = None,
    muscle_command: str = "",
    dry_run: bool = False,
    progress: PipelineProgress | None = None,
) -> PipelineResult:
    if progress:
        progress.start("align", "Checking MUSCLE and marker FASTA")
    marker_names = list(marker_map(markers))
    muscle = detect_muscle(muscle_command)
    if not muscle and not dry_run:
        raise RuntimeError("MUSCLE was not found. Install MUSCLE or pass --dry-run for validation only.")

    duplicate_rows: list[dict[str, object]] = []
    rejected_rows: list[dict[str, object]] = []
    qc_rows: list[dict[str, object]] = []
    written = 0
    for marker_index, marker in enumerate(marker_names, start=1):
        records, duplicates = collect_marker_records(marker, genbank_dir, user_dir)
        if progress and records:
            progress.update("align", f"{marker}: {len(records)} sequences")
        duplicate_rows.extend(duplicates)
        combined_path = output_dir / "intermediate" / "combined_by_marker" / f"{marker}.fasta"
        raw_path = output_dir / "intermediate" / "raw_muscle_alignment_by_marker" / f"{marker}.fasta"
        aligned_path = output_dir / "aligned_by_marker" / f"{marker}.fasta"
        if not records:
            rejected_rows.append({"marker": marker, "sample_id": "", "reason": "no_input_sequences"})
            qc_rows.append({"marker": marker, "records": 0, "alignment_length": "", "status": "skipped_no_input"})
            if progress:
                progress.set_progress(
                    "align",
                    marker_index * 100 / len(marker_names),
                    f"{marker_index}/{len(marker_names)} markers",
                )
            continue
        valid = OrderedDict((header, sequence.upper().replace(" ", "")) for header, sequence in records.items() if sequence.strip())
        if len(valid) != len(records):
            rejected_rows.append({"marker": marker, "sample_id": "", "reason": "empty_sequence"})
        write_fasta(combined_path, valid)
        if dry_run:
            qc_rows.append({"marker": marker, "records": len(valid), "alignment_length": "", "status": "validated"})
            if progress:
                progress.set_progress(
                    "align",
                    marker_index * 100 / len(marker_names),
                    f"{marker_index}/{len(marker_names)} markers",
                )
            continue
        if len(valid) == 1:
            write_fasta(raw_path, valid)
        else:
            run_muscle(combined_path, raw_path, muscle or "muscle")
        stats = terminal_process_alignment(raw_path, aligned_path)
        qc_rows.append({"marker": marker, **stats, "status": stats["length_status"]})
        written += 1
        if progress:
            progress.set_progress(
                "align",
                marker_index * 100 / len(marker_names),
                f"{marker_index}/{len(marker_names)} markers",
            )

    write_tsv(output_dir / "qc" / "duplicate_records.tsv", duplicate_rows, ["marker", "sample_id", "kept_source", "discarded_header"])
    write_tsv(output_dir / "qc" / "rejected_before_alignment.tsv", rejected_rows, ["marker", "sample_id", "reason"])
    write_tsv(output_dir / "qc" / "all_markers_qc.tsv", qc_rows, ["marker", "records", "alignment_length", "status"])
    write_tsv(
        output_dir / "parameters.tsv",
        [{"parameter": "markers", "value": ",".join(marker_names)}, {"parameter": "muscle", "value": muscle or ""}, {"parameter": "dry_run", "value": dry_run}],
        ["parameter", "value"],
    )
    summary = [
        "Cladistica fasta_aligner_by_marker",
        "",
        f"Markers requested: {', '.join(marker_names)}",
        f"Marker alignments written: {written}",
        f"Dry run: {dry_run}",
    ]
    (output_dir / "alignment_report.txt").write_text("\n".join(summary) + "\n", encoding="utf-8")
    if progress:
        progress.succeed("align", f"{written} marker alignments")
    return PipelineResult(str(output_dir), records_written=written)
