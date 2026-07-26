from __future__ import annotations

from collections import OrderedDict
from pathlib import Path

from .config import marker_map
from .io import read_fasta, sample_id_from_header, write_fasta, write_tsv
from .models import PipelineResult
from .progress import PipelineProgress


def marker_input_path(input_dir: Path, marker: str) -> Path | None:
    candidates = [
        input_dir / f"{marker}.fasta",
        input_dir / f"{marker}.aligned.fasta",
        input_dir / f"{marker}.masked.fasta",
        input_dir / f"{marker}.fa",
        input_dir / f"{marker}.fas",
    ]
    return next((path for path in candidates if path.exists()), None)


def read_aligned_marker(path: Path) -> tuple[OrderedDict[str, str], int]:
    raw = read_fasta(path)
    records: OrderedDict[str, str] = OrderedDict()
    for header, sequence in raw.items():
        sample_id = sample_id_from_header(header)
        if sample_id not in records:
            records[sample_id] = sequence.upper()
    lengths = {len(sequence) for sequence in records.values()}
    if len(lengths) > 1:
        raise ValueError(f"Alignment has inconsistent lengths: {path}")
    return records, next(iter(lengths)) if lengths else 0


def concatenate_alignments(
    *,
    input_dir: Path,
    output_dir: Path,
    markers: list[str] | None = None,
    strict: bool = False,
    progress: PipelineProgress | None = None,
) -> PipelineResult:
    if progress:
        progress.start("concat", "Reading marker alignments")
    marker_names = list(marker_map(markers))
    marker_records: dict[str, OrderedDict[str, str]] = {}
    marker_lengths: dict[str, int] = {}
    input_rows: list[dict[str, object]] = []
    for marker_index, marker in enumerate(marker_names, start=1):
        path = marker_input_path(input_dir, marker)
        if not path:
            if strict:
                raise FileNotFoundError(f"Missing aligned FASTA for marker: {marker}")
            input_rows.append({"marker": marker, "path": "", "records": 0, "alignment_length": "", "status": "missing_skipped"})
            if progress:
                progress.set_progress(
                    "concat",
                    marker_index * 70 / len(marker_names),
                    f"{marker_index}/{len(marker_names)} markers read",
                )
            continue
        records, length = read_aligned_marker(path)
        marker_records[marker] = records
        marker_lengths[marker] = length
        if progress:
            progress.update("concat", f"{marker}: {len(records)} samples x {length} bp")
            progress.set_progress(
                "concat",
                marker_index * 70 / len(marker_names),
                f"{marker_index}/{len(marker_names)} markers read",
            )
        input_rows.append({"marker": marker, "path": str(path), "records": len(records), "alignment_length": length, "status": "loaded"})

    if not marker_records:
        raise ValueError("No aligned marker FASTA files were found.")

    samples = sorted({sample for records in marker_records.values() for sample in records})
    if progress:
        progress.set_progress("concat", 80, f"{len(samples)} samples", approximate=True)
    concatenated: OrderedDict[str, str] = OrderedDict()
    presence_rows: list[dict[str, object]] = []
    for sample in samples:
        chunks: list[str] = []
        row: dict[str, object] = {"sample_id": sample}
        for marker in marker_records:
            sequence = marker_records[marker].get(sample)
            if sequence is None:
                chunks.append("?" * marker_lengths[marker])
                row[marker] = "missing"
            else:
                chunks.append(sequence)
                row[marker] = "present"
        concatenated[sample] = "".join(chunks)
        presence_rows.append(row)

    partitions: list[dict[str, object]] = []
    cursor = 1
    for marker in marker_records:
        length = marker_lengths[marker]
        partitions.append({"marker": marker, "start": cursor, "end": cursor + length - 1, "length": length})
        cursor += length

    write_fasta(output_dir / "concatenated_cpDNA.fasta", concatenated)
    write_tsv(output_dir / "partition_coordinates.tsv", partitions, ["marker", "start", "end", "length"])
    write_tsv(output_dir / "input_report.tsv", input_rows, ["marker", "path", "records", "alignment_length", "status"])
    write_tsv(output_dir / "region_presence.tsv", presence_rows, ["sample_id", *list(marker_records)])
    sample_rows = [
        {
            "sample_id": sample,
            "alignment_length": len(sequence),
            "regions_present": sum(1 for marker in marker_records if marker_records[marker].get(sample) is not None),
        }
        for sample, sequence in concatenated.items()
    ]
    write_tsv(output_dir / "concatenated_sample_report.tsv", sample_rows, ["sample_id", "alignment_length", "regions_present"])
    nexus_lines = ["begin sets;"]
    for partition in partitions:
        nexus_lines.append(f"  charset {partition['marker']} = {partition['start']}-{partition['end']};")
    nexus_lines.append("end;")
    (output_dir / "partition_coordinates.nexus_block.txt").write_text("\n".join(nexus_lines) + "\n", encoding="utf-8")
    summary = [
        "Cladistica alignment_concatenator",
        "",
        f"Samples: {len(samples)}",
        f"Regions concatenated: {len(marker_records)}",
        f"Alignment length: {sum(marker_lengths.values())}",
    ]
    (output_dir / "analysis_report.txt").write_text("\n".join(summary) + "\n", encoding="utf-8")
    if progress:
        progress.succeed(
            "concat",
            f"{len(samples)} samples x {sum(marker_lengths.values())} bp",
        )
    return PipelineResult(str(output_dir), records_written=len(samples))
