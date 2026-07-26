from __future__ import annotations

from pathlib import Path

from .config import marker_map, safe_token
from .extraction import extract_marker_sequence, sequence_has_ambiguous_bases, split_accessions
from .io import read_table, write_fasta
from .models import PipelineResult
from .ncbi import fetch_one_genbank
from .progress import PipelineProgress


def sample_header(row: dict[str, str], marker: str, accession: str) -> str:
    sample = row.get("sample_id") or f"{safe_token(row.get('taxon', 'taxon'))}_{safe_token(row.get('voucher_or_isolate', accession))}"
    metadata = [
        f"taxon={row.get('taxon', '')}",
        f"marker={marker}",
        f"accession={accession}",
    ]
    voucher = row.get("voucher_or_isolate", "")
    if voucher:
        metadata.append(f"voucher={voucher}")
    return f"{sample}|" + "|".join(item.replace("\t", " ") for item in metadata)


def validate_selected_accessions(rows: list[dict[str, str]], marker_names: list[str]) -> None:
    if not rows:
        raise ValueError("The selected accession table has no data rows.")
    seen_samples: set[str] = set()
    for row_number, row in enumerate(rows, start=2):
        taxon = row.get("taxon", "").strip()
        sample_id = row.get("sample_id", "").strip()
        if not taxon and not sample_id:
            raise ValueError(f"Row {row_number} needs taxon or sample_id.")
        identity = sample_id or taxon
        if identity in seen_samples:
            raise ValueError(f"Duplicate sample_id/taxon in accession_selected.csv: {identity}")
        seen_samples.add(identity)
        accession_count = 0
        for marker in marker_names:
            accessions = split_accessions(row.get(marker, ""))
            if len(accessions) > 1:
                raise ValueError(
                    f"Row {row_number}, marker {marker} contains multiple accessions. "
                    "Keep one accession per marker for each selected sample."
                )
            accession_count += len(accessions)
        if accession_count == 0:
            raise ValueError(f"Row {row_number} has no accession in the requested marker columns.")


def download_fasta_from_accession_table(
    *,
    accession_table: Path,
    output_dir: Path,
    email: str,
    api_key: str = "",
    markers: list[str] | None = None,
    skip_existing: bool = True,
    progress: PipelineProgress | None = None,
) -> PipelineResult:
    if progress:
        progress.start("download", "Validating accession_selected.csv")
    output_dir.mkdir(parents=True, exist_ok=True)
    marker_defs = marker_map(markers)
    rows, fieldnames = read_table(accession_table)
    marker_names = [name for name in marker_defs if name in fieldnames]
    if not marker_names:
        raise ValueError("No requested marker columns were found in the accession table.")
    validate_selected_accessions(rows, marker_names)
    total_downloads = sum(
        len(split_accessions(row.get(marker_name, "")))
        for row in rows
        for marker_name in marker_names
    )
    completed_downloads = 0

    fasta_by_marker: dict[str, dict[str, str]] = {marker: {} for marker in marker_names}
    log_rows: list[dict[str, object]] = []
    warnings: list[str] = []

    for row in rows:
        for marker_name in marker_names:
            for accession in split_accessions(row.get(marker_name, "")):
                if progress:
                    progress.update("download", f"{marker_name}: {accession}")
                existing_records = fasta_by_marker[marker_name]
                header = sample_header(row, marker_name, accession)
                if skip_existing and header in existing_records:
                    completed_downloads += 1
                    if progress:
                        progress.set_progress(
                            "download",
                            completed_downloads * 100 / max(1, total_downloads),
                            f"{completed_downloads}/{total_downloads} sequences",
                        )
                    continue
                try:
                    record = fetch_one_genbank(accession, email=email, api_key=api_key)
                    extracted = extract_marker_sequence(record, marker_defs[marker_name])
                    existing_records[header] = extracted.sequence
                    log_rows.append(
                        {
                            "taxon": row.get("taxon", ""),
                            "sample_id": row.get("sample_id", ""),
                            "marker": marker_name,
                            "accession": accession,
                            "status": "downloaded",
                            "method": extracted.method,
                            "note": extracted.note,
                        }
                    )
                    log_rows[-1].update(
                        {
                            "length": len(extracted.sequence),
                            "has_ambiguous_bases": sequence_has_ambiguous_bases(extracted.sequence),
                            "marker_kind": extracted.marker_kind,
                            "coding_qc_status": extracted.coding_qc_status,
                            "coding_qc_note": extracted.coding_qc_note,
                        }
                    )
                except Exception as exc:
                    message = str(exc)
                    warnings.append(f"{marker_name}:{accession}: {message}")
                    log_rows.append(
                        {
                            "taxon": row.get("taxon", ""),
                            "sample_id": row.get("sample_id", ""),
                            "marker": marker_name,
                            "accession": accession,
                            "status": "failed",
                            "method": "",
                            "note": message,
                        }
                    )
                finally:
                    completed_downloads += 1
                    if progress:
                        progress.set_progress(
                            "download",
                            completed_downloads * 100 / max(1, total_downloads),
                            f"{completed_downloads}/{total_downloads} sequences",
                        )
    for marker_name, records in fasta_by_marker.items():
        if records:
            write_fasta(output_dir / f"{marker_name}.fasta", records)

    log_fields = [
        "taxon",
        "sample_id",
        "marker",
        "accession",
        "status",
        "method",
        "length",
        "has_ambiguous_bases",
        "marker_kind",
        "coding_qc_status",
        "coding_qc_note",
        "note",
    ]
    log_lines = [
        "Cladistica FASTA download log",
        f"Input: {accession_table}",
        f"Markers: {','.join(marker_names)}",
        "",
        "\t".join(log_fields),
    ]
    for row in log_rows:
        log_lines.append("\t".join(str(row.get(field, "")).replace("\t", " ").replace("\n", " ") for field in log_fields))
    (output_dir / "run.log").write_text("\n".join(log_lines) + "\n", encoding="utf-8")
    summary = [
        "Cladistica fasta_downloader",
        "",
        f"Input rows: {len(rows)}",
        f"Markers requested: {', '.join(marker_names)}",
        f"Sequences written: {sum(len(records) for records in fasta_by_marker.values())}",
        f"Warnings: {len(warnings)}",
    ]
    (output_dir / "summly.txt").write_text("\n".join(summary) + "\n", encoding="utf-8")
    if progress:
        detail = f"{sum(len(records) for records in fasta_by_marker.values())} sequences"
        if warnings:
            detail += f", {len(warnings)} warnings"
        progress.succeed("download", detail)
    return PipelineResult(str(output_dir), records_written=sum(len(records) for records in fasta_by_marker.values()), warnings=warnings)
