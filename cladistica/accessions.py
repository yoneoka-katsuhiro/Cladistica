from __future__ import annotations

import re
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable

from .config import UNCERTAIN_NAME_PATTERNS, Marker, marker_map, normalize_marker_text, safe_token
from .extraction import extract_marker_sequence
from .io import write_csv, write_fasta
from .models import AccessionCandidate, MarkerHit, PipelineResult


AMBIGUOUS_DNA = set("RYSWKMBDHVN?")
PRIMARY_ACCESSION_FIELDS = ("accessions", "accession", "id", "name")


@dataclass
class RecordLike:
    accession: str
    organism: str
    description: str
    sequence: str
    voucher_or_isolate: str
    publication_status: str
    publication_score: int
    pubmed_ids: set[str]
    reference_titles: set[str]
    query_name: str
    marker_names: set[str]


@dataclass
class ExtractedMarkerRecord:
    accession: str
    organism: str
    marker: str
    sequence: str
    extraction_method: str
    extraction_note: str = ""
    marker_kind: str = ""
    feature_type: str = ""
    codon_start: int | None = None
    transl_table: str = ""
    is_partial: bool = False
    is_pseudo: bool = False
    non_triplet: bool = False
    internal_stop: bool = False
    ambiguous_codons: int = 0
    evaluable_codons: int = 0
    clean_codon_fraction: float = 0.0
    coding_qc_status: str = ""
    coding_qc_note: str = ""
    voucher_or_isolate: str = ""
    publication_status: str = "unknown"
    publication_score: int = 0
    pubmed_ids: set[str] = field(default_factory=set)
    reference_titles: set[str] = field(default_factory=set)
    query_name: str = ""


def is_uncertain_taxon(name: str) -> bool:
    text = f" {name.strip().lower()} "
    return any(re.search(pattern, text, flags=re.IGNORECASE) for pattern in UNCERTAIN_NAME_PATTERNS)


def normalize_taxon_name(name: str) -> str:
    return re.sub(r"\s+", " ", name.strip())


def taxon_from_organism(organism: str) -> str:
    parts = normalize_taxon_name(organism).split()
    if len(parts) < 2:
        return normalize_taxon_name(organism)
    return " ".join(parts[:2])


def normalize_sample_id(taxon: str, voucher_or_isolate: str, accession: str) -> str:
    if voucher_or_isolate:
        return f"{safe_token(taxon)}_{safe_token(voucher_or_isolate)}"
    return f"{safe_token(taxon)}_{safe_token(accession)}"


def ambiguous_count(sequence: str) -> int:
    return sum(1 for base in sequence.upper() if base in AMBIGUOUS_DNA)


def publication_score_from_references(references: Iterable[object]) -> tuple[str, int, set[str], set[str]]:
    pubmed_ids: set[str] = set()
    titles: set[str] = set()
    status = "unknown"
    score = 0
    for reference in references:
        title = str(getattr(reference, "title", "") or "").strip()
        journal = str(getattr(reference, "journal", "") or "").strip()
        pubmed_id = str(getattr(reference, "pubmed_id", "") or "").strip()
        if title:
            titles.add(title)
        if pubmed_id:
            pubmed_ids.add(pubmed_id)
            status = "pubmed_indexed"
            score = max(score, 3)
            continue
        journal_lower = journal.lower()
        if journal and not any(token in journal_lower for token in ("unpublished", "submitted")):
            status = "published_or_accepted"
            score = max(score, 2)
        elif "submitted" in journal_lower:
            status = "submitted_to_genbank"
            score = max(score, 1)
        elif "unpublished" in journal_lower and score == 0:
            status = "unpublished"
    return status, score, pubmed_ids, titles


def voucher_from_record(record: object) -> str:
    qualifiers: list[str] = []
    for feature in getattr(record, "features", []) or []:
        data = getattr(feature, "qualifiers", {}) or {}
        for key in ("specimen_voucher", "isolate", "clone", "cultivar", "strain"):
            for value in data.get(key, []) or []:
                if value and str(value).strip():
                    qualifiers.append(str(value).strip())
    if qualifiers:
        return qualifiers[0]
    annotations = getattr(record, "annotations", {}) or {}
    for key in ("isolate", "specimen_voucher", "clone"):
        value = str(annotations.get(key, "") or "").strip()
        if value:
            return value
    return ""


def accession_from_record(record: object) -> str:
    annotations = getattr(record, "annotations", {}) or {}
    accessions = annotations.get("accessions")
    if isinstance(accessions, list) and accessions:
        return str(accessions[0])
    for field in PRIMARY_ACCESSION_FIELDS:
        value = getattr(record, field, "")
        if value:
            return str(value).split()[0]
    return ""


def marker_names_from_record(record: object, markers: dict[str, Marker]) -> set[str]:
    text_parts = [
        str(getattr(record, "name", "") or ""),
        str(getattr(record, "id", "") or ""),
        str(getattr(record, "description", "") or ""),
    ]
    for feature in getattr(record, "features", []) or []:
        qualifiers = getattr(feature, "qualifiers", {}) or {}
        for values in qualifiers.values():
            text_parts.extend(str(value) for value in values)
    text = normalize_marker_text(" ".join(text_parts))
    found: set[str] = set()
    for marker in markers.values():
        aliases = (marker.name, *marker.aliases)
        if any(normalize_marker_text(alias) in text for alias in aliases):
            found.add(marker.name)
    return found


def record_to_recordlike(record: object, query_name: str, markers: dict[str, Marker]) -> RecordLike:
    annotations = getattr(record, "annotations", {}) or {}
    references = annotations.get("references", []) or []
    status, score, pubmed_ids, titles = publication_score_from_references(references)
    return RecordLike(
        accession=accession_from_record(record),
        organism=str(annotations.get("organism", "") or getattr(record, "description", "")).strip(),
        description=str(getattr(record, "description", "") or ""),
        sequence=str(getattr(record, "seq", "") or ""),
        voucher_or_isolate=voucher_from_record(record),
        publication_status=status,
        publication_score=score,
        pubmed_ids=pubmed_ids,
        reference_titles=titles,
        query_name=query_name,
        marker_names=marker_names_from_record(record, markers),
    )


def record_to_extracted_marker_records(
    record: object,
    query_name: str,
    markers: dict[str, Marker],
) -> tuple[list[ExtractedMarkerRecord], list[dict[str, object]]]:
    annotations = getattr(record, "annotations", {}) or {}
    organism = str(annotations.get("organism", "") or getattr(record, "description", "")).strip()
    accession = accession_from_record(record)
    taxon = taxon_from_organism(organism)
    if not taxon or is_uncertain_taxon(taxon) or is_uncertain_taxon(organism):
        return [], [
            {
                "taxon": organism or taxon,
                "accession": accession,
                "reason": "uncertain_or_incomplete_taxon_name",
                "query_name": query_name,
            }
        ]

    detected_markers = marker_names_from_record(record, markers)
    if not detected_markers:
        return [], [
            {
                "taxon": taxon,
                "accession": accession,
                "reason": "no_requested_marker_detected",
                "query_name": query_name,
            }
        ]

    references = annotations.get("references", []) or []
    status, score, pubmed_ids, titles = publication_score_from_references(references)
    voucher = voucher_from_record(record)
    extracted_records: list[ExtractedMarkerRecord] = []
    rejected: list[dict[str, object]] = []
    for marker_name in sorted(detected_markers):
        if marker_name not in markers:
            continue
        try:
            extracted = extract_marker_sequence(record, markers[marker_name])
            extracted_records.append(
                ExtractedMarkerRecord(
                    accession=accession,
                    organism=organism,
                    marker=marker_name,
                    sequence=extracted.sequence,
                    extraction_method=extracted.method,
                    extraction_note=extracted.note,
                    marker_kind=extracted.marker_kind,
                    feature_type=extracted.feature_type,
                    codon_start=extracted.codon_start,
                    transl_table=extracted.transl_table,
                    is_partial=extracted.is_partial,
                    is_pseudo=extracted.is_pseudo,
                    non_triplet=extracted.non_triplet,
                    internal_stop=extracted.internal_stop,
                    ambiguous_codons=extracted.ambiguous_codons,
                    evaluable_codons=extracted.evaluable_codons,
                    clean_codon_fraction=extracted.clean_codon_fraction,
                    coding_qc_status=extracted.coding_qc_status,
                    coding_qc_note=extracted.coding_qc_note,
                    voucher_or_isolate=voucher,
                    publication_status=status,
                    publication_score=score,
                    pubmed_ids=set(pubmed_ids),
                    reference_titles=set(titles),
                    query_name=query_name,
                )
            )
        except Exception as exc:
            rejected.append(
                {
                    "taxon": taxon,
                    "accession": accession,
                    "reason": f"{marker_name}_extraction_failed: {exc}",
                    "query_name": query_name,
                }
            )
    if not extracted_records and not rejected:
        rejected.append(
            {
                "taxon": taxon,
                "accession": accession,
                "reason": "no_requested_marker_extracted",
                "query_name": query_name,
            }
        )
    return extracted_records, rejected


def candidates_from_records(records: Iterable[RecordLike], markers: dict[str, Marker]) -> tuple[list[AccessionCandidate], list[dict[str, object]]]:
    extracted_records: list[ExtractedMarkerRecord] = []
    rejected: list[dict[str, object]] = []
    for record in records:
        taxon = taxon_from_organism(record.organism)
        if not taxon or is_uncertain_taxon(taxon) or is_uncertain_taxon(record.organism):
            rejected.append(
                {
                    "taxon": record.organism or taxon,
                    "accession": record.accession,
                    "reason": "uncertain_or_incomplete_taxon_name",
                    "query_name": record.query_name,
                }
            )
            continue
        if not record.marker_names:
            rejected.append(
                {
                    "taxon": taxon,
                    "accession": record.accession,
                    "reason": "no_requested_marker_detected",
                    "query_name": record.query_name,
                }
            )
            continue
        for marker in sorted(record.marker_names):
            if marker not in markers:
                continue
            extracted_records.append(
                ExtractedMarkerRecord(
                    accession=record.accession,
                    organism=record.organism,
                    marker=marker,
                    sequence=record.sequence,
                    extraction_method="legacy_record_sequence",
                    extraction_note="selection used the RecordLike sequence supplied by caller",
                    voucher_or_isolate=record.voucher_or_isolate,
                    publication_status=record.publication_status,
                    publication_score=record.publication_score,
                    pubmed_ids=set(record.pubmed_ids),
                    reference_titles=set(record.reference_titles),
                    query_name=record.query_name,
                )
            )
    candidates, extracted_rejected = candidates_from_extracted_records(extracted_records, markers)
    return candidates, [*rejected, *extracted_rejected]


def candidates_from_extracted_records(
    records: Iterable[ExtractedMarkerRecord],
    markers: dict[str, Marker],
) -> tuple[list[AccessionCandidate], list[dict[str, object]]]:
    grouped: dict[tuple[str, str], AccessionCandidate] = {}
    rejected: list[dict[str, object]] = []
    for record in records:
        taxon = taxon_from_organism(record.organism)
        if not taxon or is_uncertain_taxon(taxon) or is_uncertain_taxon(record.organism):
            rejected.append(
                {
                    "taxon": record.organism or taxon,
                    "accession": record.accession,
                    "reason": "uncertain_or_incomplete_taxon_name",
                    "query_name": record.query_name,
                }
            )
            continue
        if record.marker not in markers:
            rejected.append(
                {
                    "taxon": taxon,
                    "accession": record.accession,
                    "reason": f"marker_not_requested: {record.marker}",
                    "query_name": record.query_name,
                }
            )
            continue
        if not record.sequence:
            rejected.append(
                {
                    "taxon": taxon,
                    "accession": record.accession,
                    "reason": f"{record.marker}_empty_sequence",
                    "query_name": record.query_name,
                }
            )
            continue
        marker_def = markers[record.marker]
        sample_id = normalize_sample_id(taxon, record.voucher_or_isolate, record.accession)
        key = (taxon, sample_id)
        candidate = grouped.get(key)
        if candidate is None:
            candidate = AccessionCandidate(
                taxon=taxon,
                sample_id=sample_id,
                voucher_or_isolate=record.voucher_or_isolate,
                publication_status=record.publication_status,
                publication_score=record.publication_score,
                query_name=record.query_name,
            )
            grouped[key] = candidate
        candidate.publication_score = max(candidate.publication_score, record.publication_score)
        if record.publication_score >= candidate.publication_score:
            candidate.publication_status = record.publication_status
        candidate.pubmed_ids.update(record.pubmed_ids)
        candidate.reference_titles.update(record.reference_titles)
        candidate.source_record_ids.add(record.accession)

        hit = MarkerHit(
            marker=record.marker,
            accession=record.accession,
            length=len(record.sequence),
            ambiguous_bases=ambiguous_count(record.sequence),
            marker_kind=record.marker_kind or marker_def.kind,
            source_record_id=record.accession,
            feature_note=record.extraction_note,
            extraction_method=record.extraction_method,
            feature_type=record.feature_type,
            codon_start=record.codon_start,
            transl_table=record.transl_table,
            is_partial=record.is_partial,
            is_pseudo=record.is_pseudo,
            non_triplet=record.non_triplet,
            internal_stop=record.internal_stop,
            ambiguous_codons=record.ambiguous_codons,
            evaluable_codons=record.evaluable_codons,
            clean_codon_fraction=record.clean_codon_fraction,
            coding_qc_status=record.coding_qc_status,
            coding_qc_note=record.coding_qc_note,
            sequence=record.sequence,
        )
        previous = candidate.marker_hits.get(record.marker)
        if previous is None or marker_hit_sort_key(hit) > marker_hit_sort_key(previous):
            candidate.marker_hits[record.marker] = hit

    return list(grouped.values()), rejected


def coding_qc_rank(status: str) -> int:
    return {"ok": 3, "not_applicable": 3, "warning": 2, "": 2, "fail": 0}.get(status, 1)


def marker_hit_sort_key(hit: MarkerHit) -> tuple[int, int, int, str]:
    return (1 if hit.is_clean else 0, coding_qc_rank(hit.coding_qc_status), hit.length, hit.accession)


def select_representatives(
    candidates: list[AccessionCandidate],
    markers: list[str],
    max_samples_per_taxon: int = 1,
    trusted_voucher_keywords: list[str] | None = None,
    include_trusted_extra: bool = False,
) -> list[AccessionCandidate]:
    if max_samples_per_taxon <= 0:
        raise ValueError("max_samples_per_taxon must be greater than 0")
    by_taxon: dict[str, list[AccessionCandidate]] = defaultdict(list)
    for candidate in candidates:
        by_taxon[candidate.taxon].append(candidate)
    selected: list[AccessionCandidate] = []
    trusted = [keyword.lower() for keyword in trusted_voucher_keywords or [] if keyword.strip()]
    for taxon in sorted(by_taxon):
        ranked = sorted(by_taxon[taxon], key=lambda item: item.sort_key(markers), reverse=True)
        keep = ranked[:max_samples_per_taxon]
        if include_trusted_extra and trusted:
            for candidate in ranked[max_samples_per_taxon:]:
                text = candidate.voucher_or_isolate.lower()
                if any(keyword in text for keyword in trusted):
                    keep.append(candidate)
        selected.extend(keep)
    return selected


def candidate_table_rows(candidates: list[AccessionCandidate], markers: list[str]) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for candidate in sorted(candidates, key=lambda item: (item.taxon, item.sample_id)):
        row: dict[str, object] = {
            "taxon": candidate.taxon,
            "sample_id": candidate.sample_id,
            "voucher_or_isolate": candidate.voucher_or_isolate,
            "marker_coverage": candidate.marker_count(markers),
            "clean_marker_count": candidate.clean_marker_count(markers),
            "ambiguous_bases": candidate.ambiguous_bases(markers),
            "total_marker_length": candidate.total_length(markers),
            "publication_status": candidate.publication_status,
            "publication_score": candidate.publication_score,
            "pubmed_ids": ";".join(sorted(candidate.pubmed_ids)),
            "source_record_ids": ";".join(sorted(candidate.source_record_ids)),
            "marker_lengths": ";".join(
                f"{marker}:{candidate.marker_hits[marker].length}"
                for marker in markers
                if marker in candidate.marker_hits
            ),
            "marker_ambiguous_bases": ";".join(
                f"{marker}:{candidate.marker_hits[marker].ambiguous_bases}"
                for marker in markers
                if marker in candidate.marker_hits
            ),
            "marker_qc": ";".join(
                f"{marker}:{candidate.marker_hits[marker].coding_qc_status}"
                for marker in markers
                if marker in candidate.marker_hits
            ),
            "marker_qc_notes": ";".join(
                f"{marker}:{candidate.marker_hits[marker].coding_qc_note}"
                for marker in markers
                if marker in candidate.marker_hits and candidate.marker_hits[marker].coding_qc_note
            ),
        }
        row.update(candidate.selected_accessions(markers))
        rows.append(row)
    return rows


def fasta_header(candidate: AccessionCandidate, marker: str, hit: MarkerHit) -> str:
    metadata = [
        f"taxon={candidate.taxon}",
        f"marker={marker}",
        f"accession={hit.accession}",
    ]
    if candidate.voucher_or_isolate:
        metadata.append(f"voucher={candidate.voucher_or_isolate}")
    if hit.extraction_method:
        metadata.append(f"extraction={hit.extraction_method}")
    return f"{candidate.sample_id}|" + "|".join(item.replace("\t", " ") for item in metadata)


def write_marker_fastas(output_dir: Path, candidates: list[AccessionCandidate], markers: list[str], directory_name: str) -> None:
    for marker in markers:
        records: dict[str, str] = {}
        for candidate in sorted(candidates, key=lambda item: (item.taxon, item.sample_id)):
            hit = candidate.marker_hits.get(marker)
            if hit and hit.sequence:
                records[fasta_header(candidate, marker, hit)] = hit.sequence
        if records:
            write_fasta(output_dir / directory_name / f"{marker}.fasta", records)


def write_accession_outputs(
    *,
    output_dir: Path,
    selected: list[AccessionCandidate],
    all_candidates: list[AccessionCandidate],
    rejected: list[dict[str, object]],
    markers: list[str],
    parameters: dict[str, object],
) -> PipelineResult:
    selected_fieldnames = [
        "taxon",
        "sample_id",
        "voucher_or_isolate",
        *markers,
        "marker_coverage",
        "clean_marker_count",
        "ambiguous_bases",
        "total_marker_length",
        "publication_status",
        "publication_score",
        "pubmed_ids",
        "source_record_ids",
        "marker_lengths",
        "marker_ambiguous_bases",
        "marker_qc",
        "marker_qc_notes",
        "selection_status",
    ]
    selected_keys = {(candidate.taxon, candidate.sample_id) for candidate in selected}
    all_rows = candidate_table_rows(all_candidates, markers)
    for row in all_rows:
        key = (str(row["taxon"]), str(row["sample_id"]))
        row["selection_status"] = "selected" if key in selected_keys else "not_selected"
        row["rejection_reason"] = ""
        row["query_name"] = ""

    seen_rejected: set[tuple[str, str, str, str]] = set()
    for rejected_row in rejected:
        key = (
            str(rejected_row.get("taxon", "")),
            str(rejected_row.get("accession", "")),
            str(rejected_row.get("reason", "")),
            str(rejected_row.get("query_name", "")),
        )
        if key in seen_rejected:
            continue
        seen_rejected.add(key)
        all_rows.append(
            {
                "taxon": key[0],
                "source_record_ids": key[1],
                "selection_status": "rejected",
                "rejection_reason": key[2],
                "query_name": key[3],
            }
        )

    all_fieldnames = [*selected_fieldnames, "rejection_reason", "query_name"]
    selected_rows = candidate_table_rows(selected, markers)
    for row in selected_rows:
        row["selection_status"] = "selected"
    write_csv(output_dir / "accession_all.csv", all_rows, all_fieldnames)
    write_csv(output_dir / "accession_selected.csv", selected_rows, selected_fieldnames)
    write_marker_fastas(output_dir, selected, markers, "fasta_by_marker")

    summary = [
        "Cladistica accession selection",
        "",
        "Default selection policy:",
        "1. Retrieve GenBank records and extract requested marker sequences first.",
        "2. Reject uncertain or incomplete taxon names and failed marker extractions.",
        "3. Group extracted sequences by taxon and voucher/isolate-derived sample ID.",
        "4. Select one sample per taxon by extracted marker coverage, clean/QC-passing marker count, publication evidence, and extracted sequence length.",
        "",
        f"Selected samples: {len(selected)}",
        f"Candidate samples: {len(all_candidates)}",
        f"Rejected records: {len(rejected)}",
    ]
    (output_dir / "summly.txt").write_text("\n".join(summary) + "\n", encoding="utf-8")
    log_lines = [
        "Cladistica accession log",
        "",
        "Parameters:",
        *[f"{key}: {value}" for key, value in parameters.items()],
        "",
        "Rejected records:",
    ]
    log_lines.extend(
        f"{row.get('accession', '')}\t{row.get('taxon', '')}\t{row.get('reason', '')}\t{row.get('query_name', '')}"
        for row in rejected
    )
    (output_dir / "run.log").write_text("\n".join(log_lines) + "\n", encoding="utf-8")
    return PipelineResult(str(output_dir), records_written=len(selected))


def build_accession_table_from_records(
    *,
    records: Iterable[RecordLike],
    output_dir: Path,
    markers: list[str] | None = None,
    max_samples_per_taxon: int = 1,
    trusted_voucher_keywords: list[str] | None = None,
    include_trusted_extra: bool = False,
    parameters: dict[str, object] | None = None,
) -> PipelineResult:
    marker_defs = marker_map(markers)
    marker_names = list(marker_defs)
    candidates, rejected = candidates_from_records(records, marker_defs)
    selected = select_representatives(
        candidates,
        marker_names,
        max_samples_per_taxon=max_samples_per_taxon,
        trusted_voucher_keywords=trusted_voucher_keywords,
        include_trusted_extra=include_trusted_extra,
    )
    merged_parameters = {
        "markers": ",".join(marker_names),
        "max_samples_per_taxon": max_samples_per_taxon,
        "include_trusted_extra": include_trusted_extra,
    }
    merged_parameters.update(parameters or {})
    return write_accession_outputs(
        output_dir=output_dir,
        selected=selected,
        all_candidates=candidates,
        rejected=rejected,
        markers=marker_names,
        parameters=merged_parameters,
    )


def build_accession_table_from_genbank_records(
    *,
    records: Iterable[tuple[str, object]],
    output_dir: Path,
    markers: list[str] | None = None,
    max_samples_per_taxon: int = 1,
    trusted_voucher_keywords: list[str] | None = None,
    include_trusted_extra: bool = False,
    parameters: dict[str, object] | None = None,
) -> PipelineResult:
    marker_defs = marker_map(markers)
    marker_names = list(marker_defs)
    extracted_records: list[ExtractedMarkerRecord] = []
    rejected: list[dict[str, object]] = []
    for query_name, record in records:
        extracted, record_rejected = record_to_extracted_marker_records(record, query_name, marker_defs)
        extracted_records.extend(extracted)
        rejected.extend(record_rejected)

    candidates, candidate_rejected = candidates_from_extracted_records(extracted_records, marker_defs)
    rejected.extend(candidate_rejected)
    selected = select_representatives(
        candidates,
        marker_names,
        max_samples_per_taxon=max_samples_per_taxon,
        trusted_voucher_keywords=trusted_voucher_keywords,
        include_trusted_extra=include_trusted_extra,
    )
    merged_parameters = {
        "markers": ",".join(marker_names),
        "max_samples_per_taxon": max_samples_per_taxon,
        "include_trusted_extra": include_trusted_extra,
        "selection_input": "extracted_marker_sequences",
    }
    merged_parameters.update(parameters or {})
    return write_accession_outputs(
        output_dir=output_dir,
        selected=selected,
        all_candidates=candidates,
        rejected=rejected,
        markers=marker_names,
        parameters=merged_parameters,
    )
