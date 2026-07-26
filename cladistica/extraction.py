from __future__ import annotations

import re
from dataclasses import dataclass

from .config import Marker, normalize_marker_text


@dataclass
class ExtractedSequence:
    sequence: str
    method: str
    note: str = ""
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
    coding_qc_status: str = "not_applicable"
    coding_qc_note: str = ""


def split_accessions(value: str) -> list[str]:
    accessions = [item.strip() for item in re.split(r"[;,|]\s*", value or "") if item.strip()]
    return [re.sub(r"\s+", "", item) for item in accessions]


def sequence_has_ambiguous_bases(sequence: str) -> bool:
    return any(base.upper() in set("RYSWKMBDHVN?") for base in sequence)


def qualifier_values(feature: object, key: str) -> list[str]:
    values = (getattr(feature, "qualifiers", {}) or {}).get(key, []) or []
    return [str(value).strip() for value in values if str(value).strip()]


def first_qualifier(feature: object, key: str) -> str:
    values = qualifier_values(feature, key)
    return values[0] if values else ""


def qualifier_int(feature: object, key: str) -> int | None:
    value = first_qualifier(feature, key)
    if not value:
        return None
    try:
        return int(value)
    except ValueError:
        return None


def feature_type(feature: object) -> str:
    return str(getattr(feature, "type", "") or "").lower()


def feature_is_partial(feature: object) -> bool:
    location_text = str(getattr(feature, "location", "") or "")
    if "<" in location_text or ">" in location_text:
        return True
    return bool(qualifier_values(feature, "partial"))


def feature_is_pseudo(feature: object) -> bool:
    qualifiers = getattr(feature, "qualifiers", {}) or {}
    return any(key.lower() in {"pseudo", "pseudogene"} for key in qualifiers)


def stop_codons_for_table(transl_table: str) -> set[str]:
    table = str(transl_table or "11")
    try:
        from Bio.Data import CodonTable  # type: ignore

        codon_table = CodonTable.unambiguous_dna_by_id[int(table)]
        return {codon.upper() for codon in codon_table.stop_codons}
    except Exception:
        return {"TAA", "TAG", "TGA"}


def coding_qc(
    sequence: str,
    *,
    marker: Marker,
    feature: object | None = None,
    codon_start: int | None = None,
    transl_table: str = "",
) -> dict[str, object]:
    if marker.kind != "cds":
        return {
            "marker_kind": marker.kind,
            "feature_type": feature_type(feature) if feature is not None else "",
            "codon_start": codon_start,
            "transl_table": transl_table,
            "is_partial": feature_is_partial(feature) if feature is not None else False,
            "is_pseudo": False,
            "non_triplet": False,
            "internal_stop": False,
            "ambiguous_codons": 0,
            "evaluable_codons": 0,
            "clean_codon_fraction": 0.0,
            "coding_qc_status": "not_applicable",
            "coding_qc_note": "noncoding marker",
        }

    codon_start_value = codon_start or (qualifier_int(feature, "codon_start") if feature is not None else None) or 1
    transl_table_value = transl_table or (first_qualifier(feature, "transl_table") if feature is not None else "") or "11"
    offset = max(codon_start_value - 1, 0)
    coding_seq = sequence.upper()[offset:]
    tail_nt = len(coding_seq) % 3
    stop_codons = stop_codons_for_table(transl_table_value)
    ambiguous_codons = 0
    evaluable_codons = 0
    clean_codons = 0
    stop_indices: list[int] = []
    evaluable_indices: list[int] = []
    for index in range(0, len(coding_seq) - 2, 3):
        codon = coding_seq[index : index + 3]
        if any(base in "-?." for base in codon):
            continue
        codon_index = index // 3
        evaluable_indices.append(codon_index)
        evaluable_codons += 1
        if any(base not in {"A", "C", "G", "T"} for base in codon):
            ambiguous_codons += 1
            continue
        if codon in stop_codons:
            stop_indices.append(codon_index)
            continue
        clean_codons += 1
    terminal_index = evaluable_indices[-1] if evaluable_indices else -1
    internal_stop = any(index != terminal_index for index in stop_indices) if terminal_index >= 0 else False
    is_pseudo = feature_is_pseudo(feature) if feature is not None else False
    is_partial = feature_is_partial(feature) if feature is not None else True
    length_warning = len(sequence) < marker.expected_min_length * 0.25
    notes: list[str] = []
    status = "ok"
    if is_pseudo:
        status = "fail"
        notes.append("pseudo_or_pseudogene_feature")
    if internal_stop:
        status = "fail"
        notes.append("internal_stop")
    if tail_nt:
        if status != "fail":
            status = "warning"
        notes.append(f"partial_codon_tail={tail_nt}")
    if length_warning:
        if status != "fail":
            status = "warning"
        notes.append("short_relative_to_expected_marker_length")
    if is_partial:
        notes.append("partial_cds_allowed")
    return {
        "marker_kind": marker.kind,
        "feature_type": feature_type(feature) if feature is not None else "",
        "codon_start": codon_start_value,
        "transl_table": transl_table_value,
        "is_partial": is_partial,
        "is_pseudo": is_pseudo,
        "non_triplet": bool(tail_nt),
        "internal_stop": internal_stop,
        "ambiguous_codons": ambiguous_codons,
        "evaluable_codons": evaluable_codons,
        "clean_codon_fraction": clean_codons / evaluable_codons if evaluable_codons else 0.0,
        "coding_qc_status": status,
        "coding_qc_note": ";".join(notes) if notes else "coding fragment accepted; start/terminal stop not required",
    }


def build_extracted_sequence(
    *,
    sequence: str,
    method: str,
    note: str,
    marker: Marker,
    feature: object | None = None,
) -> ExtractedSequence:
    qc = coding_qc(sequence, marker=marker, feature=feature)
    return ExtractedSequence(sequence=sequence, method=method, note=note, **qc)


def feature_text(feature: object) -> str:
    qualifiers = getattr(feature, "qualifiers", {}) or {}
    parts: list[str] = [str(getattr(feature, "type", "") or "")]
    for values in qualifiers.values():
        parts.extend(str(value) for value in values)
    return normalize_marker_text(" ".join(parts))


def feature_matches_marker(feature: object, marker: Marker) -> bool:
    aliases = [normalize_marker_text(marker.name), *[normalize_marker_text(alias) for alias in marker.aliases]]
    text = feature_text(feature)
    return any(alias and alias in text for alias in aliases)


def extract_best_feature(record: object, marker: Marker) -> ExtractedSequence | None:
    allowed_types = {"cds"} if marker.kind == "cds" else {"gene", "trna", "misc_feature"}
    fallback_types = {"gene", "misc_feature"} if marker.kind == "cds" else set()
    candidates: list[tuple[int, int, int, ExtractedSequence]] = []
    for feature in getattr(record, "features", []) or []:
        ftype = feature_type(feature)
        if ftype not in allowed_types | fallback_types:
            continue
        if not feature_matches_marker(feature, marker):
            continue
        seq = str(feature.extract(getattr(record, "seq"))).upper()
        if not seq:
            continue
        extracted = build_extracted_sequence(
            sequence=seq,
            method=f"annotated_{ftype}_feature",
            note=str(getattr(feature, "location", "")),
            marker=marker,
            feature=feature,
        )
        type_rank = 2 if ftype == "cds" else 1 if ftype in allowed_types else 0
        qc_rank = {"ok": 3, "not_applicable": 3, "warning": 2, "fail": 0}.get(extracted.coding_qc_status, 1)
        candidates.append((qc_rank, type_rank, len(seq), extracted))
    if not candidates:
        return None
    return max(candidates, key=lambda item: item[:3])[3]


def extract_intergenic_region(record: object, marker: Marker) -> ExtractedSequence | None:
    if marker.name == "trnL-F":
        left_aliases = ("trnl", "trnleuaa", "trnluaa")
        right_aliases = ("trnf", "trnfgaa")
    elif marker.name == "rps4-trnS":
        left_aliases = ("rps4",)
        right_aliases = ("trns", "trnsgga", "trnsuga")
    else:
        return None
    left_positions: list[int] = []
    right_positions: list[int] = []
    for feature in getattr(record, "features", []) or []:
        text = feature_text(feature)
        location = getattr(feature, "location", None)
        if location is None:
            continue
        try:
            start = int(location.start)
            end = int(location.end)
        except Exception:
            continue
        if any(alias in text for alias in left_aliases):
            left_positions.append(end)
        if any(alias in text for alias in right_aliases):
            right_positions.append(start)
    candidates: list[tuple[int, int]] = []
    for left in left_positions:
        for right in right_positions:
            if left < right:
                candidates.append((left, right))
    if not candidates:
        return None
    start, end = min(candidates, key=lambda item: item[1] - item[0])
    sequence = str(getattr(record, "seq")[start:end]).upper()
    if sequence:
        return build_extracted_sequence(sequence=sequence, method="flanking_features", note=f"{start + 1}..{end}", marker=marker)
    return None


def extract_marker_sequence(record: object, marker: Marker) -> ExtractedSequence:
    sequence = str(getattr(record, "seq", "") or "").upper()
    if not sequence:
        raise ValueError("empty GenBank sequence")
    feature = extract_best_feature(record, marker)
    if feature:
        return feature
    if marker.kind != "cds":
        intergenic = extract_intergenic_region(record, marker)
        if intergenic:
            return intergenic
    if len(sequence) <= marker.long_record_threshold:
        return build_extracted_sequence(sequence=sequence, method="whole_record_short", note=f"length={len(sequence)}", marker=marker)
    raise ValueError(f"could not extract {marker.name} from long GenBank record")
