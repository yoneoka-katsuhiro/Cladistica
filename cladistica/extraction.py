from __future__ import annotations

import re
from dataclasses import dataclass

from .config import KNOWN_MARKERS, Marker, normalize_marker_text, resolve_marker


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
    marker_kind = marker.kind
    if marker_kind == "auto":
        marker_kind = "cds" if feature is not None and feature_type(feature) == "cds" else "noncoding"
    if marker_kind != "cds":
        return {
            "marker_kind": marker_kind,
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
            "coding_qc_note": "not a CDS feature",
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
        "marker_kind": marker_kind,
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


def marker_aliases(marker: Marker) -> set[str]:
    aliases: set[str] = set()
    for value in (marker.name, *marker.aliases):
        normalized = normalize_marker_text(value)
        if normalized:
            aliases.add(normalized)
    return aliases


def exact_identifier_match(value: str, aliases: set[str]) -> bool:
    normalized = normalize_marker_text(value)
    if normalized in aliases:
        return True
    for alias in aliases:
        if re.fullmatch(r"trn[a-z]+", alias) and normalized.startswith(alias):
            anticodon = normalized[len(alias) :]
            if len(anticodon) == 3 and set(anticodon) <= set("acgtu"):
                return True
    return False


def marker_in_text(value: str, marker: Marker) -> bool:
    for alias in (marker.name, *marker.aliases):
        words = re.findall(r"[A-Za-z0-9]+", alias)
        if words:
            pattern = r"(?<![A-Za-z0-9])" + r"[^A-Za-z0-9]*".join(map(re.escape, words)) + r"(?![A-Za-z0-9])"
            if re.search(pattern, value, re.IGNORECASE):
                return True
    return False


def feature_matches_identifier(feature: object, identifier: str) -> bool:
    """Match an annotated gene/RNA identifier without partial gene-name hits."""
    return feature_matches_marker(feature, resolve_marker(identifier))


def feature_matches_marker(feature: object, marker: Marker) -> bool:
    aliases = marker_aliases(marker)
    qualifiers = getattr(feature, "qualifiers", {}) or {}
    for key in ("gene", "standard_name", "label", "locus_tag"):
        for value in qualifiers.get(key, []) or []:
            if exact_identifier_match(str(value), aliases):
                return True
    for key in ("product", "note", "function"):
        for value in qualifiers.get(key, []) or []:
            if marker_in_text(str(value), marker):
                return True
    return False


def feature_is_specific_to_marker(feature: object, marker: Marker) -> bool:
    """Reject broad features whose note merely lists the requested subregion."""
    aliases = marker_aliases(marker)
    qualifiers = getattr(feature, "qualifiers", {}) or {}
    if marker.intron_of and feature_type(feature) == "intron":
        return feature_matches_identifier(feature, marker.intron_of) and (
            marker.intron_index is None or qualifier_int(feature, "number") == marker.intron_index
        )
    for key in ("gene", "standard_name", "label", "locus_tag"):
        for value in qualifiers.get(key, []) or []:
            if exact_identifier_match(str(value), aliases):
                return True

    for value in qualifiers.get("product", []) or []:
        normalized = normalize_marker_text(str(value))
        if normalized in aliases:
            return True
    for key in ("note", "function"):
        for value in qualifiers.get(key, []) or []:
            normalized = normalize_marker_text(str(value))
            if normalized in aliases or normalized.removeprefix("contains") in aliases:
                return True
    return False


def record_matches_marker(record: object, marker: Marker) -> bool:
    header = (
        " ".join(
            (
                str(getattr(record, "name", "") or ""),
                str(getattr(record, "id", "") or ""),
                str(getattr(record, "description", "") or ""),
            )
        )
    )
    if marker_in_text(header, marker):
        return True

    features = getattr(record, "features", []) or []
    if any(feature_matches_marker(feature, marker) for feature in features):
        return True
    if marker.flanks:
        left, right = marker.flanks
        return any(feature_matches_identifier(feature, left) for feature in features) and any(
            feature_matches_identifier(feature, right) for feature in features
        )
    if marker.intron_of:
        return any(feature_matches_identifier(feature, marker.intron_of) for feature in features)
    return False


def record_has_composite_marker_title(record: object, requested: Marker) -> bool:
    header = (
        " ".join(
            (
                str(getattr(record, "name", "") or ""),
                str(getattr(record, "id", "") or ""),
                str(getattr(record, "description", "") or ""),
            )
        )
    )
    mentioned: set[str] = set()
    for marker in KNOWN_MARKERS:
        if marker_in_text(header, marker):
            mentioned.add(marker.name)
    return bool(mentioned - {requested.name})


def feature_type_groups(marker: Marker) -> tuple[set[str], set[str]]:
    if marker.kind == "cds":
        return {"cds"}, {"gene", "misc_feature"}
    if marker.kind == "rrna":
        return {"rrna"}, {"gene", "misc_rna", "misc_feature"}
    if marker.kind == "trna":
        return {"trna"}, {"gene", "misc_rna", "misc_feature"}
    if marker.kind == "noncoding":
        return {"intron", "misc_feature", "ncrna", "misc_rna"}, {"gene", "rrna", "trna"}
    return {"cds", "rrna", "trna", "ncrna", "misc_rna"}, {"gene", "misc_feature"}


def extract_best_feature(record: object, marker: Marker) -> ExtractedSequence | None:
    if marker.intron_of:
        return None
    allowed_types, fallback_types = feature_type_groups(marker)
    candidates: list[tuple[int, int, int, ExtractedSequence]] = []
    for feature in getattr(record, "features", []) or []:
        ftype = feature_type(feature)
        if ftype not in allowed_types | fallback_types:
            continue
        if not feature_is_specific_to_marker(feature, marker):
            continue
        seq = str(feature.extract(getattr(record, "seq"))).upper()
        if not seq:
            continue
        if len(seq) > marker.long_record_threshold:
            continue
        extracted = build_extracted_sequence(
            sequence=seq,
            method=f"annotated_{ftype}_feature",
            note=str(getattr(feature, "location", "")),
            marker=marker,
            feature=feature,
        )
        type_rank = 2 if ftype in allowed_types else 1
        qc_rank = {"ok": 3, "not_applicable": 3, "warning": 2, "fail": 0}.get(extracted.coding_qc_status, 1)
        candidates.append((type_rank, qc_rank, len(seq), extracted))
    if not candidates:
        return None
    return max(candidates, key=lambda item: item[:3])[3]


def location_coordinates(location: object | None) -> tuple[int, int] | None:
    if location is None:
        return None
    try:
        return int(location.start), int(location.end)
    except Exception:
        return None


def feature_coordinates(feature: object) -> tuple[int, int] | None:
    return location_coordinates(getattr(feature, "location", None))


def extract_intergenic_region(record: object, marker: Marker) -> ExtractedSequence | None:
    if not marker.flanks:
        return None
    left_name, right_name = marker.flanks
    features = getattr(record, "features", []) or []
    left_positions = [
        coordinates
        for feature in features
        if feature_type(feature) in {"gene", "cds", "trna", "rrna"}
        if feature_matches_identifier(feature, left_name)
        if (coordinates := feature_coordinates(feature)) is not None
    ]
    right_positions = [
        coordinates
        for feature in features
        if feature_type(feature) in {"gene", "cds", "trna", "rrna"}
        if feature_matches_identifier(feature, right_name)
        if (coordinates := feature_coordinates(feature)) is not None
    ]

    record_sequence = str(getattr(record, "seq", "") or "")
    record_length = len(record_sequence)
    circular = str((getattr(record, "annotations", {}) or {}).get("topology", "")).lower() == "circular"
    candidates: list[tuple[int, str, int, int]] = []
    for left_start, left_end in left_positions:
        for right_start, right_end in right_positions:
            if not (0 <= left_start < left_end <= record_length and 0 <= right_start < right_end <= record_length):
                continue
            if left_start < right_end and right_start < left_end:
                continue
            if left_end <= right_start:
                candidates.append((right_start - left_end, "linear", left_end, right_start))
            if right_end <= left_start:
                candidates.append((left_start - right_end, "linear", right_end, left_start))
            if circular:
                candidates.append((record_length - left_end + right_start, "circular", left_end, right_start))
                candidates.append((record_length - right_end + left_start, "circular", right_end, left_start))
    candidates = [candidate for candidate in candidates if 0 < candidate[0] <= marker.long_record_threshold]
    if not candidates:
        return None

    _, mode, start, end = min(candidates, key=lambda item: item[0])
    if mode == "linear":
        sequence = record_sequence[start:end].upper()
        coordinates = f"{start + 1}..{end}"
    else:
        sequence = (record_sequence[start:] + record_sequence[:end]).upper()
        coordinates = f"{start + 1}..{record_length},1..{end}"
    if not sequence:
        return None
    return build_extracted_sequence(
        sequence=sequence,
        method="flanking_features",
        note=f"{left_name}/{right_name};{coordinates};{mode}",
        marker=marker,
    )


def extract_intron_region(record: object, marker: Marker) -> ExtractedSequence | None:
    if not marker.intron_of:
        return None
    annotated: dict[tuple[int, int], object] = {}
    candidates: set[tuple[int, int, int]] = set()
    for feature in getattr(record, "features", []) or []:
        if not feature_matches_identifier(feature, marker.intron_of):
            continue
        if feature_type(feature) == "intron":
            coordinates = feature_coordinates(feature)
            if coordinates and (marker.intron_index is None or qualifier_int(feature, "number") == marker.intron_index):
                annotated[coordinates] = feature
            continue
        if feature_type(feature) not in {"cds", "trna", "rrna"}:
            continue
        if "trans_splicing" in (getattr(feature, "qualifiers", {}) or {}):
            continue
        location = getattr(feature, "location", None)
        parts = list(getattr(location, "parts", []) or [])
        if any(getattr(part, "ref", None) for part in parts):
            continue
        coordinates = sorted(
            (pair for part in parts if (pair := location_coordinates(part)) is not None),
            key=lambda pair: pair[0],
        )
        gaps = []
        for left, right in zip(coordinates, coordinates[1:]):
            if left[1] < right[0]:
                gaps.append((left[1], right[0]))
        if getattr(location, "strand", None) == -1:
            gaps.reverse()
        for index, (start, end) in enumerate(gaps, 1):
            if marker.intron_index is None or index == marker.intron_index:
                candidates.add((start, end, getattr(location, "strand", None) or 1))
    if annotated:
        if len(annotated) != 1:
            raise ValueError(f"multiple annotated introns match {marker.name}; specify an intron number")
        feature = next(iter(annotated.values()))
        sequence = str(feature.extract(getattr(record, "seq"))).upper()
        if not sequence or len(sequence) > marker.long_record_threshold:
            return None
        return build_extracted_sequence(sequence=sequence, method="annotated_intron_feature", note=str(feature.location), marker=marker, feature=feature)
    if not candidates:
        return None
    if len(candidates) != 1:
        raise ValueError(f"multiple introns match {marker.name}; specify an intron number")
    start, end, strand = next(iter(candidates))
    sequence = str(getattr(record, "seq", "") or "")[start:end].upper()
    if strand == -1:
        from Bio.Seq import Seq

        sequence = str(Seq(sequence).reverse_complement())
    if not sequence or len(sequence) > marker.long_record_threshold:
        return None
    return build_extracted_sequence(
        sequence=sequence,
        method="compound_feature_intron",
        note=f"{marker.intron_of};{start + 1}..{end}",
        marker=marker,
    )


def extract_marker_sequence(record: object, marker: Marker) -> ExtractedSequence:
    sequence = str(getattr(record, "seq", "") or "").upper()
    if not sequence:
        raise ValueError("empty GenBank sequence")
    feature = extract_best_feature(record, marker)
    if feature:
        return feature
    if marker.flanks:
        intergenic = extract_intergenic_region(record, marker)
        if intergenic:
            return intergenic
    if marker.intron_of:
        intron = extract_intron_region(record, marker)
        if intron:
            return intron
    if len(sequence) <= marker.long_record_threshold:
        if marker.intron_of or marker.flanks:
            title = str(getattr(record, "description", "") or "")
            if not marker_in_text(title, marker):
                raise ValueError(f"no boundaries annotated for {marker.name}")
        if not marker.flanks and not marker.intron_of and record_has_composite_marker_title(record, marker):
            raise ValueError(
                f"record contains multiple marker regions but no feature specific to {marker.name}"
            )
        if not record_matches_marker(record, marker):
            raise ValueError(f"record annotation does not identify {marker.name}")
        return build_extracted_sequence(sequence=sequence, method="whole_record_short", note=f"length={len(sequence)}", marker=marker)
    raise ValueError(f"could not extract {marker.name} from long GenBank record")
