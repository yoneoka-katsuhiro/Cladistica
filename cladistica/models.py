from __future__ import annotations

from dataclasses import dataclass, field


@dataclass
class MarkerHit:
    marker: str
    accession: str
    length: int = 0
    ambiguous_bases: int = 0
    marker_kind: str = ""
    source_record_id: str = ""
    feature_note: str = ""
    extraction_method: str = ""
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
    sequence: str = ""

    @property
    def is_clean(self) -> bool:
        return self.ambiguous_bases == 0 and self.coding_qc_status not in {"fail"}


@dataclass
class AccessionCandidate:
    taxon: str
    sample_id: str
    voucher_or_isolate: str = ""
    marker_hits: dict[str, MarkerHit] = field(default_factory=dict)
    publication_status: str = "unknown"
    publication_score: int = 0
    pubmed_ids: set[str] = field(default_factory=set)
    reference_titles: set[str] = field(default_factory=set)
    source_record_ids: set[str] = field(default_factory=set)
    query_name: str = ""

    def marker_count(self, markers: list[str]) -> int:
        return sum(1 for marker in markers if marker in self.marker_hits)

    def clean_marker_count(self, markers: list[str]) -> int:
        return sum(1 for marker in markers if marker in self.marker_hits and self.marker_hits[marker].is_clean)

    def total_length(self, markers: list[str]) -> int:
        return sum(self.marker_hits[marker].length for marker in markers if marker in self.marker_hits)

    def ambiguous_bases(self, markers: list[str]) -> int:
        return sum(self.marker_hits[marker].ambiguous_bases for marker in markers if marker in self.marker_hits)

    def selected_accessions(self, markers: list[str]) -> dict[str, str]:
        return {marker: self.marker_hits[marker].accession if marker in self.marker_hits else "" for marker in markers}

    def sort_key(self, markers: list[str]) -> tuple:
        return (
            self.marker_count(markers),
            self.clean_marker_count(markers),
            self.publication_score,
            self.total_length(markers),
            -self.ambiguous_bases(markers),
            self.taxon.lower(),
            self.sample_id.lower(),
        )


@dataclass
class PipelineResult:
    output_dir: str
    records_written: int = 0
    warnings: list[str] = field(default_factory=list)
