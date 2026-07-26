from __future__ import annotations

import re
from dataclasses import dataclass


@dataclass(frozen=True)
class Marker:
    name: str
    aliases: tuple[str, ...]
    expected_min_length: int
    long_record_threshold: int
    kind: str = "cds"


DEFAULT_MARKERS: tuple[Marker, ...] = (
    Marker("rbcL", ("rbcl", "rbc-l", "ribulose-1,5-bisphosphate carboxylase/oxygenase large subunit"), 900, 3000),
    Marker("matK", ("matk", "maturase k"), 700, 3000),
    Marker("atpA", ("atpa", "atp-a", "ATP synthase CF1 alpha subunit", "ATPase alpha subunit"), 900, 3000),
    Marker("atpB", ("atpb", "atp-b", "atp synthase cf1 beta subunit"), 900, 3000),
    Marker("ndhF", ("ndhf", "ndh-f", "NADH dehydrogenase subunit F"), 1200, 4000),
    Marker("rps4", ("rps-4", "ribosomal protein S4"), 450, 2500),
    Marker("rps16", ("rps-16", "ribosomal protein S16"), 450, 2500),
    Marker("psaA", ("psaa", "psa-a", "photosystem I P700 apoprotein A1"), 1200, 4000),
    Marker("psaB", ("psab", "psa-b", "photosystem I P700 apoprotein A2"), 1200, 4000),
    Marker("psbA", ("psba", "psb-a", "photosystem II protein D1"), 700, 3000),
    Marker("psbB", ("psbb", "psb-b", "photosystem II CP47 protein"), 1200, 4000),
    Marker("psbC", ("psbc", "psb-c", "photosystem II CP43 protein"), 1200, 4000),
    Marker("psbD", ("psbd", "psb-d", "photosystem II protein D2"), 900, 3500),
    Marker("ycf1", ("ycf-1", "hypothetical chloroplast RF1", "hypothetical protein ycf1"), 2000, 7000),
    Marker("ycf2", ("ycf-2", "hypothetical chloroplast RF2", "hypothetical protein ycf2"), 3000, 9000),
    Marker("trnL-F", ("trnl-f", "trnl trnf", "trnl-trnf", "trnl intron", "trnl-trnf intergenic spacer", "trnl-f intergenic spacer"), 250, 5000, "noncoding"),
    Marker("rps4-trnS", ("rps4-trns", "rps4-trns intergenic spacer", "rps4 trns", "rps4-trns spacer"), 250, 5000, "noncoding"),
)

UNCERTAIN_NAME_PATTERNS: tuple[str, ...] = (
    r"\bsp\.?\b",
    r"\bcf\.?\b",
    r"\baff\.?\b",
    r"\bindet\.?\b",
    r"\bunidentified\b",
    r"\bunclassified\b",
    r"\buncultured\b",
    r"\benvironmental sample\b",
    r"\bmetagenome\b",
    r"\bhybrid\b",
    r"\bcultivar\b",
    r"\bx\s+[A-Z_a-z-]+",
)


def marker_map(markers: list[str] | tuple[str, ...] | None = None) -> dict[str, Marker]:
    by_name = {marker.name: marker for marker in DEFAULT_MARKERS}
    if not markers:
        return by_name
    selected: dict[str, Marker] = {}
    for name in markers:
        if name in by_name:
            selected[name] = by_name[name]
            continue
        key = normalize_marker_text(name)
        match = next(
            (
                marker
                for marker in DEFAULT_MARKERS
                if key in {normalize_marker_text(marker.name), *[normalize_marker_text(alias) for alias in marker.aliases]}
            ),
            None,
        )
        if not match:
            raise ValueError(f"Unknown marker: {name}")
        selected[match.name] = match
    return selected


def normalize_marker_text(value: str) -> str:
    return re.sub(r"[^a-z0-9]+", "", value.lower())


def safe_token(value: str) -> str:
    token = re.sub(r"[^A-Za-z0-9._-]+", "_", value.strip())
    token = re.sub(r"_+", "_", token).strip("._-")
    return token or "unnamed"
