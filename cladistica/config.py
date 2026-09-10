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
    flanks: tuple[str, str] | None = None
    intron_of: str | None = None
    intron_index: int | None = None


DEFAULT_MARKERS: tuple[Marker, ...] = (
    Marker("rbcL", ("rbcl", "rbc-l", "ribulose-1,5-bisphosphate carboxylase/oxygenase large subunit"), 900, 3000),
    Marker("matK", ("matk", "maturase k"), 700, 3000),
    Marker("atpA", ("atpa", "atp-a", "ATP synthase CF1 alpha subunit", "ATPase alpha subunit"), 900, 3000),
    Marker("atpB", ("atpb", "atp-b", "ATP synthase CF1 beta subunit"), 900, 3000),
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
    Marker(
        "trnL-F",
        ("trnl-f", "trnl trnf", "trnl-trnf", "trnl-trnf intergenic spacer", "trnl-f intergenic spacer"),
        250,
        5000,
        "noncoding",
        ("trnL", "trnF"),
    ),
    Marker(
        "rps4-trnS",
        ("rps4-trns", "rps4-trns intergenic spacer", "rps4 trns", "rps4-trns spacer"),
        250,
        5000,
        "noncoding",
        ("rps4", "trnS"),
    ),
)


# Common plant loci get stable names and types. The dynamic resolver below
# handles valid markers outside this curated registry.
KNOWN_MARKERS: tuple[Marker, ...] = DEFAULT_MARKERS + (
    Marker("accD", ("accd", "acetyl-CoA carboxylase beta subunit"), 900, 4000),
    Marker("rpoB", ("rpob", "DNA-directed RNA polymerase beta subunit"), 500, 3500),
    Marker("rpoC1", ("rpoc1", "DNA-directed RNA polymerase beta' subunit"), 500, 3500),
    Marker("ndhJ", ("ndhj", "NADH dehydrogenase subunit J"), 300, 2500),
    Marker("trnK", ("trnk", "trnK-UUU", "tRNA-Lys", "transfer RNA-Lys"), 50, 500, "trna"),
    Marker("5.8S", ("5.8s", "5.8S rRNA", "5.8S rDNA", "5.8S ribosomal RNA", "rrn5.8S"), 100, 2500, "rrna"),
    Marker("18S", ("18s", "18S rRNA", "18S ribosomal RNA", "small subunit ribosomal RNA"), 1000, 4000, "rrna"),
    Marker("26S", ("26s", "26S rRNA", "26S ribosomal RNA"), 1500, 5000, "rrna"),
    Marker("28S", ("28s", "28S rRNA", "28S ribosomal RNA", "large subunit ribosomal RNA"), 1500, 5000, "rrna"),
    Marker("ITS", ("nrITS", "ribosomal internal transcribed spacer", "internal transcribed spacer"), 200, 3000, "noncoding"),
    Marker("ITS1", ("nrITS1", "internal transcribed spacer 1"), 100, 2000, "noncoding"),
    Marker("ITS2", ("nrITS2", "internal transcribed spacer 2"), 100, 2000, "noncoding"),
    Marker("ETS", ("nrETS", "external transcribed spacer"), 100, 3000, "noncoding"),
    Marker("trnH-psbA", ("psbA-trnH", "trnH psbA spacer", "trnH-psbA intergenic spacer"), 100, 5000, "noncoding", ("trnH", "psbA")),
    Marker("atpF-atpH", ("atpF atpH spacer", "atpF-atpH intergenic spacer"), 100, 5000, "noncoding", ("atpF", "atpH")),
    Marker("psbK-psbI", ("psbK psbI spacer", "psbK-psbI intergenic spacer"), 100, 5000, "noncoding", ("psbK", "psbI")),
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


INTERGENIC_PATTERN = re.compile(
    r"^\s*([A-Za-z][A-Za-z0-9_.]*(?:\([A-Za-z0-9]+\))?)\s*[-\u2010-\u2015]\s*"
    r"([A-Za-z][A-Za-z0-9_.]*(?:\([A-Za-z0-9]+\))?)\s*"
    r"(?:(?:intergenic\s+)?(?:spacer|region|igs))?\s*$",
    flags=re.IGNORECASE,
)

INTRON_PATTERN = re.compile(
    r"^\s*([A-Za-z][A-Za-z0-9_.]*(?:\([A-Za-z0-9]+\))?)\s+"
    r"intron(?:\s+(\d+))?\s*$",
    flags=re.IGNORECASE,
)

CODING_GENE_PATTERN = re.compile(
    r"^(?:accd|atp[a-z0-9]+|ccsa|cema|ch[lbc][a-z0-9]+|clpp\d*|"
    r"cob|cox\d+|infa|matk|matr|nad\d+[a-z]*|ndh[a-z0-9]+|"
    r"pet[a-z0-9]+|psa[a-z0-9]+|psb[a-z0-9]+|rbcl|rpl\d+[a-z]*|"
    r"rpo[a-z0-9]+|rps\d+[a-z]*|ycf\d+[a-z]*)$",
    flags=re.IGNORECASE,
)


def normalize_marker_text(value: str) -> str:
    return re.sub(r"[^a-z0-9]+", "", value.lower())


def clean_marker_name(value: str) -> str:
    if any(ord(character) < 32 for character in str(value)):
        raise ValueError(f"Marker name contains a control character: {value!r}")
    name = " ".join(str(value).strip().split())
    name = re.sub(r"[\u2010-\u2015]", "-", name)
    if not name:
        raise ValueError("Marker name cannot be empty")
    if not normalize_marker_text(name) or any(char in name for char in '\\";|'):
        raise ValueError(f"Invalid marker identifier: {value!r}")
    if any(ord(character) < 32 for character in name):
        raise ValueError(f"Marker name contains a control character: {value!r}")
    return name


def _known_marker(name: str) -> Marker | None:
    key = normalize_marker_text(name)
    for marker in KNOWN_MARKERS:
        marker_keys = {
            normalize_marker_text(marker.name),
            *(normalize_marker_text(alias) for alias in marker.aliases),
        }
        if key in marker_keys:
            return marker
    return None


def _expand_flank_shorthand(left: str, right: str) -> tuple[str, str]:
    if left.lower().startswith("trn") and len(right) == 1:
        right = "trn" + right.upper()
    return left, right


def infer_marker_kind(name: str) -> str:
    normalized = normalize_marker_text(name)
    if "ribosomalrna" in normalized or "rrna" in normalized:
        return "rrna"
    if re.match(r"^(?:58|18|26|28)s(?:ribosomalrna|rrna)?$", normalized):
        return "rrna"
    if normalized.startswith("trn") or "transferrna" in normalized:
        return "trna"
    if any(token in normalized for token in ("intergenic", "spacer", "intron")):
        return "noncoding"
    if normalized.startswith(("its", "nrits", "ets", "nrets")) or "transcribedspacer" in normalized:
        return "noncoding"
    if CODING_GENE_PATTERN.fullmatch(normalized):
        return "cds"
    return "auto"


def resolve_marker(value: str) -> Marker:
    """Resolve a curated alias or construct an annotation-driven marker."""
    name = clean_marker_name(value)
    known = _known_marker(name)
    if known:
        return known
    if "/" in name:
        raise ValueError(f"Invalid marker identifier: {value!r}")

    intron = INTRON_PATTERN.fullmatch(name)
    if intron:
        gene, index_text = intron.groups()
        intron_index = int(index_text) if index_text else None
        if intron_index is not None and intron_index < 1:
            raise ValueError("Intron numbers start at 1")
        canonical = f"{gene} intron" + (f" {intron_index}" if intron_index else "")
        return Marker(
            canonical,
            (name,),
            0,
            10_000,
            "noncoding",
            intron_of=gene,
            intron_index=intron_index,
        )

    if re.fullmatch(r"trn[A-Za-z](?:[-(][ACGTUacgtu]{3}\)?)?", name, re.I):
        return Marker(name, (name,), 0, 500, "trna")

    interval = INTERGENIC_PATTERN.fullmatch(name)
    if interval:
        left, right = _expand_flank_shorthand(*interval.groups())
        canonical = f"{left}-{right}"
        return Marker(
            canonical,
            (name, f"{left} {right} intergenic spacer"),
            0,
            10_000,
            "noncoding",
            (left, right),
        )

    kind = infer_marker_kind(name)
    return Marker(name, (name,), 0, 10_000, kind)


def marker_map(markers: list[str] | tuple[str, ...] | None = None) -> dict[str, Marker]:
    selected_markers = DEFAULT_MARKERS if not markers else tuple(
        resolve_marker(name) for name in markers
    )
    selected: dict[str, Marker] = {}
    normalized_names: dict[str, str] = {}
    for marker in selected_markers:
        key = normalize_marker_text(marker.name)
        previous = normalized_names.get(key)
        if previous and previous != marker.name:
            raise ValueError(
                f"Marker names resolve to the same identifier: {previous!r} and {marker.name!r}"
            )
        normalized_names[key] = marker.name
        selected[marker.name] = marker
    return selected


def marker_query_terms(marker: Marker) -> tuple[str, ...]:
    """Return a short set of useful Entrez phrases for one marker."""
    terms = [marker.name]
    if marker.flanks:
        terms.extend(alias for alias in marker.aliases[:1] if len(alias) >= 4)
        terms.extend(marker.flanks)
    elif marker.intron_of:
        terms.extend(marker.aliases)
        terms.append(marker.intron_of)
    else:
        terms.extend(alias for alias in marker.aliases if len(alias) >= 4)
    unique: list[str] = []
    seen: set[str] = set()
    for term in terms:
        cleaned = str(term).replace('"', "").strip()
        key = normalize_marker_text(cleaned)
        if cleaned and key not in seen:
            seen.add(key)
            unique.append(cleaned)
        if len(unique) == 4:
            break
    return tuple(unique)


def safe_token(value: str) -> str:
    token = re.sub(r"[^A-Za-z0-9._-]+", "_", value.strip())
    token = re.sub(r"_+", "_", token).strip("._-")
    return token or "unnamed"
