from __future__ import annotations

import time
from pathlib import Path
from typing import Iterator

from .config import Marker, marker_map


def require_biopython():
    try:
        from Bio import Entrez, SeqIO  # type: ignore
    except ImportError as exc:
        raise RuntimeError(
            "Biopython is required for GenBank access. Install dependencies with "
            "`python3 -m pip install -r requirements.txt`."
        ) from exc
    return Entrez, SeqIO


def configure_entrez(email: str, api_key: str = "") -> None:
    Entrez, _ = require_biopython()
    if not email or "@" not in email:
        raise ValueError("A valid NCBI email is required. Use --email or NCBI_EMAIL.")
    Entrez.email = email
    if api_key:
        Entrez.api_key = api_key


def marker_query(markers: dict[str, Marker]) -> str:
    # Keep this intentionally short. NCBI ESearch can reject very long OR
    # expressions with multiple quoted aliases, especially for broad genera.
    tokens: list[str] = ["chloroplast", "plastid", "plastome"]
    for marker in markers.values():
        tokens.append(marker.name)
    unique = []
    seen = set()
    for token in tokens:
        normalized = token.lower()
        if normalized not in seen:
            seen.add(normalized)
            unique.append(token)
    return " OR ".join(f'"{token}"[All Fields]' for token in unique)


def marker_search_terms(query_name: str, markers: dict[str, Marker]) -> list[tuple[str, str]]:
    taxon = query_name.strip()
    base = f'"{taxon}"[Organism]'
    excluded = 'NOT "whole genome shotgun"[All Fields]'
    terms = [
        (
            "plastome",
            f'{base} AND ("chloroplast"[All Fields] OR "plastid"[All Fields] OR "plastome"[All Fields]) {excluded}',
        )
    ]
    for marker in markers.values():
        terms.append((marker.name, f'{base} AND "{marker.name}"[All Fields] {excluded}'))
    return terms


def build_search_term(query_name: str, markers: dict[str, Marker]) -> str:
    taxon = query_name.strip()
    return (
        f'"{taxon}"[Organism] '
        f"AND ({marker_query(markers)}) "
        'NOT "whole genome shotgun"[All Fields]'
    )


def search_term_ids(term: str, limit: int | None, delay_seconds: float, page_size: int = 500) -> list[str]:
    Entrez, _ = require_biopython()
    ids: list[str] = []
    retstart = 0
    while limit is None or len(ids) < limit:
        request_size = page_size if limit is None else min(page_size, limit - len(ids))
        with Entrez.esearch(
            db="nucleotide",
            term=term,
            retstart=retstart,
            retmax=request_size,
            sort="relevance",
        ) as handle:
            result = Entrez.read(handle)
        page = [str(value) for value in result.get("IdList", [])]
        ids.extend(page)
        total = int(str(result.get("Count", len(ids))))
        retstart += len(page)
        if not page or retstart >= total:
            break
        time.sleep(delay_seconds)
    return ids


def search_ids(query_name: str, markers: dict[str, Marker], retmax: int, delay_seconds: float) -> list[str]:
    seen: set[str] = set()
    ids: list[str] = []
    terms = marker_search_terms(query_name, markers)
    unlimited = retmax <= 0
    plastome_retmax = None if unlimited else max(20, min(retmax, max(20, retmax // 2)))
    marker_retmax = None if unlimited else max(5, min(25, max(5, retmax // max(len(terms), 1))))
    for label, term in terms:
        term_retmax = plastome_retmax if label == "plastome" else marker_retmax
        for accession_id in search_term_ids(term, term_retmax, delay_seconds):
            if accession_id not in seen:
                seen.add(accession_id)
                ids.append(accession_id)
        time.sleep(delay_seconds)
    if not unlimited:
        return ids[:retmax]
    return ids


def fetch_genbank_records(ids: list[str], batch_size: int, delay_seconds: float) -> Iterator[object]:
    Entrez, SeqIO = require_biopython()
    for start in range(0, len(ids), batch_size):
        batch = ids[start : start + batch_size]
        if not batch:
            continue
        with Entrez.efetch(db="nucleotide", id=",".join(batch), rettype="gb", retmode="text") as handle:
            yield from SeqIO.parse(handle, "genbank")
        time.sleep(delay_seconds)


def query_genbank_seqrecords(
    *,
    genus: str,
    outgroups: list[str],
    email: str,
    api_key: str = "",
    markers: list[str] | None = None,
    retmax: int = 500,
    batch_size: int = 50,
    delay_seconds: float = 0.34,
    log_path: Path | None = None,
) -> list[tuple[str, object]]:
    configure_entrez(email, api_key)
    marker_defs = marker_map(markers)
    query_names = [genus, *outgroups]
    rows: list[str] = ["query_name\tids_found\trecords_parsed"]
    records: list[tuple[str, object]] = []
    for query_name in query_names:
        ids = search_ids(query_name, marker_defs, retmax=retmax, delay_seconds=delay_seconds)
        parsed = 0
        for record in fetch_genbank_records(ids, batch_size=batch_size, delay_seconds=delay_seconds):
            records.append((query_name, record))
            parsed += 1
        rows.append(f"{query_name}\t{len(ids)}\t{parsed}")
    if log_path:
        log_path.parent.mkdir(parents=True, exist_ok=True)
        log_path.write_text("\n".join(rows) + "\n", encoding="utf-8")
    return records


def fetch_one_genbank(accession: str, email: str, api_key: str = "") -> object:
    configure_entrez(email, api_key)
    Entrez, SeqIO = require_biopython()
    with Entrez.efetch(db="nucleotide", id=accession, rettype="gb", retmode="text") as handle:
        return SeqIO.read(handle, "genbank")
