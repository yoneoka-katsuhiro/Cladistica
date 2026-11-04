"""Offline parser round-trips for diagnosing an installed/frozen package."""
from __future__ import annotations

import io
from pathlib import Path
import tempfile

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from .io import read_table, write_xlsx


def self_check() -> None:
    record = SeqRecord(Seq("ACGT"), id="DEMO.1", name="DEMO", description="Synthetic installation check")
    record.annotations["molecule_type"] = "DNA"
    buffer = io.StringIO()
    SeqIO.write(record, buffer, "genbank")
    buffer.seek(0)
    parsed = SeqIO.read(buffer, "genbank")
    if parsed.id != record.id or parsed.seq != record.seq:
        raise RuntimeError("GenBank parser round-trip failed.")
    print("Biopython GenBank read/write: OK")
    with tempfile.TemporaryDirectory(prefix="cladistica_check_") as temporary:
        path = Path(temporary) / "table.xlsx"
        expected = [{"sample_id": "demo", "rbcL": "DEMO.1"}]
        if not write_xlsx(path, expected, ["sample_id", "rbcL"]):
            raise RuntimeError("openpyxl is unavailable.")
        rows, columns = read_table(path)
        if rows != expected or columns != ["sample_id", "rbcL"]:
            raise RuntimeError("Excel parser round-trip failed.")
    print("openpyxl Excel read/write: OK")
    print("Offline installation check passed. External phylogenetic programs were not tested.")
