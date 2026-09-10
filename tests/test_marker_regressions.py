from __future__ import annotations

import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from Bio.Seq import Seq
from Bio.SeqFeature import CompoundLocation, SeqFeature, SimpleLocation
from Bio.SeqRecord import SeqRecord

from cladistica.config import resolve_marker
from cladistica.download import download_fasta_from_accession_table
from cladistica.extraction import extract_marker_sequence, record_matches_marker
from cladistica.io import read_fasta


def feature(gene: str, start: int, end: int, kind: str = "gene", **qualifiers: object) -> SeqFeature:
    return SeqFeature(SimpleLocation(start, end), type=kind, qualifiers={"gene": [gene], **qualifiers})


class MarkerRegressionTests(unittest.TestCase):
    def test_linear_record_never_wraps_around_origin(self) -> None:
        record = SeqRecord(Seq("A" * 20 + "C" * 60 + "T" * 20))
        record.features = [feature("rbcL", 5, 20), feature("accD", 80, 95)]
        extracted = extract_marker_sequence(record, resolve_marker("rbcL-accD"))
        self.assertEqual(extracted.sequence, "C" * 60)

    def test_circular_record_can_extract_origin_spacer(self) -> None:
        record = SeqRecord(Seq("A" * 20 + "C" * 60 + "T" * 20), annotations={"topology": "circular"})
        record.features = [feature("rbcL", 5, 20), feature("accD", 80, 95)]
        extracted = extract_marker_sequence(record, resolve_marker("rbcL-accD"))
        self.assertEqual(extracted.sequence, "T" * 5 + "A" * 5)

    def test_numbered_annotated_intron_selects_requested_number(self) -> None:
        record = SeqRecord(Seq("AAAACCCCGGGGTTTT"))
        record.features = [
            feature("clpP", 0, 4, "intron", number=["1"]),
            feature("clpP", 4, 8, "intron", number=["2"]),
        ]
        self.assertEqual(extract_marker_sequence(record, resolve_marker("clpP intron 2")).sequence, "CCCC")
        with self.assertRaisesRegex(ValueError, "multiple annotated introns"):
            extract_marker_sequence(record, resolve_marker("clpP intron"))

    def test_reverse_strand_intron_numbering_and_orientation(self) -> None:
        record = SeqRecord(Seq("AAAACCCCGGGGTTTTAAAA"))
        location = CompoundLocation([
            SimpleLocation(16, 20, strand=-1),
            SimpleLocation(8, 12, strand=-1),
            SimpleLocation(0, 4, strand=-1),
        ])
        record.features = [SeqFeature(location, type="CDS", qualifiers={"gene": ["clpP"]})]
        self.assertEqual(extract_marker_sequence(record, resolve_marker("clpP intron 1")).sequence, "AAAA")
        self.assertEqual(extract_marker_sequence(record, resolve_marker("clpP intron 2")).sequence, "GGGG")

    def test_missing_intron_is_not_substituted_with_whole_gene(self) -> None:
        record = SeqRecord(Seq("A" * 100), description="rps16 gene")
        record.features = [feature("rps16", 0, 100, "CDS")]
        with self.assertRaisesRegex(ValueError, "no boundaries"):
            extract_marker_sequence(record, resolve_marker("rps16 intron"))

    def test_short_names_and_gene_prefixes_have_token_boundaries(self) -> None:
        for requested, title, expected in [("ITS2", "ITS2 complete sequence", True), ("5.8S", "5.8S gene", True), ("rps4", "rps40 gene", False), ("ITS", "ITS2 gene", False)]:
            with self.subTest(marker=requested, title=title):
                record = SeqRecord(Seq("A" * 100), description=title)
                self.assertEqual(record_matches_marker(record, resolve_marker(requested)), expected)

    def test_tRNA_anticodon_is_not_interpreted_as_a_spacer(self) -> None:
        marker = resolve_marker("trnS-GGA")
        self.assertEqual(marker.kind, "trna")
        self.assertIsNone(marker.flanks)

    def test_invalid_marker_identifiers_are_rejected(self) -> None:
        for name in ("../gene", "gene\nname", "///", "rps16 intron 0"):
            with self.subTest(name=name), self.assertRaises(ValueError):
                resolve_marker(name)

    def test_selected_csv_can_use_an_rRNA_product_alias(self) -> None:
        record = SeqRecord(Seq("A" * 160))
        record.features = [SeqFeature(SimpleLocation(0, 160), type="rRNA", qualifiers={"product": ["5.8S ribosomal RNA"]})]
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            table = root / "selected.csv"
            table.write_text("taxon,sample_id,5.8S ribosomal RNA\nTest species,A,TEST1\n")
            with patch("cladistica.download.fetch_one_genbank", return_value=record):
                result = download_fasta_from_accession_table(accession_table=table, output_dir=root / "out", email="test@example.org", markers=["5.8S ribosomal RNA"])
            self.assertEqual(result.records_written, 1)
            self.assertEqual(list(read_fasta(root / "out" / "5.8S.fasta").values()), ["A" * 160])
