from __future__ import annotations

import unittest

from cladistica.config import marker_map
from cladistica.extraction import extract_marker_sequence


class FakeLocation:
    def __init__(self, start: int, end: int, text: str = "") -> None:
        self.start = start
        self.end = end
        self.text = text or f"[{start}:{end}](+)"

    def __str__(self) -> str:
        return self.text


class FakeFeature:
    def __init__(self, feature_type: str, sequence: str, qualifiers: dict[str, list[str]], location_text: str = "") -> None:
        self.type = feature_type
        self.sequence = sequence
        self.qualifiers = qualifiers
        self.location = FakeLocation(0, len(sequence), location_text)

    def extract(self, _sequence: str) -> str:
        return self.sequence


class FakeRecord:
    def __init__(self, sequence: str, features: list[FakeFeature]) -> None:
        self.seq = sequence
        self.features = features


class ExtractionTests(unittest.TestCase):
    def test_partial_cds_without_start_or_terminal_stop_is_allowed(self) -> None:
        marker = marker_map(["rbcL"])["rbcL"]
        partial_sequence = "GCT" * 100
        feature = FakeFeature(
            "CDS",
            partial_sequence,
            {"gene": ["rbcL"], "codon_start": ["1"], "transl_table": ["11"]},
            "<1..>300",
        )
        extracted = extract_marker_sequence(FakeRecord("N" * 100, [feature]), marker)

        self.assertEqual(extracted.method, "annotated_cds_feature")
        self.assertEqual(extracted.coding_qc_status, "ok")
        self.assertTrue(extracted.is_partial)
        self.assertFalse(extracted.internal_stop)

    def test_cds_internal_stop_is_flagged(self) -> None:
        marker = marker_map(["rbcL"])["rbcL"]
        feature = FakeFeature("CDS", "GCTTAAGCT", {"gene": ["rbcL"], "transl_table": ["11"]})
        extracted = extract_marker_sequence(FakeRecord("N" * 100, [feature]), marker)

        self.assertEqual(extracted.coding_qc_status, "fail")
        self.assertTrue(extracted.internal_stop)

    def test_noncoding_marker_skips_codon_qc(self) -> None:
        marker = marker_map(["trnL-F"])["trnL-F"]
        feature = FakeFeature("misc_feature", "ATGTAATGA", {"note": ["trnL-F intergenic spacer"]})
        extracted = extract_marker_sequence(FakeRecord("N" * 100, [feature]), marker)

        self.assertEqual(extracted.marker_kind, "noncoding")
        self.assertEqual(extracted.coding_qc_status, "not_applicable")
        self.assertFalse(extracted.internal_stop)

    def test_cds_feature_is_preferred_over_gene_feature(self) -> None:
        marker = marker_map(["rbcL"])["rbcL"]
        gene = FakeFeature("gene", "ATGGGGGGGTAA", {"gene": ["rbcL"]})
        cds = FakeFeature("CDS", "GCTGCTGCT", {"gene": ["rbcL"]})
        extracted = extract_marker_sequence(FakeRecord("N" * 100, [gene, cds]), marker)

        self.assertEqual(extracted.sequence, "GCTGCTGCT")
        self.assertEqual(extracted.method, "annotated_cds_feature")


if __name__ == "__main__":
    unittest.main()
