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


class FakeCompoundLocation:
    def __init__(self, parts: list[FakeLocation]) -> None:
        self.parts = parts
        self.start = parts[0].start
        self.end = parts[-1].end

    def __str__(self) -> str:
        return "join(" + ",".join(str(part) for part in self.parts) + ")"


class FakeFeature:
    def __init__(
        self,
        feature_type: str,
        sequence: str,
        qualifiers: dict[str, list[str]],
        location_text: str = "",
        start: int = 0,
        end: int | None = None,
    ) -> None:
        self.type = feature_type
        self.sequence = sequence
        self.qualifiers = qualifiers
        self.location = FakeLocation(start, end if end is not None else start + len(sequence), location_text)

    def extract(self, _sequence: str) -> str:
        return self.sequence


class FakeRecord:
    def __init__(self, sequence: str, features: list[FakeFeature], description: str = "") -> None:
        self.seq = sequence
        self.features = features
        self.name = "record"
        self.id = "A1"
        self.description = description


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

    def test_reported_cds_trna_and_rrna_features_are_extracted_by_type(self) -> None:
        cases = (
            ("accD", "CDS", {"gene": ["accD"]}, "cds", "ok"),
            ("trnK", "tRNA", {"gene": ["trnK-UUU"], "product": ["tRNA-Lys"]}, "trna", "not_applicable"),
            ("5.8S ribosomal RNA", "rRNA", {"product": ["5.8S ribosomal RNA"]}, "rrna", "not_applicable"),
        )
        for requested, feature_type, qualifiers, expected_kind, expected_qc in cases:
            with self.subTest(marker=requested):
                marker = next(iter(marker_map([requested]).values()))
                feature = FakeFeature(feature_type, "GCT" * 100, qualifiers)
                extracted = extract_marker_sequence(FakeRecord("N" * 1000, [feature]), marker)

                self.assertEqual(extracted.marker_kind, expected_kind)
                self.assertEqual(extracted.coding_qc_status, expected_qc)
                self.assertEqual(extracted.method, f"annotated_{feature_type.lower()}_feature")

    def test_intergenic_marker_is_extracted_between_annotated_flanks(self) -> None:
        marker = marker_map(["rbcL-accD"])["rbcL-accD"]
        sequence = "A" * 10 + "CGTACG" + "T" * 10
        features = [
            FakeFeature("CDS", "A" * 10, {"gene": ["rbcL"]}, start=0, end=10),
            FakeFeature("CDS", "T" * 10, {"gene": ["accD"]}, start=16, end=26),
        ]
        extracted = extract_marker_sequence(FakeRecord(sequence, features), marker)

        self.assertEqual(extracted.sequence, "CGTACG")
        self.assertEqual(extracted.method, "flanking_features")
        self.assertIn("rbcL/accD", extracted.note)

    def test_short_marker_specific_interval_record_is_accepted_whole(self) -> None:
        marker = marker_map(["rbcL-accD"])["rbcL-accD"]
        record = FakeRecord("ACGT" * 100, [], "rbcL-accD intergenic spacer")

        extracted = extract_marker_sequence(record, marker)

        self.assertEqual(extracted.method, "whole_record_short")
        self.assertEqual(len(extracted.sequence), 400)

    def test_trna_flanks_match_standard_anticodon_suffixes(self) -> None:
        marker = marker_map(["trnL-F"])["trnL-F"]
        sequence = "A" * 8 + "CGTA" + "T" * 8
        features = [
            FakeFeature("tRNA", "A" * 8, {"gene": ["trnL-UAA"]}, start=0, end=8),
            FakeFeature("tRNA", "T" * 8, {"gene": ["trnF-GAA"]}, start=12, end=20),
        ]
        extracted = extract_marker_sequence(FakeRecord(sequence, features), marker)

        self.assertEqual(extracted.sequence, "CGTA")
        self.assertEqual(extracted.method, "flanking_features")

    def test_intron_is_derived_from_compound_annotated_gene(self) -> None:
        marker = marker_map(["rps16 intron"])["rps16 intron"]
        sequence = "A" * 10 + "CGTACG" + "T" * 10
        feature = FakeFeature("CDS", sequence, {"gene": ["rps16"]})
        feature.location = FakeCompoundLocation([FakeLocation(0, 10), FakeLocation(16, 26)])

        extracted = extract_marker_sequence(FakeRecord(sequence, [feature]), marker)

        self.assertEqual(extracted.sequence, "CGTACG")
        self.assertEqual(extracted.method, "compound_feature_intron")

    def test_long_trnk_intron_is_not_mixed_with_mature_trna(self) -> None:
        marker = marker_map(["trnK"])["trnK"]
        feature = FakeFeature("gene", "A" * 2400, {"gene": ["trnK"]})
        record = FakeRecord("A" * 2400, [feature], "trnK gene and matK gene")

        with self.assertRaisesRegex(ValueError, "could not extract trnK"):
            extract_marker_sequence(record, marker)

    def test_unregistered_cds_marker_uses_codon_qc_after_annotation_match(self) -> None:
        marker = marker_map(["PHYC"])["PHYC"]
        feature = FakeFeature("CDS", "GCTGCTGCT", {"gene": ["PHYC"]})
        extracted = extract_marker_sequence(FakeRecord("N" * 100, [feature]), marker)

        self.assertEqual(extracted.marker_kind, "cds")
        self.assertEqual(extracted.coding_qc_status, "ok")

    def test_broad_multiregion_misc_rna_is_not_mislabeled_as_58s(self) -> None:
        marker = marker_map(["5.8S"])["5.8S"]
        feature = FakeFeature(
            "misc_RNA",
            "A" * 6000,
            {
                "note": [
                    "contains small subunit ribosomal RNA, internal transcribed spacer 1, "
                    "5.8S ribosomal RNA, internal transcribed spacer 2, and large subunit ribosomal RNA"
                ]
            },
        )
        record = FakeRecord(
            "A" * 6000,
            [feature],
            "18S ribosomal RNA; ITS1, 5.8S ribosomal RNA, ITS2; and 28S ribosomal RNA",
        )

        with self.assertRaisesRegex(ValueError, "could not extract 5.8S"):
            extract_marker_sequence(record, marker)

    def test_short_composite_record_without_subfeature_is_rejected(self) -> None:
        marker = marker_map(["5.8S"])["5.8S"]
        record = FakeRecord(
            "A" * 700,
            [],
            "ITS1, 5.8S ribosomal RNA, and ITS2 complete sequence",
        )

        with self.assertRaisesRegex(ValueError, "multiple marker regions"):
            extract_marker_sequence(record, marker)


if __name__ == "__main__":
    unittest.main()
