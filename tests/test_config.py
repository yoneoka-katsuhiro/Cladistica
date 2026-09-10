from __future__ import annotations

import unittest

from cladistica.config import marker_map, resolve_marker


class MarkerResolutionTests(unittest.TestCase):
    def test_reported_markers_resolve_to_scientific_feature_types(self) -> None:
        markers = marker_map(["rbcL-accD", "accD", "trnK", "5.8S ribosomal RNA"])

        self.assertEqual(markers["rbcL-accD"].flanks, ("rbcL", "accD"))
        self.assertEqual(markers["rbcL-accD"].kind, "noncoding")
        self.assertEqual(markers["accD"].kind, "cds")
        self.assertEqual(markers["trnK"].kind, "trna")
        self.assertEqual(markers["5.8S"].kind, "rrna")

    def test_unicode_dash_is_canonicalized_for_intergenic_marker(self) -> None:
        marker = resolve_marker("rbcL–accD")

        self.assertEqual(marker.name, "rbcL-accD")
        self.assertEqual(marker.flanks, ("rbcL", "accD"))

    def test_trna_anticodons_are_supported_in_interval_names(self) -> None:
        marker = resolve_marker("trnL(UAA)-trnF(GAA)")

        self.assertEqual(marker.flanks, ("trnL(UAA)", "trnF(GAA)"))

    def test_any_gene_intron_is_resolved_as_a_noncoding_marker(self) -> None:
        marker = resolve_marker("rps16 intron 1")

        self.assertEqual(marker.kind, "noncoding")
        self.assertEqual(marker.intron_of, "rps16")
        self.assertEqual(marker.intron_index, 1)

    def test_unregistered_marker_is_accepted_for_annotation_driven_detection(self) -> None:
        marker = resolve_marker("PHYC")

        self.assertEqual(marker.name, "PHYC")
        self.assertEqual(marker.kind, "auto")

    def test_blank_marker_is_rejected(self) -> None:
        with self.assertRaisesRegex(ValueError, "cannot be empty"):
            resolve_marker("  ")


if __name__ == "__main__":
    unittest.main()
