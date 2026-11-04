from __future__ import annotations

import unittest

from cladistica.download import sample_header, validate_selected_accessions
from cladistica.io import sample_id_from_header


class DownloadIdentityTests(unittest.TestCase):
    def test_taxon_only_row_keeps_one_sample_across_markers(self) -> None:
        row = {"taxon": "Taxon alba", "rbcL": "AA1", "matK": "BB1"}
        first = sample_id_from_header(sample_header(row, "rbcL", "AA1"))
        second = sample_id_from_header(sample_header(row, "matK", "BB1"))
        self.assertEqual(first, second)
        self.assertEqual(first, "Taxon_alba")

    def test_explicit_voucher_stays_in_fallback_identity(self) -> None:
        row = {"taxon": "Taxon alba", "voucher_or_isolate": "K 12"}
        self.assertEqual(sample_id_from_header(sample_header(row, "rbcL", "AA1")), "Taxon_alba_K_12")

    def test_whitespace_is_trimmed_from_explicit_sample_id(self) -> None:
        self.assertEqual(sample_header({"sample_id": "  A  "}, "rbcL", "AA1").split("|")[0], "A")

    def test_normalized_taxon_collisions_are_rejected_before_download(self) -> None:
        rows = [{"taxon": "Taxon alba", "rbcL": "AA1"}, {"taxon": "Taxon_alba", "rbcL": "BB1"}]
        with self.assertRaisesRegex(ValueError, "Duplicate FASTA sample ID"):
            validate_selected_accessions(rows, ["rbcL"])

    def test_pipe_in_sample_id_is_rejected(self) -> None:
        with self.assertRaisesRegex(ValueError, "Invalid sample_id"):
            validate_selected_accessions([{"sample_id": "A|B", "rbcL": "AA1"}], ["rbcL"])

    def test_linebreak_in_sample_id_is_rejected(self) -> None:
        with self.assertRaisesRegex(ValueError, "Invalid sample_id"):
            validate_selected_accessions([{"sample_id": "A\nB", "rbcL": "AA1"}], ["rbcL"])

    def test_metadata_cannot_inject_a_new_fasta_record(self) -> None:
        header = sample_header({"sample_id": "A", "taxon": "Taxon\n>injected"}, "rbcL", "AA1")
        self.assertNotIn("\n", header)


if __name__ == "__main__":
    unittest.main()
