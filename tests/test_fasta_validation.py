from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from cladistica.concat import concatenate_alignments, read_aligned_marker
from cladistica.io import read_fasta


class FastaValidationTests(unittest.TestCase):
    def assert_invalid_fasta(self, text: str, message: str) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "input.fasta"
            path.write_text(text, encoding="utf-8")
            with self.assertRaisesRegex(ValueError, message):
                read_fasta(path)

    def test_duplicate_headers_do_not_silently_replace_sequences(self) -> None:
        self.assert_invalid_fasta(">A\nAAAA\n>A\nCCCC\n", "Duplicate FASTA header")

    def test_sequence_before_header_is_not_discarded(self) -> None:
        self.assert_invalid_fasta("AAAA\n>A\nCCCC\n", "before.*header")

    def test_empty_header_is_rejected(self) -> None:
        self.assert_invalid_fasta(">\nAAAA\n", "Empty FASTA header")

    def test_empty_sequence_is_rejected(self) -> None:
        self.assert_invalid_fasta(">A\n>B\nAAAA\n", "Empty FASTA sequence")

    def test_final_empty_sequence_is_rejected(self) -> None:
        self.assert_invalid_fasta(">A\nAAAA\n>B\n", "Empty FASTA sequence")

    def test_duplicate_sample_ids_in_alignment_are_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "rbcL.fasta"
            path.write_text(">A|source=one\nAAAA\n>A|source=two\nCCCC\n")
            with self.assertRaisesRegex(ValueError, "Duplicate sample ID"):
                read_aligned_marker(path)

    def test_empty_sample_id_in_alignment_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "rbcL.fasta"
            path.write_text(">|source=one\nAAAA\n")
            with self.assertRaisesRegex(ValueError, "Empty sample ID"):
                read_aligned_marker(path)

    def test_empty_alignment_does_not_generate_invalid_partition(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root / "rbcL.fasta").write_text("")
            with self.assertRaisesRegex(ValueError, "Empty alignment"):
                concatenate_alignments(input_dir=root, output_dir=root / "out", markers=["rbcL"])
            self.assertFalse((root / "out" / "partition_coordinates.tsv").exists())

    def test_wrapped_sequences_and_metadata_headers_are_preserved(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "rbcL.fasta"
            path.write_text(">A|taxon=Taxon A\nAA\nCC\n\n>B\nAANN\n")
            self.assertEqual(read_fasta(path), {"A|taxon=Taxon A": "AACC", "B": "AANN"})


if __name__ == "__main__":
    unittest.main()
