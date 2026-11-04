from pathlib import Path
import tempfile
import unittest

from cladistica.align import collect_marker_records


class AlignmentIdentityTests(unittest.TestCase):
    def test_duplicate_sample_ids_within_one_source_are_rejected(self):
        for source in ("genbank", "user"):
            with self.subTest(source=source), tempfile.TemporaryDirectory() as temporary:
                root = Path(temporary)
                folder = root / source
                folder.mkdir()
                (folder / "rbcL.fasta").write_text(">A|record=one\nAAAA\n>A|record=two\nCCCC\n")
                with self.assertRaisesRegex(ValueError, "Duplicate sample ID"):
                    collect_marker_records("rbcL", root / "genbank", root / "user")

    def test_empty_sample_id_is_rejected_before_alignment(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            (root / "rbcL.fasta").write_text(">|record=one\nAAAA\n")
            with self.assertRaisesRegex(ValueError, "Empty sample ID"):
                collect_marker_records("rbcL", root, root / "unused")

    def test_user_override_of_genbank_is_preserved_and_reported(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            for source, sequence in (("genbank", "AAAA"), ("user", "CCCC")):
                (root / source).mkdir()
                (root / source / "rbcL.fasta").write_text(f">A|source={source}\n{sequence}\n")
            records, conflicts = collect_marker_records("rbcL", root / "genbank", root / "user")
            self.assertEqual(records, {"A|source=user": "CCCC"})
            self.assertEqual(conflicts, [{"marker": "rbcL", "sample_id": "A", "kept_source": "user", "discarded_header": "A|source=genbank"}])


if __name__ == "__main__":
    unittest.main()
