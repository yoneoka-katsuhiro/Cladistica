from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from cladistica.accessions import (
    ExtractedMarkerRecord,
    RecordLike,
    build_accession_table_from_records,
    candidates_from_extracted_records,
    candidates_from_records,
    select_representatives,
)
from cladistica.config import marker_map
from cladistica.io import read_csv


class AccessionSelectionTests(unittest.TestCase):
    def record(
        self,
        accession: str,
        organism: str,
        voucher: str,
        markers: set[str],
        sequence: str,
        publication_score: int = 0,
    ) -> RecordLike:
        return RecordLike(
            accession=accession,
            organism=organism,
            description="",
            sequence=sequence,
            voucher_or_isolate=voucher,
            publication_status="published_or_accepted" if publication_score else "unpublished",
            publication_score=publication_score,
            pubmed_ids=set(),
            reference_titles=set(),
            query_name="Hymenasplenium",
            marker_names=markers,
        )

    def test_selects_one_best_sample_per_taxon_by_coverage_first(self) -> None:
        markers = ["rbcL", "matK"]
        candidates, rejected = candidates_from_records(
            [
                self.record("A1", "Hymenasplenium excisum", "voucherA", {"rbcL", "matK"}, "ACGTNN", 0),
                self.record("B1", "Hymenasplenium excisum", "voucherB", {"rbcL"}, "ACGTAC", 3),
                self.record("C1", "Hymenasplenium unilaterale", "voucherC", {"rbcL"}, "ACGTAC", 0),
            ],
            marker_map(markers),
        )
        selected = select_representatives(candidates, markers)
        self.assertEqual(len(rejected), 0)
        self.assertEqual(len(selected), 2)
        by_taxon = {candidate.taxon: candidate for candidate in selected}
        self.assertEqual(by_taxon["Hymenasplenium excisum"].voucher_or_isolate, "voucherA")

    def test_rejects_uncertain_taxa(self) -> None:
        candidates, rejected = candidates_from_records(
            [
                self.record("A1", "Hymenasplenium sp.", "voucherA", {"rbcL"}, "ACGT"),
                self.record("B1", "Hymenasplenium excisum", "voucherB", {"rbcL"}, "ACGT"),
            ],
            marker_map(["rbcL"]),
        )
        self.assertEqual(len(candidates), 1)
        self.assertEqual(len(rejected), 1)
        self.assertEqual(rejected[0]["reason"], "uncertain_or_incomplete_taxon_name")

    def test_selects_by_extracted_marker_sequence_quality(self) -> None:
        markers = ["rbcL", "matK"]
        records = [
            ExtractedMarkerRecord(
                accession="A1",
                organism="Hymenasplenium excisum",
                marker="rbcL",
                sequence="ACGTNN",
                extraction_method="annotated_feature",
                voucher_or_isolate="voucherA",
            ),
            ExtractedMarkerRecord(
                accession="A2",
                organism="Hymenasplenium excisum",
                marker="matK",
                sequence="ACGTNN",
                extraction_method="annotated_feature",
                voucher_or_isolate="voucherA",
            ),
            ExtractedMarkerRecord(
                accession="B1",
                organism="Hymenasplenium excisum",
                marker="rbcL",
                sequence="ACGTAC",
                extraction_method="whole_record_short",
                voucher_or_isolate="voucherB",
            ),
            ExtractedMarkerRecord(
                accession="B2",
                organism="Hymenasplenium excisum",
                marker="matK",
                sequence="ACGTAC",
                extraction_method="whole_record_short",
                voucher_or_isolate="voucherB",
            ),
        ]
        candidates, rejected = candidates_from_extracted_records(records, marker_map(markers))
        selected = select_representatives(candidates, markers)
        self.assertEqual(rejected, [])
        self.assertEqual(selected[0].voucher_or_isolate, "voucherB")

    def test_writes_minimal_accession_outputs(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            out = Path(tmp)
            result = build_accession_table_from_records(
                records=[
                    self.record("A1", "Hymenasplenium excisum", "voucherA", {"rbcL"}, "ACGT"),
                    self.record("B1", "Hymenasplenium sp.", "voucherB", {"rbcL"}, "ACGT"),
                ],
                output_dir=out,
                markers=["rbcL"],
            )
            rows, fieldnames = read_csv(out / "accession_selected.csv")
            all_rows, _ = read_csv(out / "accession_all.csv")
            self.assertEqual(result.records_written, 1)
            self.assertIn("sample_id", fieldnames)
            self.assertEqual(rows[0]["rbcL"], "A1")
            self.assertEqual({row["selection_status"] for row in all_rows}, {"selected", "rejected"})
            self.assertEqual(next(row for row in all_rows if row["selection_status"] == "rejected")["source_record_ids"], "B1")
            self.assertTrue((out / "fasta_by_marker" / "rbcL.fasta").exists())
            self.assertTrue((out / "summly.txt").exists())


if __name__ == "__main__":
    unittest.main()
