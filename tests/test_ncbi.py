from __future__ import annotations

import unittest
from unittest.mock import patch

from cladistica.config import marker_map
from cladistica.ncbi import search_ids, search_term_ids


class FakeHandle:
    def __init__(self, payload: dict[str, object]) -> None:
        self.payload = payload

    def __enter__(self) -> FakeHandle:
        return self

    def __exit__(self, exc_type: object, exc: object, traceback: object) -> None:
        return None


class FakeEntrez:
    values = [str(index) for index in range(650)]

    @classmethod
    def esearch(
        cls,
        *,
        db: str,
        term: str,
        retstart: int,
        retmax: int,
        sort: str,
    ) -> FakeHandle:
        page = cls.values[retstart : retstart + retmax]
        return FakeHandle({"Count": str(len(cls.values)), "IdList": page})

    @staticmethod
    def read(handle: FakeHandle) -> dict[str, object]:
        return handle.payload


class NcbiSearchTests(unittest.TestCase):
    def test_search_term_ids_pages_beyond_first_500_records(self) -> None:
        with patch("cladistica.ncbi.require_biopython", return_value=(FakeEntrez, object())):
            ids = search_term_ids("test", limit=None, delay_seconds=0)
        self.assertEqual(len(ids), 650)
        self.assertEqual(ids[-1], "649")

    def test_unlimited_search_deduplicates_overlapping_marker_queries(self) -> None:
        with patch("cladistica.ncbi.require_biopython", return_value=(FakeEntrez, object())):
            ids = search_ids("Testgenus", marker_map(["rbcL"]), retmax=0, delay_seconds=0)
        self.assertEqual(len(ids), 650)

    def test_limited_search_respects_global_retmax(self) -> None:
        with patch("cladistica.ncbi.require_biopython", return_value=(FakeEntrez, object())):
            ids = search_ids("Testgenus", marker_map(["rbcL"]), retmax=10, delay_seconds=0)
        self.assertEqual(len(ids), 10)


if __name__ == "__main__":
    unittest.main()
