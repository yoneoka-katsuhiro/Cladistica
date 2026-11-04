"""XLSX readers must release handles even when a workbook stays referenced."""
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from openpyxl import Workbook, load_workbook

from cladistica.io import read_xlsx


class ExcelResourceTests(unittest.TestCase):
    def check_file_released(self, *, empty=False, broken=False):
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "input.xlsx"
            workbook = Workbook()
            if not empty:
                workbook.active.append(["sample_id", "rbcL"])
                workbook.active.append(["specimen", "ABC.1"])
            workbook.save(path)
            references = []

            def open_and_retain(source, **kwargs):
                result = load_workbook(source, **kwargs)
                references.append((source, result))
                if broken:
                    def failing_rows(**kwargs):
                        yield ("sample_id",)
                        raise ValueError("invalid sheet content")
                    result.active.iter_rows = failing_rows
                return result

            with patch("openpyxl.load_workbook", side_effect=open_and_retain):
                if broken:
                    with self.assertRaisesRegex(ValueError, "invalid sheet content"):
                        read_xlsx(path)
                else:
                    rows, columns = read_xlsx(path)
                    self.assertEqual(rows, [] if empty else [{"sample_id": "specimen", "rbcL": "ABC.1"}])
                    self.assertEqual(columns, [] if empty else ["sample_id", "rbcL"])
            # Retain the Workbook so GC cannot accidentally mask missing close().
            self.assertIsNone(references[0][1]._archive.fp)
            self.assertTrue(references[0][0].closed)
            moved = path.with_name("moved.xlsx")
            path.rename(moved)  # Windows rejects this if any XLSX handle leaked.
            moved.unlink()

    def test_populated_workbook_releases_file(self):
        self.check_file_released()

    def test_empty_workbook_releases_file(self):
        self.check_file_released(empty=True)

    def test_iteration_failure_releases_file(self):
        self.check_file_released(broken=True)


if __name__ == "__main__":
    unittest.main()
