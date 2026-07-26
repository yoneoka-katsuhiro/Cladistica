from __future__ import annotations

import tempfile
import unittest
from datetime import date
from pathlib import Path

from cladistica.workflow import default_run_dir


class WorkflowTests(unittest.TestCase):
    def test_default_run_dir_uses_date_and_increment(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            stamp = date.today().strftime("%Y%m%d")
            first = default_run_dir(root)
            self.assertEqual(first, root / "output" / f"{stamp}_1")
            first.mkdir(parents=True)
            self.assertEqual(default_run_dir(root), root / "output" / f"{stamp}_2")


if __name__ == "__main__":
    unittest.main()
