"""Tests for SQLite-backed parameter catalog persistence."""
from __future__ import annotations

import sqlite3
import tempfile
import unittest
from pathlib import Path

from tools.oneflow_gui.parameter_document import ParameterDocument
from tools.oneflow_gui.parameter_catalog_store import ParameterCatalogStore


class ParameterCatalogStoreTest(unittest.TestCase):
    def test_indexes_typed_values_and_metadata_without_changing_source(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary) / "case"
            script = root / "script"
            script.mkdir(parents=True)
            path = script / "cfd.txt"
            original = b"int maxSteps = 50; // keep comment\r\nReal maxTime = 0.2;\r\n"
            path.write_bytes(original)
            store = ParameterCatalogStore(Path(temporary) / "db" / "catalog.sqlite3")
            document = ParameterDocument.load(path)
            self.assertTrue(store.sync_document(root, document))
            self.assertFalse(store.sync_document(root, document))
            self.assertEqual(path.read_bytes(), original)
            files = store.list_project_files(root)
            self.assertEqual(len(files), 1)
            self.assertEqual(files[0]["parameter_count"], 2)
            values = store.list_values(root, "script/cfd.txt")
            self.assertEqual(values[0]["name"], "maxSteps")
            self.assertEqual(values[0]["value"], 50)
            self.assertEqual(values[1]["value"], 0.2)

    def test_changed_file_adds_snapshot_and_refreshes_current_values(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary) / "case"
            script = root / "script"
            script.mkdir(parents=True)
            path = script / "cfd.txt"
            path.write_text("int maxSteps = 50;\n", encoding="utf-8")
            store = ParameterCatalogStore(Path(temporary) / "catalog.sqlite3")
            self.assertTrue(store.sync_document(root, ParameterDocument.load(path)))
            path.write_text("int maxSteps = 100;\n", encoding="utf-8")
            self.assertTrue(store.sync_document(root, ParameterDocument.load(path)))
            self.assertEqual(len(store.list_history(root, "script/cfd.txt")), 2)
            self.assertEqual(store.list_values(root, "script/cfd.txt")[0]["value"], 100)

    def test_reverting_to_an_existing_content_hash_reuses_snapshot(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary) / "case"
            script = root / "script"
            script.mkdir(parents=True)
            path = script / "cfd.txt"
            database = Path(temporary) / "catalog.sqlite3"
            store = ParameterCatalogStore(database)
            for value in (50, 100, 50):
                path.write_text(f"int maxSteps = {value};\n", encoding="utf-8")
                self.assertTrue(store.sync_document(root, ParameterDocument.load(path)))
            self.assertEqual(len(store.list_history(root, "script/cfd.txt")), 2)
            self.assertEqual(store.list_values(root, "script/cfd.txt")[0]["value"], 50)

    def test_rejects_file_outside_project_root(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary) / "case"
            root.mkdir()
            outside = Path(temporary) / "outside.txt"
            outside.write_text("int maxSteps = 50;\n", encoding="utf-8")
            store = ParameterCatalogStore(Path(temporary) / "catalog.sqlite3")
            with self.assertRaises(ValueError):
                store.sync_document(root, ParameterDocument.load(outside))

    def test_database_integrity_and_schema_metadata(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary) / "case"
            script = root / "script"
            script.mkdir(parents=True)
            path = script / "cfd.txt"
            path.write_text("int maxSteps = 50;\n", encoding="utf-8")
            database = Path(temporary) / "catalog.sqlite3"
            ParameterCatalogStore(database).sync_document(root, ParameterDocument.load(path))
            with sqlite3.connect(database) as db:
                integrity = db.execute("PRAGMA integrity_check").fetchone()[0]
                definition = db.execute(
                    "SELECT value_kind, label FROM parameter_definitions "
                    "WHERE file_name='cfd.txt' AND parameter_key='maxSteps'"
                ).fetchone()
            self.assertEqual(integrity, "ok")
            self.assertEqual(definition, ("integer", "Maximum iterations"))


if __name__ == "__main__":
    unittest.main()
