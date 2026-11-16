"""Tests for user-owned non-STL geometry asset registration."""

from __future__ import annotations

import hashlib
import json
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

from tools.oneflow_gui.geometry_asset_store import import_geometry_asset, load_user_geometry_assets


class GeometryAssetStoreTest(unittest.TestCase):
    def test_step_import_is_copied_and_registered_by_verified_checksum(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source = root / "part.step"
            source.write_text("ISO-10303-21;\nEND-ISO-10303-21;\n", encoding="ascii")
            user_root = root / "user-data"
            with patch(
                "tools.oneflow_gui.geometry_asset_store.UserDataPaths.discover",
                return_value=SimpleNamespace(root=user_root),
            ):
                asset = import_geometry_asset(source)
                copied = user_root / asset["representations"][0]["path"]
                self.assertTrue(copied.is_file())
                digest = hashlib.sha256(source.read_bytes()).hexdigest()
                self.assertIsNone(asset["dimension"])
                self.assertEqual(asset["provenance"]["source_sha256"], digest)
                self.assertEqual(asset["representations"][0]["sha256"], digest)
                self.assertEqual(asset["representations"][0]["role"], "cad_brep")
                self.assertEqual(asset["metadata"]["registration_status"], "registered_not_parsed")
                self.assertEqual(load_user_geometry_assets(), [asset])

    def test_implicit_json_import_is_not_claimed_as_parsed_geometry(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source = root / "field.implicit.json"
            source.write_text(json.dumps({"field": "signed_distance"}), encoding="utf-8")
            with patch(
                "tools.oneflow_gui.geometry_asset_store.UserDataPaths.discover",
                return_value=SimpleNamespace(root=root / "user"),
            ):
                asset = import_geometry_asset(source)
            self.assertEqual(asset["definition"]["kind"], "implicit")
            self.assertEqual(asset["representations"][0]["role"], "implicit_field")
            self.assertEqual(asset["metadata"]["registration_status"], "registered_not_parsed")

    def test_unknown_extension_is_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            source = Path(directory) / "model.unknown"
            source.write_text("data", encoding="ascii")
            with self.assertRaisesRegex(ValueError, "Unsupported geometry import format"):
                import_geometry_asset(source)


if __name__ == "__main__":
    unittest.main()
