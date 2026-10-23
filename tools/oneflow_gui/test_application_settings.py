"""Tests for persistent OneFLOW application settings."""

from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

from tools.oneflow_gui.application_settings import (
    ApplicationSettings,
    DEFAULT_REPOSITORY_ROOT,
    resolve_repository_root,
)


class ApplicationSettingsTest(unittest.TestCase):
    def test_missing_settings_use_default_repository_root(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            settings = ApplicationSettings.load(Path(directory) / "settings.json")
            self.assertEqual(settings.repository_root, DEFAULT_REPOSITORY_ROOT)

    def test_save_and_load_repository_root(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            settings_file = root / "user-data" / "settings.json"
            repository_root = root / "checkout"
            repository_root.mkdir()

            ApplicationSettings(repository_root).save(settings_file)
            settings = ApplicationSettings.load(settings_file)

            self.assertEqual(settings.repository_root, repository_root.resolve())
            self.assertEqual(
                json.loads(settings_file.read_text(encoding="utf-8")),
                {"repository_root": str(repository_root.resolve())},
            )

    def test_missing_configured_checkout_falls_back_to_launch_checkout(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            launch_root = Path(directory)
            resolved = resolve_repository_root(
                ApplicationSettings(DEFAULT_REPOSITORY_ROOT),
                launch_root,
            )
            self.assertEqual(resolved, launch_root.resolve())

    def test_missing_custom_checkout_falls_back_to_launch_checkout(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            missing = root / "missing-checkout"
            launch_root = root / "launch"
            launch_root.mkdir()
            resolved = resolve_repository_root(ApplicationSettings(missing), launch_root)
            self.assertEqual(resolved, launch_root.resolve())

    def test_existing_configured_checkout_takes_precedence(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            configured = root / "configured"
            launch_root = root / "launch"
            configured.mkdir()
            launch_root.mkdir()

            resolved = resolve_repository_root(
                ApplicationSettings(configured),
                launch_root,
            )
            self.assertEqual(resolved, configured.resolve())


if __name__ == "__main__":
    unittest.main()
