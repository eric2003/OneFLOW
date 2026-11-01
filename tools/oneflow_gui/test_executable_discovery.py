"""Platform-specific executable discovery tests."""

from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from tools.oneflow_gui.executable_discovery import executable_candidates


class ExecutableDiscoveryTest(unittest.TestCase):
    def test_posix_candidates_use_unsuffixed_build_output_and_path_fallback(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory) / "repo"
            executable = Path(directory) / "bin" / "OneFLOW"
            executable.parent.mkdir(parents=True)
            executable.write_text("", encoding="utf-8")
            executable.chmod(0o755)

            lookup_calls: list[tuple[str, str | None]] = []

            def find_on_path(name: str, search_path: str | None) -> str:
                lookup_calls.append((name, search_path))
                return str(executable)

            candidates = executable_candidates(
                root,
                {"PATH": str(executable.parent)},
                platform_name="posix",
                path_lookup=find_on_path,
            )

            self.assertEqual(candidates[0], (root / "build" / "bin" / "OneFLOW").resolve())
            self.assertEqual(candidates[1], (root / "build" / "bin" / "Debug" / "OneFLOW").resolve())
            self.assertEqual(candidates[2], (root / "build" / "bin" / "Release" / "OneFLOW").resolve())
            self.assertIn(executable.resolve(), candidates)
            self.assertEqual(lookup_calls, [("OneFLOW", str(executable.parent))])
            self.assertFalse(any(path.name.endswith(".exe") for path in candidates))

    def test_windows_candidates_use_configuration_specific_executable(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            candidates = executable_candidates(
                root,
                {"ProgramFiles": str(root / "Program Files")},
                platform_name="nt",
            )

            self.assertEqual(
                candidates[:2],
                [
                    (root / "build" / "bin" / "Debug" / "OneFLOW.exe").resolve(),
                    (root / "build" / "bin" / "Release" / "OneFLOW.exe").resolve(),
                ],
            )
            self.assertFalse(any(path.name == "OneFLOW" for path in candidates))


if __name__ == "__main__":
    unittest.main()
