#!/usr/bin/env python3
"""Check repository path references into the documentation tree."""
from __future__ import annotations
import re
import subprocess
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
ALLOWLIST_PATH = REPO_ROOT / "scripts" / "doc-path-allowlist.txt"
# Detect the documentation root from the checked-out tree during root migration.
CANONICAL_ROOT = "docs" if (REPO_ROOT / "docs" / "source").is_dir() else "doc"
LEGACY_ROOT = "doc" if CANONICAL_ROOT == "docs" else "docs"

TEXT_SUFFIXES = {
    ".md", ".rst", ".py", ".json", ".yml", ".yaml", ".h", ".hpp",
    ".c", ".cc", ".cpp", ".txt", ".sh", ".bat", ".cmake",
}
EXTRA_NAMES = {"Makefile", "make.bat", "Dockerfile"}
REF_PATTERN = re.compile(
    r"(?<![\w./:-])(?P<root>docs?)/(?P<tail>[A-Za-z0-9_][A-Za-z0-9_./-]*)"
)
CD_PATTERN = re.compile(r"\bcd\s+(?P<dir>[^\s;&|()<>]+)")
BARE_ROOT_PATTERN = re.compile(
    rf"(?<![\w./:-]){LEGACY_ROOT}/(?![A-Za-z0-9_])"
)
TRAILING_PUNCTUATION = ".,;:)]}'\"*" + chr(96)


def repository_files() -> list[Path]:
    try:
        listed = subprocess.run(
            ["git", "-C", str(REPO_ROOT), "ls-files", "-z"],
            check=True, capture_output=True, text=True,
        ).stdout
    except (OSError, subprocess.CalledProcessError):
        return [
            path for path in REPO_ROOT.rglob("*")
            if path.is_file() and ".git" not in path.parts
        ]
    return [REPO_ROOT / entry for entry in listed.split("\0") if entry]


def is_in_scope(path: Path) -> bool:
    try:
        relative = path.relative_to(REPO_ROOT).as_posix()
    except ValueError:
        return False
    if relative.startswith("ThirdParty/"):
        return False
    if "/build/" in relative or relative.startswith("build/"):
        return False
    if any(
        relative.startswith(f"{CANONICAL_ROOT}/{directory}/")
        for directory in ("reports", "plans")
    ):
        return False
    if relative in (
        "scripts/check-doc-paths.py",
        "scripts/test_check_doc_paths.py",
        "scripts/doc-path-allowlist.txt",
    ):
        return False
    return path.suffix in TEXT_SUFFIXES or path.name in EXTRA_NAMES


def load_allowlist() -> dict[str, str]:
    entries: dict[str, str] = {}
    if not ALLOWLIST_PATH.is_file():
        return entries
    for raw in ALLOWLIST_PATH.read_text(encoding="utf-8").splitlines():
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        reference, separator, reason = line.partition(" ")
        if not separator or not reason.strip():
            raise SystemExit(
                f"{ALLOWLIST_PATH.name}: entry needs '<ref> <reason>': {line}"
            )
        entries[reference] = reason.strip()
    return entries


def iter_references(text: str):
    """Yield (line_number, reference, is_prose_only) for each reference."""
    for number, line in enumerate(text.splitlines(), start=1):
        seen: set[str] = set()
        for match in REF_PATTERN.finditer(line):
            reference = (match.group("root") + "/" + match.group("tail")).rstrip(
                TRAILING_PUNCTUATION
            )
            if reference.endswith("/"):
                reference = reference.rstrip("/")
            if len(reference.split("/")) < 2:
                continue
            seen.add(reference)
            yield number, reference, False

        # Avoid reporting a path twice as both a full reference and bare prose.
        has_legacy_reference = any(
            reference.startswith(f"{LEGACY_ROOT}/") for reference in seen
        )
        for match in CD_PATTERN.finditer(line):
            directory = match.group("dir").strip(TRAILING_PUNCTUATION)
            if directory in (LEGACY_ROOT, f"./{LEGACY_ROOT}") and not has_legacy_reference:
                yield number, f"{LEGACY_ROOT}/", False
        if BARE_ROOT_PATTERN.search(line) and not has_legacy_reference:
            yield number, f"{LEGACY_ROOT}/", True


def main() -> int:
    allowlist = load_allowlist()
    legacy: list[str] = []
    broken: list[str] = []
    exempt: list[str] = []
    prose: list[str] = []
    for path in repository_files():
        if not is_in_scope(path):
            continue
        try:
            text = path.read_text(encoding="utf-8")
        except (OSError, UnicodeDecodeError):
            continue
        relative = path.relative_to(REPO_ROOT).as_posix()
        for number, reference, prose_only in iter_references(text):
            entry = f"{relative}:{number} {reference}"
            if reference in allowlist:
                exempt.append(entry)
                continue
            build_output = f"{CANONICAL_ROOT}/build"
            if reference == build_output or reference.startswith(build_output + "/"):
                continue
            if reference.startswith(f"{LEGACY_ROOT}/"):
                (prose if prose_only else legacy).append(entry)
            elif not (REPO_ROOT / reference).exists():
                broken.append(entry)

    for label, items in (
        (f"legacy {LEGACY_ROOT}/", legacy),
        (f"missing {CANONICAL_ROOT}/", broken),
    ):
        if items:
            print(f"{label} references ({len(items)}):")
            for item in items:
                print(f"  {item}")
            print()
    if prose:
        print(
            f"bare {LEGACY_ROOT}/ mentions ({len(prose)}) - review, not a failure. "
            "Historical mentions are valid; current-layout statements should use "
            f"{CANONICAL_ROOT}/:"
        )
        for item in prose:
            print(f"  {item}")
        print()
    if exempt:
        print(f"exempt references ({len(exempt)}):")
        for item in exempt:
            print(f"  {item}")
        print()
    total = len(legacy) + len(broken)
    if total:
        print(f"FAIL: {total} documentation path reference(s) must be fixed.")
        return 1
    print(f"OK: documentation path references are consistent with {CANONICAL_ROOT}/.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
