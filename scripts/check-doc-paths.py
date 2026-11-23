#!/usr/bin/env python3
"""Check repository path references into the documentation tree.

The documentation root is ``docs/``.  A reference to ``doc/...`` is a leftover
from the historical ``doc/`` root and must not survive, because a reader or a
tool that follows it silently hits a missing path.  A reference to ``docs/...``
must resolve to an existing file or directory.

Scope and exclusions follow the document retention rules:

* ``docs/reports/`` holds published reports that are kept as historical
  records; their internal path mentions describe the tree at publication time.
* ``docs/plans/`` holds dated plans and handoffs; some describe the upstream
  ``doc/`` layout on purpose.
* ``https://.../blob/.../doc/...`` URLs point at another branch, not this
  tree, and are skipped by the reference matcher.
* ``docs/build/`` is generated output and never committed.

Refs that are intentionally exempt must be declared in
``scripts/doc-path-allowlist.txt`` with a reason, so an exemption is visible in
review instead of being implicit in the code.
"""

from __future__ import annotations

import re
import subprocess
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
ALLOWLIST_PATH = REPO_ROOT / "scripts" / "doc-path-allowlist.txt"

TEXT_SUFFIXES = {
    ".md",
    ".rst",
    ".py",
    ".json",
    ".yml",
    ".yaml",
    ".h",
    ".hpp",
    ".c",
    ".cc",
    ".cpp",
    ".txt",
    ".sh",
    ".bat",
    ".cmake",
}
EXTRA_NAMES = {"Makefile", "make.bat", "Dockerfile"}

# A reference must not be glued to a preceding path separator, word character
# or colon; that keeps URL bodies such as ".../blob/hexin/doc/plans/x.md" out of
# scope while still matching "See doc/source/foo.rst".
REF_PATTERN = re.compile(r"(?<![\w./:-])(?P<root>docs?)/(?P<tail>[A-Za-z0-9_][A-Za-z0-9_./-]*)")

# A shell `cd <dir>` names a directory without a trailing path segment, so
# REF_PATTERN alone misses the exact defect this check was written for:
# `cd doc` left behind after the tree moved to `docs/`.  Match the command form
# separately rather than loosening REF_PATTERN, which would then flag prose.
CD_PATTERN = re.compile(r"\bcd\s+(?P<dir>[^\s;&|()<>]+)")

# Prose also names the root on its own ("the repository's `doc/` directory").
# This is reported as a warning rather than a failure: recording a rename
# legitimately quotes the old name, and a checker cannot tell a stale policy
# statement from a deliberate historical mention.  Only the trailing-slash form
# is matched, because the bare word `doc` appears inside unrelated words and
# identifiers.
BARE_ROOT_PATTERN = re.compile(r"(?<![\w./:-])doc/(?![A-Za-z0-9_])")

TRAILING_PUNCTUATION = ".,;:)]}`'\"*"


def repository_files() -> list[Path]:
    try:
        listed = subprocess.run(
            ["git", "-C", str(REPO_ROOT), "ls-files", "-z"],
            check=True,
            capture_output=True,
            text=True,
        ).stdout
    except (OSError, subprocess.CalledProcessError):
        return [
            path
            for path in REPO_ROOT.rglob("*")
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
    if relative.startswith("docs/reports/"):
        return False
    if relative.startswith("docs/plans/"):
        return False
    # This file quotes the old form in its comments and pattern definitions;
    # checking it would report itself.  The allowlist names exempt references
    # by definition, so it is not a document either.
    if relative in (
        "scripts/check-doc-paths.py",
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
    """Yield ``(line_number, reference, is_prose_only)`` for each reference."""
    for number, line in enumerate(text.splitlines(), start=1):
        seen: set[str] = set()
        for match in REF_PATTERN.finditer(line):
            reference = (match.group("root") + "/" + match.group("tail")).rstrip(
                TRAILING_PUNCTUATION
            )
            if reference.endswith("/"):
                keep = reference.rstrip("/")
                reference = keep
            if len(reference.split("/")) < 2:
                continue
            seen.add(reference)
            yield number, reference, False
        for match in CD_PATTERN.finditer(line):
            directory = match.group("dir").strip(TRAILING_PUNCTUATION)
            if directory in ("doc", "./doc") and "doc/" not in seen:
                yield number, "doc/", False
        if BARE_ROOT_PATTERN.search(line) and "doc/" not in seen:
            yield number, "doc/", True


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
            if reference == "docs/build" or reference.startswith("docs/build/"):
                continue
            if reference.startswith("doc/"):
                (prose if prose_only else legacy).append(entry)
            elif not (REPO_ROOT / reference).exists():
                broken.append(entry)

    for label, items in (("legacy doc/", legacy), ("missing docs/", broken)):
        if items:
            print(f"{label} references ({len(items)}):")
            for item in items:
                print(f"  {item}")
            print()

    if prose:
        print(
            f"bare doc/ mentions ({len(prose)}) - review, not a failure. "
            "A historical mention is fine; a policy statement describing the "
            "current root must say docs/:"
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
    print("OK: documentation path references are consistent with docs/.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
