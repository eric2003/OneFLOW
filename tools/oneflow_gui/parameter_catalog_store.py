"""SQLite index and immutable snapshots for parsed OneFLOW parameters."""
from __future__ import annotations

import hashlib
import json
import sqlite3
from datetime import datetime, timezone
from pathlib import Path

from .parameter_document import ParameterDocument, parameter_specs_for_file


_SCHEMA = """
CREATE TABLE IF NOT EXISTS projects (
    id INTEGER PRIMARY KEY,
    root TEXT NOT NULL UNIQUE,
    created_at TEXT NOT NULL,
    updated_at TEXT NOT NULL
);
CREATE TABLE IF NOT EXISTS parameter_files (
    id INTEGER PRIMARY KEY,
    project_id INTEGER NOT NULL REFERENCES projects(id) ON DELETE CASCADE,
    relative_path TEXT NOT NULL,
    encoding TEXT NOT NULL,
    sha256 TEXT NOT NULL,
    indexed_at TEXT NOT NULL,
    UNIQUE(project_id, relative_path)
);
CREATE TABLE IF NOT EXISTS parameter_values (
    file_id INTEGER NOT NULL REFERENCES parameter_files(id) ON DELETE CASCADE,
    ordinal INTEGER NOT NULL,
    name TEXT NOT NULL,
    declared_type TEXT NOT NULL,
    raw_value TEXT NOT NULL,
    parsed_json TEXT NOT NULL,
    PRIMARY KEY(file_id, ordinal)
);
CREATE INDEX IF NOT EXISTS idx_parameter_values_name ON parameter_values(name);
CREATE TABLE IF NOT EXISTS parameter_snapshots (
    id INTEGER PRIMARY KEY,
    file_id INTEGER NOT NULL REFERENCES parameter_files(id) ON DELETE CASCADE,
    sha256 TEXT NOT NULL,
    encoding TEXT NOT NULL,
    captured_at TEXT NOT NULL,
    UNIQUE(file_id, sha256)
);
CREATE TABLE IF NOT EXISTS parameter_snapshot_values (
    snapshot_id INTEGER NOT NULL REFERENCES parameter_snapshots(id) ON DELETE CASCADE,
    ordinal INTEGER NOT NULL,
    name TEXT NOT NULL,
    declared_type TEXT NOT NULL,
    raw_value TEXT NOT NULL,
    parsed_json TEXT NOT NULL,
    PRIMARY KEY(snapshot_id, ordinal)
);
CREATE TABLE IF NOT EXISTS parameter_definitions (
    file_name TEXT NOT NULL,
    parameter_key TEXT NOT NULL,
    label TEXT NOT NULL,
    description TEXT NOT NULL,
    value_kind TEXT NOT NULL,
    group_name TEXT NOT NULL,
    choices_json TEXT NOT NULL,
    minimum REAL,
    maximum REAL,
    decimals INTEGER NOT NULL,
    PRIMARY KEY(file_name, parameter_key)
);
"""


def _now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


class ParameterCatalogStore:
    """Persist searchable parameter snapshots without replacing source files."""

    def __init__(self, database_path: str | Path):
        self.database_path = Path(database_path).expanduser().resolve()

    def _connect(self) -> sqlite3.Connection:
        self.database_path.parent.mkdir(parents=True, exist_ok=True)
        connection = sqlite3.connect(self.database_path, timeout=10.0)
        connection.execute("PRAGMA foreign_keys = ON")
        connection.executescript(_SCHEMA)
        return connection

    def sync_document(self, project_root: str | Path, document: ParameterDocument) -> bool:
        """Index one file; return True only when a new content snapshot is recorded."""
        root = Path(project_root).expanduser().resolve()
        try:
            relative = document.path.resolve().relative_to(root).as_posix()
        except ValueError as error:
            raise ValueError("Parameter files must be inside their project root.") from error
        digest = hashlib.sha256(document.original_bytes).hexdigest()
        now = _now()
        with self._connect() as db:
            row = db.execute("SELECT id FROM projects WHERE root = ?", (str(root),)).fetchone()
            if row is None:
                project_id = db.execute(
                    "INSERT INTO projects(root, created_at, updated_at) VALUES (?, ?, ?)",
                    (str(root), now, now),
                ).lastrowid
            else:
                project_id = int(row[0])
                db.execute("UPDATE projects SET updated_at = ? WHERE id = ?", (now, project_id))
            row = db.execute(
                "SELECT id, sha256 FROM parameter_files WHERE project_id = ? AND relative_path = ?",
                (project_id, relative),
            ).fetchone()
            if row is None:
                file_id = db.execute(
                    "INSERT INTO parameter_files(project_id, relative_path, encoding, sha256, indexed_at) "
                    "VALUES (?, ?, ?, ?, ?)",
                    (project_id, relative, document.encoding, digest, now),
                ).lastrowid
                previous_digest = None
            else:
                file_id, previous_digest = int(row[0]), str(row[1])
                db.execute(
                    "UPDATE parameter_files SET encoding = ?, sha256 = ?, indexed_at = ? WHERE id = ?",
                    (document.encoding, digest, now, file_id),
                )
            if previous_digest == digest:
                return False
            db.execute("DELETE FROM parameter_values WHERE file_id = ?", (file_id,))
            snapshot_row = db.execute(
                "SELECT id FROM parameter_snapshots WHERE file_id = ? AND sha256 = ?",
                (file_id, digest),
            ).fetchone()
            is_new_snapshot = snapshot_row is None
            if is_new_snapshot:
                db.execute(
                    "INSERT INTO parameter_snapshots(file_id, sha256, encoding, captured_at) "
                    "VALUES (?, ?, ?, ?)", (file_id, digest, document.encoding, now)
                )
                snapshot_row = db.execute(
                    "SELECT id FROM parameter_snapshots WHERE file_id = ? AND sha256 = ?",
                    (file_id, digest),
                ).fetchone()
            snapshot_id = int(snapshot_row[0])
            for ordinal, item in enumerate(document.declarations):
                values = (
                    ordinal, item.name, item.declared_type, item.raw_value,
                    json.dumps(item.value, ensure_ascii=True, allow_nan=False),
                )
                db.execute(
                    "INSERT INTO parameter_values(file_id, ordinal, name, declared_type, raw_value, parsed_json) "
                    "VALUES (?, ?, ?, ?, ?, ?)", (file_id, *values)
                )
                if is_new_snapshot:
                    db.execute(
                        "INSERT INTO parameter_snapshot_values(snapshot_id, ordinal, name, declared_type, raw_value, parsed_json) "
                        "VALUES (?, ?, ?, ?, ?, ?)", (snapshot_id, *values)
                    )
            for spec in parameter_specs_for_file(document.path):
                choices = [{"label": item.label, "value": item.value} for item in spec.choices]
                db.execute(
                    "INSERT INTO parameter_definitions VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?) "
                    "ON CONFLICT(file_name, parameter_key) DO UPDATE SET "
                    "label=excluded.label, description=excluded.description, value_kind=excluded.value_kind, "
                    "group_name=excluded.group_name, choices_json=excluded.choices_json, "
                    "minimum=excluded.minimum, maximum=excluded.maximum, decimals=excluded.decimals",
                    (document.path.name.casefold(), spec.key, spec.label, spec.description,
                     spec.value_kind, spec.group, json.dumps(choices, ensure_ascii=True),
                     spec.minimum, spec.maximum, spec.decimals),
                )
        return True

    def list_project_files(self, project_root: str | Path) -> list[dict[str, object]]:
        root = Path(project_root).expanduser().resolve()
        with self._connect() as db:
            rows = db.execute(
                "SELECT f.relative_path, f.encoding, f.sha256, f.indexed_at, "
                "(SELECT COUNT(*) FROM parameter_values v WHERE v.file_id=f.id) "
                "FROM parameter_files f JOIN projects p ON p.id=f.project_id "
                "WHERE p.root=? ORDER BY f.relative_path", (str(root),)
            ).fetchall()
        return [{"relative_path": r[0], "encoding": r[1], "sha256": r[2],
                 "indexed_at": r[3], "parameter_count": r[4]} for r in rows]

    def list_values(self, project_root: str | Path, relative_path: str) -> list[dict[str, object]]:
        root = Path(project_root).expanduser().resolve()
        with self._connect() as db:
            rows = db.execute(
                "SELECT v.name, v.declared_type, v.raw_value, v.parsed_json, v.ordinal "
                "FROM parameter_values v JOIN parameter_files f ON f.id=v.file_id "
                "JOIN projects p ON p.id=f.project_id WHERE p.root=? AND f.relative_path=? "
                "ORDER BY v.ordinal", (str(root), relative_path.replace("\\", "/"))
            ).fetchall()
        return [{"name": r[0], "declared_type": r[1], "raw_value": r[2],
                 "value": json.loads(r[3]), "ordinal": r[4]} for r in rows]

    def list_history(self, project_root: str | Path, relative_path: str) -> list[dict[str, str]]:
        root = Path(project_root).expanduser().resolve()
        with self._connect() as db:
            rows = db.execute(
                "SELECT s.sha256, s.encoding, s.captured_at FROM parameter_snapshots s "
                "JOIN parameter_files f ON f.id=s.file_id JOIN projects p ON p.id=f.project_id "
                "WHERE p.root=? AND f.relative_path=? ORDER BY s.id DESC",
                (str(root), relative_path.replace("\\", "/"))
            ).fetchall()
        return [{"sha256": r[0], "encoding": r[1], "captured_at": r[2]} for r in rows]
