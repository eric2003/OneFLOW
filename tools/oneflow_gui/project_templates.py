"""Project templates and safe project creation outside the source checkout."""

from __future__ import annotations

import json
import re
import shutil
from dataclasses import dataclass
from pathlib import Path

from .project_repository import ProjectDefinition, open_project


@dataclass(frozen=True, slots=True)
class ProjectTemplate:
    """A user-facing project template backed by repository content."""

    template_id: str
    title: str
    description: str
    source_relative_path: str | None


def available_project_templates() -> tuple[ProjectTemplate, ...]:
    """Return the initial templates supported by the current project structure."""
    return (
        ProjectTemplate(
            "blank",
            "Blank project",
            "Create an empty project directory with OneFLOW project metadata.",
            None,
        ),
        ProjectTemplate(
            "mesh_workflow",
            "Mesh conversion",
            "Start from the PLOT3D-to-OneFLOW mesh conversion case.",
            "grid/mesh_workflow",
        ),
        ProjectTemplate(
            "minimal_rectangle",
            "Mesh generation",
            "Start from the minimal rectangle mesh generation case.",
            "grid/minimal_rectangle",
        ),
        ProjectTemplate(
            "rae2822_multicase",
            "RAE2822 grid and CFD workflow",
            "Create a grid conversion case and a dependent Roe-SA flow case.",
            None,
        ),
    )


def _copy_tree(source: Path, destination: Path) -> None:
    if not source.is_dir():
        raise FileNotFoundError(f"Project template directory not found: {source}")
    shutil.copytree(
        source,
        destination,
        ignore=shutil.ignore_patterns(".git", "__pycache__"),
    )


def _create_rae2822_workflow(repository: Path, destination: Path, name: str) -> None:
    """Create a workspace for the existing grid-to-flow dependency graph."""
    _copy_tree(repository / "grid" / "mesh_workflow", destination / "rae2822grid")
    _copy_tree(repository / "test" / "cfd" / "rae2822_roe_sa", destination / "rae2822_roe_sa")

    example = repository / "tools" / "oneflow_gui" / "examples" / "rae2822_multicase.json"
    if not example.is_file():
        raise FileNotFoundError(f"Workflow template file not found: {example}")
    shutil.copy2(example, destination / "workflow.json")
    (destination / ".oneflow-project.json").write_text(
        json.dumps(
            {
                "format": "oneflow.project",
                "version": 1,
                "name": name,
                "template": "rae2822_multicase",
            },
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )
    (destination / "README.md").write_text(
        f"# {name}\n\n"
        "This workspace contains two linked OneFLOW cases:\n\n"
        "1. rae2822grid converts the RAE2822 PLOT3D mesh into "
        "results/rae2822_2d.ofl.\n"
        "2. rae2822_roe_sa consumes that generated grid and runs the Roe-SA "
        "flow case.\n\n"
        "The dependency and artifact paths are defined in workflow.json. "
        "Run the workflow from this workspace so generated results stay here.\n",
        encoding="utf-8",
    )


def create_project(
    repository_root: str | Path,
    template_id: str,
    project_name: str,
    parent_directory: str | Path,
) -> ProjectDefinition:
    """Create a project from a template without writing into the source checkout."""
    repository = Path(repository_root).expanduser().resolve()
    parent = Path(parent_directory).expanduser().resolve()
    name = project_name.strip()
    if not name or name in {".", ".."} or not re.fullmatch(r"[^<>:\"/\\|?*]+", name):
        raise ValueError("Project name must be a valid single directory name.")

    template = next(
        (item for item in available_project_templates() if item.template_id == template_id),
        None,
    )
    if template is None:
        raise ValueError(f"Unknown project template: {template_id}")

    try:
        parent.relative_to(repository)
    except ValueError:
        pass
    else:
        raise ValueError("New projects must be created outside the OneFLOW source repository.")

    destination = (parent / name).resolve()
    if destination == repository or repository in destination.parents:
        raise ValueError("New projects must be created outside the OneFLOW source repository.")
    parent.mkdir(parents=True, exist_ok=True)
    if destination.exists():
        raise FileExistsError(f"Project directory already exists: {destination}")

    if template.template_id == "blank":
        destination.mkdir()
        metadata = {
            "format": "oneflow.project",
            "version": 1,
            "name": name,
            "template": template.template_id,
        }
        (destination / ".oneflow-project.json").write_text(
            json.dumps(metadata, indent=2) + "\n", encoding="utf-8"
        )
        (destination / "README.md").write_text(
            f"# {name}\n\nOneFLOW project created from the Blank project template.\n",
            encoding="utf-8",
        )
    elif template.template_id == "rae2822_multicase":
        destination.mkdir()
        _create_rae2822_workflow(repository, destination, name)
    else:
        source = (repository / str(template.source_relative_path)).resolve()
        _copy_tree(source, destination)

    return open_project(destination)
