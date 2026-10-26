"""Convert portable .ofw definitions into existing workflow drafts."""

from __future__ import annotations

from pathlib import Path

from .workflow_repository import WorkflowDefinition
from .workspace_repository import WorkspaceDefinition


def resolve_workspace_case_directory(
    case_path: str, workspace_path: str | Path
) -> Path:
    """Resolve a case path relative to the containing .ofw file."""
    path = Path(case_path).expanduser()
    base_directory = Path(workspace_path).expanduser().resolve().parent
    return path.resolve() if path.is_absolute() else (base_directory / path).resolve()


def store_workspace_case_path(case_path: str, workspace_path: str | Path) -> str:
    """Keep paths inside a workspace relative for portable serialization."""
    workspace_file = Path(workspace_path).expanduser().resolve()
    resolved = resolve_workspace_case_directory(case_path, workspace_file)
    try:
        return resolved.relative_to(workspace_file.parent).as_posix()
    except ValueError:
        return str(resolved)


def build_workflow_definition_from_workspace(
    workspace: WorkspaceDefinition,
    workspace_path: str | Path,
) -> WorkflowDefinition:
    """Translate a saved workspace into the existing workflow editor format.

    Single, sequential, and independent workspaces use OneFLOW's existing
    single-process multi-directory contract. Dependency workspaces preserve
    their edges, but .ofw currently has no artifact binding schema.
    """
    source_path = Path(workspace_path).expanduser().resolve()
    base_directory = source_path.parent

    if workspace.execution_mode == "dependency_workflow":
        cases = [
            {
                "id": case.case_id,
                "directory": case.path,
                "depends_on": list(case.depends_on),
                "inputs": [item.to_dict() for item in case.inputs],
                "outputs": [item.to_dict() for item in case.outputs],
            }
            for case in workspace.cases
        ]
        step = {
            "type": "case_graph",
            "base_directory": str(base_directory),
            "cases": cases,
        }
    elif workspace.execution_mode in {
        "single_case", "sequential_cases", "independent_cases"
    }:
        step = {
            "type": "multi_directory",
            "directories": [str(resolve_workspace_case_directory(case.path, source_path)) for case in workspace.cases],
        }
    else:
        raise ValueError(f"Unsupported workspace execution mode: {workspace.execution_mode!r}.")

    return WorkflowDefinition(
        workflow_id=f"workspace-{workspace.workspace_id}",
        name=workspace.name,
        description=workspace.description,
        steps=(step,),
    )
