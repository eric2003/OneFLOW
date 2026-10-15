# OneFLOW Desktop Case Runner

The desktop application provides a small, extensible GUI for selecting cases
from the shared catalog and running them through the normal OneFLOW executable.

## Launch

From the repository root, using Python 3.10 or newer:

```powershell
python -m pip install -r tools/requirements-gui.txt
python -m tools.oneflow_gui
```

The executable field stays as a single compact line. Existing executables are
detected in this order: repository Debug build, repository Release build,
`C:/dev/OneFLOW/bin/OneFLOW.exe`, then conventional Windows installation
folders under `ProgramFiles`, `ProgramFiles(x86)`, and `LOCALAPPDATA`.
The first detected path is selected automatically; use the dropdown to switch
between detected copies, type a path directly, or use **Browse...** for an
unlisted location. No full-disk search is performed.

Choose an external run workspace, select one or more cases, and click
**Run selected cases**. The queue runs cases sequentially
and displays process output and per-case status. **Stop queue** terminates the
active process and prevents pending cases from launching.

Each case is copied to a unique directory under the selected workspace before
execution. Generated files therefore do not overwrite the checked-in case
inputs. The application invokes the same command shape used by the existing
workflow runner: `OneFLOW.exe r <case-directory>`.

## Architecture direction

The long-term product direction is documented in both [Chinese](../../docs/ONEFLOW_WORKBENCH_ARCHITECTURE.zh-CN.md) and [English](../../docs/ONEFLOW_WORKBENCH_ARCHITECTURE.en.md). These memos describe the intended engineering-workbench boundaries and roadmap; they do not imply that all listed capabilities are implemented yet.

## Structure

- `models.py`: immutable case and result models plus the run-state enum.
- `case_repository.py`: catalog loading and shared validation.
- `executable_discovery.py`: ordered lookup of conventional executable paths.
- `runner.py`: workspace isolation, process launch, output streaming, and
  cancellation; it has no Qt dependency.
- `workflow_execution.py`: shell-free construction of one-process,
  multi-directory command arguments, preserving the user-specified order.
- `main_window.py`: Qt widgets and the sequential queue worker.
- `test_services.py`: standard-library tests for the non-UI services.
- `test_workflow_execution.py`: tests for argument order, spaces, and invalid input.
- `workflow_repository.py`: versioned workflow-definition and run-history persistence, independent of Qt.
- `test_workflow_repository.py`: persistence tests for version checks, CRUD behavior, and history retention.

The GUI is an optional tool and is not a runtime dependency of OneFLOW. Core
case discovery and process execution remain separable from the GUI framework.

## Service tests

From the repository root:

```powershell
python -m unittest tools.oneflow_gui.test_services
python -m unittest tools.oneflow_gui.test_workflow_repository
python -m unittest tools.oneflow_gui.test_workflow_execution
python tools/test_mesh_cases.py
```

The first version reports process exit status and logs. It does not yet claim
that a successful process exit proves mesh geometry is valid; result validators
will be integrated as a separate workflow stage.
