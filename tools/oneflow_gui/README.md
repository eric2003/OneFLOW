# OneFLOW Desktop Case Runner

The desktop application provides a small, extensible GUI for selecting cases
from the shared catalog and running them through the normal OneFLOW executable.

## Launch

From the repository root, using Python 3.10 or newer:

```powershell
python -m pip install -r tools/requirements-gui.txt
python -m tools.oneflow_gui
```

Choose the built `OneFLOW.exe`, choose an external run workspace, select one or
more cases, and click **Run selected cases**. The queue runs cases sequentially
and displays process output and per-case status. **Stop queue** terminates the
active process and prevents pending cases from launching.

Each case is copied to a unique directory under the selected workspace before
execution. Generated files therefore do not overwrite the checked-in case
inputs. The application invokes the same command shape used by the existing
workflow runner: `OneFLOW.exe r <case-directory>`.

## Structure

- `models.py`: immutable case and result models plus the run-state enum.
- `case_repository.py`: catalog loading and shared validation.
- `runner.py`: workspace isolation, process launch, output streaming, and
  cancellation; it has no Qt dependency.
- `main_window.py`: Qt widgets and the sequential queue worker.
- `test_services.py`: standard-library tests for the non-UI services.

The GUI is an optional tool and is not a runtime dependency of OneFLOW. Core
case discovery and process execution remain separable from the GUI framework.

## Service tests

From the repository root:

```powershell
python -m unittest tools.oneflow_gui.test_services
python tools/test_mesh_cases.py
```

The first version reports process exit status and logs. It does not yet claim
that a successful process exit proves mesh geometry is valid; result validators
will be integrated as a separate workflow stage.
