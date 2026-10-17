# OneFLOW Case Workbench

The Tkinter workbench is a lightweight desktop front end for launching the
existing OneFLOW executable and running case directories in a deliberate order.
It uses only the Python standard library.

## Start the interface

From this directory, run:

```powershell
py -3 tkgui.py
```

On Linux, use:

```bash
python3 tkgui.py
```

Select the built `OneFLOW.exe` (or the platform's OneFLOW executable), keep
mode `r` for a normal run, and add case directories in dependency order. Use
**Run all** to execute the entire list, or select rows and use **Run selected**.
Each case is launched as its own process:

```text
OneFLOW <mode> <absolute-case-directory>
```

The current command-line execution path selects one case per process. The
workbench therefore runs multiple cases sequentially rather than passing all
case paths to one process and assuming the solver executes them all. A nonzero
exit code stops the sequence by default, so dependent cases are not run after a
failed prerequisite.

Run plans can be saved to or loaded from versioned JSON files. The list order
is significant: move cases up or down to make grid-generation cases precede
the solver cases that consume their output.

## Current scope

- Select the executable and execution mode.
- Maintain an ordered list of case directories.
- Run all cases or only selected cases sequentially.
- Stream stdout and stderr into the log panel.
- Stop the active process and optionally stop after the first failed case.
- Save and reload a run plan.

This interface intentionally does not edit solver input files or infer
dependencies automatically. Those can be added as separate layers once the
case metadata and dependency contract are stable.

## Tests

Run the standard-library unit tests from this directory:

```powershell
py -3 -m unittest discover -s . -p "test_*.py"
```
