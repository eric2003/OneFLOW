#!/usr/bin/env python3
"""Tkinter front end for running OneFLOW cases in a controlled sequence.

Start this file with the Python interpreter. It uses only the Python standard
library and does not need to be built by CMake.
"""

from __future__ import annotations

from pathlib import Path
from queue import Empty, Queue
import subprocess
import sys
import tkinter as tk
from tkinter import filedialog, messagebox, ttk

from oneflow_runner import OneFlowRunner, RunPlan, load_run_plan, save_run_plan


REPOSITORY_ROOT = Path(__file__).resolve().parents[3]
if str(REPOSITORY_ROOT) not in sys.path:
    sys.path.insert(0, str(REPOSITORY_ROOT))

from tools.oneflow_gui.executable_discovery import executable_candidates


class OneFlowWorkbench:
    """Small case workbench with an explicit ordered run list and live output."""

    def __init__(self, root: tk.Tk) -> None:
        self.root = root
        self.root.title("OneFLOW Case Workbench")
        self.root.geometry("1100x760")
        self.root.minsize(850, 580)

        self.events: Queue[dict] = Queue()
        self.runner = OneFlowRunner(self.events.put)
        self.case_directories: list[Path] = []
        self.active = False

        candidates = executable_candidates(REPOSITORY_ROOT)
        build_candidates = [
            REPOSITORY_ROOT / "build" / "bin" / "Debug" / "OneFLOW.exe",
            REPOSITORY_ROOT / "build" / "bin" / "Release" / "OneFLOW.exe",
        ]
        executable_options = list(dict.fromkeys(
            str(path) for path in [*candidates, *build_candidates]
        ))
        default_executable = str(candidates[0]) if candidates else ""
        self.executable_var = tk.StringVar(value=default_executable)
        self.executable_options = executable_options
        self.mode_var = tk.StringVar(value="r")
        self.stop_on_error_var = tk.BooleanVar(value=True)
        self.status_var = tk.StringVar(value="Ready")

        self._build_interface()
        self.root.protocol("WM_DELETE_WINDOW", self._on_close)
        self.root.after(100, self._poll_events)

    def _build_interface(self) -> None:
        """Create the interface in separate configuration, cases, and log areas."""
        self.root.columnconfigure(0, weight=1)
        self.root.rowconfigure(2, weight=1)

        executable_frame = ttk.LabelFrame(
            self.root, text="Execution", padding=(12, 10)
        )
        executable_frame.grid(row=0, column=0, sticky="ew", padx=12, pady=(12, 6))
        executable_frame.columnconfigure(1, weight=1)

        ttk.Label(executable_frame, text="OneFLOW executable").grid(
            row=0, column=0, sticky="w", padx=(0, 8), pady=4
        )
        self.executable_combo = ttk.Combobox(
            executable_frame,
            textvariable=self.executable_var,
            values=self.executable_options,
            state="normal",
        )
        self.executable_combo.grid(row=0, column=1, sticky="ew", pady=4)
        ttk.Button(
            executable_frame, text="Browse...", command=self._choose_executable
        ).grid(row=0, column=2, padx=(8, 0), pady=4)

        ttk.Label(executable_frame, text="Mode").grid(
            row=1, column=0, sticky="w", padx=(0, 8), pady=4
        )
        ttk.Entry(executable_frame, textvariable=self.mode_var, width=8).grid(
            row=1, column=1, sticky="w", pady=4
        )
        ttk.Checkbutton(
            executable_frame,
            text="Stop after the first failed case",
            variable=self.stop_on_error_var,
        ).grid(row=1, column=2, sticky="w", padx=(8, 0), pady=4)

        case_frame = ttk.LabelFrame(
            self.root, text="Ordered case list", padding=(12, 10)
        )
        case_frame.grid(row=1, column=0, sticky="ew", padx=12, pady=6)
        case_frame.columnconfigure(0, weight=1)

        self.case_list = tk.Listbox(
            case_frame,
            height=7,
            selectmode=tk.EXTENDED,
            activestyle="underline",
            exportselection=False,
        )
        self.case_list.grid(row=0, column=0, rowspan=2, sticky="nsew")
        case_scrollbar = ttk.Scrollbar(
            case_frame, orient=tk.VERTICAL, command=self.case_list.yview
        )
        case_scrollbar.grid(row=0, column=1, rowspan=2, sticky="ns")
        self.case_list.configure(yscrollcommand=case_scrollbar.set)

        case_buttons = ttk.Frame(case_frame)
        case_buttons.grid(row=0, column=2, sticky="n", padx=(10, 0))
        ttk.Button(case_buttons, text="Add case...", command=self._add_case).grid(
            row=0, column=0, sticky="ew", pady=2
        )
        ttk.Button(case_buttons, text="Remove selected", command=self._remove_cases).grid(
            row=1, column=0, sticky="ew", pady=2
        )
        ttk.Button(case_buttons, text="Move up", command=lambda: self._move_case(-1)).grid(
            row=2, column=0, sticky="ew", pady=2
        )
        ttk.Button(case_buttons, text="Move down", command=lambda: self._move_case(1)).grid(
            row=3, column=0, sticky="ew", pady=2
        )
        ttk.Separator(case_buttons, orient=tk.HORIZONTAL).grid(
            row=4, column=0, sticky="ew", pady=7
        )
        ttk.Button(case_buttons, text="Load plan...", command=self._load_plan).grid(
            row=5, column=0, sticky="ew", pady=2
        )
        ttk.Button(case_buttons, text="Save plan...", command=self._save_plan).grid(
            row=6, column=0, sticky="ew", pady=2
        )

        log_frame = ttk.LabelFrame(self.root, text="Run output", padding=(12, 10))
        log_frame.grid(row=2, column=0, sticky="nsew", padx=12, pady=6)
        log_frame.columnconfigure(0, weight=1)
        log_frame.rowconfigure(0, weight=1)

        self.log_view = tk.Text(
            log_frame,
            wrap=tk.NONE,
            state=tk.DISABLED,
            font=("Consolas", 10),
            undo=False,
        )
        self.log_view.grid(row=0, column=0, sticky="nsew")
        log_y = ttk.Scrollbar(
            log_frame, orient=tk.VERTICAL, command=self.log_view.yview
        )
        log_y.grid(row=0, column=1, sticky="ns")
        log_x = ttk.Scrollbar(
            log_frame, orient=tk.HORIZONTAL, command=self.log_view.xview
        )
        log_x.grid(row=1, column=0, sticky="ew")
        self.log_view.configure(yscrollcommand=log_y.set, xscrollcommand=log_x.set)

        footer = ttk.Frame(self.root, padding=(12, 6, 12, 12))
        footer.grid(row=3, column=0, sticky="ew")
        footer.columnconfigure(0, weight=1)
        ttk.Label(footer, textvariable=self.status_var).grid(
            row=0, column=0, sticky="w"
        )
        self.run_selected_button = ttk.Button(
            footer, text="Run selected", command=self._run_selected
        )
        self.run_selected_button.grid(row=0, column=1, padx=(8, 0))
        self.run_all_button = ttk.Button(
            footer, text="Run all", command=self._run_all
        )
        self.run_all_button.grid(row=0, column=2, padx=(8, 0))
        self.stop_button = ttk.Button(
            footer, text="Stop", command=self._stop_run, state=tk.DISABLED
        )
        self.stop_button.grid(row=0, column=3, padx=(8, 0))

    def _choose_executable(self) -> None:
        executable = filedialog.askopenfilename(
            title="Select the OneFLOW executable",
            filetypes=[
                ("OneFLOW executable", "OneFLOW.exe"),
                ("Executable files", "*.exe"),
                ("All files", "*"),
            ],
        )
        if executable:
            self.executable_var.set(str(Path(executable).resolve()))

    def _add_case(self) -> None:
        directory = filedialog.askdirectory(title="Select a OneFLOW case directory")
        if not directory:
            return
        case_directory = Path(directory).resolve()
        if case_directory in self.case_directories:
            messagebox.showinfo("Case already added", str(case_directory))
            return
        self.case_directories.append(case_directory)
        self._refresh_case_list()

    def _remove_cases(self) -> None:
        selected = list(self.case_list.curselection())
        for index in reversed(selected):
            del self.case_directories[index]
        self._refresh_case_list()

    def _move_case(self, direction: int) -> None:
        """Move one selected case so the dependency order remains unambiguous."""
        selected = list(self.case_list.curselection())
        if len(selected) != 1:
            messagebox.showinfo("Select one case", "Select exactly one case to reorder.")
            return

        old_index = selected[0]
        new_index = old_index + direction
        if not 0 <= new_index < len(self.case_directories):
            return

        self.case_directories[old_index], self.case_directories[new_index] = (
            self.case_directories[new_index],
            self.case_directories[old_index],
        )
        self._refresh_case_list()
        self.case_list.selection_set(new_index)
        self.case_list.activate(new_index)

    def _refresh_case_list(self) -> None:
        self.case_list.delete(0, tk.END)
        for index, case_directory in enumerate(self.case_directories, start=1):
            self.case_list.insert(tk.END, f"{index:02d}. {case_directory}")

    def _make_plan(self, cases: tuple[Path, ...]) -> RunPlan:
        executable_text = self.executable_var.get().strip()
        if not executable_text:
            raise ValueError("Select the OneFLOW executable first.")

        return RunPlan(
            executable=Path(executable_text).expanduser().resolve(),
            mode=self.mode_var.get().strip(),
            case_directories=cases,
            stop_on_error=self.stop_on_error_var.get(),
        )

    def _run_selected(self) -> None:
        selected_indices = self.case_list.curselection()
        if not selected_indices:
            messagebox.showinfo(
                "No cases selected", "Select one or more cases, or use Run all."
            )
            return
        selected_cases = tuple(self.case_directories[index] for index in selected_indices)
        self._start_run(selected_cases)

    def _run_all(self) -> None:
        self._start_run(tuple(self.case_directories))

    def _start_run(self, cases: tuple[Path, ...]) -> None:
        try:
            plan = self._make_plan(cases)
            plan.validate()
            self.runner.start(plan)
        except (OSError, ValueError, RuntimeError) as error:
            messagebox.showerror("Cannot start run", str(error))
            return

        self.active = True
        self._set_running_state(True)
        self._append_log(f"Starting run with {len(cases)} case(s).\n")

    def _stop_run(self) -> None:
        if self.active:
            self.status_var.set("Stopping current case...")
            self.runner.stop()

    def _set_running_state(self, running: bool) -> None:
        state = tk.DISABLED if running else tk.NORMAL
        self.run_selected_button.configure(state=state)
        self.run_all_button.configure(state=state)
        self.stop_button.configure(state=tk.NORMAL if running else tk.DISABLED)

    def _append_log(self, text: str) -> None:
        self.log_view.configure(state=tk.NORMAL)
        self.log_view.insert(tk.END, text)
        self.log_view.see(tk.END)
        self.log_view.configure(state=tk.DISABLED)

    def _poll_events(self) -> None:
        """Consume worker events on the Tk event loop; never touch widgets in a worker."""
        try:
            while True:
                event = self.events.get_nowait()
                event_type = event["type"]

                if event_type == "case_started":
                    command = subprocess.list2cmdline(event["command"])
                    self.status_var.set(
                        f"Case {event['index'] + 1}/{event['total']}: {event['case']}"
                    )
                    self._append_log(f"\n$ {command}\n")
                elif event_type == "log":
                    self._append_log(event["line"])
                elif event_type == "case_finished":
                    self._append_log(
                        f"\nCase exit code: {event['return_code']}\n"
                    )
                elif event_type == "error":
                    self._append_log(f"ERROR: {event.get('message', '')}\n")
                elif event_type == "finished":
                    self.active = False
                    self._set_running_state(False)
                    if event["stopped"]:
                        self.status_var.set("Run cancelled.")
                    elif event["return_code"] == 0:
                        self.status_var.set("All requested cases completed successfully.")
                    else:
                        self.status_var.set(
                            f"Run failed with exit code {event['return_code']}."
                        )
                    self._append_log(
                        f"Run finished with exit code {event['return_code']}.\n"
                    )
        except Empty:
            pass
        finally:
            self.root.after(100, self._poll_events)

    def _save_plan(self) -> None:
        try:
            plan = self._make_plan(tuple(self.case_directories))
        except ValueError as error:
            messagebox.showerror("Cannot save plan", str(error))
            return

        filename = filedialog.asksaveasfilename(
            title="Save OneFLOW run plan",
            defaultextension=".json",
            filetypes=[("JSON run plan", "*.json"), ("All files", "*")],
        )
        if not filename:
            return
        try:
            save_run_plan(Path(filename), plan)
        except OSError as error:
            messagebox.showerror("Cannot save plan", str(error))
            return
        self.status_var.set(f"Saved run plan: {filename}")

    def _load_plan(self) -> None:
        filename = filedialog.askopenfilename(
            title="Load OneFLOW run plan",
            filetypes=[("JSON run plan", "*.json"), ("All files", "*")],
        )
        if not filename:
            return
        try:
            plan = load_run_plan(Path(filename))
        except (OSError, ValueError, KeyError, TypeError) as error:
            messagebox.showerror("Cannot load plan", str(error))
            return

        self.executable_var.set(str(plan.executable))
        self.mode_var.set(plan.mode)
        self.stop_on_error_var.set(plan.stop_on_error)
        self.case_directories = list(plan.case_directories)
        self._refresh_case_list()
        self.status_var.set(f"Loaded run plan: {filename}")

    def _on_close(self) -> None:
        if self.active:
            should_close = messagebox.askyesno(
                "Run in progress",
                "A case is still running. Stop it and close the workbench?",
            )
            if not should_close:
                return
            self.runner.stop()
        self.root.destroy()


def main() -> None:
    root = tk.Tk()
    OneFlowWorkbench(root)
    root.mainloop()


if __name__ == "__main__":
    main()
