#!/usr/bin/env python3
"""Local graphical front end for the shared OneFLOW regression test runner.

Run from any directory with: python path/to/test/gui.py
The GUI delegates execution and assertions to test.py; it does not implement
a second test engine.
"""

from __future__ import annotations

import os
import queue
import subprocess
import sys
import threading
import tkinter as tk
from pathlib import Path
from tkinter import filedialog, messagebox, ttk


TEST_DIR = Path(__file__).resolve().parent
RUNNER = TEST_DIR / "test.py"


class OneFLOWTestGUI:
    def __init__(self, root: tk.Tk) -> None:
        self.root = root
        self.root.title("OneFLOW Regression Test Runner")
        self.root.geometry("1020x720")
        self.root.minsize(760, 520)
        self.messages: queue.Queue = queue.Queue()
        self.process: subprocess.Popen | None = None
        self.worker: threading.Thread | None = None

        suite_names = sorted(path.name for path in (TEST_DIR / "suites").glob("*.txt"))
        default_suite = "ci-serial.txt" if "ci-serial.txt" in suite_names else (
            suite_names[0] if suite_names else ""
        )
        default_exe = os.environ.get("ONEFLOW_EXE", "")
        default_mpi = "mpiexec -n 1" if os.name == "nt" else "mpirun -np 1"

        self.suite_var = tk.StringVar(value=default_suite)
        self.mpi_var = tk.StringVar(value=default_mpi)
        self.exe_var = tk.StringVar(value=default_exe)
        self.residual_var = tk.StringVar(value="")
        self.tolerance_var = tk.StringVar(value="1e-8")
        self.timeout_var = tk.StringVar(value="")

        self._build_layout(suite_names)
        self.root.after(100, self._poll_messages)
        self.root.protocol("WM_DELETE_WINDOW", self._on_close)

    def _build_layout(self, suite_names: list[str]) -> None:
        outer = ttk.Frame(self.root, padding=14)
        outer.pack(fill=tk.BOTH, expand=True)
        ttk.Label(
            outer,
            text="OneFLOW Regression Tests",
            font=("Segoe UI", 16, "bold"),
        ).pack(anchor=tk.W)
        ttk.Label(
            outer,
            text="Cross-platform GUI · uses the same suite files and test runner as CI",
        ).pack(anchor=tk.W, pady=(2, 12))

        config = ttk.LabelFrame(outer, text="Run configuration", padding=10)
        config.pack(fill=tk.X)
        config.columnconfigure(1, weight=1)

        ttk.Label(config, text="Suite").grid(row=0, column=0, sticky=tk.W, padx=(0, 8), pady=4)
        suite_box = ttk.Combobox(
            config, textvariable=self.suite_var, values=suite_names, state="readonly"
        )
        suite_box.grid(row=0, column=1, columnspan=2, sticky=tk.EW, pady=4)

        ttk.Label(config, text="MPI command").grid(row=1, column=0, sticky=tk.W, padx=(0, 8), pady=4)
        ttk.Entry(config, textvariable=self.mpi_var).grid(row=1, column=1, columnspan=2, sticky=tk.EW, pady=4)

        ttk.Label(config, text="OneFLOW executable").grid(row=2, column=0, sticky=tk.W, padx=(0, 8), pady=4)
        ttk.Entry(config, textvariable=self.exe_var).grid(row=2, column=1, sticky=tk.EW, pady=4)
        ttk.Button(config, text="Browse…", command=self._browse_executable).grid(row=2, column=2, padx=(8, 0), pady=4)

        ttk.Label(config, text="Residual baseline (optional)").grid(row=3, column=0, sticky=tk.W, padx=(0, 8), pady=4)
        ttk.Entry(config, textvariable=self.residual_var).grid(row=3, column=1, sticky=tk.EW, pady=4)
        ttk.Button(config, text="Browse…", command=self._browse_residual).grid(row=3, column=2, padx=(8, 0), pady=4)

        options = ttk.Frame(config)
        options.grid(row=4, column=0, columnspan=3, sticky=tk.EW, pady=(6, 0))
        ttk.Label(options, text="Residual tolerance").pack(side=tk.LEFT)
        ttk.Entry(options, textvariable=self.tolerance_var, width=12).pack(side=tk.LEFT, padx=(6, 18))
        ttk.Label(options, text="Per-case timeout (seconds; blank = unlimited)").pack(side=tk.LEFT)
        ttk.Entry(options, textvariable=self.timeout_var, width=12).pack(side=tk.LEFT, padx=6)

        actions = ttk.Frame(outer)
        actions.pack(fill=tk.X, pady=10)
        self.start_button = ttk.Button(actions, text="Run tests", command=self._start)
        self.start_button.pack(side=tk.LEFT)
        self.stop_button = ttk.Button(actions, text="Stop run", command=self._stop, state=tk.DISABLED)
        self.stop_button.pack(side=tk.LEFT, padx=(8, 0))
        self.status_var = tk.StringVar(value="Ready")
        ttk.Label(actions, textvariable=self.status_var).pack(side=tk.RIGHT)

        log_frame = ttk.LabelFrame(outer, text="Live test log and report", padding=6)
        log_frame.pack(fill=tk.BOTH, expand=True)
        self.log = tk.Text(log_frame, wrap=tk.NONE, state=tk.DISABLED, font=("Consolas", 10))
        y_scroll = ttk.Scrollbar(log_frame, orient=tk.VERTICAL, command=self.log.yview)
        x_scroll = ttk.Scrollbar(log_frame, orient=tk.HORIZONTAL, command=self.log.xview)
        self.log.configure(yscrollcommand=y_scroll.set, xscrollcommand=x_scroll.set)
        self.log.grid(row=0, column=0, sticky=tk.NSEW)
        y_scroll.grid(row=0, column=1, sticky=tk.NS)
        x_scroll.grid(row=1, column=0, sticky=tk.EW)
        log_frame.rowconfigure(0, weight=1)
        log_frame.columnconfigure(0, weight=1)

    def _browse_executable(self) -> None:
        path = filedialog.askopenfilename(title="Select OneFLOW executable")
        if path:
            self.exe_var.set(path)

    def _browse_residual(self) -> None:
        path = filedialog.askopenfilename(
            title="Select residual baseline JSON",
            filetypes=(("JSON files", "*.json"), ("All files", "*.*")),
        )
        if path:
            self.residual_var.set(path)

    def _append_log(self, text: str) -> None:
        self.log.configure(state=tk.NORMAL)
        self.log.insert(tk.END, text)
        self.log.see(tk.END)
        self.log.configure(state=tk.DISABLED)

    def _start(self) -> None:
        if self.process is not None:
            return
        suite = self.suite_var.get().strip()
        mpi = self.mpi_var.get().strip()
        executable = self.exe_var.get().strip()
        if not suite or not mpi or not executable:
            messagebox.showerror(
                "Missing configuration",
                "Select a suite and provide both the MPI command and OneFLOW executable.",
            )
            return
        if not (TEST_DIR / "suites" / suite).is_file():
            messagebox.showerror("Invalid suite", f"Suite file not found: {suite}")
            return

        # test.py currently accepts the MPI launcher as a shell command string.
        if any(char in executable for char in (" ", "\t")) and not (
            executable.startswith('"') and executable.endswith('"')
        ):
            executable = f'"{executable}"'
        command = [sys.executable, str(RUNNER), mpi, executable, f"suites/{suite}"]
        residual = self.residual_var.get().strip()
        tolerance = self.tolerance_var.get().strip()
        if residual:
            command.append(residual)
            if tolerance:
                command.append(tolerance)
        timeout = self.timeout_var.get().strip()
        if timeout:
            try:
                if float(timeout) <= 0:
                    raise ValueError
            except ValueError:
                messagebox.showerror("Invalid timeout", "Timeout must be a positive number.")
                return
            command.extend(["--timeout", timeout])

        self._append_log("\n" + "=" * 72 + "\n")
        self._append_log("Command: " + subprocess.list2cmdline(command) + "\n")
        self._append_log("=" * 72 + "\n")
        self.status_var.set("Running…")
        self.start_button.configure(state=tk.DISABLED)
        self.stop_button.configure(state=tk.NORMAL)
        self.process = subprocess.Popen(
            command,
            cwd=str(TEST_DIR),
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            encoding="utf-8",
            errors="replace",
            bufsize=1,
        )
        self.worker = threading.Thread(target=self._read_output, daemon=True)
        self.worker.start()

    def _read_output(self) -> None:
        process = self.process
        if process is None or process.stdout is None:
            return
        for line in process.stdout:
            self.messages.put(("line", line))
        return_code = process.wait()
        self.messages.put(("finished", return_code))

    def _poll_messages(self) -> None:
        try:
            while True:
                kind, value = self.messages.get_nowait()
                if kind == "line":
                    self._append_log(value)
                elif kind == "finished":
                    self.process = None
                    self.start_button.configure(state=tk.NORMAL)
                    self.stop_button.configure(state=tk.DISABLED)
                    if value == 0:
                        self.status_var.set("Passed")
                        self._append_log("\nGUI status: test runner exited successfully.\n")
                    else:
                        self.status_var.set(f"Failed (exit code {value})")
                        self._append_log(f"\nGUI status: test runner exited with code {value}.\n")
        except queue.Empty:
            pass
        self.root.after(100, self._poll_messages)

    def _stop(self) -> None:
        if self.process is not None and self.process.poll() is None:
            self.process.terminate()
            self.status_var.set("Stopping…")

    def _on_close(self) -> None:
        if self.process is not None and self.process.poll() is None:
            if not messagebox.askyesno("Exit", "A test run is active. Stop it and exit?"):
                return
            self.process.terminate()
        self.root.destroy()


def main() -> None:
    root = tk.Tk()
    OneFLOWTestGUI(root)
    root.mainloop()


if __name__ == "__main__":
    main()
