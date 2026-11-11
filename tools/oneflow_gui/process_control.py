"""Cross-platform process-group creation and process-tree termination."""

from __future__ import annotations

import ctypes
import os
import signal
import subprocess
from collections.abc import Callable
from typing import Any


CommandRunner = Callable[..., Any]


def process_group_options(*, platform_name: str | None = None) -> dict[str, object]:
    """Return Popen options that isolate a command and its descendants."""
    platform_value = (platform_name or os.name).lower()
    if platform_value in {"nt", "windows"}:
        return {"creationflags": getattr(subprocess, "CREATE_NEW_PROCESS_GROUP", 0x00000200)}
    return {"start_new_session": True}


def _process_job(process: subprocess.Popen[object]) -> tuple[Any, Any] | None:
    """Read only a job explicitly attached by this module.

    Popen and test doubles may synthesize arbitrary attributes from getattr;
    inspecting the instance dictionary avoids treating a Mock child as a
    real (kernel32, HANDLE) pair.
    """
    job = vars(process).get("_oneflow_job")
    if (
        isinstance(job, tuple)
        and len(job) == 2
        and hasattr(job[0], "CloseHandle")
    ):
        return job
    return None


def attach_process_job(
    process: subprocess.Popen[object],
    *,
    platform_name: str | None = None,
) -> bool:
    """Attach a Windows process to a kill-on-close Job Object when available.

    Descendants inherit the job automatically. The Job Object also guarantees
    that remaining descendants are terminated when the workflow process exits.
    """
    platform_value = (platform_name or os.name).lower()
    if platform_value not in {"nt", "windows"}:
        return False

    from ctypes import wintypes

    class BasicLimitInformation(ctypes.Structure):
        _fields_ = [
            ("PerProcessUserTimeLimit", ctypes.c_longlong),
            ("PerJobUserTimeLimit", ctypes.c_longlong),
            ("LimitFlags", wintypes.DWORD),
            ("MinimumWorkingSetSize", ctypes.c_size_t),
            ("MaximumWorkingSetSize", ctypes.c_size_t),
            ("ActiveProcessLimit", wintypes.DWORD),
            ("Affinity", ctypes.c_size_t),
            ("PriorityClass", wintypes.DWORD),
            ("SchedulingClass", wintypes.DWORD),
        ]

    class IoCounters(ctypes.Structure):
        _fields_ = [
            ("ReadOperationCount", ctypes.c_ulonglong),
            ("WriteOperationCount", ctypes.c_ulonglong),
            ("OtherOperationCount", ctypes.c_ulonglong),
            ("ReadTransferCount", ctypes.c_ulonglong),
            ("WriteTransferCount", ctypes.c_ulonglong),
            ("OtherTransferCount", ctypes.c_ulonglong),
        ]

    class ExtendedLimitInformation(ctypes.Structure):
        _fields_ = [
            ("BasicLimitInformation", BasicLimitInformation),
            ("IoInfo", IoCounters),
            ("ProcessMemoryLimit", ctypes.c_size_t),
            ("JobMemoryLimit", ctypes.c_size_t),
            ("PeakProcessMemoryUsed", ctypes.c_size_t),
            ("PeakJobMemoryUsed", ctypes.c_size_t),
        ]

    try:
        kernel32 = ctypes.WinDLL("kernel32", use_last_error=True)
        kernel32.CreateJobObjectW.argtypes = [wintypes.LPVOID, wintypes.LPCWSTR]
        kernel32.CreateJobObjectW.restype = wintypes.HANDLE
        kernel32.SetInformationJobObject.argtypes = [
            wintypes.HANDLE, ctypes.c_int, wintypes.LPVOID, wintypes.DWORD
        ]
        kernel32.SetInformationJobObject.restype = wintypes.BOOL
        kernel32.AssignProcessToJobObject.argtypes = [wintypes.HANDLE, wintypes.HANDLE]
        kernel32.AssignProcessToJobObject.restype = wintypes.BOOL
        kernel32.TerminateJobObject.argtypes = [wintypes.HANDLE, wintypes.UINT]
        kernel32.TerminateJobObject.restype = wintypes.BOOL
        kernel32.CloseHandle.argtypes = [wintypes.HANDLE]
        kernel32.CloseHandle.restype = wintypes.BOOL

        job_handle = kernel32.CreateJobObjectW(None, None)
        if not job_handle:
            return False

        information = ExtendedLimitInformation()
        information.BasicLimitInformation.LimitFlags = 0x00002000  # KILL_ON_JOB_CLOSE
        configured = kernel32.SetInformationJobObject(
            job_handle,
            9,  # JobObjectExtendedLimitInformation
            ctypes.byref(information),
            ctypes.sizeof(information),
        )
        assigned = configured and kernel32.AssignProcessToJobObject(
            job_handle, wintypes.HANDLE(int(process._handle))
        )
        if not assigned:
            kernel32.CloseHandle(job_handle)
            return False

        process._oneflow_job = (kernel32, job_handle)
        return True
    except (AttributeError, OSError, TypeError, ValueError):
        return False


def close_process_job(process: subprocess.Popen[object]) -> bool:
    """Close an attached Job Object, terminating any remaining descendants."""
    job = _process_job(process)
    if job is None:
        return False
    kernel32, job_handle = job
    try:
        kernel32.CloseHandle(job_handle)
    finally:
        delattr(process, "_oneflow_job")
    return True



def _terminate_windows_descendants(parent_pid: int) -> None:
    """Terminate descendants still alive after the root process was stopped.

    taskkill /T and Job Object assignment can both miss a child created during
    process startup. A Toolhelp snapshot taken after stopping the root closes
    that race: descendants retain their parent PID even after the parent exits.
    """
    if os.name != "nt":
        return

    from ctypes import wintypes

    class ProcessEntry32W(ctypes.Structure):
        _fields_ = [
            ("dwSize", wintypes.DWORD),
            ("cntUsage", wintypes.DWORD),
            ("th32ProcessID", wintypes.DWORD),
            ("th32DefaultHeapID", ctypes.c_size_t),
            ("th32ModuleID", wintypes.DWORD),
            ("cntThreads", wintypes.DWORD),
            ("th32ParentProcessID", wintypes.DWORD),
            ("pcPriClassBase", ctypes.c_long),
            ("dwFlags", wintypes.DWORD),
            ("szExeFile", wintypes.WCHAR * 260),
        ]

    kernel32 = ctypes.WinDLL("kernel32", use_last_error=True)
    kernel32.CreateToolhelp32Snapshot.argtypes = [wintypes.DWORD, wintypes.DWORD]
    kernel32.CreateToolhelp32Snapshot.restype = wintypes.HANDLE
    kernel32.Process32FirstW.argtypes = [
        wintypes.HANDLE, ctypes.POINTER(ProcessEntry32W)
    ]
    kernel32.Process32FirstW.restype = wintypes.BOOL
    kernel32.Process32NextW.argtypes = [
        wintypes.HANDLE, ctypes.POINTER(ProcessEntry32W)
    ]
    kernel32.Process32NextW.restype = wintypes.BOOL
    kernel32.OpenProcess.argtypes = [
        wintypes.DWORD, wintypes.BOOL, wintypes.DWORD
    ]
    kernel32.OpenProcess.restype = wintypes.HANDLE
    kernel32.TerminateProcess.argtypes = [wintypes.HANDLE, wintypes.UINT]
    kernel32.TerminateProcess.restype = wintypes.BOOL
    kernel32.CloseHandle.argtypes = [wintypes.HANDLE]
    kernel32.CloseHandle.restype = wintypes.BOOL

    snapshot = kernel32.CreateToolhelp32Snapshot(0x00000002, 0)
    if snapshot == wintypes.HANDLE(-1).value:
        return

    try:
        children: dict[int, list[int]] = {}
        entry = ProcessEntry32W()
        entry.dwSize = ctypes.sizeof(entry)
        if kernel32.Process32FirstW(snapshot, ctypes.byref(entry)):
            while True:
                children.setdefault(int(entry.th32ParentProcessID), []).append(
                    int(entry.th32ProcessID)
                )
                if not kernel32.Process32NextW(snapshot, ctypes.byref(entry)):
                    break
    finally:
        kernel32.CloseHandle(snapshot)

    descendants: list[int] = []
    pending = list(children.get(parent_pid, ()))
    seen = {parent_pid}
    while pending:
        process_id = pending.pop()
        if process_id in seen:
            continue
        seen.add(process_id)
        descendants.append(process_id)
        pending.extend(children.get(process_id, ()))

    # Terminate deepest descendants first so they cannot keep creating children.
    for process_id in reversed(descendants):
        handle = kernel32.OpenProcess(0x0001, False, process_id)
        if handle:
            try:
                kernel32.TerminateProcess(handle, 1)
            finally:
                kernel32.CloseHandle(handle)

def terminate_process_tree(
    process: subprocess.Popen[object],
    *,
    platform_name: str | None = None,
    command_runner: CommandRunner | None = None,
) -> bool:
    """Force-stop a process and its descendants; return whether it was running."""
    if process.poll() is not None:
        return False

    platform_value = (platform_name or os.name).lower()
    if platform_value in {"nt", "windows"}:
        # A child can be spawned in the short interval between Popen returning
        # and attach_process_job assigning the parent to its Job Object. Use
        # taskkill while the parent is still alive so Windows enumerates its
        # current process tree, including any child that escaped that race.
        runner = command_runner or subprocess.run
        try:
            result = runner(
                ["taskkill", "/PID", str(process.pid), "/T", "/F"],
                stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL,
                check=False,
            )
        except OSError:
            result = None
        if result is not None and result.returncode == 0:
            # taskkill may report success while a child created during startup
            # remains alive. The root is now stopped, so its descendants can
            # be enumerated reliably by their recorded parent-process IDs.
            _terminate_windows_descendants(process.pid)
            return True

        # Keep the Job Object as a second termination mechanism when taskkill
        # is unavailable or cannot terminate the process tree.
        job = _process_job(process)
        if job is not None and job[0].TerminateJobObject(job[1], 1):
            _terminate_windows_descendants(process.pid)
            return True
        if process.poll() is None:
            process.kill()
        _terminate_windows_descendants(process.pid)
        return True

    try:
        os.killpg(process.pid, signal.SIGKILL)
    except ProcessLookupError:
        pass
    except OSError:
        process.kill()
    return True


def terminate_and_reap(
    process: subprocess.Popen[object],
    *,
    timeout: float = 3.0,
    platform_name: str | None = None,
    command_runner: CommandRunner | None = None,
) -> bool:
    """Terminate a running process tree and wait for the direct child to exit."""
    if process.poll() is not None:
        return False

    terminate_process_tree(
        process,
        platform_name=platform_name,
        command_runner=command_runner,
    )
    try:
        process.wait(timeout=timeout)
    except subprocess.TimeoutExpired:
        terminate_process_tree(
            process,
            platform_name=platform_name,
            command_runner=command_runner,
        )
        if process.poll() is None:
            process.kill()
        process.wait()
    return True
