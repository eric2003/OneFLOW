"""Launch the OneFLOW desktop application as a Python module."""

from .main_window import run_desktop_app


if __name__ == "__main__":
    raise SystemExit(run_desktop_app())
