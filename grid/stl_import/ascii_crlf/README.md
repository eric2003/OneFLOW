# STL ASCII CRLF - native OneFLOW Case

This is a native OneFLOW Grid Case for the same two-triangle planar surface
used by `stl_ascii_basic`, with the STL text fixture stored using CRLF line
endings.

Run it through the standard application entry point:

    OneFLOW.exe r grid/stl/ascii_crlf

The Case intentionally has no `script/control.txt`. The runtime discovers
`dimension.txt`, `grid.txt`, and `simu.txt` from the script directory.
The C++ STL reader should produce the same four welded nodes and two TRI_3
cells as the LF fixture in `ascii_basic`, then write `results/model.cgns`.

This is a line-ending compatibility Case, not a separate parser or conversion
implementation.
