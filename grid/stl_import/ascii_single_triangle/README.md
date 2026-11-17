# STL ASCII Single Triangle - native OneFLOW Case

This is a native C++ OneFLOW grid case containing one valid ASCII STL facet. It exercises the minimum non-empty surface input through the production STL-to-CGNS path.

Run it from the repository root:

    OneFLOW.exe r grid/stl_import/ascii_single_triangle

The case uses the standard parameter files under `script/`; generated CGNS and STL read-report files are written under `results/` in the staged workspace. The expected mesh contains three distinct vertices and one TRI_3 element.
