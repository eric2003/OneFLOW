"""Architecture contract for the shared unstructured CGNS writer."""

from __future__ import annotations

import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]


class UnstructuredCgnsWriterArchitectureTest(unittest.TestCase):
    def test_gmsh_builds_common_mesh_data_and_calls_common_writer(self) -> None:
        source = (ROOT / "codes/geometry/src/GridConversion.cpp").read_text(
            encoding="utf-8"
        )
        start = source.index("UnstructuredMeshData BuildUnstructuredMeshFromGmsh")
        end = source.index("void ConvertGmshToCgns", start)
        adapter = source[start:end]
        self.assertIn("UnstructuredMeshData grid;", adapter)
        self.assertIn("UnstructuredCgnsZoneWriter::Write", adapter)
        self.assertNotIn("Su2Grid grid", adapter)
        self.assertNotIn("FillSU2CgnsZone", adapter)

    def test_su2_reader_delegates_zone_population_to_common_writer(self) -> None:
        source = (ROOT / "codes/geometry/src/Su2Grid.cpp").read_text(
            encoding="utf-8"
        )
        start = source.index("void Su2Grid::FillSU2CgnsZone")
        end = source.index("void Su2ToOneFlowGrid", start)
        adapter = source[start:end]
        self.assertIn("UnstructuredCgnsZoneWriter::Write", adapter)

    def test_common_writer_has_no_source_format_dependency(self) -> None:
        header = (
            ROOT / "codes/geometry/include/UnstructuredCgnsZoneWriter.h"
        ).read_text(encoding="utf-8")
        source = (
            ROOT / "codes/geometry/src/UnstructuredCgnsZoneWriter.cpp"
        ).read_text(encoding="utf-8")
        self.assertNotIn("Su2Grid.h", header)
        self.assertNotIn("Su2Grid.h", source)
        self.assertIn("cgnsZone.cgnsZoneType = Unstructured", source)
        self.assertIn("mesh.mmark", source)


if __name__ == "__main__":
    unittest.main()
