#include <gtest/gtest.h>
#include <stdexcept>

#include "CellTopo.h"
#include "FaceTopo.h"
#include "UnsGrid.h"
#include "Boundary.h"
#include "CellMesh.h"
#include "BcRecord.h"

TEST(FaceTopoConnectivityTest, ReordersAllFaceIndexedDataTogether)
{
    ONEFLOW::FaceTopo faces;
    faces.fTypes = ONEFLOW::IntField{10, 20, 30};
    faces.faces = ONEFLOW::LinkField{{10, 11}, {20, 21}, {30, 31}};
    faces.GetLeftCells() = ONEFLOW::IntField{1, 0, 2};
    faces.GetRightCells() = ONEFLOW::IntField{2, ONEFLOW::INVALID_INDEX, 3};
    faces.lPosition = ONEFLOW::IntField{11, 22, 33};
    faces.rPosition = ONEFLOW::IntField{111, 222, 333};
    faces.GetFaceFlags() = ONEFLOW::IntField{1, 2, 3};
    faces.c2f = ONEFLOW::LinkField{{0}, {1}};

    faces.ReorderLink();

    EXPECT_EQ(faces.fTypes, (ONEFLOW::IntField{20, 10, 30}));
    EXPECT_EQ(faces.faces,
        (ONEFLOW::LinkField{{20, 21}, {10, 11}, {30, 31}}));
    EXPECT_EQ(faces.GetLeftCells(), (ONEFLOW::IntField{0, 1, 2}));
    EXPECT_EQ(
        faces.GetRightCells(),
        (ONEFLOW::IntField{ONEFLOW::INVALID_INDEX, 2, 3}));
    EXPECT_EQ(faces.lPosition, (ONEFLOW::IntField{22, 11, 33}));
    EXPECT_EQ(faces.rPosition, (ONEFLOW::IntField{222, 111, 333}));
    EXPECT_EQ(faces.GetFaceFlags(), (ONEFLOW::IntField{2, 1, 3}));
    EXPECT_TRUE(faces.c2f.empty());
}

TEST(FaceTopoConnectivityTest, InvalidatesCellToFaceAdjacencyAfterTopologyUpdate)
{
    ONEFLOW::FaceTopo faces;
    faces.c2f = ONEFLOW::LinkField{{0, 2}, {1, 2}};

    faces.UpdateOtherTopologyTerm();

    EXPECT_TRUE(faces.c2f.empty());
}

TEST(CellTopoConnectivityTest, BuildsCellToFaceAdjacencyForBoundaryAndInternalFaces)
{
    ONEFLOW::CellTopo cells;
    cells.Alloc(2);

    ONEFLOW::FaceTopo faces;
    faces.SetNBFaces(2);
    faces.fTypes.resize(3);
    faces.GetLeftCells() = ONEFLOW::IntField{0, 1, 0};
    faces.GetRightCells() = ONEFLOW::IntField{2, 3, 1};

    ASSERT_NO_THROW(cells.CalcC2f(faces));

    const auto& c2f = cells.c2f;
    ASSERT_EQ(c2f.size(), 4u);
    EXPECT_EQ(c2f[0].size(), 2u);
    EXPECT_EQ(c2f[0][0], 0);
    EXPECT_EQ(c2f[0][1], 2);
    EXPECT_EQ(c2f[1].size(), 2u);
    EXPECT_EQ(c2f[1][0], 1);
    EXPECT_EQ(c2f[1][1], 2);
    EXPECT_TRUE(c2f[2].empty());
    EXPECT_TRUE(c2f[3].empty());
}

TEST(CellTopoConnectivityTest, RejectsInvalidCellReferenceWithoutCachingPartialAdjacency)
{
    ONEFLOW::CellTopo cells;
    cells.Alloc(1);

    ONEFLOW::FaceTopo faces;
    faces.fTypes.resize(1);
    faces.GetLeftCells() = ONEFLOW::IntField{1};
    faces.GetRightCells() = ONEFLOW::IntField{0};

    EXPECT_THROW(cells.CalcC2f(faces), std::out_of_range);
    EXPECT_EQ(cells.c2f.size(), 0u);
}

TEST(CellTopoConnectivityTest, RejectsSelfConnectedInternalFace)
{
    ONEFLOW::CellTopo cells;
    cells.Alloc(1);

    ONEFLOW::FaceTopo faces;
    faces.fTypes.resize(1);
    faces.GetLeftCells() = ONEFLOW::IntField{0};
    faces.GetRightCells() = ONEFLOW::IntField{0};

    EXPECT_THROW(cells.CalcC2f(faces), std::logic_error);
    EXPECT_EQ(cells.c2f.size(), 0u);
}

TEST(CellTopoConnectivityTest, RejectsMismatchedFaceCellArrays)
{
    ONEFLOW::CellTopo cells;
    cells.Alloc(1);

    ONEFLOW::FaceTopo faces;
    faces.fTypes.resize(1);
    faces.GetLeftCells() = ONEFLOW::IntField{0};

    EXPECT_THROW(cells.CalcC2f(faces), std::logic_error);
    EXPECT_EQ(cells.c2f.size(), 0u);
}


TEST(CellTopoConnectivityTest, RebuildsAdjacencyAfterFaceConnectivityChanges)
{
    ONEFLOW::CellTopo cells;
    cells.Alloc(3);

    ONEFLOW::FaceTopo faces;
    faces.SetNBFaces(0);
    faces.fTypes.resize(2);
    faces.GetLeftCells() = ONEFLOW::IntField{0, 1};
    faces.GetRightCells() = ONEFLOW::IntField{1, 2};

    ASSERT_NO_THROW(cells.CalcC2f(faces));
    ASSERT_EQ(cells.c2f.size(), 3u);
    ASSERT_EQ(cells.c2f[0], (ONEFLOW::IntField{0}));
    ASSERT_EQ(cells.c2f[1], (ONEFLOW::IntField{0, 1}));
    ASSERT_EQ(cells.c2f[2], (ONEFLOW::IntField{1}));

    faces.GetLeftCells()[1] = 0;
    faces.GetRightCells()[1] = 2;

    ASSERT_NO_THROW(cells.CalcC2f(faces));
    EXPECT_EQ(cells.c2f[0], (ONEFLOW::IntField{0, 1}));
    EXPECT_EQ(cells.c2f[1], (ONEFLOW::IntField{0}));
    EXPECT_EQ(cells.c2f[2], (ONEFLOW::IntField{1}));
}

TEST(FaceTopoConnectivityTest, RebuildsAdjacencyAfterFaceConnectivityChanges)
{
    ONEFLOW::UnsGrid grid;
    grid.Init();
    grid.nCells = 3;

    ONEFLOW::CellTopo & cells = grid.GetCellMesh().GetCellTopo();
    cells.Alloc(3);

    ONEFLOW::FaceTopo & faces = grid.GetFaceTopo();
    faces.SetNBFaces(0);
    faces.fTypes.resize(2);
    faces.GetLeftCells() = ONEFLOW::IntField{0, 1};
    faces.GetRightCells() = ONEFLOW::IntField{1, 2};

    ASSERT_NO_THROW(cells.CalcC2C(faces));
    ASSERT_EQ(cells.c2c.size(), 3u);
    EXPECT_EQ(cells.c2c[0], (ONEFLOW::IntField{1}));
    EXPECT_EQ(cells.c2c[1], (ONEFLOW::IntField{0, 2}));
    EXPECT_EQ(cells.c2c[2], (ONEFLOW::IntField{1}));

    faces.GetLeftCells()[1] = 0;
    faces.GetRightCells()[1] = 2;

    ASSERT_NO_THROW(cells.CalcC2C(faces));
    EXPECT_EQ(cells.c2c[0], (ONEFLOW::IntField{1, 2}));
    EXPECT_EQ(cells.c2c[1], (ONEFLOW::IntField{0}));
    EXPECT_EQ(cells.c2c[2], (ONEFLOW::IntField{0}));
}

TEST(FaceTopoConnectivityTest, BuildsReciprocalInternalNeighborsAndKeepsInterfaceGhostIds)
{
    ONEFLOW::UnsGrid grid;
    grid.Init();
    grid.nCells = 2;

    ONEFLOW::FaceTopo & faces = grid.GetFaceTopo();
    faces.SetNBFaces(2);
    faces.fTypes.resize(3);
    faces.GetBcRecord().bcType = {
        ONEFLOW::BC::INTERFACE,
        ONEFLOW::BC::SOLID_SURFACE
    };
    faces.GetBcRecord().bcNameId.resize(2);
    faces.GetLeftCells() = ONEFLOW::IntField{0, 1, 0};
    faces.GetRightCells() = ONEFLOW::IntField{2, 3, 1};

    ONEFLOW::LinkField c2c;
    ASSERT_NO_THROW(faces.CalcC2C(c2c));

    ASSERT_EQ(c2c.size(), 2u);
    ASSERT_EQ(c2c[0].size(), 2u);
    EXPECT_EQ(c2c[0][0], 2);
    EXPECT_EQ(c2c[0][1], 1);
    ASSERT_EQ(c2c[1].size(), 1u);
    EXPECT_EQ(c2c[1][0], 0);
}

TEST(FaceTopoConnectivityTest, RejectsInvalidInternalCellIdsWithoutCachingPartialAdjacency)
{
    ONEFLOW::UnsGrid grid;
    grid.Init();
    grid.nCells = 2;

    ONEFLOW::FaceTopo & faces = grid.GetFaceTopo();
    faces.SetNBFaces(0);
    faces.fTypes.resize(1);
    faces.GetLeftCells() = ONEFLOW::IntField{0};
    faces.GetRightCells() = ONEFLOW::IntField{2};

    ONEFLOW::LinkField c2c;
    EXPECT_THROW(faces.CalcC2C(c2c), std::out_of_range);
    EXPECT_TRUE(c2c.empty());

    faces.GetRightCells()[0] = 0;
    EXPECT_THROW(faces.CalcC2C(c2c), std::logic_error);
    EXPECT_TRUE(c2c.empty());
}

TEST(CellMeshConnectivityTest, RecalculatesCellSpanAfterBoundaryConnectivityChanges)
{
    ONEFLOW::UnsGrid grid;
    grid.Init();
    grid.nCells = 2;
    grid.nFaces = 3;
    grid.nBFaces = 2;

    ONEFLOW::CellMesh & cellMesh = grid.GetCellMesh();
    cellMesh.GetCellTopo().Alloc(2);

    ONEFLOW::FaceTopo & faces = grid.GetFaceTopo();
    faces.SetNBFaces(2);
    faces.fTypes.resize(3);
    faces.GetLeftCells() = ONEFLOW::IntField{0, 1, 0};
    faces.GetRightCells() = ONEFLOW::IntField{2, 3, 1};

    cellMesh.xcc = ONEFLOW::RealField{0.0, 1.0, 10.0, 4.0};
    cellMesh.ycc = ONEFLOW::RealField(4, 0.0);
    cellMesh.zcc = ONEFLOW::RealField(4, 0.0);

    cellMesh.CalcCellSpan(grid);
    ASSERT_EQ(cellMesh.span.size(), 2u);
    EXPECT_DOUBLE_EQ(cellMesh.span[0], 10.0);
    EXPECT_DOUBLE_EQ(cellMesh.span[1], 3.0);

    faces.GetLeftCells()[0] = 1;
    cellMesh.CalcCellSpan(grid);

    EXPECT_DOUBLE_EQ(cellMesh.span[0], 1.0);
    EXPECT_DOUBLE_EQ(cellMesh.span[1], 9.0);
}

TEST(CellTopoConnectivityTest, CellCountApiMatchesCompatibilityWrapper)
{
    // Intent: GetCellCount is preferred; GetNumberOfCells is a compatibility
    // alias that must report the same value after Alloc and after PushElement.
    ONEFLOW::CellTopo cellTopo;
    EXPECT_EQ(cellTopo.GetCellCount(), static_cast<ONEFLOW::HXSize_t>(0));
    EXPECT_EQ(cellTopo.GetNumberOfCells(), cellTopo.GetCellCount());

    cellTopo.Alloc(4);
    EXPECT_EQ(cellTopo.GetCellCount(), static_cast<ONEFLOW::HXSize_t>(4));
    EXPECT_EQ(cellTopo.GetNumberOfCells(), cellTopo.GetCellCount());
}
