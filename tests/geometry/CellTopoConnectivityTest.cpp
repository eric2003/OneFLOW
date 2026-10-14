#include <gtest/gtest.h>

#include "CellTopo.h"
#include "FaceTopo.h"
#include "UnsGrid.h"
#include "Boundary.h"

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
