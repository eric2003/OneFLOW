#include <gtest/gtest.h>
#include <stdexcept>

#include "ScalarGrid.h"
#include "Boundary.h"
#include "ScalarIFace.h"

TEST(ScalarGridTopologyTest, CalcC2CRejectsInvalidInternalCellReferences)
{
    ONEFLOW::ScalarGrid grid;
    grid.eTypes.AddData(1);
    grid.faces.Resize(1);
    grid.lc.AddData(0);
    grid.rc.AddData(1);

    ONEFLOW::EList adjacency;
    EXPECT_THROW(grid.CalcC2C(adjacency), std::runtime_error);
    EXPECT_EQ(adjacency.GetNElements(), 0u);
}

TEST(ScalarGridTopologyTest, CalcC2CRejectsSelfConnectedInternalFace)
{
    ONEFLOW::ScalarGrid grid;
    grid.eTypes.AddData(1);
    grid.faces.Resize(1);
    grid.lc.AddData(0);
    grid.rc.AddData(0);

    ONEFLOW::EList adjacency;
    EXPECT_THROW(grid.CalcC2C(adjacency), std::runtime_error);
    EXPECT_EQ(adjacency.GetNElements(), 0u);
}

TEST(ScalarGridTopologyTest, CalcC2CBuildsSymmetricInternalAdjacency)
{
    ONEFLOW::ScalarGrid grid;
    grid.eTypes.AddData(1);
    grid.eTypes.AddData(1);
    grid.faces.Resize(1);
    grid.lc.AddData(0);
    grid.rc.AddData(1);

    ONEFLOW::EList adjacency;
    grid.CalcC2C(adjacency);

    ASSERT_EQ(adjacency.GetNElements(), 2u);
    ASSERT_EQ(adjacency[0].size(), 1u);
    ASSERT_EQ(adjacency[1].size(), 1u);
    EXPECT_EQ(adjacency[0][0], 1);
    EXPECT_EQ(adjacency[1][0], 0);
}

TEST(ScalarGridTopologyTest, CalcInterfaceToBcFaceRejectsInconsistentFaceArrays)
{
    ONEFLOW::ScalarGrid grid;
    grid.faces.Resize(1);
    grid.lc.AddData(0);
    grid.rc.AddData(1);
    grid.bcTypes.AddData(ONEFLOW::BC::INTERFACE);
    grid.scalarIFace->interface_to_bcface = {42};

    EXPECT_THROW(grid.CalcInterfaceToBcFace(), std::runtime_error);
    ASSERT_EQ(grid.scalarIFace->interface_to_bcface.size(), 1u);
    EXPECT_EQ(grid.scalarIFace->interface_to_bcface[0], 42);
}

TEST(ScalarGridTopologyTest, CalcInterfaceToBcFaceRejectsMissingInterfaceTopology)
{
    ONEFLOW::ScalarGrid grid;
    grid.faces.Resize(1);
    grid.lc.AddData(0);
    grid.rc.AddData(1);
    grid.fBcTypes.AddData(ONEFLOW::BC::INTERFACE);
    grid.bcTypes.AddData(ONEFLOW::BC::INTERFACE);
    grid.scalarIFace.reset();

    EXPECT_THROW(grid.CalcInterfaceToBcFace(), std::logic_error);
}
