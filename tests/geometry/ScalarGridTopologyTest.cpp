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

TEST(ScalarGridTopologyTest, CalcTopologyRejectsNullBoundaryConditionWithoutResettingTopology)
{
    ONEFLOW::ScalarGrid grid;
    grid.xn.Resize(2);
    grid.yn.Resize(2);
    grid.zn.Resize(2);
    grid.faces.Resize(1);
    grid.lc.AddData(37);
    grid.scalarBccos->bccos.push_back(nullptr);

    EXPECT_THROW(grid.CalcTopology(), std::runtime_error);
    ASSERT_EQ(grid.faces.GetNElements(), 1u);
    ASSERT_EQ(grid.lc.GetNElements(), 1u);
    EXPECT_EQ(grid.lc[0], 37);
}

TEST(ScalarGridTopologyTest, CalcTopologyRejectsInvalidBoundaryNodeWithoutResettingTopology)
{
    ONEFLOW::ScalarGrid grid;
    grid.xn.Resize(2);
    grid.yn.Resize(2);
    grid.zn.Resize(2);
    grid.faces.Resize(1);
    grid.lc.AddData(41);

    auto boundaryCondition = std::make_unique< ONEFLOW::ScalarBcco >();
    boundaryCondition->AddBcPoint(2);
    grid.scalarBccos->AddBcco(std::move(boundaryCondition));

    EXPECT_THROW(grid.CalcTopology(), std::runtime_error);
    ASSERT_EQ(grid.faces.GetNElements(), 1u);
    ASSERT_EQ(grid.lc.GetNElements(), 1u);
    EXPECT_EQ(grid.lc[0], 41);
}

TEST(ScalarGridTopologyTest, CalcTopologyRejectsMismatchedCoordinateArraysWithoutResettingTopology)
{
    ONEFLOW::ScalarGrid grid;
    grid.xn.Resize(2);
    grid.yn.Resize(1);
    grid.zn.Resize(2);
    grid.faces.Resize(1);
    grid.lc.AddData(31);

    EXPECT_THROW(grid.CalcTopology(), std::runtime_error);
    ASSERT_EQ(grid.faces.GetNElements(), 1u);
    ASSERT_EQ(grid.lc.GetNElements(), 1u);
    EXPECT_EQ(grid.lc[0], 31);
}

TEST(ScalarGridTopologyTest, CalcTopologyRejectsMismatchedCellArraysWithoutResettingTopology)
{
    ONEFLOW::ScalarGrid grid;
    grid.eTypes.AddData(ONEFLOW::BAR_2);
    grid.faces.Resize(1);
    grid.lc.AddData(17);

    EXPECT_THROW(grid.CalcTopology(), std::runtime_error);
    ASSERT_EQ(grid.faces.GetNElements(), 1u);
    ASSERT_EQ(grid.lc.GetNElements(), 1u);
    EXPECT_EQ(grid.lc[0], 17);
}

TEST(ScalarGridTopologyTest, CalcTopologyRejectsOutOfRangeNodeBeforeResettingTopology)
{
    ONEFLOW::ScalarGrid grid;
    grid.xn.Resize(2);
    grid.elements.AddElem(std::vector<int>{0, 2});
    grid.eTypes.AddData(ONEFLOW::BAR_2);
    grid.faces.Resize(1);
    grid.lc.AddData(17);

    EXPECT_THROW(grid.CalcTopology(), std::runtime_error);
    ASSERT_EQ(grid.faces.GetNElements(), 1u);
    ASSERT_EQ(grid.lc.GetNElements(), 1u);
    EXPECT_EQ(grid.lc[0], 17);
}

TEST(ScalarGridTopologyTest, CalcTopologyRejectsNonManifoldFaceWithoutResettingTopology)
{
    ONEFLOW::ScalarGrid grid;
    grid.xn.Resize(4);
    grid.elements.AddElem(std::vector<int>{0, 1});
    grid.elements.AddElem(std::vector<int>{0, 2});
    grid.elements.AddElem(std::vector<int>{0, 3});
    grid.eTypes.AddData(ONEFLOW::BAR_2);
    grid.eTypes.AddData(ONEFLOW::BAR_2);
    grid.eTypes.AddData(ONEFLOW::BAR_2);
    grid.faces.Resize(1);
    grid.lc.AddData(23);

    EXPECT_THROW(grid.CalcTopology(), std::runtime_error);
    ASSERT_EQ(grid.faces.GetNElements(), 1u);
    ASSERT_EQ(grid.lc.GetNElements(), 1u);
    EXPECT_EQ(grid.lc[0], 23);
}

TEST(ScalarGridTopologyTest, CalcTopologyRejectsDuplicateCellNodesWithoutResettingTopology)
{
    ONEFLOW::ScalarGrid grid;
    grid.xn.Resize(2);
    grid.elements.AddElem(std::vector<int>{0, 0});
    grid.eTypes.AddData(ONEFLOW::BAR_2);
    grid.faces.Resize(1);
    grid.lc.AddData(29);

    EXPECT_THROW(grid.CalcTopology(), std::runtime_error);
    ASSERT_EQ(grid.faces.GetNElements(), 1u);
    ASSERT_EQ(grid.lc.GetNElements(), 1u);
    EXPECT_EQ(grid.lc[0], 29);
}
