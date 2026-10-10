#include <gtest/gtest.h>
#include <stdexcept>

#include "ScalarIFace.h"

TEST(ScalarIFaceLifecycleTest, AddInterfaceKeepsMappingsConsistent)
{
    ONEFLOW::ScalarIFace interfaceData;

    interfaceData.AddInterface(100, 2, 7);
    interfaceData.AddInterface(105, 3, 9);

    ASSERT_EQ(interfaceData.iglobalfaces.size(), 2u);
    ASSERT_EQ(interfaceData.zones.size(), 2u);
    ASSERT_EQ(interfaceData.cells.size(), 2u);
    EXPECT_EQ(interfaceData.iglobalfaces[0], 100);
    EXPECT_EQ(interfaceData.zones[0], 2);
    EXPECT_EQ(interfaceData.cells[0], 7);
    EXPECT_EQ(interfaceData.GetLocalInterfaceId(105), 1);
    EXPECT_EQ(interfaceData.local_to_global_interfaces.at(0), 100);
    EXPECT_EQ(interfaceData.local_to_global_interfaces.at(1), 105);
}

TEST(ScalarIFaceLifecycleTest, RejectsDuplicateInterfaceWithoutChangingState)
{
    ONEFLOW::ScalarIFace interfaceData;
    interfaceData.AddInterface(100, 2, 7);

    EXPECT_THROW(interfaceData.AddInterface(100, 4, 8), std::invalid_argument);

    ASSERT_EQ(interfaceData.iglobalfaces.size(), 1u);
    ASSERT_EQ(interfaceData.zones.size(), 1u);
    ASSERT_EQ(interfaceData.cells.size(), 1u);
    EXPECT_EQ(interfaceData.iglobalfaces[0], 100);
    EXPECT_EQ(interfaceData.zones[0], 2);
    EXPECT_EQ(interfaceData.cells[0], 7);
    EXPECT_EQ(interfaceData.global_to_local_interfaces.at(100), 0);
    EXPECT_EQ(interfaceData.local_to_global_interfaces.at(0), 100);
}

TEST(ScalarIFaceLifecycleTest, RejectsInvalidTopologyWithoutReplacingNeighborGroups)
{
    ONEFLOW::ScalarIFace interfaceData;
    interfaceData.AddInterface(30, 3, 300);
    interfaceData.ReconstructNeighbor();
    ASSERT_EQ(interfaceData.data.size(), 1u);

    // Public legacy arrays can be populated outside AddInterface; reject invalid
    // state before replacing the previously reconstructed neighbor groups.
    interfaceData.iglobalfaces.push_back(30);
    interfaceData.zones.push_back(4);
    interfaceData.cells.push_back(400);

    EXPECT_THROW(interfaceData.ReconstructNeighbor(), std::runtime_error);
    ASSERT_EQ(interfaceData.data.size(), 1u);
    EXPECT_EQ(interfaceData.data[0].zonej, 3);
    ASSERT_EQ(interfaceData.data[0].iglobalfaces.size(), 1u);
    EXPECT_EQ(interfaceData.data[0].iglobalfaces[0], 30);
}

TEST(ScalarIFaceLifecycleTest, ReconstructNeighborIsGroupedOrderedAndRepeatable)
{
    ONEFLOW::ScalarIFace interfaceData;
    interfaceData.AddInterface(30, 3, 300);
    interfaceData.AddInterface(20, 2, 200);
    interfaceData.AddInterface(31, 3, 301);

    interfaceData.ReconstructNeighbor();

    ASSERT_EQ(interfaceData.data.size(), 2u);
    EXPECT_EQ(interfaceData.data[0].zonej, 2);
    ASSERT_EQ(interfaceData.data[0].iglobalfaces.size(), 1u);
    EXPECT_EQ(interfaceData.data[0].iglobalfaces[0], 20);
    EXPECT_EQ(interfaceData.data[0].ifaces[0], 1);
    EXPECT_EQ(interfaceData.data[0].cells[0], 200);

    EXPECT_EQ(interfaceData.data[1].zonej, 3);
    ASSERT_EQ(interfaceData.data[1].iglobalfaces.size(), 2u);
    EXPECT_EQ(interfaceData.data[1].iglobalfaces[0], 30);
    EXPECT_EQ(interfaceData.data[1].iglobalfaces[1], 31);
    EXPECT_EQ(interfaceData.data[1].ifaces[0], 0);
    EXPECT_EQ(interfaceData.data[1].ifaces[1], 2);
    EXPECT_EQ(interfaceData.data[1].cells[0], 300);
    EXPECT_EQ(interfaceData.data[1].cells[1], 301);

    interfaceData.ReconstructNeighbor();

    ASSERT_EQ(interfaceData.data.size(), 2u);
    EXPECT_EQ(interfaceData.data[0].iglobalfaces.size(), 1u);
    EXPECT_EQ(interfaceData.data[1].iglobalfaces.size(), 2u);
}

TEST(ScalarIFaceLifecycleTest, RejectsNegativeNeighborZoneWithoutReplacingNeighborGroups)
{
    ONEFLOW::ScalarIFace interfaceData;
    interfaceData.AddInterface(30, 3, 300);
    interfaceData.ReconstructNeighbor();
    ASSERT_EQ(interfaceData.data.size(), 1u);

    interfaceData.zones[0] = -1;

    EXPECT_THROW(interfaceData.ReconstructNeighbor(), std::runtime_error);
    ASSERT_EQ(interfaceData.data.size(), 1u);
    EXPECT_EQ(interfaceData.data[0].zonej, 3);
    ASSERT_EQ(interfaceData.data[0].iglobalfaces.size(), 1u);
    EXPECT_EQ(interfaceData.data[0].iglobalfaces[0], 30);
}
