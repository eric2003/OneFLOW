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


TEST(ScalarIFaceLifecycleTest, ReadOnlyQueriesCanUseConstInterface)
{
    ONEFLOW::ScalarIFace interfaceData;
    interfaceData.AddInterface(100, 2, 7);
    interfaceData.ReconstructNeighbor();

    const ONEFLOW::ScalarIFace & readOnlyInterface = interfaceData;
    EXPECT_EQ(readOnlyInterface.GetNIFaces(), 1);
    EXPECT_EQ(readOnlyInterface.GetLocalInterfaceId(100), 0);
    EXPECT_EQ(readOnlyInterface.FindINeibor(2), 0);
    EXPECT_EQ(readOnlyInterface.FindINeibor(9), -1);
}

TEST(ScalarIFaceLifecycleTest, CalculatesLocalIdsFromReadOnlyGlobalIds)
{
    ONEFLOW::ScalarIFace interfaceData;
    interfaceData.AddInterface(100, 2, 7);
    interfaceData.ReconstructNeighbor();

    const std::vector<int> globalFaces{100};
    std::vector<int> localFaces;
    interfaceData.CalcLocalInterfaceId(2, globalFaces, localFaces);

    ASSERT_EQ(localFaces.size(), 1u);
    EXPECT_EQ(localFaces[0], 0);
    ASSERT_EQ(interfaceData.data[0].recv_ifaces.size(), 1u);
    EXPECT_EQ(interfaceData.data[0].recv_ifaces[0], 0);
}

TEST(ScalarIFaceLifecycleTest, RejectsNonReciprocalMapsBeforeAddingInterface)
{
    ONEFLOW::ScalarIFace interfaceData;
    interfaceData.AddInterface(100, 2, 7);

    // Equal map sizes are not sufficient: the two directions must describe
    // the same local/global interface identity.
    interfaceData.local_to_global_interfaces[0] = 999;

    EXPECT_THROW(interfaceData.AddInterface(105, 3, 9), std::logic_error);

    ASSERT_EQ(interfaceData.iglobalfaces.size(), 1u);
    EXPECT_EQ(interfaceData.iglobalfaces[0], 100);
    EXPECT_EQ(interfaceData.zones.size(), 1u);
    EXPECT_EQ(interfaceData.cells.size(), 1u);
    EXPECT_EQ(interfaceData.global_to_local_interfaces.size(), 1u);
    EXPECT_EQ(interfaceData.local_to_global_interfaces.at(0), 999);
    EXPECT_EQ(interfaceData.global_to_local_interfaces.count(105), 0u);
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

TEST(ScalarIFaceLifecycleTest, RejectsSerializedArraysWhenInterfaceCountIsZero)
{
    ONEFLOW::ScalarIFace interfaceData;
    interfaceData.target_interfaces.push_back(4);

    // Validation must run before writing to the output buffer.
    EXPECT_THROW(interfaceData.WriteInterfaceTopology(nullptr), std::logic_error);
}

TEST(ScalarIFaceLifecycleTest, ReconstructNeighborRejectsCorruptedIdentityMapsWithoutReplacingGroups)
{
    ONEFLOW::ScalarIFace interfaceData;
    interfaceData.AddInterface(100, 2, 7);
    interfaceData.ReconstructNeighbor();
    ASSERT_EQ(interfaceData.data.size(), 1u);

    // Simulate a legacy direct mutation that bypasses AddInterface().
    interfaceData.global_to_local_interfaces[100] = 4;

    EXPECT_THROW(interfaceData.ReconstructNeighbor(), std::logic_error);

    ASSERT_EQ(interfaceData.data.size(), 1u);
    EXPECT_EQ(interfaceData.data[0].zonej, 2);
    ASSERT_EQ(interfaceData.data[0].ifaces.size(), 1u);
    EXPECT_EQ(interfaceData.data[0].ifaces[0], 0);
}

TEST(ScalarIFaceLifecycleTest, CalcLocalInterfaceIdRejectsCorruptedMapsWithoutChangingOutputs)
{
    ONEFLOW::ScalarIFace interfaceData;
    interfaceData.AddInterface(100, 2, 7);
    interfaceData.ReconstructNeighbor();

    std::vector<int> localFaces{77};
    interfaceData.data[0].recv_ifaces = {42};
    interfaceData.global_to_local_interfaces[100] = 4;

    EXPECT_THROW(interfaceData.CalcLocalInterfaceId(2, std::vector<int>{100}, localFaces), std::logic_error);

    ASSERT_EQ(localFaces.size(), 1u);
    EXPECT_EQ(localFaces[0], 77);
    ASSERT_EQ(interfaceData.data[0].recv_ifaces.size(), 1u);
    EXPECT_EQ(interfaceData.data[0].recv_ifaces[0], 42);
}

TEST(ScalarIFaceLifecycleTest, GetLocalInterfaceIdRejectsCorruptedIdentityMaps)
{
    ONEFLOW::ScalarIFace interfaceData;
    interfaceData.AddInterface(100, 2, 7);

    // Public legacy maps can be modified without going through AddInterface().
    interfaceData.local_to_global_interfaces[0] = 999;

    EXPECT_THROW(interfaceData.GetLocalInterfaceId(100), std::logic_error);
}
