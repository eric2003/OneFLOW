#include <gtest/gtest.h>
#include <memory>
#include <stdexcept>

#include "FaceTopo.h"
#include "BcRecord.h"
#include "Grid.h"
#include "InterFace.h"
#include "IFaceLink.h"
#include "GridHandles.h"
#include "Boundary.h"
#include "FaceSearch.h"
#include "NodeMesh.h"

namespace
{
class TestGrid final : public ONEFLOW::Grid
{
public:
    std::unique_ptr< ONEFLOW::Grid > Clone() const override
    {
        return std::make_unique< TestGrid >();
    }
};

struct MinimalZoneFixture
{
    ONEFLOW::Grids grids;

    MinimalZoneFixture()
    {
        grids.push_back(std::make_unique<TestGrid>());
        grids[0]->id = 0;
        grids[0]->nCells = 1;
        // Intentionally no nodeMesh: IFaceLink must construct without building
        // the PointLocator ADT so topology-only rebuild tests stay lightweight.
    }

    ONEFLOW::Grid & grid() { return *grids[0]; }
};

// Boundary BC type values match Boundary.cpp assignments.
// Any negative type is treated as interface by BC::IsInterfaceBc.
constexpr int kSolidSurfaceBc = 2;
constexpr int kInterfaceBc = -1;

void SeedTwoFaceTopology(ONEFLOW::FaceTopo & faceTopo, int boundaryBcType)
{
    faceTopo.GetFaceTypes().clear();
    faceTopo.GetFaces().clear();
    faceTopo.GetLeftCells().clear();
    faceTopo.GetRightCells().clear();

    faceTopo.GetFaceTypes().push_back(2);
    faceTopo.GetFaceTypes().push_back(2);
    faceTopo.GetFaces().push_back(ONEFLOW::IntField{0, 1});
    faceTopo.GetFaces().push_back(ONEFLOW::IntField{1, 2});
    faceTopo.GetLeftCells().push_back(0);
    faceTopo.GetRightCells().push_back(-1);
    faceTopo.GetLeftCells().push_back(0);
    faceTopo.GetRightCells().push_back(1);

    faceTopo.SetBoundaryFaceCount(1);
    faceTopo.GetBcRecord().bcType[0] = boundaryBcType;
    faceTopo.GetBcRecord().bcNameId[0] = 0;
}
}

TEST(FaceTopoLifecycleTest, GetGridRejectsUnboundGrid)
{
    ONEFLOW::FaceTopo faceTopo;
    const ONEFLOW::FaceTopo & constFaceTopo = faceTopo;

    EXPECT_THROW(faceTopo.GetGrid(), std::logic_error);
    EXPECT_THROW(constFaceTopo.GetGrid(), std::logic_error);
    int value = 0;
    EXPECT_THROW(faceTopo.GetSId(0, 0, value), std::logic_error);
    EXPECT_THROW(faceTopo.GetTId(0, 0, value), std::logic_error);
}

TEST(FaceTopoLifecycleTest, SideAccessRejectsMissingInterfaceData)
{
    TestGrid grid;
    ONEFLOW::FaceTopo faceTopo;
    faceTopo.BindGrid(grid);

    int value = 0;
    EXPECT_THROW(faceTopo.GetSId(0, 0, value), std::logic_error);
    EXPECT_THROW(faceTopo.GetTId(0, 0, value), std::logic_error);
}

TEST(FaceTopoLifecycleTest, SideAccessRejectsInvalidInterfaceFaceId)
{
    TestGrid grid;
    grid.interFace = std::make_unique< ONEFLOW::InterFace >();

    ONEFLOW::FaceTopo faceTopo;
    faceTopo.BindGrid(grid);

    int value = 0;
    EXPECT_THROW(faceTopo.GetSId(0, 0, value), std::out_of_range);
    EXPECT_THROW(faceTopo.GetTId(0, 0, value), std::out_of_range);
}

TEST(FaceTopoLifecycleTest, SideAccessUsesInterfaceToBoundaryMapping)
{
    TestGrid grid;
    grid.interFace = std::make_unique< ONEFLOW::InterFace >();
    grid.interFace->i2b.push_back(1);

    ONEFLOW::FaceTopo faceTopo;
    faceTopo.BindGrid(grid);
    faceTopo.GetLeftCells().push_back(10);
    faceTopo.GetLeftCells().push_back(20);
    faceTopo.GetRightCells().push_back(30);
    faceTopo.GetRightCells().push_back(40);

    int sourceCell = -1;
    int targetCell = -1;
    EXPECT_TRUE(faceTopo.GetSId(0, 0, sourceCell));
    EXPECT_TRUE(faceTopo.GetTId(0, 0, targetCell));
    EXPECT_EQ(sourceCell, 20);
    EXPECT_EQ(targetCell, 40);
}

TEST(FaceTopoLifecycleTest, SideAccessRejectsInvalidBoundaryFaceMapping)
{
    TestGrid grid;
    grid.interFace = std::make_unique< ONEFLOW::InterFace >();
    grid.interFace->i2b.push_back(2);

    ONEFLOW::FaceTopo faceTopo;
    faceTopo.BindGrid(grid);
    faceTopo.GetLeftCells().push_back(10);
    faceTopo.GetLeftCells().push_back(20);
    faceTopo.GetRightCells().push_back(30);
    faceTopo.GetRightCells().push_back(40);

    int value = -1;
    EXPECT_THROW(faceTopo.GetSId(0, 0, value), std::out_of_range);
    EXPECT_THROW(faceTopo.GetTId(0, 0, value), std::out_of_range);
}

TEST(FaceTopoLifecycleTest, FaceCountApisMatchCompatibilityWrappers)
{
    // Intent: GetFaceCount / GetBoundaryFaceCount are preferred; GetNFaces /
    // GetNBFaces / SetNBFaces remain pure aliases so existing call sites stay
    // valid while new code can use clearer names.
    ONEFLOW::FaceTopo faceTopo;

    faceTopo.GetFaceTypes().resize(5);
    EXPECT_EQ(faceTopo.GetFaceCount(), static_cast<ONEFLOW::HXSize_t>(5));
    EXPECT_EQ(faceTopo.GetNFaces(), faceTopo.GetFaceCount());

    // Boundary face count is owned by BcRecord; setting either API must agree.
    faceTopo.SetBoundaryFaceCount(2);
    EXPECT_EQ(faceTopo.GetBoundaryFaceCount(), static_cast<ONEFLOW::HXSize_t>(2));
    EXPECT_EQ(faceTopo.GetNBFaces(), faceTopo.GetBoundaryFaceCount());

    faceTopo.SetNBFaces(3);
    EXPECT_EQ(faceTopo.GetBoundaryFaceCount(), static_cast<ONEFLOW::HXSize_t>(3));
    EXPECT_EQ(faceTopo.GetNBFaces(), faceTopo.GetBoundaryFaceCount());
}

TEST(FaceTopoLifecycleTest, FreshFaceTopoHasNoPendingRebuild)
{
    // Intent: rebuild staging buffers start empty.
    ONEFLOW::FaceTopo faceTopo;
    EXPECT_FALSE(faceTopo.HasPendingTopologyRebuild());
}

TEST(FaceTopoLifecycleTest, UpdateOtherTopologyTermClearsRebuildBuffers)
{
    // Intent: UpdateOtherTopologyTerm always leaves staging buffers empty after
    // commit. Callers must fill rebuildFaces / rebuildLeftCells /
    // rebuildRightCells (via SetNewFace2Node / SetNewFace2Cell) before commit;
    // committing with empty staging replaces live topology with empty containers.
    TestGrid grid;
    ONEFLOW::FaceTopo faceTopo;
    faceTopo.BindGrid(grid);

    faceTopo.GetFaceTypes().push_back(2);
    faceTopo.GetFaces().push_back(ONEFLOW::IntField{0, 1});
    faceTopo.GetLeftCells().push_back(0);
    faceTopo.GetRightCells().push_back(1);

    EXPECT_FALSE(faceTopo.HasPendingTopologyRebuild());

    // Commit with empty staging (no prior SetNewFace2* call).
    faceTopo.UpdateOtherTopologyTerm();

    EXPECT_FALSE(faceTopo.HasPendingTopologyRebuild());
    EXPECT_TRUE(faceTopo.GetFaces().empty());
    EXPECT_TRUE(faceTopo.GetLeftCells().empty());
    EXPECT_TRUE(faceTopo.GetRightCells().empty());
}

TEST(FaceTopoLifecycleTest, NonInterfaceRebuildPipelineCommitsAndClearsStaging)
{
    // Intent: for non-interface boundary faces, SetNewFace2Node/SetNewFace2Cell
    // stage a 1:1 copy of live topology; UpdateOtherTopologyTerm commits staging
    // and leaves no pending rebuild. Boundary right cells become iFace + nCells.
    MinimalZoneFixture fixture;
    ONEFLOW::IFaceLink iFaceLink(fixture.grids);
    ONEFLOW::FaceTopo faceTopo;
    faceTopo.BindGrid(fixture.grid());
    SeedTwoFaceTopology(faceTopo, kSolidSurfaceBc);

    ASSERT_EQ(faceTopo.GetFaceCount(), static_cast<ONEFLOW::HXSize_t>(2));
    ASSERT_EQ(faceTopo.GetBoundaryFaceCount(), static_cast<ONEFLOW::HXSize_t>(1));
    ASSERT_FALSE(faceTopo.HasPendingTopologyRebuild());

    ASSERT_NO_THROW(faceTopo.SetNewFace2Node(iFaceLink));
    ASSERT_NO_THROW(faceTopo.SetNewFace2Cell(iFaceLink));
    EXPECT_TRUE(faceTopo.HasPendingTopologyRebuild());

    EXPECT_EQ(faceTopo.GetFaces().size(), static_cast<std::size_t>(2));
    EXPECT_EQ(faceTopo.GetLeftCells()[0], 0);
    EXPECT_EQ(faceTopo.GetRightCells()[0], -1);

    ASSERT_NO_THROW(faceTopo.UpdateOtherTopologyTerm());
    EXPECT_FALSE(faceTopo.HasPendingTopologyRebuild());

    ASSERT_EQ(faceTopo.GetFaces().size(), static_cast<std::size_t>(2));
    EXPECT_EQ(faceTopo.GetFaces()[0][0], 0);
    EXPECT_EQ(faceTopo.GetFaces()[0][1], 1);
    EXPECT_EQ(faceTopo.GetFaces()[1][0], 1);
    EXPECT_EQ(faceTopo.GetFaces()[1][1], 2);

    ASSERT_EQ(faceTopo.GetLeftCells().size(), static_cast<std::size_t>(2));
    ASSERT_EQ(faceTopo.GetRightCells().size(), static_cast<std::size_t>(2));
    EXPECT_EQ(faceTopo.GetLeftCells()[0], 0);
    EXPECT_EQ(faceTopo.GetRightCells()[0], fixture.grid().nCells);
    EXPECT_EQ(faceTopo.GetLeftCells()[1], 0);
    EXPECT_EQ(faceTopo.GetRightCells()[1], 1);
}

TEST(FaceTopoLifecycleTest, InterfaceFaceWithEmptyChildrenStillRebuildsOneToOne)
{
    // Intent: interface BC with zero child faces still goes through the
    // IFaceLink lookup path, then stages a 1:1 copy.
    MinimalZoneFixture fixture;
    ONEFLOW::IFaceLink iFaceLink(fixture.grids);

    ASSERT_EQ(static_cast<int>(iFaceLink.l2g.size()), 1);
    iFaceLink.l2g[0].push_back(0);
    iFaceLink.face_search->interfaceSplit.childFaceIdsOfParent.resize(1);

    ONEFLOW::FaceTopo faceTopo;
    faceTopo.BindGrid(fixture.grid());
    SeedTwoFaceTopology(faceTopo, kInterfaceBc);

    ASSERT_NO_THROW(faceTopo.ModifyFaceNodeId(iFaceLink));
    EXPECT_TRUE(faceTopo.HasPendingTopologyRebuild());

    ASSERT_NO_THROW(faceTopo.UpdateOtherTopologyTerm());
    EXPECT_FALSE(faceTopo.HasPendingTopologyRebuild());

    ASSERT_EQ(faceTopo.GetFaces().size(), static_cast<std::size_t>(2));
    EXPECT_EQ(faceTopo.GetFaces()[0][0], 0);
    EXPECT_EQ(faceTopo.GetFaces()[0][1], 1);
    EXPECT_EQ(faceTopo.GetLeftCells()[0], 0);
    EXPECT_EQ(faceTopo.GetRightCells()[0], fixture.grid().nCells);
}

TEST(FaceTopoLifecycleTest, InterfaceRebuildRejectsMissingFaceSearchMapping)
{
    // Intent: interface BC requires a valid IFaceLink mapping; a default
    // IFaceLink without l2g entries must fail before leaving a pending rebuild.
    MinimalZoneFixture fixture;
    ONEFLOW::IFaceLink iFaceLink(fixture.grids);

    ONEFLOW::FaceTopo faceTopo;
    faceTopo.BindGrid(fixture.grid());
    SeedTwoFaceTopology(faceTopo, kInterfaceBc);

    EXPECT_THROW(faceTopo.SetNewFace2Node(iFaceLink), std::out_of_range);
    EXPECT_FALSE(faceTopo.HasPendingTopologyRebuild());
}

TEST(FaceTopoLifecycleTest, IFaceLinkConstructsWithoutNodeMesh)
{
    // Intent: constructing IFaceLink must not require nodeMesh / PointLocator ADT.
    // Topology rebuild paths only need face_search and l2g.
    MinimalZoneFixture fixture;
    ASSERT_EQ(fixture.grid().nodeMesh, nullptr);

    ONEFLOW::IFaceLink iFaceLink(fixture.grids);
    EXPECT_NE(iFaceLink.face_search, nullptr);
    EXPECT_NE(iFaceLink.point_search, nullptr);
    EXPECT_FALSE(iFaceLink.IsPointLocatorReady());
    EXPECT_EQ(static_cast<int>(iFaceLink.l2g.size()), 1);
}

TEST(FaceTopoLifecycleTest, ResolveChildFaceNodeIds_RelativeFlagsMapThroughParent)
{
    // Intent: flag==1 means childNodeIds[i] indexes into the parent face node list.
    const ONEFLOW::IntField parent{10, 20, 30, 40};
    const ONEFLOW::IntField childIds{0, 2, 3};
    const ONEFLOW::IntField flags{1, 1, 1};

    const ONEFLOW::IntField resolved = ONEFLOW::ResolveChildFaceNodeIds(
        parent, childIds, flags, false, 0);

    ASSERT_EQ(resolved.size(), static_cast<std::size_t>(3));
    EXPECT_EQ(resolved[0], 10);
    EXPECT_EQ(resolved[1], 30);
    EXPECT_EQ(resolved[2], 40);
}

TEST(FaceTopoLifecycleTest, ResolveChildFaceNodeIds_AbsoluteFlagsUseMeshIds)
{
    // Intent: non-relative flags pass absolute mesh node IDs through unchanged.
    const ONEFLOW::IntField parent{10, 20};
    const ONEFLOW::IntField childIds{5, 7};
    const ONEFLOW::IntField flags{0, 0};

    const ONEFLOW::IntField resolved = ONEFLOW::ResolveChildFaceNodeIds(
        parent, childIds, flags, true, 10);

    ASSERT_EQ(resolved.size(), static_cast<std::size_t>(2));
    EXPECT_EQ(resolved[0], 5);
    EXPECT_EQ(resolved[1], 7);
}

TEST(FaceTopoLifecycleTest, ResolveChildFaceNodeIds_MixedFlags)
{
    const ONEFLOW::IntField parent{100, 200, 300};
    const ONEFLOW::IntField childIds{1, 9, 0};
    const ONEFLOW::IntField flags{1, 0, 1};

    const ONEFLOW::IntField resolved = ONEFLOW::ResolveChildFaceNodeIds(
        parent, childIds, flags, true, 20);

    ASSERT_EQ(resolved.size(), static_cast<std::size_t>(3));
    EXPECT_EQ(resolved[0], 200);
    EXPECT_EQ(resolved[1], 9);
    EXPECT_EQ(resolved[2], 100);
}

TEST(FaceTopoLifecycleTest, ResolveChildFaceNodeIds_RejectsMismatchedSizes)
{
    const ONEFLOW::IntField parent{1, 2};
    const ONEFLOW::IntField childIds{0, 1};
    const ONEFLOW::IntField flags{1};

    EXPECT_THROW(
        ONEFLOW::ResolveChildFaceNodeIds(parent, childIds, flags, false, 0),
        std::logic_error);
}

TEST(FaceTopoLifecycleTest, ResolveChildFaceNodeIds_RejectsOutOfRangeRelative)
{
    const ONEFLOW::IntField parent{1, 2};
    const ONEFLOW::IntField childIds{2};
    const ONEFLOW::IntField flags{1};

    EXPECT_THROW(
        ONEFLOW::ResolveChildFaceNodeIds(parent, childIds, flags, false, 0),
        std::out_of_range);
}

TEST(FaceTopoLifecycleTest, ResolveChildFaceNodeIds_RejectsAbsoluteWithoutMesh)
{
    const ONEFLOW::IntField parent{1, 2};
    const ONEFLOW::IntField childIds{0};
    const ONEFLOW::IntField flags{0};

    EXPECT_THROW(
        ONEFLOW::ResolveChildFaceNodeIds(parent, childIds, flags, false, 0),
        std::out_of_range);
}

TEST(FaceTopoLifecycleTest, InterfaceParentFaceSplitsIntoTwoRelativeChildren)
{
    // Intent: one interface parent face with two relative children becomes two
    // staged faces whose nodes are resolved through the parent face node list.
    // This is the minimal multi-child rebuild path enabled by ResolveChildFaceNodeIds.
    MinimalZoneFixture fixture;
    ONEFLOW::IFaceLink iFaceLink(fixture.grids);

    // local interface face 0 -> global face 0; two children 0 and 1.
    iFaceLink.l2g[0].push_back(0);
    iFaceLink.face_search->interfaceSplit.childFaceIdsOfParent.resize(1);
    iFaceLink.face_search->interfaceSplit.childFaceIdsOfParent[0] = ONEFLOW::IntField{0, 1};
    iFaceLink.face_search->interfaceSplit.childFaceNodeIds.resize(2);
    iFaceLink.face_search->interfaceSplit.childFaceNodeFlags.resize(2);
    // child 0: parent local corners 0 and 1
    iFaceLink.face_search->interfaceSplit.childFaceNodeIds[0] = ONEFLOW::IntField{0, 1};
    iFaceLink.face_search->interfaceSplit.childFaceNodeFlags[0] = ONEFLOW::IntField{1, 1};
    // child 1: parent local corners 2 and 3
    iFaceLink.face_search->interfaceSplit.childFaceNodeIds[1] = ONEFLOW::IntField{2, 3};
    iFaceLink.face_search->interfaceSplit.childFaceNodeFlags[1] = ONEFLOW::IntField{1, 1};

    ONEFLOW::FaceTopo faceTopo;
    faceTopo.BindGrid(fixture.grid());

    // Single interface boundary face; parent nodes 10,20,30,40.
    faceTopo.GetFaceTypes().push_back(2);
    faceTopo.GetFaces().push_back(ONEFLOW::IntField{10, 20, 30, 40});
    faceTopo.GetLeftCells().push_back(7);
    faceTopo.GetRightCells().push_back(-1);
    faceTopo.SetBoundaryFaceCount(1);
    faceTopo.GetBcRecord().bcType[0] = kInterfaceBc;
    faceTopo.GetBcRecord().bcNameId[0] = 0;

    ASSERT_NO_THROW(faceTopo.ModifyFaceNodeId(iFaceLink));
    EXPECT_TRUE(faceTopo.HasPendingTopologyRebuild());

    ASSERT_NO_THROW(faceTopo.UpdateOtherTopologyTerm());
    EXPECT_FALSE(faceTopo.HasPendingTopologyRebuild());

    // Two faces after split; no internal faces were present.
    ASSERT_EQ(faceTopo.GetFaces().size(), static_cast<std::size_t>(2));
    ASSERT_EQ(faceTopo.GetFaces()[0].size(), static_cast<std::size_t>(2));
    ASSERT_EQ(faceTopo.GetFaces()[1].size(), static_cast<std::size_t>(2));
    EXPECT_EQ(faceTopo.GetFaces()[0][0], 10);
    EXPECT_EQ(faceTopo.GetFaces()[0][1], 20);
    EXPECT_EQ(faceTopo.GetFaces()[1][0], 30);
    EXPECT_EQ(faceTopo.GetFaces()[1][1], 40);

    // Left cell of the parent is duplicated onto each child; right cells are
    // iFace + nCells for the new boundary faces.
    ASSERT_EQ(faceTopo.GetLeftCells().size(), static_cast<std::size_t>(2));
    ASSERT_EQ(faceTopo.GetRightCells().size(), static_cast<std::size_t>(2));
    EXPECT_EQ(faceTopo.GetLeftCells()[0], 7);
    EXPECT_EQ(faceTopo.GetLeftCells()[1], 7);
    EXPECT_EQ(faceTopo.GetRightCells()[0], 0 + fixture.grid().nCells);
    EXPECT_EQ(faceTopo.GetRightCells()[1], 1 + fixture.grid().nCells);
}

TEST(FaceTopoLifecycleTest, GetInterfaceChildFaceIds_ReturnsEmptyChildren)
{
    // Intent: valid l2g + empty childFaceIdsOfParent slot returns an empty child list.
    MinimalZoneFixture fixture;
    ONEFLOW::IFaceLink iFaceLink(fixture.grids);
    iFaceLink.l2g[0].push_back(0);
    iFaceLink.face_search->interfaceSplit.childFaceIdsOfParent.resize(1);

    const ONEFLOW::IntField & children =
        ONEFLOW::GetInterfaceChildFaceIds(iFaceLink, 0, 0);
    EXPECT_TRUE(children.empty());
}

TEST(FaceTopoLifecycleTest, GetInterfaceChildFaceIds_RejectsMissingLocalMapping)
{
    MinimalZoneFixture fixture;
    ONEFLOW::IFaceLink iFaceLink(fixture.grids);

    EXPECT_THROW(
        ONEFLOW::GetInterfaceChildFaceIds(iFaceLink, 0, 0),
        std::out_of_range);
}

TEST(FaceTopoLifecycleTest, GetInterfaceChildFaceIds_RejectsMismatchedChildNodeTables)
{
    MinimalZoneFixture fixture;
    ONEFLOW::IFaceLink iFaceLink(fixture.grids);
    iFaceLink.l2g[0].push_back(0);
    iFaceLink.face_search->interfaceSplit.childFaceIdsOfParent.resize(1);
    iFaceLink.face_search->interfaceSplit.childFaceIdsOfParent[0] = ONEFLOW::IntField{0};
    iFaceLink.face_search->interfaceSplit.childFaceNodeIds.resize(1);
    iFaceLink.face_search->interfaceSplit.childFaceNodeFlags.resize(1);
    iFaceLink.face_search->interfaceSplit.childFaceNodeIds[0] = ONEFLOW::IntField{0, 1};
    iFaceLink.face_search->interfaceSplit.childFaceNodeFlags[0] = ONEFLOW::IntField{1}; // size mismatch

    EXPECT_THROW(
        ONEFLOW::GetInterfaceChildFaceIds(iFaceLink, 0, 0),
        std::logic_error);
}

TEST(FaceTopoLifecycleTest, InterfaceChildWithAbsoluteNodeIdsNeedsNodeMesh)
{
    // Intent: absolute child node IDs (flag != 1) require a nodeMesh; without
    // it SetNewFace2Node fails before commit. With a mesh, absolute IDs pass
    // through unchanged.
    MinimalZoneFixture fixture;
    fixture.grid().nodeMesh = std::make_unique<ONEFLOW::NodeMesh>();
    fixture.grid().nodeMesh->AddPoint(0.0, 0.0, 0.0); // id 0
    fixture.grid().nodeMesh->AddPoint(1.0, 0.0, 0.0); // id 1
    fixture.grid().nodeMesh->AddPoint(1.0, 1.0, 0.0); // id 2
    fixture.grid().nodeMesh->CalcMinMaxBox();

    ONEFLOW::IFaceLink iFaceLink(fixture.grids);
    iFaceLink.l2g[0].push_back(0);
    iFaceLink.face_search->interfaceSplit.childFaceIdsOfParent.resize(1);
    iFaceLink.face_search->interfaceSplit.childFaceIdsOfParent[0] = ONEFLOW::IntField{0};
    iFaceLink.face_search->interfaceSplit.childFaceNodeIds.resize(1);
    iFaceLink.face_search->interfaceSplit.childFaceNodeFlags.resize(1);
    // Absolute mesh node IDs 1 and 2 (not relative to parent).
    iFaceLink.face_search->interfaceSplit.childFaceNodeIds[0] = ONEFLOW::IntField{1, 2};
    iFaceLink.face_search->interfaceSplit.childFaceNodeFlags[0] = ONEFLOW::IntField{0, 0};

    ONEFLOW::FaceTopo faceTopo;
    faceTopo.BindGrid(fixture.grid());
    faceTopo.GetFaceTypes().push_back(2);
    faceTopo.GetFaces().push_back(ONEFLOW::IntField{10, 20});
    faceTopo.GetLeftCells().push_back(3);
    faceTopo.GetRightCells().push_back(-1);
    faceTopo.SetBoundaryFaceCount(1);
    faceTopo.GetBcRecord().bcType[0] = kInterfaceBc;
    faceTopo.GetBcRecord().bcNameId[0] = 0;

    ASSERT_NO_THROW(faceTopo.ModifyFaceNodeId(iFaceLink));
    ASSERT_NO_THROW(faceTopo.UpdateOtherTopologyTerm());

    ASSERT_EQ(faceTopo.GetFaces().size(), static_cast<std::size_t>(1));
    EXPECT_EQ(faceTopo.GetFaces()[0][0], 1);
    EXPECT_EQ(faceTopo.GetFaces()[0][1], 2);
    EXPECT_EQ(faceTopo.GetLeftCells()[0], 3);
}

TEST(FaceTopoLifecycleTest, InterfaceChildAbsoluteIdsRejectMissingNodeMesh)
{
    MinimalZoneFixture fixture;
    ASSERT_EQ(fixture.grid().nodeMesh, nullptr);

    ONEFLOW::IFaceLink iFaceLink(fixture.grids);
    iFaceLink.l2g[0].push_back(0);
    iFaceLink.face_search->interfaceSplit.childFaceIdsOfParent.resize(1);
    iFaceLink.face_search->interfaceSplit.childFaceIdsOfParent[0] = ONEFLOW::IntField{0};
    iFaceLink.face_search->interfaceSplit.childFaceNodeIds.resize(1);
    iFaceLink.face_search->interfaceSplit.childFaceNodeFlags.resize(1);
    iFaceLink.face_search->interfaceSplit.childFaceNodeIds[0] = ONEFLOW::IntField{1};
    iFaceLink.face_search->interfaceSplit.childFaceNodeFlags[0] = ONEFLOW::IntField{0};

    ONEFLOW::FaceTopo faceTopo;
    faceTopo.BindGrid(fixture.grid());
    faceTopo.GetFaceTypes().push_back(2);
    faceTopo.GetFaces().push_back(ONEFLOW::IntField{10, 20});
    faceTopo.GetLeftCells().push_back(0);
    faceTopo.GetRightCells().push_back(-1);
    faceTopo.SetBoundaryFaceCount(1);
    faceTopo.GetBcRecord().bcType[0] = kInterfaceBc;
    faceTopo.GetBcRecord().bcNameId[0] = 0;

    EXPECT_THROW(faceTopo.SetNewFace2Node(iFaceLink), std::out_of_range);
    EXPECT_FALSE(faceTopo.HasPendingTopologyRebuild());
}
