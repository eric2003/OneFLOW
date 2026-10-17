#include <gtest/gtest.h>
#include <memory>
#include <stdexcept>

#include "FaceTopo.h"
#include "Grid.h"
#include "InterFace.h"
#include "IFaceLink.h"
#include "GridHandles.h"
#include "Boundary.h"
#include "BcRecord.h"
#include "FaceSearch.h"


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


namespace
{
// Minimal single-zone grids collection required by IFaceLink's constructor.
struct MinimalZoneFixture
{
    ONEFLOW::Grids grids;

    MinimalZoneFixture()
    {
        grids.push_back(std::make_unique<TestGrid>());
        grids[0]->id = 0;
        grids[0]->nCells = 1;
    }

    ONEFLOW::Grid & grid() { return *grids[0]; }
};

void SeedTwoFaceTopology(ONEFLOW::FaceTopo & faceTopo, int boundaryBcType)
{
    // Face 0: boundary, Face 1: internal.
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

TEST(FaceTopoLifecycleTest, NonInterfaceRebuildPipelineCommitsAndClearsStaging)
{
    // Intent: for non-interface boundary faces, SetNewFace2Node/SetNewFace2Cell
    // stage a 1:1 copy of live topology; ModifyFaceNodeId is the combined entry;
    // UpdateOtherTopologyTerm commits staging and leaves no pending rebuild.
    MinimalZoneFixture fixture;
    ONEFLOW::IFaceLink iFaceLink(fixture.grids);
    ONEFLOW::FaceTopo faceTopo;
    faceTopo.BindGrid(fixture.grid());
    SeedTwoFaceTopology(faceTopo, ONEFLOW::BC::SOLID_SURFACE);

    EXPECT_FALSE(faceTopo.HasPendingTopologyRebuild());

    faceTopo.ModifyFaceNodeId(iFaceLink);
    EXPECT_TRUE(faceTopo.HasPendingTopologyRebuild());

    // Before commit, live topology is still the seeded content.
    EXPECT_EQ(faceTopo.GetFaceCount(), static_cast<ONEFLOW::HXSize_t>(2));
    EXPECT_EQ(faceTopo.GetFaces()[0].size(), static_cast<std::size_t>(2));
    EXPECT_EQ(faceTopo.GetLeftCells()[0], 0);

    faceTopo.UpdateOtherTopologyTerm();
    EXPECT_FALSE(faceTopo.HasPendingTopologyRebuild());

    // After commit: faces/cells come from staging (1:1 copy).
    EXPECT_EQ(faceTopo.GetFaces().size(), static_cast<std::size_t>(2));
    EXPECT_EQ(faceTopo.GetFaces()[0][0], 0);
    EXPECT_EQ(faceTopo.GetFaces()[0][1], 1);
    EXPECT_EQ(faceTopo.GetFaces()[1][0], 1);
    EXPECT_EQ(faceTopo.GetFaces()[1][1], 2);

    // Boundary face right cell is rewritten as iFace + nCells.
    EXPECT_EQ(faceTopo.GetLeftCells().size(), static_cast<std::size_t>(2));
    EXPECT_EQ(faceTopo.GetRightCells().size(), static_cast<std::size_t>(2));
    EXPECT_EQ(faceTopo.GetLeftCells()[0], 0);
    EXPECT_EQ(faceTopo.GetRightCells()[0], 0 + fixture.grid().nCells);
    EXPECT_EQ(faceTopo.GetLeftCells()[1], 0);
    EXPECT_EQ(faceTopo.GetRightCells()[1], 1);
}

TEST(FaceTopoLifecycleTest, InterfaceFaceWithEmptyChildrenStillRebuildsOneToOne)
{
    // Intent: interface BC with zero child faces still goes through the
    // IFaceLink lookup path, then stages a 1:1 copy. The minimal IFaceLink
    // fixture only needs l2g and an empty cFaceId slot.
    MinimalZoneFixture fixture;
    ONEFLOW::IFaceLink iFaceLink(fixture.grids);

    // local interface face 0 -> global face 0; no children.
    iFaceLink.l2g[0].push_back(0);
    iFaceLink.face_search->cFaceId.resize(1); // cFaceId[0] empty

    ONEFLOW::FaceTopo faceTopo;
    faceTopo.BindGrid(fixture.grid());
    SeedTwoFaceTopology(faceTopo, ONEFLOW::BC::INTERFACE);

    faceTopo.ModifyFaceNodeId(iFaceLink);
    EXPECT_TRUE(faceTopo.HasPendingTopologyRebuild());

    faceTopo.UpdateOtherTopologyTerm();
    EXPECT_FALSE(faceTopo.HasPendingTopologyRebuild());

    EXPECT_EQ(faceTopo.GetFaces().size(), static_cast<std::size_t>(2));
    EXPECT_EQ(faceTopo.GetFaces()[0][0], 0);
    EXPECT_EQ(faceTopo.GetFaces()[0][1], 1);
    EXPECT_EQ(faceTopo.GetLeftCells()[0], 0);
    EXPECT_EQ(faceTopo.GetRightCells()[0], 0 + fixture.grid().nCells);
}

TEST(FaceTopoLifecycleTest, InterfaceRebuildRejectsMissingFaceSearchMapping)
{
    // Intent: interface BC requires a valid IFaceLink mapping; a default
    // IFaceLink without l2g entries must fail before mutating live topology.
    MinimalZoneFixture fixture;
    ONEFLOW::IFaceLink iFaceLink(fixture.grids);

    ONEFLOW::FaceTopo faceTopo;
    faceTopo.BindGrid(fixture.grid());
    SeedTwoFaceTopology(faceTopo, ONEFLOW::BC::INTERFACE);

    EXPECT_THROW(faceTopo.SetNewFace2Node(iFaceLink), std::out_of_range);
    EXPECT_FALSE(faceTopo.HasPendingTopologyRebuild());
}
