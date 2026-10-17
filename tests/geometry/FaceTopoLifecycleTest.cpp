#include <gtest/gtest.h>
#include <memory>
#include <stdexcept>

#include "FaceTopo.h"
#include "Grid.h"
#include "InterFace.h"

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
