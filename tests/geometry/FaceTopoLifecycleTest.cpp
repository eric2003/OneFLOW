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
