// tests/PointLocatorLifecycleTest.cpp
#include <gtest/gtest.h>
#include "PointLocator.h"
#include "Grid.h"
#include "NodeMesh.h"
#include <memory>
#include "Constant.h"

// Verify that PointLocator can be safely instantiated and destroyed on the stack.
// The unique_ptr<AdtTree> member should be correctly initialized to nullptr.
TEST(PointLocatorLifecycleTest, DefaultConstructionIsSafe)
{
    EXPECT_NO_THROW({
        ONEFLOW::PointLocator locator;
        // locator goes out of scope here; destructor runs safely
        // even though coorTree is nullptr (unique_ptr handles this correctly)
    });
}

// Verify that Initialize correctly creates the ADT tree and cleanup is automatic.
TEST(PointLocatorLifecycleTest, InitializeAndCleanupIsSafe)
{
    EXPECT_NO_THROW({
        ONEFLOW::PointLocator locator;
        ONEFLOW::RealField pmin(3, 0.0);
        ONEFLOW::RealField pmax(3, 10.0);
        locator.Initialize(pmin, pmax, 1.0e-6);
        
        // Add some points to exercise the tree
        locator.AddPoint(1.0, 1.0, 1.0);
        locator.AddPoint(5.0, 5.0, 5.0);
        locator.AddPoint(9.0, 9.0, 9.0);
        
        // locator goes out of scope here.
        // unique_ptr automatically deletes the AdtTree and all its nodes.
        // No memory leak, no manual delete needed.
    });
}

// Verify that re-initializing does not leak the previous tree.
TEST(PointLocatorLifecycleTest, ReInitializeDoesNotLeak)
{
    EXPECT_NO_THROW({
        ONEFLOW::PointLocator locator;
        ONEFLOW::RealField pmin(3, 0.0);
        ONEFLOW::RealField pmax(3, 10.0);
        
        // First initialization
        locator.Initialize(pmin, pmax, 1.0e-6);
        locator.AddPoint(1.0, 1.0, 1.0);
        
        // Second initialization: unique_ptr's operator= automatically 
        // deletes the old tree before taking ownership of the new one.
        // This was a memory leak in the legacy raw pointer version!
        locator.Initialize(pmin, pmax, 1.0e-6);
        EXPECT_EQ(locator.GetNPoint(), 0);
        EXPECT_EQ(locator.FindPoint(1.0, 1.0, 1.0), ONEFLOW::INVALID_INDEX);
        EXPECT_EQ(locator.AddPoint(2.0, 2.0, 2.0), 0);
        EXPECT_EQ(locator.GetNPoint(), 1);
    });
}
namespace
{
class PointLocatorTestGrid : public ONEFLOW::Grid
{
public:
    std::unique_ptr< ONEFLOW::Grid > Clone() const override
    {
        return std::make_unique< PointLocatorTestGrid >();
    }

    void GetMinMaxDistance( ONEFLOW::Real & minDistance,
                            ONEFLOW::Real & maxDistance ) const override
    {
        minDistance = 0.1;
        maxDistance = 1.0;
    }
};
} // namespace

TEST( PointLocatorLifecycleTest, ExplicitToleranceControlsFinalMatchPredicate )
{
    ONEFLOW::PointLocator locator;
    const ONEFLOW::RealField pmin{ 0.0, 0.0, 0.0 };
    const ONEFLOW::RealField pmax{ 10.0, 10.0, 10.0 };
    locator.Initialize( pmin, pmax, 0.5 );

    EXPECT_EQ( locator.AddPoint( 1.0, 1.0, 1.0 ), 0 );
    EXPECT_EQ( locator.FindPoint( 1.49, 1.0, 1.0 ), 0 );
    EXPECT_EQ( locator.FindPoint( 1.5, 1.0, 1.0 ), ONEFLOW::INVALID_INDEX );
}

TEST(PointLocatorLifecycleTest, SpecialToleranceExpandsAdtBounds)
{
    PointLocatorTestGrid grid;
    grid.Init();
    grid.nodeMesh->SetCoordinates(
        ONEFLOW::RealField{ 0.0, 1.0 },
        ONEFLOW::RealField{ 0.0, 1.0 },
        ONEFLOW::RealField{ 0.0, 1.0 } );

    ONEFLOW::PointLocator locator;
    const double explicitTolerance = 0.5;
    locator.InitializeSpecial( &grid, explicitTolerance );

    // This point is outside the grid's raw bounds but inside the explicit
    // tolerance-expanded ADT bounds.
    EXPECT_EQ( locator.AddPoint( 1.4, 0.5, 0.5 ), 0 );
    EXPECT_EQ( locator.FindPoint( 1.45, 0.5, 0.5 ), 0 );
}
