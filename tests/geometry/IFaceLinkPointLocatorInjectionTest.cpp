// tests/geometry/IFaceLinkPointLocatorInjectionTest.cpp
// Locks PointLocator injection on IFaceLink without grids / ADT build.

#include <gtest/gtest.h>
#include <memory>
#include <stdexcept>
#include <vector>

#include "IFaceLink.h"
#include "FaceSearch.h"
#include "PointLocator.h"
#include "GridHandles.h"
#include "HXDefine.h"

using namespace ONEFLOW;

namespace {

// Minimal double: identity map of synthetic points, no ADT, no nodeMesh.
class RecordingPointLocator : public PointLocator
{
public:
    int FindPoint( Real xm, Real ym, Real zm ) override
    {
        ++ findCalls;
        for ( std::size_t i = 0; i < xs.size(); ++ i )
        {
            if ( xs[ i ] == xm && ys[ i ] == ym && zs[ i ] == zm )
            {
                return static_cast< int >( i );
            }
        }
        return -1;
    }

    int AddPoint( Real xm, Real ym, Real zm ) override
    {
        ++ addCalls;
        const int id = static_cast< int >( xs.size() );
        xs.push_back( xm );
        ys.push_back( ym );
        zs.push_back( zm );
        return id;
    }

    void GetPoint( int id, Real & xm, Real & ym, Real & zm ) const override
    {
        xm = xs.at( static_cast< std::size_t >( id ) );
        ym = ys.at( static_cast< std::size_t >( id ) );
        zm = zs.at( static_cast< std::size_t >( id ) );
    }

    int GetNPoint() const override { return static_cast< int >( xs.size() ); }

    int findCalls = 0;
    int addCalls = 0;
    std::vector< Real > xs, ys, zs;
};

} // namespace

TEST( IFaceLinkPointLocatorInjectionTest, ReplaceRejectsNull )
{
    Grids grids;
    IFaceLink link( grids );
    EXPECT_THROW(
        link.ReplacePointLocator( std::unique_ptr< PointLocator >{}, true ),
        std::invalid_argument );
}

TEST( IFaceLinkPointLocatorInjectionTest, InjectedDoubleIsReadyWithoutGrids )
{
    Grids grids;
    IFaceLink link( grids );
    EXPECT_FALSE( link.IsPointLocatorReady() );

    auto fake = std::make_unique< RecordingPointLocator >();
    RecordingPointLocator * raw = fake.get();
    link.ReplacePointLocator( std::move( fake ), /*markReady=*/true );

    EXPECT_TRUE( link.IsPointLocatorReady() );
    // Ready path must not rebuild ADT from empty grids.
    EXPECT_NO_THROW( link.EnsurePointLocatorReady() );

    const int id = link.GetPointLocator().AddPoint( 1.0, 2.0, 3.0 );
    EXPECT_EQ( id, 0 );
    EXPECT_EQ( raw->addCalls, 1 );
    EXPECT_EQ( link.GetPointLocator().FindPoint( 1.0, 2.0, 3.0 ), 0 );
    EXPECT_EQ( raw->findCalls, 1 );
    const PointLocator & readOnlyLocator = link.GetPointLocator();
    EXPECT_EQ( readOnlyLocator.GetNPoint(), 1 );
    EXPECT_DOUBLE_EQ( readOnlyLocator.GetTol(), 1.0e-6 );

    Real x = 0, y = 0, z = 0;
    readOnlyLocator.GetPoint( 0, x, y, z );
    EXPECT_EQ( x, 1.0 );
    EXPECT_EQ( y, 2.0 );
    EXPECT_EQ( z, 3.0 );
}

TEST( IFaceLinkPointLocatorInjectionTest, MarkReadyFalseKeepsLocatorUnready )
{
    Grids grids;
    IFaceLink link( grids );
    auto fake = std::make_unique< RecordingPointLocator >();
    RecordingPointLocator * raw = fake.get();
    link.ReplacePointLocator( std::move( fake ), /*markReady=*/false );
    EXPECT_FALSE( link.IsPointLocatorReady() );
    // Object is still reachable; readiness is independent of injection.
    EXPECT_EQ( &link.GetPointLocator(), static_cast< PointLocator * >( raw ) );
}

TEST( IFaceLinkPointLocatorInjectionTest, SplitLineWithoutRegisteredChildrenDoesNotPublishMapping )
{
    Grids grids;
    IFaceLink link( grids );

    auto fake = std::make_unique< RecordingPointLocator >();
    fake->AddPoint( 0.0, 0.0, 0.0 );
    fake->AddPoint( 2.0, 0.0, 0.0 );
    fake->AddPoint( 1.0, 0.0, 0.0 );
    link.ReplacePointLocator( std::move( fake ), /*markReady=*/true );

    link.EnsureInterfaceZoneCount( 1 );
    link.EnsureLocalInterfaceFaceCount( 0, 1 );
    link.RegisterLocalInterfaceFace( 0, /*isNew=*/true, 0, 0 );

    FaceSearch faceSearch;
    EXPECT_EQ( faceSearch.AddFace( IntField{ 0, 1 } ), 0 );
    faceSearch.interfaceSplit.ResetForFaceCount( faceSearch.Size() );

    // The midpoint exists geometrically, but neither child edge is registered.
    // SplitLine must return before looking up or mutating the zone grid.
    EXPECT_NO_THROW( faceSearch.SplitLine( link, 0 ) );

    EXPECT_TRUE( faceSearch.interfaceSplit.childFaceIdsOfParent[ 0 ].empty() );
    EXPECT_TRUE( faceSearch.interfaceSplit.childFaceNodeIdsByParent[ 0 ].empty() );
    EXPECT_TRUE( faceSearch.interfaceSplit.childFaceNodeFlagsByParent[ 0 ].empty() );
}

