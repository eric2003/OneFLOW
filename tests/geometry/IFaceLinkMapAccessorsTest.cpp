// tests/geometry/IFaceLinkMapAccessorsTest.cpp
// Locks dual-map accessors on IFaceLink without full mesh scan.
// Fields are private; seed through Ensure*/Register*.

#include <gtest/gtest.h>
#include <type_traits>
#include <utility>

#include "IFaceLink.h"
#include "GridHandles.h"
#include "StrGrid.h"
#include "HXDefine.h"

using namespace ONEFLOW;

static_assert(
    std::is_same_v<
        decltype( GridAt( std::declval< const ConstGridViews & >(), std::size_t{} ) ),
        const Grid * >,
    "ConstGridViews must expose const Grid objects" );

TEST( IFaceLinkMapAccessorsTest, ConstAccessorsReflectEmptyTables )
{
    Grids grids;
    IFaceLink link( grids );

    EXPECT_TRUE( link.GetLocalToGlobalMap().empty() );
    EXPECT_TRUE( link.GetGlobalToZoneIds().empty() );
    EXPECT_TRUE( link.GetGlobalToLocalIds().empty() );
}

TEST( IFaceLinkMapAccessorsTest, AccessorsSeeRegisteredTables )
{
    Grids grids;
    IFaceLink link( grids );
    link.EnsureInterfaceZoneCount( 1 );
    link.EnsureLocalInterfaceFaceCount( 0, 1 );
    link.RegisterLocalInterfaceFace( /*gIid=*/0, /*isNew=*/true, /*zid=*/0, /*localFaceId=*/0 );

    ASSERT_EQ( link.GetLocalToGlobalMap().size(), static_cast< std::size_t >( 1 ) );
    EXPECT_EQ( link.GetLocalToGlobalMap()[ 0 ], ( IntField{ 0 } ) );
    EXPECT_EQ( link.GetGlobalToZoneIds()[ 0 ], ( IntField{ 0 } ) );
    EXPECT_EQ( link.GetGlobalToLocalIds()[ 0 ], ( IntField{ 0 } ) );
}


TEST( IFaceLinkMapAccessorsTest, NewInverseMapAccessorsAliasStagingTables )
{
    Grids grids;
    IFaceLink link( grids );
    link.MutableGlobalToZoneIdsNew().push_back( IntField{ 1 } );
    link.MutableGlobalToLocalIdsNew().push_back( IntField{ 2 } );

    EXPECT_EQ( link.GetGlobalToZoneIdsNew()[ 0 ], ( IntField{ 1 } ) );
    EXPECT_EQ( link.GetGlobalToLocalIdsNew()[ 0 ], ( IntField{ 2 } ) );
    EXPECT_EQ( &link.MutableGlobalToZoneIdsNew(), &link.GetGlobalToZoneIdsNew() );
}

TEST( IFaceLinkMapAccessorsTest, ConstGetGridReturnsReadOnlyZoneView )
{
    Grids grids;
    grids.push_back( std::make_unique< StrGrid >() );
    IFaceLink link( grids );
    const IFaceLink & constLink = link;

    EXPECT_EQ( &constLink.GetGrid( 0 ), grids[ 0 ].get() );
}

TEST( GridHandlesTest, ConstViewPreservesReadOnlyGridType )
{
    Grids grids;
    grids.push_back( std::make_unique< StrGrid >() );

    const ConstGridViews views = AsConstGridViews( grids );

    ASSERT_EQ( GridsSize( views ), 1 );
    EXPECT_EQ( GridAt( views, 0 ), grids[ 0 ].get() );
}
