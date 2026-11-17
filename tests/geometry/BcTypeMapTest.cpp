// tests/geometry/BcTypeMapTest.cpp
// Locks CGNS -> OneFLOW BC type mapping without GridElem or a full mesh.

#include <gtest/gtest.h>
#include <iostream>
#include <sstream>
#include <string>

#include "Boundary.h"
#include "HXCgns.h"
#include "HXDefine.h"

using namespace ONEFLOW;

TEST( BcTypeMapTest, MapFromCgnsInPlaceConvertsKnownTypes )
{
    BcTypeMap map;
    map.Init();

    IntField types{
        static_cast< int >( CGNS_ENUMV( BCWall ) ),
        static_cast< int >( CGNS_ENUMV( BCInflow ) ),
        static_cast< int >( CGNS_ENUMV( BCFarfield ) ),
        static_cast< int >( CGNS_ENUMV( BCTypeNull ) )
    };

    map.MapFromCgnsInPlace( types );

    EXPECT_EQ( types[ 0 ], BC::SOLID_SURFACE );
    EXPECT_EQ( types[ 1 ], BC::INFLOW );
    EXPECT_EQ( types[ 2 ], BC::FARFIELD );
    EXPECT_EQ( types[ 3 ], BC::INTERFACE );
}

TEST( BcTypeMapTest, UnknownCgnsTypeMapsToGeneric2 )
{
    BcTypeMap map;
    map.Init();

    // Arbitrary id not present in the table.
    EXPECT_EQ( map.Cgns2OneFlow( 999999 ), BC::GENERIC_2 );
}

TEST( BcTypeMapTest, MapFromCgnsInPlaceEmptyIsNoOp )
{
    BcTypeMap map;
    map.Init();
    IntField empty;
    EXPECT_NO_THROW( map.MapFromCgnsInPlace( empty ) );
    EXPECT_TRUE( empty.empty() );
}

TEST( BcTypeMapTest, RoundTripCommonOneFlowTypes )
{
    BcTypeMap map;
    map.Init();

    EXPECT_EQ( map.Cgns2OneFlow( map.OneFlow2Cgns( BC::INFLOW ) ), BC::INFLOW );
    EXPECT_EQ( map.Cgns2OneFlow( map.OneFlow2Cgns( BC::OUTFLOW ) ), BC::OUTFLOW );
    EXPECT_EQ( map.Cgns2OneFlow( map.OneFlow2Cgns( BC::SOLID_SURFACE ) ), BC::SOLID_SURFACE );
}

TEST( BcTypeMapTest, GmshPhysicalNamesMapToStandardCgnsTypes )
{
    EXPECT_EQ( GmshPhysicalNameToCgnsBcType( "Inlet" ),
               static_cast< int >( CGNS_ENUMV( BCInflow ) ) );
    EXPECT_EQ( GmshPhysicalNameToCgnsBcType( "Velocity-Inlet" ),
               static_cast< int >( CGNS_ENUMV( BCInflow ) ) );
    EXPECT_EQ( GmshPhysicalNameToCgnsBcType( "Outlet" ),
               static_cast< int >( CGNS_ENUMV( BCOutflow ) ) );
    EXPECT_EQ( GmshPhysicalNameToCgnsBcType( "LowerWall" ),
               static_cast< int >( CGNS_ENUMV( BCWall ) ) );
    EXPECT_EQ( GmshPhysicalNameToCgnsBcType( "NoSlipWall" ),
               static_cast< int >( CGNS_ENUMV( BCWallViscous ) ) );
    EXPECT_EQ( GmshPhysicalNameToCgnsBcType( "SlipWall" ),
               static_cast< int >( CGNS_ENUMV( BCWallInviscid ) ) );
    EXPECT_EQ( GmshPhysicalNameToCgnsBcType( "Symmetry" ),
               static_cast< int >( CGNS_ENUMV( BCSymmetryPlane ) ) );
    EXPECT_EQ( GmshPhysicalNameToCgnsBcType( "Farfield" ),
               static_cast< int >( CGNS_ENUMV( BCFarfield ) ) );
    EXPECT_EQ( GmshPhysicalNameToCgnsBcType( "UnclassifiedPatch" ),
               static_cast< int >( CGNS_ENUMV( BCTypeUserDefined ) ) );

    EXPECT_EQ( GmshPhysicalNameToCgnsBcType( "xmin_inflow" ),
               static_cast< int >( CGNS_ENUMV( BCInflow ) ) );
    EXPECT_EQ( GmshPhysicalNameToCgnsBcType( "xmax_outflow" ),
               static_cast< int >( CGNS_ENUMV( BCOutflow ) ) );
    EXPECT_EQ( GmshPhysicalNameToCgnsBcType( "ymin_wall" ),
               static_cast< int >( CGNS_ENUMV( BCWall ) ) );
    EXPECT_EQ( GmshPhysicalNameToCgnsBcType( "ymax_farfield" ),
               static_cast< int >( CGNS_ENUMV( BCFarfield ) ) );
    EXPECT_EQ( GmshPhysicalNameToCgnsBcType( "zmin_symmetry" ),
               static_cast< int >( CGNS_ENUMV( BCSymmetryPlane ) ) );
    EXPECT_EQ( GmshPhysicalNameToCgnsBcType( "zmax_symmetry" ),
               static_cast< int >( CGNS_ENUMV( BCSymmetryPlane ) ) );
}

TEST( BcTypeMapTest, PrintSummaryListsUniqueTypesAndOnlyUsedMappings )
{
    BcTypeMap map;
    map.Init();

    IntField cgnsTypes{
        static_cast< int >( CGNS_ENUMV( BCWall ) ),
        static_cast< int >( CGNS_ENUMV( BCWall ) ),
        static_cast< int >( CGNS_ENUMV( BCInflow ) ),
        static_cast< int >( CGNS_ENUMV( BCFarfield ) )
    };

    std::ostringstream output;
    std::streambuf * oldBuffer = std::cout.rdbuf( output.rdbuf() );
    map.PrintSummary( cgnsTypes );
    std::cout.rdbuf( oldBuffer );

    const std::string summary = output.str();
    EXPECT_NE( summary.find( "OneFLOW types (3)" ), std::string::npos );
    EXPECT_NE( summary.find( "SOLID_SURFACE" ), std::string::npos );
    EXPECT_NE( summary.find( "INFLOW" ), std::string::npos );
    EXPECT_NE( summary.find( "FARFIELD" ), std::string::npos );
    EXPECT_EQ( summary.find( "BCOutflow" ), std::string::npos );
}
