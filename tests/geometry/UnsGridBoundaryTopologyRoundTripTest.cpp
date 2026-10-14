#include <gtest/gtest.h>

#include "BcRecord.h"
#include "Boundary.h"
#include "DataBook.h"
#include "DataBaseIO.h"
#include "FaceTopo.h"
#include "InterFace.h"
#include "UnsGrid.h"

namespace
{
void ConfigureBoundaryTopology( ONEFLOW::UnsGrid & grid )
{
    auto & faceTopo = grid.GetFaceTopo();
    faceTopo.SetNBFaces( 4 );

    auto & bcRecord = faceTopo.GetBcRecord();
    bcRecord.bcType = {
        ONEFLOW::BC::INTERFACE,
        ONEFLOW::BC::SOLID_SURFACE,
        ONEFLOW::BC::INTERFACE,
        ONEFLOW::BC::INTERFACE
    };
    bcRecord.bcNameId = { 10, 20, 30, 40 };

    grid.nIFaces = 3;
    grid.interFace->Set( 3, &grid );
    grid.interFace->zoneId = { 1, 1, 2 };
    grid.interFace->localInterfaceId = { 4, 7, 0 };
    grid.interFace->i2b = { 0, 2, 3 };
}
}

TEST( UnsGridBoundaryTopologyRoundTripTest, DataBookPreservesBoundaryAndInterfaceMappings )
{
    ONEFLOW::UnsGrid source;
    source.Init();
    ConfigureBoundaryTopology( source );

    ONEFLOW::DataBook dataBook( 16 );
    ASSERT_NO_THROW( source.WriteBoundaryTopology( &dataBook ) );

    dataBook.MoveToBegin();

    ONEFLOW::UnsGrid restored;
    restored.Init();
    ASSERT_NO_THROW( restored.ReadBoundaryTopology( &dataBook ) );

    const auto & sourceBc = source.GetFaceTopo().GetBcRecord();
    const auto & restoredBc = restored.GetFaceTopo().GetBcRecord();
    ASSERT_EQ( restored.GetFaceTopo().GetNBFaces(), 4u );
    EXPECT_EQ( restoredBc.bcType, sourceBc.bcType );
    EXPECT_EQ( restoredBc.bcNameId, sourceBc.bcNameId );
    EXPECT_EQ( restored.nIFaces, 3 );
    ASSERT_NE( restored.interFace, nullptr );
    EXPECT_EQ( restored.interFace->zoneId, source.interFace->zoneId );
    EXPECT_EQ( restored.interFace->localInterfaceId, source.interFace->localInterfaceId );
    EXPECT_EQ( restored.interFace->i2b, ( ONEFLOW::IntField{ 0, 2, 3 } ) );

    for ( int interfaceId = 0; interfaceId < restored.nIFaces; ++ interfaceId )
    {
        const int boundaryId = restored.interFace->i2b[ interfaceId ];
        ASSERT_GE( boundaryId, 0 );
        ASSERT_LT( boundaryId, static_cast< int >( restoredBc.bcType.size() ) );
        EXPECT_TRUE( ONEFLOW::BC::IsInterfaceBc( restoredBc.bcType[ boundaryId ] ) );
    }
}

TEST( UnsGridBoundaryTopologyRoundTripTest, SplitInterfaceRecordsRegenerateBoundaryMap )
{
    ONEFLOW::BcRecord bcRecord;
    bcRecord.bcType = {
        ONEFLOW::BC::INTERFACE,
        ONEFLOW::BC::INTERFACE,
        ONEFLOW::BC::SOLID_SURFACE,
        ONEFLOW::BC::INTERFACE,
        ONEFLOW::BC::INTERFACE
    };
    bcRecord.bcNameId = { 10, 10, 20, 30, 30 };

    ONEFLOW::InterFace interFace;
    interFace.Set( 4 );
    bcRecord.GenerateI2B( interFace );

    EXPECT_EQ( interFace.i2b, ( ONEFLOW::IntField{ 0, 1, 3, 4 } ) );
    for ( int interfaceId = 0; interfaceId < interFace.nIFaces; ++ interfaceId )
    {
        const int boundaryId = interFace.i2b[ interfaceId ];
        ASSERT_GE( boundaryId, 0 );
        ASSERT_LT( boundaryId, static_cast< int >( bcRecord.bcType.size() ) );
        EXPECT_TRUE( ONEFLOW::BC::IsInterfaceBc( bcRecord.bcType[ boundaryId ] ) );
    }

    EXPECT_EQ( bcRecord.bcNameId[ interFace.i2b[ 0 ] ], 10 );
    EXPECT_EQ( bcRecord.bcNameId[ interFace.i2b[ 1 ] ], 10 );
    EXPECT_EQ( bcRecord.bcNameId[ interFace.i2b[ 2 ] ], 30 );
    EXPECT_EQ( bcRecord.bcNameId[ interFace.i2b[ 3 ] ], 30 );
}
