#include <gtest/gtest.h>

#include "InterFace.h"
#include "Parallel.h"
#include "SlipFace.h"

#include <memory>
#include <stdexcept>

namespace
{

TEST( InterfaceNeighborMappingTest, InterfaceUsesZoneMapInsteadOfLoopIndex )
{
    ONEFLOW::InterFace interFace;
    interFace.interFacePairs.resize( 2 );
    interFace.interFacePairs[ 0 ] =
        std::make_unique< ONEFLOW::InterfacePair >();
    interFace.interFacePairs[ 1 ] =
        std::make_unique< ONEFLOW::InterfacePair >();
    interFace.interFacePairs[ 0 ]->nzid = 7;
    interFace.interFacePairs[ 1 ]->nzid = 3;
    interFace.z2n = { { 7, 0 }, { 3, 1 } };

    EXPECT_EQ( interFace.GetNeighborIndex( 3 ), 1 );
    EXPECT_EQ( interFace.GetNeighborIndex( 7 ), 0 );

    ONEFLOW::IntField sendIds = { 11, 13 };
    interFace.SetSendId( 3, sendIds );
    EXPECT_EQ( interFace.interFacePairs[ 1 ]->idsend, sendIds );
    EXPECT_TRUE( interFace.interFacePairs[ 0 ]->idsend.empty() );
}

TEST( InterfaceNeighborMappingTest, InterfaceMissingLookupDoesNotMutateMap )
{
    ONEFLOW::InterFace interFace;
    const auto originalSize = interFace.z2n.size();

    EXPECT_THROW( interFace.GetNeighborIndex( 4 ), std::out_of_range );
    EXPECT_EQ( interFace.z2n.size(), originalSize );
}

TEST( InterfaceNeighborMappingTest, InterfaceRejectsCorruptedMap )
{
    ONEFLOW::InterFace interFace;
    interFace.interFacePairs.resize( 1 );
    interFace.interFacePairs[ 0 ] =
        std::make_unique< ONEFLOW::InterfacePair >();
    interFace.interFacePairs[ 0 ]->nzid = 2;

    interFace.z2n[ 2 ] = 5;
    EXPECT_THROW( interFace.GetNeighborIndex( 2 ), std::logic_error );

    interFace.z2n[ 2 ] = 0;
    interFace.interFacePairs[ 0 ]->nzid = 3;
    EXPECT_THROW( interFace.GetNeighborIndex( 2 ), std::logic_error );
}

TEST( InterfaceNeighborMappingTest, SlipFaceUsesZoneMapAndAllowsEmptyIds )
{
    ONEFLOW::SlipFace slipFace;
    slipFace.slipfacePairs.resize( 2 );
    slipFace.slipfacePairs[ 0 ] =
        std::make_unique< ONEFLOW::SlipfacePair >();
    slipFace.slipfacePairs[ 1 ] =
        std::make_unique< ONEFLOW::SlipfacePair >();
    slipFace.slipfacePairs[ 0 ]->nzid = 9;
    slipFace.slipfacePairs[ 1 ]->nzid = 4;
    slipFace.z2n = { { 9, 0 }, { 4, 1 } };

    EXPECT_EQ( slipFace.GetNeighborIndex( 4 ), 1 );

    ONEFLOW::IntField emptyIds;
    EXPECT_NO_THROW( slipFace.SetSendId( 4, emptyIds ) );
    EXPECT_TRUE( slipFace.slipfacePairs[ 1 ]->idsend.empty() );
}

TEST( InterfaceNeighborMappingTest, ZeroLengthSwapIsANoOp )
{
    EXPECT_NO_THROW(
        ONEFLOW::HXSwapData(
            static_cast< int * >( nullptr ),
            0,
            0,
            1 ) );
}

}
