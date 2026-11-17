#include <gtest/gtest.h>

#include "DistributedDomainPlan.h"

#include <stdexcept>
#include <vector>

namespace
{

using ONEFLOW::BoundaryExchange;
using ONEFLOW::DistributedDomainPlan;

TEST( DistributedDomainPlanTest, BuildsDeterministicPairSchedule )
{
    const std::vector< int > owners = { 0, 1, 2, 3 };
    const std::vector< std::vector< int > > neighbors = {
        { 2, 1, 1 },
        { 3, 0 },
        { 0 },
        { 1 }
    };

    const DistributedDomainPlan plan =
        DistributedDomainPlan::Build( owners, neighbors, 4 );

    EXPECT_EQ(
        plan.BoundaryExchanges(),
        ( std::vector< BoundaryExchange >{
            { 0, 1 },
            { 1, 0 },
            { 0, 2 },
            { 2, 0 },
            { 1, 3 },
            { 3, 1 }
        } ) );
}

TEST( DistributedDomainPlanTest, OneSidedPairStillCreatesBothDirections )
{
    const DistributedDomainPlan plan = DistributedDomainPlan::Build(
        { 0, 1 },
        { { 1 }, {} },
        2 );

    EXPECT_EQ(
        plan.BoundaryExchanges(),
        ( std::vector< BoundaryExchange >{
            { 0, 1 },
            { 1, 0 }
        } ) );
}

TEST( DistributedDomainPlanTest, NeighborOrderDoesNotChangeSchedule )
{
    const DistributedDomainPlan first = DistributedDomainPlan::Build(
        { 0, 1, 2 },
        { { 2, 1, 1 }, { 2, 0 }, { 1, 0 } },
        3 );
    const DistributedDomainPlan reordered = DistributedDomainPlan::Build(
        { 0, 1, 2 },
        { { 1, 2 }, { 0, 2 }, { 0, 1 } },
        3 );

    EXPECT_EQ(
        first.BoundaryExchanges(),
        reordered.BoundaryExchanges() );
}

TEST( DistributedDomainPlanTest, SelfNeighborIsSkipped )
{
    const DistributedDomainPlan plan = DistributedDomainPlan::Build(
        { 0, 1 },
        { { 1, 0, 0 }, { 0 } },
        2 );

    EXPECT_EQ(
        plan.BoundaryExchanges(),
        ( std::vector< BoundaryExchange >{
            { 0, 1 },
            { 1, 0 }
        } ) );
}

TEST( DistributedDomainPlanTest, StoresUniqueZoneOwnership )
{
    const DistributedDomainPlan plan = DistributedDomainPlan::Build(
        { 2, 0, 1 },
        { {}, {}, {} },
        3 );

    EXPECT_EQ( plan.ZoneCount(), 3 );
    EXPECT_EQ( plan.OwnerRanks(), ( std::vector< int >{ 2, 0, 1 } ) );
    EXPECT_EQ( plan.OwnerRank( 0 ), 2 );
    EXPECT_EQ( plan.OwnerRank( 1 ), 0 );
    EXPECT_EQ( plan.OwnerRank( 2 ), 1 );
    EXPECT_THROW( plan.OwnerRank( -1 ), std::out_of_range );
    EXPECT_THROW( plan.OwnerRank( 3 ), std::out_of_range );
}

TEST( DistributedDomainPlanTest, RejectsInvalidOwnershipAndTopology )
{
    EXPECT_THROW(
        DistributedDomainPlan::Build( { 0 }, { {} }, 0 ),
        std::invalid_argument );
    EXPECT_THROW(
        DistributedDomainPlan::Build( { 0 }, {}, 1 ),
        std::invalid_argument );
    EXPECT_THROW(
        DistributedDomainPlan::Build( { -1 }, { {} }, 1 ),
        std::out_of_range );
    EXPECT_THROW(
        DistributedDomainPlan::Build( { 1 }, { {} }, 1 ),
        std::out_of_range );
    EXPECT_THROW(
        DistributedDomainPlan::Build( { 0 }, { { -1 } }, 1 ),
        std::out_of_range );
    EXPECT_THROW(
        DistributedDomainPlan::Build( { 0 }, { { 1 } }, 1 ),
        std::out_of_range );
}

}
