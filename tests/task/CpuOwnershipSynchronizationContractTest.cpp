#include <gtest/gtest.h>

#include "BoundaryExchangeContract.h"
#include "DistributedDomainPlan.h"

#include <stdexcept>
#include <vector>

namespace
{

using ONEFLOW::BoundaryEndpointMap;
using ONEFLOW::BoundaryEndpointRole;
using ONEFLOW::BoundaryExchangeContract;
using ONEFLOW::BoundaryMemorySpace;
using ONEFLOW::BoundaryBufferView;
using ONEFLOW::BoundaryTransferMode;
using ONEFLOW::DistributedDomainPlan;
using ONEFLOW::PackBoundaryHost;
using ONEFLOW::Real;
using ONEFLOW::UnpackBoundaryHost;
using ONEFLOW::ValidateBoundaryExchangeBuffers;

BoundaryExchangeContract MakeCpuInterfaceContract()
{
    BoundaryExchangeContract contract;
    contract.sourceZone = 0;
    contract.targetZone = 1;
    contract.sourceRank = 0;
    contract.targetRank = 1;
    contract.nEquations = 2;
    contract.source = BoundaryEndpointMap{
        4, { 0, 2 }, 3, BoundaryEndpointRole::OwnedCell };
    contract.target = BoundaryEndpointMap{
        5, { 3, 4 }, 3, BoundaryEndpointRole::GhostCell };
    return contract;
}

TEST( CpuOwnershipSynchronizationContractTest,
      ScheduleSkipsSelfNeighborsAndKeepsOwnership )
{
    const DistributedDomainPlan plan = DistributedDomainPlan::Build(
        { 0, 0, 1 },
        { { 0, 1, 2, 2 }, { 0, 1, 2 }, { 0, 2 } },
        2 );

    EXPECT_EQ( plan.OwnerRanks(), ( std::vector< int >{ 0, 0, 1 } ) );
    EXPECT_EQ(
        plan.BoundaryExchanges(),
        ( std::vector< ONEFLOW::BoundaryExchange >{
            { 0, 1 }, { 1, 0 },
            { 0, 2 }, { 2, 0 },
            { 1, 2 }, { 2, 1 }
        } ) );
    for ( const auto & exchange : plan.BoundaryExchanges() )
        EXPECT_NE( exchange.sourceZone, exchange.targetZone );
}

TEST( CpuOwnershipSynchronizationContractTest,
      EmptyTopologyProducesNoMessages )
{
    const DistributedDomainPlan plan = DistributedDomainPlan::Build(
        { 0, 1 }, { { 0 }, { 1 } }, 2 );

    EXPECT_TRUE( plan.BoundaryExchanges().empty() );
    EXPECT_EQ( plan.ZoneCount(), 2 );
    EXPECT_EQ( plan.OwnerRank( 0 ), 0 );
    EXPECT_EQ( plan.OwnerRank( 1 ), 1 );
}

TEST( CpuOwnershipSynchronizationContractTest,
      HostExchangeCopiesOwnedValuesIntoGhostCells )
{
    const BoundaryExchangeContract contract = MakeCpuInterfaceContract();
    ASSERT_NO_THROW( contract.Validate() );

    const std::vector< Real > source0 = { 10, 11, 12, 13 };
    const std::vector< Real > source1 = { 20, 21, 22, 23 };
    std::vector< Real > target0( 5, -1 );
    std::vector< Real > target1( 5, -1 );
    std::vector< Real > packed;

    PackBoundaryHost(
        contract.source,
        contract.nEquations,
        { source0.data(), source1.data() },
        packed );
    ASSERT_EQ( packed, ( std::vector< Real >{ 10, 12, 20, 22 } ) );

    BoundaryBufferView send{
        packed.data(), packed.size(), BoundaryMemorySpace::Host };
    BoundaryBufferView receive{
        packed.data(), packed.size(), BoundaryMemorySpace::Host };
    ASSERT_NO_THROW(
        ValidateBoundaryExchangeBuffers(
            contract, send, receive, BoundaryTransferMode::HostStaged ) );

    UnpackBoundaryHost(
        contract.target,
        contract.nEquations,
        packed,
        { target0.data(), target1.data() } );

    EXPECT_EQ( target0, ( std::vector< Real >{ -1, -1, -1, 10, 12 } ) );
    EXPECT_EQ( target1, ( std::vector< Real >{ -1, -1, -1, 20, 22 } ) );
}

TEST( CpuOwnershipSynchronizationContractTest,
      CrossRankReverseDirectionUsesItsOwnEndpointRoles )
{
    BoundaryExchangeContract forward = MakeCpuInterfaceContract();
    ASSERT_NO_THROW( forward.Validate() );
    EXPECT_EQ( forward.sourceRank, 0 );
    EXPECT_EQ( forward.targetRank, 1 );

    BoundaryExchangeContract reverse;
    reverse.sourceZone = forward.targetZone;
    reverse.targetZone = forward.sourceZone;
    reverse.sourceRank = forward.targetRank;
    reverse.targetRank = forward.sourceRank;
    reverse.nEquations = forward.nEquations;
    reverse.source = BoundaryEndpointMap{
        5, { 0, 2 }, 3, BoundaryEndpointRole::OwnedCell };
    reverse.target = BoundaryEndpointMap{
        5, { 3, 4 }, 3, BoundaryEndpointRole::GhostCell };

    ASSERT_NO_THROW( reverse.Validate() );
    EXPECT_EQ( reverse.InterfaceCount(), forward.InterfaceCount() );
    EXPECT_EQ( reverse.ValueCount(), forward.ValueCount() );
    EXPECT_NE( reverse.sourceZone, reverse.targetZone );
    EXPECT_NE( reverse.sourceRank, reverse.targetRank );
}

TEST( CpuOwnershipSynchronizationContractTest,
      EmptyInterfaceLeavesGhostStateUnchanged )
{
    BoundaryExchangeContract contract = MakeCpuInterfaceContract();
    contract.source.cellIds.clear();
    contract.target.cellIds.clear();
    ASSERT_NO_THROW( contract.Validate() );

    std::vector< Real > target0 = { 7, 8, 9, 101, 102 };
    std::vector< Real > target1 = { 17, 18, 19, 201, 202 };
    const std::vector< Real > before0 = target0;
    const std::vector< Real > before1 = target1;
    std::vector< Real > packed = { 99 };

    BoundaryBufferView empty;
    ASSERT_NO_THROW(
        ValidateBoundaryExchangeBuffers(
            contract, empty, empty, BoundaryTransferMode::HostStaged ) );
    PackBoundaryHost(
        contract.source, contract.nEquations,
        { target0.data(), target1.data() }, packed );
    EXPECT_TRUE( packed.empty() );
    UnpackBoundaryHost(
        contract.target, contract.nEquations, packed,
        { target0.data(), target1.data() } );

    EXPECT_EQ( target0, before0 );
    EXPECT_EQ( target1, before1 );
}

TEST( CpuOwnershipSynchronizationContractTest,
      InterfaceRequiresOwnedSourceAndGhostTarget )
{
    BoundaryExchangeContract contract = MakeCpuInterfaceContract();

    contract.source.cellIds.clear();
    contract.source.role = BoundaryEndpointRole::GhostCell;
    EXPECT_THROW( contract.Validate(), std::invalid_argument );

    contract = MakeCpuInterfaceContract();
    contract.target.cellIds.clear();
    contract.target.role = BoundaryEndpointRole::OwnedCell;
    EXPECT_THROW( contract.Validate(), std::invalid_argument );

    contract = MakeCpuInterfaceContract();
    contract.source.cellIds[ 0 ] = contract.source.ownedCellCount;
    EXPECT_THROW( contract.Validate(), std::out_of_range );

    contract = MakeCpuInterfaceContract();
    contract.target.cellIds[ 0 ] = contract.target.ownedCellCount - 1;
    EXPECT_THROW( contract.Validate(), std::out_of_range );
}

} // namespace
