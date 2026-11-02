#include <gtest/gtest.h>

#include "BoundaryExchangeContract.h"

#include <stdexcept>
#include <vector>

namespace
{

using ONEFLOW::BoundaryFieldExchangePlan;
using ONEFLOW::BoundaryEndpointMap;
using ONEFLOW::BoundaryEndpointRole;
using ONEFLOW::BoundaryExchangeContract;
using ONEFLOW::BoundaryBufferView;
using ONEFLOW::BoundaryMemorySpace;
using ONEFLOW::BoundaryTransferMode;
using ONEFLOW::PackBoundaryHost;
using ONEFLOW::Real;
using ONEFLOW::UnpackBoundaryHost;
using ONEFLOW::ValidateBoundaryBuffer;
using ONEFLOW::ValidateBoundaryExchangeBuffers;

BoundaryFieldExchangePlan MakePlan()
{
    BoundaryFieldExchangePlan plan;
    plan.sourceZone = 2;
    plan.targetZone = 7;
    plan.sourceRank = 0;
    plan.targetRank = 1;
    plan.ghostLevel = 1;
    plan.nEquations = 2;
    plan.sourceCellCount = 5;
    plan.targetCellCount = 6;
    plan.sendCellIds = { 4, 1 };
    plan.receiveCellIds = { 3, 5 };
    return plan;
}

TEST( BoundaryExchangeContractTest, PacksAndUnpacksEquationMajorValues )
{
    const BoundaryFieldExchangePlan plan = MakePlan();
    const std::vector< Real > source0 = { 10, 11, 12, 13, 14 };
    const std::vector< Real > source1 = { 20, 21, 22, 23, 24 };
    std::vector< Real > packed;

    PackBoundaryHost(
        plan,
        { source0.data(), source1.data() },
        packed );

    EXPECT_EQ(
        packed,
        ( std::vector< Real >{ 14, 11, 24, 21 } ) );

    std::vector< Real > target0( 6, -1 );
    std::vector< Real > target1( 6, -1 );
    UnpackBoundaryHost(
        plan,
        packed,
        { target0.data(), target1.data() } );

    EXPECT_EQ( target0[ 3 ], 14 );
    EXPECT_EQ( target0[ 5 ], 11 );
    EXPECT_EQ( target1[ 3 ], 24 );
    EXPECT_EQ( target1[ 5 ], 21 );
    EXPECT_EQ( target0[ 0 ], -1 );
    EXPECT_EQ( target1[ 4 ], -1 );
}

TEST( BoundaryExchangeContractTest, SupportsZeroLengthExchange )
{
    BoundaryFieldExchangePlan plan = MakePlan();
    plan.sendCellIds.clear();
    plan.receiveCellIds.clear();
    std::vector< Real > packed = { 1 };
    const std::vector< Real > source0( 5, 0 );
    const std::vector< Real > source1( 5, 0 );

    PackBoundaryHost(
        plan,
        { source0.data(), source1.data() },
        packed );

    EXPECT_TRUE( packed.empty() );
    std::vector< Real > target0( 6, -1 );
    std::vector< Real > target1( 6, -1 );
    EXPECT_NO_THROW(
        UnpackBoundaryHost(
            plan,
            packed,
            { target0.data(), target1.data() } ) );
}

TEST( BoundaryExchangeContractTest, RejectsInconsistentMappings )
{
    BoundaryFieldExchangePlan plan = MakePlan();
    plan.receiveCellIds.pop_back();
    EXPECT_THROW( plan.Validate(), std::invalid_argument );

    plan = MakePlan();
    plan.sendCellIds[ 0 ] = plan.sourceCellCount;
    EXPECT_THROW( plan.Validate(), std::out_of_range );

    plan = MakePlan();
    plan.receiveCellIds[ 0 ] = -1;
    EXPECT_THROW( plan.Validate(), std::out_of_range );
}

TEST( BoundaryExchangeContractTest, ValidatesOwnedAndGhostEndpointRoles )
{
    BoundaryEndpointMap send;
    send.cellCount = 6;
    send.ownedCellCount = 4;
    send.role = BoundaryEndpointRole::OwnedCell;
    send.cellIds = { 0, 3 };
    EXPECT_NO_THROW( send.Validate() );

    BoundaryEndpointMap receive;
    receive.cellCount = 6;
    receive.ownedCellCount = 4;
    receive.role = BoundaryEndpointRole::GhostCell;
    receive.cellIds = { 4, 5 };
    EXPECT_NO_THROW( receive.Validate() );

    send.cellIds[ 1 ] = 4;
    EXPECT_THROW( send.Validate(), std::out_of_range );

    receive.cellIds[ 0 ] = 2;
    EXPECT_THROW( receive.Validate(), std::out_of_range );
}

TEST( BoundaryExchangeContractTest, RejectsIncompleteEndpointOwnership )
{
    BoundaryEndpointMap endpoint;
    endpoint.cellCount = 4;
    endpoint.ownedCellCount = 5;
    endpoint.role = BoundaryEndpointRole::OwnedCell;
    EXPECT_THROW( endpoint.Validate(), std::invalid_argument );

    endpoint.ownedCellCount = -1;
    EXPECT_THROW( endpoint.Validate(), std::invalid_argument );
}

TEST( BoundaryExchangeContractTest, ValidatesOwnedToGhostExchangeBuffers )
{
    BoundaryExchangeContract contract;
    contract.sourceZone = 2;
    contract.targetZone = 7;
    contract.sourceRank = 0;
    contract.targetRank = 1;
    contract.nEquations = 2;
    contract.source.cellCount = 6;
    contract.source.ownedCellCount = 4;
    contract.source.role = BoundaryEndpointRole::OwnedCell;
    contract.source.cellIds = { 0, 3 };
    contract.target.cellCount = 7;
    contract.target.ownedCellCount = 5;
    contract.target.role = BoundaryEndpointRole::GhostCell;
    contract.target.cellIds = { 5, 6 };

    Real sendValues[ 4 ] = {};
    Real receiveValues[ 4 ] = {};
    BoundaryBufferView send{
        sendValues, 4, BoundaryMemorySpace::Host };
    BoundaryBufferView receive{
        receiveValues, 4, BoundaryMemorySpace::Host };

    EXPECT_NO_THROW(
        ValidateBoundaryExchangeBuffers(
            contract, send, receive,
            BoundaryTransferMode::HostStaged ) );
    EXPECT_EQ( contract.InterfaceCount(), 2u );
    EXPECT_EQ( contract.ValueCount(), 4u );
}

TEST( BoundaryExchangeContractTest, RejectsWrongMemorySpaceOrEndpointRole )
{
    BoundaryExchangeContract contract;
    contract.sourceZone = 0;
    contract.targetZone = 1;
    contract.sourceRank = 0;
    contract.targetRank = 0;
    contract.nEquations = 1;
    contract.source.cellCount = 3;
    contract.source.ownedCellCount = 2;
    contract.source.role = BoundaryEndpointRole::OwnedCell;
    contract.source.cellIds = { 1 };
    contract.target.cellCount = 3;
    contract.target.ownedCellCount = 2;
    contract.target.role = BoundaryEndpointRole::GhostCell;
    contract.target.cellIds = { 2 };

    Real values[ 1 ] = {};
    BoundaryBufferView deviceBuffer{
        values, 1, BoundaryMemorySpace::Device };
    EXPECT_THROW(
        ValidateBoundaryBuffer(
            contract.source, contract.nEquations, deviceBuffer,
            BoundaryTransferMode::HostStaged ),
        std::invalid_argument );

    contract.target.role = BoundaryEndpointRole::OwnedCell;
    contract.target.cellIds[ 0 ] = 1;
    EXPECT_THROW( contract.Validate(), std::invalid_argument );
}

TEST( BoundaryExchangeContractTest, AllowsEmptyInterfaceWithEmptyBuffers )
{
    BoundaryExchangeContract contract;
    contract.sourceZone = 0;
    contract.targetZone = 1;
    contract.sourceRank = 0;
    contract.targetRank = 0;
    contract.nEquations = 5;
    contract.source.cellCount = 2;
    contract.source.ownedCellCount = 2;
    contract.source.role = BoundaryEndpointRole::OwnedCell;
    contract.target.cellCount = 3;
    contract.target.ownedCellCount = 2;
    contract.target.role = BoundaryEndpointRole::GhostCell;

    BoundaryBufferView emptySend;
    BoundaryBufferView emptyReceive;
    EXPECT_NO_THROW(
        ValidateBoundaryExchangeBuffers(
            contract, emptySend, emptyReceive,
            BoundaryTransferMode::HostStaged ) );
}

}
