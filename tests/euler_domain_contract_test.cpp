#include "EulerDomain.h"
#include "AccelViews.h"

#include <gtest/gtest.h>

#include <stdexcept>

namespace
{

using namespace ONEFLOW;

EulerDomainProblem Problem()
{
    EulerDomainProblem result;
    result.nCells = 32;
    result.dt = 0.0005;
    result.dx = 1.0 / result.nCells;
    return result;
}

EulerDomainProblem FiveEquationProblem()
{
    EulerDomainProblem result = Problem();
    result.nEquations = 5;
    return result;
}

TEST( EulerDomainContract, ValidatesBackendNeutralProblemAndViews )
{
    const EulerDomainProblem problem = Problem();
    Real values[ 3 * 32 ] = {};
    EulerDomainConstFieldView input{ 32, 3, values };
    EulerDomainFieldView output{ 32, 3, values };

    EXPECT_NO_THROW( ValidateEulerDomainProblem( problem ) );
    EXPECT_NO_THROW( ValidateEulerDomainField( problem, input ) );
    EXPECT_NO_THROW( ValidateEulerDomainField( problem, output ) );
}

TEST( EulerDomainContract, ValidatesFiveEquationProblem )
{
    const EulerDomainProblem problem = FiveEquationProblem();
    Real values[ 5 * 32 ] = {};
    EulerDomainConstFieldView input{ 32, 5, values };

    EXPECT_NO_THROW( ValidateEulerDomainProblem( problem ) );
    EXPECT_NO_THROW( ValidateEulerDomainField( problem, input ) );
}

TEST( EulerDomainContract, RejectsInternalFieldShapeMismatch )
{
    const EulerDomainProblem problem = Problem();
    Real values[ 3 * 32 ] = {};
    EulerDomainConstFieldView wrongCells{ 31, 3, values };
    EulerDomainConstFieldView wrongEquations{ 32, 5, values };
    EulerDomainConstFieldView nullValues{ 32, 3, nullptr };

    EXPECT_THROW(
        ValidateEulerDomainField( problem, wrongCells ), std::invalid_argument );
    EXPECT_THROW(
        ValidateEulerDomainField( problem, wrongEquations ), std::invalid_argument );
    EXPECT_THROW(
        ValidateEulerDomainField( problem, nullValues ), std::invalid_argument );
}

TEST( EulerDomainContract, StateKeySeparatesExecutionOwnership )
{
    const EulerDomainStateKey cpu{ 0, 4, 0, AccelBackendKind::CPU };
    const EulerDomainStateKey same{ 0, 4, 0, AccelBackendKind::CPU };
    const EulerDomainStateKey otherBackend{ 0, 4, 0, AccelBackendKind::HIP };
    const EulerDomainStateKey otherZone{ 0, 5, 0, AccelBackendKind::CPU };

    EXPECT_TRUE( cpu == same );
    EXPECT_FALSE( cpu == otherBackend );
    EXPECT_FALSE( cpu == otherZone );
}

TEST( AccelViewsContract, DescribesLayoutsAreaPolicyAndCapabilities )
{
    const SolverDomainCapabilities threeEq{
        3, 2, 1, true, true, true, true };
    const SolverDomainCapabilities fiveEq{
        5, 4, 2, true, true, true, true };
    const SolverDomainCapabilities fourEq{
        4, 0, 0, true, true, false, false };

    EXPECT_TRUE( threeEq.SupportsEulerState() );
    EXPECT_TRUE( fiveEq.SupportsEulerState() );
    EXPECT_FALSE( fourEq.SupportsEulerState() );

    Real values[ 3 * 4 ] = {};
    const SolverConstFieldView field{
        4, 3, values, FieldLayout::EquationMajor,
        FieldRepresentation::Primitive };
    EXPECT_NO_THROW( ValidateSolverFieldView( field ) );

    FaceGeometryView geometry;
    geometry.nFaces = 4;
    geometry.faceArea = values;
    geometry.areaPolicy = FaceAreaPolicy::BackendMultiplies;
    EXPECT_NO_THROW( ValidateFaceGeometryView( geometry ) );
}

TEST( AccelViewsContract, RejectsUnsupportedStateAndGeometryContracts )
{
    Real values[ 4 ] = {};
    int cells[ 4 ] = {};
    FaceStateView state;
    state.nFaces = 4;
    state.nEquations = 3;
    state.qLeft = values;
    state.qRight = values;
    state.layout = FieldLayout::EntityMajor;
    EXPECT_THROW( ValidateFaceStateView( state ), std::invalid_argument );

    FaceConnectivityView connectivity;
    connectivity.nFaces = 4;
    connectivity.nBoundaryFaces = 5;
    connectivity.leftCell = nullptr;
    connectivity.rightCell = cells;
    EXPECT_THROW( ValidateFaceConnectivityView( connectivity ),
        std::invalid_argument );

    FaceGeometryView geometry;
    geometry.nFaces = 4;
    EXPECT_THROW( ValidateFaceGeometryView( geometry ), std::invalid_argument );
}

} // namespace
