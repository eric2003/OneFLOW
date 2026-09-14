#include "EulerDomain.h"

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

} // namespace
