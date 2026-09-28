#include <gtest/gtest.h>
#include <stdexcept>
#include "Solver.h"

namespace
{
    class StubSolver : public ONEFLOW::Solver
    {
    public:
        ONEFLOW::Solver * Clone() const override
        {
            return new StubSolver( *this );
        }
    };
}

TEST( SolverTest, SafeCloneOnUnregisteredTypeThrows )
{
    EXPECT_THROW(
        ONEFLOW::Solver::SafeClone( "SolverTest_NeverRegistered" ),
        std::runtime_error
    );
}

TEST( SolverTest, RegisterThenSafeCloneReturnsANewInstance )
{
    ONEFLOW::Solver::Register( "SolverTest_TypeA", new StubSolver() );

    ONEFLOW::Solver * cloned = ONEFLOW::Solver::SafeClone( "SolverTest_TypeA" );

    ASSERT_NE( cloned, nullptr );
    delete cloned;
}