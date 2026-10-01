#include <gtest/gtest.h>
#include <memory>
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
        ONEFLOW::Solver::SafeCloneUnique( "SolverTest_NeverRegistered" ),
        std::runtime_error
    );
}

TEST( SolverTest, RegisterThenSafeCloneReturnsANewInstance )
{
    ONEFLOW::Solver::Register( "SolverTest_TypeA", std::make_unique< StubSolver >() );

    std::unique_ptr< ONEFLOW::Solver > cloned =
        ONEFLOW::Solver::SafeCloneUnique( "SolverTest_TypeA" );

    ASSERT_NE( cloned, nullptr );
}
