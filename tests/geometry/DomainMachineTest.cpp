#include "DomainMachine.h"

#include <gtest/gtest.h>

#include <stdexcept>

namespace ONEFLOW
{

TEST( DomainMachineTest, PreservesBoundaryIds )
{
    DomainMachine machine;

    machine.SetBcType( 4, 44 );
    machine.SetBcType( 2, 22 );

    EXPECT_EQ( machine.GetBcType( 2 ), 22 );
    EXPECT_EQ( machine.GetBcType( 4 ), 44 );
}

TEST( DomainMachineTest, RejectsInvalidBoundaryId )
{
    DomainMachine machine;

    EXPECT_THROW( machine.SetBcType( 0, 1 ), std::invalid_argument );
}

} // namespace ONEFLOW
