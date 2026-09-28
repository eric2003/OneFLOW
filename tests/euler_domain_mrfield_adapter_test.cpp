#include "EulerDomainMrFieldAdapter.h"

#include <gtest/gtest.h>

#include <stdexcept>

namespace
{

using namespace ONEFLOW;

TEST( EulerDomainMrFieldAdapter, PacksInternalCellsEquationMajor )
{
    MRField field( 5, 4 );
    for ( int equation = 0; equation < 5; ++ equation )
    {
        for ( int cell = 0; cell < 4; ++ cell )
        {
            field[ equation ][ cell ] = 100.0 * equation + cell;
        }
    }

    EulerDomainMrFieldSnapshot snapshot( field, 3 );
    EulerDomainConstFieldView view = snapshot.View();

    ASSERT_EQ( view.nCells, 3 );
    ASSERT_EQ( view.nEquations, 5 );
    ASSERT_NE( view.values, nullptr );
    for ( int equation = 0; equation < 5; ++ equation )
    {
        for ( int cell = 0; cell < 3; ++ cell )
        {
            EXPECT_EQ( view.values[ equation * 3 + cell ],
                100.0 * equation + cell );
        }
    }
}

TEST( EulerDomainMrFieldAdapter, RejectsUnsupportedShape )
{
    MRField fourEquationField( 4, 3 );
    MRField shortField( 5, 2 );

    EXPECT_THROW(
        EulerDomainMrFieldSnapshot( fourEquationField, 3 ),
        std::invalid_argument );
    EXPECT_THROW(
        EulerDomainMrFieldSnapshot( shortField, 3 ),
        std::invalid_argument );
}

} // namespace
