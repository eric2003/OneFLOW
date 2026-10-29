#include <gtest/gtest.h>
#include <type_traits>

#include "Marray.h"

namespace
{
using ONEFLOW::HXSize_t;
using ONEFLOW::HXVector;
using ONEFLOW::Marray;

using MutableEquationAccess = HXVector< int > & ( Marray< int >::* )( HXSize_t );
using ConstEquationAccess = const HXVector< int > & ( Marray< int >::* )( HXSize_t ) const;

static_assert( std::is_same_v<
    decltype( static_cast< MutableEquationAccess >( &Marray< int >::operator[] ) ),
    MutableEquationAccess > );
static_assert( std::is_same_v<
    decltype( static_cast< ConstEquationAccess >( &Marray< int >::operator[] ) ),
    ConstEquationAccess > );

TEST(MarrayTest, ConstructorPreservesBothSizeTypeExtents)
{
    const HXSize_t equationCount = 2;
    const HXSize_t valueCount = 3;

    Marray< int > values( equationCount, valueCount );

    EXPECT_EQ(values.GetNEqu(), equationCount);
    ASSERT_EQ(values[ 0 ].size(), valueCount);
    ASSERT_EQ(values[ 1 ].size(), valueCount);
}
}
