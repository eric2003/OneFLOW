#include <gtest/gtest.h>
#include <type_traits>
#include <utility>

#include "Marray.h"

namespace
{
using ONEFLOW::HXSize_t;
using ONEFLOW::HXVector;
using ONEFLOW::Marray;

static_assert( std::is_same_v<
    decltype( std::declval< Marray< int > & >()[ std::declval< HXSize_t >() ] ),
    HXVector< int > & > );
static_assert( std::is_same_v<
    decltype( std::declval< const Marray< int > & >()[ std::declval< HXSize_t >() ] ),
    const HXVector< int > & > );

TEST( MarrayTest, ConstructorPreservesBothSizeTypeExtents )
{
    const HXSize_t equationCount = 2;
    const HXSize_t valueCount = 3;

    Marray< int > values( equationCount, valueCount );

    EXPECT_EQ( values.GetNEqu(), equationCount );
    ASSERT_EQ( values[ 0 ].size(), valueCount );
    ASSERT_EQ( values[ 1 ].size(), valueCount );
}

TEST( MarrayTest, MutableAndConstEquationAccessPreserveReferenceSemantics )
{
    Marray< int > values( 2, 3 );
    values[ 1 ][ 2 ] = 42;

    const Marray< int > & readOnlyValues = values;

    EXPECT_EQ( readOnlyValues[ 1 ][ 2 ], 42 );
}
}
