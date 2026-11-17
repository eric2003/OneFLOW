#include <gtest/gtest.h>

#include "DataBaseType.h"

using namespace ONEFLOW;

TEST( DataBaseTypeInitTest, AccessorsInitializeMappingsOnFirstUse )
{
    ASSERT_FALSE( DataBaseType::init_flag );

    EXPECT_EQ( DataBaseType::GetName( HX_INT ), "int" );
    EXPECT_EQ( DataBaseType::GetIndex( "string" ), HX_STRING );
}
