#include <gtest/gtest.h>
#include "ZoneSharedFields.h"
#include <cstring>

using namespace ONEFLOW;

TEST( ZoneSharedFields, NameLiteralsAreStable )
{
    EXPECT_STREQ( ZoneSharedFieldNames::ViscosityLaminar, "visl" );
    EXPECT_STREQ( ZoneSharedFieldNames::ViscosityTurbulent, "vist" );
    EXPECT_STREQ( ZoneSharedFieldNames::TimeStep, "timestep" );
    EXPECT_STREQ( ZoneSharedFieldNames::Gamma, "gama" );
}

TEST( ZoneSharedFields, NamesAreNotSolverPrefixed )
{
    EXPECT_EQ( std::strstr( ZoneSharedFieldNames::ViscosityLaminar, "turb" ), nullptr );
    EXPECT_EQ( std::strstr( ZoneSharedFieldNames::TimeStep, "turb" ), nullptr );
}

TEST( ZoneSharedFields, DefaultViewPointersAreNull )
{
    ZoneSharedFieldView view;
    EXPECT_EQ( view.visl, nullptr );
    EXPECT_EQ( view.vist, nullptr );
    EXPECT_EQ( view.timestep, nullptr );
    EXPECT_EQ( view.gama, nullptr );
}
