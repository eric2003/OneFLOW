#include <gtest/gtest.h>
#include "ZoneSharedFields.h"
#include "IFieldLookup.h"
#include <cstring>
#include <string>
#include <vector>

using namespace ONEFLOW;

namespace
{
    class FakeFieldLookup : public IFieldLookup
    {
    public:
        std::vector< std::string > requested;
        MRField vislStub;
        MRField vistStub;
        MRField timestepStub;
        MRField gamaStub;

        MRField * GetMRField( Grid * /*grid*/, const std::string & name ) override
        {
            requested.push_back( name );
            if ( name == ZoneSharedFieldNames::ViscosityLaminar ) return &vislStub;
            if ( name == ZoneSharedFieldNames::ViscosityTurbulent ) return &vistStub;
            if ( name == ZoneSharedFieldNames::TimeStep ) return &timestepStub;
            if ( name == ZoneSharedFieldNames::Gamma ) return &gamaStub;
            return nullptr;
        }
    };
}

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

TEST( ZoneSharedFields, Bind_UsesLookupPortForAllSharedNames )
{
    FakeFieldLookup fake;
    ZoneSharedFieldView view;
    view.Bind( nullptr, fake );

    ASSERT_EQ( fake.requested.size(), 4u );
    EXPECT_EQ( fake.requested[ 0 ], "visl" );
    EXPECT_EQ( fake.requested[ 1 ], "vist" );
    EXPECT_EQ( fake.requested[ 2 ], "timestep" );
    EXPECT_EQ( fake.requested[ 3 ], "gama" );

    EXPECT_EQ( view.visl, &fake.vislStub );
    EXPECT_EQ( view.vist, &fake.vistStub );
    EXPECT_EQ( view.timestep, &fake.timestepStub );
    EXPECT_EQ( view.gama, &fake.gamaStub );
}
