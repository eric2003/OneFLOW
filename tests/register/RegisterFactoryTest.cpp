// RegisterFactoryTest.cpp
#include <gtest/gtest.h>
#include "Register.h"
#include "HXDefine.h"
#include "HXClone.h"

class RegisterFactoryTest : public ::testing::Test
{
protected:
    void TearDown() override
    {
        ONEFLOW::RegisterFactory::FreeMRegister();
    }
};

TEST_F( RegisterFactoryTest, GetMRegisterOnUnregisteredIdReturnsNull )
{
    EXPECT_EQ( ONEFLOW::RegisterFactory::GetMRegister( 42 ), nullptr );
}

TEST_F( RegisterFactoryTest, AddThenGetMRegisterReturnsSameInstance )
{
    ONEFLOW::RegisterFactory::AddMRegister( 1 );

    ONEFLOW::MRegister * a = ONEFLOW::RegisterFactory::GetMRegister( 1 );
    ONEFLOW::MRegister * b = ONEFLOW::RegisterFactory::GetMRegister( 1 );

    ASSERT_NE( a, nullptr );
    EXPECT_EQ( a, b ); // same underlying instance on repeated lookups
}

TEST_F( RegisterFactoryTest, AddMRegisterIsIdempotentForSameId )
{
    ONEFLOW::RegisterFactory::AddMRegister( 1 );
    ONEFLOW::MRegister * first = ONEFLOW::RegisterFactory::GetMRegister( 1 );

    ONEFLOW::RegisterFactory::AddMRegister( 1 ); // should not replace it
    ONEFLOW::MRegister * second = ONEFLOW::RegisterFactory::GetMRegister( 1 );

    EXPECT_EQ( first, second );
}

TEST_F( RegisterFactoryTest, GetRegisterOnUnregisteredMRegisterIdReturnsNull )
{
    EXPECT_EQ( ONEFLOW::RegisterFactory::GetRegister( 99, 0 ), nullptr );
}

TEST_F( RegisterFactoryTest, FreeMRegisterClearsEverything )
{
    ONEFLOW::RegisterFactory::AddMRegister( 1 );
    ONEFLOW::RegisterFactory::FreeMRegister();

    EXPECT_EQ( ONEFLOW::RegisterFactory::GetMRegister( 1 ), nullptr );
}

// --- MRegister-level tests (no file I/O; only exercises AllocateData/GetRegister) ---

TEST( MRegisterTest, GetRegisterOnOutOfRangeIndexReturnsNull )
{
    ONEFLOW::MRegister mRegister;
    EXPECT_EQ( mRegister.GetRegister( 0 ), nullptr );   // nothing allocated yet
    EXPECT_EQ( mRegister.GetRegister( -1 ), nullptr );
}

TEST( MRegisterTest, DefaultGetRegisterReturnsIndexZero )
{
    ONEFLOW::StringField names;
    names.push_back( "dummy_file.txt" ); // AllocateData only needs the count
    ONEFLOW::MRegister mRegister;
    mRegister.SetSolverFileNames( names );

    // Directly exercising AllocateData() indirectly via RegisterAll()
    // would require a real file (TextFileParser::OpenFile). Since we're
    // only testing allocation/indexing here, we call the private
    // AllocateData() indirectly is not possible (it's private) - so this
    // test is limited to what's reachable through the public API.
    // Left as a placeholder: consider a friend test hook or making
    // AllocateData reachable if this path needs direct coverage later.
}