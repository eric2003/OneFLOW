#include <gtest/gtest.h>
#include <string>
#include <memory>

// OneFLOW headers
#include "DataBase.h"
#include "DataField.h"
#include "DataPointer.h"
#include "DataBaseType.h"

using namespace ONEFLOW;

// Simple POD type for testing field data
struct DummyField
{
    int    id    = 0;
    double value = 0.0;
};

// ============================================================================
// Test Fixture
// ============================================================================
class DataFieldTest : public ::testing::Test
{
protected:
    void SetUp() override
    {
        DataBaseType::Init();
        // Use the global database for now (consistent with DataPara tests)
        db_ = GetGlobalDataBase();
    }

    void TearDown() override
    {
        // Clean up any leftover fields to avoid pollution between tests
        // (In the future we will inject a fresh DataBase)
    }

    DataBase* db_ = nullptr;
};

// ----------------------------------------------------------------------------
// 1. Create and retrieve a field pointer successfully
// ----------------------------------------------------------------------------
TEST_F(DataFieldTest, CreateAndRetrievePointer)
{
    // Create a DataPointer that owns a DummyField
    auto* rawField = new DummyField{ 42, 3.14159 };
    auto* wrap = new DataPointer<DummyField>( rawField );

    CreateFieldPointer( db_, wrap, "dummy_field" );

    // Retrieve via GetFieldPointer
    DummyField* ptr = GetFieldPointer<DummyField>( db_, "dummy_field" );

    ASSERT_NE( ptr, nullptr );
    EXPECT_EQ( ptr->id, 42 );
    EXPECT_DOUBLE_EQ( ptr->value, 3.14159 );
}

// ----------------------------------------------------------------------------
// 2. GetFieldReference should give a usable reference
// ----------------------------------------------------------------------------
TEST_F(DataFieldTest, GetFieldReference)
{
    auto* rawField = new DummyField{ 7, 2.718 };
    auto* wrap = new DataPointer<DummyField>( rawField );

    CreateFieldPointer( db_, wrap, "ref_field" );

    DummyField& ref = GetFieldReference<DummyField>( db_, "ref_field" );

    EXPECT_EQ( ref.id, 7 );
    EXPECT_DOUBLE_EQ( ref.value, 2.718 );

    // Modify through reference and verify
    ref.id = 99;
    DummyField* ptr = GetFieldPointer<DummyField>( db_, "ref_field" );
    EXPECT_EQ( ptr->id, 99 );
}

// ----------------------------------------------------------------------------
// 3. Query a non-existent field -> should return nullptr
// ----------------------------------------------------------------------------
TEST_F(DataFieldTest, GetNonExistentFieldReturnsNull)
{
    DummyField* ptr = GetFieldPointer<DummyField>( db_, "this_field_does_not_exist" );
    EXPECT_EQ( ptr, nullptr );
}

// ----------------------------------------------------------------------------
// 4. Delete a field and verify it is gone
// ----------------------------------------------------------------------------
TEST_F(DataFieldTest, DeleteField)
{
    auto* rawField = new DummyField{ 123, 1.0 };
    auto* wrap = new DataPointer<DummyField>( rawField );

    CreateFieldPointer( db_, wrap, "to_delete_field" );

    // Confirm it exists
    ASSERT_NE( GetFieldPointer<DummyField>( db_, "to_delete_field" ), nullptr );

    // Delete
    db_->dataField->DeleteDataF( "to_delete_field" );

    // Should now be gone
    EXPECT_EQ( GetFieldPointer<DummyField>( db_, "to_delete_field" ), nullptr );
}

// ----------------------------------------------------------------------------
// 5. Destructor safety (critical test)
//    If the old reinterpret_cast bug still exists, AddressSanitizer
//    or the test runner will crash / report error here.
// ----------------------------------------------------------------------------
TEST_F(DataFieldTest, DestructorSafety)
{
    // Create a local DataBase to control lifetime precisely
    DataBase localDb;

    {
        auto* rawField = new DummyField{ 1, 1.0 };
        auto* wrap = new DataPointer<DummyField>( rawField );

        CreateFieldPointer( &localDb, wrap, "temp_field" );

        DummyField* p = GetFieldPointer<DummyField>( &localDb, "temp_field" );
        ASSERT_NE( p, nullptr );
        EXPECT_EQ( p->id, 1 );
    } // wrap and DataF are owned by localDb.dataField

      // When localDb goes out of scope, DataField destructor runs.
      // This test passes if and only if there is no crash / double-free /
      // invalid delete.
}

// ----------------------------------------------------------------------------
// 6. Multiple fields coexist
// ----------------------------------------------------------------------------
TEST_F(DataFieldTest, MultipleFieldsCoexist)
{
    auto* f1 = new DummyField{ 10, 1.1 };
    auto* f2 = new DummyField{ 20, 2.2 };

    CreateFieldPointer( db_, new DataPointer<DummyField>( f1 ), "field_a" );
    CreateFieldPointer( db_, new DataPointer<DummyField>( f2 ), "field_b" );

    DummyField* p1 = GetFieldPointer<DummyField>( db_, "field_a" );
    DummyField* p2 = GetFieldPointer<DummyField>( db_, "field_b" );

    ASSERT_NE( p1, nullptr );
    ASSERT_NE( p2, nullptr );

    EXPECT_EQ( p1->id, 10 );
    EXPECT_EQ( p2->id, 20 );
    EXPECT_NE( p1, p2 );
}