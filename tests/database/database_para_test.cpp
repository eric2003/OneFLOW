#include <gtest/gtest.h>

#include <stdexcept>
#include <string>

// OneFLOW headers
#include "DataBase.h"
#include "DataBaseType.h"
#include "DataPara.h"
#include "DataObject.h"

using namespace ONEFLOW;


// ============================================================================
// Test Fixture
// ============================================================================

class DataParaTest : public ::testing::Test
{
protected:

    void SetUp() override
    {
        DataBaseType::Init();

        // Use the global database for consistency with existing database tests.
        db_ = GetGlobalDataBase();
    }

    void TearDown() override
    {
        // The global database is currently reused between tests.
        // Test keys should therefore be unique to avoid accidental collisions.
    }

    DataBase* db_ = nullptr;
};


// ----------------------------------------------------------------------------
// 1. Update an existing key with the same type and size
// ----------------------------------------------------------------------------

TEST_F(DataParaTest, UpdateExistingKeyWithSameTypeAndSize)
{
    // Create the initial integer value.
    SetDataInt( "test_value", 10 );

    ASSERT_NE(
        GetDataPointer< int >( "test_value" ),
        nullptr
    );

    EXPECT_EQ(
        GetDataValue< int >( "test_value" ),
        10
    );

    // Update the existing entry with the same type and size.
    SetDataInt( "test_value", 20 );

    // The existing value should be replaced successfully.
    EXPECT_EQ(
        GetDataValue< int >( "test_value" ),
        20
    );
}


// ----------------------------------------------------------------------------
// 2. Reject an update when the type is different
// ----------------------------------------------------------------------------

TEST_F(DataParaTest, RejectUpdateWithDifferentType)
{
    // Create the original integer entry.
    SetDataInt( "test_value", 10 );

    ASSERT_NE(
        GetDataPointer< int >( "test_value" ),
        nullptr
    );

    EXPECT_EQ(
        GetDataValue< int >( "test_value" ),
        10
    );

    // Construct an update entry with a different data type.
    DataEntry* dataEntry = new DataEntry();

    dataEntry->name = "test_value";
    dataEntry->type = HX_REAL;
    dataEntry->size = 1;

    Real value = 20.0;

    TDataObject< Real >* dataObject =
        new TDataObject< Real >( 1 );

    dataObject->CopyValue( &value );

    dataEntry->data = dataObject;

    // Updating an existing entry with a different type must fail.
    EXPECT_THROW(
        db_->dataPara->UpdateDataPointer( dataEntry ),
        std::runtime_error
    );

    // The original value must remain unchanged.
    EXPECT_EQ(
        GetDataValue< int >( "test_value" ),
        10
    );
}


// ----------------------------------------------------------------------------
// 3. Reject an update when the size is different
// ----------------------------------------------------------------------------

TEST_F(DataParaTest, RejectUpdateWithDifferentSize)
{
    // Create the original integer array with two elements.
    int initialValues[2] = { 10, 20 };

    SetData(
        "test_value",
        initialValues,
        HX_INT,
        2
    );

    DataEntry* existing =
        db_->dataPara->GetDataPointer( "test_value" );

    ASSERT_NE( existing, nullptr );

    EXPECT_EQ( existing->type, HX_INT );
    EXPECT_EQ( existing->size, 2 );

    // Construct an update entry with the same type but a different size.
    DataEntry* dataEntry = new DataEntry();

    dataEntry->name = "test_value";
    dataEntry->type = HX_INT;
    dataEntry->size = 3;

    int newValues[3] = { 30, 40, 50 };

    TDataObject< int >* dataObject =
        new TDataObject< int >( 3 );

    dataObject->CopyValue( newValues );

    dataEntry->data = dataObject;

    // Updating an existing entry with a different size must fail.
    EXPECT_THROW(
        db_->dataPara->UpdateDataPointer( dataEntry ),
        std::runtime_error
    );

    // Verify that the original entry was not modified.
    existing =
        db_->dataPara->GetDataPointer( "test_value" );

    ASSERT_NE( existing, nullptr );

    EXPECT_EQ( existing->type, HX_INT );
    EXPECT_EQ( existing->size, 2 );

    // Verify that the original data is still intact.
    int* values =
        GetDataPointer< int >( "test_value" );

    ASSERT_NE( values, nullptr );

    EXPECT_EQ( values[0], 10 );
    EXPECT_EQ( values[1], 20 );
}