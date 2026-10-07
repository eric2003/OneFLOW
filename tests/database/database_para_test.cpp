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
    Real value = 20.0;
    auto dataObject = std::make_unique<TDataObject<Real>>( 1 );
    dataObject->CopyValue( &value, 1 );
    auto dataEntry = std::make_unique<DataEntry>(
        "test_value", HX_REAL, 1, std::move( dataObject ) );

    // Updating an existing entry with a different type must fail.
    EXPECT_THROW(
        db_->GetDataPara()->SetDataEntry( std::move( dataEntry ) ),
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
        db_->GetDataPara()->FindDataEntry( "test_value" );

    ASSERT_NE( existing, nullptr );

    EXPECT_EQ( existing->GetType(), HX_INT );
    EXPECT_EQ( existing->GetSize(), 2 );

    // Construct an update entry with the same type but a different size.
    int newValues[3] = { 30, 40, 50 };
    auto dataObject = std::make_unique<TDataObject< int >>( 3 );
    dataObject->CopyValue( newValues, 3 );
    auto dataEntry = std::make_unique<DataEntry>(
        "test_value", HX_INT, 3, std::move( dataObject ) );

    // Updating an existing entry with a different size must fail.
    EXPECT_THROW(
        db_->GetDataPara()->SetDataEntry( std::move( dataEntry ) ),
        std::runtime_error
    );

    // Verify that the original entry was not modified.
    existing =
        db_->GetDataPara()->FindDataEntry( "test_value" );

    ASSERT_NE( existing, nullptr );

    EXPECT_EQ( existing->GetType(), HX_INT );
    EXPECT_EQ( existing->GetSize(), 2 );

    // Verify that the original data is still intact.
    int* values =
        GetDataPointer< int >( "test_value" );

    ASSERT_NE( values, nullptr );

    EXPECT_EQ( values[0], 10 );
    EXPECT_EQ( values[1], 20 );
}

// ----------------------------------------------------------------------------
// 4. Clear all entries without destroying the DataPara object
// ----------------------------------------------------------------------------

TEST(DataParaTestStandalone, ClearReleasesAllEntries)
{
    DataPara dataPara;

    int value = 42;
    auto dataObject = std::make_unique<TDataObject< int >>( 1 );
    dataObject->CopyValue( &value, 1 );
    auto dataEntry = std::make_unique<DataEntry>(
        "clear_value", HX_INT, 1, std::move( dataObject ) );

    dataPara.SetDataEntry( std::move( dataEntry ) );

    ASSERT_NE(
        dataPara.FindDataEntry( "clear_value" ),
        nullptr
    );

    dataPara.Clear();

    EXPECT_EQ(
        dataPara.FindDataEntry( "clear_value" ),
        nullptr
    );
}

TEST(DataBookOwnership, UniquePtrAndRawViewShareObject)
{
    auto ownedBook = std::make_unique<DataBook>();
    DataBook * dataBook = ownedBook.get();

    ASSERT_EQ( dataBook, ownedBook.get() );

    // Simulate compress: write through the raw view
    int n = 42;
    ONEFLOW::HXWrite( dataBook, n );

    dataBook->MoveToBegin();
    int out = 0;
    ONEFLOW::HXRead( dataBook, out );
    EXPECT_EQ( out, 42 );

    // Rebinding the local pointer must not change ownership
    DataBook * other = nullptr;
    dataBook = other;
    EXPECT_NE( ownedBook.get(), nullptr );  // still owns original
}

TEST(DataEntryInvariant, DataObjectIsRequired)
{
    int value = 42;
    auto dataObject = std::make_unique<TDataObject< int >>( 1 );
    dataObject->CopyValue( &value, 1 );
    DataEntry dataEntry( "required_data", HX_INT, 1, std::move( dataObject ) );

    DataObject & object = dataEntry.GetDataObject();
    EXPECT_EQ( object.GetVoidPointer(), dataEntry.GetDataObject().GetVoidPointer() );
    EXPECT_EQ( GetDataValue< int >( &object ), 42 );
}

TEST(DataEntryInvariant, NullDataObjectIsRejected)
{
    EXPECT_THROW(
        DataEntry( "invalid_data", HX_INT, 1, nullptr ),
        std::invalid_argument
    );
}
