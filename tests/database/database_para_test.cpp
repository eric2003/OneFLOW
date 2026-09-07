#include <gtest/gtest.h>
#include <string>
#include <stdexcept>

// OneFLOW headers (ensure include path is correct in CMake)
#include "DataBase.h"
#include "DataPara.h"
#include "DataBaseType.h"
#include "HXType.h"          // for Real

using namespace ONEFLOW;

// ============================================================================
// Test Fixture: ensures clean global DataBase for every test
// ============================================================================
class DataParaTest : public ::testing::Test {
protected:
    void SetUp() override {
        // Make sure type system is initialized
        DataBaseType::Init();

        // Clear any previous data (important because of global singleton)
        DataBase* db = GetGlobalDataBase();
        // For short-term we keep the global, but we can delete all entries
        // In future we will inject a fresh DataBase instance
    }

    void TearDown() override {
        // Optional: clean up after each test to avoid pollution
    }
};

// ----------------------------------------------------------------------------
// 1. Basic Set / Get for int
// ----------------------------------------------------------------------------
TEST_F(DataParaTest, SetAndGetInt) {
    int value = 42;
    SetDataInt("test_int", value);

    int retrieved = GetDataValue<int>("test_int");
    EXPECT_EQ(retrieved, 42);
}

// ----------------------------------------------------------------------------
// 2. Basic Set / Get for Real
// ----------------------------------------------------------------------------
TEST_F(DataParaTest, SetAndGetReal) {
    Real value = 3.14159;
    SetDataReal("test_real", value);

    Real retrieved = GetDataValue<Real>("test_real");
    EXPECT_DOUBLE_EQ(retrieved, 3.14159);
}

// ----------------------------------------------------------------------------
// 3. Fixed SetDataString (this is the bug we fixed in step 1)
// ----------------------------------------------------------------------------
TEST_F(DataParaTest, SetAndGetString) {
    std::string value = "OneFLOW_Database_Test";
    SetDataString("test_string", value);

    std::string retrieved = GetDataValue<std::string>("test_string");
    EXPECT_EQ(retrieved, "OneFLOW_Database_Test");
}

// ----------------------------------------------------------------------------
// 4. Update existing key (Copy semantics of UpdateDataPointer)
// ----------------------------------------------------------------------------
TEST_F(DataParaTest, UpdateExistingKey) {
    int v1 = 10;
    SetDataInt("update_key", v1);
    EXPECT_EQ(GetDataValue<int>("update_key"), 10);

    int v2 = 99;
    SetDataInt("update_key", v2);   // should overwrite
    EXPECT_EQ(GetDataValue<int>("update_key"), 99);
}

// ----------------------------------------------------------------------------
// 5. Multiple types coexist
// ----------------------------------------------------------------------------
TEST_F(DataParaTest, MultipleTypesCoexist) {
    SetDataInt("multi_int", 7);
    SetDataReal("multi_real", 2.718);
    SetDataString("multi_str", "hello");

    EXPECT_EQ(GetDataValue<int>("multi_int"), 7);
    EXPECT_DOUBLE_EQ(GetDataValue<Real>("multi_real"), 2.718);
    EXPECT_EQ(GetDataValue<std::string>("multi_str"), "hello");
}

// ----------------------------------------------------------------------------
// 6. Delete key
// ----------------------------------------------------------------------------
TEST_F(DataParaTest, DeleteKey) {
    SetDataInt("to_delete", 123);
    EXPECT_EQ(GetDataValue<int>("to_delete"), 123);

    DataBase* db = GetGlobalDataBase();
    db->dataPara->DeleteDataPointer("to_delete");

    // After delete, Get should fail (we will test the exception in another file)
    EXPECT_THROW(GetDataValue<int>("to_delete"), std::runtime_error);
}