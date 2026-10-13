#include <gtest/gtest.h>
#include <stdexcept>

#include "ElementHome.h"
#include "HXCgns.h"

TEST(ElementHomeTest, InitializeBuildsRegistryWithStableEnumIndices)
{
    ONEFLOW::ElementHome::Initialize();

    EXPECT_GT(ONEFLOW::ElementHome::GetNumberOfUnitElements(), ONEFLOW::HEXA_8);
    EXPECT_TRUE(ONEFLOW::ElementHome::HasUnitElement(ONEFLOW::NODE));
    EXPECT_TRUE(ONEFLOW::ElementHome::HasUnitElement(ONEFLOW::BAR_2));
    EXPECT_TRUE(ONEFLOW::ElementHome::HasUnitElement(ONEFLOW::HEXA_8));

    EXPECT_EQ(ONEFLOW::GetElementNodeNumbers(ONEFLOW::NODE), 1);
    EXPECT_EQ(ONEFLOW::GetElementNodeNumbers(ONEFLOW::BAR_2), 2);
    EXPECT_EQ(ONEFLOW::GetElementNodeNumbers(ONEFLOW::HEXA_8), 8);

    EXPECT_EQ(ONEFLOW::GetElementDimension(ONEFLOW::NODE), 0);
    EXPECT_EQ(ONEFLOW::GetElementDimension(ONEFLOW::BAR_2), 1);
    EXPECT_EQ(ONEFLOW::GetElementDimension(ONEFLOW::HEXA_8), 3);
}

TEST(ElementHomeTest, ReadOnlyDefinitionAccessorPreservesRegistryTopology)
{
    ONEFLOW::ElementHome::Initialize();

    const ONEFLOW::UnitElement & element =
        ONEFLOW::ElementHome::GetUnitElementDefinition(ONEFLOW::HEXA_8);

    EXPECT_EQ(element.GetElementType(), ONEFLOW::HEXA_8);
    EXPECT_EQ(element.GetElementNodeNumbers(ONEFLOW::HEXA_8), 8);
    EXPECT_EQ(element.GetElementFaceNumber(), 6);
}

TEST(ElementHomeTest, RepeatedInitializeRebuildsEquivalentRegistry)
{
    ONEFLOW::ElementHome::Initialize();
    const int initialCount = ONEFLOW::ElementHome::GetNumberOfUnitElements();
    ASSERT_TRUE(ONEFLOW::ElementHome::HasUnitElement(ONEFLOW::HEXA_8));
    const int initialNodeCount =
        ONEFLOW::GetElementNodeNumbers(ONEFLOW::HEXA_8);

    ONEFLOW::ElementHome::Initialize();

    EXPECT_EQ(ONEFLOW::ElementHome::GetNumberOfUnitElements(), initialCount);
    ASSERT_TRUE(ONEFLOW::ElementHome::HasUnitElement(ONEFLOW::HEXA_8));
    EXPECT_EQ(ONEFLOW::GetElementNodeNumbers(ONEFLOW::HEXA_8), initialNodeCount);
}

TEST(ElementHomeTest, InvalidElementIndicesAreRejectedConsistently)
{
    EXPECT_FALSE(ONEFLOW::ElementHome::HasUnitElement(-1));
    EXPECT_FALSE(
        ONEFLOW::ElementHome::HasUnitElement(
            ONEFLOW::ElementHome::GetNumberOfUnitElements()));

    EXPECT_THROW(ONEFLOW::ElementHome::GetUnitElement(-1), std::out_of_range);
    EXPECT_THROW(
        ONEFLOW::GetElementNodeNumbers(
            ONEFLOW::ElementHome::GetNumberOfUnitElements()),
        std::out_of_range);
}

TEST(ElementHomeTest, FreeClearsRegistryAndInitializeRestoresIt)
{
    ONEFLOW::ElementHome::Initialize();
    ASSERT_TRUE(ONEFLOW::ElementHome::HasUnitElement(ONEFLOW::HEXA_8));

    ONEFLOW::ElementHome::Free();

    EXPECT_EQ(ONEFLOW::ElementHome::GetNumberOfUnitElements(), 0);
    EXPECT_FALSE(ONEFLOW::ElementHome::HasUnitElement(ONEFLOW::HEXA_8));
    EXPECT_THROW(ONEFLOW::ElementHome::GetUnitElement(ONEFLOW::HEXA_8),
                 std::out_of_range);

    ONEFLOW::ElementHome::Initialize();

    EXPECT_TRUE(ONEFLOW::ElementHome::HasUnitElement(ONEFLOW::HEXA_8));
    EXPECT_EQ(ONEFLOW::GetElementNodeNumbers(ONEFLOW::HEXA_8), 8);
}

TEST(UnitElementTest, FailedInitializePreservesExistingTopology)
{
    ONEFLOW::UnitElement element;
    element.Initialize(ONEFLOW::QUAD_4);

    const int initialElementType = element.GetElementType();
    const int initialFaceCount = element.GetElementFaceNumber();
    const int initialNodeCount = element.GetElementNodeNumbers(ONEFLOW::QUAD_4);

    EXPECT_THROW(element.Initialize(-1), std::invalid_argument);

    EXPECT_EQ(element.GetElementType(), initialElementType);
    EXPECT_EQ(element.GetElementFaceNumber(), initialFaceCount);
    EXPECT_EQ(element.GetElementNodeNumbers(ONEFLOW::QUAD_4), initialNodeCount);
    ASSERT_EQ(element.GetElementFaceNumber(), 4);
    EXPECT_EQ(element.GetElementFace(0).size(), 2);
}
