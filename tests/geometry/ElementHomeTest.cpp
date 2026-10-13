#include <gtest/gtest.h>
#include <stdexcept>
#include <type_traits>
#include <utility>

#include "ElementHome.h"
#include "HXCgns.h"

static_assert(
    std::is_const_v<std::remove_reference_t<
        decltype(ONEFLOW::ElementHome::GetUnitElementDefinition(ONEFLOW::HEXA_8))>>,
    "ElementHome must expose unit element definitions as read-only");

using ConstElementFaceResult = decltype(
    std::declval< const ONEFLOW::UnitElement & >().GetElementFace(0));
using ConstPhysicsFaceResult = decltype(
    std::declval< const ONEFLOW::UnitElement & >().GetElementPhysicsFace(0));
using ConstMiddlePointResult = decltype(
    std::declval< const ONEFLOW::UnitElement & >()
        .GetRelatedPointListForMiddlePointCalculation(0));
using ConstChildNodeResult = decltype(
    std::declval< const ONEFLOW::UnitElement & >()
        .GetChildElementRelativeNodeIndex(0));

static_assert(
    std::is_const_v<std::remove_reference_t<ConstElementFaceResult>>,
    "Reading a face from a const unit element must not expose mutable topology");
static_assert(
    std::is_const_v<std::remove_reference_t<ConstPhysicsFaceResult>>,
    "Reading a physics face from a const unit element must not expose mutable topology");
static_assert(
    std::is_const_v<std::remove_reference_t<ConstMiddlePointResult>>,
    "Reading middle-point relations from a const unit element must be read-only");
static_assert(
    std::is_const_v<std::remove_reference_t<ConstChildNodeResult>>,
    "Reading child-element node mappings from a const unit element must be read-only");

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

    EXPECT_THROW(
        ONEFLOW::ElementHome::GetUnitElementDefinition(-1),
        std::out_of_range);
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
    EXPECT_THROW(
        ONEFLOW::ElementHome::GetUnitElementDefinition(ONEFLOW::HEXA_8),
        std::out_of_range);

    ONEFLOW::ElementHome::Initialize();

    EXPECT_TRUE(ONEFLOW::ElementHome::HasUnitElement(ONEFLOW::HEXA_8));
    EXPECT_EQ(ONEFLOW::GetElementNodeNumbers(ONEFLOW::HEXA_8), 8);
}

TEST(UnitElementTest, DefaultConstructionStartsUninitialized)
{
    const ONEFLOW::UnitElement element;

    EXPECT_EQ(element.GetElementType(), -1);
    EXPECT_EQ(element.GetElementFaceNumber(), 0);
    EXPECT_EQ(element.GetElementPhysicsFaceNumber(), 0);
    EXPECT_EQ(element.GetChildElementNumbers(), 0);
}

TEST(UnitElementTest, CorrectlySpelledMiddlePointAccessorIsAvailable)
{
    ONEFLOW::UnitElement element;
    element.Initialize(ONEFLOW::QUAD_4);

    EXPECT_THROW(
        element.GetRelatedPointListForMiddlePointCalculation(-1),
        std::out_of_range);

    const ONEFLOW::UnitElement & constElement = element;
    EXPECT_NO_THROW(
        constElement.GetRelatedPointListForMiddlePointCalculation(0));
    EXPECT_TRUE(
        constElement.GetRelatedPointListForMiddlePointCalculation(0).empty());
    EXPECT_THROW(
        constElement.GetRelatedPointListForMiddlePointCalculation(4),
        std::out_of_range);
}

TEST(UnitElementTest, TopologyAccessorsRejectOutOfRangeIndices)
{
    ONEFLOW::UnitElement element;
    element.Initialize(ONEFLOW::QUAD_4);

    EXPECT_THROW(element.GetElementFace(-1), std::out_of_range);
    EXPECT_THROW(element.GetElementFace(4), std::out_of_range);
    EXPECT_THROW(element.GetFaceType(-1), std::out_of_range);
    EXPECT_THROW(element.GetFaceType(4), std::out_of_range);

    EXPECT_THROW(element.GetRelatedPointListForMiddlePointCalculation(-1),
                 std::out_of_range);
    EXPECT_THROW(element.GetRelatedPointListForMiddlePointCalculation(4),
                 std::out_of_range);

    EXPECT_THROW(element.GetElementPhysicsFace(-1), std::out_of_range);
    EXPECT_THROW(element.GetElementPhysicsFace(0), std::out_of_range);
    EXPECT_THROW(element.GetElementPhysicsFaceType(-1), std::out_of_range);
    EXPECT_THROW(element.GetElementPhysicsFaceType(0), std::out_of_range);

    EXPECT_THROW(element.GetChildElementRelativeNodeIndex(-1),
                 std::out_of_range);
    EXPECT_THROW(element.GetChildElementRelativeNodeIndex(0),
                 std::out_of_range);
    EXPECT_THROW(element.GetChildElementType(-1), std::out_of_range);
    EXPECT_THROW(element.GetChildElementType(0), std::out_of_range);
}

TEST(UnitElementTest, PushMiddlePointRejectsInvalidLocalNodeIndicesWithoutMutation)
{
    ONEFLOW::UnitElement element;
    element.Initialize(ONEFLOW::BAR_3);

    EXPECT_THROW(element.PushMiddlePoint(0, 1, 2), std::out_of_range);
    EXPECT_THROW(element.PushMiddlePoint(4, 1, 2), std::out_of_range);
    EXPECT_THROW(element.PushMiddlePoint(3, 0, 2), std::out_of_range);
    EXPECT_THROW(element.PushMiddlePoint(3, 1, 4), std::out_of_range);

    EXPECT_EQ(element.GetRelatedPointListForMiddlePointCalculation(2).size(), 2);
    EXPECT_EQ(element.GetRelatedPointListForMiddlePointCalculation(0).size(), 0);
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


TEST(UnitElementTest, TopologyBuildersRejectInvalidLocalNodeIndicesWithoutMutation)
{
    ONEFLOW::UnitElement element;
    element.Initialize(ONEFLOW::QUAD_4);

    const int initialFaceCount = element.GetElementFaceNumber();
    const int initialChildCount = element.GetChildElementNumbers();
    const int initialPhysicsFaceCount = element.GetElementPhysicsFaceNumber();

    EXPECT_THROW(element.PushElementFace(ONEFLOW::BAR_2, 0, 2), std::out_of_range);
    EXPECT_THROW(element.PushChildElement(ONEFLOW::TRI_3, 1, 0, 3), std::out_of_range);
    EXPECT_THROW(element.PushCompositeFace(ONEFLOW::BAR_3, 1, 2, 5), std::out_of_range);

    EXPECT_EQ(element.GetElementFaceNumber(), initialFaceCount);
    EXPECT_EQ(element.GetChildElementNumbers(), initialChildCount);
    EXPECT_EQ(element.GetElementPhysicsFaceNumber(), initialPhysicsFaceCount);
}
