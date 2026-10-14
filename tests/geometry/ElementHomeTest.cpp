#include <gtest/gtest.h>
#include <cstddef>
#include <stdexcept>
#include <type_traits>
#include <utility>

#include "ElementHome.h"
#include "HXCgns.h"

static_assert(
    std::is_const_v<std::remove_reference_t<
        decltype(ONEFLOW::ElementHome::GetUnitElementDefinition(ONEFLOW::HEXA_8))>>,
    "ElementHome must expose unit element definitions as read-only");

using MutableElementFaceResult = decltype(
    std::declval< ONEFLOW::UnitElement & >().GetElementFace(0));
using MutableCompositeFaceResult = decltype(
    std::declval< ONEFLOW::UnitElement & >().GetCompositeFace(0));
using MutableMiddlePointResult = decltype(
    std::declval< ONEFLOW::UnitElement & >()
        .GetRelatedPointListForMiddlePointCalculation(0));
using MutableChildNodeResult = decltype(
    std::declval< ONEFLOW::UnitElement & >()
        .GetChildElementRelativeNodeIndex(0));

using ConstElementFaceResult = decltype(
    std::declval< const ONEFLOW::UnitElement & >().GetElementFace(0));
using ConstCompositeFaceResult = decltype(
    std::declval< const ONEFLOW::UnitElement & >().GetCompositeFace(0));
using ConstMiddlePointResult = decltype(
    std::declval< const ONEFLOW::UnitElement & >()
        .GetRelatedPointListForMiddlePointCalculation(0));
using ConstChildNodeResult = decltype(
    std::declval< const ONEFLOW::UnitElement & >()
        .GetChildElementRelativeNodeIndex(0));

static_assert(
    std::is_const_v<std::remove_reference_t<MutableElementFaceResult>>,
    "A mutable unit element must not expose mutable face topology");
static_assert(
    std::is_const_v<std::remove_reference_t<MutableCompositeFaceResult>>,
    "A mutable unit element must not expose mutable composite-face topology");
static_assert(
    std::is_const_v<std::remove_reference_t<MutableMiddlePointResult>>,
    "A mutable unit element must not expose mutable midpoint relations");
static_assert(
    std::is_const_v<std::remove_reference_t<MutableChildNodeResult>>,
    "A mutable unit element must not expose mutable child mappings");

static_assert(
    std::is_const_v<std::remove_reference_t<ConstElementFaceResult>>,
    "Reading a face from a const unit element must not expose mutable topology");
static_assert(
    std::is_const_v<std::remove_reference_t<ConstCompositeFaceResult>>,
    "Reading a composite face from a const unit element must not expose mutable topology");
static_assert(
    std::is_const_v<std::remove_reference_t<ConstMiddlePointResult>>,
    "Reading middle-point relations from a const unit element must be read-only");
static_assert(
    std::is_const_v<std::remove_reference_t<ConstChildNodeResult>>,
    "Reading child-element node mappings from a const unit element must be read-only");

TEST(ElementHomeTest, InitializeBuildsRegistryWithStableEnumIndices)
{
    ONEFLOW::ElementHome::Initialize();

    EXPECT_GT(ONEFLOW::ElementHome::GetUnitElementSlotCount(), ONEFLOW::HEXA_8);
    EXPECT_EQ(
        ONEFLOW::ElementHome::GetNumberOfUnitElements(),
        ONEFLOW::ElementHome::GetUnitElementSlotCount());
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
    EXPECT_EQ(element.GetElementFaceCount(), 6);
    EXPECT_EQ(element.GetElementFaceNumber(), element.GetElementFaceCount());
}

TEST(UnitElementTest, StandardElementDefinitionsHaveExpectedTopologyCounts)
{
    struct ExpectedTopology
    {
        int elementType;
        int faceCount;
        int compositeFaceCount;
        int childElementCount;
    };

    const ExpectedTopology expected[] = {
        {ONEFLOW::NODE, 1, 0, 0},
        {ONEFLOW::BAR_2, 2, 0, 0},
        {ONEFLOW::BAR_3, 2, 0, 2},
        {ONEFLOW::TRI_3, 3, 0, 0},
        {ONEFLOW::TRI_6, 6, 3, 4},
        {ONEFLOW::QUAD_4, 4, 0, 0},
        {ONEFLOW::QUAD_8, 8, 4, 0},
        {ONEFLOW::TETRA_4, 4, 0, 0},
        {ONEFLOW::HEXA_8, 6, 0, 0},
    };

    ONEFLOW::ElementHome::Initialize();

    for (const ExpectedTopology & item : expected)
    {
        SCOPED_TRACE(item.elementType);
        ASSERT_TRUE(ONEFLOW::ElementHome::HasUnitElement(item.elementType));

        const ONEFLOW::UnitElement & element =
            ONEFLOW::ElementHome::GetUnitElementDefinition(item.elementType);

        EXPECT_EQ(element.GetElementFaceCount(), item.faceCount);
        EXPECT_EQ(element.GetCompositeFaceCount(), item.compositeFaceCount);
        EXPECT_EQ(element.GetChildElementCount(), item.childElementCount);
    }
}

TEST(ElementHomeTest, RegisteredTopologyIndicesMatchElementDefinitions)
{
    ONEFLOW::ElementHome::Initialize();

    for (int elementType = 0;
         elementType < ONEFLOW::ElementHome::GetUnitElementSlotCount();
         ++elementType)
    {
        if (!ONEFLOW::ElementHome::HasUnitElement(elementType))
        {
            continue;
        }

        const ONEFLOW::UnitElement & element =
            ONEFLOW::ElementHome::GetUnitElementDefinition(elementType);
        const int nodeCount = ONEFLOW::GetElementNodeNumbers(elementType);

        ASSERT_EQ(element.GetElementType(), elementType);
        ASSERT_GT(nodeCount, 0);

        for (int iFace = 0; iFace < element.GetElementFaceCount(); ++iFace)
        {
            const auto & face = element.GetElementFace(iFace);
            const int faceType = element.GetFaceType(iFace);

            ASSERT_EQ(
                face.size(),
                static_cast<std::size_t>(
                    ONEFLOW::GetElementNodeNumbers(faceType)))
                << "element type " << elementType << ", face " << iFace;

            for (int nodeIndex : face)
            {
                EXPECT_GE(nodeIndex, 0)
                    << "element type " << elementType << ", face " << iFace;
                EXPECT_LT(nodeIndex, nodeCount)
                    << "element type " << elementType << ", face " << iFace;
            }
        }

        for (int iFace = 0; iFace < element.GetCompositeFaceCount(); ++iFace)
        {
            const auto & face = element.GetCompositeFace(iFace);
            const int faceType = element.GetCompositeFaceType(iFace);

            ASSERT_EQ(
                face.size(),
                static_cast<std::size_t>(
                    ONEFLOW::GetElementNodeNumbers(faceType)))
                << "element type " << elementType
                << ", composite face " << iFace;

            for (int nodeIndex : face)
            {
                EXPECT_GE(nodeIndex, 0)
                    << "element type " << elementType
                    << ", composite face " << iFace;
                EXPECT_LT(nodeIndex, nodeCount)
                    << "element type " << elementType
                    << ", composite face " << iFace;
            }
        }

        for (int iChild = 0; iChild < element.GetChildElementCount(); ++iChild)
        {
            const auto & childNodes =
                element.GetChildElementRelativeNodeIndex(iChild);
            const int childType = element.GetChildElementType(iChild);

            ASSERT_EQ(
                childNodes.size(),
                static_cast<std::size_t>(
                    ONEFLOW::GetElementNodeNumbers(childType)))
                << "element type " << elementType << ", child " << iChild;

            for (int nodeIndex : childNodes)
            {
                EXPECT_GE(nodeIndex, 0)
                    << "element type " << elementType << ", child " << iChild;
                EXPECT_LT(nodeIndex, nodeCount)
                    << "element type " << elementType << ", child " << iChild;
            }
        }

        for (int iNode = 0; iNode < nodeCount; ++iNode)
        {
            const auto & relatedNodes =
                element.GetRelatedPointListForMiddlePointCalculation(iNode);

            for (int nodeIndex : relatedNodes)
            {
                EXPECT_GE(nodeIndex, 0)
                    << "element type " << elementType
                    << ", midpoint node " << iNode;
                EXPECT_LT(nodeIndex, nodeCount)
                    << "element type " << elementType
                    << ", midpoint node " << iNode;
            }
        }
    }
}

TEST(ElementHomeTest, RepeatedInitializeRebuildsEquivalentRegistry)
{
    ONEFLOW::ElementHome::Initialize();
    const int initialSlotCount = ONEFLOW::ElementHome::GetUnitElementSlotCount();
    ASSERT_TRUE(ONEFLOW::ElementHome::HasUnitElement(ONEFLOW::HEXA_8));
    const int initialNodeCount =
        ONEFLOW::GetElementNodeNumbers(ONEFLOW::HEXA_8);

    ONEFLOW::ElementHome::Initialize();

    EXPECT_EQ(ONEFLOW::ElementHome::GetUnitElementSlotCount(), initialSlotCount);
    ASSERT_TRUE(ONEFLOW::ElementHome::HasUnitElement(ONEFLOW::HEXA_8));
    EXPECT_EQ(ONEFLOW::GetElementNodeNumbers(ONEFLOW::HEXA_8), initialNodeCount);
}

TEST(ElementHomeTest, InvalidElementIndicesAreRejectedConsistently)
{
    EXPECT_FALSE(ONEFLOW::ElementHome::HasUnitElement(-1));
    EXPECT_FALSE(
        ONEFLOW::ElementHome::HasUnitElement(
            ONEFLOW::ElementHome::GetUnitElementSlotCount()));

    EXPECT_THROW(
        ONEFLOW::ElementHome::GetUnitElementDefinition(-1),
        std::out_of_range);
    EXPECT_THROW(
        ONEFLOW::GetElementNodeNumbers(
            ONEFLOW::ElementHome::GetUnitElementSlotCount()),
        std::out_of_range);
}

TEST(ElementHomeTest, FreeClearsRegistryAndInitializeRestoresIt)
{
    ONEFLOW::ElementHome::Initialize();
    ASSERT_TRUE(ONEFLOW::ElementHome::HasUnitElement(ONEFLOW::HEXA_8));

    ONEFLOW::ElementHome::Free();

    EXPECT_EQ(ONEFLOW::ElementHome::GetUnitElementSlotCount(), 0);
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
    EXPECT_EQ(element.GetElementFaceCount(), 0);
    EXPECT_EQ(element.GetCompositeFaceCount(), 0);
    EXPECT_EQ(element.GetChildElementCount(), 0);
}

TEST(UnitElementTest, ReinitializeReplacesThePreviousTopologyDefinition)
{
    ONEFLOW::UnitElement element;
    element.Initialize(ONEFLOW::QUAD_4);

    element.Initialize(ONEFLOW::BAR_2);

    EXPECT_EQ(element.GetElementType(), ONEFLOW::BAR_2);
    EXPECT_EQ(element.GetElementFaceCount(), 2);
    EXPECT_EQ(element.GetChildElementCount(), 0);
    EXPECT_EQ(element.GetCompositeFaceCount(), 0);
    EXPECT_EQ(
        element.GetRelatedPointListForMiddlePointCalculation(1).size(), 0);
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

    EXPECT_THROW(element.GetCompositeFace(-1), std::out_of_range);
    EXPECT_THROW(element.GetCompositeFace(0), std::out_of_range);
    EXPECT_THROW(element.GetCompositeFaceType(-1), std::out_of_range);
    EXPECT_THROW(element.GetCompositeFaceType(0), std::out_of_range);

    EXPECT_THROW(element.GetChildElementRelativeNodeIndex(-1),
                 std::out_of_range);
    EXPECT_THROW(element.GetChildElementRelativeNodeIndex(0),
                 std::out_of_range);
    EXPECT_THROW(element.GetChildElementType(-1), std::out_of_range);
    EXPECT_THROW(element.GetChildElementType(0), std::out_of_range);
}

TEST(UnitElementTest, FailedInitializePreservesExistingTopology)
{
    ONEFLOW::UnitElement element;
    element.Initialize(ONEFLOW::QUAD_4);

    const int initialElementType = element.GetElementType();
    const int initialFaceCount = element.GetElementFaceCount();
    const int initialNodeCount = element.GetElementNodeNumbers(ONEFLOW::QUAD_4);

    EXPECT_THROW(element.Initialize(-1), std::invalid_argument);

    EXPECT_EQ(element.GetElementType(), initialElementType);
    EXPECT_EQ(element.GetElementFaceCount(), initialFaceCount);
    EXPECT_EQ(element.GetElementNodeNumbers(ONEFLOW::QUAD_4), initialNodeCount);
    ASSERT_EQ(element.GetElementFaceCount(), 4);
    EXPECT_EQ(element.GetElementFace(0).size(), 2);
}

