#include <gtest/gtest.h>

#include "FaceSearch.h"
#include "Constant.h"

TEST(FaceSearchBasicTest, AddFaceAssignsStableIdsAndDeduplicates)
{
    // Intent: faces with the same node set (order-independent via lookup) share
    // one ID; distinct node sets get distinct IDs. Size tracks unique faces.
    ONEFLOW::FaceSearchBasic search;

    const int id0 = search.AddFace(ONEFLOW::IntField{0, 1, 2});
    const int id1 = search.AddFace(ONEFLOW::IntField{3, 4, 5});
    EXPECT_EQ(id0, 0);
    EXPECT_EQ(id1, 1);
    EXPECT_EQ(search.Size(), static_cast<std::size_t>(2));

    // Exact same node sequence is treated as the same face.
    const int id0Again = search.AddFace(ONEFLOW::IntField{0, 1, 2});
    EXPECT_EQ(id0Again, id0);
    EXPECT_EQ(search.Size(), static_cast<std::size_t>(2));
}

TEST(FaceSearchBasicTest, FindFaceReturnsInvalidWhenMissing)
{
    ONEFLOW::FaceSearchBasic search;
    search.AddFace(ONEFLOW::IntField{1, 2});

    EXPECT_EQ(search.FindFace(ONEFLOW::IntField{1, 2}), 0);
    EXPECT_EQ(search.FindFace(ONEFLOW::IntField{9, 8}), ONEFLOW::INVALID_INDEX);
}

TEST(FaceSearchBasicTest, ClearRemovesAllFaces)
{
    ONEFLOW::FaceSearchBasic search;
    search.AddFace(ONEFLOW::IntField{0, 1});
    search.AddFace(ONEFLOW::IntField{2, 3});
    EXPECT_EQ(search.Size(), static_cast<std::size_t>(2));

    search.Clear();
    EXPECT_EQ(search.Size(), static_cast<std::size_t>(0));
    EXPECT_EQ(search.FindFace(ONEFLOW::IntField{0, 1}), ONEFLOW::INVALID_INDEX);
}

TEST(FaceSearchBasicTest, FaceArrayPreservesOriginalNodeOrder)
{
    // Intent: faceArray stores the original node order of each added face so
    // downstream topology can recover orientation, not only the lookup key.
    ONEFLOW::FaceSearchBasic search;
    search.AddFace(ONEFLOW::IntField{5, 6, 7, 8});

    ASSERT_EQ(search.faceArray.size(), static_cast<std::size_t>(1));
    ASSERT_EQ(search.faceArray[0].size(), static_cast<std::size_t>(4));
    EXPECT_EQ(search.faceArray[0][0], 5);
    EXPECT_EQ(search.faceArray[0][1], 6);
    EXPECT_EQ(search.faceArray[0][2], 7);
    EXPECT_EQ(search.faceArray[0][3], 8);
}

TEST(FaceSearchBasicTest, InterfaceFaceSplitTable_ResetForParentFaceCount)
{
    // Intent: ResetForParentFaceCount prepares one empty child-id list per parent
    // face and parallel node/flag tables of the same parent capacity.
    ONEFLOW::InterfaceFaceSplitTable table;
    table.childFaceIdsOfParent.push_back(ONEFLOW::IntField{1});
    table.ResetForParentFaceCount(3);

    ASSERT_EQ(table.childFaceIdsOfParent.size(), static_cast<std::size_t>(3));
    ASSERT_EQ(table.childFaceNodeIds.size(), static_cast<std::size_t>(3));
    ASSERT_EQ(table.childFaceNodeFlags.size(), static_cast<std::size_t>(3));
    EXPECT_TRUE(table.childFaceIdsOfParent[0].empty());
    EXPECT_TRUE(table.childFaceIdsOfParent[1].empty());
    EXPECT_TRUE(table.childFaceIdsOfParent[2].empty());
}

TEST(FaceSearchBasicTest, InterfaceFaceSplitTable_ClearEmptiesAllTables)
{
    ONEFLOW::InterfaceFaceSplitTable table;
    table.ResetForParentFaceCount(2);
    table.childFaceIdsOfParent[0].push_back(9);
    table.childFaceNodeIds[0].push_back(1);
    table.childFaceNodeFlags[0].push_back(1);

    table.Clear();
    EXPECT_TRUE(table.childFaceIdsOfParent.empty());
    EXPECT_TRUE(table.childFaceNodeIds.empty());
    EXPECT_TRUE(table.childFaceNodeFlags.empty());
}

TEST(FaceSearchBasicTest, FaceSearchKeepsSplitTableSeparateFromRegistry)
{
    // Intent: FaceSearch owns both a face registry (base) and an interfaceSplit
    // table; clearing the registry Size does not require split table to be empty
    // and vice versa -- they are independent responsibilities.
    ONEFLOW::FaceSearch search;
    search.AddFace(ONEFLOW::IntField{0, 1, 2});
    search.interfaceSplit.ResetForParentFaceCount(1);
    search.interfaceSplit.childFaceIdsOfParent[0].push_back(0);

    EXPECT_EQ(search.Size(), static_cast<std::size_t>(1));
    EXPECT_EQ(search.interfaceSplit.childFaceIdsOfParent[0].size(), static_cast<std::size_t>(1));

    search.Clear();
    EXPECT_EQ(search.Size(), static_cast<std::size_t>(0));
    // Split table is independent; Clear() on registry does not touch it.
    EXPECT_EQ(search.interfaceSplit.childFaceIdsOfParent[0].size(), static_cast<std::size_t>(1));
}
