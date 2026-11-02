// test_data_page_book.cpp
// Unit test for DataPage and DataBook, using Google Test
#include <gtest/gtest.h>
#include <stdexcept>
#include <string>
#include <fstream>
#include <cstring>

// OneFLOW headers
#include "DataPage.h"
#include "DataBook.h"
#include "HXType.h"

#ifdef ONEFLOW_TEST_MPI
#include "Parallel.h"
#endif

using namespace ONEFLOW;

// ============================================================================
// Test Fixture for DataPage & DataBook
// ============================================================================
class DataPageBookTest : public ::testing::Test
{
protected:
    void SetUp() override
    {
        // Executed before each TEST_F
    }

    void TearDown() override
    {
        // Executed after each TEST_F
    }

    // Helper: create temp file name for file IO test
    std::string GetTempFileName()
    {
        return "tmp_datapage_book_test.bin";
    }
};

// ----------------------------------------------------------------------------
// DataPage Basic Function Test
// ----------------------------------------------------------------------------
TEST_F(DataPageBookTest, DataPage_Create_Destroy_Empty)
{
    DataPage page;
    EXPECT_EQ(page.size(), 0U);
    EXPECT_EQ(page.data(), nullptr);
    EXPECT_EQ(page.CurrentPtr(), nullptr);
}

TEST_F(DataPageBookTest, DataPage_Resize_Basic)
{
    DataPage page;
    page.ReSize(1024);
    EXPECT_EQ(page.size(), 1024U);
    EXPECT_NE(page.data(), nullptr);
    EXPECT_NE(page.CurrentPtr(), nullptr);
}

TEST_F(DataPageBookTest, DataPage_Write_Read_SmallBuffer)
{
    DataPage page;
    const HXSize_t bufSize = 64;
    page.ReSize(bufSize);

    char src[bufSize];
    for (HXSize_t i = 0; i < bufSize; ++i)
    {
        src[i] = static_cast<char>(i & 0xff);
    }

    page.MoveToBegin();
    page.Write(src, bufSize);

    page.MoveToBegin();
    char dst[bufSize];
    page.Read(dst, bufSize);

    EXPECT_EQ(std::memcmp(src, dst, bufSize), 0);
}

TEST_F(DataPageBookTest, DataPage_Write_Read_RandomPosition)
{
    DataPage page;
    page.ReSize(256);

    int srcVal = 0x12345678;
    page.Write(&srcVal, 100, sizeof(int));

    int dstVal = 0;
    page.Read(&dstVal, 100, sizeof(int));
    EXPECT_EQ(srcVal, dstVal);
}

TEST_F(DataPageBookTest, DataPage_MoveToEnd_PositionEqualSize)
{
    // Important boundary: position == size() should be valid for append
    DataPage page;
    page.ReSize(128);
    page.MoveToEnd();
    // currPos equals size, should not trigger fatal error
    EXPECT_EQ( page.CurrentPtr(), page.PtrAt( 128 ) );
}

TEST_F(DataPageBookTest, DataPage_ToString)
{
    DataPage page;
    std::string srcStr = "HelloOneFlowDataPage";
    page.ReSize(srcStr.size());
    page.MoveToBegin();
    page.Write(const_cast<char*>(srcStr.c_str()), srcStr.size());

    std::string outStr;
    page.ToString(outStr);
    EXPECT_EQ(outStr, srcStr);
}

TEST_F(DataPageBookTest, DataPage_FileIO)
{
    DataPage pageWrite;
    const std::string tmpFile = GetTempFileName();
    const HXSize_t testSize = 128;
    pageWrite.ReSize(testSize);

    char src[testSize];
    for (HXSize_t i = 0; i < testSize; ++i)
    {
        src[i] = static_cast<char>(i);
    }
    pageWrite.MoveToBegin();
    pageWrite.Write(src, testSize);

    {
        std::fstream fout(tmpFile, std::ios::binary | std::ios::out);
        ASSERT_TRUE(fout.is_open());
        pageWrite.WriteFile(fout);
        fout.close();
    }

    DataPage pageRead;
    pageRead.ReSize(testSize);
    {
        std::fstream fin(tmpFile, std::ios::binary | std::ios::in);
        ASSERT_TRUE(fin.is_open());
        pageRead.ReadFile(fin);
        fin.close();
    }

    char dst[testSize];
    pageRead.MoveToBegin();
    pageRead.Read(dst, testSize);
    EXPECT_EQ(std::memcmp(src, dst, testSize), 0);

    // clean temp file
    std::remove(tmpFile.c_str());
}

// Fatal() now throws std::runtime_error instead of terminating the
// process, so the out-of-range path can be exercised directly in gtest.
TEST_F(DataPageBookTest, DataPage_Write_OutOfRangePosition_Throws)
{
    DataPage page;
    page.ReSize(10);

    char buf[4] = {0};
    // Update expected exception type to std::out_of_range.
    EXPECT_THROW(page.Write(buf, 11, sizeof(buf)), std::out_of_range);
}

TEST_F(DataPageBookTest, DataPage_Read_OutOfRangePosition_Throws)
{
    DataPage page;
    page.ReSize(10);

    char buf[4] = {0};
    // Update expected exception type to std::out_of_range.
    EXPECT_THROW(page.Read(buf, 11, sizeof(buf)), std::out_of_range);
}
// ----------------------------------------------------------------------------
// DataBook Basic Test
// ----------------------------------------------------------------------------
TEST_F(DataPageBookTest, DataBook_Create_Empty)
{
    DataBook book;
    EXPECT_EQ(book.size(), 0LL);
}

TEST_F(DataPageBookTest, DataBook_Write_Read_Small_NoCrossPage)
{
    DataBook book;
    const HXOffset_t bufSize = 1024;
    book.Reserve(bufSize);

    char src[bufSize];
    for (HXOffset_t i = 0; i < bufSize; ++i)
    {
        src[i] = static_cast<char>(i & 0xff);
    }

    book.MoveToBegin();
    book.Write(src, bufSize);

    book.MoveToBegin();
    char dst[bufSize];
    book.Read(dst, bufSize);

    EXPECT_EQ(std::memcmp(src, dst, static_cast<std::size_t>(bufSize)), 0);
}



TEST_F(DataPageBookTest, DataBook_Write_Read_SingleLargePage)
{
    DataBook book;
    // Force cross‑page: override maxUnitSize for test, small page size
    // Note: original member maxUnitSize is not public,
    // If you add setter void SetMaxUnitSize(HXOffset_t s), enable below.
    // book.maxUnitSize = 100;

    // Test data cross page boundary, e.g total size 250, page size=100 -> 3 pages
    const HXOffset_t totalSize = 250;
    book.Reserve(totalSize);

    char src[250];
    for (HXOffset_t i = 0; i < totalSize; ++i)
    {
        src[i] = static_cast<char>(i);
    }

    book.MoveToBegin();
    book.Write(src, totalSize);

    book.MoveToBegin();
    char dst[250];
    book.Read(dst, totalSize);

    EXPECT_EQ(std::memcmp(src, dst, static_cast<std::size_t>(totalSize)), 0);
}

TEST_F(DataPageBookTest, DataBook_TrulyCrossMultiplePages)
{
    // Use a tiny page size (10 bytes) so 250 bytes of data must span
    // multiple real pages (25 pages), unlike the default 1GB page size
    // where "cross-page" tests never actually leave page 0.
    DataBook book(10);

    const HXOffset_t totalSize = 250;
    char src[250];
    for (HXOffset_t i = 0; i < totalSize; ++i)
    {
        src[i] = static_cast<char>(i);
    }

    book.MoveToBegin();
    book.Write(src, totalSize);

    // Sanity check: with a 10-byte page size and 250 bytes written,
    // the data must have actually spanned multiple pages (25 pages),
    // not silently stayed within a single page like the old test did.
    EXPECT_GT(book.GetPageCount(), 1U); 

    book.MoveToBegin();
    char dst[250];
    book.Read(dst, totalSize);

    EXPECT_EQ(std::memcmp(src, dst, totalSize), 0);
}

TEST_F(DataPageBookTest, DataBook_SmallUnitSize_StringRoundTrip)
{
    // Force a string write/read to cross several tiny pages internally.
    DataBook book(8);
    std::string testStr = "OneFLOW-CFD-DataBook-CrossPage-StringTest";

    book.WriteString(testStr);
    book.MoveToBegin();

    std::string outStr;
    book.ReadString(outStr);

    EXPECT_EQ(outStr, testStr);
}

TEST_F(DataPageBookTest, DataBook_InvalidUnitSize_Throws)
{
    EXPECT_THROW(DataBook book(0), std::invalid_argument);
    EXPECT_THROW(DataBook book(-1), std::invalid_argument);
}

TEST_F(DataPageBookTest, DataBook_Append)
{
    DataBook book;
    char block1[] = {0x11, 0x22, 0x33};
    char block2[] = {0x44, 0x55};

    book.Append(block1, 3);
    book.Append(block2, 2);

    EXPECT_EQ(book.size(), 5LL);

    char out[5];
    book.MoveToBegin();
    book.Read(out,5);
    EXPECT_EQ(out[0], 0x11);
    EXPECT_EQ(out[2], 0x33);
    EXPECT_EQ(out[3], 0x44);
    EXPECT_EQ(out[4], 0x55);
}

TEST_F(DataPageBookTest, DataBook_String_Write_Read)
{
    DataBook book;
    std::string testStr = "OneFLOW-CFD-DataBook-StringTest-0123456789";

    book.MoveToBegin();
    book.WriteString(testStr);

    book.MoveToBegin();
    std::string outStr;
    book.ReadString(outStr);

    EXPECT_EQ(outStr, testStr);
}


TEST_F(DataPageBookTest, DataBook_AppendString)
{
    DataBook book;
    std::string s1 = "AAA";
    std::string s2 = "BBB";
    book.AppendString(s1);
    book.AppendString(s2);

    book.MoveToBegin();
    std::string r1, r2;
    book.ReadString(r1);
    book.ReadString(r2);

    EXPECT_EQ(r1, "AAA");
    EXPECT_EQ(r2, "BBB");
}

TEST_F(DataPageBookTest, DataBook_MoveBegin_MoveEnd)
{
    DataBook book;
    // NOTE: do NOT pre-allocate space here. MoveToEnd() moves the cursor to
    // size(); if we pre-size the book, "end" is no longer position 0, so
    // writing at the end and reading from the begin would compare different
    // offsets. This test is only meaningful starting from an empty book.
    int val = 0xABCD;
    book.MoveToEnd();              // book is empty -> equivalent to MoveToBegin()
    book.Write(&val, sizeof(int));

    book.MoveToBegin();
    int outVal = 0;
    book.Read(&outVal, sizeof(int));
    EXPECT_EQ(val, outVal);
}

TEST_F(DataPageBookTest, DataBook_ToString)
{
    DataBook book;
    char data[] = {'X','Y','Z'};
    book.Append(data, 3);

    std::string strOut;
    book.ToString(strOut);
    EXPECT_EQ(strOut.size(), 3U);
    EXPECT_EQ(strOut[0], 'X');
    EXPECT_EQ(strOut[2], 'Z');
}

TEST_F(DataPageBookTest, DataBook_ReSize_Zero)
{
    DataBook book;
    char buf[16] = {0};
    book.Append(buf,16);
    EXPECT_GT(book.size(), 0LL);

    book.Resize(0);
    EXPECT_EQ(book.size(), 0LL);
}

TEST_F(DataPageBookTest, DataBook_FileIO_Write_Read)
{
    DataBook bookWrite;
    std::string s1("DataBookFileTest‑Part1");
    std::string s2("DataBookFileTest‑Part2");
    bookWrite.WriteString(s1);
    bookWrite.WriteString(s2);

    const std::string tmpFile = GetTempFileName();
    {
        std::fstream fout(tmpFile, std::ios::binary | std::ios::out);
        ASSERT_TRUE(fout.is_open());
        bookWrite.WriteFile(fout);
        fout.close();
    }

    DataBook bookRead;
    {
        std::fstream fin(tmpFile, std::ios::binary | std::ios::in);
        ASSERT_TRUE(fin.is_open());
        bookRead.ReadFile(fin);
        fin.close();
    }

    bookRead.MoveToBegin();
    std::string r1, r2;
    bookRead.ReadString(r1);
    bookRead.ReadString(r2);

    EXPECT_EQ(r1, s1);
    EXPECT_EQ(r2, s2);
    std::remove(tmpFile.c_str());
}

TEST_F(DataPageBookTest, DataBook_String_WithEmbeddedNull_RoundTrip)
{
    DataBook book;
    // Old implementation truncated at the first '\0' when doing cs = data;
    // New implementation must preserve embedded null bytes.
    std::string s("abc\0def", 7);   // explicit length, contains an embedded '\0'

    book.WriteString(s);
    book.MoveToBegin();

    std::string out;
    book.ReadString(out);

    EXPECT_EQ(out.size(), 7U);
    EXPECT_EQ(out, s);
}

TEST_F(DataPageBookTest, DataBook_EmptyString_RoundTrip)
{
    DataBook book;
    std::string empty;

    book.WriteString(empty);
    book.MoveToBegin();

    std::string out;
    book.ReadString(out);

    EXPECT_EQ(out, empty);
}
TEST_F(DataPageBookTest, DataBook_ManyTinyPages_NoStackOverflow)
{
    // Regression test: DataBook::Read/Write used to be recursive, causing
    // ~5000 stack frames here (unitSize=1, 5000 bytes) and taking ~0.5s.
    // After converting to an iterative implementation, this should run in
    // a few milliseconds with no risk of stack overflow for larger inputs.
    DataBook book(1);
    std::vector<char> src(5000);
    for (size_t i = 0; i < src.size(); ++i) src[i] = static_cast<char>(i);

    book.MoveToBegin();
    book.Write(src.data(), static_cast<HXOffset_t>(src.size()));

    book.MoveToBegin();
    std::vector<char> dst(5000);
    book.Read(dst.data(), static_cast<HXOffset_t>(dst.size()));

    EXPECT_EQ(src, dst);
}

TEST_F(DataPageBookTest, DataBook_SetPageCount_GrowAndShrink)
{
    // Use tiny page size so we can force multiple pages easily.
    DataBook book(10);   // maxUnitSize = 10

    // 1. Grow: write enough data to create several pages
    const HXOffset_t total = 35;   // 4 pages (10+10+10+5)
    std::vector<char> src(total, 0xAB);
    book.MoveToBegin();
    book.Write(src.data(), total);

    EXPECT_EQ(book.GetPageCount(), 4U);
    EXPECT_EQ(book.size(), total);

    // 2. Shrink via Resize (this calls SetPageCount internally)
    book.Resize(15);                 // should become 2 pages (10+5)
    EXPECT_EQ(book.GetPageCount(), 2U);
    EXPECT_EQ(book.size(), 15LL);

    // 3. Shrink to zero
    book.Resize(0);
    EXPECT_EQ(book.GetPageCount(), 0U);  // or 1U depending on your policy;
    // current implementation leaves 0 pages after Resize(0)
    EXPECT_EQ(book.size(), 0LL);

    // 4. Grow again from empty
    book.Resize(25);                 // 3 pages
    EXPECT_EQ(book.GetPageCount(), 3U);
    EXPECT_EQ(book.size(), 25LL);
}

TEST_F(DataPageBookTest, DataBook_SetPageCount_PreservesExistingDataWhenGrowing)
{
    // Optional but useful: verify that growing does not destroy already-written data.
    DataBook book(16);
    std::string s = "HelloResize";
    book.WriteString(s);

    HXSize_t oldPages = book.GetPageCount();
    book.Resize(book.size() + 100);   // force growth

    EXPECT_GT(book.GetPageCount(), oldPages);

    book.MoveToBegin();
    std::string out;
    book.ReadString(out);
    EXPECT_EQ(out, s);   // original content must still be readable
}

#include <type_traits>

TEST(DataBookDesign, MoveSemantics_RuleOfZero)
{
    // After removing the vestigial empty destructor, DataBook should get
    // move ctor/assignment for free (cheap pointer-swap of unique_ptrs),
    // while copy remains implicitly deleted because of the unique_ptr member.
    EXPECT_TRUE ( std::is_move_constructible<DataBook>::value );
    EXPECT_TRUE ( std::is_move_assignable<DataBook>::value );
    EXPECT_FALSE( std::is_copy_constructible<DataBook>::value );
    EXPECT_FALSE( std::is_copy_assignable<DataBook>::value );
}

TEST(DataBookDesign, Move_IsCheapPointerTransfer_NotDeepCopy)
{
    // Behavioral proof that moving does NOT reallocate/copy page buffers:
    // the underlying DataPage address should stay identical after the move.
    DataBook src;
    std::string s = "MoveShouldNotCopyBuffer";
    src.WriteString(s);

    char * originalPageAddr = src.GetPage(0)->data();

    DataBook dst( std::move(src) );

    EXPECT_EQ( dst.GetPage(0)->data(), originalPageAddr );

    dst.MoveToBegin();
    std::string out;
    dst.ReadString(out);
    EXPECT_EQ(out, s);
}

TEST(DataPageDesign, CopyDisabled_MoveEnabled)
{
    EXPECT_FALSE( std::is_copy_constructible<DataPage>::value );
    EXPECT_FALSE( std::is_copy_assignable<DataPage>::value );
    EXPECT_TRUE ( std::is_move_constructible<DataPage>::value );
    EXPECT_TRUE ( std::is_move_assignable<DataPage>::value );
}

TEST(DataPageDesign, Move_IsCheap_BufferAddressUnchanged)
{
    // Proves the move is a real O(1) pointer transfer, not a hidden deep
    // copy: the underlying heap buffer address must survive the move.
    DataPage src;
    src.ReSize(4096);
    src.MoveToBegin();
    char fill[4096];
    memset(fill, 0x42, sizeof(fill));
    src.Write(fill, sizeof(fill));

    char * originalAddr = src.data();

    DataPage dst( std::move(src) );

    EXPECT_EQ( dst.data(), originalAddr );
    EXPECT_EQ( dst.size(), 4096U );
}

// Test case to verify that Append() correctly handles data that crosses DataPage boundaries.
TEST( DataBookTest, Append_CrossPageBoundary_Success )
{
    // 1. Initialize DataBook with a tiny unitSize (10 bytes per page) to force multiple pages easily.
    constexpr HXOffset_t tinyUnitSize = 10;
    DataBook book( tinyUnitSize );

    // 2. Prepare string data that exceeds tinyUnitSize (25 bytes).
    std::string appendData = "1234567890123456789012345"; // 25 bytes
    ASSERT_EQ( appendData.size(), 25 );

    // 3. Perform Append. Under the old bug, this would call GetCurrentPage()->Write() directly 
    // and throw std::out_of_range or crash. With the fix, it delegates to DataBook::Write().
    EXPECT_NO_THROW( book.Append( appendData.data(), static_cast<HXOffset_t>( appendData.size() ) ) );

    // 4. Verify total size and allocated page count.
    EXPECT_EQ( book.size(), 25 );
    // 25 bytes with unitSize=10 should allocate 3 pages (10 + 10 + 5).
    EXPECT_EQ( book.GetPageCount(), 3 );

    // 5. Read back the data from the beginning to verify memory integrity across page boundaries.
    book.MoveToBegin();
    std::string readBuffer( 25, '\0' );
    book.Read( readBuffer.data(), 25 );

    EXPECT_EQ( readBuffer, appendData );
}

// ----------------------------------------------------------------------------
// MPI related test, only compiled when ONEFLOW_TEST_MPI defined
// Must launch with mpirun -np N ./test_binary
// ----------------------------------------------------------------------------
#ifdef ONEFLOW_TEST_MPI
TEST_F(DataPageBookTest, DataBook_MPI_SendRecv)
{
    int pid = Parallel::pid;
    int np = Parallel::np;
    ASSERT_GE(np,2);

    const int sendRank = 0;
    const int recvRank = 1;
    const int tag = 1001;

    if(pid == sendRank)
    {
        DataBook bookSend;
        std::string sendStr("MPI-SendRecv-TestCase-001");
        bookSend.WriteString(sendStr);
        bookSend.Send(recvRank, tag);
    }
    else if(pid == recvRank)
    {
        DataBook bookRecv;
        bookRecv.Recv(sendRank, tag);
        bookRecv.MoveToBegin();
        std::string outStr;
        bookRecv.ReadString(outStr);
        EXPECT_EQ(outStr, std::string("MPI-SendRecv-TestCase-001"));
    }
}

TEST_F(DataPageBookTest, DataBook_MPI_Bcast)
{
    int pid = Parallel::pid;
    int root = 0;
    const int tag = 1002;
    const std::string refStr("Bcast-Test-0001");

    DataBook book;
    if(pid == root)
    {
        book.WriteString(refStr);
    }

    book.Bcast(root);

    book.MoveToBegin();
    std::string out;
    book.ReadString(out);
    EXPECT_EQ(out, refStr);
}
#endif

// ----------------------------------------------------------------------------
// main for gtest
// ----------------------------------------------------------------------------
int main(int argc, char **argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}

