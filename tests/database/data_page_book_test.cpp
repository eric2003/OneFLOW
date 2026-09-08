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
    EXPECT_EQ(page.GetSize(), 0U);
    EXPECT_EQ(page.GetBeginDataPointer(), nullptr);
    EXPECT_EQ(page.GetCurrentDataPointer(), nullptr);
}

TEST_F(DataPageBookTest, DataPage_Resize_Basic)
{
    DataPage page;
    page.ReSize(1024);
    EXPECT_EQ(page.GetSize(), 1024U);
    EXPECT_NE(page.GetBeginDataPointer(), nullptr);
    EXPECT_NE(page.GetCurrentDataPointer(), nullptr);
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
    page.Write(&srcVal, sizeof(int), 100);

    int dstVal = 0;
    page.Read(&dstVal, sizeof(int), 100);
    EXPECT_EQ(srcVal, dstVal);
}

TEST_F(DataPageBookTest, DataPage_MoveToEnd_PositionEqualSize)
{
    // Important boundary: position == GetSize() should be valid for append
    DataPage page;
    page.ReSize(128);
    page.MoveToEnd();
    // currPos equals size, should not trigger fatal error
    EXPECT_EQ(page.GetCurrentDataPointer(), page.GetDataPointer(128));
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

// Note: Out-of-bound memcpy / Fatal crash test cannot run in gtest,
// because Fatal() terminates whole process. Need manual test.

// ----------------------------------------------------------------------------
// DataBook Basic Test
// ----------------------------------------------------------------------------
TEST_F(DataPageBookTest, DataBook_Create_Empty)
{
    DataBook book;
    EXPECT_EQ(book.GetSize(), 0LL);
}

TEST_F(DataPageBookTest, DataBook_Write_Read_Small_NoCrossPage)
{
    DataBook book;
    const HXLongLong_t bufSize = 1024;
    book.SecureAbsoluteSpace(bufSize);

    char src[bufSize];
    for (HXLongLong_t i = 0; i < bufSize; ++i)
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

TEST_F(DataPageBookTest, DataBook_Write_Read_CrossPageBoundary)
{
    DataBook book;
    // Force cross‑page: override maxUnitSize for test, small page size
    // Note: original member maxUnitSize is not public,
    // If you add setter void SetMaxUnitSize(HXLongLong_t s), enable below.
    // book.maxUnitSize = 100;

    // Test data cross page boundary, e.g total size 250, page size=100 -> 3 pages
    const HXLongLong_t totalSize = 250;
    book.SecureAbsoluteSpace(totalSize);

    char src[250];
    for (HXLongLong_t i = 0; i < totalSize; ++i)
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

TEST_F(DataPageBookTest, DataBook_Append)
{
    DataBook book;
    char block1[] = {0x11, 0x22, 0x33};
    char block2[] = {0x44, 0x55};

    book.Append(block1, 3);
    book.Append(block2, 2);

    EXPECT_EQ(book.GetSize(), 5LL);

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
    // GetSize(); if we pre-size the book, "end" is no longer position 0, so
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
    EXPECT_GT(book.GetSize(), 0LL);

    book.ReSize(0);
    EXPECT_EQ(book.GetSize(), 0LL);
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

