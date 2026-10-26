// tests/geometry/ExpandBcRecordSplitTest.cpp
// Locks BC array expansion when interface faces split into children.

#include <gtest/gtest.h>
#include <stdexcept>

#include "FaceTopo.h"
#include "Boundary.h"
#include "HXDefine.h"

using namespace ONEFLOW;

TEST( ExpandBcRecordSplitTest, NonInterfaceFacesEmitOnce )
{
    const IntField bcType{ BC::SOLID_SURFACE, BC::INFLOW };
    const IntField bcNameId{ 10, 20 };
    const IntField childCounts; // no interface faces
    IntField outType, outName;

    ExpandBcRecordForSplitInterfaces(
        bcType, bcNameId, childCounts, outType, outName );

    EXPECT_EQ( outType, bcType );
    EXPECT_EQ( outName, bcNameId );
}

TEST( ExpandBcRecordSplitTest, UnsplitInterfaceEmitsOnce )
{
    const IntField bcType{ BC::SOLID_SURFACE, BC::INTERFACE, BC::FARFIELD };
    const IntField bcNameId{ 1, 2, 3 };
    const IntField childCounts{ 0 }; // interface local 0 unsplit
    IntField outType, outName;

    ExpandBcRecordForSplitInterfaces(
        bcType, bcNameId, childCounts, outType, outName );

    EXPECT_EQ( outType, ( IntField{ BC::SOLID_SURFACE, BC::INTERFACE, BC::FARFIELD } ) );
    EXPECT_EQ( outName, ( IntField{ 1, 2, 3 } ) );
}

TEST( ExpandBcRecordSplitTest, SplitInterfaceEmitsPerChild )
{
    const IntField bcType{ BC::INTERFACE, BC::SOLID_SURFACE };
    const IntField bcNameId{ 7, 8 };
    const IntField childCounts{ 3 };
    IntField outType, outName;

    ExpandBcRecordForSplitInterfaces(
        bcType, bcNameId, childCounts, outType, outName );

    EXPECT_EQ( outType,
               ( IntField{ BC::INTERFACE, BC::INTERFACE, BC::INTERFACE, BC::SOLID_SURFACE } ) );
    EXPECT_EQ( outName, ( IntField{ 7, 7, 7, 8 } ) );
}

TEST( ExpandBcRecordSplitTest, NegativeInterfaceCodeCountsAsInterface )
{
    const IntField bcType{ -5, BC::INFLOW };
    const IntField bcNameId{ 1, 2 };
    const IntField childCounts{ 2 };
    IntField outType, outName;

    ExpandBcRecordForSplitInterfaces(
        bcType, bcNameId, childCounts, outType, outName );

    EXPECT_EQ( outType, ( IntField{ -5, -5, BC::INFLOW } ) );
    EXPECT_EQ( outName, ( IntField{ 1, 1, 2 } ) );
}

TEST( ExpandBcRecordSplitTest, RejectsChildCountShortfall )
{
    const IntField bcType{ BC::INTERFACE };
    const IntField bcNameId{ 1 };
    const IntField childCounts; // empty but one interface face
    IntField outType, outName;

    EXPECT_THROW(
        ExpandBcRecordForSplitInterfaces(
            bcType, bcNameId, childCounts, outType, outName ),
        std::out_of_range );
}
