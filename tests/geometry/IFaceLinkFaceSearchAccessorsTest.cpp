// tests/geometry/IFaceLinkFaceSearchAccessorsTest.cpp
// Locks private face_search exposure through accessors.

#include <gtest/gtest.h>

#include "IFaceLink.h"
#include "FaceSearch.h"
#include "GridHandles.h"
#include "HXDefine.h"

using namespace ONEFLOW;

TEST( IFaceLinkFaceSearchAccessorsTest, SplitTableAliasesFaceSearch )
{
    Grids grids;
    IFaceLink link( grids );

    EXPECT_EQ( &link.MutableInterfaceSplit(), &link.GetFaceSearch().interfaceSplit );
    EXPECT_EQ( &link.GetInterfaceSplit(), &link.GetFaceSearch().interfaceSplit );
}

TEST( IFaceLinkFaceSearchAccessorsTest, MutableSplitAcceptsParentResize )
{
    Grids grids;
    IFaceLink link( grids );
    link.MutableInterfaceSplit().ResetForFaceCount( 2 );
    EXPECT_EQ( link.GetInterfaceSplit().childFaceIdsOfParent.size(),
               static_cast< std::size_t >( 2 ) );
}
