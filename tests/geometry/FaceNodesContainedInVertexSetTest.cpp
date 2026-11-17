// tests/geometry/FaceNodesContainedInVertexSetTest.cpp
// Locks BC face matching predicate: all face nodes in BC vertex set.

#include <gtest/gtest.h>

#include "FaceSolver.h"
#include "HXDefine.h"

using namespace ONEFLOW;

TEST( FaceNodesContainedInVertexSetTest, AllNodesPresent )
{
    IntSet bc;
    bc.insert( 1 );
    bc.insert( 2 );
    bc.insert( 3 );
    bc.insert( 9 );
    EXPECT_TRUE( FaceNodesContainedInVertexSet( bc, IntField{ 1, 2, 3 } ) );
}

TEST( FaceNodesContainedInVertexSetTest, MissingNodeFails )
{
    IntSet bc;
    bc.insert( 1 );
    bc.insert( 2 );
    EXPECT_FALSE( FaceNodesContainedInVertexSet( bc, IntField{ 1, 2, 3 } ) );
}

TEST( FaceNodesContainedInVertexSetTest, EmptyFaceIsTrue )
{
    IntSet bc;
    bc.insert( 1 );
    EXPECT_TRUE( FaceNodesContainedInVertexSet( bc, IntField{} ) );
}

TEST( FaceNodesContainedInVertexSetTest, OrderIrrelevant )
{
    IntSet bc;
    bc.insert( 3 );
    bc.insert( 1 );
    bc.insert( 2 );
    EXPECT_TRUE( FaceNodesContainedInVertexSet( bc, IntField{ 2, 1, 3 } ) );
}
