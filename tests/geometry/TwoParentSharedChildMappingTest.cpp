// tests/geometry/TwoParentSharedChildMappingTest.cpp
// Locks parent-specific node mappings when two parents share one global child.
// See face_stage_pipeline.rst Open design pressures / parent-relative mapping.

#include <gtest/gtest.h>
#include <stdexcept>

#include "FaceSearch.h"
#include "FaceTopo.h"
#include "HXDefine.h"

using namespace ONEFLOW;

namespace {

// Parent A: [10,20,30,40]  child triangle uses local corners 1,2,3 -> [20,30,40]
// Parent B: [20,30,40,50]  same undirected child ring maps as local 0,1,2 -> [20,30,40]
// Global parent ids 0 and 1; global child id 2.

void SeedTwoParentsOneChild( InterfaceFaceSplitTable & split )
{
    split.ResetForFaceCount( 3 ); // parents 0,1 and child slot 2

    // Parent 0 owns child 2 once
    split.childFaceIdsOfParent[ 0 ] = IntField{ 2 };
    // Parent 1 owns the same global child 2 once
    split.childFaceIdsOfParent[ 1 ] = IntField{ 2 };

    // Different parent-local encodings for the same global child identity
    split.SetChildMapping(
        /*parentFaceId=*/0,
        /*childFaceId=*/2,
        IntField{ 1, 2, 3 },
        IntField{ 1, 1, 1 } ); // relative to parent A ring

    split.SetChildMapping(
        /*parentFaceId=*/1,
        /*childFaceId=*/2,
        IntField{ 0, 1, 2 },
        IntField{ 1, 1, 1 } ); // relative to parent B ring
}

} // namespace

TEST( TwoParentSharedChildMappingTest, ParentSpecificRowsStayIndependent )
{
    InterfaceFaceSplitTable split;
    SeedTwoParentsOneChild( split );

    ASSERT_EQ( split.childFaceNodeIdsByParent[ 0 ].size(), static_cast< std::size_t >( 1 ) );
    ASSERT_EQ( split.childFaceNodeIdsByParent[ 1 ].size(), static_cast< std::size_t >( 1 ) );

    EXPECT_EQ( split.childFaceNodeIdsByParent[ 0 ][ 0 ], ( IntField{ 1, 2, 3 } ) );
    EXPECT_EQ( split.childFaceNodeIdsByParent[ 1 ][ 0 ], ( IntField{ 0, 1, 2 } ) );
    EXPECT_EQ( split.childFaceNodeFlagsByParent[ 0 ][ 0 ], ( IntField{ 1, 1, 1 } ) );
    EXPECT_EQ( split.childFaceNodeFlagsByParent[ 1 ][ 0 ], ( IntField{ 1, 1, 1 } ) );
}

TEST( TwoParentSharedChildMappingTest, ResolveUsesEachParentRing )
{
    InterfaceFaceSplitTable split;
    SeedTwoParentsOneChild( split );

    const IntField parentA{ 10, 20, 30, 40 };
    const IntField parentB{ 20, 30, 40, 50 };

    const IntField resolvedA = ResolveChildFaceNodeIds(
        parentA,
        split.childFaceNodeIdsByParent[ 0 ][ 0 ],
        split.childFaceNodeFlagsByParent[ 0 ][ 0 ],
        /*hasNodeMesh=*/false,
        /*nodeMeshCount=*/0 );
    const IntField resolvedB = ResolveChildFaceNodeIds(
        parentB,
        split.childFaceNodeIdsByParent[ 1 ][ 0 ],
        split.childFaceNodeFlagsByParent[ 1 ][ 0 ],
        false,
        0 );

    EXPECT_EQ( resolvedA, ( IntField{ 20, 30, 40 } ) );
    EXPECT_EQ( resolvedB, ( IntField{ 20, 30, 40 } ) );
}

TEST( TwoParentSharedChildMappingTest, GlobalSeedFilledOnceDoesNotReplaceParentRows )
{
    InterfaceFaceSplitTable split;
    SeedTwoParentsOneChild( split );

    // First SetChildMapping seeds global childFaceNodeIds[2] from parent 0.
    // Parent 1 still has its own ByParent row; global seed is not the authority
    // for rebuild when ByParent rows exist.
    EXPECT_EQ( split.childFaceNodeIds[ 2 ], ( IntField{ 1, 2, 3 } ) );
    EXPECT_EQ( split.childFaceNodeIdsByParent[ 1 ][ 0 ], ( IntField{ 0, 1, 2 } ) );
}
