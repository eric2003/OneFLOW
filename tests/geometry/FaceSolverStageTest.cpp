// tests/geometry/FaceSolverStageTest.cpp
// Locks StageNewFace / ClaimStagedFaceAsLeft / AttachStagedFaceAsRight
// invariants for the polygon-style path without a full CGNS scan.

#include <gtest/gtest.h>
#include <stdexcept>

#include "Constant.h"
#include "FaceSolver.h"
#include "FaceTopo.h"
#include "BcRecord.h"

using namespace ONEFLOW;

TEST( FaceSolverStageTest, StageNewFaceGrowsAllParallelArrays )
{
    FaceSolver solver;
    solver.StageNewFace( IntField{ 1, 2, 3 }, /*faceType=*/5, /*leftCellId=*/0 );
    solver.StageNewFace( IntField{ 4, 5, 6, 7 }, /*faceType=*/7, INVALID_INDEX );

    ASSERT_TRUE( solver.StagedBcMatchesFaceCount() );

    const FaceTopo & topo = solver.GetFaceTopo();
    EXPECT_EQ( topo.GetFaces().size(), static_cast< std::size_t >( 2 ) );
    EXPECT_EQ( topo.GetFaceFlags().size(), static_cast< std::size_t >( 2 ) );
    EXPECT_EQ( topo.GetLeftCells().size(), static_cast< std::size_t >( 2 ) );
    EXPECT_EQ( topo.GetRightCells().size(), static_cast< std::size_t >( 2 ) );
    EXPECT_EQ( topo.lPosition.size(), static_cast< std::size_t >( 2 ) );
    EXPECT_EQ( topo.rPosition.size(), static_cast< std::size_t >( 2 ) );
    EXPECT_EQ( solver.GetStagedBcSlotCount(), static_cast< std::size_t >( 2 ) );
    EXPECT_EQ( solver.GetStagedBcSlotCount(), static_cast< std::size_t >( 2 ) );
    EXPECT_EQ( solver.GetChildFidSlotCount(), static_cast< std::size_t >( 2 ) );
    EXPECT_EQ( topo.GetFaceTypes().size(), static_cast< std::size_t >( 2 ) );

    EXPECT_EQ( topo.GetFaceFlags()[ 0 ], FaceSolver::kFaceUnclaimed );
    EXPECT_EQ( topo.GetFaceFlags()[ 1 ], FaceSolver::kFaceUnclaimed );
    EXPECT_EQ( topo.GetLeftCells()[ 1 ], INVALID_INDEX );
    EXPECT_EQ( solver.GetStagedBcType( 0 ), INVALID_INDEX );
    EXPECT_EQ( solver.GetStagedBcNameId( 1 ), INVALID_INDEX );
}

TEST( FaceSolverStageTest, ClaimStagedFaceAsLeftAssignsOwnership )
{
    FaceSolver solver;
    solver.StageNewFace( IntField{ 1, 2, 3, 4 }, /*faceType=*/7, INVALID_INDEX );
    ASSERT_TRUE( solver.StagedBcMatchesFaceCount() );

    solver.ClaimStagedFaceAsLeft( /*faceId=*/0, /*cellId=*/9 );

    const FaceTopo & topo = solver.GetFaceTopo();
    EXPECT_EQ( topo.GetFaceFlags()[ 0 ], FaceSolver::kFaceClaimed );
    EXPECT_EQ( topo.GetLeftCells()[ 0 ], 9 );
    EXPECT_EQ( topo.GetRightCells()[ 0 ], INVALID_INDEX );
    EXPECT_EQ( solver.GetStagedBcType( 0 ), INVALID_INDEX );
    EXPECT_EQ( solver.GetStagedBcNameId( 0 ), INVALID_INDEX );
    EXPECT_TRUE( solver.StagedBcMatchesFaceCount() );
}

// Claim must not depend on EnsureParallelFaceArraysSized: StageNewFace alone
// is enough for a consistent claim.
TEST( FaceSolverStageTest, ClaimDoesNotRequireEnsureParallelFaceArraysSized )
{
    FaceSolver solver;
    solver.StageNewFace( IntField{ 10, 11, 12 }, /*faceType=*/5, INVALID_INDEX );
    EXPECT_TRUE( solver.StagedBcMatchesFaceCount() );
    EXPECT_NO_THROW( solver.ClaimStagedFaceAsLeft( 0, 4 ) );
    EXPECT_EQ( solver.GetFaceTopo().GetLeftCells()[ 0 ], 4 );
}

TEST( FaceSolverStageTest, AttachStagedFaceAsRightSetsNeighbor )
{
    FaceSolver solver;
    solver.StageNewFace( IntField{ 1, 2, 3 }, /*faceType=*/5, INVALID_INDEX );
    solver.ClaimStagedFaceAsLeft( 0, /*cellId=*/3 );
    solver.AttachStagedFaceAsRight( 0, /*cellId=*/8 );

    const FaceTopo & topo = solver.GetFaceTopo();
    EXPECT_EQ( topo.GetLeftCells()[ 0 ], 3 );
    EXPECT_EQ( topo.GetRightCells()[ 0 ], 8 );
    EXPECT_EQ( topo.GetFaceFlags()[ 0 ], FaceSolver::kFaceClaimed );
}

TEST( FaceSolverStageTest, ClaimStagedFaceAsLeftRejectsDoubleClaim )
{
    FaceSolver solver;
    solver.StageNewFace( IntField{ 1, 2, 3 }, /*faceType=*/5, INVALID_INDEX );
    solver.ClaimStagedFaceAsLeft( 0, 1 );
    EXPECT_THROW( solver.ClaimStagedFaceAsLeft( 0, 2 ), std::logic_error );
}

TEST( FaceSolverStageTest, AttachStagedFaceAsRightRejectsUnclaimed )
{
    FaceSolver solver;
    solver.StageNewFace( IntField{ 1, 2, 3 }, /*faceType=*/5, INVALID_INDEX );
    EXPECT_THROW( solver.AttachStagedFaceAsRight( 0, 1 ), std::logic_error );
}

TEST( FaceSolverStageTest, SetStagedFaceBcRejectsOutOfRange )
{
    FaceSolver solver;
    EXPECT_THROW( solver.SetStagedFaceBc( 0, 1, 1 ), std::out_of_range );
}

// CSR polygon staging without CgnsSection: two faces, one duplicated.
TEST( FaceSolverStageTest, StagePolygonFacesFromConnectivityDedupsAndAligns )
{
    FaceSolver solver;
    // face0: nodes 1,2,3  face1: nodes 4,5,6  face2: same as face0
    CgIntField offsets{ 0, 3, 6, 9 };
    CgIntField connectivity{ 1, 2, 3, 4, 5, 6, 1, 2, 3 };
    const int faceType = 7;

    solver.StagePolygonFacesFromConnectivity( 3, offsets, connectivity, faceType );

    EXPECT_TRUE( solver.StagedBcMatchesFaceCount() );
    EXPECT_EQ( solver.GetFaceTopo().GetFaces().size(), static_cast< std::size_t >( 2 ) );
    EXPECT_EQ( solver.GetFaceTopo().GetFaceTypes()[ 0 ], faceType );
    EXPECT_EQ( solver.GetFaceTopo().GetLeftCells()[ 0 ], INVALID_INDEX );
    EXPECT_EQ( solver.GetFaceTopo().GetFaceFlags()[ 0 ], FaceSolver::kFaceUnclaimed );

    solver.ClaimStagedFaceAsLeft( 0, /*cellId=*/10 );
    solver.ClaimStagedFaceAsLeft( 1, /*cellId=*/11 );
    EXPECT_EQ( solver.GetFaceTopo().GetLeftCells()[ 0 ], 10 );
    EXPECT_EQ( solver.GetFaceTopo().GetLeftCells()[ 1 ], 11 );
}

TEST( FaceSolverStageTest, StagePolygonFacesFromConnectivityRejectsBadOffsets )
{
    FaceSolver solver;
    CgIntField offsets{ 0, 2 }; // claims 1 element but end index 2 past connectivity
    CgIntField connectivity{ 1 };
    EXPECT_THROW(
        solver.StagePolygonFacesFromConnectivity( 1, offsets, connectivity, 7 ),
        std::runtime_error );
}

TEST( FaceSolverStageTest, StagePolygonFacesFromConnectivityRejectsIncompleteOffsetTable )
{
    FaceSolver solver;
    CgIntField offsets{ 0 }; // need nElement+1 entries
    CgIntField connectivity{};
    EXPECT_THROW(
        solver.StagePolygonFacesFromConnectivity( 1, offsets, connectivity, 7 ),
        std::runtime_error );
}

// Two polyhedra share one face: first claims left, second attaches right.
TEST( FaceSolverStageTest, ClaimPolyhedronCellsAssignsLeftThenRight )
{
    FaceSolver solver;
    // Three faces staged as unclaimed.
    CgIntField faceOffsets{ 0, 3, 6, 9 };
    CgIntField faceConn{ 1, 2, 3, 4, 5, 6, 7, 8, 9 };
    solver.StagePolygonFacesFromConnectivity( 3, faceOffsets, faceConn, /*faceType=*/7 );

    // Cell 0 uses faces 0,1; cell 1 uses faces 1,2 (shared face 1).
    CgIntField cellOffsets{ 0, 2, 4 };
    CgIntField cellFaceRefs{ 0, 1, 1, 2 };
    solver.ClaimPolyhedronCellsFromConnectivity( 2, cellOffsets, cellFaceRefs );

    const FaceTopo & topo = solver.GetFaceTopo();
    EXPECT_EQ( topo.GetLeftCells()[ 0 ], 0 );
    EXPECT_EQ( topo.GetRightCells()[ 0 ], INVALID_INDEX );
    EXPECT_EQ( topo.GetLeftCells()[ 1 ], 0 );
    EXPECT_EQ( topo.GetRightCells()[ 1 ], 1 );
    EXPECT_EQ( topo.GetLeftCells()[ 2 ], 1 );
    EXPECT_EQ( topo.GetRightCells()[ 2 ], INVALID_INDEX );
    EXPECT_EQ( topo.GetFaceFlags()[ 1 ], FaceSolver::kFaceClaimed );
}

TEST( FaceSolverStageTest, ClaimPolyhedronCellsRejectsThirdCellOnSameFace )
{
    FaceSolver solver;
    CgIntField faceOffsets{ 0, 3 };
    CgIntField faceConn{ 1, 2, 3 };
    solver.StagePolygonFacesFromConnectivity( 1, faceOffsets, faceConn, 7 );

    // Three cells all reference face 0.
    CgIntField cellOffsets{ 0, 1, 2, 3 };
    CgIntField cellFaceRefs{ 0, 0, 0 };
    EXPECT_THROW(
        solver.ClaimPolyhedronCellsFromConnectivity( 3, cellOffsets, cellFaceRefs ),
        std::runtime_error );
}

TEST( FaceSolverStageTest, ClaimPolyhedronCellsRejectsDuplicateFaceInOneCell )
{
    FaceSolver solver;
    CgIntField faceOffsets{ 0, 3, 6 };
    CgIntField faceConn{ 1, 2, 3, 4, 5, 6 };
    solver.StagePolygonFacesFromConnectivity( 2, faceOffsets, faceConn, 7 );

    CgIntField cellOffsets{ 0, 2 };
    CgIntField cellFaceRefs{ 0, 0 }; // same face twice in one polyhedron
    EXPECT_THROW(
        solver.ClaimPolyhedronCellsFromConnectivity( 1, cellOffsets, cellFaceRefs ),
        std::runtime_error );
}

TEST( FaceSolverStageTest, ClaimPolyhedronCellsRejectsOutOfRangeFaceId )
{
    FaceSolver solver;
    CgIntField faceOffsets{ 0, 3 };
    CgIntField faceConn{ 1, 2, 3 };
    solver.StagePolygonFacesFromConnectivity( 1, faceOffsets, faceConn, 7 );

    CgIntField cellOffsets{ 0, 1 };
    CgIntField cellFaceRefs{ 5 }; // only face 0 exists
    EXPECT_THROW(
        solver.ClaimPolyhedronCellsFromConnectivity( 1, cellOffsets, cellFaceRefs ),
        std::runtime_error );
}

TEST( FaceSolverStageTest, AssignBcToMatchingFacesBoundaryOnly )
{
    FaceSolver solver;
    // Face 0: boundary (right INVALID). Face 1: interior (right set).
    solver.StageNewFace( IntField{ 1, 2, 3 }, /*faceType=*/5, /*leftCellId=*/0 );
    solver.StageNewFace( IntField{ 2, 3, 4 }, /*faceType=*/5, /*leftCellId=*/0 );
    solver.GetFaceTopo().GetRightCells()[ 1 ] = 1;

    IntSet bcVertex{ 1, 2, 3, 4 };
    const int assigned = solver.AssignBcToMatchingFaces(
        bcVertex, /*bcType=*/10, /*bcNameId=*/3, /*boundaryFacesOnly=*/true );

    EXPECT_EQ( assigned, 1 );
    EXPECT_EQ( solver.GetStagedBcType( 0 ), 10 );
    EXPECT_EQ( solver.GetStagedBcNameId( 0 ), 3 );
    EXPECT_EQ( solver.GetStagedBcType( 1 ), INVALID_INDEX );
}

TEST( FaceSolverStageTest, AssignBcToMatchingFacesIgnoresAlreadySetBc )
{
    FaceSolver solver;
    solver.StageNewFace( IntField{ 1, 2, 3 }, 5, 0 );
    solver.SetStagedFaceBc( 0, 99, 1 );

    IntSet bcVertex{ 1, 2, 3 };
    const int assigned = solver.AssignBcToMatchingFaces(
        bcVertex, /*bcType=*/10, /*bcNameId=*/2, /*boundaryFacesOnly=*/true );

    EXPECT_EQ( assigned, 0 );
    EXPECT_EQ( solver.GetStagedBcType( 0 ), 99 );
}

TEST( FaceSolverStageTest, AssignBcToMatchingFacesDetailModeMatchesInterior )
{
    FaceSolver solver;
    solver.StageNewFace( IntField{ 1, 2, 3 }, 5, 0 );
    solver.GetFaceTopo().GetRightCells()[ 0 ] = 1; // interior

    IntSet bcVertex{ 1, 2, 3 };
    const int assigned = solver.AssignBcToMatchingFaces(
        bcVertex, /*bcType=*/7, /*bcNameId=*/4, /*boundaryFacesOnly=*/false );

    EXPECT_EQ( assigned, 1 );
    EXPECT_EQ( solver.GetStagedBcType( 0 ), 7 );
}

TEST( FaceSolverStageTest, BindCellToExistingFaceFillsLeftThenRight )
{
    FaceSolver solver;
    solver.StageNewFace( IntField{ 1, 2, 3 }, 5, INVALID_INDEX );

    solver.BindCellToExistingFace( 0, /*cellId=*/3 );
    EXPECT_EQ( solver.GetFaceTopo().GetLeftCells()[ 0 ], 3 );
    EXPECT_EQ( solver.GetFaceTopo().GetRightCells()[ 0 ], INVALID_INDEX );

    solver.BindCellToExistingFace( 0, /*cellId=*/8 );
    EXPECT_EQ( solver.GetFaceTopo().GetLeftCells()[ 0 ], 3 );
    EXPECT_EQ( solver.GetFaceTopo().GetRightCells()[ 0 ], 8 );
}

TEST( FaceSolverStageTest, BindCellToExistingFaceIgnoresSameLeftCell )
{
    FaceSolver solver;
    solver.StageNewFace( IntField{ 1, 2, 3 }, 5, /*leftCellId=*/3 );

    solver.BindCellToExistingFace( 0, /*cellId=*/3 );
    EXPECT_EQ( solver.GetFaceTopo().GetLeftCells()[ 0 ], 3 );
    EXPECT_EQ( solver.GetFaceTopo().GetRightCells()[ 0 ], INVALID_INDEX );
}

TEST( FaceSolverStageTest, BindCellToExistingFaceRejectsOutOfRange )
{
    FaceSolver solver;
    EXPECT_THROW( solver.BindCellToExistingFace( 0, 1 ), std::out_of_range );
}

TEST( FaceSolverStageTest, StagedBcMatchesFaceCountRequiresChildFid )
{
    FaceSolver solver;
    solver.StageNewFace( IntField{ 1, 2, 3 }, 5, 0 );
    ASSERT_TRUE( solver.StagedBcMatchesFaceCount() );

    solver.DropLastChildFidEntryForTest();
    EXPECT_FALSE( solver.StagedBcMatchesFaceCount() );
}

TEST( FaceSolverStageTest, AppendBoundaryBcToRecordCopiesBoundaryOnly )
{
    FaceSolver solver;
    solver.StageNewFace( IntField{ 1, 2, 3 }, 5, 0 ); // boundary: right INVALID
    solver.StageNewFace( IntField{ 2, 3, 4 }, 5, 0 );
    solver.GetFaceTopo().GetRightCells()[ 1 ] = 1; // interior
    solver.SetStagedFaceBc( 0, 20, 5 );
    solver.SetStagedFaceBc( 1, 30, 6 );

    BcRecord record;
    const int nBFaces = solver.AppendBoundaryBcToRecord( record );

    EXPECT_EQ( nBFaces, 1 );
    ASSERT_EQ( record.bcType.size(), static_cast< std::size_t >( 1 ) );
    EXPECT_EQ( record.bcType[ 0 ], 20 );
    EXPECT_EQ( record.bcNameId[ 0 ], 5 );
}

TEST( FaceSolverStageTest, AppendBoundaryBcToRecordRequiresAlignment )
{
    FaceSolver solver;
    solver.StageNewFace( IntField{ 1, 2, 3 }, 5, 0 );
    solver.DropLastChildFidEntryForTest();

    BcRecord record;
    EXPECT_THROW( solver.AppendBoundaryBcToRecord( record ), std::logic_error );
}

TEST( FaceSolverStageTest, GetStagedBcRejectsOutOfRange )
{
    FaceSolver solver;
    EXPECT_THROW( solver.GetStagedBcType( 0 ), std::out_of_range );
    EXPECT_THROW( solver.GetStagedBcNameId( 0 ), std::out_of_range );

    solver.StageNewFace( IntField{ 1, 2, 3 }, 5, 0 );
    EXPECT_EQ( solver.GetStagedBcType( 0 ), INVALID_INDEX );
    EXPECT_THROW( solver.GetStagedBcType( 1 ), std::out_of_range );
}

TEST( FaceSolverStageTest, GetChildFaceIdsStartsEmptyAndRejectsOutOfRange )
{
    FaceSolver solver;
    EXPECT_THROW( solver.GetChildFaceIds( 0 ), std::out_of_range );

    solver.StageNewFace( IntField{ 1, 2, 3 }, 5, 0 );
    EXPECT_EQ( solver.GetChildFidSlotCount(), static_cast< std::size_t >( 1 ) );
    EXPECT_TRUE( solver.GetChildFaceIds( 0 ).empty() );
}

TEST( FaceSolverStageTest, CountFacesWithoutChildrenAllUnsplitAfterStage )
{
    FaceSolver solver;
    solver.StageNewFace( IntField{ 1, 2, 3 }, 5, 0 );
    solver.StageNewFace( IntField{ 2, 3, 4 }, 5, 1 );
    EXPECT_EQ( solver.CountFacesWithoutChildren(), 2 );
    EXPECT_EQ( solver.GetNSimpleFace(), 2 );
}

TEST( FaceSolverStageTest, CountFacesWithoutChildrenRejectsMisalignedChildFid )
{
    FaceSolver solver;
    solver.StageNewFace( IntField{ 1, 2, 3 }, 5, 0 );
    solver.DropLastChildFidEntryForTest();
    EXPECT_THROW( solver.CountFacesWithoutChildren(), std::logic_error );
}
