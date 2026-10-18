#include <gtest/gtest.h>
#include <memory>

#include "BcRecord.h"
#include "Boundary.h"
#include "DataBook.h"
#include "DataBaseIO.h"
#include "FaceTopo.h"
#include "InterFace.h"
#include "IFaceLink.h"
#include "GridHandles.h"
#include "UnsGrid.h"
#include "NodeMesh.h"
#include "FaceSearch.h"
#include "Tolerance.h"


namespace
{
class ScopedPointTolerance
{
public:
    ScopedPointTolerance()
        : previousTolerance( ONEFLOW::Tolerance::tol ),
          previousFlag( ONEFLOW::Tolerance::flag )
    {
        ONEFLOW::Tolerance::tol = 1.0e-8;
        ONEFLOW::Tolerance::flag = true;
    }

    ~ScopedPointTolerance()
    {
        ONEFLOW::Tolerance::tol = previousTolerance;
        ONEFLOW::Tolerance::flag = previousFlag;
    }

private:
    ONEFLOW::Real previousTolerance;
    bool previousFlag;
};

void ConfigureBoundaryTopology( ONEFLOW::UnsGrid & grid )
{
    auto & faceTopo = grid.GetFaceTopo();
    faceTopo.SetNBFaces( 4 );

    auto & bcRecord = faceTopo.GetBcRecord();
    bcRecord.bcType = {
        ONEFLOW::BC::INTERFACE,
        ONEFLOW::BC::SOLID_SURFACE,
        ONEFLOW::BC::INTERFACE,
        ONEFLOW::BC::INTERFACE
    };
    bcRecord.bcNameId = { 10, 20, 30, 40 };

    grid.nIFaces = 3;
    grid.interFace->Set( 3, &grid );
    grid.interFace->zoneId = { 1, 1, 2 };
    grid.interFace->localInterfaceId = { 4, 7, 0 };
    grid.interFace->i2b = { 0, 2, 3 };
}
}

TEST( UnsGridBoundaryTopologyRoundTripTest, DataBookPreservesBoundaryAndInterfaceMappings )
{
    ONEFLOW::UnsGrid source;
    source.Init();
    source.nFaces = 4;
    ConfigureBoundaryTopology( source );

    ONEFLOW::DataBook dataBook( 16 );
    ASSERT_NO_THROW( source.WriteBoundaryTopology( &dataBook ) );

    dataBook.MoveToBegin();

    ONEFLOW::UnsGrid restored;
    restored.Init();
    restored.nFaces = 4;
    ASSERT_NO_THROW( restored.ReadBoundaryTopology( &dataBook ) );

    const auto & sourceBc = source.GetFaceTopo().GetBcRecord();
    const auto & restoredBc = restored.GetFaceTopo().GetBcRecord();
    ASSERT_EQ( restored.GetFaceTopo().GetNBFaces(), 4u );
    EXPECT_EQ( restoredBc.bcType, sourceBc.bcType );
    EXPECT_EQ( restoredBc.bcNameId, sourceBc.bcNameId );
    EXPECT_EQ( restored.nIFaces, 3 );
    ASSERT_NE( restored.interFace, nullptr );
    EXPECT_EQ( restored.interFace->zoneId, source.interFace->zoneId );
    EXPECT_EQ( restored.interFace->localInterfaceId, source.interFace->localInterfaceId );
    EXPECT_EQ( restored.interFace->i2b, ( ONEFLOW::IntField{ 0, 2, 3 } ) );

    for ( int interfaceId = 0; interfaceId < restored.nIFaces; ++ interfaceId )
    {
        const int boundaryId = restored.interFace->i2b[ interfaceId ];
        ASSERT_GE( boundaryId, 0 );
        ASSERT_LT( boundaryId, static_cast< int >( restoredBc.bcType.size() ) );
        EXPECT_TRUE( ONEFLOW::BC::IsInterfaceBc( restoredBc.bcType[ boundaryId ] ) );
    }
}

TEST( UnsGridBoundaryTopologyRoundTripTest, SplitInterfaceRecordsRegenerateBoundaryMap )
{
    ONEFLOW::BcRecord bcRecord;
    bcRecord.bcType = {
        ONEFLOW::BC::INTERFACE,
        ONEFLOW::BC::INTERFACE,
        ONEFLOW::BC::SOLID_SURFACE,
        ONEFLOW::BC::INTERFACE,
        ONEFLOW::BC::INTERFACE
    };
    bcRecord.bcNameId = { 10, 10, 20, 30, 30 };

    ONEFLOW::InterFace interFace;
    interFace.Set( 4 );
    bcRecord.GenerateI2B( interFace );

    EXPECT_EQ( interFace.i2b, ( ONEFLOW::IntField{ 0, 1, 3, 4 } ) );
    for ( int interfaceId = 0; interfaceId < interFace.nIFaces; ++ interfaceId )
    {
        const int boundaryId = interFace.i2b[ interfaceId ];
        ASSERT_GE( boundaryId, 0 );
        ASSERT_LT( boundaryId, static_cast< int >( bcRecord.bcType.size() ) );
        EXPECT_TRUE( ONEFLOW::BC::IsInterfaceBc( bcRecord.bcType[ boundaryId ] ) );
    }

    EXPECT_EQ( bcRecord.bcNameId[ interFace.i2b[ 0 ] ], 10 );
    EXPECT_EQ( bcRecord.bcNameId[ interFace.i2b[ 1 ] ], 10 );
    EXPECT_EQ( bcRecord.bcNameId[ interFace.i2b[ 2 ] ], 30 );
    EXPECT_EQ( bcRecord.bcNameId[ interFace.i2b[ 3 ] ], 30 );
}

TEST( UnsGridBoundaryTopologyRoundTripTest, InterfaceSplitKeepsLocalGlobalAndBoundaryMappingsAligned )
{
    ScopedPointTolerance pointTolerance;
    ONEFLOW::Grids grids;
    auto ownedGrid = std::make_unique< ONEFLOW::UnsGrid >();
    ownedGrid->Init();
    ownedGrid->id = 0;
    ownedGrid->nCells = 1;
    ownedGrid->nFaces = 3;
    ownedGrid->nNodes = 2;
    ownedGrid->nodeMesh->CreateNodes( 2 );
    ownedGrid->nodeMesh->xN[ 0 ] = 0.0;
    ownedGrid->nodeMesh->xN[ 1 ] = 1.0;
    ownedGrid->nodeMesh->CalcMinMaxBox();
    ownedGrid->GetFaceTopo().GetFaceTypes().resize( 3 );
    ownedGrid->GetFaceTopo().GetFaces() = {
        { 0, 1 }, { 0, 1 }, { 0, 1 }
    };
    ownedGrid->GetFaceTopo().SetNBFaces( 3 );
    ownedGrid->GetFaceTopo().GetBcRecord().bcType = {
        ONEFLOW::BC::INTERFACE,
        ONEFLOW::BC::SOLID_SURFACE,
        ONEFLOW::BC::INTERFACE
    };
    ownedGrid->GetFaceTopo().GetBcRecord().bcNameId = { 10, 20, 30 };
    ownedGrid->nIFaces = 2;
    ownedGrid->interFace->Set( 2, ownedGrid.get() );
    grids.push_back( std::move( ownedGrid ) );
    ONEFLOW::UnsGrid & grid = static_cast< ONEFLOW::UnsGrid & >( * grids[ 0 ] );

    ONEFLOW::IFaceLink link( grids );
    link.l2g[ 0 ] = { 0, 1 };
    link.gI2Zid = { { 0, 1 }, { 0, 1 } };
    link.g2l = { { 0, 0 }, { 1, 1 } };
    link.InitNewLgMapping();

    link.face_search->interfaceSplit.childFaceIdsOfParent.resize( 4 );
    link.face_search->interfaceSplit.childFaceIdsOfParent[ 0 ] = { 2, 3 };
    link.face_search->interfaceSplit.childFaceNodeIds.resize( 4 );
    link.face_search->interfaceSplit.childFaceNodeFlags.resize( 4 );
    link.face_search->interfaceSplit.childFaceNodeIds[ 2 ] = { 0 };
    link.face_search->interfaceSplit.childFaceNodeFlags[ 2 ] = { 1 };
    link.face_search->interfaceSplit.childFaceNodeIds[ 3 ] = { 0 };
    link.face_search->interfaceSplit.childFaceNodeFlags[ 3 ] = { 1 };

    grid.GetFaceTopo().ModifyBoundaryInformation( link );

    ASSERT_EQ( link.l2gNew[ 0 ], ( ONEFLOW::IntField{ 2, 3, 1 } ) );
    link.UpdateLgMapping();
    grid.UpdateOtherTopologyTerm( link );

    const auto & bcRecord = grid.GetFaceTopo().GetBcRecord();
    ASSERT_EQ( bcRecord.bcType.size(), 4u );
    EXPECT_EQ( bcRecord.bcType, ( ONEFLOW::IntField{
        ONEFLOW::BC::INTERFACE,
        ONEFLOW::BC::INTERFACE,
        ONEFLOW::BC::SOLID_SURFACE,
        ONEFLOW::BC::INTERFACE
    } ) );
    EXPECT_EQ( bcRecord.bcNameId, ( ONEFLOW::IntField{ 10, 10, 20, 30 } ) );
    EXPECT_EQ( grid.nIFaces, 3 );
    EXPECT_EQ( grid.interFace->i2b, ( ONEFLOW::IntField{ 0, 1, 3 } ) );
    EXPECT_EQ( link.l2g[ 0 ], ( ONEFLOW::IntField{ 2, 3, 1 } ) );

    for ( int interfaceId = 0; interfaceId < grid.nIFaces; ++ interfaceId )
    {
        const int boundaryId = grid.interFace->i2b[ interfaceId ];
        ASSERT_GE( boundaryId, 0 );
        ASSERT_LT( boundaryId, static_cast< int >( bcRecord.bcType.size() ) );
        EXPECT_TRUE( ONEFLOW::BC::IsInterfaceBc( bcRecord.bcType[ boundaryId ] ) );
    }
}
