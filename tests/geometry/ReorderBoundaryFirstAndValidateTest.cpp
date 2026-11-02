// tests/geometry/ReorderBoundaryFirstAndValidateTest.cpp
// Locks ReorderLink + ValidateBoundaryBcLayout as one post-stage transaction.

#include <gtest/gtest.h>
#include <stdexcept>

#include "Constant.h"
#include "FaceTopo.h"
#include "BcRecord.h"
#include "HXDefine.h"

using namespace ONEFLOW;

namespace {

void PushFace( FaceTopo & topo, IntField nodes, int left, int right, int type = 5 )
{
    topo.GetFaces().push_back( std::move( nodes ) );
    topo.GetFaceTypes().push_back( type );
    topo.GetLeftCells().push_back( left );
    topo.GetRightCells().push_back( right );
}

} // namespace

TEST( ReorderBoundaryFirstAndValidateTest, BoundaryPrefixMatchesBcRecord )
{
    FaceTopo topo;
    // Discovery order: interior, boundary, boundary
    PushFace( topo, IntField{ 0, 1 }, 0, 1 );
    PushFace( topo, IntField{ 1, 2 }, 0, INVALID_INDEX );
    PushFace( topo, IntField{ 2, 3 }, 1, INVALID_INDEX );

    // BcRecord filled in discovery order of boundary faces only
    topo.GetBcRecord().bcType = IntField{ 10, 20 };
    topo.GetBcRecord().bcNameId = IntField{ 1, 2 };

    EXPECT_NO_THROW( ReorderBoundaryFirstAndValidate( topo ) );

    const IntField & right = topo.GetRightCells();
    ASSERT_EQ( right.size(), static_cast< std::size_t >( 3 ) );
    EXPECT_EQ( right[ 0 ], INVALID_INDEX );
    EXPECT_EQ( right[ 1 ], INVALID_INDEX );
    EXPECT_NE( right[ 2 ], INVALID_INDEX );

    EXPECT_EQ( topo.GetBcRecord().bcType.size(), static_cast< std::size_t >( 2 ) );
}

TEST( ReorderBoundaryFirstAndValidateTest, RejectsBcRecordSizeMismatch )
{
    FaceTopo topo;
    PushFace( topo, IntField{ 0, 1 }, 0, INVALID_INDEX );
    PushFace( topo, IntField{ 1, 2 }, 0, 1 );

    topo.GetBcRecord().bcType = IntField{ 10, 20 }; // two BC rows, one boundary face
    topo.GetBcRecord().bcNameId = IntField{ 0, 1 };

    EXPECT_THROW( ReorderBoundaryFirstAndValidate( topo ), std::logic_error );
}
