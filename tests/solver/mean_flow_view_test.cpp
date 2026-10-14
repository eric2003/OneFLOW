#include <gtest/gtest.h>
#include "MeanFlowView.h"
#include "SolverDef.h"

using namespace ONEFLOW;

// ---------------------------------------------------------------------------
// Mean-flow coupling is always NS_SOLVER (loose coupling contract).
// Turbulence must not bind mean-flow via ambient SolverState.
// ---------------------------------------------------------------------------

TEST( MeanFlowView, SolverTypeIsNsSolver )
{
    EXPECT_EQ( MeanFlowView::SolverType(), NS_SOLVER );
    EXPECT_NE( MeanFlowView::SolverType(), TURB_SOLVER );
    EXPECT_NE( MeanFlowView::SolverType(), INC_NS_SOLVER );
}

TEST( MeanFlowView, DefaultPointersAreNull )
{
    MeanFlowView view;
    EXPECT_EQ( view.q, nullptr );
    EXPECT_EQ( view.dqdx, nullptr );
    EXPECT_EQ( view.dqdy, nullptr );
    EXPECT_EQ( view.dqdz, nullptr );
}
