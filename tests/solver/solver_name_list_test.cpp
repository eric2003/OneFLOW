#include <gtest/gtest.h>
#include "SolverNamePolicy.h"
#include "SolverNameList.h"
#include "SolverMap.h"
#include "GridState.h"  // UMESH / SMESH if needed
#include <string>
#include <vector>

using namespace ONEFLOW;

TEST( SolverNameClassSeam, LoadFromBaseNamesSkipsFileAndExpands )
{
    SolverNameClass::Reset();

    StringField base;
    base.push_back( "NsSolver" );
    base.push_back( "TurbSolver" );

    SolverNameClass::LoadFromBaseNames( base );

    StringField & uns = SolverNameClass::GetSolverNames( ONEFLOW::UMESH );
    StringField & str = SolverNameClass::GetSolverNames( ONEFLOW::SMESH );

    ASSERT_EQ( uns.size(), 2u );
    ASSERT_EQ( str.size(), 2u );
    EXPECT_EQ( uns[ 0 ], "UNsSolver" );
    EXPECT_EQ( uns[ 1 ], "UTurbSolver" );
    EXPECT_EQ( str[ 0 ], "SNsSolver" );
    EXPECT_EQ( str[ 1 ], "STurbSolver" );

    // Init must not wipe injected names (flag already true)
    SolverNameClass::Init();
    EXPECT_EQ( SolverNameClass::GetSolverNames( ONEFLOW::UMESH ).size(), 2u );

    SolverNameClass::Reset();
}