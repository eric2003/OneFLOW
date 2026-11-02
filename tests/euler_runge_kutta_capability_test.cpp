#include "EulerRungeKuttaCapability.h"

#include <gtest/gtest.h>

namespace
{

using namespace ONEFLOW;

EulerRungeKuttaCapabilityRequest SupportedRequest()
{
    EulerRungeKuttaCapabilityRequest request;
    request.localZoneCount = 1;
    request.gridLevel = 0;
    request.gridCount = 1;
    request.nEquations = 5;
    request.inviscidScheme = EULER_RK_LAX_FRIEDRICHS_SCHEME;
    request.hasViscousTerms = false;
    request.hasSourceTerms = false;
    request.hasLimiter = false;
    request.hasInterfaceExchange = false;
    request.backendSupportsAdvance = true;
    return request;
}

TEST( EulerRungeKuttaCapability, AcceptsOnlyFullySupportedRequest )
{
    const auto decision = EvaluateEulerRungeKuttaCapability(
        SupportedRequest() );
    EXPECT_TRUE( decision.enabled );
    EXPECT_EQ(
        decision.reason,
        EulerRungeKuttaCapabilityReason::Supported );
    EXPECT_STREQ(
        EulerRungeKuttaCapabilityReasonName( decision.reason ), "supported" );
}

TEST( EulerRungeKuttaCapability, RejectsUnsupportedRequestWithReason )
{
    struct Case
    {
        EulerRungeKuttaCapabilityReason reason;
        void ( * mutate )( EulerRungeKuttaCapabilityRequest & );
    };

    const Case cases[] = {
        { EulerRungeKuttaCapabilityReason::UnsupportedSolver,
          []( auto & request ) { request.solverType = CFD_SOLVER; } },
        { EulerRungeKuttaCapabilityReason::MultipleLocalZones,
          []( auto & request ) { request.localZoneCount = 2; } },
        { EulerRungeKuttaCapabilityReason::NonFinestGrid,
          []( auto & request ) { request.gridLevel = 1; } },
        { EulerRungeKuttaCapabilityReason::MultipleGridLevels,
          []( auto & request ) { request.gridCount = 2; } },
        { EulerRungeKuttaCapabilityReason::UnsupportedEquationCount,
          []( auto & request ) { request.nEquations = 3; } },
        { EulerRungeKuttaCapabilityReason::UnsupportedInviscidScheme,
          []( auto & request ) { request.inviscidScheme = 1; } },
        { EulerRungeKuttaCapabilityReason::NonExplicitTimeIntegral,
          []( auto & request ) {
              request.timeIntegral = EulerRungeKuttaTimeIntegral::Lusgs;
          } },
        { EulerRungeKuttaCapabilityReason::ViscousTermsEnabled,
          []( auto & request ) { request.hasViscousTerms = true; } },
        { EulerRungeKuttaCapabilityReason::SourceTermsEnabled,
          []( auto & request ) { request.hasSourceTerms = true; } },
        { EulerRungeKuttaCapabilityReason::LimiterEnabled,
          []( auto & request ) { request.hasLimiter = true; } },
        { EulerRungeKuttaCapabilityReason::InterfaceExchangeEnabled,
          []( auto & request ) { request.hasInterfaceExchange = true; } },
        { EulerRungeKuttaCapabilityReason::BackendAdvanceUnsupported,
          []( auto & request ) { request.backendSupportsAdvance = false; } },
    };

    for ( const auto & testCase : cases )
    {
        auto request = SupportedRequest();
        testCase.mutate( request );
        const auto decision = EvaluateEulerRungeKuttaCapability( request );
        EXPECT_FALSE( decision.enabled );
        EXPECT_EQ( decision.reason, testCase.reason );
        EXPECT_STRNE(
            EulerRungeKuttaCapabilityReasonName( decision.reason ),
            "supported" );
    }
}

} // namespace
