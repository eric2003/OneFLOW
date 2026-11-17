#include "EulerInvFluxCapability.h"

#include <gtest/gtest.h>

namespace
{

using namespace ONEFLOW;

EulerInvFluxCapabilityRequest SupportedRequest()
{
    EulerInvFluxCapabilityRequest request;
    request.hipBuildEnabled = true;
    request.runtimeInitialized = true;
    request.hipBackendSelected = true;
    request.supportedSolver = true;
    request.localZoneCount = 1;
    request.gridLevel = 0;
    request.gridCount = 1;
    request.nEquations = 5;
    request.limiterEquations = 5;
    request.laxFriedrichsScheme = true;
    return request;
}

TEST( EulerInvFluxCapability, AcceptsOnlyFullySupportedRequest )
{
    const auto decision = EvaluateEulerInvFluxCapability(
        SupportedRequest() );
    EXPECT_TRUE( decision.enabled );
    EXPECT_EQ( decision.reason, EulerInvFluxCapabilityReason::Supported );
    EXPECT_STREQ(
        EulerInvFluxCapabilityReasonName( decision.reason ), "supported" );
}

TEST( EulerInvFluxCapability, RejectsEveryUnsupportedRequestWithReason )
{
    struct Case
    {
        EulerInvFluxCapabilityReason reason;
        void ( * mutate )( EulerInvFluxCapabilityRequest & );
    };

    const Case cases[] = {
        { EulerInvFluxCapabilityReason::HipNotBuilt,
          []( auto & request ) { request.hipBuildEnabled = false; } },
        { EulerInvFluxCapabilityReason::RuntimeNotInitialized,
          []( auto & request ) { request.runtimeInitialized = false; } },
        { EulerInvFluxCapabilityReason::HipBackendNotSelected,
          []( auto & request ) { request.hipBackendSelected = false; } },
        { EulerInvFluxCapabilityReason::UnsupportedSolver,
          []( auto & request ) { request.supportedSolver = false; } },
        { EulerInvFluxCapabilityReason::MultipleLocalZones,
          []( auto & request ) { request.localZoneCount = 2; } },
        { EulerInvFluxCapabilityReason::NonFinestGrid,
          []( auto & request ) { request.gridLevel = 1; } },
        { EulerInvFluxCapabilityReason::MultipleGridLevels,
          []( auto & request ) { request.gridCount = 2; } },
        { EulerInvFluxCapabilityReason::UnsupportedEquationCount,
          []( auto & request ) { request.nEquations = 3; } },
        { EulerInvFluxCapabilityReason::UnsupportedLimiterEquationCount,
          []( auto & request ) { request.limiterEquations = 3; } },
        { EulerInvFluxCapabilityReason::UnsupportedInviscidScheme,
          []( auto & request ) { request.laxFriedrichsScheme = false; } },
        { EulerInvFluxCapabilityReason::InterfaceExchangeEnabled,
          []( auto & request ) { request.hasInterfaceExchange = true; } },
    };

    for ( const auto & testCase : cases )
    {
        auto request = SupportedRequest();
        testCase.mutate( request );
        const auto decision = EvaluateEulerInvFluxCapability( request );
        EXPECT_FALSE( decision.enabled );
        EXPECT_EQ( decision.reason, testCase.reason );
        EXPECT_STRNE(
            EulerInvFluxCapabilityReasonName( decision.reason ),
            "supported" );
    }
}

TEST( EulerInvFluxCapability, ReportsFirstBlockingReasonDeterministically )
{
    auto request = SupportedRequest();
    request.runtimeInitialized = false;
    request.hipBackendSelected = false;
    request.nEquations = 3;

    const auto decision = EvaluateEulerInvFluxCapability( request );
    EXPECT_FALSE( decision.enabled );
    EXPECT_EQ(
        decision.reason,
        EulerInvFluxCapabilityReason::RuntimeNotInitialized );
}

TEST( EulerInvFluxCapability, AcceptsSupportedInterfaceExchange )
{
    auto request = SupportedRequest();
    request.hasInterfaceExchange = true;
    request.backendSupportsInterfaceExchange = true;

    const auto decision = EvaluateEulerInvFluxCapability( request );
    EXPECT_TRUE( decision.enabled );
    EXPECT_EQ(
        decision.reason,
        EulerInvFluxCapabilityReason::Supported );
}

TEST( EulerInvFluxCapability, AcceptsMultipleLocalZonesWithSupportedExchange )
{
    auto request = SupportedRequest();
    request.localZoneCount = 4;
    request.hasInterfaceExchange = true;
    request.backendSupportsInterfaceExchange = true;

    const auto decision = EvaluateEulerInvFluxCapability( request );
    EXPECT_TRUE( decision.enabled );
    EXPECT_EQ(
        decision.reason,
        EulerInvFluxCapabilityReason::Supported );
}

} // namespace
