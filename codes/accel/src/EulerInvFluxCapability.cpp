#include "EulerInvFluxCapability.h"

BeginNameSpace( ONEFLOW )

EulerInvFluxCapabilityDecision EvaluateEulerInvFluxCapability(
    const EulerInvFluxCapabilityRequest & request )
{
    EulerInvFluxCapabilityDecision decision;

    if ( ! request.hipBuildEnabled )
    {
        decision.reason = EulerInvFluxCapabilityReason::HipNotBuilt;
        return decision;
    }
    if ( ! request.runtimeInitialized )
    {
        decision.reason = EulerInvFluxCapabilityReason::RuntimeNotInitialized;
        return decision;
    }
    if ( ! request.hipBackendSelected )
    {
        decision.reason = EulerInvFluxCapabilityReason::HipBackendNotSelected;
        return decision;
    }
    if ( ! request.supportedSolver )
    {
        decision.reason = EulerInvFluxCapabilityReason::UnsupportedSolver;
        return decision;
    }
    if ( request.localZoneCount <= 0 )
    {
        decision.reason = EulerInvFluxCapabilityReason::MultipleLocalZones;
        return decision;
    }
    if ( request.gridLevel != 0 )
    {
        decision.reason = EulerInvFluxCapabilityReason::NonFinestGrid;
        return decision;
    }
    if ( request.gridCount != 1 )
    {
        decision.reason = EulerInvFluxCapabilityReason::MultipleGridLevels;
        return decision;
    }
    if ( request.nEquations != 5 )
    {
        decision.reason = EulerInvFluxCapabilityReason::UnsupportedEquationCount;
        return decision;
    }
    if ( request.limiterEquations != 5 )
    {
        decision.reason =
            EulerInvFluxCapabilityReason::UnsupportedLimiterEquationCount;
        return decision;
    }
    if ( ! request.laxFriedrichsScheme )
    {
        decision.reason =
            EulerInvFluxCapabilityReason::UnsupportedInviscidScheme;
        return decision;
    }
    if ( request.hasInterfaceExchange
         && ! request.backendSupportsInterfaceExchange )
    {
        decision.reason =
            EulerInvFluxCapabilityReason::InterfaceExchangeEnabled;
        return decision;
    }
    if ( request.localZoneCount > 1
         && ( ! request.hasInterfaceExchange
              || ! request.backendSupportsInterfaceExchange ) )
    {
        decision.reason = EulerInvFluxCapabilityReason::MultipleLocalZones;
        return decision;
    }

    decision.enabled = true;
    decision.reason = EulerInvFluxCapabilityReason::Supported;
    return decision;
}

const char * EulerInvFluxCapabilityReasonName(
    EulerInvFluxCapabilityReason reason )
{
    switch ( reason )
    {
    case EulerInvFluxCapabilityReason::Supported:
        return "supported";
    case EulerInvFluxCapabilityReason::HipNotBuilt:
        return "hip_not_built";
    case EulerInvFluxCapabilityReason::RuntimeNotInitialized:
        return "runtime_not_initialized";
    case EulerInvFluxCapabilityReason::HipBackendNotSelected:
        return "hip_backend_not_selected";
    case EulerInvFluxCapabilityReason::UnsupportedSolver:
        return "unsupported_solver";
    case EulerInvFluxCapabilityReason::MultipleLocalZones:
        return "multiple_local_zones";
    case EulerInvFluxCapabilityReason::NonFinestGrid:
        return "non_finest_grid";
    case EulerInvFluxCapabilityReason::MultipleGridLevels:
        return "multiple_grid_levels";
    case EulerInvFluxCapabilityReason::UnsupportedEquationCount:
        return "unsupported_equation_count";
    case EulerInvFluxCapabilityReason::UnsupportedLimiterEquationCount:
        return "unsupported_limiter_equation_count";
    case EulerInvFluxCapabilityReason::UnsupportedInviscidScheme:
        return "unsupported_inviscid_scheme";
    case EulerInvFluxCapabilityReason::InterfaceExchangeEnabled:
        return "interface_exchange_enabled";
    }
    return "unknown";
}

EndNameSpace
