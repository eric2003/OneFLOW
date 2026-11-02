#include "EulerRungeKuttaCapability.h"

BeginNameSpace( ONEFLOW )

EulerRungeKuttaCapabilityDecision EvaluateEulerRungeKuttaCapability(
    const EulerRungeKuttaCapabilityRequest & request )
{
    EulerRungeKuttaCapabilityDecision decision;

    if ( request.solverType != NS_SOLVER )
    {
        decision.reason = EulerRungeKuttaCapabilityReason::UnsupportedSolver;
        return decision;
    }
    if ( request.localZoneCount != 1 )
    {
        decision.reason = EulerRungeKuttaCapabilityReason::MultipleLocalZones;
        return decision;
    }
    if ( request.gridLevel != 0 )
    {
        decision.reason = EulerRungeKuttaCapabilityReason::NonFinestGrid;
        return decision;
    }
    if ( request.gridCount != 1 )
    {
        decision.reason = EulerRungeKuttaCapabilityReason::MultipleGridLevels;
        return decision;
    }
    if ( request.nEquations != 5 )
    {
        decision.reason = EulerRungeKuttaCapabilityReason::UnsupportedEquationCount;
        return decision;
    }
    if ( request.inviscidScheme != EULER_RK_LAX_FRIEDRICHS_SCHEME )
    {
        decision.reason = EulerRungeKuttaCapabilityReason::UnsupportedInviscidScheme;
        return decision;
    }
    if ( request.timeIntegral != EulerRungeKuttaTimeIntegral::ExplicitMultiStage )
    {
        decision.reason = EulerRungeKuttaCapabilityReason::NonExplicitTimeIntegral;
        return decision;
    }
    if ( request.hasViscousTerms )
    {
        decision.reason = EulerRungeKuttaCapabilityReason::ViscousTermsEnabled;
        return decision;
    }
    if ( request.hasSourceTerms )
    {
        decision.reason = EulerRungeKuttaCapabilityReason::SourceTermsEnabled;
        return decision;
    }
    if ( request.hasLimiter )
    {
        decision.reason = EulerRungeKuttaCapabilityReason::LimiterEnabled;
        return decision;
    }
    if ( request.hasInterfaceExchange )
    {
        decision.reason = EulerRungeKuttaCapabilityReason::InterfaceExchangeEnabled;
        return decision;
    }
    if ( ! request.backendSupportsAdvance )
    {
        decision.reason = EulerRungeKuttaCapabilityReason::BackendAdvanceUnsupported;
        return decision;
    }

    decision.enabled = true;
    decision.reason = EulerRungeKuttaCapabilityReason::Supported;
    return decision;
}

const char * EulerRungeKuttaCapabilityReasonName(
    EulerRungeKuttaCapabilityReason reason )
{
    switch ( reason )
    {
    case EulerRungeKuttaCapabilityReason::Supported:
        return "supported";
    case EulerRungeKuttaCapabilityReason::UnsupportedSolver:
        return "unsupported_solver";
    case EulerRungeKuttaCapabilityReason::MultipleLocalZones:
        return "multiple_local_zones";
    case EulerRungeKuttaCapabilityReason::NonFinestGrid:
        return "non_finest_grid";
    case EulerRungeKuttaCapabilityReason::MultipleGridLevels:
        return "multiple_grid_levels";
    case EulerRungeKuttaCapabilityReason::UnsupportedEquationCount:
        return "unsupported_equation_count";
    case EulerRungeKuttaCapabilityReason::UnsupportedInviscidScheme:
        return "unsupported_inviscid_scheme";
    case EulerRungeKuttaCapabilityReason::NonExplicitTimeIntegral:
        return "non_explicit_time_integral";
    case EulerRungeKuttaCapabilityReason::ViscousTermsEnabled:
        return "viscous_terms_enabled";
    case EulerRungeKuttaCapabilityReason::SourceTermsEnabled:
        return "source_terms_enabled";
    case EulerRungeKuttaCapabilityReason::LimiterEnabled:
        return "limiter_enabled";
    case EulerRungeKuttaCapabilityReason::InterfaceExchangeEnabled:
        return "interface_exchange_enabled";
    case EulerRungeKuttaCapabilityReason::BackendAdvanceUnsupported:
        return "backend_advance_unsupported";
    }
    return "unknown";
}

EndNameSpace
