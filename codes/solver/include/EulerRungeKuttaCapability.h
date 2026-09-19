#pragma once

#include "SolverDef.h"

BeginNameSpace( ONEFLOW )

enum class EulerRungeKuttaTimeIntegral
{
    ExplicitMultiStage = 1,
    Lusgs = 2,
    Simple = 3
};

enum class EulerRungeKuttaCapabilityReason
{
    Supported,
    UnsupportedSolver,
    MultipleLocalZones,
    NonFinestGrid,
    MultipleGridLevels,
    UnsupportedEquationCount,
    UnsupportedInviscidScheme,
    NonExplicitTimeIntegral,
    ViscousTermsEnabled,
    SourceTermsEnabled,
    LimiterEnabled,
    InterfaceExchangeEnabled,
    BackendAdvanceUnsupported
};

constexpr int EULER_RK_LAX_FRIEDRICHS_SCHEME = 5;

// A solver-side request used before entering the Euler RK fast path. It is
// deliberately independent of MRField and backend implementation details.
struct EulerRungeKuttaCapabilityRequest
{
    int solverType = NS_SOLVER;
    int localZoneCount = 0;
    int gridLevel = -1;
    int gridCount = 0;
    int nEquations = 0;
    int inviscidScheme = 0;
    EulerRungeKuttaTimeIntegral timeIntegral =
        EulerRungeKuttaTimeIntegral::ExplicitMultiStage;
    bool hasViscousTerms = true;
    bool hasSourceTerms = true;
    bool hasLimiter = true;
    bool hasInterfaceExchange = false;
    bool backendSupportsAdvance = false;
};

struct EulerRungeKuttaCapabilityDecision
{
    bool enabled = false;
    EulerRungeKuttaCapabilityReason reason =
        EulerRungeKuttaCapabilityReason::UnsupportedSolver;
};

EulerRungeKuttaCapabilityDecision EvaluateEulerRungeKuttaCapability(
    const EulerRungeKuttaCapabilityRequest & request );

const char * EulerRungeKuttaCapabilityReasonName(
    EulerRungeKuttaCapabilityReason reason );

EndNameSpace
