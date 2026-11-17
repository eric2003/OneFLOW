#pragma once

#include "NamespaceMacros.h"

BeginNameSpace( ONEFLOW )

enum class EulerInvFluxCapabilityReason
{
    Supported,
    HipNotBuilt,
    RuntimeNotInitialized,
    HipBackendNotSelected,
    UnsupportedSolver,
    MultipleLocalZones,
    NonFinestGrid,
    MultipleGridLevels,
    UnsupportedEquationCount,
    UnsupportedLimiterEquationCount,
    UnsupportedInviscidScheme,
    InterfaceExchangeEnabled
};

// Solver-side facts required before the production UNs flux path may enter
// HipFluxBackend. Keeping this request independent of global solver state makes
// every rejection reason deterministic and unit-testable.
struct EulerInvFluxCapabilityRequest
{
    bool hipBuildEnabled = false;
    bool runtimeInitialized = false;
    bool hipBackendSelected = false;
    bool supportedSolver = false;
    int localZoneCount = 0;
    int gridLevel = -1;
    int gridCount = 0;
    int nEquations = 0;
    int limiterEquations = 0;
    bool laxFriedrichsScheme = false;
    // Multi-zone partitions can have one local zone per rank while still
    // requiring cross-rank interface exchange. Keep this explicit so the HIP
    // path cannot infer support from localZoneCount alone.
    bool hasInterfaceExchange = false;
    bool backendSupportsInterfaceExchange = false;
};

struct EulerInvFluxCapabilityDecision
{
    bool enabled = false;
    EulerInvFluxCapabilityReason reason =
        EulerInvFluxCapabilityReason::HipNotBuilt;
};

EulerInvFluxCapabilityDecision EvaluateEulerInvFluxCapability(
    const EulerInvFluxCapabilityRequest & request );

const char * EulerInvFluxCapabilityReasonName(
    EulerInvFluxCapabilityReason reason );

EndNameSpace
