#pragma once

#include "AccelViews.h"

BeginNameSpace( ONEFLOW )

class FluxBackend;

// For one face, primitive values are [rho, u, p] for 3 equations and
// [rho, u, v, w, p] for 5 equations. Face batches use equation-major input
// and output: [rho, rho*u, rho*E] or [rho, rho*u, rho*v, rho*w, rho*E].
void PrimitiveToConserved(
    const Real * primitive, int nEquations, Real gamma, Real * conserved );

struct PrimitiveFaceStateView
{
    int nFaces = 0;
    int nEquations = 0;
    const Real * primitiveLeft = nullptr;
    const Real * primitiveRight = nullptr;
    const Real * primitiveLeftComponents[ 5 ] = {};
    const Real * primitiveRightComponents[ 5 ] = {};
    const Real * xNormal = nullptr;
    const Real * yNormal = nullptr;
    const Real * zNormal = nullptr;
    const Real * meshVelocityNormal = nullptr;
    const Real * faceArea = nullptr;
    Real gamma = 1.4;
};

void PackPrimitiveFaceStates(
    const Real * primitiveLeft,
    const Real * primitiveRight,
    int nFaces,
    int nEquations,
    Real gamma,
    Real * conservedLeft,
    Real * conservedRight );

class EulerCpuAdapter
{
public:
    // The primitive-to-conserved conversion and equation-major pack are
    // execution-space independent. Production callers inject CPU or HIP;
    // keeping this seam here prevents the two paths from drifting.
    void CalcInvFlux(
        const PrimitiveFaceStateView & state,
        FaceFluxView & flux,
        FluxBackend & backend,
        int scheme = 0 ) const;

    void CalcInvFlux(
        const PrimitiveFaceStateView & state,
        FaceFluxView & flux,
        int scheme = 0 ) const;

    void AddFaceFlux(
        const FaceFluxView & flux,
        const FaceConnectivityView & connectivity,
        ResidualView & residual ) const;
};

EndNameSpace
