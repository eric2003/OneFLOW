/*---------------------------------------------------------------------------*\
    OneFLOW - LargeScale Multiphysics Scientific Simulation Environment
    Copyright (C) 2017-2026 He Xin and the OneFLOW contributors.
-------------------------------------------------------------------------------
License
    This file is part of OneFLOW.
\*---------------------------------------------------------------------------*/

#include "CpuFluxBackend.h"

#include <cmath>
#include <stdexcept>

namespace ONEFLOW
{
namespace
{

// --- Scalar convection (nEquations == 1) ---

void CalcScalarInvFlux( const FaceStateView & state, FaceFluxView & flux )
{
    for ( int face = 0; face < state.nFaces; ++ face )
    {
        const Real normalVelocity = state.xNormal[ face ];
        const Real positive = 0.5 * ( normalVelocity + std::abs( normalVelocity ) );
        const Real negative = 0.5 * ( normalVelocity - std::abs( normalVelocity ) );
        flux.values[ face ] =
            ( state.qLeft[ face ] * positive + state.qRight[ face ] * negative )
            * state.faceArea[ face ];
    }
}

bool IsBoundaryFace( const FaceConnectivityView & connectivity, int face )
{
    return connectivity.boundaryMask != nullptr
        ? connectivity.boundaryMask[ face ] != 0
        : face < connectivity.nBoundaryFaces;
}

void AddScalarFaceFlux(
    const FaceFluxView & flux,
    const FaceConnectivityView & connectivity,
    ResidualView & residual )
{
    for ( int face = 0; face < connectivity.nFaces; ++ face )
    {
        const int left = connectivity.leftCell[ face ];
        residual.values[ left ] -= flux.values[ face ];
        if ( ! IsBoundaryFace( connectivity, face ) )
        {
            residual.values[ connectivity.rightCell[ face ] ] += flux.values[ face ];
        }
    }
}

// --- Euler Rusanov / Lax-Friedrichs (nEquations >= 3) ---
//
// Input:  conserved variables Q = [rho, rho*u, rho*v, rho*w, rho*E]
//         at each face, equation-major layout.
// Output: numerical flux F_num = 0.5*(F_L + F_R) - 0.5*|lambda_max|*(Q_R - Q_L)
//
// Physical flux F(Q):
//   F[0] = rho*u_n
//   F[1] = rho*u_n*u + p*n_x
//   F[2] = rho*u_n*v + p*n_y
//   F[3] = rho*u_n*w + p*n_z
//   F[4] = rho*u_n*H + p*v_mesh   (H = total enthalpy)
//
// Wave speed: |lambda_max| = |u_n - v_mesh| + c

inline Real Square( Real x ) { return x * x; }

void ComputePhysicalFlux(
    Real density, Real u, Real v, Real w, Real pressure,
    Real totalEnthalpy,
    Real nx, Real ny, Real nz, Real vfn,
    Real * physFlux, int nEq )
{
    const Real normalVelocity = u * nx + v * ny + w * nz;
    const Real massFlux = density * normalVelocity;

    physFlux[ 0 ] = massFlux;
    physFlux[ 1 ] = massFlux * u + pressure * nx;
    physFlux[ 2 ] = massFlux * v + pressure * ny;
    if ( nEq >= 4 )
    {
        physFlux[ 3 ] = massFlux * w + pressure * nz;
        physFlux[ 4 ] = massFlux * totalEnthalpy + pressure * vfn;
    }
    else
    {
        // 1D Euler: 3 equations
        physFlux[ 2 ] = massFlux * totalEnthalpy + pressure * vfn;
    }
}

Real ComputeMaxWaveSpeed(
    Real density, Real u, Real v, Real w, Real pressure,
    Real nx, Real ny, Real nz, Real vfn, Real gamma )
{
    const Real normalVelocity = u * nx + v * ny + w * nz;
    const Real soundSpeed = std::sqrt( std::abs( gamma * pressure / density ) );
    return std::abs( normalVelocity - vfn ) + soundSpeed;
}

void CalcEulerRusanovFlux( const FaceStateView & state, FaceFluxView & flux )
{
    const int nEq = state.nEquations;
    const int nFaces = state.nFaces;
    const Real gamma = state.gamma;
    const Real gamm1 = gamma - 1.0;

    for ( int face = 0; face < nFaces; ++ face )
    {
        // Extract left conserved variables
        const Real rhoL  = state.qLeft[ 0 * nFaces + face ];
        const Real rhouL = state.qLeft[ 1 * nFaces + face ];
        const Real rhoEL = state.qLeft[ ( nEq - 1 ) * nFaces + face ];

        Real uL, vL = 0.0, wL = 0.0;
        uL = rhouL / rhoL;
        Real keL = 0.5 * Square( uL );
        if ( nEq >= 4 )
        {
            const Real rhovL = state.qLeft[ 2 * nFaces + face ];
            const Real rhowL = state.qLeft[ 3 * nFaces + face ];
            vL = rhovL / rhoL;
            wL = rhowL / rhoL;
            keL = 0.5 * ( Square( uL ) + Square( vL ) + Square( wL ) );
        }
        const Real pL = gamm1 * ( rhoEL - keL * rhoL );
        const Real hL = ( rhoEL + pL ) / rhoL;

        // Extract right conserved variables
        const Real rhoR  = state.qRight[ 0 * nFaces + face ];
        const Real rhouR = state.qRight[ 1 * nFaces + face ];
        const Real rhoER = state.qRight[ ( nEq - 1 ) * nFaces + face ];

        Real uR, vR = 0.0, wR = 0.0;
        uR = rhouR / rhoR;
        Real keR = 0.5 * Square( uR );
        if ( nEq >= 4 )
        {
            const Real rhovR = state.qRight[ 2 * nFaces + face ];
            const Real rhowR = state.qRight[ 3 * nFaces + face ];
            vR = rhovR / rhoR;
            wR = rhowR / rhoR;
            keR = 0.5 * ( Square( uR ) + Square( vR ) + Square( wR ) );
        }
        const Real pR = gamm1 * ( rhoER - keR * rhoR );
        const Real hR = ( rhoER + pR ) / rhoR;

        // Face normal and mesh velocity
        const Real nx = state.xNormal ? state.xNormal[ face ] : 1.0;
        const Real ny = state.yNormal ? state.yNormal[ face ] : 0.0;
        const Real nz = state.zNormal ? state.zNormal[ face ] : 0.0;
        const Real vfn = state.meshVelocityNormal
            ? state.meshVelocityNormal[ face ] : 0.0;

        // Physical fluxes
        Real fL[ 5 ], fR[ 5 ];
        ComputePhysicalFlux( rhoL, uL, vL, wL, pL, hL,
            nx, ny, nz, vfn, fL, nEq );
        ComputePhysicalFlux( rhoR, uR, vR, wR, pR, hR,
            nx, ny, nz, vfn, fR, nEq );

        // Maximum wave speed
        const Real lambdaL = ComputeMaxWaveSpeed(
            rhoL, uL, vL, wL, pL, nx, ny, nz, vfn, gamma );
        const Real lambdaR = ComputeMaxWaveSpeed(
            rhoR, uR, vR, wR, pR, nx, ny, nz, vfn, gamma );
        const Real lambdaMax = lambdaL > lambdaR ? lambdaL : lambdaR;

        // Face area
        const Real area = state.faceArea ? state.faceArea[ face ] : 1.0;

        // Rusanov flux: F = 0.5*(F_L + F_R) - 0.5*|lambda|*(Q_R - Q_L)
        for ( int eq = 0; eq < nEq; ++ eq )
        {
            const int idx = eq * nFaces + face;
            flux.values[ idx ] = area * (
                0.5 * ( fL[ eq ] + fR[ eq ] )
                - 0.5 * lambdaMax * ( state.qRight[ idx ] - state.qLeft[ idx ] ) );
        }
    }
}

void AddEulerFaceFlux(
    const FaceFluxView & flux,
    const FaceConnectivityView & connectivity,
    ResidualView & residual )
{
    const int nEq = flux.nEquations;
    const int nFaces = flux.nFaces;

    for ( int face = 0; face < nFaces; ++ face )
    {
        const int left = connectivity.leftCell[ face ];
        for ( int eq = 0; eq < nEq; ++ eq )
        {
            const Real value = flux.values[ eq * nFaces + face ];
            residual.values[ eq * residual.nCells + left ] -= value;
            if ( ! IsBoundaryFace( connectivity, face ) )
            {
                const int right = connectivity.rightCell[ face ];
                residual.values[ eq * residual.nCells + right ] += value;
            }
        }
    }
}

} // namespace

void CpuFluxBackend::CalcInvFlux(
    const FaceStateView & state,
    FaceFluxView & flux,
    int )
{
    if ( state.nFaces != flux.nFaces || state.nEquations != flux.nEquations )
    {
        throw std::invalid_argument( "CPU flux view dimensions are inconsistent." );
    }
    ValidateFaceStateView( state );

    if ( state.nEquations == 1 )
    {
        CalcScalarInvFlux( state, flux );
    }
    else if ( state.nEquations >= 3 )
    {
        CalcEulerRusanovFlux( state, flux );
    }
    else
    {
        throw std::invalid_argument(
            "CPU flux backend: unsupported number of equations." );
    }
}

void CpuFluxBackend::AddFaceFlux(
    const FaceFluxView & flux,
    const FaceConnectivityView & connectivity,
    ResidualView & residual )
{
    if ( flux.nFaces != connectivity.nFaces
         || flux.nEquations != residual.nEquations )
    {
        throw std::invalid_argument( "CPU residual view dimensions are inconsistent." );
    }
    ValidateFaceConnectivityView( connectivity );

    if ( flux.nEquations == 1 )
    {
        AddScalarFaceFlux( flux, connectivity, residual );
    }
    else
    {
        AddEulerFaceFlux( flux, connectivity, residual );
    }
}

}
