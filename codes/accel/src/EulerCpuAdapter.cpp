#include "EulerCpuAdapter.h"
#include "CpuFluxBackend.h"

#include <cmath>
#include <stdexcept>
#include <vector>

BeginNameSpace( ONEFLOW )

namespace
{

void CheckPrimitive( const Real * primitive, int nEquations, Real gamma )
{
    if ( primitive == nullptr || ( nEquations != 3 && nEquations != 5 )
         || gamma <= 1.0 )
    {
        throw std::invalid_argument( "invalid Euler primitive state" );
    }
    for ( int i = 0; i < nEquations; ++ i )
    {
        if ( ! std::isfinite( primitive[ i ] ) )
        {
            throw std::invalid_argument( "non-finite Euler primitive state" );
        }
    }
    const int pressureIndex = nEquations == 3 ? 2 : 4;
    if ( primitive[ 0 ] <= 0.0 || primitive[ pressureIndex ] <= 0.0 )
    {
        throw std::invalid_argument( "non-positive Euler primitive state" );
    }
}

}

void PrimitiveToConserved(
    const Real * primitive, int nEquations, Real gamma, Real * conserved )
{
    CheckPrimitive( primitive, nEquations, gamma );
    if ( conserved == nullptr )
    {
        throw std::invalid_argument( "null Euler conserved output" );
    }

    const Real density = primitive[ 0 ];
    const Real u = primitive[ 1 ];
    const Real v = nEquations == 5 ? primitive[ 2 ] : 0.0;
    const Real w = nEquations == 5 ? primitive[ 3 ] : 0.0;
    const Real pressure = primitive[ nEquations == 3 ? 2 : 4 ];
    const Real kineticEnergy = 0.5 * density * ( u * u + v * v + w * w );

    conserved[ 0 ] = density;
    conserved[ 1 ] = density * u;
    if ( nEquations == 3 )
    {
        conserved[ 2 ] = pressure / ( gamma - 1.0 ) + kineticEnergy;
    }
    else
    {
        conserved[ 2 ] = density * v;
        conserved[ 3 ] = density * w;
        conserved[ 4 ] = pressure / ( gamma - 1.0 ) + kineticEnergy;
    }
}

void PackPrimitiveFaceStates(
    const Real * primitiveLeft,
    const Real * primitiveRight,
    int nFaces,
    int nEquations,
    Real gamma,
    Real * conservedLeft,
    Real * conservedRight )
{
    if ( nFaces <= 0 || primitiveLeft == nullptr || primitiveRight == nullptr
         || conservedLeft == nullptr || conservedRight == nullptr )
    {
        throw std::invalid_argument( "invalid Euler face state pack request" );
    }
    Real primitive[ 5 ] = {};
    Real conserved[ 5 ] = {};
    for ( int face = 0; face < nFaces; ++ face )
    {
        for ( int equation = 0; equation < nEquations; ++ equation )
        {
            primitive[ equation ] = primitiveLeft[ equation * nFaces + face ];
        }
        PrimitiveToConserved( primitive, nEquations, gamma, conserved );
        for ( int equation = 0; equation < nEquations; ++ equation )
        {
            conservedLeft[ equation * nFaces + face ] = conserved[ equation ];
        }

        for ( int equation = 0; equation < nEquations; ++ equation )
        {
            primitive[ equation ] = primitiveRight[ equation * nFaces + face ];
        }
        PrimitiveToConserved( primitive, nEquations, gamma, conserved );
        for ( int equation = 0; equation < nEquations; ++ equation )
        {
            conservedRight[ equation * nFaces + face ] = conserved[ equation ];
        }
    }
}

void EulerCpuAdapter::CalcInvFlux(
    const PrimitiveFaceStateView & primitiveState,
    FaceFluxView & flux,
    int scheme ) const
{
    if ( primitiveState.nFaces <= 0
         || primitiveState.nFaces != flux.nFaces
         || primitiveState.nEquations != flux.nEquations
         || primitiveState.faceArea == nullptr
         || flux.values == nullptr )
    {
        throw std::invalid_argument( "invalid Euler CPU adapter request" );
    }

    std::vector< Real > conservedLeft(
        primitiveState.nFaces * primitiveState.nEquations );
    std::vector< Real > conservedRight(
        primitiveState.nFaces * primitiveState.nEquations );
    PackPrimitiveFaceStates(
        primitiveState.primitiveLeft, primitiveState.primitiveRight,
        primitiveState.nFaces, primitiveState.nEquations, primitiveState.gamma,
        conservedLeft.data(), conservedRight.data() );

    FaceStateView conservedState;
    conservedState.nFaces = primitiveState.nFaces;
    conservedState.nEquations = primitiveState.nEquations;
    conservedState.qLeft = conservedLeft.data();
    conservedState.qRight = conservedRight.data();
    conservedState.xNormal = primitiveState.xNormal;
    conservedState.yNormal = primitiveState.yNormal;
    conservedState.zNormal = primitiveState.zNormal;
    conservedState.meshVelocityNormal = primitiveState.meshVelocityNormal;
    conservedState.faceArea = primitiveState.faceArea;
    conservedState.gamma = primitiveState.gamma;

    CpuFluxBackend backend;
    backend.CalcInvFlux( conservedState, flux, scheme );
}

void EulerCpuAdapter::AddFaceFlux(
    const FaceFluxView & flux,
    const FaceConnectivityView & connectivity,
    ResidualView & residual ) const
{
    CpuFluxBackend backend;
    backend.AddFaceFlux( flux, connectivity, residual );
}

EndNameSpace
