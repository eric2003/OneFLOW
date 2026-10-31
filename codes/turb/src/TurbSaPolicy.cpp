#include "TurbSaPolicy.h"
#include "HXMath.h"
#include <cmath>

BeginNameSpace( ONEFLOW )

Real SaChi( Real nuTilde, Real invLaminarKinematicViscosity )
{
    return nuTilde * invLaminarKinematicViscosity + SMALL;
}

Real SaFv1( Real chi, Real cv1Cubed )
{
    const Real chi3 = POWER3( chi );
    return chi3 / ( chi3 + cv1Cubed );
}

Real SaFv2( Real chi, Real fv1 )
{
    return one - chi / ( one + chi * fv1 );
}

Real SaSBar( Real rs, Real fv2 )
{
    return rs * fv2;
}

Real SaModifiedStrainRate( Real omega, Real sBar, Real c2, Real c3 )
{
    if ( sBar >= - c2 * omega )
    {
        return omega + sBar;
    }
    const Real term1 = ( c2 * c2 * omega + c3 * sBar );
    const Real term2 = ( c3 - two * c2 ) * omega - sBar;
    return omega + omega * term1 / term2;
}

Real SaRRatio( Real rs, Real modifiedStrain )
{
    return MIN( rs / ( modifiedStrain + SMALL ), ten );
}

Real SaG( Real r, Real cw2 )
{
    const Real r5 = std::pow( r, 5 );
    return r + cw2 * ( r5 * r - r );
}

Real SaFw( Real g, Real cw3ToSix )
{
    const Real g6 = std::pow( g, 6 );
    const Real or6 = one / six;
    return g * std::pow( ( one + cw3ToSix ) / ( g6 + cw3ToSix ), or6 );
}

Real SaFt2( Real chi, Real ct3, Real ct4, bool enabled )
{
    if ( ! enabled )
    {
        return 0.0;
    }
    return ct3 * std::exp( - ct4 * SQR( chi ) );
}

Real SaProduction( Real cb1, Real ft2, Real modifiedStrain, Real nuTilde )
{
    return cb1 * ( one - ft2 ) * modifiedStrain * nuTilde;
}

Real SaDestruction( Real cw1k, Real fw, Real cb1, Real ft2, Real nuTildeRs )
{
    return ( cw1k * fw - cb1 * ft2 ) * nuTildeRs;
}

Real SaDiffusion( Real cb2OverSigma, Real gradNuSquared, Real invReynolds )
{
    return cb2OverSigma * gradNuSquared * invReynolds;
}

Real SaLimitedDiffusion( Real diffusion, Real rprod, Real production )
{
    return MIN( rprod * production, diffusion );
}

Real SaSourceResidual( Real production, Real destruction, Real diffusion )
{
    return production - destruction + diffusion;
}

Real SaNegativePart( Real x )
{
    return - half * ( x - ABS( x ) );
}

Real SaEddyViscosity( Real density, Real nuTilde, Real fv1 )
{
    return density * nuTilde * fv1;
}

Real SaCappedEddyViscosity(
    Real eddyViscosity, Real laminarViscosity, Real maxViscosityRatio )
{
    return MIN( laminarViscosity * maxViscosityRatio, eddyViscosity );
}


Real SaStrainProductionClassical( Real work1, Real work2, Real work3 )
{
    const Real stp = two * work2 + work3 - two3rd * work1 * work1;
    return std::sqrt( ABS( stp ) );
}

Real SaStrainProductionRotationalV4( Real omega, Real strainProductionTensor )
{
    const Real root = std::sqrt( ABS( strainProductionTensor ) );
    return omega - two * MAX( zero, omega - root );
}

Real SaStrainProductionRotationalV5( Real omega, Real absStrainProductionTensor )
{
    const Real stp = absStrainProductionTensor;
    return omega * ( one - MAX( zero,
        ( omega * omega - stp ) / ( omega * omega + stp + SMALL ) ) );
}


Real SaWorkingVariableUpperBound(
    Real maxViscosityRatio, Real laminarViscosity,
    Real density, Real fv1 )
{
    return maxViscosityRatio * laminarViscosity / ( density * fv1 );
}

Real SaNonNegativeWorkingVariable(
    Real workingVariable, Real fallbackValue )
{
    if ( workingVariable < 0.0 )
    {
        return fallbackValue;
    }
    return workingVariable;
}


Real SaUnderRelaxedUpdate( Real qOld, Real deltaQ, Real relaxCoeff )
{
    const Real qCandidate = qOld + deltaQ;
    return ( one - relaxCoeff ) * qOld + relaxCoeff * qCandidate;
}

EndNameSpace
