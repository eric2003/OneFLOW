#include "TurbSstPolicy.h"
#include "HXMath.h"
#include <cmath>

BeginNameSpace( ONEFLOW )

Real SstProdkBoussinesq(
    Real eddyViscosity, Real strainRateSquared, Real divergence,
    Real invReynolds, Real density, Real kineticEnergy )
{
    return eddyViscosity * ( strainRateSquared - two3rd * SQR( divergence ) ) * invReynolds
         - two3rd * density * kineticEnergy * divergence;
}

Real SstProdkStrainOnly(
    Real eddyViscosity, Real strainRateSquared, Real invReynolds )
{
    return eddyViscosity * strainRateSquared * invReynolds;
}

Real SstProdkVorticity(
    Real eddyViscosity, Real vorticitySquared,
    Real invReynolds, Real density, Real kineticEnergy, Real divergence )
{
    return eddyViscosity * vorticitySquared * invReynolds
         - two3rd * density * kineticEnergy * divergence;
}

Real SstDisskStandard(
    Real fBetas, Real betas, Real density, Real kineticEnergy, Real omega )
{
    return fBetas * betas * density * kineticEnergy * omega;
}

Real SstDisskDes(
    Real fBetas, Real density, Real kineticEnergy, Real lengthScale )
{
    return fBetas * density * kineticEnergy * std::sqrt( kineticEnergy )
         / lengthScale;
}

Real SstLimitProdk( Real production, Real pklim, Real dissipation )
{
    return MIN( production, pklim * dissipation );
}

Real SstBlend( Real f1, Real value1, Real value2 )
{
    return f1 * value1 + ( one - f1 ) * value2;
}

Real SstProdwMenter(
    Real alphaw, Real density, Real eddyViscosity,
    Real prodk, Real reynolds )
{
    return alphaw * density / ( eddyViscosity + SMALL ) * prodk * reynolds;
}

Real SstProdwWilcox(
    Real alphaw, Real omega, Real kineticEnergy, Real prodk )
{
    return alphaw * omega / ( kineticEnergy + SMALL ) * prodk;
}

Real SstDissw( Real fBeta, Real beta, Real density, Real omega )
{
    return fBeta * beta * density * SQR( omega );
}

Real SstCrossDiffMenter(
    Real f1, Real density, Real sigw2, Real crossTerm, Real omega )
{
    return two * ( one - f1 ) * density * sigw2 * crossTerm / ( omega + SMALL );
}

Real SstCrossDiffWilcox(
    Real sigd, Real density, Real crossTerm, Real omega )
{
    return sigd * density * MAX( crossTerm, zero ) / ( omega + SMALL );
}

Real SstCrossDiffScalar(
    Real density, Real sigw2, Real crossTerm, Real omega )
{
    return two * density * sigw2 * crossTerm / ( omega + SMALL );
}

Real SstCdkwWithFloor( Real crossDiff, Real floorValue )
{
    return MAX( crossDiff, floorValue );
}

Real SstBlendingArg1(
    Real termSqrtKe, Real termViscous, Real termCrossDiff )
{
    return MIN( MAX( termSqrtKe, termViscous ), termCrossDiff );
}

Real SstF1FromArg1( Real arg1 )
{
    return std::tanh( POWER4( arg1 ) );
}

Real SstGradientDotProduct(
    Real dkdx, Real dkdy, Real dkdz,
    Real dwdx, Real dwdy, Real dwdz )
{
    return dkdx * dwdx + dkdy * dwdy + dkdz * dwdz;
}


Real WilcoxXkFromCrossTerm( Real crossTerm, Real omega )
{
    return crossTerm / ( std::pow( omega, 3 ) + SMALL );
}

Real WilcoxFkFromXk( Real xk )
{
    if ( xk <= 0.0 )
    {
        return 1.0;
    }
    const Real xk2 = SQR( xk );
    return ( one + 680.0 * xk2 ) / ( one + 400.0 * xk2 );
}

Real WilcoxFw1998( Real xw )
{
    return ( one + 70.0 * xw ) / ( one + 80.0 * xw );
}

Real WilcoxFw2006( Real xw )
{
    return ( one + 85.0 * xw ) / ( one + 100.0 * xw );
}

Real WilcoxVortexStretchingInvariant(
    Real w12, Real w13, Real w23,
    Real s11, Real s22, Real s33, Real s12, Real s13, Real s23 )
{
    Real xw = - ( w12 * w12 * ( s11 + s22 ) + w13 * w13 * ( s11 + s33 ) + w23 * w23 * ( s22 + s33 ) )
         + two * ( - w13 * w23 * s12 + w12 * w23 * s13 - w12 * w13 * s23 );
    return ABS( xw );
}

EndNameSpace
