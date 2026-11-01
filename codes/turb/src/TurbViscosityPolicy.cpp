#include "TurbViscosityPolicy.h"
#include "HXMath.h"

BeginNameSpace( ONEFLOW )

Real TurbCappedEddyViscosity(
    Real eddyViscosity, Real laminarViscosity, Real maxViscosityRatio )
{
    return MIN( laminarViscosity * maxViscosityRatio, eddyViscosity );
}

bool TurbDensityIsInvalid( Real density )
{
    return density < 0.0 || NotANumber( density );
}

EndNameSpace
