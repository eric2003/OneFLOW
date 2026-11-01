#include "TurbViscosityPolicy.h"
#include "HXMath.h"

BeginNameSpace( ONEFLOW )

Real TurbCappedEddyViscosity(
    Real eddyViscosity, Real laminarViscosity, Real maxViscosityRatio )
{
    return MIN( laminarViscosity * maxViscosityRatio, eddyViscosity );
}

EndNameSpace
