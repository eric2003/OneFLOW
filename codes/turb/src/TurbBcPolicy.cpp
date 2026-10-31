#include "TurbBcPolicy.h"
#include "HXMath.h"

BeginNameSpace( ONEFLOW )

Real TurbWallSaImage( Real interiorWorkingVariable )
{
    return - interiorWorkingVariable;
}

Real TurbWallOmegaMenter(
    Real laminarViscosity, Real density, Real beta1,
    Real wallDistanceSquared, Real reynolds )
{
    return 60.0 * laminarViscosity
         / ( density * beta1 * wallDistanceSquared * reynolds );
}

Real TurbWallKineticEnergyImage( Real interiorKineticEnergy )
{
    return - interiorKineticEnergy;
}

EndNameSpace
