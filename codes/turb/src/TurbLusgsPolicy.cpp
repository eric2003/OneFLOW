#include "TurbLusgsPolicy.h"
#include "HXMath.h"

BeginNameSpace( ONEFLOW )

Real TurbLusgsMatrixTimesDq( Real matrixDiag, Real deltaQ )
{
    return matrixDiag * deltaQ;
}

Real TurbLusgsSolveByRadius( Real rightHandSide, Real spectralRadius )
{
    return rightHandSide / spectralRadius;
}

Real TurbLusgsDifference( Real a, Real b )
{
    return a - b;
}

Real TurbLusgsSweepNormContribution( Real dqNew, Real dqOld )
{
    return SQR( dqNew - dqOld );
}

EndNameSpace
