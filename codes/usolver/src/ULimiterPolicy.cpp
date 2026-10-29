/*---------------------------------------------------------------------------*\
    OneFLOW - LargeScale Multiphysics Scientific Simulation Environment
    Copyright (C) 2017-2026 He Xin and the OneFLOW contributors.
-------------------------------------------------------------------------------
License
    This file is part of OneFLOW.

    OneFLOW is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OneFLOW is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OneFLOW.  If not, see <http://www.gnu.org/licenses/>.
\*---------------------------------------------------------------------------*/

#include "ULimiterPolicy.h"
#include "HXMath.h"

BeginNameSpace( ONEFLOW )

Real BarthFunction( Real dq, Real minv, Real maxv, Real dot )
{
    if ( dot < 0.0 ) 
    {
        return 0.0;
    }
    else if ( dq > maxv )
    {
        return maxv / dq;
    }
    else if ( dq < minv )
    {
        return minv / dq;
    }

    return 1.0;
}

Real Vencat( Real x, Real y, Real eps )
{
    Real x2 = x * x;
    Real xy = x * y;
    Real y2 = y * y;
    Real cc = 2.0;
    //Real cc = 0.1;
    return ( ( x2 + eps + cc * xy ) / ( x2 + cc * y2 + xy + eps ) );
}

Real VencatEric( Real x, Real y, Real eps )
{
    Real x2 = x * x;
    Real xy = x * y;
    Real y2 = y * y;
    return ( ( x2 + eps + 1.0 * xy ) / ( x2 + 2.0 * y2 + xy + eps ) );
}

Real VencatFunction( Real dq, Real minv, Real maxv, Real eps )
{
    if ( dq > SMALL )
    {
        return Vencat( maxv, dq, eps );
    }
    else if ( dq < - SMALL )
    {
        return Vencat( minv, dq, eps );
    }

    return 1.0;
}

Real VencatFunctionEric( Real dq, Real minv, Real maxv, Real eps )
{
    if ( dq > SMALL )
    {
        return VencatEric( maxv, dq, eps );
    }
    else if ( dq < - SMALL )
    {
        return VencatEric( minv, dq, eps );
    }

    return 1.0;
}

EndNameSpace
