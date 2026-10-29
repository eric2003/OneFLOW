#pragma once

#include "HXDefine.h"

BeginNameSpace( ONEFLOW )

Real BarthFunction( Real dq, Real minv, Real maxv, Real dot );
Real VencatFunction( Real dq, Real minv, Real maxv, Real eps );
Real VencatFunctionEric( Real dq, Real minv, Real maxv, Real eps );
Real Vencat( Real x, Real y, Real eps );
Real VencatEric( Real x, Real y, Real eps );

EndNameSpace
