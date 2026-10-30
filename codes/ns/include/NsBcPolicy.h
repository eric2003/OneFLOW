#pragma once

#include "HXDefine.h"

BeginNameSpace( ONEFLOW )

// @doc-algorithm boundary-velocity
// @doc-stage symmetry-reflect
// @doc-pure
// @doc-formula v' = v - 2 n (n·v - v_fn)
void ReflectVelocityRelativeToMovingFace(
    Real u, Real v, Real w,
    Real xfn, Real yfn, Real zfn, Real vfn,
    Real & uReflected, Real & vReflected, Real & wReflected );

// @doc-algorithm boundary-velocity
// @doc-stage wall-image
// @doc-pure
// @doc-formula v' = -v + 2 v_wall
void ImageVelocityAboutWall(
    Real u, Real v, Real w,
    Real uWall, Real vWall, Real wWall,
    Real & uImage, Real & vImage, Real & wImage );

// @doc-algorithm boundary-velocity
// @doc-stage face-average
// @doc-pure
Real ArithmeticAverage( Real a, Real b );

EndNameSpace
