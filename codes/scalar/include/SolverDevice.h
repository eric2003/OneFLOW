/*---------------------------------------------------------------------------*\
    OneFLOW - LargeScale Multiphysics Scientific Simulation Environment
    Copyright (C) 2017-2022 He Xin and the OneFLOW contributors.
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


#pragma once

#include "Configure.h"
#include "HXType.h"

BeginNameSpace( ONEFLOW )

#ifdef ENABLE_CUDA
void SetFaceValueCuda(Real *aface, Real *bcell, int *id, int nFaces, int nTCells);
void MyCalcInvFluxCuda(Real *qf1, Real *qf2, Real *invflux, Real *xfn, Real *yfn, Real *zfn, Real *area, int nFaces);
void MyAddF2CFieldCuda(Real *fField, Real *cField, int *lc, int * rc, int nBFaces, int nFaces, int nTCells);
void MyZoneTimeIntergralCuda(Real *res, Real *vol, Real dt, int nCells);
void MyZoneUpdateCuda(Real *q, Real *res, int nCells);

void addWithCuda(int *a, int *b, int *c, unsigned int nElems);
void addRealWithCuda(Real *a, Real *b, Real *c, unsigned int nElems);
void addRealSwapWithCuda(Real *a, Real *b, int * id, Real *c, unsigned int nElems);
void setRealSwapWithCuda(Real *a, int * id, Real *c, unsigned int nElems);
void setRealSwapWithCudaNew(Real *a, Real *b, int * id,  unsigned int nElems);
void setRealSwapWithCudaNewRealProblem(Real *a, Real *b, int * id, unsigned int nFaces, unsigned int nCells);
#endif

EndNameSpace
