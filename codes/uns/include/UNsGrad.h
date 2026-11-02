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


#pragma once
#include "Grad.h"

BeginNameSpace( ONEFLOW )

class UNsGrad : public Grad
{
public:
    UNsGrad();
    ~UNsGrad();
public:
    void Init() override;
    // downloadToHost / mirrorBoundaryGradient are capability decisions, not
    // environment guesses.  A device reconstruction reads the resident
    // gradient, so the host mirror is only needed where the host interface
    // exchange touches it: the boundary-adjacent owned cells (interface send
    // values) and the ghost cells that Grad::SwapBcGrad fills from the receive
    // buffer.  The one-round exchange lag of the CPU oracle is preserved
    // because the swap still runs on the host.
    // packedExchange: run the device pack/unpack interface exchange that fills
    //                 the resident ghost gradient (the CPU lag contract is kept
    //                 because the send/recv storages are exchanged by the
    //                 post-process interface round trip).
    // mirrorToHost:   also copy the boundary-adjacent owned gradient into the
    //                 host fields; only aero forces, the optional boundary
    //                 store and field dumps read it.
    void CalcGradHip( bool downloadToHost, bool packedExchange,
                      bool mirrorToHost );
};

class UTGrad : public Grad
{
public:
    UTGrad();
    ~UTGrad();
public:
    void Init() override;
};

extern UNsGrad uns_grad;
extern UTGrad ut_grad;

EndNameSpace
