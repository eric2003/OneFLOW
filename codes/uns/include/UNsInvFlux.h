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
#include "NsInvFlux.h"
#include "AccelViews.h"
#include <vector>
#include <memory>

BeginNameSpace( ONEFLOW )

class UNsFField;
class Limiter;
class LimField;
class FluxBackend;
class UnsGrid;
class BcRecord;

class UNsInvFlux : public NsInvFlux
{
public:
    UNsInvFlux ();
    ~UNsInvFlux();
public:
    void CalcFlux();
    void Alloc();
    void DeAlloc();
    void CalcInvFlux();
    void CalcInvFluxCpuBatch();
    void CalcInvFluxHipBatch();
    void CalcAndAddInvFluxHipBatch();
    void CalcInvFluxBatch( FluxBackend & backend );
    bool UseCpuBatchAdapter() const;
    bool UseHipBatchAdapter() const;
    bool UseHipGradient() const;
    bool UseHipDeviceReconstruction() const;
    bool UseHipDeviceStateUpdate() const;
    void CalcInvFace();
    void CalcLimiter();
    void AddInvFlux();
    void PrepareFaceValue();
    void UpdateFaceInvFlux();
    void DumpInvFluxTrace();
    void DumpInvFluxStageTrace();
public:
    void GetQlQrField();
    void ReconstructFaceValueField();
    void BoundaryQlQrFixField();
public:
    std::unique_ptr<Limiter> limiter;
    LimField * limf;
    std::unique_ptr<MRField> invflux;
    // Boundary reconstruction metadata is topology-owned and remains valid
    // across RK stages.  Keep the host containers stable so the HIP backend
    // can detect and skip redundant operation metadata uploads.
    std::vector< ReconstructionBoundaryOperation > hipBoundaryOperation;
    std::vector< BoundaryGhostUpdatePolicy > hipBoundaryGhostPolicy;
    std::vector< Real > hipBoundaryQ;
    const UnsGrid * hipBoundaryGrid = nullptr;
    const BcRecord * hipBoundaryRecord = nullptr;
    int hipBoundaryFaces = -1;
    int hipBoundaryBoundaryFaces = -1;
    int hipBoundaryEquations = -1;
    bool hipBoundaryHasSolid = false;
};

EndNameSpace
