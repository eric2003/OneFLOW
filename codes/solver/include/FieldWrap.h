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
#include "HXDefine.h"
#include "UnsteadyFieldNames.h"
#include "HXArray.h"
#include <memory>

BeginNameSpace( ONEFLOW )

const int FLOW_RHS_ORDER = 0;
const int RHS_FLOW_ORDER = 1;

class FieldWrap
{
public:
    FieldWrap() = default;
    ~FieldWrap() = default;
protected:
    // When owning: owned holds the field and view == owned.get().
    // When non-owning: owned is null and view points at an external field.
    std::unique_ptr<MRField> owned;
    MRField * view = nullptr;
public:
    MRField * GetUnsField();
    const MRField * GetConstUnsField() const;
    void SetUnsField( MRField * unsField, bool deleteFlag = false );
    void SetOwnedField( std::unique_ptr<MRField> field );
};

class Grid;

class FieldHome
{
public:
    FieldHome() = default;
    ~FieldHome() = default;
public:
    // Factory: allocate owned MRField + FieldWrap (solver plugs in via solverType/level).
    static std::unique_ptr<FieldWrap> CreateField();
    static std::unique_ptr<FieldWrap> CreateField( int solverType );
    static std::unique_ptr<FieldWrap> CreateField( int solverType, int level );
    // Non-owning view of an existing named field on the current grid (wrapper still unique_ptr).
    static std::unique_ptr<FieldWrap> GetFieldWrap( const std::string & fieldName );
public:
    static void SetField( const std::string & fieldName, Real value );
    static void SetField( int fieldId, const std::string & fieldName, int orderFlag );
    static void SetUnsField( int fieldId, const std::string & fieldName, int orderFlag );
public:
    static FieldWrap * GetFieldWrap( int fieldId );
    static const FieldWrap * GetConstFieldWrap( int fieldId );
    static MRField * GetUnsField( int fieldId );
    static const MRField * GetConstUnsField( int fieldId );
    static MRField * GetUnsField( Grid * grid, const std::string & fieldName );
    static const MRField * GetConstUnsField( const Grid * grid, const std::string & fieldName );
    // No-arg overloads use SolverState::solverType (current solver only).
    // solverType overloads are for cross-solver binders (e.g. UNsField::Init under TURB).

    static const std::string & GetResidualName();
    static const std::string & GetResidualName( int solverType );
    static MRField * GetResidualField( Grid * grid );
    static MRField * GetResidualField( Grid * grid, int solverType );
    static const MRField * GetConstResidualField( const Grid * grid );
    static const MRField * GetConstResidualField( const Grid * grid, int solverType );

    static const std::string & GetSolutionName();
    static const std::string & GetSolutionName( int solverType );
    static MRField * GetSolutionField( Grid * grid );
    static MRField * GetSolutionField( Grid * grid, int solverType );
    static const MRField * GetConstSolutionField( const Grid * grid );

    static const std::string & GetGradName( int index );
    static const std::string & GetGradName( int solverType, int index );
    static MRField * GetGradField( Grid * grid, int index );
    static MRField * GetGradField( Grid * grid, int solverType, int index );
    static const MRField * GetConstGradField( const Grid * grid, int index );

    static const std::string & GetUpdateDeltaName();
    static const std::string & GetUpdateDeltaName( int solverType );
    static MRField * GetUpdateDeltaField( Grid * grid );
    static MRField * GetUpdateDeltaField( Grid * grid, int solverType );

    static const std::string & GetBcGradName( int index );
    static const std::string & GetBcGradName( int solverType, int index );
    static MRField * GetBcGradField( Grid * grid, int index );
    static MRField * GetBcGradField( Grid * grid, int solverType, int index );

    // Dual-time stage solution: stageSolutionNames[0]=q1, [1]=q2 (or turbq1/2).
    static const std::string & GetStageSolutionName( int index );
    static const std::string & GetStageSolutionName( int solverType, int index );
    static MRField * GetStageSolutionField( Grid * grid, int index );
    static MRField * GetStageSolutionField( Grid * grid, int solverType, int index );

    // Dual-time stage residual: stageResidualNames[0]=res1, [1]=res2 (may be empty).
    static const std::string & GetStageResidualName( int index );
    static const std::string & GetStageResidualName( int solverType, int index );
    static MRField * GetStageResidualField( Grid * grid, int index );
    static MRField * GetStageResidualField( Grid * grid, int solverType, int index );

    // Build Current/Previous/Old flow+residual name lists from SolverInfo
    // (aligned with system/*/alloc/unsteady.txt).
    static UnsteadyFieldNames BuildUnsteadyFieldNames( int solverType );

    static void GetSourceTargetField( Grid * grid, int fieldId, const std::string & fieldName, MRField *& sField, MRField *& tField, int orderFlag );
};

// ---------------------------------------------------------------------------
// Field naming convention (prefix strategy)
// ---------------------------------------------------------------------------
// Driven by SolverInfo for the *current* solver (SolverState::solverType):
//
//   Role              NS / INS              Turbulence
//   ---------------   -------------------   -------------------
//   solution          q                     turbq
//   residual          res                   turbres
//   grad (d*/dx,y,z)  dqdx, dqdy, dqdz      turbdqdx, turbdqdy, turbdqdz
//   update delta      dq                    turbdq
//   BC face grad      bcdqdx/y/z            turbbcdqdx/y/z
//
// Role separation (multi-solver):
//   - "Current solver logistics"  -> Get*(grid) or Get*(grid, SolverState::solverType)
//   - "Foreign solver read"       -> Get*(grid, foreignSolverType)  e.g. turb reads NS
//   - "Dual-use binder"           -> must pin identity (UNsField::Init -> NS_SOLVER)
//   - "Shared physics fields"     -> literal names (visl, vist, timestep) until owned
//   Prefer typed views (UNsField / UTurbField) over ad-hoc GetFieldPointer in kernels.
//
// Rules:
//   1. Current-solver: FieldHome::Get*(grid) via SolverState.
//   2. Cross-solver / dual-use: FieldHome::Get*(grid, solverType).
//      UNsField::Init is also called from UTurbInvFlux under TURB - must use
//      NS_SOLVER overloads so unsf.* stay on mean-flow fields.
//   3. NsLimField (nslimiter under TURB): literal NS names only in Init/BcQlQrFix.
//   4. Stage q1/q2/res1/res2: FieldHome::GetStageSolution* / GetStageResidual*.
//   5. limiter / bc_q / visl / vist / timestep / tempr: literal names (shared or LimField).
// ---------------------------------------------------------------------------

void SetField( FieldWrap * fieldWrap, Real value );
void SetField( MRField * field, Real value );
void SetField( RealField & field, Real value );
void SetField( MRField * field1, MRField * field2 );
void NegField( MRField * field1, MRField * field2 );


EndNameSpace
