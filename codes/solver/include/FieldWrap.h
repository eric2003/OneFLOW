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
    // Named residual on the zone grid (NS/INS: "res", Turb: "turbres" via SolverInfo::residualName).
    static const std::string & GetResidualName();
    static MRField * GetResidualField( Grid * grid );
    static const MRField * GetConstResidualField( const Grid * grid );
    // Primary solution field name: SolverInfo::gradString[0] (NS: "q", Turb: "turbq").
    // Prefer this over hard-coded "q" for current-solver logistics. Do not use for
    // cross-solver coupling (e.g. turb reading NS mean flow still uses literal "q").
    static const std::string & GetSolutionName();
    static MRField * GetSolutionField( Grid * grid );
    static const MRField * GetConstSolutionField( const Grid * grid );
    // Gradient names from SolverInfo::gradString: [0]=solution, [1]=dx, [2]=dy, [3]=dz.
    // index 0 == GetSolutionName(); out-of-range returns empty string / null field.
    static const std::string & GetGradName( int index );
    static MRField * GetGradField( Grid * grid, int index );
    static const MRField * GetConstGradField( const Grid * grid, int index );
    // Update-delta database name: SolverInfo::implicitString[1] (NS: "dq", Turb: "turbdq").
    // Binds the named field only. MULTI_STAGE may later alias the *Com dq pointer to residual
    // storage; that alias is independent of this name lookup.
    static const std::string & GetUpdateDeltaName();
    static MRField * GetUpdateDeltaField( Grid * grid );
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
// Rules:
//   1. Current-solver logistics (UNsCom, UTurbCom, U*Grad, residual tasks):
//      use FieldHome::GetSolutionField / GetResidualField / GetGradField.
//   2. Cross-solver coupling (turb reads NS mean flow, aero force, visualize
//      primitive NS components with IDX::*): keep literal names ("q", "dqdx", ...).
//   3. Derived fields not in SolverInfo (q1/q2, limiter, visl, tempr, ...):
//      stay literal until a dedicated registry exists.
//   4. Update delta: use GetUpdateDeltaName/GetUpdateDeltaField for the database
//      field named in implicitString[1]. MULTI_STAGE still rebinds *Com.dq to the
//      residual pointer afterward; GetUpdateField (LUSGS) also uses implicitString.
//      ULhs scales residual via GetResidualField (not the "dq" database field).
// ---------------------------------------------------------------------------

void SetField( FieldWrap * fieldWrap, Real value );
void SetField( MRField * field, Real value );
void SetField( RealField & field, Real value );
void SetField( MRField * field1, MRField * field2 );
void NegField( MRField * field1, MRField * field2 );


EndNameSpace
