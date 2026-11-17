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
#include "HXVector.h"
#include <memory>
#include <vector>
BeginNameSpace( ONEFLOW )

class FieldWrap;
class BasicBgField
{
public:
    BasicBgField() = default;
    ~BasicBgField();  // out-of-line: Free()
public:
    using Field3DType = HXVector< HXVector< HXVector< std::unique_ptr<FieldWrap> > > >;
public:
    Field3DType data;
public:
    void Init();
    void Free();
};

// Multi-solver / multi-level field slots (logistics for algorithm plugins).
// Layout: data[zone] -> BasicBgField::data[solverIndex][fid][gridLevel]
//   fid: 0 = FLOW  (solver working buffer; not the same object as named grid field "q")
//   fid: 1 = RHS   (solver residual buffer accessed via FieldHome::GetUnsField(fieldId))
// Each slot is FieldHome::CreateField() -> unique_ptr<FieldWrap> owning MRField.
// GetFieldWrap returns a non-owning raw view (.get()) for solvers to read/write in place.
//
// Relation to named residual fields on UnsGrid:
//   SolverInfo::residualName is the zone-local database field name:
//     NS / INS -> "res"
//     Turbulence -> "turbres"
//   Never hard-code "res"/"turbres" at call sites; use FieldHome::GetResidualField /
//   GetConstResidualField / GetResidualName (driven by SolverState::solverType).
//   LoadResiduals: NegField(named residual, BgField RHS slot) - fid=1 -> residualName field.
//   ResidualTask aggregates the named residual (FieldHome::GetConstResidualField).
//   Do not confuse BgField fid=1 storage with the string keys "res" / "turbres" on the grid DataBase.
class BgField
{
public:
    BgField() = default;
    ~BgField() = default;
protected:
    static HXVector< std::unique_ptr< BasicBgField > > data;
public:
    static bool flag;
    static void Init();
    static void Free();
    // Non-owning: lifetime owned by BgField::data[zid].
    static FieldWrap * GetFieldWrap( int zid, int solverIndex, int fid, int gl );
    // Read-only view for diagnostics / pure-read algorithms.
    static const FieldWrap * GetConstFieldWrap( int zid, int solverIndex, int fid, int gl );
};


EndNameSpace
