/*---------------------------------------------------------------------------*\
OneFLOW - LargeScale Multiphysics Scientific Simulation Environment
Copyright (C) 2017-2019 He Xin and the OneFLOW contributors.
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

BeginNameSpace( ONEFLOW )

class MolecularProperty;
class ReactionRate;
class Stoichiometric;
class BlotterCurve;
class Thermodynamic;

class FileIO;
class DataBook;

//通用气体常数，单位为J/( Mol * K )
const double rjmk = 8.31434;

class Chemical
{
public:
    Chemical();
    ~Chemical();
public:
    MolecularProperty * moleProp;
    ReactionRate * reactionRate;
    Stoichiometric * stoichiometric;
    BlotterCurve * blotterCurve;
    Thermodynamic * thermodynamic;
    int nSpecies, nReaction;
public:
    //working variables
    RealField xi_s;
    RealField cs_s;
    RealField vis_s;
    RealField vis_phi_s;
    RealField cp_s;
    RealField dim_cp_s;
    RealField hint_s;
    RealField work_s;
public:
    void Init();
    void ComputeRefPara();
    void ComputeRefGasInfo();
    void ComputeRefGama();
    void ComputeRefSoundSpeed();
    void ComputeRefMachAndVel();
    void ComputeStateCoef();
    void ComputeStateCoefNs();
    void ComputeStateCoefChemical();
    void ComputeRefPrim();
    void ComputeRefReynolds();
    void ComputeDimRefViscosity();
    void ComputeDimRefViscosityNs();
    void ComputeDimRefViscosityChemical();
    void ComputeSutherlandConstant();
    void ComputeMoleFractionByMassFraction( RealField & massFrac, RealField & moleFrac );
    void ComputeDimSpeciesViscosity( Real tm, RealField & vis_s_dim );
    void ComputeMixtureCoefByWilkeFormula( RealField & moleFrac, RealField & var, RealField & phi );
    void ComputeMixtureByWilkeFormula( RealField & moleFrac, RealField & mixs, RealField & phi, Real & mixture );
    void ComputeMixtureByWilkeFormula( RealField & moleFrac, RealField & mixs, Real & mixture );
    void SetAirInformationByDataBase();
    void NormalizeAirInfo();
    void ComputeRefMolecularInfo();
    void ComputeRefMolecularInfoAir();
    void ComputeRefMolecularInfoChem();
public:
    void Alloc();
    void DeAlloc();
    void InitRefPara();
    void InitGasModel();
    void ReadGasModel();
    void Init( int nSpecies, int nReaction );
    void ReadChemical( FileIO * ioFile );
    void InitWorkingSpace();
    void AllocWorkingSpace();
public:
    void Read( DataBook * dataBook );
    void Write( DataBook * dataBook );
    void CompressData( DataBook *& dataBook );
    void DecompressData( DataBook * dataBook );
public:
    void ComputeDimCps( Real tm, RealField & dim_cps );
    void ComputeMixtureByMassFraction( RealField & cs, RealField & var, Real & mixture );
};

extern Chemical chem;

void ChemicalCompressData( DataBook *& dataBook );
void ChemicalDecompressData( DataBook * dataBook );

EndNameSpace