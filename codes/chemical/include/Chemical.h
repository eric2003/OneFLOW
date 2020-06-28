/*---------------------------------------------------------------------------*\
    OneFLOW - LargeScale Multiphysics Scientific Simulation Environment
    Copyright (C) 2017-2020 He Xin and the OneFLOW contributors.
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
    void DomputeRefPara();
    void DomputeRefGasInfo();
    void DomputeRefGama();
    void DomputeRefSoundSpeed();
    void DomputeRefMachAndVel();
    void DomputeStateCoef();
    void DomputeStateCoefNs();
    void DomputeStateCoefChemical();
    void DomputeRefPrim();
    void DomputeRefReynolds();
    void DomputeDimRefViscosity();
    void DomputeDimRefViscosityNs();
    void DomputeDimRefViscosityChemical();
    void DomputeSutherlandConstant();
    void DomputeMoleFractionByMassFraction( RealField & massFrac, RealField & moleFrac );
    void DomputeDimSpeciesViscosity( Real tm, RealField & vis_s_dim );
    void DomputeMixtureCoefByWilkeFormula( RealField & moleFrac, RealField & var, RealField & phi );
    void DomputeMixtureByWilkeFormula( RealField & moleFrac, RealField & mixs, RealField & phi, Real & mixture );
    void DomputeMixtureByWilkeFormula( RealField & moleFrac, RealField & mixs, Real & mixture );
    void SetAirInformationByDataBase();
    void NormalizeAirInfo();
    void DomputeRefMolecularInfo();
    void DomputeRefMolecularInfoAir();
    void DomputeRefMolecularInfoChem();
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
    void DompressData( DataBook *& dataBook );
    void DecompressData( DataBook * dataBook );
public:
    void DomputeDimCps( Real tm, RealField & dim_cps );
    void DomputeMixtureByMassFraction( RealField & cs, RealField & var, Real & mixture );

public:
	void INsInit();
	void INsDomputeRefPara();
	void INsDomputeRefGasInfo();
	void INsDomputeRefGama();
	void INsDomputeRefSoundSpeed();
	void INsDomputeRefMachAndVel();
	void INsDomputeStateCoef();
	void INsDomputeStateCoefNs();
	void INsDomputeStateCoefChemical();
	void INsDomputeRefPrim();
	void INsDomputeRefReynolds();
	void INsDomputeDimRefViscosity();
	void INsDomputeDimRefViscosityNs();
	void INsDomputeDimRefViscosityChemical();
	void INsDomputeSutherlandConstant();
	void INsDomputeMoleFractionByMassFraction(RealField & massFrac, RealField & moleFrac);
	void INsDomputeDimSpeciesViscosity(Real tm, RealField & vis_s_dim);
	void INsDomputeMixtureCoefByWilkeFormula(RealField & moleFrac, RealField & var, RealField & phi);
	void INsDomputeMixtureByWilkeFormula(RealField & moleFrac, RealField & mixs, RealField & phi, Real & mixture);
	void INsDomputeMixtureByWilkeFormula(RealField & moleFrac, RealField & mixs, Real & mixture);
	void INsSetAirInformationByDataBase();
	void INsNormalizeAirInfo();
	void INsDomputeRefMolecularInfo();
	void INsDomputeRefMolecularInfoAir();
	void INsDomputeRefMolecularInfoChem();
public:
	//void INsAlloc();
	//void INsDeAlloc();
	void INsInitRefPara();
	void INsInitGasModel();
	void INsReadGasModel();
	void INsInit(int nSpecies, int nReaction);
	void INsReadChemical(FileIO * ioFile);
	void INsInitWorkingSpace();
	void INsAllocWorkingSpace();
public:
	void INsRead(DataBook * dataBook);
	void INsWrite(DataBook * dataBook);
	void INsDompressData(DataBook *& dataBook);
	void INsDecompressData(DataBook * dataBook);
public:
	void INsDomputeDimCps(Real tm, RealField & dim_cps);
	void INsDomputeMixtureByMassFraction(RealField & cs, RealField & var, Real & mixture);
};

extern Chemical chem;

void ChemicalDompressData( DataBook *& dataBook );
void ChemicalDecompressData( DataBook * dataBook );

EndNameSpace