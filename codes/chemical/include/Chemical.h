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
    void CalcRefPara();
    void CalcRefGasInfo();
    void CalcRefGama();
    void CalcRefSoundSpeed();
    void CalcRefMachAndVel();
    void CalcStateCoef();
    void CalcStateCoefNs();
    void CalcStateCoefChemical();
    void CalcRefPrim();
    void CalcRefReynolds();
    void CalcDimRefViscosity();
    void CalcDimRefViscosityNs();
    void CalcDimRefViscosityChemical();
    void CalcSutherlandConstant();
    void CalcMoleFractionByMassFraction( RealField & massFrac, RealField & moleFrac );
    void CalcDimSpeciesViscosity( Real tm, RealField & vis_s_dim );
    void CalcMixtureCoefByWilkeFormula( RealField & moleFrac, RealField & var, RealField & phi );
    void CalcMixtureByWilkeFormula( RealField & moleFrac, RealField & mixs, RealField & phi, Real & mixture );
    void CalcMixtureByWilkeFormula( RealField & moleFrac, RealField & mixs, Real & mixture );
    void SetAirInformationByDataBase();
    void NormalizeAirInfo();
    void CalcRefMolecularInfo();
    void CalcRefMolecularInfoAir();
    void CalcRefMolecularInfoChem();
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
    void CalcDimCps( Real tm, RealField & dim_cps );
    void CalcMixtureByMassFraction( RealField & cs, RealField & var, Real & mixture );

public:
	void INsInit();
	void INsCalcRefPara();
	void INsCalcRefGasInfo();
	void INsCalcRefGama();
	void INsCalcRefSoundSpeed();
	void INsCalcRefMachAndVel();
	void INsCalcStateCoef();
	void INsCalcStateCoefNs();
	void INsCalcStateCoefChemical();
	void INsCalcRefPrim();
	void INsCalcRefReynolds();
	void INsCalcDimRefViscosity();
	void INsCalcDimRefViscosityNs();
	void INsCalcDimRefViscosityChemical();
	void INsCalcSutherlandConstant();
	void INsCalcMoleFractionByMassFraction(RealField & massFrac, RealField & moleFrac);
	void INsCalcDimSpeciesViscosity(Real tm, RealField & vis_s_dim);
	void INsCalcMixtureCoefByWilkeFormula(RealField & moleFrac, RealField & var, RealField & phi);
	void INsCalcMixtureByWilkeFormula(RealField & moleFrac, RealField & mixs, RealField & phi, Real & mixture);
	void INsCalcMixtureByWilkeFormula(RealField & moleFrac, RealField & mixs, Real & mixture);
	void INsSetAirInformationByDataBase();
	void INsNormalizeAirInfo();
	void INsCalcRefMolecularInfo();
	void INsCalcRefMolecularInfoAir();
	void INsCalcRefMolecularInfoChem();
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
	void INsCompressData(DataBook *& dataBook);
	void INsDecompressData(DataBook * dataBook);
public:
	void INsCalcDimCps(Real tm, RealField & dim_cps);
	void INsCalcMixtureByMassFraction(RealField & cs, RealField & var, Real & mixture);
};

extern Chemical chem;

void ChemicalCompressData( DataBook *& dataBook );
void ChemicalDecompressData( DataBook * dataBook );

EndNameSpace