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
#include "Chemical.h"
#include "MolecularProperty.h"
#include "ReactionRate.h"
#include "Stoichiometric.h"
#include "BlotterCurve.h"
#include "Thermodynamic.h"
#include "NsCom.h"
#include "Parallel.h"
#include "FileIO.h"
#include "DataBook.h"
#include "Atmosphere.h"
#include "HXMath.h"
#include "NsIdx.h"
#include "Ctrl.h"
#include "INsCom.h"
#include "INsIdx.h"

BeginNameSpace( ONEFLOW )

Chemical chem;

Chemical::Chemical()
{
    Alloc();
}

Chemical::~Chemical()
{
    DeAlloc();
}

void Chemical::Alloc()
{
    moleProp = new MolecularProperty();
    reactionRate = new ReactionRate();
    stoichiometric = new Stoichiometric();
    blotterCurve = new BlotterCurve();
    thermodynamic = new Thermodynamic();
}

void Chemical::DeAlloc()
{
    delete moleProp;
    delete reactionRate;
    delete stoichiometric;
    delete blotterCurve;
    delete thermodynamic;
}

void Chemical::InitGasModel()
{
    if ( nscom.chemModel <= 0 ) return;

    if ( Parallel::pid == Parallel::serverid )
    {
        ReadGasModel();
    }
    ONEFLOW::HXBcast( ChemicalCompressData, ChemicalDecompressData, Parallel::GetServerid() );
}

void Chemical::InitWorkingSpace()
{
    nscom.nBEqu = nscom.nEqu;
    if ( nscom.chemModel == 0 )
    {
        nscom.nSpecies = 0;
    }
    else
    {
        nscom.nSpecies = nSpecies;
    }

    nscom.nTEqu = nscom.nBEqu + MAX( nSpecies - 1, 0 );

    if ( nscom.chemModel > 0 )
    {
        AllocWorkingSpace();
    }
}

void Chemical::AllocWorkingSpace()
{
    xi_s.resize( nSpecies );
    cs_s.resize( nSpecies );
    vis_s.resize( nSpecies );
    vis_phi_s.resize( nSpecies );
    cp_s.resize( nSpecies );
    dim_cp_s.resize( nSpecies );
    hint_s.resize( nSpecies );
    work_s.resize( nSpecies );
}

void Chemical::ReadGasModel()
{
    FileIO ioFile;

    ioFile.OpenPrjFile( nscom.gasModelFile, ios_base::in );
    string separator = " =\r\n#$,;\"'";
    ioFile.SetDefaultSeparator( separator );

    ioFile.SkipLines( 2 );

    ioFile.ReadNextNonEmptyLine();

    nSpecies  = ioFile.ReadNextDigit< int >();
    nReaction = ioFile.ReadNextDigit< int >();
    Init( nSpecies, nReaction );
    ReadChemical( & ioFile );

    ioFile.CloseFile();
}

void Chemical::Init( int nSpecies, int nReaction )
{
    moleProp->Init( nSpecies );
    reactionRate->Init( nReaction );
    stoichiometric->Init( nReaction, nSpecies );
    blotterCurve->Init( nSpecies );
    thermodynamic->Init( nSpecies );
}

void Chemical::ReadChemical( FileIO * ioFile )
{
    moleProp->Read( ioFile );
    reactionRate->Read( ioFile );
    thermodynamic->Read( ioFile );
    blotterCurve->Read( ioFile );
    stoichiometric->Read( ioFile );
}

void Chemical::Read( DataBook * dataBook )
{
    moleProp->Read( dataBook );
    reactionRate->Read( dataBook );
    thermodynamic->Read( dataBook );
    blotterCurve->Read( dataBook );
    stoichiometric->Read( dataBook );
}

void Chemical::Write( DataBook * dataBook )
{
    moleProp->Write( dataBook );
    reactionRate->Write( dataBook );
    thermodynamic->Write( dataBook );
    blotterCurve->Write( dataBook );
    stoichiometric->Write( dataBook );
}

void Chemical::Init()
{
    InitRefPara();

    CalcRefPara();
}

void Chemical::InitRefPara()
{
    InitGasModel();

    InitWorkingSpace();
}

void Chemical::CalcRefPara()
{
    CalcRefMolecularInfo();

    CalcRefGasInfo();

    CalcRefGama();

    CalcRefSoundSpeed();

    CalcRefMachAndVel();

    CalcStateCoef();

    CalcRefPrim();

    CalcRefReynolds();
}

void Chemical::CalcRefMolecularInfo()
{
    if ( nscom.chemModel > 0 )
    {
        CalcRefMolecularInfoChem();
    }
    else
    {
        CalcRefMolecularInfoAir();
    }
}

void Chemical::CalcRefMolecularInfoAir()
{
    //dimensional average molecular weight
    Real mair = 28.75481;
    Real coef = 1.0e-3;
    nscom.dim_amw  = mair * coef;
    nscom.amw = 1.0;
}

void Chemical::CalcRefMolecularInfoChem()
{
    moleProp->CalcProperty();

    nscom.dim_amw = moleProp->dim_amw;
    nscom.amw = moleProp->amw;
}

void Chemical::CalcRefGama()
{
    if ( nscom.chemModel <= 0 ) return;
    Real t_ref = 1.0;
    CalcDimCps( t_ref, dim_cp_s );

    RealField & mfrac = moleProp->mfrac;
    Real dimCp, dimCv;
    CalcMixtureByMassFraction( mfrac, dim_cp_s, dimCp );

    Real dim_oamw = 1.0 / nscom.dim_amw;

    dimCv = dimCp - rjmk * dim_oamw;
    nscom.gama_ref = dimCp / dimCv;
}

void Chemical::CalcRefSoundSpeed()
{
    Real dim_oamw = 1.0 / nscom.dim_amw;

    nscom.cref_dim = sqrt( nscom.gama_ref * rjmk * dim_oamw * nscom.tref_dim );
}

void Chemical::CalcRefMachAndVel()
{
    if ( nscom.machStrategy == 0 )
    {
        nscom.vref_dim = nscom.cref_dim * nscom.mach_ref;
    }
    else
    {
        nscom.mach_ref = nscom.vref_dim / nscom.cref_dim;
    }
}

void Chemical::CalcStateCoef()
{
    if ( nscom.chemModel > 0 )
    {
        CalcStateCoefChemical();
    }
    else
    {
        CalcStateCoefNs();
    }
}

void Chemical::CalcStateCoefNs()
{
    nscom.statecoef = 1.0 / ( nscom.gama_ref * SQR( nscom.mach_ref ) );
}

void Chemical::CalcStateCoefChemical()
{
    //这个和NS的计算应该是等价的
    nscom.statecoef = rjmk * nscom.tref_dim / ( SQR( nscom.vref_dim ) * nscom.dim_amw );
}

void Chemical::CalcRefPrim()
{
    nscom.dref = 1.0;
    nscom.tref = 1.0;
    nscom.vref = 1.0;
    nscom.pref = nscom.statecoef * nscom.dref * nscom.tref / nscom.amw;
        
    nscom.inflow.resize( nscom.nTEqu );
    nscom.refns.resize( nscom.nTEqu );

    nscom.refns[ IDX::IR ] = nscom.dref;
    nscom.refns[ IDX::IU ] = 1.0 * cos( nscom.aoa ) * cos( nscom.aos );
    nscom.refns[ IDX::IV ] = 1.0 * sin( nscom.aoa ) * cos( nscom.aos );
    nscom.refns[ IDX::IW ] = 1.0                    * sin( nscom.aos );
    nscom.refns[ IDX::IP ] = nscom.pref;

    if ( nscom.chemModel > 0 )
    {
        for ( int iSpecies = 0; iSpecies < nSpecies; ++ iSpecies )
        {
            nscom.refns[ nscom.nBEqu + iSpecies ] = moleProp->mfrac[ iSpecies ];
        }
    }

    for ( int iEqu = 0; iEqu < nscom.nTEqu; ++ iEqu )
    {
        nscom.inflow[ iEqu ] = nscom.refns[ iEqu ];
    }

    if ( ctrl.inflowType == 1 )
    {
        nscom.inflow[ IDX::IU ] = zero;
        nscom.inflow[ IDX::IV ] = zero;
        nscom.inflow[ IDX::IW ] = zero;
    }
}

void Chemical::CalcRefReynolds()
{
    this->CalcDimRefViscosity();

    if ( nscom.chemModel > 0 || nscom.gasInfoStrategy == 0 )
    {
        nscom.reynolds = nscom.dref_dim * nscom.vref_dim * nscom.reylref_dim / nscom.visref_dim;
    }

    nscom.oreynolds = 1.0 / nscom.reynolds;
}

void Chemical::CalcDimRefViscosity()
{
    if ( nscom.chemModel > 0 )
    {
        CalcDimRefViscosityChemical();
    }
    else
    {
        CalcDimRefViscosityNs();
    }

    CalcSutherlandConstant();
}

void Chemical::CalcDimRefViscosityNs()
{
    Real t0 = 273.15; // zero degree celsius temperature(K)
    Real c = 110.4; //dimensional temperature constant in sutherland formula
    Real mu0 = 1.715e-5; //dimensional air viscosity at zero celsius degree

    Real t = nscom.tref_dim;
    Real coef = ( t0 + c ) / ( t + c );
    nscom.visref_dim = coef * pow( t / t0, 1.5 ) * mu0;
}

void Chemical::CalcDimRefViscosityChemical()
{
    CalcMoleFractionByMassFraction( moleProp->mfrac, xi_s );
    CalcDimSpeciesViscosity( nscom.tref, vis_s );
    CalcMixtureByWilkeFormula( xi_s, vis_s, nscom.visref_dim );
}

void Chemical::CalcSutherlandConstant()
{
    nscom.csuth_dim = 110.4;
    nscom.csuth = nscom.csuth_dim / nscom.tref_dim;
}

void Chemical::CalcMoleFractionByMassFraction( RealField & massFrac, RealField & moleFrac )
{
    Real xiSum = 0.0;
    for ( int iSpecies = 0; iSpecies < nSpecies; ++ iSpecies )
    {
        xiSum += massFrac[ iSpecies ] * moleProp->omw[ iSpecies ];
    }
    Real oxiSum = 1.0 / xiSum;

    for ( int iSpecies = 0; iSpecies < nSpecies; ++ iSpecies )
    {
        moleFrac[ iSpecies ] = massFrac[ iSpecies ] * moleProp->omw[ iSpecies ] * oxiSum;
    }
}

void Chemical::CalcDimSpeciesViscosity( Real tm, RealField & vis_s_dim )
{
    Real t1, lt1, lt2, lt3, lt4, tmp;
    t1 = tm * nscom.tref_dim;

    lt1 = log( t1 );
    lt2 = lt1 * lt1;
    lt3 = lt1 * lt2;
    lt4 = lt1 * lt3;

    for ( int iSpecies = 0; iSpecies < nSpecies; ++ iSpecies )
    {
        tmp = blotterCurve->a[ iSpecies ] * lt4
            + blotterCurve->b[ iSpecies ] * lt3
            + blotterCurve->c[ iSpecies ] * lt2
            + blotterCurve->d[ iSpecies ] * lt1
            + blotterCurve->e[ iSpecies ];
        vis_s_dim[ iSpecies ] = 0.10 * exp( tmp ) + SMALL;
    }
    return;
}

void Chemical::CalcMixtureCoefByWilkeFormula( RealField & moleFrac, RealField & var, RealField & phi )
{
    RealField & mw = moleProp->mw;

    for ( int is = 0; is < nSpecies; ++ is )
    {
        phi[ is ] = 0.0;
        for ( int js = 0; js < nSpecies; ++ js )
        {
            Real tmp1 = 1.0 + sqrt( var[ is ] / var[ js ] ) * pow( mw[ js ] / mw[ is ], 0.25 );
            Real tmp2 = sqrt( 8.0 * ( 1.0 + mw[ is ] / mw[ js ] ) );
            phi[ is ] += moleFrac[ js ] * SQR( tmp1 ) / tmp2;
        }
    }
}

void Chemical::CalcMixtureByWilkeFormula( RealField & moleFrac, RealField & mixs, RealField & phi, Real & mixture )
{
    mixture = 0.0;
    for ( int is = 0; is < nSpecies; ++ is )
    {
        mixture += moleFrac[ is ] * mixs[ is ] / phi[ is ];
    }
}

void Chemical::CalcMixtureByWilkeFormula( RealField & moleFrac, RealField & mixs, Real & mixture )
{
    RealField & phi = vis_phi_s;
    CalcMixtureCoefByWilkeFormula( moleFrac, mixs, phi );
    CalcMixtureByWilkeFormula( moleFrac, mixs, phi, mixture );
}

void Chemical::CalcRefGasInfo()
{
    //compute reference pressure and density
    if ( nscom.gasInfoStrategy == 0 )
    {
        SetAirInformationByDataBase();
    }
    else if ( nscom.gasInfoStrategy == 1 )
    {
        Real odim_amw = 1.0 / nscom.dim_amw;
        nscom.pref_dim = nscom.dref_dim * rjmk * odim_amw * nscom.tref_dim;
    }
    else if ( nscom.gasInfoStrategy == 2 )
    {
        Real odim_amw = 1.0 / nscom.dim_amw;
        nscom.dref_dim = nscom.pref_dim / ( rjmk * odim_amw * nscom.tref_dim );
    }
}

void Chemical::SetAirInformationByDataBase()
{
    atmosphere.Init();
    atmosphere.GetAirPara( nscom.elevation );
    nscom.tref_dim = atmosphere.tm;
    nscom.pref_dim = atmosphere.pm;
    nscom.dref_dim = atmosphere.rm;
    nscom.cref_dim = atmosphere.cm;

    NormalizeAirInfo();
}

void Chemical::NormalizeAirInfo()
{
    Real odim_amw = 1.0 / nscom.dim_amw;

    nscom.dref_dim = nscom.pref_dim / ( rjmk * odim_amw * nscom.tref_dim );
}

void Chemical::CalcDimCps( Real tm, RealField & dim_cps )
{
    Real t1, t2, t3, t4;

    t1 = tm * nscom.tref_dim;
    t2 = t1 * t1;
    t3 = t1 * t2;
    t4 = t1 * t3;

    int it;
    thermodynamic->GetTRangeId( t1, it );

    RealField & dim_omw = moleProp->dim_omw;

    for ( int iSpecies = 0; iSpecies < nSpecies; ++ iSpecies )
    {
        RealField & polyCoef = thermodynamic->GetPolyCoef( iSpecies, it );

        dim_cps[ iSpecies ] = polyCoef[ 0 ] +
                              polyCoef[ 1 ] * t1 +
                              polyCoef[ 2 ] * t2 +
                              polyCoef[ 3 ] * t3 +
                              polyCoef[ 4 ] * t4;
        dim_cps[ iSpecies ] *= rjmk * dim_omw[ iSpecies ];
    }
}

void Chemical::CalcMixtureByMassFraction( RealField & cs, RealField & var, Real & mixture )
{
    mixture = 0.0;
    for ( int iSpecies = 0; iSpecies < nSpecies; ++ iSpecies )
    {
        mixture += cs[ iSpecies ] * var[ iSpecies ];
    }
}

void Chemical::CompressData( DataBook *& dataBook )
{
    HXAppend( dataBook, nSpecies );
    HXAppend( dataBook, nReaction );

    Write( dataBook );
}

void Chemical::DecompressData( DataBook * dataBook )
{
    dataBook->MoveToBegin();

    HXRead( dataBook, nSpecies );
    HXRead( dataBook, nReaction );

    Init( nSpecies, nReaction );
    Read( dataBook );
}

void ChemicalCompressData( DataBook *& dataBook )
{
    chem.CompressData( dataBook );
}

void ChemicalDecompressData( DataBook * dataBook )
{
    chem.DecompressData( dataBook );
}




//Chemical::Chemical()
//{
	//INsAlloc();
//}

//Chemical::~Chemical()
//{
	//INsDeAlloc();
//}

//void Chemical::INsAlloc()
//{
	//moleProp = new MolecularProperty();
	//reactionRate = new ReactionRate();
	//stoichiometric = new Stoichiometric();
	//blotterCurve = new BlotterCurve();
	//thermodynamic = new Thermodynamic();
//}

//void Chemical::INsDeAlloc()
//{
	//delete moleProp;
	//delete reactionRate;
	//delete stoichiometric;
	//delete blotterCurve;
	//delete thermodynamic;
//}

void Chemical::INsInitGasModel()
{
	if (inscom.chemModel <= 0) return;

	if (Parallel::pid == Parallel::serverid)
	{
		INsReadGasModel();
	}
	ONEFLOW::HXBcast(ChemicalCompressData, ChemicalDecompressData, Parallel::GetServerid());
}

void Chemical::INsInitWorkingSpace()
{
	inscom.nBEqu = inscom.nEqu;
	if (inscom.chemModel == 0)
	{
		inscom.nSpecies = 0;
	}
	else
	{
		inscom.nSpecies = nSpecies;
	}

	inscom.nTEqu = inscom.nBEqu + MAX(nSpecies - 1, 0);

	if (inscom.chemModel > 0)
	{
		INsAllocWorkingSpace();
	}
}

void Chemical::INsAllocWorkingSpace()
{
	xi_s.resize(nSpecies);
	cs_s.resize(nSpecies);
	vis_s.resize(nSpecies);
	vis_phi_s.resize(nSpecies);
	cp_s.resize(nSpecies);
	dim_cp_s.resize(nSpecies);
	hint_s.resize(nSpecies);
	work_s.resize(nSpecies);
}

void Chemical::INsReadGasModel()
{
	FileIO ioFile;

	ioFile.OpenPrjFile(inscom.gasModelFile, ios_base::in);
	string separator = " =\r\n#$,;\"'";
	ioFile.SetDefaultSeparator(separator);

	ioFile.SkipLines(2);

	ioFile.ReadNextNonEmptyLine();

	nSpecies = ioFile.ReadNextDigit< int >();
	nReaction = ioFile.ReadNextDigit< int >();
	Init(nSpecies, nReaction);
	ReadChemical(&ioFile);

	ioFile.CloseFile();
}

void Chemical::INsInit(int nSpecies, int nReaction)
{
	moleProp->Init(nSpecies);
	reactionRate->Init(nReaction);
	stoichiometric->Init(nReaction, nSpecies);
	blotterCurve->Init(nSpecies);
	thermodynamic->Init(nSpecies);
}

void Chemical::INsReadChemical(FileIO * ioFile)
{
	moleProp->Read(ioFile);
	reactionRate->Read(ioFile);
	thermodynamic->Read(ioFile);
	blotterCurve->Read(ioFile);
	stoichiometric->Read(ioFile);
}

void Chemical::INsRead(DataBook * dataBook)
{
	moleProp->Read(dataBook);
	reactionRate->Read(dataBook);
	thermodynamic->Read(dataBook);
	blotterCurve->Read(dataBook);
	stoichiometric->Read(dataBook);
}

void Chemical::INsWrite(DataBook * dataBook)
{
	moleProp->Write(dataBook);
	reactionRate->Write(dataBook);
	thermodynamic->Write(dataBook);
	blotterCurve->Write(dataBook);
	stoichiometric->Write(dataBook);
}

void Chemical::INsInit()
{
	INsInitRefPara();

	INsCalcRefPara();
}

void Chemical::INsInitRefPara()
{
	INsInitGasModel();

	INsInitWorkingSpace();
}

void Chemical::INsCalcRefPara()
{
	INsCalcRefMolecularInfo();

	INsCalcRefGasInfo();

	INsCalcRefGama();

	INsCalcRefSoundSpeed();

	INsCalcRefMachAndVel();

	INsCalcStateCoef();

	INsCalcRefPrim();

	INsCalcRefReynolds();
}

void Chemical::INsCalcRefMolecularInfo()
{
	if (inscom.chemModel > 0)
	{
		INsCalcRefMolecularInfoChem();
	}
	else
	{
		INsCalcRefMolecularInfoAir();
	}
}

void Chemical::INsCalcRefMolecularInfoAir()
{
	//dimensional average molecular weight
	Real mair = 28.75481;
	Real coef = 1.0e-3;
	inscom.dim_amw = mair * coef;
	inscom.amw = 1.0;
}

void Chemical::INsCalcRefMolecularInfoChem()
{
	moleProp->CalcProperty();

	inscom.dim_amw = moleProp->dim_amw;
	inscom.amw = moleProp->amw;
}

void Chemical::INsCalcRefGama()
{
	if (inscom.chemModel <= 0) return;
	Real t_ref = 1.0;
	INsCalcDimCps(t_ref, dim_cp_s);

	RealField & mfrac = moleProp->mfrac;
	Real dimCp, dimCv;
	INsCalcMixtureByMassFraction(mfrac, dim_cp_s, dimCp);

	Real dim_oamw = 1.0 / inscom.dim_amw;

	dimCv = dimCp - rjmk * dim_oamw;
	inscom.gama_ref = dimCp / dimCv;
}

void Chemical::INsCalcRefSoundSpeed()
{
	Real dim_oamw = 1.0 / inscom.dim_amw;

	inscom.cref_dim = sqrt(inscom.gama_ref * rjmk * dim_oamw * inscom.tref_dim);
}

void Chemical::INsCalcRefMachAndVel()
{
	if (inscom.machStrategy == 0)
	{
		inscom.vref_dim = inscom.cref_dim * inscom.mach_ref;
	}
	else
	{
		inscom.mach_ref = inscom.vref_dim / inscom.cref_dim;
	}
}

void Chemical::INsCalcStateCoef()
{
	if (inscom.chemModel > 0)
	{
		INsCalcStateCoefChemical();
	}
	else
	{
		INsCalcStateCoefNs();
	}
}

void Chemical::INsCalcStateCoefNs()
{
	inscom.statecoef = 1.0 / (inscom.gama_ref * SQR(inscom.mach_ref));
}

void Chemical::INsCalcStateCoefChemical()
{
	//这个和NS的计算应该是等价的
	inscom.statecoef = rjmk * inscom.tref_dim / (SQR(inscom.vref_dim) * inscom.dim_amw);
}

void Chemical::INsCalcRefPrim()
{
	inscom.dref = 1.0;
	inscom.tref = 1.0;
	inscom.vref = 1.0;
	inscom.pref = inscom.statecoef * inscom.dref * inscom.tref / inscom.amw;

	inscom.inflow.resize(inscom.nTEqu);
	inscom.refns.resize(inscom.nTEqu);

	inscom.refns[IIDX::IIR] = inscom.dref;
	inscom.refns[IIDX::IIU] = 1.0 * cos(inscom.aoa) * cos(inscom.aos);
	inscom.refns[IIDX::IIV] = 1.0 * sin(inscom.aoa) * cos(inscom.aos);
	inscom.refns[IIDX::IIW] = 1.0                    * sin(inscom.aos);
	inscom.refns[IIDX::IIP] = inscom.pref;

	if (inscom.chemModel > 0)
	{
		for (int iSpecies = 0; iSpecies < nSpecies; ++iSpecies)
		{
			inscom.refns[inscom.nBEqu + iSpecies] = moleProp->mfrac[iSpecies];
		}
	}

	for (int iEqu = 0; iEqu < inscom.nTEqu; ++iEqu)
	{
		inscom.inflow[iEqu] = inscom.refns[iEqu];
	}

	if (ctrl.inflowType == 1)
	{
		inscom.inflow[IIDX::IIU] = zero;
		inscom.inflow[IIDX::IIV] = zero;
		inscom.inflow[IIDX::IIW] = zero;
	}
}

void Chemical::INsCalcRefReynolds()
{
	this->INsCalcDimRefViscosity();

	if (inscom.chemModel > 0 || inscom.gasInfoStrategy == 0)
	{
		inscom.reynolds = inscom.dref_dim * inscom.vref_dim * inscom.reylref_dim / inscom.visref_dim;
	}

	inscom.oreynolds = 1.0 / inscom.reynolds;
}

void Chemical::INsCalcDimRefViscosity()
{
	if (inscom.chemModel > 0)
	{
		INsCalcDimRefViscosityChemical();
	}
	else
	{
		INsCalcDimRefViscosityNs();
	}

	INsCalcSutherlandConstant();
}

void Chemical::INsCalcDimRefViscosityNs()
{
	Real t0 = 273.15; // zero degree celsius temperature(K)
	Real c = 110.4; //dimensional temperature constant in sutherland formula
	Real mu0 = 1.715e-5; //dimensional air viscosity at zero celsius degree

	Real t = inscom.tref_dim;
	Real coef = (t0 + c) / (t + c);
	inscom.visref_dim = coef * pow(t / t0, 1.5) * mu0;
}

void Chemical::INsCalcDimRefViscosityChemical()
{
	INsCalcMoleFractionByMassFraction(moleProp->mfrac, xi_s);
	INsCalcDimSpeciesViscosity(inscom.tref, vis_s);
	INsCalcMixtureByWilkeFormula(xi_s, vis_s, inscom.visref_dim);
}

void Chemical::INsCalcSutherlandConstant()
{
	inscom.csuth_dim = 110.4;
	inscom.csuth = inscom.csuth_dim / inscom.tref_dim;
}

void Chemical::INsCalcMoleFractionByMassFraction(RealField & massFrac, RealField & moleFrac)
{
	Real xiSum = 0.0;
	for (int iSpecies = 0; iSpecies < nSpecies; ++iSpecies)
	{
		xiSum += massFrac[iSpecies] * moleProp->omw[iSpecies];
	}
	Real oxiSum = 1.0 / xiSum;

	for (int iSpecies = 0; iSpecies < nSpecies; ++iSpecies)
	{
		moleFrac[iSpecies] = massFrac[iSpecies] * moleProp->omw[iSpecies] * oxiSum;
	}
}

void Chemical::INsCalcDimSpeciesViscosity(Real tm, RealField & vis_s_dim)
{
	Real t1, lt1, lt2, lt3, lt4, tmp;
	t1 = tm * inscom.tref_dim;

	lt1 = log(t1);
	lt2 = lt1 * lt1;
	lt3 = lt1 * lt2;
	lt4 = lt1 * lt3;

	for (int iSpecies = 0; iSpecies < nSpecies; ++iSpecies)
	{
		tmp = blotterCurve->a[iSpecies] * lt4
			+ blotterCurve->b[iSpecies] * lt3
			+ blotterCurve->c[iSpecies] * lt2
			+ blotterCurve->d[iSpecies] * lt1
			+ blotterCurve->e[iSpecies];
		vis_s_dim[iSpecies] = 0.10 * exp(tmp) + SMALL;
	}
	return;
}

void Chemical::INsCalcMixtureCoefByWilkeFormula(RealField & moleFrac, RealField & var, RealField & phi)
{
	RealField & mw = moleProp->mw;

	for (int is = 0; is < nSpecies; ++is)
	{
		phi[is] = 0.0;
		for (int js = 0; js < nSpecies; ++js)
		{
			Real tmp1 = 1.0 + sqrt(var[is] / var[js]) * pow(mw[js] / mw[is], 0.25);
			Real tmp2 = sqrt(8.0 * (1.0 + mw[is] / mw[js]));
			phi[is] += moleFrac[js] * SQR(tmp1) / tmp2;
		}
	}
}

void Chemical::INsCalcMixtureByWilkeFormula(RealField & moleFrac, RealField & mixs, RealField & phi, Real & mixture)
{
	mixture = 0.0;
	for (int is = 0; is < nSpecies; ++is)
	{
		mixture += moleFrac[is] * mixs[is] / phi[is];
	}
}

void Chemical::INsCalcMixtureByWilkeFormula(RealField & moleFrac, RealField & mixs, Real & mixture)
{
	RealField & phi = vis_phi_s;
	INsCalcMixtureCoefByWilkeFormula(moleFrac, mixs, phi);
	INsCalcMixtureByWilkeFormula(moleFrac, mixs, phi, mixture);
}

void Chemical::INsCalcRefGasInfo()
{
	//compute reference pressure and density
	if (inscom.gasInfoStrategy == 0)
	{
		INsSetAirInformationByDataBase();
	}
	else if (inscom.gasInfoStrategy == 1)
	{
		Real odim_amw = 1.0 / inscom.dim_amw;
		inscom.pref_dim = inscom.dref_dim * rjmk * odim_amw * inscom.tref_dim;
	}
	else if (inscom.gasInfoStrategy == 2)
	{
		Real odim_amw = 1.0 / inscom.dim_amw;
		inscom.dref_dim = inscom.pref_dim / (rjmk * odim_amw * inscom.tref_dim);
	}
}

void Chemical::INsSetAirInformationByDataBase()
{
	atmosphere.Init();
	atmosphere.GetAirPara(inscom.elevation);
	inscom.tref_dim = atmosphere.tm;
	inscom.pref_dim = atmosphere.pm;
	inscom.dref_dim = atmosphere.rm;
	inscom.cref_dim = atmosphere.cm;

	INsNormalizeAirInfo();
}

void Chemical::INsNormalizeAirInfo()
{
	Real odim_amw = 1.0 / inscom.dim_amw;

	inscom.dref_dim = inscom.pref_dim / (rjmk * odim_amw * inscom.tref_dim);
}

void Chemical::INsCalcDimCps(Real tm, RealField & dim_cps)
{
	Real t1, t2, t3, t4;

	t1 = tm * inscom.tref_dim;
	t2 = t1 * t1;
	t3 = t1 * t2;
	t4 = t1 * t3;

	int it;
	thermodynamic->GetTRangeId(t1, it);

	RealField & dim_omw = moleProp->dim_omw;

	for (int iSpecies = 0; iSpecies < nSpecies; ++iSpecies)
	{
		RealField & polyCoef = thermodynamic->GetPolyCoef(iSpecies, it);

		dim_cps[iSpecies] = polyCoef[0] +
			polyCoef[1] * t1 +
			polyCoef[2] * t2 +
			polyCoef[3] * t3 +
			polyCoef[4] * t4;
		dim_cps[iSpecies] *= rjmk * dim_omw[iSpecies];
	}
}

void Chemical::INsCalcMixtureByMassFraction(RealField & cs, RealField & var, Real & mixture)
{
	mixture = 0.0;
	for (int iSpecies = 0; iSpecies < nSpecies; ++iSpecies)
	{
		mixture += cs[iSpecies] * var[iSpecies];
	}
}

void Chemical::INsCompressData(DataBook *& dataBook)
{
	HXAppend(dataBook, nSpecies);
	HXAppend(dataBook, nReaction);

	INsWrite(dataBook);
}

void Chemical::INsDecompressData(DataBook * dataBook)
{
	dataBook->MoveToBegin();

	HXRead(dataBook, nSpecies);
	HXRead(dataBook, nReaction);

	INsInit(nSpecies, nReaction);
	INsRead(dataBook);
}

//void INsChemicalCompressData(DataBook *& dataBook)
//{
	//chem.INsCompressData(dataBook);
//}

//void INsChemicalDecompressData(DataBook * dataBook)
//{
//	chem.INsDecompressData(dataBook);
//}

EndNameSpace