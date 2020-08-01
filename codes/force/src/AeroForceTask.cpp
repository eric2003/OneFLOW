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
#include "AeroForceTask.h"
#include "AeroForce.h"
#include "ActionState.h"
#include "Zone.h"
#include "ZoneState.h"
#include "Prj.h"
#include "OStream.h"
#include "UnsGrid.h"
#include "BcRecord.h"
#include "FaceTopo.h"
#include "Boundary.h"
#include "FaceMesh.h"
#include "CellMesh.h"
#include "DataBase.h"
#include "NsIdx.h"
#include "HXMath.h"
#include "NsCtrl.h"
#include "Stress.h"
#include "VisGrad.h"
#include "NsCom.h"
#include "Parallel.h"
#include "Iteration.h"
#include "FileUtil.h"
#include "FileIO.h"
#include "UNsCom.h"

#include "Ctrl.h"
#include "INsIdx.h"
#include "INsCtrl.h"
#include "INsCom.h"
#include "UINsCom.h"
#include "SolverDef.h"
#include <sstream>
#include <iomanip>
#include <iostream>
using namespace std;

BeginNameSpace( ONEFLOW )

AerodynamicForceTask::AerodynamicForceTask()
{
    ;
}

AerodynamicForceTask::~AerodynamicForceTask()
{
    ;
}

void AerodynamicForceTask::Run()
{
    this->Init();

    for ( int zId = 0; zId < ZoneState::nZones; ++ zId )
    {
        if (  ! ZoneState::IsValidZone( zId ) ) continue;
        ZoneState::zid = zId;

        this->CalcForce();
    }

    aeroForceInfo.CollectForce();
    this->Dump();
}

void AerodynamicForceTask::Init()
{
    aeroForceInfo.Init();

    this->fileName = GetDataValue< string >( "aeroFile" );
}

void AerodynamicForceTask::Dump()
{
    if ( Parallel::pid != Parallel::serverid ) return;

    aeroForceInfo.CalcCoef();

    ostringstream oss;

    int wordWidth = 16;
    oss << setiosflags( ios::right );
    oss << setiosflags( ios::scientific );

    fstream file;
    OpenPrjFile( file, fileName, ios_base::out | ios_base::app );

    if ( IsEmpty( file ) )
    {
        StringField title;
        title.push_back( "Title=\"Aerodynamic Force\"" );
        title.push_back( "Variables=" );
        title.push_back( "\"iter\"" );
        title.push_back( "\"sub-iter\"" );
        title.push_back( "\"time\"" );    
        title.push_back( "\"CL\"" );    
        title.push_back( "\"CD\"" );    
        title.push_back( "\"CD_PR\"" );    
        title.push_back( "\"CD_SF\"" );    
        title.push_back( "\"CDL2\"" );    
        title.push_back( "\"Xcp\"" );    
        title.push_back( "\"Fx\"" );    
        title.push_back( "\"Fy\"" );
        title.push_back( "\"Fz\"" );
        title.push_back( "\"Cmx\"" );
        title.push_back( "\"Cmy\"" );
        title.push_back( "\"Cmz\"" );

        for ( UInt iTitle = 0; iTitle < title.size(); ++ iTitle )
        {
            oss << title[ iTitle ] << endl;
        }
    }

    oss << Iteration::outerSteps << "    ";
    oss << Iteration::innerSteps << "    ";
    oss << setprecision( 6 ) << ctrl.currTime << "    ";
    oss << setprecision( 4 );
    oss << aeroForceInfo.cl << "    ";
    oss << aeroForceInfo.cd << "    ";
    oss << aeroForceInfo.cd_pres << "    ";
    oss << aeroForceInfo.cd_vis << "    ";
    oss << aeroForceInfo.cdl << "    ";
    oss << aeroForceInfo.pres_center << "    ";
    oss << aeroForceInfo.cf.x << "    ";
    oss << aeroForceInfo.cf.y << "    ";
    oss << aeroForceInfo.cf.z << "    ";
    oss << aeroForceInfo.cmom.x << "    ";
    oss << aeroForceInfo.cmom.y << "    ";
    oss << aeroForceInfo.cmom.z << "    ";
    oss << endl;

    file << oss.str();

    CloseFile( file );
}

void AerodynamicForceTask::CalcForce()
{
    int idump_pres = 0;
    CalcAeroForce( idump_pres );
}

int GetNSolidCell( UnsGrid * grid )
{
    BcRecord * bcRecord = grid->faceTopo->bcManager->bcRecord;
    bcRecord->CreateBcRegion();

    BcInfo * bcInfo = bcRecord->bcInfo;

    int nRegion = bcInfo->bcType.size();

    int nSolidCell = 0;
    for ( int ir = 0; ir < nRegion; ++ ir )
    {
        int bcType = bcInfo->bcType[ ir ];
        if ( bcType != BC::SOLID_SURFACE ) continue;
        bcInfo->bcFace.size();

        int nBCFace = bcInfo->bcFace[ ir ].size();
        nSolidCell += nBCFace;
    }

    return nSolidCell;
}

void CalcAeroForce(int idump_pres)
{
	UnsGrid * grid = Zone::GetUnsGrid();
	BcRecord * bcRecord = grid->faceTopo->bcManager->bcRecord;
	bcRecord->CreateBcRegion();

	BcInfo * bcInfo = bcRecord->bcInfo;

	int nRegion = bcInfo->bcType.size();

	RealField & xcc = grid->cellMesh->xcc;
	RealField & ycc = grid->cellMesh->ycc;
	RealField & zcc = grid->cellMesh->zcc;
	RealField & vol = grid->cellMesh->vol;

	RealField & xfn = grid->faceMesh->xfn;
	RealField & yfn = grid->faceMesh->yfn;
	RealField & zfn = grid->faceMesh->zfn;

	RealField & xfc = grid->faceMesh->xfc;
	RealField & yfc = grid->faceMesh->yfc;
	RealField & zfc = grid->faceMesh->zfc;

	RealField & vfx = grid->faceMesh->vfx;
	RealField & vfy = grid->faceMesh->vfy;
	RealField & vfz = grid->faceMesh->vfz;

	RealField & area = grid->faceMesh->area;

	MRField * q = GetFieldPointer< MRField >(grid, "q");
	MRField * visl = GetFieldPointer< MRField >(grid, "visl");

	MRField * bcdqdx = GetFieldPointer< MRField >(grid, "bcdqdx");
	MRField * bcdqdy = GetFieldPointer< MRField >(grid, "bcdqdy");
	MRField * bcdqdz = GetFieldPointer< MRField >(grid, "bcdqdz");

	int startStrategy = ONEFLOW::GetDataValue< int >("startStrategy");
	if (startStrategy == 2)
	{
		RealField & dpdx = (*bcdqdx)[IIDX::IIP];
		RealField & dpdy = (*bcdqdy)[IIDX::IIP];
		RealField & dpdz = (*bcdqdz)[IIDX::IIP];
	}
	else
	{
		RealField & dpdx = (*bcdqdx)[IDX::IP];
		RealField & dpdy = (*bcdqdy)[IDX::IP];
		RealField & dpdz = (*bcdqdz)[IDX::IP];
	}

	stress.rey = GetDataValue< Real >("reynolds");
	stress.orey = 1.0 / stress.rey;

	int nSolidCell = GetNSolidCell(grid);
	if (nSolidCell == 0) return;

	if (idump_pres == 1)
	{
		StrIO.ClearAll();

		StringField title;
		title.push_back("title=\"THE FLOW FIELD OF ONEFLOW\"");
		title.push_back("variables=");
		title.push_back("\"x\"");
		title.push_back("\"y\"");
		title.push_back("\"z\"");
		title.push_back("\"-cp\"");
		title.push_back("\"cf\"");
		for (UInt i = 0; i < title.size(); ++i)
		{
			StrIO << title[i] << "\n";
		}

		StrIO << "Zone  i = " << nSolidCell << " \n";
	}

	for (int ir = 0; ir < nRegion; ++ir)
	{
		int bcType = bcInfo->bcType[ir];
		if (bcType != BC::SOLID_SURFACE) continue;
		bcInfo->bcFace.size();
		int nBCFace = bcInfo->bcFace[ir].size();

		AeroForce aeroForce;
		for (int iBCFace = 0; iBCFace < nBCFace; ++iBCFace)
		{
			int fId = bcInfo->bcFace[ir][iBCFace];
			int lc = grid->faceTopo->lCell[fId];
			int rc = grid->faceTopo->rCell[fId];
			stress.area = area[fId];
			stress.fnx = xfn[fId];
			stress.fny = yfn[fId];
			stress.fnz = zfn[fId];
			stress.fanx = xfn[fId] * area[fId];
			stress.fany = yfn[fId] * area[fId];
			stress.fanz = zfn[fId] * area[fId];

			Real dx = xfc[fId] - xcc[lc];
			Real dy = yfc[fId] - ycc[lc];
			Real dz = zfc[fId] - zcc[lc];

			//pressure drag


			if (startStrategy == 2)
			{
				Real wp = (*q)[IIDX::IIP][lc] + (*bcdqdx)[IIDX::IIP][fId] * dx + (*bcdqdy)[IIDX::IIP][fId] * dy + (*bcdqdz)[IIDX::IIP][fId] * dz;
				if (wp < 0.0) wp = (*q)[IIDX::IIP][lc];
				Real pref = inscom.inflow[IIDX::IIP];
				Real cp = two * (wp - pref);

				Real dpx = stress.fanx * cp;
				Real dpy = stress.fany * cp;
				Real dpz = stress.fanz * cp;

				aeroForce.pres.x = dpx;
				aeroForce.pres.y = dpy;
				aeroForce.pres.z = dpz;

				if (vis_model.vismodel > 0)
				{
					stress.dudx = (*bcdqdx)[IIDX::IIU][fId];
					stress.dudy = (*bcdqdy)[IIDX::IIU][fId];
					stress.dudz = (*bcdqdz)[IIDX::IIU][fId];

					stress.dvdx = (*bcdqdx)[IIDX::IIV][fId];
					stress.dvdy = (*bcdqdy)[IIDX::IIV][fId];
					stress.dvdz = (*bcdqdz)[IIDX::IIV][fId];

					stress.dwdx = (*bcdqdx)[IIDX::IIW][fId];
					stress.dwdy = (*bcdqdy)[IIDX::IIW][fId];
					stress.dwdz = (*bcdqdz)[IIDX::IIW][fId];

					//gradient correction
					dx = xcc[rc] - xcc[lc];
					dy = ycc[rc] - ycc[lc];
					dz = zcc[rc] - zcc[lc];

					Real ods = 1.0 / DIST(dx, dy, dz);
					dx *= ods;
					dy *= ods;
					dz *= ods;

					CorrectGrad((*q)[IIDX::IIU][lc], (*q)[IIDX::IIU][rc], stress.dudx, stress.dudy, stress.dudz, dx, dy, dz, ods);
					CorrectGrad((*q)[IIDX::IIV][lc], (*q)[IIDX::IIV][rc], stress.dvdx, stress.dvdy, stress.dvdz, dx, dy, dz, ods);
					CorrectGrad((*q)[IIDX::IIW][lc], (*q)[IIDX::IIW][rc], stress.dwdx, stress.dwdy, stress.dwdz, dx, dy, dz, ods);

					stress.viscosity = (*visl)[0][lc];

					stress.CalcForce(&aeroForce.vis);
				}
				aeroForce.SumForce();

				aeroForce.vfx = vfx[fId];
				aeroForce.vfy = vfy[fId];
				aeroForce.vfz = vfz[fId];

				aeroForce.CalcPower();

				Real xc = xfc[fId];
				Real yc = yfc[fId];
				Real zc = zfc[fId];

				aeroForce.CalcMoment(xc, yc, zc);
				aeroForceInfo.totalForce.AddForce(&aeroForce);

				if (idump_pres == 1)
				{
					Real cf = aeroCom.CalcCF(&aeroForce.vis, area[fId]);

					int wordWidth = 20;
					StrIO << setiosflags(ios::left);
					StrIO << setiosflags(ios::scientific);
					StrIO << setprecision(10);
					StrIO << setw(wordWidth) << xc;
					StrIO << setw(wordWidth) << yc;
					StrIO << setw(wordWidth) << zc;
					StrIO << setw(wordWidth) << -cp;
					StrIO << setw(wordWidth) << cf;
					StrIO << endl;
				}
			}
			else
			{
				Real wp = (*q)[IDX::IP][lc] + (*bcdqdx)[IDX::IP][fId] * dx + (*bcdqdy)[IDX::IP][fId] * dy + (*bcdqdz)[IDX::IP][fId] * dz;
				if (wp < 0.0) wp = (*q)[IDX::IP][lc];
				Real pref = nscom.inflow[IDX::IP];
				Real cp = two * (wp - pref);

				Real dpx = stress.fanx * cp;
				Real dpy = stress.fany * cp;
				Real dpz = stress.fanz * cp;

				aeroForce.pres.x = dpx;
				aeroForce.pres.y = dpy;
				aeroForce.pres.z = dpz;

				if (vis_model.vismodel > 0)
				{
					stress.dudx = (*bcdqdx)[IDX::IU][fId];
					stress.dudy = (*bcdqdy)[IDX::IU][fId];
					stress.dudz = (*bcdqdz)[IDX::IU][fId];

					stress.dvdx = (*bcdqdx)[IDX::IV][fId];
					stress.dvdy = (*bcdqdy)[IDX::IV][fId];
					stress.dvdz = (*bcdqdz)[IDX::IV][fId];

					stress.dwdx = (*bcdqdx)[IDX::IW][fId];
					stress.dwdy = (*bcdqdy)[IDX::IW][fId];
					stress.dwdz = (*bcdqdz)[IDX::IW][fId];

					//gradient correction
					dx = xcc[rc] - xcc[lc];
					dy = ycc[rc] - ycc[lc];
					dz = zcc[rc] - zcc[lc];

					Real ods = 1.0 / DIST(dx, dy, dz);
					dx *= ods;
					dy *= ods;
					dz *= ods;

					CorrectGrad((*q)[IDX::IU][lc], (*q)[IDX::IU][rc], stress.dudx, stress.dudy, stress.dudz, dx, dy, dz, ods);
					CorrectGrad((*q)[IDX::IV][lc], (*q)[IDX::IV][rc], stress.dvdx, stress.dvdy, stress.dvdz, dx, dy, dz, ods);
					CorrectGrad((*q)[IDX::IW][lc], (*q)[IDX::IW][rc], stress.dwdx, stress.dwdy, stress.dwdz, dx, dy, dz, ods);

					stress.viscosity = (*visl)[0][lc];

					stress.CalcForce(&aeroForce.vis);
				}
				aeroForce.SumForce();

				aeroForce.vfx = vfx[fId];
				aeroForce.vfy = vfy[fId];
				aeroForce.vfz = vfz[fId];

				aeroForce.CalcPower();

				Real xc = xfc[fId];
				Real yc = yfc[fId];
				Real zc = zfc[fId];

				aeroForce.CalcMoment(xc, yc, zc);
				aeroForceInfo.totalForce.AddForce(&aeroForce);

				if (idump_pres == 1)
				{
					Real cf = aeroCom.CalcCF(&aeroForce.vis, area[fId]);

					int wordWidth = 20;
					StrIO << setiosflags(ios::left);
					StrIO << setiosflags(ios::scientific);
					StrIO << setprecision(10);
					StrIO << setw(wordWidth) << xc;
					StrIO << setw(wordWidth) << yc;
					StrIO << setw(wordWidth) << zc;
					StrIO << setw(wordWidth) << -cp;
					StrIO << setw(wordWidth) << cf;
					StrIO << endl;
				}
			}
		}
		if (idump_pres == 1)
		{
			ToDataBook(ActionState::dataBook, StrIO);
		}
	}
}

EndNameSpace