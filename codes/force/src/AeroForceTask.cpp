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
#include "AeroForceTask.h"
#include "AeroForce.h"
#include "ActionState.h"
#include "Zone.h"
#include "ZoneState.h"
#include "Prj.h"
#include "OStream.h"
#include "UnsGrid.h"
#include "GridUtils.h"
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
#include "FieldWrap.h"
#include "NsCom.h"
#include "NsRhs.h"
#include "HipFluxBackend.h"
#include "HipProductionPolicy.h"
#include "Parallel.h"
#include "Iteration.h"
#include <algorithm>
#include <vector>
#include "TextFileParser.h"
#include "UNsCom.h"

#include "Ctrl.h"
#include "INsIdx.h"
#include "INsCtrl.h"
#include "INsCom.h"
#include "UINsCom.h"
#include "SolverDef.h"

#include <stdexcept>
#include <sstream>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <cstdlib>


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

    this->fileName = GetDataValue< std::string >( "aeroFile" );
}

void AerodynamicForceTask::Dump()
{
    if ( Parallel::pid != Parallel::serverid ) return;

    aeroForceInfo.CalcCoef();

    std::ostringstream oss;

    int wordWidth = 16;
    oss << std::setiosflags( std::ios::right );
    oss << std::setiosflags( std::ios::scientific );

    std::fstream file;
    Prj::OpenPrjFile( file, fileName, std::ios_base::out | std::ios_base::app );

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

        for ( HXSize_t iTitle = 0; iTitle < title.size(); ++ iTitle )
        {
            oss << title[ iTitle ] << std::endl;
        }
    }

    oss << Iteration::outerSteps << "    ";
    oss << Iteration::innerSteps << "    ";
    oss << std::setprecision( 6 ) << ctrl.currTime << "    ";
    oss << std::setprecision( 4 );
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
    oss << std::endl;

    file << oss.str();

	Prj::CloseFile( file );
}

void AerodynamicForceTask::CalcForce()
{
    int idump_pres = 0;
    CalcAeroForce( idump_pres );
}

void CalcAeroForce(int idump_pres)
{
	UnsGrid * grid = Zone::GetUnsGrid();
	BcRecord * bcRecord = &grid->GetFaceTopo().GetBcRecord();
	bcRecord->EnsureBcTypeRegion();

	BcInfo * bcInfo = bcRecord->GetBcInfo();

	int nRegion = bcInfo->bcType.size();

	RealField & xcc = grid->GetCellMesh().xcc;
	RealField & ycc = grid->GetCellMesh().ycc;
	RealField & zcc = grid->GetCellMesh().zcc;
	RealField & vol = grid->GetCellMesh().vol;

	RealField & xfn = grid->GetFaceMesh().xfn;
	RealField & yfn = grid->GetFaceMesh().yfn;
	RealField & zfn = grid->GetFaceMesh().zfn;

	RealField & xfc = grid->GetFaceMesh().xfc;
	RealField & yfc = grid->GetFaceMesh().yfc;
	RealField & zfc = grid->GetFaceMesh().zfc;

	RealField & vfx = grid->GetFaceMesh().vfx;
	RealField & vfy = grid->GetFaceMesh().vfy;
	RealField & vfz = grid->GetFaceMesh().vfz;

	RealField & area = grid->GetFaceMesh().area;

	// Literal "q": aero force is always NS mean-flow (IDX/IIDX), not current-solver
	// solution name. See FieldWrap.h naming convention (cross-solver / NS-only path).
	MRField * q = GetFieldPointer< MRField >(grid, "q");
	MRField * visl = GetFieldPointer< MRField >(grid, "visl");

	// Literal NS cell gradients. The aero force used to read the face-indexed
	// copy bcdqdx/y/z, which Grad::StoreBcGrad filled every RK stage by
	// gathering dqdx/y/z at the left cell of each boundary face. Read that same
	// cell value directly instead: identical numbers, no per-stage face-sized
	// copy, and the HIP path already mirrors exactly these boundary-adjacent
	// cells to the host.
	MRField * dqdx = FieldHome::GetGradField( grid, NS_SOLVER, 1 );
	MRField * dqdy = FieldHome::GetGradField( grid, NS_SOLVER, 2 );
	MRField * dqdz = FieldHome::GetGradField( grid, NS_SOLVER, 3 );
	if ( dqdx == nullptr || dqdy == nullptr || dqdz == nullptr )
	{
		throw std::runtime_error(
			"aero force requires the NS cell gradient fields" );
	}

	int startStrategy = ctrl.startStrategy;

	stress.rey = GetDataValue< Real >("reynolds");
	stress.orey = 1.0 / stress.rey;

	int nSolidCells = GetNumberOfSolidCells( *grid );
	if (nSolidCells == 0) return;

#ifdef ONEFLOW_ENABLE_HIP
	// Holistic aero-input mirror: every host field the force loop reads is
	// refreshed here, in one place, so q / gradient / right-cell inputs can
	// never disagree at the read point.  No-op on the CPU path or when the
	// device does not own the state.
	if ( HipDeviceAuthoritativeStateEnabled() )
	{
		const auto & leftCells = grid->GetFaceTopo().GetLeftCells();
		const auto & rightCells = grid->GetFaceTopo().GetRightCells();

		// The pressure loop reads q[IP][lc]; with vismodel>0 the viscous
		// CorrectGrad also reads q[IU/IV/IW][rc] on the right cell (a ghost
		// at solid boundaries).  Mirror one sorted owner+right cell set for
		// both q and gradients.  includeGhost=true lets the scatter cover rc.
		std::vector< int > aeroCells;
		aeroCells.reserve( static_cast< std::size_t >( grid->nBFaces ) * 2 );
		for ( int face = 0; face < grid->nBFaces; ++ face )
		{
			aeroCells.push_back( leftCells[ face ] );
			aeroCells.push_back( rightCells[ face ] );
		}
		std::sort( aeroCells.begin(), aeroCells.end() );
		aeroCells.erase( std::unique( aeroCells.begin(), aeroCells.end() ),
		                 aeroCells.end() );

		Real * qComp[ 5 ] = {};
		for ( int equation = 0; equation < nscom.nEqu; ++ equation )
			qComp[ equation ] = ( *q )[ equation ].data();
		HipFluxBackend::Shared().MirrorPrimitiveCellsToHost(
			grid, nscom.nEqu, aeroCells.data(),
			static_cast< int >( aeroCells.size() ), qComp, true );

		// The gradient mirror is a hard contract (it throws when no resident
		// gradient exists). On flagship configs where the device gradient
		// stage is not active, the host gradients are already authoritative,
		// so skip the mirror instead of failing.
		if ( HipFluxBackend::Shared().HasResidentGradient(
				 grid, nscom.nEqu ) )
		{
			Real * gx[ 5 ] = {};
			Real * gy[ 5 ] = {};
			Real * gz[ 5 ] = {};
			for ( int equation = 0; equation < nscom.nEqu; ++ equation )
			{
				gx[ equation ] = ( *dqdx )[ equation ].data();
				gy[ equation ] = ( *dqdy )[ equation ].data();
				gz[ equation ] = ( *dqdz )[ equation ].data();
			}
			HipFluxBackend::Shared().MirrorGradientCellsToHost(
				grid, nscom.nEqu, aeroCells, gx, gy, gz );
		}
	}
	// Also refresh the host q mirror + derived gamma/temperature fields for
	// non-force host readers; kept for parity with the pre-existing path.
	NsMirrorHostBoundaryState();
#endif

	if (idump_pres == 1)
	{
		OStream &logger = OStream::Instance();
		logger.ClearAll();

		StringField title;
		title.push_back("title=\"THE FLOW FIELD OF ONEFLOW\"");
		title.push_back("variables=");
		title.push_back("\"x\"");
		title.push_back("\"y\"");
		title.push_back("\"z\"");
		title.push_back("\"-cp\"");
		title.push_back("\"cf\"");
		for (HXSize_t i = 0; i < title.size(); ++i)
		{
			logger << title[i] << "\n";
		}

		logger << "Zone  i = " << nSolidCells << " \n";
	}

	for (int ir = 0; ir < nRegion; ++ir)
	{
		int bcType = bcInfo->bcType[ir];
		if (bcType != BC::SOLID_SURFACE) continue;
		int nBCFace = bcInfo->bcFace[ir].size();

		AeroForce aeroForce;
		for (int iBCFace = 0; iBCFace < nBCFace; ++iBCFace)
		{
			int fId = bcInfo->bcFace[ir][iBCFace];
			int lc = grid->GetFaceTopo().GetLeftCells()[fId];
			int rc = grid->GetFaceTopo().GetRightCells()[fId];
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
				Real wp = (*q)[IIDX::IIP][lc] + (*dqdx)[IIDX::IIP][lc] * dx + (*dqdy)[IIDX::IIP][lc] * dy + (*dqdz)[IIDX::IIP][lc] * dz;
				if (wp < 0.0) wp = (*q)[IIDX::IIP][lc];
				Real pref = nscom.inflow[IIDX::IIP];
				Real cp = two * (wp - pref);

				Real dpx = stress.fanx * cp;
				Real dpy = stress.fany * cp;
				Real dpz = stress.fanz * cp;

				aeroForce.pres.x = dpx;
				aeroForce.pres.y = dpy;
				aeroForce.pres.z = dpz;

				if (vis_model.vismodel > 0)
				{
					stress.dudx = (*dqdx)[IIDX::IIU][lc];
					stress.dudy = (*dqdy)[IIDX::IIU][lc];
					stress.dudz = (*dqdz)[IIDX::IIU][lc];

					stress.dvdx = (*dqdx)[IIDX::IIV][lc];
					stress.dvdy = (*dqdy)[IIDX::IIV][lc];
					stress.dvdz = (*dqdz)[IIDX::IIV][lc];

					stress.dwdx = (*dqdx)[IIDX::IIW][lc];
					stress.dwdy = (*dqdy)[IIDX::IIW][lc];
					stress.dwdz = (*dqdz)[IIDX::IIW][lc];

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
					OStream &logger = OStream::Instance();
					logger << std::setiosflags(std::ios::left);
					logger << std::setiosflags(std::ios::scientific);
					logger << std::setprecision(10);
					logger << std::setw(wordWidth) << xc;
					logger << std::setw(wordWidth) << yc;
					logger << std::setw(wordWidth) << zc;
					logger << std::setw(wordWidth) << -cp;
					logger << std::setw(wordWidth) << cf;
					logger << std::endl;
				}
			}
			else
			{
				// Optional host aero-input trace (off by default): pins down whether
				// the 21% DCU-vs-CPU aero gap comes from stale host q/gradient mirror.
				if ( std::getenv( "ONEFLOW_AERO_INPUT_PROBE" ) != nullptr )
				{
					static int aeroProbeCount = 0;
					if ( aeroProbeCount < 8 )
					{
						++aeroProbeCount;
						const Real qip = (*q)[IDX::IP][lc];
						const Real gpx = (*dqdx)[IDX::IP][lc];
						const Real gpy = (*dqdy)[IDX::IP][lc];
						const Real gpz = (*dqdz)[IDX::IP][lc];
						const Real wpTrace = qip + gpx * dx + gpy * dy + gpz * dz;
						std::cout << std::setprecision( 17 );
						std::cout << "AERO_INPUT step=" << Iteration::outerSteps
						          << " fId=" << fId << " lc=" << lc
						          << " qip=" << qip
						          << " gpx=" << gpx << " gpy=" << gpy << " gpz=" << gpz
						          << " dx=" << dx << " dy=" << dy << " dz=" << dz
						          << " wp=" << wpTrace
						          << " fanx=" << stress.fanx << " fany=" << stress.fany
						          << " fanz=" << stress.fanz
						          << std::endl;
					}
				}
				Real wp = (*q)[IDX::IP][lc] + (*dqdx)[IDX::IP][lc] * dx + (*dqdy)[IDX::IP][lc] * dy + (*dqdz)[IDX::IP][lc] * dz;
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
					stress.dudx = (*dqdx)[IDX::IU][lc];
					stress.dudy = (*dqdy)[IDX::IU][lc];
					stress.dudz = (*dqdz)[IDX::IU][lc];

					stress.dvdx = (*dqdx)[IDX::IV][lc];
					stress.dvdy = (*dqdy)[IDX::IV][lc];
					stress.dvdz = (*dqdz)[IDX::IV][lc];

					stress.dwdx = (*dqdx)[IDX::IW][lc];
					stress.dwdy = (*dqdy)[IDX::IW][lc];
					stress.dwdz = (*dqdz)[IDX::IW][lc];

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

					{
						const char * wf =
							std::getenv( "ONEFLOW_WALL_FORCE_PROBE" );
						static std::ofstream wout;
						if ( wf != nullptr && ! wout.is_open() )
						{
							wout.open( wf, std::ios::out );
						}
						if ( wout.is_open() )
						{
							wout << fId << ' ' << lc << ' ' << rc;
							wout << ' ' << (*q)[IDX::IU][lc] << ' ' << (*q)[IDX::IU][rc];
							wout << ' ' << (*q)[IDX::IV][lc] << ' ' << (*q)[IDX::IV][rc];
							wout << ' ' << (*q)[IDX::IW][lc] << ' ' << (*q)[IDX::IW][rc];
							wout << ' ' << stress.dudx << ' ' << stress.dudy << ' ' << stress.dudz;
							wout << ' ' << stress.dvdx << ' ' << stress.dvdy << ' ' << stress.dvdz;
							wout << ' ' << stress.dwdx << ' ' << stress.dwdy << ' ' << stress.dwdz;
							wout << ' ' << stress.viscosity;
							wout << ' ' << stress.fanx << ' ' << stress.fany << ' ' << stress.fanz;
							wout << ' ' << stress.area << ' ' << stress.orey;
							wout << '\n';
						}
					}

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
					OStream &logger = OStream::Instance();
					logger << std::setiosflags(std::ios::left);
					logger << std::setiosflags(std::ios::scientific);
					logger << std::setprecision(10);
					logger << std::setw(wordWidth) << xc;
					logger << std::setw(wordWidth) << yc;
					logger << std::setw(wordWidth) << zc;
					logger << std::setw(wordWidth) << -cp;
					logger << std::setw(wordWidth) << cf;
					logger << std::endl;
				}
			}
		}
		if (idump_pres == 1)
		{
			OStream &logger = OStream::Instance();
			ToDataBook(ActionState::dataBook, logger);
		}
	}
}

EndNameSpace
