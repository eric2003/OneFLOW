/*---------------------------------------------------------------------------*\
    OneFLOW - LargeScale Multiphysics Scientific Simulation Environment
    Copyright (C) 2017-2021 He Xin and the OneFLOW contributors.
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

#include "ScalarField.h"
#include "FieldPara.h"
#include "Prj.h"
#include "Constant.h"
#include "HXMath.h"
#include "HXCgns.h"
#include "FileUtil.h"
#include <iostream>
#include <fstream>
#include <vector>
using namespace std;

BeginNameSpace( ONEFLOW )

ScalarField::ScalarField()
{
	;
}

ScalarField::~ScalarField()
{
	;
}

void ScalarField::InitFlowField( ScalarGrid * grid )
{
    this->grid = grid;

    int nTCells = grid->GetNTCells();
    q.Resize( nTCells );
    res.Resize( nTCells );

    for ( int iCell = 0; iCell < nTCells; ++ iCell )
    {
        Real xm = grid->xcc[ iCell ];
        q[ iCell ] = this->ScalarFun( xm );
    }
    res = 0;
    qInf = q[ 0 ]; // farfield value
    int nFaces = grid->GetNFaces();

    qL.Resize( nFaces );
    qR.Resize( nFaces );

    invflux.Resize( nFaces );
}

Real ScalarField::ScalarFun( Real xm )
{
    return SquareFun( xm );
}

Real ScalarField::SquareFun( Real xm )
{
    if ( xm >= 0.5 && xm <= 1.0 )
    {
        return 2.0;
    }
    return 1.0;
}

void ScalarField::ToTecplot( RealList & varlist, string const & fileName )
{
    int nCells = grid->GetNCells();

    fstream file;
    OpenPrjFile( file, fileName, ios_base::out );

    int nSize = nCells;
    file << "TITLE = " << "\"OneFLOW X-Y Plot\"" << "\n";
    file << "VARIABLES = " << "\"x\", " << "\"u\" " << "\n";
    file << "ZONE T = " << "\"Scalar Results\"," << " I = " << nSize << "\n";

    for ( int iCell = 0; iCell < nCells; ++ iCell )
    {
        Real xm = grid->xcc[ iCell ];
        Real fm = varlist[ iCell ];
        file << xm << " " << fm << "\n";
    }

    CloseFile( file );
}

void ScalarField::SolveFlowField( FieldPara * para )
{
    this->para = para;
    for ( int n = 0; n < para->nt; ++ n )
    {
        //cout << " iStep = " << n << " nStep = " << para->nt << "\n";

        this->SolveOneStep();
    }

    this->Visual();
}

void ScalarField::Visual()
{
    Real time = para->dt * para->nt;
    RealList theory;
    Theory( time, theory );
    ToTecplot( q, "OneFLOW.plt" );
    ToTecplot( theory, "theory.plt" );
}

void ScalarField::Theory( Real time, RealList & theory )
{
    int nCells = grid->GetNCells();
    theory.Resize( nCells );

    Real xs = para->c * time;

    for ( int iCell = 0; iCell < nCells; ++ iCell )
    {
        Real xm = grid->xcc[ iCell ];
        Real xm_new = xm - xs;
        theory[ iCell ] = this->ScalarFun( xm_new );
    }
}

void ScalarField::SolveOneStep()
{
    this->Boundary();
    this->GetQLQR();
    this->CalcInvFlux();
    this->UpdateResidual();
    this->TimeIntergral();
    this->Update();
}

void ScalarField::Update()
{
    int nCells = grid->GetNCells();

    for ( int iCell = 0; iCell < nCells; ++ iCell )
    {
        q[ iCell ] += res[ iCell ];
    }
}

void ScalarField::TimeIntergral()
{
    int nCells = grid->GetNCells();
    for ( int iCell = 0; iCell < nCells; ++ iCell )
    {
        Real ovol = 1.0 / grid->vol[ iCell ];
        Real coef = para->dt * ovol;
        res[ iCell ] *= coef;
    }
}

void ScalarField::GetQLQR()
{
    int nFaces = grid->GetNFaces();
    for ( int iFace = 0; iFace < nFaces; ++ iFace )
    {
        int lc = grid->lc[ iFace ];
        int rc = grid->rc[ iFace ];

        qL[ iFace ] = q[ lc ];
        qR[ iFace ] = q[ rc ];
    }
}

void ScalarField::CalcInvFlux()
{
    int nFaces = grid->GetNFaces();
    Real vxl = 1.0;
    Real vyl = 0.0;
    Real vzl = 0.0;

    Real vxr = 1.0;
    Real vyr = 0.0;
    Real vzr = 0.0;
    for ( int iFace = 0; iFace < nFaces; ++ iFace )
    {
        Real q_L = qL[ iFace ];
        Real q_R = qR[ iFace ];

        Real vnl  = grid->xfn[ iFace ] * vxl + grid->yfn[ iFace ] * vyl + grid->zfn[ iFace ] * vzl;
        Real vnr  = grid->xfn[ iFace ] * vxr + grid->yfn[ iFace ] * vyr + grid->zfn[ iFace ] * vzr;

        Real eigenL = vnl;
        Real eigenR = vnr;

        eigenL = half * ( eigenL + ABS( eigenL ) );
        eigenR = half * ( eigenR - ABS( eigenR ) );

        Real fL = q_L * eigenL;
        Real fR = q_R * eigenR;
        Real fM = fL + fR;

        Real area = grid->area[ iFace ];
        invflux[ iFace ] = fM * area;
    }
}

void ScalarField::Boundary()
{
    int nBFaces = grid->GetNBFaces();
    for ( int iFace = 0; iFace < nBFaces; ++ iFace )
    {
        int bcType = grid->bcTypes[ iFace ];
        int lc = grid->lc[ iFace ];
        int rc = grid->rc[ iFace ];
        if ( bcType == ONEFLOW::BCInflow )
        {
            this->q[ rc ] = qInf;
        }
        else if ( bcType == ONEFLOW::BCOutflow )
        {
            this->q[ rc ] = this->q[ lc ];
        }

    }
}

void ScalarField::UpdateResidual()
{
    this->res = 0;
    this->AddF2CField( this->res, this->invflux );
}

void ScalarField::AddF2CField( RealList & cField, RealList & fField )
{
    int nFaces = grid->GetNFaces();
    int nBFaces = grid->GetNBFaces();
    for ( int iFace = 0; iFace < nBFaces; ++ iFace )
    {
        int lc = grid->lc[ iFace ];
        cField[ lc ] -= fField[ iFace ];
    }
    
    for ( int iFace = nBFaces; iFace < nFaces; ++ iFace )
    {
        int lc = grid->lc[ iFace ];
        int rc = grid->rc[ iFace ];

        cField[ lc ] -= fField[ iFace ];
        cField[ rc ] += fField[ iFace ];
    }
}


EndNameSpace