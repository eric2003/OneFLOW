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

#include "ElemFeature.h"
#include "FaceSolver.h"
#include "ElementHome.h"
#include <iostream>
using namespace std;

BeginNameSpace( ONEFLOW )

int EState::SPLIT    = 0;
int EState::NOCHANGE = 1;
int EState::MERGE    = 2;
int EState::DELETE   = 3;
int EState::HIDDEN   = 4;
int EState::ON       = 5;
int EState::BC_FACE  = 6;
int EState::H_FACE   = 7;
int EState::G_FACE   = 8;

EState::EState()
{
    ;
}

EState::~EState()
{
    ;
}

ElemFeature::ElemFeature()
{
    this->eType = new IntField();
}

ElemFeature::~ElemFeature()
{
    delete this->eType;
}

void ElemFeature::ScanElements()
{
    int nElement = this->eType->size();

    cout << " nElement = " << nElement << endl;

    int nIo = 200000;

    int iCount = 0;
    for ( int eId = 0; eId < nElement; ++ eId )
    {
        int eType = ( * this->eType )[ eId ];
        ++ iCount;
        if ( ! ONEFLOW::IsBasicVolumeElementType( eType ) )
        {
            continue;
        }
        //cout << "eId = " << eId << "\n";

        if ( eId % nIo == 0 ) 
        {
            cout << " eId = " << eId << " Total Element Number = " << nElement << endl;
        }

        this->face_solver->ScanElementFace( this->eNodeId[ eId ], eType, eId );
    }

    cout << " ScanElements Face number = " << this->face_solver->GetNSimpleFace() << endl;
}


EndNameSpace