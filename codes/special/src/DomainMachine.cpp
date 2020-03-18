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

#include "DomainMachine.h"
#include "FileIO.h"
#include "HXMath.h"
#include <iostream>
using namespace std;

BeginNameSpace( ONEFLOW )

DomainMachine domain_Machine;

DomainMachine::DomainMachine()
{
}

DomainMachine::~DomainMachine()
{
}

void DomainMachine::AddBcType( FileIO * ioFile )
{
    int id = ioFile->ReadNextDigit< int >();
    int bctype = ioFile->ReadNextDigit< int >();
    bctypeList.push_back( bctype );
}


int DomainMachine::GetBcType( int id )
{
    int ida = ABS( id ) - 1;
    return this->bctypeList[ ida ];
}




EndNameSpace