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
#include "Stop.h"
#include <iostream>
using namespace std;

BeginNameSpace( ONEFLOW )

void StopProgramFunction( const string & stopInformation, const string & fileName, const int & fileLine, const string & dateName, const string & timeName )
{
    cout << endl;
    cout << "++++++++++++++++++Stop Information  +++++++++++++++++++++++++++++\n";
    cout <<  stopInformation << endl;
    cout << " The stop filename is : " << fileName << endl;
    cout << " at line " << fileLine << endl;
    cout << " compiled on " << dateName << " at " << timeName << endl;
    cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
    exit( 0 );
}


EndNameSpace
