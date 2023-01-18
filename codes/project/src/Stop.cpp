/*---------------------------------------------------------------------------*\
    OneFLOW - LargeScale Multiphysics Scientific Simulation Environment
    Copyright (C) 2017-2023 He Xin and the OneFLOW contributors.
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

BeginNameSpace( ONEFLOW )

void StopProgramFunction( const std::string & stopInformation, const std::string & fileName, const int & fileLine, const std::string & dateName, const std::string & timeName )
{
    std::cout << std::endl;
    std::cout << "++++++++++++++++++Stop Information  +++++++++++++++++++++++++++++\n";
    std::cout <<  stopInformation << std::endl;
    std::cout << " The stop filename is : " << fileName << std::endl;
    std::cout << " at line " << fileLine << std::endl;
    std::cout << " compiled on " << dateName << " at " << timeName << std::endl;
    std::cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
    exit( 0 );
}


EndNameSpace
