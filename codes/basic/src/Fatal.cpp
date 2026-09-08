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
#include "Fatal.h"
#include <iostream>
#include <string>
#include <stdexcept>

BeginNameSpace( ONEFLOW )

void StopProgramFunction( const std::string & stopInformation,
                          const std::string & fileName,
                          const int & fileLine,
                          const std::string & dateName,
                          const std::string & timeName )
{
    std::ostringstream oss;
    oss << "\n++++++++++++++++++ Fatal Error +++++++++++++++++++++++++++++\n"
        << stopInformation << "\n"
        << " File   : " << fileName << "\n"
        << " Line   : " << fileLine << "\n"
        << " Compiled on " << dateName << " at " << timeName << "\n"
        << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";

    throw std::runtime_error( oss.str() );
}

EndNameSpace