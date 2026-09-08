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
#pragma once
#include "NamespaceMacros.h"
#include <string>
#include <stdexcept>
#include <sstream>

BeginNameSpace( ONEFLOW )

// ---------------------------------------------------------------------
// Recommended modern name
// ---------------------------------------------------------------------
#define Fatal( _Expression ) \
    do { \
        std::ostringstream _oss; \
        _oss << "\n++++++++++++++++++ Fatal Error +++++++++++++++++++++++++++++\n" \
             << (_Expression) << "\n" \
             << " File   : " << __FILE__ << "\n" \
             << " Line   : " << __LINE__ << "\n" \
             << " Compiled on " << __DATE__ << " at " << __TIME__ << "\n" \
             << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n"; \
        throw std::runtime_error( _oss.str() ); \
    } while (0)

// Optional: keep the old function declaration if other places call it directly
void StopProgramFunction( const std::string & stopInformation,
                          const std::string & fileName,
                          const int & fileLine,
                          const std::string & dateName,
                          const std::string & timeName );

EndNameSpace
