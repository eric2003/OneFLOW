/*---------------------------------------------------------------------------*\
    OneFLOW - LargeScale Multiphysics Scientific Simulation Environment
    Copyright (C) 2017-2024 He Xin and the OneFLOW contributors.
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
#include "HXDefine.h"
#include "HXCgns.h"
#include <string>

BeginNameSpace( ONEFLOW )

#ifdef ENABLE_CGNS

class CgnsZone;

class CgnsBcLink
{
public:
    CgnsBcLink( CgnsZone * cgnsZone );
    ~CgnsBcLink();
public:
    void ConvertToInnerDataStandard();
public:
    CgInt          nConnPoints;
    CgIntField     connPoint;

    CgInt          nConnDonorPoints;
    CgIntField     connDonorPoint;

    ZoneType_t     donorZoneType;
    PointSetType_t donorPointSetType;
    DataType_t     donorDataType;

    std::string connName;
    std::string donorZoneName;

    CgnsZone * cgnsZone;
    int bcId;
};

#endif

EndNameSpace
