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
#include "UnstructuredMeshData.h"
#include <string>
#include <map>
#include <set>

BeginNameSpace( ONEFLOW )

class GridMediator;
struct GridConfig;
class CgnsZone;

class Su2Bc
{
public:
    Su2Bc();
    ~Su2Bc();
public:
    std::set<std::string> bcList;
    std::map<std::string, std::string> bcMap;
    std::map<std::string, int> bcNameToValueMap;
    void Init();
    void AddBc(const std::string &geoName, const std::string &bcName);
    void Process(const StringField &markerBCNameList, const StringField &markerNameList);
    std::string GetBcName(const std::string &geoName) const;
    int GetCgnsBcType(const std::string &geoName) const;
};

// SU2-specific reader and compatibility adapter over the shared unstructured
// value model. CGNS zone population belongs to the common writer.
class Su2Grid : public UnstructuredMeshData
{
public:
    Su2Grid();
    ~Su2Grid() = default;
public:
    void ReadSu2Grid(GridMediator *gridMediator);
    void ReadSu2GridAscii(const std::string &fileName, const std::string &caseDir);
    void Su2ToOneFlowGrid(const GridConfig &config, const std::string &caseDir);
    void MarkBoundary(const std::string &su2cfgFile, const std::string &caseDir);
    void FillSU2CgnsZone(CgnsZone &cgnsZone, bool convertToInnerDataStandard = true);
public:
    IntField vtkmap;
    IntField elemId;
    Su2Bc su2Bc;
    int nZone{1};
};

void Su2ToOneFlowGrid(Su2Grid &su2Grid);

EndNameSpace
