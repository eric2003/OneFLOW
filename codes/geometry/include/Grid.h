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
#include "Constant.h"
#include "HXDefine.h"
#include <vector>
#include <string>
#include <memory>
#include <stdexcept>


BeginNameSpace( ONEFLOW )

#define IMPLEMENT_GRID_CLONE( TYPE ) \
std::unique_ptr< Grid > Clone() const override { return std::make_unique< TYPE >(); }

#define REGISTER_GRID( TYPE ) \
    Grid * TYPE ## _myClass = \
        Grid::Register( #TYPE, std::make_unique< TYPE >() );

class DataBook;
class NodeMesh;
class InterFace;
class SlipFace;
class DataBase;
class IFaceLink;

class Grid
{
public:
    Grid();
    virtual ~Grid();
public:
    virtual std::unique_ptr< Grid > Clone() const = 0;
public:
    // Preferred: exclusive ownership of a registered grid prototype clone.
    static std::unique_ptr< Grid > SafeCloneUnique( const std::string & type );
    static Grid * Register( const std::string & type, std::unique_ptr< Grid > clone );
public:
    std::string name;
    int dimension;
    int type, level;
    int id, localId;
    int nNodes;
    int nFaces, nCells;
    int nBFaces;
    int nIFaces;
    int volBcType;
    std::unique_ptr< NodeMesh > nodeMesh;
    std::unique_ptr< InterFace > interFace;
    std::unique_ptr< SlipFace > slipFace;
    std::unique_ptr< DataBase > dataBase;
public:
    DataBase & GetDataBase()
    {
        if ( ! dataBase )
        {
            throw std::logic_error( "Grid DataBase is not initialized" );
        }
        return *dataBase;
    }

    const DataBase & GetDataBase() const
    {
        if ( ! dataBase )
        {
            throw std::logic_error( "Grid DataBase is not initialized" );
        }
        return *dataBase;
    }

    [[nodiscard]] DataBase * TryGetDataBase() noexcept { return dataBase.get(); }
    [[nodiscard]] const DataBase * TryGetDataBase() const noexcept { return dataBase.get(); }
public:
    void BasicInit();
    void Free();
    virtual void Init();
public:
    bool IsOneD();
    bool IsTwoD();
    bool IsThreeD();
public:
    virtual void ReadGrid ( std::fstream & file ) {};
    virtual void WriteGrid( std::fstream & file ) {};
    virtual void Decode( DataBook * databook ){};
    virtual void Encode( DataBook * databook ){};
    virtual void ReadGrid( DataBook * databook ){};
    virtual void WriteGrid( DataBook * databook ){};
    virtual void ModifyBcType( int bcType1, int bcType2 ) {};
    virtual void GenerateLgMapping( IFaceLink * iFaceLink ){};
    virtual void ReGenerateLgMapping( IFaceLink * iFaceLink ){};
    virtual void UpdateOtherTopologyTerm( IFaceLink * iFaceLink ){};
public:
    virtual void GetMinMaxDistance( Real & dismin, Real & dismax ) {};
    virtual void CalcMetrics() {};
};

EndNameSpace
