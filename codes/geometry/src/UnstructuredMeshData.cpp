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

#include "UnstructuredMeshData.h"
#include "GridTypes.h"
#include "HXStd.h"
#include "ElementHome.h"

BeginNameSpace( ONEFLOW )

UnstructuredMeshData::UnstructuredMeshData()
    : ndim(0), nPoin(0), nElem(0)
{
}

VTK_CgnsMap vtk_CgnsMap;

int VTK_TYPE::VERTEX        = 1;
int VTK_TYPE::LINE          = 3;
int VTK_TYPE::TRIANGLE      = 5;
int VTK_TYPE::QUADRILATERAL = 9;
int VTK_TYPE::TETRAHEDRON   = 10;
int VTK_TYPE::HEXAHEDRON    = 12;
int VTK_TYPE::PRISM         = 13;
int VTK_TYPE::PYRAMID       = 14;

VTK_CgnsMap::VTK_CgnsMap()
{
    this->Init();
}

VTK_CgnsMap::~VTK_CgnsMap()
{
    ;
}

void VTK_CgnsMap::Init()
{
    using IntPair = std::pair< int, int >;

    vtk2Cgns.insert( IntPair( VTK_TYPE::VERTEX       , NODE    ) );
    vtk2Cgns.insert( IntPair( VTK_TYPE::LINE         , BAR_2   ) );
    vtk2Cgns.insert( IntPair( VTK_TYPE::TRIANGLE     , TRI_3   ) );
    vtk2Cgns.insert( IntPair( VTK_TYPE::QUADRILATERAL, QUAD_4  ) );
    vtk2Cgns.insert( IntPair( VTK_TYPE::TETRAHEDRON  , TETRA_4 ) );
    vtk2Cgns.insert( IntPair( VTK_TYPE::HEXAHEDRON   , HEXA_8  ) );
    vtk2Cgns.insert( IntPair( VTK_TYPE::PRISM        , PENTA_6 ) );
    vtk2Cgns.insert( IntPair( VTK_TYPE::PYRAMID      , PYRA_5 ) );
}

void SecMarkerManager::Alloc( int nType )
{
    this->nType = nType;
    data.clear();
    data.resize( nType );
}

int SecMarkerManager::CalcTotalElem() const
{
    const int nType = data.size();
    int nElem = 0;
    for ( int i = 0; i < nType; ++ i )
    {
        const SecMarker * sec = & this->data[ i ];
        nElem += sec->nElem;
    }
    return nElem;
}

void MarkerManager::CreateMarkerList( int nMarker )
{
    markerList.clear();
    markerList.resize( nMarker );
}

void MarkerManager::CalcSecMarker( SecMarkerManager & secMarkerManager )
{
    const int nMarker = this->markerList.size();
    IntSet typeSet;
    for ( int iMarker = 0; iMarker < nMarker; ++ iMarker )
    {
        const Marker * marker = & this->markerList[ iMarker ];
        for ( int iElem = 0; iElem < marker->nElem; ++ iElem )
        {
            const int type = marker->eTypes[ iElem ];
            typeSet.insert( type );
        }
    }

    ONEFLOW::Set2Array( typeSet, types );

    l2g.resize( nMarker );
    for ( int iMarker = 0; iMarker < nMarker; ++ iMarker )
    {
        const Marker * marker = & this->markerList[ iMarker ];
        l2g[ iMarker ].resize( marker->nElem );
    }

    const int nType = types.size();

    secMarkerManager.Alloc( nType );

    int gid = 0;
    for ( int iType = 0; iType < nType; ++ iType )
    {
        const int eType = types[ iType ];
        SecMarker * secMarker = & secMarkerManager.data[ iType ];
        secMarker->vtk_type = eType;
        secMarker->cgns_type = vtk_CgnsMap.vtk2Cgns[ eType ];
        secMarker->name = ElementTypeName[ secMarker->cgns_type ];
        for ( int iMarker = 0; iMarker < nMarker; ++ iMarker )
        {
            const Marker * marker = & this->markerList[ iMarker ];
            for ( int iElem = 0; iElem < marker->nElem; ++ iElem )
            {
                const int type = marker->eTypes[ iElem ];
                if ( type == eType )
                {
                    secMarker->elems.push_back( marker->elems[ iElem ] );
                    this->l2g[ iMarker ][ iElem ] = gid ++;
                }
            }
        }
        secMarker->nElem = secMarker->elems.size();
    }

}

void VolumeSecManager::CalcVolSec( const UnstructuredMeshData & mesh, SecMarkerManager & secMarkerManager )
{
    IntSet typeSet;

    for ( int iElem = 0; iElem < mesh.nElem; ++ iElem )
    {
        const int eVtkType = mesh.elemVTKType[ iElem ];
        typeSet.insert( eVtkType );
    }
    
    ONEFLOW::Set2Array( typeSet, types );

    int nType = types.size();
    l2g.resize( nType );

    secMarkerManager.Alloc( nType );

    for ( int iType = 0; iType < nType; ++ iType )
    {
        int eType = types[ iType ];
        SecMarker * secMarker = & secMarkerManager.data[ iType ];
        secMarker->vtk_type = eType;
        secMarker->cgns_type = vtk_CgnsMap.vtk2Cgns[ eType ];
        secMarker->name = ElementTypeName[ secMarker->cgns_type ];
        int gid = 0;
        for ( int iElem = 0; iElem < mesh.nElem; ++ iElem )
        {
            const int e_VtkType = mesh.elemVTKType[ iElem ];
            if ( eType == e_VtkType )
            {
                secMarker->elems.push_back( mesh.elems[ iElem ] );
                this->l2g[ iType ].push_back( iElem );
            }
        }
        secMarker->nElem = secMarker->elems.size();
    }

}


EndNameSpace
