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


#include "UnstructuredCgnsZoneWriter.h"
#include "GridTypes.h"
#include "CgnsZsection.h"
#include "CgnsZbc.h"
#include "CgnsZbcBoco.h"
#include "CgnsBcBoco.h"
#include "CgnsSection.h"
#include "CgnsZone.h"
#include "CgnsCoor.h"
#include "NodeMesh.h"
#include "ElementHome.h"

BeginNameSpace( ONEFLOW )

void UnstructuredCgnsZoneWriter::Write(
    UnstructuredMeshData &mesh,
    CgnsZone &cgnsZone,
    bool convertToInnerDataStandard)
{
    // This routine always materializes an unstructured zone. Set the type
    // before creating sections and BC lists so their allocation follows the
    // unstructured CGNS path (the structured path expands two-value ranges).
    cgnsZone.cgnsZoneType = Unstructured;

    const int nNodes = mesh.xN.size();
    const int nCells = mesh.nElem;

    CgnsCoor & cgnsCoor = cgnsZone.RequireCgnsCoor();
    cgnsCoor.SetNNode( nNodes );
    cgnsCoor.SetNCell( nCells );

    NodeMesh & nodeMesh = cgnsCoor.RequireNodeMesh();

    nodeMesh.SetCoordinates( mesh.xN, mesh.yN, mesh.zN );

    cgnsZone.zoneName = "Zone1";
    cgnsZone.isize[ 0 ] = nNodes;
    cgnsZone.isize[ 1 ] = nCells;
    cgnsZone.isize[ 2 ] = 0;

    if ( ! convertToInnerDataStandard )
    {
        cgnsCoor.nCoor = 3;
        cgnsCoor.coorNameList[ 0 ] = "CoordinateX";
        cgnsCoor.coorNameList[ 1 ] = "CoordinateY";
        cgnsCoor.coorNameList[ 2 ] = "CoordinateZ";
        for ( int iCoor = 0; iCoor < 3; ++ iCoor )
        {
            cgnsCoor.typeList[ iCoor ] = RealDouble;
            cgnsCoor.nNodeList[ iCoor ] = nNodes;
            cgnsCoor.Alloc( iCoor, nNodes, RealDouble );
        }
        cgnsCoor.SetAllCoorData();
    }
    
    SecMarkerManager volSec;

    mesh.volSec.CalcVolSec(mesh, volSec);
    SecMarkerManager bcSec;
    mesh.mmark.CalcSecMarker( bcSec );

    const int nVolSec = volSec.nType;
    const int nBcSec = bcSec.nType;

    const int nSection = nVolSec + nBcSec;

    CgnsZsection & cgnsZsection = cgnsZone.RequireCgnsZsection();

    cgnsZsection.CreateCgnsSections( nSection );

    int nVolCell = volSec.CalcTotalElem();
 
    int sumElem = 0;
    for ( int iSection = 0; iSection < nSection; ++ iSection )
    {
        CgnsSection & cgnsSection = cgnsZsection.GetCgnsSection( iSection );
        const SecMarker * sec = nullptr;
        if ( iSection < nVolSec )
        {
            sec = & volSec.data[ iSection ];
        }
        else
        {
            int jSection = iSection - nVolSec;
            sec = & bcSec.data[ jSection ];
        }
            

        const int nElem = sec->nElem;
        cgnsSection.SetSectionInfo( sec->name, sec->cgns_type, sumElem + 1, sumElem + nElem );
        cgnsSection.CreateConnList();
        sumElem += nElem;

        int pos = 0;
        for ( int iElem = 0; iElem < nElem; ++ iElem )
        {
            const IntField & elem = sec->elems[ iElem ];
            const int nNodes = elem.size();
            for ( int i = 0; i < nNodes; ++ i )
            {
                cgnsSection.connList[ pos ++ ]= elem[ i ] + 1;
            }
        }

    }

    for ( int iSection = 0; iSection < nSection; ++ iSection )
    {
        CgnsSection & cgnsSection = cgnsZsection.GetCgnsSection( iSection );
        cgnsSection.SetElemPosition();
    }

    CgnsZbc & cgnsZbc = cgnsZone.RequireCgnsZbc();
    CgnsZbcBoco & cgnsZbcBoco = cgnsZbc.RequireCgnsZbcBoco();
    cgnsZbcBoco.ReadZnboco( mesh.mmark.nMarker );
    cgnsZbcBoco.CreateCgnsZbc();

    for ( int iMarker = 0; iMarker < mesh.mmark.nMarker; ++ iMarker )
    {
        const Marker * marker = & mesh.mmark.markerList[ iMarker ];
        const std::string & name = marker->name;

        CgnsBcBoco & cgnsBcBoco = cgnsZbcBoco.GetCgnsBc( iMarker );
        cgnsBcBoco.name = name;
        // Boundary elements are one topological dimension below the cells.
        // CGNS therefore locates 2D unstructured boundaries at EdgeCenter
        // and 3D unstructured boundaries at FaceCenter.
        cgnsBcBoco.gridLocation = mesh.ndim == 2 ? EdgeCenter :
                                 mesh.ndim == 3 ? FaceCenter : Vertex;
        cgnsBcBoco.nElements    = marker->nElem;
        cgnsBcBoco.bcType = static_cast< BCType_t >( marker->cgns_bcType );
        cgnsBcBoco.pointSetType = PointList;
        cgnsBcBoco.CreateCgnsBcBoco();

        for ( int iElem = 0; iElem < marker->nElem; ++ iElem )
        {
            const int elemId = mesh.mmark.l2g[ iMarker ][ iElem ];
            cgnsBcBoco.SetConnListValue( iElem, elemId + 1 + nVolCell );
        }
        
        //string bcName = GetCgnsBcName( cgnsBcBoco.bcType );
    }

    if ( convertToInnerDataStandard )
    {
        cgnsZone.ConvertToInnerDataStandard();
    }


}

EndNameSpace
