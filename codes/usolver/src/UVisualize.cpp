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

#include "UVisualize.h"
#include "NodeField.h"
#include "HXMath.h"
#include "DataBook.h"
#include "DataBase.h"
#include "UnsGrid.h"
#include "FaceTopo.h"
#include "BcRecord.h"
#include "Boundary.h"
#include "ActionState.h"
#include "Dimension.h"
#include "NsIdx.h"
#include "Zone.h"
#include "ZoneState.h"
#include "StringUtils.h"
#include "Prj.h"
#include "HXLookup.h"
#include "NodeMesh.h"
#include "NsCtrl.h"
#include <sstream>
#include <iostream>
#include <algorithm>
#include <iomanip>
#include <unordered_map>

BeginNameSpace( ONEFLOW )

VisualTool::VisualTool()
{
    ;
}

VisualTool::~VisualTool()
{
    int nSize = qNodeField.size();
    for ( int i = 0; i < nSize; ++ i )
    {
        delete qNodeField[ i ];
    }
}

void VisualTool::Init()
{
    title.push_back( "title=\"THE FLOW FIELD OF ONEFLOW\"" );
    title.push_back( "variables=" );
    title.push_back( "\"x\"" );
    title.push_back( "\"y\"" );
    title.push_back( "\"z\"" );
}

void VisualTool::AddTitle( const std::string & varName )
{
    title.push_back( AddString( "\"",  varName, "\"" ) );
}

MRField * VisualTool::AddField( const std::string & varName )
{
    this->AddTitle( varName );
    MRField * fn = InterpolateCellToNode( varName );
    qNodeField.push_back( fn );
    return fn;
}

MRField * VisualTool::AddField( RealField & qc, const std::string & varName )
{
    this->AddTitle( varName );
    MRField * fn = InterpolateCellToNode( qc );
    qNodeField.push_back( fn );
    return fn;
}

MRField * VisualTool::CreateField( const std::string & varName, int nEqu )
{
    this->AddTitle( varName );
    MRField * fn = AllocateNodeField( nEqu );
    qNodeField.push_back( fn );
    return fn;
}


BcVisual::BcVisual()
{
    ;
}

BcVisual::~BcVisual()
{
    ;
}

void BcVisual::Calc( int bcType )
{
    UnsGrid * grid = Zone::GetUnsGrid();

    this->Calcf2n( bcType );

    e2n.resize( 0 );
    lcell.resize( 0 );
    rcell.resize( 0 );

    ResolveElementEdge();
}

void BcVisual::ResolveElementEdge()
{
    int nFaces = this->f2n.size();
    int nSize = 2;
    HXLookup<int> faceLookup;

    for ( int fId = 0; fId < nFaces; ++ fId )
    {
        int nNodes = this->f2n[ fId ].size();
        for ( int iNode = 0; iNode < nNodes; ++ iNode )
        {
            int iNode0 = iNode;
            int iNode1 = ( iNode + 1 ) % nNodes;

            int ip1 = this->f2n[ fId ][ iNode0 ];
            int ip2 = this->f2n[ fId ][ iNode1 ];

            if ( ip1 == ip2 ) continue;

            IntField edgeNodeId;
            edgeNodeId.push_back( ip1 );
            edgeNodeId.push_back( ip2 );

            auto [ edgeIndex, isNew ] = faceLookup.FindOrAdd( edgeNodeId );

            if ( isNew )
            {
                this->lcell.push_back( fId );
                this->rcell.push_back( -1 );      
                this->e2n.push_back( edgeNodeId );
            }
            else
            {
                int ip = this->rcell[ edgeIndex ];
                if ( ip != -1 )
                {
                    std::cout << "Fatal Error\n";
                    std::cout << " edgeIndex = " << edgeIndex << "\n";
                    for ( int i = 0; i < this->e2n[ edgeIndex ].size(); ++ i )
                    {
                        std::cout << this->e2n[ edgeIndex ][ i ] << " ";
                    }
                }

                this->rcell[ edgeIndex ] = fId;
            }
        }
    }
}

void BcVisual::Calcf2n( int bcType )
{
    UnsGrid * grid = Zone::GetUnsGrid();
    FaceTopo * faceTopo = grid->faceTopo;
    LinkField & total_f2n = faceTopo->faces;
    BcRecord * bcRecord = faceTopo->bcManager->bcRecord;

    // 清空数据
    this->f2n.clear();
    this->l2g.clear();

    HXLookup<int> nodeLookup;

    std::unordered_map<int, int> globalNodeToIndex;  // 全局节点编号 -> 局部索引
    globalNodeToIndex.reserve(grid->nBFaces * 4);     // 预分配空间

    int nBFaces = grid->nBFaces;

    for ( int iFace = 0; iFace < nBFaces; ++ iFace )
    {
        if ( bcType != bcRecord->bcType[ iFace ] ) continue;

        const auto& faceNodes = total_f2n[iFace];
        int nNodes = faceNodes.size();

        IntField localf2n;
        localf2n.reserve(nNodes);

        for (int globalNodeId : faceNodes)
        {
            auto it = globalNodeToIndex.find(globalNodeId);
            if (it == globalNodeToIndex.end())
            {
                // 新节点：添加到全局池
                int localIndex = this->l2g.size();
                globalNodeToIndex[globalNodeId] = localIndex;
                this->l2g.push_back(globalNodeId);
                localf2n.push_back(localIndex);
            }
            else
            {
                // 已存在的节点：使用已有索引
                localf2n.push_back(it->second);
            }
        }

        this->f2n.push_back(std::move(localf2n));
    }
}

void BcVisual::Dump( std::ostringstream & oss, VisualTool * visualTool, std::string & bcTitle )
{
    UnsGrid * grid = Zone::GetUnsGrid();

    for ( HXSize_t i = 0; i < visualTool->title.size(); ++ i )
    {
        oss << visualTool->title[ i ] << std::endl;
    }
        
    int nNodes = l2g.size();
    int nFaces = e2n.size();
    int nElem = this->f2n.size();

    // output for Tecplot
    oss << "ZONE\n";

    oss << "T = " << bcTitle << std::endl;

    oss << "ZoneType = FEPolygon\n";

    oss << "Nodes    = " << nNodes << std::endl;
    oss << "Faces    = " << nFaces << std::endl;  
    oss << "Elements = " << nElem << std::endl;  
    oss << "NumConnectedBoundaryFaces = 0\n";
    oss << "TotalNumBoundaryConnections = 0\n";

    Plot::oss = & oss;
    Plot::DumpField( l2g, grid->nodeMesh->xN );
    Plot::DumpField( l2g, grid->nodeMesh->yN );
    Plot::DumpField( l2g, grid->nodeMesh->zN );

    int nVar = visualTool->qNodeField.size();        
    for ( int iVar = 0; iVar < nVar; ++ iVar )
    {
        RealField & q = ( * visualTool->qNodeField[ iVar ] )[ 0 ];
        Plot::DumpField( l2g, q );
    }

    Plot::DumpFaceNodeLink( e2n );
    Plot::DumpFaceElementLink( lcell, nElem );
    Plot::DumpFaceElementLink( rcell, nElem );
       
}

void BcVisual::DumpDebug( std::ostringstream & oss, VisualTool * visualTool, std::string & bcTitle )
{
    UnsGrid * grid = Zone::GetUnsGrid();

    int nNodes = l2g.size();
    int nFaces = e2n.size();
    int nElem = this->f2n.size();

    // output for Tecplot

    oss << " VARIALBES = ";
    oss << " \"x\" ";
    oss << " \"y\" ";
    oss << " \"z\" ";
    oss << "\n";
    //oss << " ZONE N = " << nNodes << " E = " << nElem << " F = FEPOINT, ET = TRIANGLE \n";
    //oss << " ZONE N = " << nNodes << " E = " << nElem << " F = FEPOINT, ET = QUADRILATERAL \n";
    oss << "title = \"THE FLOW FIELD OF ONEFLOW\" \n";
    oss << "VARIALBES = \"x\" \"y\" \"z\" " << "\n";

    oss << " ZONE N = " << nNodes << " E = " << nElem << " F = FEPOINT, ET = QUADRILATERAL \n";

    for ( int iNode = 0; iNode < nNodes; ++ iNode )
    {
        int id = l2g[ iNode ];
        oss << grid->nodeMesh->xN[ id ] << " ";
        oss << grid->nodeMesh->yN[ id ] << " ";
        oss << grid->nodeMesh->zN[ id ] << " ";
        oss << "\n";
    }

    for ( int iElem = 0; iElem < nElem; ++ iElem )
    {
        int p1 = this->f2n[ iElem ][ 0 ] + 1;
        int p2 = this->f2n[ iElem ][ 1 ] + 1;
        int p3 = this->f2n[ iElem ][ 2 ] + 1;
        int p4 = this->f2n[ iElem ][ 3 ] + 1;

        oss << p1 << " ";
        oss << p2 << " ";
        oss << p3 << " ";
        oss << p4 << " ";
        oss << "\n";
    }

    //DumpSeveralElement();
}

void BcVisual::DumpSeveralElement()
{
    UnsGrid * grid = Zone::GetUnsGrid();

    int nNodes = l2g.size();
    int nFaces = e2n.size();
    int nElem = this->f2n.size();

    std::fstream file;
    std::string fileName = "test.plt";
    Prj::OpenPrjFile( file, fileName, std::ios_base::out );

    file << " VARIALBES = ";
    file << " \"x\" ";
    file << " \"y\" ";
    file << " \"z\" ";
    file << "\n";

    IntField eList = {1360, 1361};
    IntField nList, nList1, nList2;

    // 使用 unordered_map 做节点去重和索引映射
    std::unordered_map<int, int> nodeToIndex;  // 节点编号 -> 去重后的索引
    IntField uniqueNodes;                      // 去重后的节点列表（按索引顺序）
    IntField nodeIndices;                      // 每个原始节点的去重索引
    IntField localNodeIds;                     // 原始节点编号列表

    int iCount = 0;
    for ( int e : eList )
    {
        const auto & faceNodes = this->f2n[ e ];
        for ( int nodeId : faceNodes )
        {
            localNodeIds.push_back( nodeId );
            auto it = nodeToIndex.find( nodeId );
            if ( it == nodeToIndex.end() )
            {
                // 新节点：分配索引
                int newIndex = uniqueNodes.size();
                nodeToIndex[ nodeId ] = newIndex;
                uniqueNodes.push_back( nodeId );
                nodeIndices.push_back( newIndex );
            }
            else
            {
                // 已存在的节点：使用已有索引
                nodeIndices.push_back( it->second );
            }
        }
    }
 
    file << " ZONE N = " << iCount << " E = " << eList.size() << " F = FEPOINT, ET = QUADRILATERAL \n";
    int width = 20;
    int pre = 20;

    for (int nodeId : uniqueNodes)
    {
        int globalId = l2g[nodeId];
        file << std::setw(width) << std::setprecision(pre) << grid->nodeMesh->xN[globalId] << " ";
        file << std::setw(width) << std::setprecision(pre) << grid->nodeMesh->yN[globalId] << " ";
        file << std::setw(width) << std::setprecision(pre) << grid->nodeMesh->zN[globalId] << " ";
        file << "\n";
    }

    // 输出节点编号（调试信息）
    for (int nodeId : localNodeIds)
    {
        int globalId = l2g[nodeId];
        file << globalId << "\n";
    }

    // 输出单元连接
    int pos = 0;
    int nodesPerFace = 4;  // QUADRILATERAL
    for (int e = 0; e < eList.size(); ++e)
    {
        for (int i = 0; i < nodesPerFace; ++i)
        {
            file << nodeIndices[pos + i] + 1 << " ";
        }
        file << "\n";
        pos += nodesPerFace;
    }
}

UVisualize::UVisualize()
{
    ;
}

UVisualize::~UVisualize()
{
    ;
}

bool UVisualize::NeedVisualField()
{
    bool flag1 = Dim::dimension == TWO_D;
    bool flag2 = Dim::dimension == THREE_D && ctrl.showfield == 1;
    return  flag1 || flag2;
}

void UVisualize::Visual()
{
    VisualTool visualTool;
    visualTool.Init();

    this->CalcNodeField( & visualTool );

    std::ostringstream oss;

    this->ShowBc( oss, & visualTool );
    //this->ShowBcDebug( oss, & visualTool );
    //this->ShowBcDebugTest( oss, & visualTool );

    if ( this->NeedVisualField() )
    {
        this->ShowField( oss, & visualTool );
    }


    ToDataBook( ActionState::dataBook, oss );
}

void UVisualize::ExtractLinkNum( LinkField & f2n, IntField & fnNumber )
{
    int nSize = f2n.size();
    fnNumber.resize( nSize );
    for ( int i = 0; i < nSize; ++ i )
    {
        fnNumber[ i ] = f2n[ i ].size();
    }
}

int UVisualize::GetTotalNumFaceNodes( LinkField & f2n )
{
    int totalNumFaceNodes = 0;
    int nSize = f2n.size();
    for ( int i = 0; i < nSize; ++ i )
    {
        totalNumFaceNodes += f2n[ i ].size();
    }
    return totalNumFaceNodes;
}

void UVisualize::ShowField( std::ostringstream & oss, VisualTool * visualTool )
{
    UnsGrid * grid = Zone::GetUnsGrid();

    FaceTopo * faceTopo = grid->faceTopo;
    LinkField & f2n = faceTopo->faces;

    int nNodes = grid->nNodes;
    int nCells = grid->nCells;
    int nFaces = grid->nFaces;

    for ( HXSize_t i = 0; i < visualTool->title.size(); ++ i )
    {
        oss << visualTool->title[ i ] << std::endl;
    }

    int totalNumFaceNodes = this->GetTotalNumFaceNodes( f2n );

    // output for Tecplot
    oss << "ZONE\n";
    if ( Dim::dimension == THREE_D )
    {
        oss << "ZoneType = FEPolyhedron\n";
    }
    else
    {
        oss << "ZoneType = FEPolygon\n";
    }
    oss << "Nodes    = " << nNodes << std::endl;
    oss << "Faces    = " << nFaces << std::endl;
    oss << "Elements = " << nCells << std::endl;
    oss << "TotalNumFaceNodes = " << totalNumFaceNodes << std::endl;
    oss << "NumConnectedBoundaryFaces = 0\n";
    oss << "TotalNumBoundaryConnections = 0\n";
    Plot::oss = & oss;
    Plot::DumpField( grid->nodeMesh->xN );
    Plot::DumpField( grid->nodeMesh->yN );
    Plot::DumpField( grid->nodeMesh->zN );
    
    int nVar = visualTool->qNodeField.size();        
    for ( int iVar = 0; iVar < nVar; ++ iVar )
    {
        RealField & q = ( * visualTool->qNodeField[ iVar ] )[ 0 ];
        Plot::DumpField( q );
    }

    if ( Dim::dimension == THREE_D )
    {
        Plot::DumpFaceNodeNumber( f2n );
    }

    Plot::DumpFaceNodeLink( f2n );
    Plot::DumpFaceElementLink( faceTopo->lCells, nCells );
    Plot::DumpFaceElementLink( faceTopo->rCells, nCells );
}

void UVisualize::ShowBc( std::ostringstream & oss, VisualTool * visualTool )
{
    if ( IsTwoD() ) return;
    UnsGrid * grid = Zone::GetUnsGrid();

    IntField bcTypeList;
    grid->faceTopo->bcManager->CalcBcType( bcTypeList );
    int nBcType = bcTypeList.size();

    for ( int iBcType = 0; iBcType < nBcType; ++ iBcType )
    {
        int bcType = bcTypeList[ iBcType ];

        if ( BC::IsInterfaceBc( bcType )  ) continue;

        std::string bcTitle = AddString( "\"",  ZoneState::zid, "BC=", bcType, "\"" );

        BcVisual bcVisual;

        bcVisual.Calc( bcType );

        bcVisual.Dump( oss, visualTool, bcTitle );
    }    
}

void UVisualize::ShowBcDebugTest( std::ostringstream & oss, VisualTool * visualTool )
{
    if ( IsTwoD() ) return;
    UnsGrid * grid = Zone::GetUnsGrid();

    IntField bcTypeList;
    grid->faceTopo->bcManager->CalcBcType( bcTypeList );
    int nBcType = bcTypeList.size();

    for ( int iBcType = 0; iBcType < nBcType; ++ iBcType )
    {
        int bcType = bcTypeList[ iBcType ];

        if ( BC::IsInterfaceBc( bcType ) ) continue;
        if ( bcType != BC::SYMMETRY ) continue;

        std::string bcTitle = AddString( "\"", ZoneState::zid, "BC=", bcType, "\"" );

        BcVisual bcVisual;

        bcVisual.Calc( bcType );

        bcVisual.DumpDebug( oss, visualTool, bcTitle );
    }
}

void UVisualize::CalcNodeField( VisualTool * visualTool )
{
    UnsGrid * grid = Zone::GetUnsGrid();
    MRField * q = ONEFLOW::GetFieldPointer< MRField >( grid, "q" );

    MRField * rn = visualTool->AddField( ( * q )[ IDX::IR ], "r" );
    MRField * un = visualTool->AddField( ( * q )[ IDX::IU ], "u" );
    MRField * vn = visualTool->AddField( ( * q )[ IDX::IV ], "v" );
    MRField * wn = visualTool->AddField( ( * q )[ IDX::IW ], "w" );
    MRField * pn = visualTool->AddField( ( * q )[ IDX::IP ], "p" );

    MRField * gaman = InterpolateCellToNode( "gama" );
    MRField * machn = visualTool->CreateField( "mach" );
    CalcMach( rn, un, vn, wn, pn, gaman, machn );
    delete gaman;

    MRField * tempr = ONEFLOW::GetFieldPointer< MRField >( grid, "tempr" );
    visualTool->AddField( ( * tempr )[ IDX::ITT ], "tempr" );

    if ( vis_model.vismodel > 0 )
    {
        visualTool->AddField( "visl" );
        visualTool->AddField( "vist" );
    }
}

void CalcMach( MRField * r, MRField * u, MRField * v, MRField * w, MRField * p, MRField * gama, MRField * mach )
{
    UnsGrid * grid = Zone::GetUnsGrid();
    int nNodes = grid->nNodes;
    for ( int iNode = 0; iNode < nNodes; ++ iNode )
    {
        Real rm = ( * r )[ 0 ][ iNode ];
        Real um = ( * u )[ 0 ][ iNode ];
        Real vm = ( * v )[ 0 ][ iNode ];
        Real wm = ( * w )[ 0 ][ iNode ];
        Real pm = ( * p )[ 0 ][ iNode ];

        Real gm = ( * gama )[ 0 ][ iNode ];
        Real v2 = SQR( um, vm, wm );
        Real c2 = gm * pm / rm;
        Real mm = sqrt( v2 / c2 );
        ( * mach )[ 0 ][ iNode ] = mm;
    }
}

EndNameSpace
