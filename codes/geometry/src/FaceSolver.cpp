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

#include "FaceSolver.h"
#include "ElementHome.h"
#include "Stop.h"
#include "FaceTopo.h"
#include <iostream>
#include <algorithm>
using namespace std;

BeginNameSpace( ONEFLOW )

FaceSolver::FaceSolver()
{
    //this->nodeId = new LinkField();
    //this->lCell = new IntField();
    //this->rCell = new IntField();
    this->faceBcKey = new IntField();
    this->faceBcType = new IntField();
    //this->faceType = new IntField();
    this->childFid = new LinkField();

    this->refFaces = new set< Mid<int> >();
    this->faceTopo = new FaceTopo();
}

FaceSolver::~FaceSolver()
{
    //delete this->nodeId;
    //delete this->lCell;
    //delete this->rCell;
    delete this->faceBcKey;
    delete this->faceBcType;
    //delete this->faceType;
    delete this->childFid;

    delete this->refFaces;
    delete this->faceTopo;
}

int FaceSolver::FindFace( Mid<int> & face )
{
    set< Mid<int> >::iterator iter = this->refFaces->find( face );
    if ( iter == this->refFaces->end() )
    {
        return ONEFLOW::INVALID_INDEX;
    }
    return iter->id;
}

bool FaceSolver::CheckBcFace( IntSet & bcVertex, IntField & nodeId )
{
    int size = nodeId.size();
    for ( int iNode = 0; iNode < size; ++ iNode )
    {
        IntSet::iterator iter = bcVertex.find( nodeId[ iNode ] );
        if ( iter == bcVertex.end() )
        {
            return false;
        }
    }
    return true;
}


void FaceSolver::ScanElementFace( CgIntField & eNodeId, int eType, int eId )
{
    UnitElement * unitElement = ElementHome::GetUnitElement( eType );

    //composite Element not to be involved in analysis !!!
    int nElemFace = unitElement->faceList.size();
    for ( int iFace = 0; iFace < nElemFace; ++ iFace )
    {
        IntField & rNodeId = unitElement->faceList[ iFace ];
        int fType = unitElement->GetFaceType( iFace );
         
        int nNode = rNodeId.size();

        IntField aNodeId;
        for ( int iNode = 0; iNode < nNode; ++ iNode )
        {
            aNodeId.push_back( eNodeId[ rNodeId[ iNode ] ] );
        }                                                              

        Mid<int> fMid( nNode, this->faceTopo->f2n.size() );
        fMid.data = aNodeId;
        std::sort( fMid.data.begin(), fMid.data.end() );
        int gFid = this->FindFace( fMid );

        if ( gFid == ONEFLOW::INVALID_INDEX )
        {
            int totalfn    = this->refFaces->size();
            int faceNumber = this->faceTopo->lCell.size();
            if ( totalfn != faceNumber )
            {
                cout << "totalfn != faceNumber " << totalfn << " " << faceNumber << endl;
                Stop("");
            }
            this->refFaces->insert( fMid );

            this->faceTopo->lCell.push_back( eId );
            this->faceTopo->rCell.push_back( ONEFLOW::INVALID_INDEX );

            this->faceBcType->push_back( ONEFLOW::INVALID_INDEX );
            this->faceBcKey->push_back( ONEFLOW::INVALID_INDEX );
            this->faceTopo->faceType.push_back( fType );

            this->faceTopo->f2n.push_back( aNodeId );
            this->childFid->resize( totalfn + 1 );
        }
        else
        {
            if ((this->faceTopo->lCell)[gFid] == ONEFLOW::INVALID_INDEX)
            {
                //说明此面虽然存在，但是由于种种原因没有处理
                (this->faceTopo->lCell)[gFid] = eId; //例如分裂过程中新加入体单元面
            }
            else
            {
                if ( (this->faceTopo->rCell)[gFid] == ONEFLOW::INVALID_INDEX )
                {
                    if ((this->faceTopo->lCell)[gFid] != eId)
                    {
                        (this->faceTopo->rCell)[gFid] = eId;
                    }
                }
            }
        }
    }
}

void FaceSolver::ScanBcFace( IntSet& bcVertex, int bcType, int bcNameId )
{
    int nBFace = 0;

    //cout << " this->faceTopo = " << this->faceTopo << "\n";
    int nFace = this->faceTopo->lCell.size();

    cout << " nFace = " << nFace << "\n";
    int nTraditionalBc = 0;
    for ( int iFace = 0; iFace < nFace; ++ iFace )
    {
        int rCell = ( this->faceTopo->rCell )[ iFace ];

        if ( rCell == ONEFLOW::INVALID_INDEX )
        {
            ++ nTraditionalBc;
        }
    }
    cout << " nTraditionalBc = " << nTraditionalBc << "\n";


    for ( int iFace = 0; iFace < nFace; ++ iFace )
    {
        if ( iFace % 200000 == 0 ) 
        {
            //cout << " iFace = " << iFace << " numberOfTotalFaces = " << nFace << endl;
        }
        int originalBcType = ( * this->faceBcType )[ iFace ];
        int rCell     = ( this->faceTopo->rCell )[ iFace ];

        if ( ( rCell          == ONEFLOW::INVALID_INDEX ) && 
             ( originalBcType == ONEFLOW::INVALID_INDEX ) )
        {
            if ( this->CheckBcFace( bcVertex, ( this->faceTopo->f2n )[ iFace ] ) )
            {
                ++ nBFace;

                ( * this->faceBcType )[ iFace ] = bcType;
                ( * this->faceBcKey  )[ iFace ] = bcNameId;
            }
        }
    }

    //cout << " nBFace = " << nBFace << endl;
    int kkk = 1;
}

void FaceSolver::ScanBcFaceDetail( IntSet& bcVertex, int bcType, int bcNameId )
{
    int nFace = this->faceTopo->lCell.size();
    cout << " nFace = " << nFace << "\n";

    int nTraditionalBc = 0;
    for ( int iFace = 0; iFace < nFace; ++ iFace )
    {
        int rCell = ( this->faceTopo->rCell )[ iFace ];

        if ( rCell == ONEFLOW::INVALID_INDEX )
        {
            ++ nTraditionalBc;
        }
    }
    cout << " nTraditionalBc = " << nTraditionalBc << "\n";

    int nBFace = 0;
    for ( int iFace = 0; iFace < nFace; ++ iFace )
    {
        if ( iFace % 200000 == 0 ) 
        {
            //cout << " iFace = " << iFace << " numberOfTotalFaces = " << nFace << endl;
        }
        int originalBcType = ( * this->faceBcType )[ iFace ];

        if ( originalBcType == ONEFLOW::INVALID_INDEX )
        {
            if ( this->CheckBcFace( bcVertex, ( this->faceTopo->f2n )[ iFace ] ) )
            {
                ++ nBFace;

                ( * this->faceBcType )[ iFace ] = bcType;
                ( * this->faceBcKey  )[ iFace ] = bcNameId;
            }
        }
    }

    cout << " nFinalBcFace = " << nBFace << " bcType = " << bcType << endl;
    int kkk = 1;
}

void FaceSolver::ScanInterfaceBc()
{
    int nFace = this->faceTopo->lCell.size();

    int bcNameId = -1;
    int nInterFace = 0;
    for ( int iFace = 0; iFace < nFace; ++ iFace )
    {
        int originalBcType = ( * this->faceBcType )[ iFace ];

        if ( originalBcType == ONEFLOW::INVALID_INDEX )
        {
            nInterFace ++;
            ( * this->faceBcType )[ iFace ] = BCTypeNull;
            ( * this->faceBcKey  )[ iFace ] = bcNameId;
        }
    }

    cout << " nInterFace = " << nInterFace << endl;
}

int FaceSolver::GetNSimpleFace()
{
    int nSimpleFace = 0;

    for ( int iFace = 0; iFace < this->faceTopo->f2n.size(); ++ iFace )
    {
        int nCFace = ( * this->childFid )[ iFace ].size();
        if ( nCFace == 0 )
        {
            ++ nSimpleFace;
        }
    }
    return nSimpleFace;
}


EndNameSpace