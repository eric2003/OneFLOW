/*---------------------------------------------------------------------------*\
    OneFLOW - LargeScale Multiphysics Scientific Simulation Environment
    Copyright (C) 2017-2019 He Xin and the OneFLOW contributors.
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
#include "StrUtil.h"
#include "Prj.h"
#include "Mid.h"
#include "NodeMesh.h"
#include "NsCtrl.h"
#include <sstream>
#include <iostream>
#include <algorithm>
#include <iomanip>
using namespace std;


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

void VisualTool::AddTitle( const string & varName )
{
    title.push_back( AddString( "\"",  varName, "\"" ) );
}

MRField * VisualTool::AddField( const string & varName )
{
    this->AddTitle( varName );
    MRField * fn = CreateNodeVar( varName );
    qNodeField.push_back( fn );
    return fn;
}

MRField * VisualTool::AddField( RealField & qc, const string & varName )
{
    this->AddTitle( varName );
    MRField * fn = CreateNodeVar( qc );
    qNodeField.push_back( fn );
    return fn;
}

MRField * VisualTool::CreateField( const string & varName, int nEqu )
{
    this->AddTitle( varName );
    MRField * fn = AllocNodeVar( nEqu );
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

void BcVisual::Calc(int bcType)
{
	UnsGrid * grid = Zone::GetUnsGrid();

	this->Calcf2n(bcType);

	e2n.resize(0);
	lcell.resize(0);
	rcell.resize(0);

	ResolveElementEdge();
}

void BcVisual::ResolveElementEdge()
{
	int nFace = this->f2n.size();
	int nSize = 2;

	set< Mid<int> > edgeSet;

	for (int fId = 0; fId < nFace; ++fId)
	{
		int nNode = this->f2n[fId].size();
		for (int iNode = 0; iNode < nNode; ++iNode)
		{
			int iNode0 = iNode;
			int iNode1 = (iNode + 1) % nNode;

			int ip1 = this->f2n[fId][iNode0];
			int ip2 = this->f2n[fId][iNode1];

			if (ip1 == ip2) continue;

			IntField eNodeId;
			eNodeId.push_back(ip1);
			eNodeId.push_back(ip2);

			IntField sortedNodeId = eNodeId;
			sort(sortedNodeId.begin(), sortedNodeId.end());

			int eIdddd = this->e2n.size();
			Mid<int> edge(nSize, eIdddd);
			edge.data = sortedNodeId;

			int  edgeIndex;
			set< Mid<int> >::iterator iter = edgeSet.find(edge);
			if (iter == edgeSet.end())
			{
				edgeIndex = -1;
			}
			else
			{
				edgeIndex = iter->id;
			}

			if (edgeIndex == -1)
			{
				edgeSet.insert(edge);
				this->lcell.push_back(fId);
				this->rcell.push_back(-1);
				this->e2n.push_back(eNodeId);
			}
			else
			{
				int ip = this->rcell[edgeIndex];
				if (ip != -1)
				{
					cout << "Fatal Error\n";
					cout << " edgeIndex = " << edgeIndex << "\n";
					for (int i = 0; i < this->e2n[edgeIndex].size(); ++i)
					{
						cout << this->e2n[edgeIndex][i] << " ";
					}
				}

				this->rcell[edgeIndex] = fId;
			}
		}
	}
}

void BcVisual::Calcf2n(int bcType)
{
	UnsGrid * grid = Zone::GetUnsGrid();
	FaceTopo * faceTopo = grid->faceTopo;
	LinkField & f2n = faceTopo->f2n;
	BcRecord * bcRecord = faceTopo->bcManager->bcRecord;

	IntField localf2n(4);
	set< HXSort< int > > sets;
	HXSort< int > data;

	this->f2n.resize(0);
	this->l2g.resize(0);

	int nBFace = grid->nBFace;

	for (int iFace = 0; iFace < nBFace; ++iFace)
	{
		if (bcType != bcRecord->bcType[iFace]) continue;

		int nNode = f2n[iFace].size();

		localf2n.resize(0);
		for (int iNode = 0; iNode < nNode; ++iNode)
		{
			int gId = f2n[iFace][iNode];

			data.value = gId;
			set< HXSort< int > >::iterator iter = sets.find(data);

			if (iter == sets.end())
			{
				this->l2g.push_back(gId);
				data.index = this->l2g.size() - 1;

				sets.insert(data);
				localf2n.push_back(data.index);
			}
			else
			{
				localf2n.push_back(iter->index);
			}

		}
		this->f2n.push_back(localf2n);
	}

	//if ( bcType == 3 && false )
	//{
	//    int le = 1360;
	//    int re = 1361;

	//    cout << "Elem id = " << le << " " << re << "\n";
	//    cout << " this->f2n.size() = " << this->f2n.size() << "\n";
	//    int nle = this->f2n[ le ].size();
	//    int nre = this->f2n[ re ].size();
	//    cout << "left elem node size =  " << nle << "\n";
	//    cout << "right elem node size =  " << nre << "\n";
	//    for ( int i = 0; i < nle; ++ i )
	//    {
	//        cout << this->f2n[ le ][ i ] << " ";
	//    }
	//    cout << "\n";
	//    for ( int i = 0; i < nre; ++ i )
	//    {
	//        cout << this->f2n[ re ][ i ] << " ";
	//    }
	//    cout << "\n";
	//}
}

void BcVisual::Dump(ostringstream & oss, VisualTool * visualTool, string & bcTitle)
{
	UnsGrid * grid = Zone::GetUnsGrid();

	for (UInt i = 0; i < visualTool->title.size(); ++i)
	{
		oss << visualTool->title[i] << endl;
	}

	int nNode = l2g.size();
	int nFace = e2n.size();
	int nElem = this->f2n.size();

	// output for Tecplot
	oss << "ZONE\n";

	oss << "T = " << bcTitle << endl;

	oss << "ZoneType = FEPolygon\n";

	oss << "Nodes    = " << nNode << endl;
	oss << "Faces    = " << nFace << endl;
	oss << "Elements = " << nElem << endl;
	oss << "NumConnectedBoundaryFaces = 0\n";
	oss << "TotalNumBoundaryConnections = 0\n";

	Plot::oss = &oss;
	Plot::DumpField(l2g, grid->nodeMesh->xN);
	Plot::DumpField(l2g, grid->nodeMesh->yN);
	Plot::DumpField(l2g, grid->nodeMesh->zN);

	int nVar = visualTool->qNodeField.size();
	for (int iVar = 0; iVar < nVar; ++iVar)
	{
		RealField & q = (*visualTool->qNodeField[iVar])[0];
		Plot::DumpField(l2g, q);
	}

	Plot::DumpFaceNodeLink(e2n);
	Plot::DumpFaceElementLink(lcell, nElem);
	Plot::DumpFaceElementLink(rcell, nElem);

}

void BcVisual::DumpDebug(ostringstream & oss, VisualTool * visualTool, string & bcTitle)
{
	UnsGrid * grid = Zone::GetUnsGrid();

	int nNode = l2g.size();
	int nFace = e2n.size();
	int nElem = this->f2n.size();

	// output for Tecplot

	oss << " VARIALBES = ";
	oss << " \"x\" ";
	oss << " \"y\" ";
	oss << " \"z\" ";
	oss << "\n";
	//oss << " ZONE N = " << nNode << " E = " << nElem << " F = FEPOINT, ET = TRIANGLE \n";
	//oss << " ZONE N = " << nNode << " E = " << nElem << " F = FEPOINT, ET = QUADRILATERAL \n";
	oss << "title = \"THE FLOW FIELD OF ONEFLOW\" \n";
	oss << "VARIALBES = \"x\" \"y\" \"z\" " << "\n";

	oss << " ZONE N = " << nNode << " E = " << nElem << " F = FEPOINT, ET = QUADRILATERAL \n";

	for (int iNode = 0; iNode < nNode; ++iNode)
	{
		int id = l2g[iNode];
		oss << grid->nodeMesh->xN[id] << " ";
		oss << grid->nodeMesh->yN[id] << " ";
		oss << grid->nodeMesh->zN[id] << " ";
		oss << "\n";
	}

	for (int iElem = 0; iElem < nElem; ++iElem)
	{
		int p1 = this->f2n[iElem][0] + 1;
		int p2 = this->f2n[iElem][1] + 1;
		int p3 = this->f2n[iElem][2] + 1;
		int p4 = this->f2n[iElem][3] + 1;

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

	int nNode = l2g.size();
	int nFace = e2n.size();
	int nElem = this->f2n.size();

	fstream file;
	string fileName = "test.plt";
	ONEFLOW::OpenPrjFile(file, fileName, ios_base::out);

	file << " VARIALBES = ";
	file << " \"x\" ";
	file << " \"y\" ";
	file << " \"z\" ";
	file << "\n";

	IntField eList;
	eList.push_back(1360);
	eList.push_back(1361);
	IntField nList, nList1, nList2;

	IntField localf2n(4);
	set< HXSort< int > > sList;
	HXSort< int > data;
	int iCount = 0;
	for (int e = 0; e < eList.size(); ++e)
	{
		int ee = eList[e];
		int nsize = this->f2n[ee].size();
		for (int in = 0; in < nsize; ++in)
		{
			int ip = this->f2n[ee][in];
			nList1.push_back(ip);
			data.value = ip;
			set< HXSort< int > >::iterator iter = sList.find(data);

			if (iter == sList.end())
			{
				data.index = iCount;
				sList.insert(data);
				nList2.push_back(ip);

				nList.push_back(iCount);
				++iCount;
			}
			else
			{
				nList.push_back(iter->index);
			}
		}
	}

	file << " ZONE N = " << iCount << " E = " << eList.size() << " F = FEPOINT, ET = QUADRILATERAL \n";
	int width = 20;
	int pre = 20;
	for (int iNode = 0; iNode < iCount; ++iNode)
	{
		int id = l2g[nList2[iNode]];
		file << setw(width) << setprecision(pre) << grid->nodeMesh->xN[id] << " ";
		file << setw(width) << setprecision(pre) << grid->nodeMesh->yN[id] << " ";
		file << setw(width) << setprecision(pre) << grid->nodeMesh->zN[id] << " ";
		file << "\n";
	}

	for (int iNode = 0; iNode < nList1.size(); ++iNode)
	{
		int id = l2g[nList1[iNode]];
		file << id << "\n";
	}


	int pos = 0;
	for (int e = 0; e < eList.size(); ++e)
	{
		int p1 = nList[pos + 0] + 1;
		int p2 = nList[pos + 1] + 1;
		int p3 = nList[pos + 2] + 1;
		int p4 = nList[pos + 3] + 1;
		pos += 4;

		file << p1 << " ";
		file << p2 << " ";
		file << p3 << " ";
		file << p4 << " ";
		file << "\n";
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

	this->CalcNodeField(&visualTool);

	ostringstream oss;

	this->ShowBc(oss, &visualTool);
	//this->ShowBcDebug( oss, & visualTool );
	//this->ShowBcDebugTest( oss, & visualTool );

	if (this->NeedVisualField())
	{
		this->ShowField(oss, &visualTool);
	}


	ToDataBook(ActionState::dataBook, oss);
}

void UVisualize::ExtractLinkNum(LinkField & f2n, IntField & fnNumber)
{
	int nSize = f2n.size();
	fnNumber.resize(nSize);
	for (int i = 0; i < nSize; ++i)
	{
		fnNumber[i] = f2n[i].size();
	}
}

int UVisualize::GetTotalNumFaceNodes(LinkField & f2n)
{
	int totalNumFaceNodes = 0;
	int nSize = f2n.size();
	for (int i = 0; i < nSize; ++i)
	{
		totalNumFaceNodes += f2n[i].size();
	}
	return totalNumFaceNodes;
}

void UVisualize::ShowField(ostringstream & oss, VisualTool * visualTool)
{
	UnsGrid * grid = Zone::GetUnsGrid();

	FaceTopo * faceTopo = grid->faceTopo;
	LinkField & f2n = faceTopo->f2n;

	int nNode = grid->nNode;
	int nCell = grid->nCell;
	int nFace = grid->nFace;

	for (UInt i = 0; i < visualTool->title.size(); ++i)
	{
		oss << visualTool->title[i] << endl;
	}

	int totalNumFaceNodes = this->GetTotalNumFaceNodes(f2n);

	// output for Tecplot
	oss << "ZONE\n";
	if (Dim::dimension == THREE_D)
	{
		oss << "ZoneType = FEPolyhedron\n";
	}
	else
	{
		oss << "ZoneType = FEPolygon\n";
	}
	oss << "Nodes    = " << nNode << endl;
	oss << "Faces    = " << nFace << endl;
	oss << "Elements = " << nCell << endl;
	oss << "TotalNumFaceNodes = " << totalNumFaceNodes << endl;
	oss << "NumConnectedBoundaryFaces = 0\n";
	oss << "TotalNumBoundaryConnections = 0\n";
	Plot::oss = &oss;
	Plot::DumpField(grid->nodeMesh->xN);
	Plot::DumpField(grid->nodeMesh->yN);
	Plot::DumpField(grid->nodeMesh->zN);

	int nVar = visualTool->qNodeField.size();
	for (int iVar = 0; iVar < nVar; ++iVar)
	{
		RealField & q = (*visualTool->qNodeField[iVar])[0];
		Plot::DumpField(q);
	}

	if (Dim::dimension == THREE_D)
	{
		Plot::DumpFaceNodeNumber(f2n);
	}

	Plot::DumpFaceNodeLink(f2n);
	Plot::DumpFaceElementLink(faceTopo->lCell, nCell);
	Plot::DumpFaceElementLink(faceTopo->rCell, nCell);
}

void UVisualize::ShowBc(ostringstream & oss, VisualTool * visualTool)
{
	if (IsTwoD()) return;
	UnsGrid * grid = Zone::GetUnsGrid();

	IntField bcTypeList;
	grid->faceTopo->bcManager->CalcBcType(bcTypeList);
	int nBcType = bcTypeList.size();

	for (int iBcType = 0; iBcType < nBcType; ++iBcType)
	{
		int bcType = bcTypeList[iBcType];

		if (BC::IsInterfaceBc(bcType)) continue;

		string bcTitle = AddString("\"", ZoneState::zid, "BC=", bcType, "\"");

		BcVisual bcVisual;

		bcVisual.Calc(bcType);

		bcVisual.Dump(oss, visualTool, bcTitle);
	}
}

void UVisualize::ShowBcDebugTest(ostringstream & oss, VisualTool * visualTool)
{
	if (IsTwoD()) return;
	UnsGrid * grid = Zone::GetUnsGrid();

	IntField bcTypeList;
	grid->faceTopo->bcManager->CalcBcType(bcTypeList);
	int nBcType = bcTypeList.size();

	for (int iBcType = 0; iBcType < nBcType; ++iBcType)
	{
		int bcType = bcTypeList[iBcType];

		if (BC::IsInterfaceBc(bcType)) continue;
		if (bcType != BC::SYMMETRY) continue;

		string bcTitle = AddString("\"", ZoneState::zid, "BC=", bcType, "\"");

		BcVisual bcVisual;

		bcVisual.Calc(bcType);

		bcVisual.DumpDebug(oss, visualTool, bcTitle);
	}
}

void UVisualize::CalcNodeField(VisualTool * visualTool)
{
	UnsGrid * grid = Zone::GetUnsGrid();
	MRField * q = ONEFLOW::GetFieldPointer< MRField >(grid, "q");

	MRField * rn = visualTool->AddField((*q)[IDX::IR], "r");
	MRField * un = visualTool->AddField((*q)[IDX::IU], "u");
	MRField * vn = visualTool->AddField((*q)[IDX::IV], "v");
	MRField * wn = visualTool->AddField((*q)[IDX::IW], "w");
	MRField * pn = visualTool->AddField((*q)[IDX::IP], "p");

	MRField * gaman = CreateNodeVar("gama");
	MRField * machn = visualTool->CreateField("mach");
	CalcMach(rn, un, vn, wn, pn, gaman, machn);
	delete gaman;

	MRField * tempr = ONEFLOW::GetFieldPointer< MRField >(grid, "tempr");
	visualTool->AddField((*tempr)[IDX::ITT], "tempr");

	if (vis_model.vismodel > 0)
	{
		visualTool->AddField("visl");
		visualTool->AddField("vist");
	}
}

void CalcMach(MRField * r, MRField * u, MRField * v, MRField * w, MRField * p, MRField * gama, MRField * mach)
{
    UnsGrid * grid = Zone::GetUnsGrid();
    int nNode = grid->nNode;
	for (int iNode = 0; iNode < nNode; ++iNode)
	{
		Real rm = (*r)[0][iNode];
		Real um = (*u)[0][iNode];
		Real vm = (*v)[0][iNode];
		Real wm = (*w)[0][iNode];
		Real pm = (*p)[0][iNode];

		Real gm = (*gama)[0][iNode];
		int startStrategy = ONEFLOW::GetDataValue< int >("startStrategy");
		if (startStrategy == 2)
		{
			(*mach)[0][iNode] = 0;
		}
		else
		{
			Real v2 = SQR(um, vm, wm);
			Real c2 = gm * pm / rm;
			Real mm = sqrt(v2 / c2);
			(*mach)[0][iNode] = mm;
		}
	}
		
}

EndNameSpace