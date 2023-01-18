/*---------------------------------------------------------------------------*\
OneFLOW - LargeScale Multiphysics Scientific Simulation Environment
Copyright (C) 2017-2023 He Xin and the OneFLOW contributors.
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
#include <vector>
#include <set>
#include <map>
#include <string>

class Grid
{
public:
    Grid();
    ~Grid();
public:
    std::vector<float> xcoor;
};

class BoundarySolver;
class InterfaceSolver;

const int BCInterface = -1;
//const int BCInflow = 1;
//const int BCOutflow = 2;

class BcTypeMap
{
public:
    BcTypeMap();
    ~BcTypeMap();
public:
    void Init();
    std::string GetBcName( int bcId );
private:
    std::map< int, std::string > bc_map;
};

class Interface_t
{
public:
    Interface_t();
    ~Interface_t();
public:
    int cell;
    int zone;
    int ghost_cell;
    int face;
};

class IFaceBasic
{
public:
    IFaceBasic();
    ~IFaceBasic();
public:
    int GetNeighborZone() const;
    int GetNeighborFace() const;
    Interface_t * GetCurrentInterface_t();
    Interface_t * GetNeighborInterface_t();
public:
    int global_face_id;
    Interface_t face_pair[ 2 ];
    int count;
    int bctype;
public:
    int nid;
};

class PhysicalBoundarySolver
{
public:
    PhysicalBoundarySolver();
    ~PhysicalBoundarySolver();
public:
    void Init();
public:
    std::vector<int> physical_faceids;
    std::vector<int> physical_ghost_cells;
    std::vector<int> physical_bctypes;
    std::vector<IFaceBasic * > physical_bc_list;
};

class BoundarySolver
{
public:
    BoundarySolver();
    ~BoundarySolver();
public:
    int GetNBFace() { return bctypes.size(); };
    int GetNIFace();
public:
    void Init( int zoneId, int nZones, int ni );
    void FillBCPoints();
    void MarkInterface();
    void PrintBcInfo();
    void SetBcType( int bc_face_id, int bcType );
    IFaceBasic * FindIFaceBasic( int id );
    void InsertVector( std::vector<int> & a, std::vector<int> & b );
public:
    int zoneId, nZones;
    int ni;
    std::map<int,IFaceBasic *> global_face_map;
    std::map<int,int> local_face_map;
    InterfaceSolver* interfaceSolver;
    PhysicalBoundarySolver * physicalSolver;
public:
    std::vector<int> bctypes;
    std::vector<int> bc_faceids;
    std::vector<int> bc_ghostcells;
};

class IData
{
public:
    IData();
    ~IData();
public:
    void Print();
    void SetSendQs( float * q );
    void GetRecvQs( float * q );
public:
    int zone_id;
    std::vector< int > send_cells; //from myzone to neighbor zone(zone_id)
    std::vector< float > send_qs; //from myzone to neighbor zone(zone_id)
    std::vector< int > recv_cells; //from neighbor zone(zone_id) to myzone
    std::vector< float > recv_qs; //from neighbor zone(zone_id) to myzone

    std::vector< int > interface_faceids;
    std::vector< int > interface_neighbor_faceids;
public:
    std::vector<IFaceBasic * > tmplist;
};

class InterfaceSolver
{
public:
    InterfaceSolver();
    ~InterfaceSolver();
public:
    void Init();
    void SortInterface();
    void SetInterface();
    void SetNeighborData();
    void ReorderCurrentData();
    void SwapData();
    void SwapData( float * q );
    void PrepareSendData( float * q );
    void ShowInfo( int iface );
    float GetBcValue( int iface );
    int GetNIFace()
    { 
        return this->interface_bctypes.size();
    }
public:
    int zone_id;
public:
    std::vector<int> interface_neighbor_zoneids;
    std::vector<int> interface_neighbor_faceids;
    std::vector<int> interface_neighbor_cells;
    std::vector<int> interface_neighbor_ghost_cells;
    std::vector<int> interface_ipos;
    std::vector<int> interface_jpos;
    std::vector<int> interface_faceids;
    std::vector<int> interface_cells;
    std::vector<int> interface_ghost_cells;
    std::vector<int> interface_bctypes;
public:
    std::set<int> neighbor_zones;
    std::vector<IFaceBasic * > interlist;
    std::vector<IData> neighbor_datas;
};
