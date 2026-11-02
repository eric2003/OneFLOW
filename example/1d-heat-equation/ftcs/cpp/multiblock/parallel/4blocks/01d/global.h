#pragma once
#include <vector>
#include <map>
#include <unordered_map>
#include <set>
#include <algorithm>
#include <string>
#include "CgnsGrid.h"

class Grid;

class BaseZone
{
public:
    //int baseId;
    std::string zone_name;
    bool operator < ( const BaseZone & rhs ) const
    {
        //if ( this->baseId != rhs.baseId )
        //{
        //    return this->baseId < rhs.baseId;
        //}

        return this->zone_name < rhs.zone_name;
    } 
};

class BaseZoneList
{
public:
    std::map<BaseZone, int> basezone_map;
public:
    void AddBaseZone( const BaseZone & baseZone )
    {
        std::map<BaseZone, int>::iterator iter;
        iter = basezone_map.find( baseZone );
        if ( iter == basezone_map.end() )
        {
            int id = basezone_map.size();
            basezone_map.insert( std::make_pair( baseZone, id ) );
        }
    }

    int FindBaseZone( const BaseZone & baseZone )
    {
        std::map<BaseZone, int>::iterator iter = basezone_map.find( baseZone );
        if ( iter != basezone_map.end() )
        {
            return iter->second;
        }
        return -1;
    }
};

class Grid
{
public:
    int zoneIndex;
    std::vector<double> x;
};

class Face
{
public:
    int zone = 0;
    int i = 0;
    int j = 1;
    int k = 1;
public:
    bool operator < ( const Face & rhs ) const;
    bool operator == ( const Face & rhs ) const;
};

class FacePair
{
public:
    Face left;
    Face right;
public:
    bool operator < ( const FacePair & rhs ) const;
    void AddPair( const Face &face1, const Face &face2);
};

class Transform;


class Interface
{
public:
    int zoneid;
    std::vector<int> zoneList;
    std::vector<int> global_faceids;
    std::vector<int> local_faceids;
    std::vector<int> ijk_ghosts;
    std::vector<int> ijk_donors;
    std::vector<int> ijk_local_donors;
    std::vector<double> data_recv;
    std::vector<double> data_send;
    std::unordered_map<int, int> global_local_face_map;
public:
    std::vector<int> neighbor_donor_zones;
    std::vector<std::vector<int>> neighbor_donorfaces;
    std::vector<std::vector<int>> sub_local_faceids;
    std::vector<int> send_to_zones;
    std::vector<std::vector<int>> donorfaces_for_send;
    std::vector<std::vector<int>> donorijk_for_send;
    std::vector<std::vector<double>> donordata_for_send;
    std::vector<int> tmp_recv_donor_zones;
    std::vector<std::vector<int>> tmp_recv_donorfaces;
public:
    std::vector<int> recv_donor_zones;
    std::vector<std::vector<double>> recv_donor_fields;
public:
    void CalcInterface( Transform * transform, std::vector<int> & start, std::vector<int> & end, int donor_zoneid );
};

class Field;
class InterFaceZone;

class Global
{
public:
    static std::vector<Grid *> grids;
    static std::vector<Field *> fields;
public:
    static std::vector<Zone *> zones;
    static BaseZoneList zone_names;
    static std::vector<Interface *> interfaces;
public:
    static std::map<Face, int> faceMap;
    static std::map<FacePair, int> facePairMap;
    static std::vector<FacePair> facePairList;
    static std::vector<std::set<int>> donor_zone_sets;
    static std::vector<std::vector<int>> donor_zones;
public:
    static int nt;
    static int cell_dim;
    static int phys_dim;
public:
    static void InsertFaceMap( const Face & face );
    static int InsertFacePairMap( const FacePair & facePair );
    static void SendGeom( int zone, int donorzone, std::vector<int> & donorfaces );
    static void ReSendGeomTest( int zone, int zone_to_send, std::vector<int> & donorfaces );
};

