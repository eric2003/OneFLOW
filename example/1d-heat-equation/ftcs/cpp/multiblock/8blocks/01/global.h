#pragma once
#include <vector>
#include <map>
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

class Transform;

class Interface
{
public:
    std::vector<int> zoneList;
    std::vector<int> faceidList;
    std::vector<int> ijk_ghosts;
    std::vector<int> ijk_donors;
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
    static int nt;
    static int cell_dim;
    static int phys_dim;
};

