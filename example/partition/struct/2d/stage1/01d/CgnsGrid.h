#pragma once
#include <iostream>
#include <vector>
#include <map>

class Grid;

class Zone;

class Grid
{
public:
    Grid();
    ~Grid();
public:
    int nZone;
    std::vector<double> x;
    std::vector<Zone *> zones;
};

class ZoneBc;

class Coor
{
public:
    Coor();
    ~Coor();
public:
    std::string coorname;
    int nNodes;
    std::vector<int> nijk;
    std::vector<char> coord;
public:
    void DumpCoor();
};

class Zone
{
public:
    Zone();
    ~Zone();
public:
    std::vector<int> nijk;
    std::vector<ZoneBc *> bccos;
    std::vector<Coor *> coors;
};

class ZoneBcPatch;
class ZoneBc
{
public:
    ZoneBc();
    ~ZoneBc();
public:
    int bcType;
    int zoneid;
    std::vector<int> pnts;
    ZoneBcPatch * patch = nullptr;
};

class ZoneBcPatch
{
public:
    ZoneBcPatch();
    ~ZoneBcPatch();
public:
    int donor_zoneid;
    std::vector<int> donor_pnts;
};

class BaseZone
{
public:
    std::string zone_name;
    bool operator < ( const BaseZone & rhs ) const
    {
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

class Global
{
public:
    Global();
    ~Global();
public:
    static int cell_dim;
    static int phys_dim;
    static std::vector<Zone *> zones;
    static BaseZoneList zone_names;
public:
    static void FreeData();
};

void ReadCgnsGrid( const std::string & filename );
void ReadCgnsGridBaseZone( const std::string & filename );