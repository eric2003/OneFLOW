#pragma once
#include <string>
#include <vector>

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
    void DumpCoorX( std::vector<double> & x );
};

class ZoneBcPatch
{
public:
    ZoneBcPatch();
    ~ZoneBcPatch();
public:
    int donor_zoneid;
    std::vector<int> transform;
    std::vector<int> donor_pnts;
};


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

class ZoneBc1To1
{
public:
    ZoneBc1To1();
    ~ZoneBc1To1();
public:
    int zoneid;
    int donor_zoneid;
    std::vector<int> pnts;
    std::vector<int> donor_pnts;
    std::vector<int> transform;
};


class Zone
{
public:
    Zone();
    ~Zone();
public:
    std::vector<int> nijk;
    std::vector<ZoneBc *> bccos;
    std::vector<ZoneBc1To1 *> bc1to1s;
    std::vector<Coor *> coors;
};

void ReadCgnsGridBaseZone( const std::string & filename );
void ReadCgnsGrid( const std::string & filename );