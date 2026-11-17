#pragma once
#include <iostream>
#include <vector>

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

class ZoneBc
{
public:
    ZoneBc();
    ~ZoneBc();
public:
    int bcType;
    std::vector<int> pnts;
};

class Global
{
public:
    Global();
    ~Global();
public:
    static std::vector<Zone *> zones;
    static void FreeData();
    static int cell_dim;
    static int phys_dim;
};

void ReadCgnsGrid( const std::string & filename );
