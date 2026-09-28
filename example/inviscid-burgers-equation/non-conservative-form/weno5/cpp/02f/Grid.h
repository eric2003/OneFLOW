#pragma once
#include <string>
#include <vector>

class Region
{
public:
    Region();
    ~Region();
public:
    std::vector<int> start;
    std::vector<int> end;
public:
    Region & operator = ( const Region & rhs );
    void SetRegion( std::vector<int> & pnts );
    void Print();
};

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

class ZoneBc
{
public:
    ZoneBc();
    ~ZoneBc();
public:
    int bcType;
    int zoneid;
    std::vector<int> pnts;
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
    int zoneIndex;
    std::vector<int> nijk;
    std::vector<ZoneBc *> bccos;
    std::vector<ZoneBc1To1 *> bc1to1s;
    std::vector<Coor *> coors;
};

class Trans
{
public:
    static int M[ 3 ][ 3 ];
    static std::vector<int> transform;
    static int sgn( int x );
    static int del( int x, int y );
    static void ZeroMatrix();
    static void CalcTransformMatrix();
};

class Transform
{
public:
    Transform();
    ~Transform();
private:
    std::vector<int> diff;
    std::vector<int> mul;
public:
    int Mt[ 3 ][ 3 ];
    std::vector<int> begin1;
    std::vector<int> begin2;
    std::vector<int> transform;
public:
    void Init();
    void MapIndex( std::vector<int> & index1, std::vector<int> & index2 );
    void Multiply( std::vector<int> & a, std::vector<int> & b );
};
