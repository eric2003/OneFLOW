#include <string>
#include <vector>
#include "cgnslib.h"

class Boco
{
public:
    Boco();
    ~Boco();
public:
    int bocoId;
    double bc_double_id;
    std::string name;
    BCType_t bcType;
    PointSetType_t pointSetType;
    GridLocation_t gridLocation;
    GridLocation_t modifiedLocation;
    GridConnectivityType_t gridConnType;  //Overset, Abutting, Abutting1to1
    cgsize_t nElements;
    DataType_t normalDataType;
    cgsize_t normalListSize;
    int normalIndex[ 3 ];
    int nDataSets;
    std::vector<cgsize_t> conn;
public:
    void Read( int fileId, int baseId, int zoneId, int bocoId );
    void ReadGridLocation( int fileId, int baseId, int zoneId );
};

class Bc
{
public:
    Bc();
    ~Bc();
public:
    void ReadBoundaries( int fileId, int baseId, int zoneId );
    void AllocateBocos();
public:
    int nBocos;
    std::vector< Boco * > bocos;
};

