#pragma once
#include <string>
#include <vector>
#include "CgnsHeader.h"

std::string GetCgnsFileTypeName( int file_type );

class Zone;
class Base;
class BBase;

Zone * GetZone( int fileId, int BaseId, int zoneId );
void ReadCgnsFile( const std::string &fileName );
BBase * GetBBase();

class BBase
{
public:
    BBase();
    ~BBase();
public:
    int nBases;
    std::vector<Base *> bases;
public:
    int fileId;
    int file_type;
    float fileVersion;
    int precision;
public:
    void AllocateBases( int nBases );
    void ReadFile( const std::string & fileName );
    Base * GetBase() { return bases[0]; }
    Base * GetBase( int baseId ) { return bases[baseId-1]; }
    Zone * GetZone();
    Zone * GetZone( int zoneId );
};

class Base
{
public:
    Base();
    ~Base();
public:
    int nZones;
    std::vector<Zone *> zones;
public:
    double double_base_id;
    int celldim;
    int phydim;
    char33 name;
    int baseId;
public:
    void AllocateZones( int nZones );
    void ReadBase( int fileId, int baseId );
    void ReadBaseUnits( int fileId, int baseId );
};

class Coor;
class Field;
class Solution;
class Section;
class Bc;
class Zone
{
public:
    Zone();
    ~Zone();
public:
    ZoneType_t zoneType;
    char33 zoneName;
    int zoneId;
public:
    cgsize_t isize[ 9 ];
    cgsize_t irmin[ 3 ], irmax[ 3 ], cellSize[ 3 ];
    cgsize_t nNodes, nCells;
    cgsize_t nFaces;
    cgsize_t totalNumFaceNodes;
    int nCoords;
    std::vector<Coor *> coors;
    int nSolutions = 0;
    std::vector<Solution *> solutions;

    int nFields;
    std::vector<Field *> fields;
    int nSections;
    std::vector<Section *> sections;
    std::vector<cgsize_t> left_elements;
    std::vector<cgsize_t> right_elements;
    Bc * bc;
public:
    void ReadZone( int fileId, int baseId, int zoneId );
    void ReadCoordinates( int fileId, int baseId, int zoneId );
    void ReadElements( int fileId, int baseId, int zoneId );
    void ReadBoundaries( int fileId, int baseId, int zoneId );
    void ReadFlowSolution( int fileId, int baseId, int zoneId );
public:
    void AllocateSolutions( int nSolutions );
public:
    void SetDimensions();
    void AllocateCoors( int nCoords );
    void DeAllocateCoors();
    void AllocateSections( int nSections );
    void DeAllocateSections();
    void DumpZone();
    void DrawZone();
    void DumpCoor( std::ostringstream & oss );
    void DumpFaceNodeNumber( std::ostringstream & oss );
    void DumpFaceNodeLink( std::ostringstream & oss );
    void DumpFaceElementLink( std::ostringstream & oss, std::vector<cgsize_t> & elementId );
    void CalcFaceElementLink();
    void DumpSectionMesh( std::ostringstream & oss );

};

class Coor
{
public:
    Coor();
    ~Coor();
public:
    DataType_t dataType;
    char33 coorName;
    void * data;
    cgsize_t nNodes;
public:
    void ReadCoor( int fileId, int baseId, int zoneId, int coordId, Zone * zone );
    void AllocateData( int nNodes );
    void DeAllocateData();
    void DumpCoor();
};

class Field
{
public:
    Field();
    ~Field();
public:
    DataType_t dataType;
    char33 fieldName;
    int fieldId = -1;
    double doubleFieldId = -1;
    void * data;
    cgsize_t nNodes;
    DataClass_t dataclass;
    int dimflag=-1; //0: nondimensional 1: dimensional
    int nExponents;
    std::vector<float> exponents;
public:
    void ReadField( int fileId, int baseId, int zoneId, int solutionId, int fieldId, Zone * zone );
    void AllocateData( int nNodes );
    void DeAllocateData();
};

class Solution
{
public:
    Solution();
    ~Solution();
public:
    void AllocateFields( int nFields );
    void ReadField( int fileId, int baseId, int zoneId, int solutionId, Zone * zone );
public:
    char33 solutionName;
    GridLocation_t location;
    double solver_double_id = -1;
    int nFields;
    std::vector<Field *> fields;
};

class Uns2D;

class Section
{
public:
    Section();
    ~Section();
public:
    ElementType_t elementType;
    char33 sectionName;
    cgsize_t startId;
    cgsize_t endId;
    cgsize_t nElements;
    int nbndry;
    int iparentflag;
    cgsize_t elementDataSize;
    int pos_shift;
    std::vector<cgsize_t> conn;
    std::vector<cgsize_t> conn_offsets;
    Zone * zone;
    std::vector<cgsize_t> left_elements;
    std::vector<cgsize_t> right_elements;
public:
    void Read( int fileId, int baseId, int zoneId, int sectionId );
    bool IsMixedSection();
    void PrintElementInfo();
    void PrintPolygonFaceInfo();
    void PrintPolyhedronElementInfo();
    int  GetNumFaceNodes();
    void DumpSectionMesh( std::ostringstream & oss );
};
