#pragma once
#include <iostream>
#include <vector>
#include <set>

class SplitZone;
class RangeRegion;

void GetDomainRegion( RangeRegion & zoneBox, int iDomain, RangeRegion & domainBox );

class DimInfo
{
public:
    DimInfo();
    ~DimInfo();
public:
    int dim;
    std::vector<int> dimList;
    std::vector<int> oriPoint;
    std::vector<int> isize;
    std::vector<int> rmin;
    std::vector<int> rmax;
public:
    int GetSize() const { return dim; }
    DimInfo & operator = ( const DimInfo & rhs );
public:
    int ComputeNumNodes();
    void ComputeISize();
    void ComputeRminmax( std::vector<int> & ref );
    void Init( std::vector<int> & dimList );
    void SetNewDimension( int axis, int newStart, int newWidth );
    void GetZoneGlobalBox( RangeRegion & zoneBox );
    void GetZoneLocalBox( RangeRegion & globalBox, RangeRegion & localBox );
};

class SplitZone;
class SplitBcPatch;

class RangeRegion
{
public:
    RangeRegion( int nSize = 1 );
    ~RangeRegion();
public:
    int regionAxis = -1;
    std::vector<int> start;
    std::vector<int> end;
public:
    RangeRegion & operator = ( const RangeRegion & rhs );
    bool operator == ( const RangeRegion & rhs ) const;
    void SetRegion( std::vector<int> &pnts );
    int ComputeRegionAxis();
};

class SplitBc
{
public:
    SplitBc();
    ~SplitBc();
public:
    int bcType;
    RangeRegion region;
public:
    SplitZone * splitZone = nullptr;
    SplitBcPatch * patch = nullptr;
    std::vector< SplitBc * > * child = nullptr;
public:
    void SetRegion( std::vector<int> &pnts );
    void SetRegion( RangeRegion &region );
    void ChangeRegionToLocalCoordinate();
public:
    void SetZoneBcFromParentBc( SplitBc * parentSplitBc );
    void CreatePatchBc( SplitBc * parentSplitBc );
};

class SplitBcPatch
{
public:
    SplitBcPatch ();
    ~SplitBcPatch();
public:
    int bcType;
    RangeRegion region;
public:
    SplitZone * splitZone;
public:
    void SetRegion( std::vector<int> &pnts );
    void SetRegion( RangeRegion &region );
};


class SplitZone
{
public:
    SplitZone();
    ~SplitZone();
public:
    DimInfo dimInfo;
    int zoneIndex, iProc;
    //int oriZoneId;
protected:
    SplitZone * parent;
public:
    std::vector< SplitBc * > splitBcList;
public:
    int GetNCell();
    void ComputeRminmax();
    void GetRootInfo( SplitZone *& root, std::vector<int> & ref );
    void SetZoneIndex( int zoneIndex ) { this->zoneIndex = zoneIndex; }
    int  GetZoneIndex() { return zoneIndex; }
    void SetProcId( int iProc ) { this->iProc = iProc; }
    void SetParent( SplitZone * parent ) { this->parent = parent; };
    std::vector< SplitBc * > & GetSplitBcList() { return this->splitBcList; };
    void CreateBcFromParent();
    void CreateSplitBc( SplitBc * parentSplitBc );
    void AddSplitBc( SplitBc * splitBc );
    void GetZoneGlobalBox( RangeRegion & zoneBox );
    void GetZoneLocalBox( RangeRegion & globalBox, RangeRegion & localBox );
    SplitZone * GetRootZone();
public:
    void Split( SplitZone *& zoneL, SplitZone *& zoneR, double nCell );
    void SetLeftDimension( std::vector<int> & dimList, int axisNeedSplit, int iSplit );
    void SetRightDimension( std::vector<int> & dimList, int axisNeedSplit, int iSplit );
private:
    int FindAxisNeedSplit();
    int FindSplitPosition( int axisNeedSplit, double nCell );
    bool IsPoleBc( int splitDirection );
    int GetCrossSection( int axis );
};

class Partition
{
public:
    Partition();
    ~Partition();
protected:
    std::set < SplitZone * > zoneStorage;
    std::set < SplitZone * > unsignedZoneGroup;
    std::vector< SplitZone * > refZone;
    std::vector< std::vector< SplitZone * > > procSplitZone;
    std::vector<double> procSpeed;
    std::vector<double> timeProc;
    double aveFloatOp, delta;
    double idealUniProcTime, eps;
    int minNode;
    int nTZone;
    int nOriZone;
    int nProc;
public:
    std::string inputName;
    std::string outName;
public:
    void Solve();
    void ReadGrid();
    void Split();
    void InitSplit();
    void BinarySplit();
    //void GreedySplit();
public:
    void AddNewZone( SplitZone * zone );
    void AddZoneToUnsignedGroup( SplitZone * zone );
    void RemoveZoneFromUnsignedGroup( SplitZone * zone );
    void AddZoneToProcess( int iProc, SplitZone * zone );
    void DumpPartitionedGridFile();
    void ModifyZoneIndex();
    void FreeAllZones();
    int GetNZones( int iProc );
    SplitZone * GetSplitZone( int iProc, int zoneId );
    void ComputeCoor( SplitZone * splitZone, int icoord, std::vector<char> & coord );
    void DumpBc( SplitZone * splitZone );
private:
    int GetMinValueIndex( double * data, int numberOfData );
    int GetMinLoadProc();
    SplitZone * GetLargestZone();
    SplitZone * Split( SplitZone * zone, double nCell );
};

