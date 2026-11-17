#pragma once
#include <iostream>
#include <vector>
#include <set>

class SplitZone;

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
};

class SplitZone
{
public:
    SplitZone();
    ~SplitZone();
public:
    DimInfo dimInfo;
    int zoneIndex, iProc;
    int oriZoneId;
protected:
    SplitZone * parent;
public:
    int GetNCell();
    void ComputeRminmax();
    void GetRootInfo( SplitZone *& root, std::vector<int> & ref );
    void SetZoneIndex( int zoneIndex ) { this->zoneIndex = zoneIndex; }
    int  GetZoneIndex() { return zoneIndex; }
    void SetProcId( int iProc ) { this->iProc = iProc; }
    void SetParent( SplitZone * parent ) { this->parent = parent; };
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
    double minNode;
    int nTZone;
    int nOriZone;
    int nProc;
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
    void FreeAllZones();
    int GetNZones( int iProc );
    SplitZone * GetSplitZone( int iProc, int zoneId );
private:
    int GetMinValueIndex( double * data, int numberOfData );
    int GetMinLoadProc();
    SplitZone * GetLargestZone();
    SplitZone * Split( SplitZone * zone, double nCell );
};

