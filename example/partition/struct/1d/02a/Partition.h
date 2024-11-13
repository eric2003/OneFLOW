#pragma once
#include <iostream>
#include <vector>
#include <set>

class SplitZone;

struct CmpSplitZone 
{
    bool operator() ( SplitZone * a, SplitZone * b ) const;
};

class Partition
{
public:
    Partition();
    ~Partition();
protected:
    std::set < SplitZone * > zoneStorage;
    std::set < SplitZone *, CmpSplitZone > unsignedZoneGroup;
    std::vector< SplitZone * > refZone;
    std::vector< std::vector< SplitZone * > > procSplitZone;
    std::vector<double> procSpeed;
    std::vector<double> timeProc;
    double aveFloatOp, delta;
    double idealUniProcTime, eps;
    int minNumberCell;
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

