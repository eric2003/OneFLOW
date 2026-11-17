#pragma once
#include <iostream>
#include <vector>
#include <set>

class Dim
{
public:
    static int dim;
};

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
    void Init( int dim );
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
    std::vector<int> start;
    std::vector<int> end;
public:
    void Init( int nSize );
    RangeRegion & operator = ( const RangeRegion & rhs );
    bool operator == ( const RangeRegion & rhs ) const;
    bool InRegion( const RangeRegion & region );
    void SetRegion( std::vector<int> &pnts );
    void Normalize();
    void Localize( std::vector<int> & oriPoint );
    void AddShift( std::vector<int> & oriPoint );
    void PrintGlobal( std::vector<int> & oriPoint );
    void Print();
};

class PhysicalSplitBc
{
public:
    PhysicalSplitBc();
    ~PhysicalSplitBc();
public:
    int bcType;
    RangeRegion region;
public:
    SplitZone * splitZone = nullptr;
public:
    void SetRegion( std::vector<int> &pnts );
    void SetRegion( RangeRegion &region );
    void ChangeRegionToLocalCoordinate();
};

class BasicSplitBc
{
public:
    BasicSplitBc();
    ~BasicSplitBc();
public:
    RangeRegion region;
    SplitZone * splitZone = nullptr;
public:
    void SetRegion( std::vector<int> & pnts );
    void SetRegion( RangeRegion & region );
};

class trans
{
public:
    static int M[ 3 ][ 3 ];
    static std::vector<int> transform;
    static int sgn( int x );
    static int del( int x, int y );
    static void ZeroMatrix();
    static void CalcTransformMatrix();
};

class InterfaceSplitBc;
class SplitZone;

class InterfaceInfo
{
public:
    InterfaceInfo();
    ~InterfaceInfo();
public:
    InterfaceSplitBc * donor_interface = nullptr;
    SplitZone * donor_zone = nullptr;
    RangeRegion donor_region;
};

class InterfaceSplitBc
{
public:
    InterfaceSplitBc();
    ~InterfaceSplitBc();
public:
    std::vector<int> transform;
    int Mt[ 3 ][ 3 ];
    RangeRegion region;
    SplitZone * zone = nullptr;
    RangeRegion donor_region;
    SplitZone * donor_zone = nullptr;
    InterfaceSplitBc * donor_bc = nullptr;
public:
    std::vector<InterfaceSplitBc *> child;
public:
    void CopyMatrix( InterfaceSplitBc * interfaceSplitBc );
    void CalcTransformMatrix();
    void CalcSubDonorRegion( RangeRegion & subRegion, RangeRegion & subDonorRegion );
    void mapindex( std::vector<int> & begin1, std::vector<int> & begin2, std::vector<int> & index1, std::vector<int> & index2 );
    void Multiply( int Mt[][ 3 ], std::vector<int> & a, std::vector<int> & b );
    void ChangeRegionToLocalCoordinate();
};

class SplitZone
{
public:
    SplitZone();
    ~SplitZone();
public:
    DimInfo dimInfo;
    int zoneIndex, iProc;
    int newIndex;
protected:
    SplitZone * parent;
    std::vector<SplitZone *> child;
public:
    std::vector< PhysicalSplitBc * > physicalSplitBcList;
    std::vector< InterfaceSplitBc * > interfaceSplitBcList;
public:
    int GetNCell();
    void ComputeRminmax();
    void GetRootInfo( SplitZone *& root, std::vector<int> & ref );
    void SetZoneIndex( int zoneIndex ) { this->zoneIndex = zoneIndex; }
    int  GetZoneIndex() { return zoneIndex; }
    void SetProcId( int iProc ) { this->iProc = iProc; }
    void SetParent( SplitZone * parent ) { this->parent = parent; };
    void SetParentAndChild( SplitZone * parent );
    InterfaceSplitBc * FindDonorInterface( InterfaceSplitBc * interfaceSplitBc, InterfaceInfo & interfaceInfo );
    void GetChildDonorRegion( RangeRegion & subDonorRegion, InterfaceInfo & interfaceInfo );
public:
    void CreateBcFromParent();
    void CreatePhysicalBcFromParent();
    void CreatePhysicalBcFromParent( PhysicalSplitBc * parentPhysicalSplitBc );
    void CreatePhysicalBc( RangeRegion & overlapRegion, int bcType );
    void CreateInterfaceBcFromParent();
    void CreateInterfaceBcFromParent( InterfaceSplitBc * parentInterfaceSplitBc );
    void CreateInterfaceBcFromParent( InterfaceSplitBc * parentInterfaceSplitBc, RangeRegion & overlapRegion );

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

