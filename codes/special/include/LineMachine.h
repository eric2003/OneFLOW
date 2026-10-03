#pragma once
#include "HXDefine.h"
#include "HXLookup.h"
#include <memory>

BeginNameSpace( ONEFLOW )

class SegmentCtrl;
class CurveInfo;
class CurveMesh;
class TextFileParser;
struct GridDistributionDefinition;

class LineMachine
{
public:
    LineMachine();
    ~LineMachine();
public:
    HXVector< std::unique_ptr< SegmentCtrl > > segmentCtrlList;
    HXVector< std::unique_ptr< CurveInfo > > curveInfoList;
    HXVector< std::unique_ptr< CurveMesh > > curveMeshList;
    IntField dimList;
    RealField ds1List, ds2List;
public:
    HXLookup<int> lineLookup;
    LinkField lineList;
public:
    void Reset();
    int AddLine( int p1, int p2 );
    void AddLine( int p1, int p2, int id );
    void AddCircle( int p1, int pc, int p2, int id );
    void AddDimension( TextFileParser & textFileParser );
    void AddDs( TextFileParser & textFileParser );
    void SetDimension( int id, int pointCount );
    void SetDistribution( const GridDistributionDefinition & distribution );
    void GenerateAllLineMesh();
    void CreateAllLineMesh();
public:
    SegmentCtrl * GetSegmentCtrl( int id ) const;
    CurveMesh * GetCurveMesh( int id ) const;
    CurveInfo * GetCurveInfo( int id ) const;
public:
    CurveMesh * GetLineMeshByTwoPoint( const int & p1, const int & p2, int & direction ) const;
    int GetLineIdByTwoPoint( const int & p1, const int & p2 ) const;
};

extern LineMachine line_Machine;

EndNameSpace
