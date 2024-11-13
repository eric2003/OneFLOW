#pragma once
#include <vector>
#include "global.h"
#include "CgnsGrid.h"
#include "Field.h"

double compute_l2norm( int ni, std::vector<double> & r );
double compute_max_error( int ni, std::vector<double> & u_error );
void DumpErrorDetails( std::vector<double> & u_error );
void DumpCsvFile( const std::string & filename, std::vector<double> & x, std::vector<double> & ue, std::vector<double> & un, std::vector<double> & uerror );

class Solver
{
public:
    void Run();
    void ReadGrid();
    void InitFields();
    void InitTopo();
    void SolveMultiZones();
    void ExchangeInterfaceValue();
public:
    void InitInterfaceField();
    void ExchangeInterfaceField();
    void PostProcess();
    void PostProcess( std::vector<double> & x, std::vector<double> & u_e, std::vector<double> & un );
    void PrintField( std::vector<double> & f );
public:
    void PostProcessNewVersion();
};


class Post
{
public:
    std::vector<int> zoneids;
    std::vector<double> x;
    std::vector<double> u_e;
    std::vector<double> un;
public:
    void Process();
    void ReorderZones();
    void GatherField();
    void DumpField();
    void DumpErrorDetails( std::vector<double> & u_error );
    void DumpCsvFile( const std::string & filename, std::vector<double> & x, std::vector<double> & ue, std::vector<double> & un, std::vector<double> & uerror );
};