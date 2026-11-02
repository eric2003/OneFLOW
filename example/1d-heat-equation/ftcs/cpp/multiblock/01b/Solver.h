#pragma once
#include <vector>
#include "global.h"
#include "CgnsGrid.h"
#include "Field.h"

double compute_l2norm( int ni, std::vector<double> & r );
double compute_max_error( int ni, std::vector<double> & u_error );

class Solver
{
public:
    void Run();
    void ReadGrid();
    void InitFields();
    void InitTopo();
    void SolveMultiZones();
public:
    void Boundary();
    void ExchangeInterfaceField();
    void PrintField( std::vector<double> & f );
public:
    void PostProcess();
};


class Post
{
public:
    std::vector<int> zoneids;
    std::vector<double> x;
    std::vector<double> u_e;
    std::vector<double> un;
    std::vector<double> u;
public:
    void Process();
    void ReorderZones();
    void GatherField();
    void DumpField();
    void DumpErrorDetails( std::vector<double> & u_error );
    void DumpCsvFile( const std::string & filename, std::vector<double> & x, std::vector<double> & ue, std::vector<double> & un, std::vector<double> & uerror );
    void DumpCsvFile( const std::string & filename, std::vector<double> & x, std::vector<double> & ue, std::vector<double> & un, std::vector<double> & u, std::vector<double> & uerror );
};