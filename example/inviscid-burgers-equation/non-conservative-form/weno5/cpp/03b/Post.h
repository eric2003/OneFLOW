#pragma once
#include <vector>
#include <string>

double compute_l2norm( int ni, std::vector<double> & r );
double compute_max_error( int ni, std::vector<double> & u_error );

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
public:
    void AddVector( std::vector<double> & a, std::vector<double> & b );
    void DumpField();
    void DumpErrorDetails( std::vector<double> & u_error );
    void DumpCsvFile( const std::string & filename, std::vector<double> & x, std::vector<double> & ue, std::vector<double> & un, std::vector<double> & uerror );
};