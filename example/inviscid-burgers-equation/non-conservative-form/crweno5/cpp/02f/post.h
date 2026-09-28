#pragma once
#include <vector>
#include <string>

class Vec1d;
void DumpCsvFile( const std::string & filename, std::vector<double> & x, std::vector<std::vector<double>> & u );
void DumpField( int iter, Vec1d & x, Vec1d & u );