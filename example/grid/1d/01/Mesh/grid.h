#include "cgnslib.h"
#include <iostream>
#include <vector>
#include <string>

class BC
{
public:
    std::string bcName;
    BCType_t bctype;
};

void write_grid_str( const std::string & gridfilename, int ni, double xstart, double xend );
void write_bc_str( const std::string & gridfilename, BC &left, BC &right );
void DumpZoneBc( int index_file, int index_base, int index_zone, cgsize_t * isize, std::vector<cgsize_t> & ipnts, BC &left, BC &right );
