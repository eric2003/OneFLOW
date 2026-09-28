#include "grid.h"

void write_grid_str( const std::string & gridfilename, int ni, double xstart, double xend )
{
    int index_file, icelldim, iphysdim, index_base;
    int index_zone, index_coord;
    char basename[ 33 ], zonename[ 33 ];

    // A++++++++B
    // o-------->i
    // x_l      x_r
    // 1        ni
    std::vector<double> x(ni);

    double x_l = xstart;
    double x_r = xend;

    double dx = ( x_r - x_l ) / ( ni - 1 );

    for ( int i = 0; i < ni; ++ i )
    {
        x[ i ] = x_l + i * dx;
    }

    if ( cg_open( gridfilename.c_str(), CG_MODE_WRITE, &index_file ) ) cg_error_exit();
    /*  create base (user can give any name) */
    strcpy_s( basename, "Base" );
    icelldim = 1;
    iphysdim = 1;

    std::vector<cgsize_t> isize( icelldim * 3, 0 );

    cg_base_write( index_file, basename, icelldim, iphysdim, &index_base );
    /*  vertex size */
    isize[ 0 ] = ni;
    /*  cell size */
    isize[ 1 ] = ni - 1;
    /*  boundary vertex size (always zero for structured grids) */
    isize[ 2 ] = 0;

    /*  define zone 1 name (user can give any name) */
    strcpy_s( zonename, "Zone1" );
    /*  create zone */
    cg_zone_write( index_file, index_base, zonename, isize.data(), CGNS_ENUMV( Structured ), &index_zone );
    /*  write grid coordinates (user must use SIDS-standard names here) */
    cg_coord_write( index_file, index_base, index_zone, CGNS_ENUMV( RealDouble ), "CoordinateX", x.data(), &index_coord );
        
    cg_close(index_file);
}

void write_bc_str( const std::string & gridfilename, BC &left, BC &right )
{
    std::printf( "\nProgram write_bc_str\n" );

    int index_file = -1;
    if ( cg_open( gridfilename.c_str(), CG_MODE_MODIFY, &index_file ) ) cg_error_exit();

    int index_base = 1;
    int icelldim = -1;
    int iphysdim = -1;
    char basename[ 33 ];
    cg_base_read( index_file, index_base, basename, &icelldim, &iphysdim );

    std::cout << "icelldim = " << icelldim << " iphysdim = " << iphysdim << "\n";

    /* get number of zones (should be 1 for our case) */
    int nzone = -1;
    cg_nzones( index_file, index_base, &nzone );
    if ( nzone != 1 )
    {
        std::printf( "\nError. This program expects 1 zones. %d read.\n", nzone );
        return;
    }

    std::vector<cgsize_t> isize( icelldim * 3, 0 );
    std::vector<cgsize_t> ipnts( 2 * icelldim );

    char zonename[ 33 ];

    int index_zone = 1;
    cg_zone_read( index_file, index_base, index_zone, zonename, isize.data() );

    DumpZoneBc( index_file, index_base, index_zone, isize.data(), ipnts, left, right );

    cg_close(index_file);
}

void DumpZoneBc( int index_file, int index_base, int index_zone, cgsize_t *isize, std::vector<cgsize_t> &ipnts, BC &left, BC &right )
{
    int ilo = 1;
    int ihi = isize[ 0 ];
    int index_bc = -1;

    /* lower point of range */
    ipnts[ 0 ] = ilo;
    /* upper point of range */
    ipnts[ 1 ] = ilo;

    cg_boco_write( index_file, index_base, index_zone, left.bcName.c_str(), left.bctype, CGNS_ENUMV( PointRange ), 2, ipnts.data(), &index_bc );

    /* lower point of range */
    ipnts[ 0 ] = ihi;
    /* upper point of range */
    ipnts[ 1 ] = ihi;

    cg_boco_write( index_file, index_base, index_zone, right.bcName.c_str(), right.bctype, CGNS_ENUMV( PointRange ), 2, ipnts.data(), &index_bc );
}


