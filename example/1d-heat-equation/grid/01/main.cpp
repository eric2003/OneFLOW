#include "cgnslib.h"
#include <iostream>
#include <vector>

void write_grid_str();
//void write_con2zn_str();

void write_grid_str()
{
    const int ni = 101;
    double x[ni];
    cgsize_t isize[9];
    int index_file,icelldim,iphysdim,index_base;
    int index_zone,index_coord;
    char basename[33],zonename[33];
    cgsize_t ipnts[ 2 ];
    int ilo,ihi;

    double x_l = -1.0;
    double x_r = 1.0;

    double dx = ( x_r - x_l ) / ( ni - 1 );

    /* create gridpoints for simple example: */
    for ( int i = 0; i < ni; ++ i )
    {
        x[ i ] = x_l + i * dx;
    }

    if (cg_open("heat1d.cgns",CG_MODE_WRITE,&index_file)) cg_error_exit();

    strcpy(basename,"Base");
    icelldim=1;
    iphysdim=3;
    cg_base_write(index_file,basename,icelldim,iphysdim,&index_base);
    /*  vertex size */
    isize[0]=ni;
    /*  cell size */
    isize[1]=ni-1;
    /*  boundary vertex size (always zero for structured grids) */
    isize[2]=0;
    isize[3]=0;
    isize[4]=0;
    isize[5]=0;
    isize[6]=0;
    isize[7]=0;
    isize[8]=0;

    ilo=1;
    ihi=isize[0];
    strcpy(zonename,"Zone 1");
    cg_zone_write(index_file,index_base,zonename,isize,CGNS_ENUMV(Structured),&index_zone);
    cg_coord_write(index_file,index_base,index_zone,CGNS_ENUMV(RealDouble),"x",x,&index_coord);

    /* lower point of range */
    ipnts[0]=ilo;
    /* upper point of range */
    ipnts[1]=ilo;

    int index_bc = -1;
    cg_boco_write(index_file,index_base,index_zone,"L",CGNS_ENUMV(BCInflow),CGNS_ENUMV(PointRange),2,ipnts,&index_bc);

    /* lower point of range */
    ipnts[0]=ihi;
    /* upper point of range */
    ipnts[1]=ihi;

    cg_boco_write(index_file,index_base,index_zone,"R",CGNS_ENUMV(BCOutflow),CGNS_ENUMV(PointRange),2,ipnts,&index_bc);

    cg_close(index_file);
}

int main(int argc, char **argv)
{
    write_grid_str();
    return 0;
}

