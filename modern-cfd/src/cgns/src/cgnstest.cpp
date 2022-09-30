/*---------------------------------------------------------------------------*\
OneFLOW - LargeScale Multiphysics Scientific Simulation Environment
Copyright (C) 2017-2022 He Xin and the OneFLOW contributors.
-------------------------------------------------------------------------------
License
This file is part of OneFLOW.

OneFLOW is free software: you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

OneFLOW is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
for more details.

You should have received a copy of the GNU General Public License
along with OneFLOW.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "cgnstest.h"
#include "Cgns_t.h"
#include "cgnslib.h"
#include "tools.h"
#include <string>
#include <vector>
#include <iostream>
#include <iomanip>

#ifdef PRJ_ENABLE_CGNS

void cgns_write_base_test()
{
    if ( cg_open( "test.cgns", CG_MODE_WRITE, &Cgns_t::file_id ) )
    {
        cg_error_exit();
    }

    Cgns_t::cell_dim=1;
    Cgns_t::phys_dim=3;

    Cgns_t::WriteBase( "baseA");
    Cgns_t::WriteBase( "baseB");

    cg_close( Cgns_t::file_id );
}

int cgnstest()
{
    int index_file = -1;
    if ( cg_open( "test.cgns", CG_MODE_WRITE, &index_file ) )
    {
        cg_error_exit();
    }
    std::string basename = "Base";

    int icelldim=3;
    int iphysdim=3;
    int index_base = -1;
    cg_base_write(index_file,basename.c_str(), icelldim, iphysdim, &index_base);
    cg_close(index_file);
    return 0;
}

void cgns_read_grid( std::vector<float> & xcoor, const std::string & filename )
{
    if ( cg_open( filename.c_str(), CG_MODE_READ, &Cgns_t::file_id) )
    {
        cg_error_exit();
    }

    //Determine the of bases in the grid
    cg_nbases( Cgns_t::file_id, &Cgns_t::nbases );
    std::cout << "   Total number of CGNS Base = " << Cgns_t::nbases << "\n";

    for ( int ibase = 0; ibase < Cgns_t::nbases; ++ ibase )
    {
        Cgns_t::char33 cgns_base_name;
        Cgns_t::base_id = ibase + 1;

        double double_base_id;
        cg_base_id( Cgns_t::file_id, Cgns_t::base_id, & double_base_id );
        std::cout << "   double_base_id = " << double_base_id << "\n";

        //Check the cell and physical dimensions of the bases.
        cg_base_read( Cgns_t::file_id, Cgns_t::base_id, cgns_base_name, & Cgns_t::cell_dim, & Cgns_t::phys_dim );
        std::cout << "   Cgns_t::base_id = " << Cgns_t::base_id << " baseName = " << cgns_base_name << "\n";
        std::cout << "   cell dim = " << Cgns_t::cell_dim << " physical dim = " << Cgns_t::phys_dim << "\n";

        //Read the number of zones in the grid.
        cg_nzones( Cgns_t::file_id, Cgns_t::base_id, & Cgns_t::nzones );

        std::cout << "** Reading CGNS Grid In Base " << Cgns_t::base_id << "\n";
        std::cout << "   Reading CGNS Family Specified BC \n";
        //this->ReadFamilySpecifiedBc();
        std::cout << "   numberOfCgnsZones       = " << Cgns_t::nzones << "\n\n";

        for ( int iZone = 0; iZone < Cgns_t::nzones; ++ iZone )
        {
            Cgns_t::zone_id = iZone + 1;
            std::cout << "==>Cgns_t::zone_id = " << Cgns_t::zone_id << " numberOfCgnsZones = " << Cgns_t::nzones << "\n";

            ZoneType_t cgnsZoneType;
            cg_zone_type( Cgns_t::file_id, Cgns_t::base_id, Cgns_t::zone_id, & cgnsZoneType );
            Cgns_t::zone_type = cgnsZoneType;
            //Check the zone type
            std::cout << "   The Zone Type is " << cg_ZoneTypeName( Cgns_t::zone_type ) << " Zone" << "\n";

            Cgns_t::char33 cgns_zone_name;

            //Determine the number of vertices and cellVolume elements in this zone
            cg_zone_read( Cgns_t::file_id, Cgns_t::base_id, Cgns_t::zone_id, cgns_zone_name, Cgns_t::isize );
            Cgns_t::SetDimensionStr( Cgns_t::isize );

            std::cout << "   CGNS Zone Name = " << cgns_zone_name << "\n";

            int ncoords = -1;
            cg_ncoords ( Cgns_t::file_id, Cgns_t::base_id, Cgns_t::zone_id, &ncoords );
            std::cout << "   ncoords = " << ncoords << "\n";
            Cgns_t::char33 coor_name;
            DataType_t dataType;
            xcoor.resize( Cgns_t::nnodes );
            for ( int icoor = 1; icoor <= ncoords; ++ icoor )
            {
                Cgns_t::coor_id = icoor;
                cg_coord_info( Cgns_t::file_id, Cgns_t::base_id, Cgns_t::zone_id, Cgns_t::coor_id, & dataType, coor_name );
                std::cout << "   coor_name = " << coor_name << " dataType = " << dataType << " dataTypeName = " << cg_DataTypeName( dataType ) << "\n";
                cg_coord_read( Cgns_t::file_id, Cgns_t::base_id, Cgns_t::zone_id, coor_name, dataType, Cgns_t::irmin, Cgns_t::irmax, xcoor.data() );
            }

            ReadFieldExample();

            cg_nbocos( Cgns_t::file_id, Cgns_t::base_id, Cgns_t::zone_id, & Cgns_t::nbccos );
            std::cout << "   nBoco        = " << Cgns_t::nbccos << std::endl;

            cgsize_t npts, ndpts;
            cgsize_t range[6], d_range[6];
            int transform[3];
            CGNS_ENUMT(GridLocation_t) location;
            CGNS_ENUMT(GridConnectivityType_t) type;
            CGNS_ENUMT(PointSetType_t) ptype, d_ptype;
            CGNS_ENUMT(ZoneType_t) d_ztype;
            CGNS_ENUMT(DataType_t) datatype;
            CGNS_ENUMT(BCType_t) bctype;

            Cgns_t::char33 cgns_bc_name;
            cgsize_t normal_list_size;
            int n_data_sets;
   
            for ( int iBoco = 0; iBoco < Cgns_t::nbccos; ++ iBoco )
            {
                Cgns_t::bc_id = iBoco + 1;
                std::cout << "\n";
                std::cout << "-->Cgns_t::bc_id  = " << Cgns_t::bc_id << " Cgns_t::nbccos = " << Cgns_t::nbccos << "\n";
                double bc_double_id = -1.0;
                cg_boco_id( Cgns_t::file_id, Cgns_t::base_id, Cgns_t::zone_id, Cgns_t::bc_id, & bc_double_id );
                std::cout << "   bc_double_id              = " << bc_double_id << "\n";
                cg_boco_info ( Cgns_t::file_id, Cgns_t::base_id, Cgns_t::zone_id, Cgns_t::bc_id, cgns_bc_name,
                    &bctype, &ptype, &npts, transform, &normal_list_size, &datatype, &n_data_sets );
                std::cout << "   CGNS Boundary Name             = " << cgns_bc_name << "\n";
                std::cout << "   CGNS Boundary Condition Name   = " << cg_BCTypeName( bctype ) << "\n";
                std::cout << "   CGNS PointSet Type Name        = " << cg_PointSetTypeName(ptype) << "\n";
                GridLocation_t bc_grid_location;
                cg_boco_gridlocation_read( Cgns_t::file_id, Cgns_t::base_id, Cgns_t::zone_id, Cgns_t::bc_id, &bc_grid_location );
                std::cout << "   CGNS Boco Grid Location Name   = " << ::cg_GridLocationName( bc_grid_location ) << "\n";
                int cgnsNormalList;
                std::vector<cgsize_t> pts( 3 * npts );
                // Read the element IDs.
                cg_boco_read( Cgns_t::file_id, Cgns_t::base_id, Cgns_t::zone_id, Cgns_t::bc_id, pts.data(), 0);
                int kkk = 1;
            }

            int kkk = 1;
        }
        int kkk = 1;
    }
    cg_close( Cgns_t::file_id );
}

void cgns_read_and_modify_grid( std::vector<float> & xcoor, const std::string & filename )
{
    if ( cg_open( filename.c_str(), CG_MODE_MODIFY, &Cgns_t::file_id) )
    {
        cg_error_exit();
    }

    //Determine the of bases in the grid
    cg_nbases( Cgns_t::file_id, &Cgns_t::nbases );
    std::cout << "   Total number of CGNS Base = " << Cgns_t::nbases << "\n";

    for ( int ibase = 0; ibase < Cgns_t::nbases; ++ ibase )
    {
        Cgns_t::char33 cgns_base_name;
        Cgns_t::base_id = ibase + 1;

        double double_base_id;
        cg_base_id( Cgns_t::file_id, Cgns_t::base_id, & double_base_id );
        std::cout << "   double_base_id = " << double_base_id << "\n";

        //Check the cell and physical dimensions of the bases.
        cg_base_read( Cgns_t::file_id, Cgns_t::base_id, cgns_base_name, & Cgns_t::cell_dim, & Cgns_t::phys_dim );
        std::cout << "   Cgns_t::base_id = " << Cgns_t::base_id << " baseName = " << cgns_base_name << "\n";
        std::cout << "   cell dim = " << Cgns_t::cell_dim << " physical dim = " << Cgns_t::phys_dim << "\n";

        //Read the number of zones in the grid.
        cg_nzones( Cgns_t::file_id, Cgns_t::base_id, & Cgns_t::nzones );

        std::cout << "** Reading CGNS Grid In Base " << Cgns_t::base_id << "\n";
        std::cout << "   Reading CGNS Family Specified BC \n";
        //this->ReadFamilySpecifiedBc();
        std::cout << "   numberOfCgnsZones       = " << Cgns_t::nzones << "\n\n";

        for ( int iZone = 0; iZone < Cgns_t::nzones; ++ iZone )
        {
            Cgns_t::zone_id = iZone + 1;
            std::cout << "==>Cgns_t::zone_id = " << Cgns_t::zone_id << " numberOfCgnsZones = " << Cgns_t::nzones << "\n";

            ZoneType_t cgnsZoneType;
            cg_zone_type( Cgns_t::file_id, Cgns_t::base_id, Cgns_t::zone_id, & cgnsZoneType );
            Cgns_t::zone_type = cgnsZoneType;
            //Check the zone type
            std::cout << "   The Zone Type is " << cg_ZoneTypeName( Cgns_t::zone_type ) << " Zone" << "\n";

            Cgns_t::char33 cgns_zone_name;

            //Determine the number of vertices and cellVolume elements in this zone
            cg_zone_read( Cgns_t::file_id, Cgns_t::base_id, Cgns_t::zone_id, cgns_zone_name, Cgns_t::isize );
            Cgns_t::SetDimensionStr( Cgns_t::isize );

            std::cout << "   CGNS Zone Name = " << cgns_zone_name << "\n";

            int ncoords = -1;
            cg_ncoords ( Cgns_t::file_id, Cgns_t::base_id, Cgns_t::zone_id, &ncoords );
            std::cout << "   ncoords = " << ncoords << "\n";
            Cgns_t::char33 coor_name;
            DataType_t dataType;
            xcoor.resize( Cgns_t::nnodes );
            for ( int icoor = 1; icoor <= ncoords; ++ icoor )
            {
                Cgns_t::coor_id = icoor;
                cg_coord_info( Cgns_t::file_id, Cgns_t::base_id, Cgns_t::zone_id, Cgns_t::coor_id, & dataType, coor_name );
                std::cout << "   coor_name = " << coor_name << " dataType = " << dataType << " dataTypeName = " << cg_DataTypeName( dataType ) << "\n";
                cg_coord_read( Cgns_t::file_id, Cgns_t::base_id, Cgns_t::zone_id, coor_name, dataType, Cgns_t::irmin, Cgns_t::irmax, xcoor.data() );
            }

            ModifyFieldExample();

            cg_nbocos( Cgns_t::file_id, Cgns_t::base_id, Cgns_t::zone_id, & Cgns_t::nbccos );
            std::cout << "   nBoco        = " << Cgns_t::nbccos << std::endl;

            cgsize_t npts, ndpts;
            cgsize_t range[6], d_range[6];
            int transform[3];
            CGNS_ENUMT(GridLocation_t) location;
            CGNS_ENUMT(GridConnectivityType_t) type;
            CGNS_ENUMT(PointSetType_t) ptype, d_ptype;
            CGNS_ENUMT(ZoneType_t) d_ztype;
            CGNS_ENUMT(DataType_t) datatype;
            CGNS_ENUMT(BCType_t) bctype;

            Cgns_t::char33 cgns_bc_name;
            cgsize_t normal_list_size;
            int n_data_sets;

            for ( int iBoco = 0; iBoco < Cgns_t::nbccos; ++ iBoco )
            {
                Cgns_t::bc_id = iBoco + 1;
                std::cout << "\n";
                std::cout << "-->Cgns_t::bc_id  = " << Cgns_t::bc_id << " Cgns_t::nbccos = " << Cgns_t::nbccos << "\n";
                double bc_double_id = -1.0;
                cg_boco_id( Cgns_t::file_id, Cgns_t::base_id, Cgns_t::zone_id, Cgns_t::bc_id, & bc_double_id );
                std::cout << "   bc_double_id              = " << bc_double_id << "\n";
                cg_boco_info ( Cgns_t::file_id, Cgns_t::base_id, Cgns_t::zone_id, Cgns_t::bc_id, cgns_bc_name,
                    &bctype, &ptype, &npts, transform, &normal_list_size, &datatype, &n_data_sets );
                std::cout << "   CGNS Boundary Name             = " << cgns_bc_name << "\n";
                std::cout << "   CGNS Boundary Condition Name   = " << cg_BCTypeName( bctype ) << "\n";
                std::cout << "   CGNS PointSet Type Name        = " << cg_PointSetTypeName(ptype) << "\n";
                GridLocation_t bc_grid_location;
                cg_boco_gridlocation_read( Cgns_t::file_id, Cgns_t::base_id, Cgns_t::zone_id, Cgns_t::bc_id, &bc_grid_location );
                std::cout << "   CGNS Boco Grid Location Name   = " << ::cg_GridLocationName( bc_grid_location ) << "\n";
                int cgnsNormalList;
                std::vector<cgsize_t> pts( 3 * npts );
                // Read the element IDs.
                cg_boco_read( Cgns_t::file_id, Cgns_t::base_id, Cgns_t::zone_id, Cgns_t::bc_id, pts.data(), 0);
                int kkk = 1;
            }

            int kkk = 1;
        }
        int kkk = 1;
    }
    cg_close( Cgns_t::file_id );
}

void cgns_dump_grid( float * xcoor, int ni, const std::string & filename )
{
    if ( cg_open( filename.c_str(), CG_MODE_WRITE, &Cgns_t::file_id) )
    {
        cg_error_exit();
    }

    Cgns_t::cell_dim = 1;
    Cgns_t::phys_dim = 3;
    cg_base_write( Cgns_t::file_id, "Base", Cgns_t::cell_dim, Cgns_t::phys_dim, &Cgns_t::base_id );
    std::vector<int> dim_array;
    cgsize_t isize[ 9 ];
    Cgns_t::SetDimArray( dim_array, ni );
    Cgns_t::SetISize( isize, dim_array );
    cg_zone_write( Cgns_t::file_id, Cgns_t::base_id, "Zone 1", isize, CGNS_ENUMV(Structured), &Cgns_t::zone_id );
    cgsize_t ipnts[ 6 ];
    {
        Cgns_t::SetStructuredIpnts( ipnts, 1 );
        cg_boco_write( Cgns_t::file_id, Cgns_t::base_id, Cgns_t::zone_id, "left", CGNS_ENUMV( BCInflow ),
            CGNS_ENUMV( PointRange ), 2, ipnts, &Cgns_t::bc_id );
    }
    {
        Cgns_t::SetStructuredIpnts( ipnts, ni );
        cg_boco_write( Cgns_t::file_id, Cgns_t::base_id, Cgns_t::zone_id, "right", CGNS_ENUMV( BCOutflow ),
            CGNS_ENUMV( PointRange ), 2, ipnts, &Cgns_t::bc_id );
    }
    cg_coord_write ( Cgns_t::file_id, Cgns_t::base_id, Cgns_t::zone_id, CGNS_ENUMV(RealSingle), "CoordinateX", xcoor, &Cgns_t::coor_id );
    WriteFieldExample( ni );
    cg_close( Cgns_t::file_id );
}

void ReadFieldExample()
{
    cg_nsols( Cgns_t::file_id, Cgns_t::base_id, Cgns_t::zone_id, & Cgns_t::nsols );
    std::cout << "   Cgns_t::nsols       = " << Cgns_t::nsols << "\n";

    for ( int isol = 0; isol < Cgns_t::nsols; ++ isol )
    {
        Cgns_t::char33 sol_name;
        GridLocation_t location;
        double double_sol_id = -1;
        Cgns_t::sols_id = isol + 1;
        cg_sol_info ( Cgns_t::file_id, Cgns_t::base_id, Cgns_t::zone_id, Cgns_t::sols_id, sol_name, &location );
        std::cout << "   Cgns_t::sols_id  = " << Cgns_t::sols_id << " sol_name = " <<  sol_name << \
            " location = " << location << " location name = " << cg_GridLocationName( location ) << "\n";
        cg_sol_id( Cgns_t::file_id, Cgns_t::base_id, Cgns_t::zone_id, Cgns_t::sols_id, & double_sol_id );
        std::cout << "   double_sol_id  = " << double_sol_id << "\n";
        cg_nfields( Cgns_t::file_id, Cgns_t::base_id, Cgns_t::zone_id, Cgns_t::sols_id, & Cgns_t::nfields );
        std::cout << "   Cgns_t::nfields  = " << Cgns_t::nfields << "\n";
        for ( int ifield = 0; ifield < Cgns_t::nfields; ++ ifield )
        {
            std::vector<float> q_field( Cgns_t::nnodes, -1.0 );
            Cgns_t::field_id = ifield + 1;
            CGNS_ENUMT( DataType_t ) dataType;
            Cgns_t::char33 field_name;
            double double_field_id = -1;
            cg_field_info( Cgns_t::file_id, Cgns_t::base_id, Cgns_t::zone_id, Cgns_t::sols_id, Cgns_t::field_id, \
                & dataType, field_name);
            std::cout << "   Cgns_t::field_id  = " << Cgns_t::field_id << " field_name = " <<  field_name << \
                " dataType = " << dataType << " dataType name = " << cg_DataTypeName( dataType ) << "\n";
            cg_field_id( Cgns_t::file_id, Cgns_t::base_id, Cgns_t::zone_id, Cgns_t::sols_id, Cgns_t::field_id, & double_field_id );
            std::cout << "   double_field_id  = " << double_field_id << "\n";
            for ( int idim = 0; idim < 3; ++ idim )
            {
                std::cout << "   Cgns_t::irmin[" << idim << "] = " << Cgns_t::irmin[ idim ];
                std::cout << "   Cgns_t::irmax[" << idim << "] = " << Cgns_t::irmax[ idim ] << "\n";
            }
            cg_field_read ( Cgns_t::file_id, Cgns_t::base_id, Cgns_t::zone_id, Cgns_t::sols_id, field_name, dataType, Cgns_t::irmin, Cgns_t::irmax, q_field.data() );
        }
    }
}

void ModifyFieldExample()
{
    cg_nsols( Cgns_t::file_id, Cgns_t::base_id, Cgns_t::zone_id, & Cgns_t::nsols );
    std::cout << "   Cgns_t::nsols       = " << Cgns_t::nsols << "\n";

    for ( int isol = 0; isol < Cgns_t::nsols; ++ isol )
    {
        Cgns_t::char33 sol_name;
        GridLocation_t location;
        double double_sol_id = -1;
        Cgns_t::sols_id = isol + 1;
        cg_sol_info ( Cgns_t::file_id, Cgns_t::base_id, Cgns_t::zone_id, Cgns_t::sols_id, sol_name, &location );
        std::cout << "   Cgns_t::sols_id  = " << Cgns_t::sols_id << " sol_name = " <<  sol_name << \
            " location = " << location << " location name = " << cg_GridLocationName( location ) << "\n";
        cg_sol_id( Cgns_t::file_id, Cgns_t::base_id, Cgns_t::zone_id, Cgns_t::sols_id, & double_sol_id );
        std::cout << "   double_sol_id  = " << double_sol_id << "\n";
        cg_nfields( Cgns_t::file_id, Cgns_t::base_id, Cgns_t::zone_id, Cgns_t::sols_id, & Cgns_t::nfields );
        std::cout << "   Cgns_t::nfields  = " << Cgns_t::nfields << "\n";
        for ( int ifield = 0; ifield < Cgns_t::nfields; ++ ifield )
        {
            std::vector<float> q_field( Cgns_t::nnodes, -1.0 );
            Cgns_t::field_id = ifield + 1;
            CGNS_ENUMT( DataType_t ) dataType;
            Cgns_t::char33 field_name;
            double double_field_id = -1;
            cg_field_info( Cgns_t::file_id, Cgns_t::base_id, Cgns_t::zone_id, Cgns_t::sols_id, Cgns_t::field_id, \
                & dataType, field_name);
            std::cout << "   Cgns_t::field_id  = " << Cgns_t::field_id << " field_name = " <<  field_name << \
                " dataType = " << dataType << " dataType name = " << cg_DataTypeName( dataType ) << "\n";
            cg_field_id( Cgns_t::file_id, Cgns_t::base_id, Cgns_t::zone_id, Cgns_t::sols_id, Cgns_t::field_id, & double_field_id );
            std::cout << "   double_field_id  = " << double_field_id << "\n";
            for ( int idim = 0; idim < 3; ++ idim )
            {
                std::cout << "   Cgns_t::irmin[" << idim << "] = " << Cgns_t::irmin[ idim ];
                std::cout << "   Cgns_t::irmax[" << idim << "] = " << Cgns_t::irmax[ idim ] << "\n";
            }
            cg_field_read ( Cgns_t::file_id, Cgns_t::base_id, Cgns_t::zone_id, Cgns_t::sols_id, field_name, dataType, Cgns_t::irmin, Cgns_t::irmax, q_field.data() );
            std::vector<float> q_tmp;
            q_tmp = q_field;
            static int ii_count_tmp = 0;
            for ( int i = 0; i < q_tmp.size(); ++ i )
            {
                q_tmp[ i ] = ii_count_tmp + i;
            }
            ii_count_tmp += 100;
            //std::cout << alignl( 10, "abcd" ) << ":" << " haha\n";
            //std::cout << alignl( 10, "abcde" ) << ":" << " haha\n";
            //std::cout << alignl( 10, "abcdefgh" ) << ":" << " haha\n";
            //std::cout << alignl( 10, "abcdefghijklnm" ) << ":" << " haha\n";

            //std::cout << alignr( 10, "abcd" ) << ":" << " haha\n";
            //std::cout << alignr( 10, "abcde" ) << ":" << " haha\n";
            //std::cout << alignr( 10, "abcdefgh" ) << ":" << " haha\n";
            //std::cout << alignr( 10, "abcdefghijklnm" ) << ":" << " haha\n";
            int ifld = -1;
            cg_field_write( Cgns_t::file_id, Cgns_t::base_id, Cgns_t::zone_id, Cgns_t::sols_id, dataType, field_name, q_tmp.data(), &ifld);
        }
    }
}


void WriteFieldExample( int ni )
{
    {
        CGNS_ENUMT(DataType_t) datatype = CGNS_ENUMV(RealSingle);
        int isol, ifld, nv;
        cg_sol_write(Cgns_t::file_id, Cgns_t::base_id, Cgns_t::zone_id, "NodeVariables", CGNS_ENUMV(Vertex), &isol );
        std::cout << "isol=" << isol << "\n";
        std::vector<float> q( ni, 0.0 );
        cg_field_write( Cgns_t::file_id, Cgns_t::base_id, Cgns_t::zone_id, isol, datatype, "q", q.data(), &ifld);
        std::cout << "ifld=" << ifld << "\n";
        std::vector<float> u( ni, 1.0 );
        cg_field_write( Cgns_t::file_id, Cgns_t::base_id, Cgns_t::zone_id, isol, datatype, "u", u.data(), &ifld);
        std::cout << "ifld=" << ifld << "\n";
    }
    {
        CGNS_ENUMT(DataType_t) datatype = CGNS_ENUMV(RealSingle);
        int isol, ifld, nv;
        cg_sol_write(Cgns_t::file_id, Cgns_t::base_id, Cgns_t::zone_id, "CFD-FlowField", CGNS_ENUMV(Vertex), &isol );
        std::cout << "isol=" << isol << "\n";
        std::vector<float> q( ni, 0.0 );
        cg_field_write( Cgns_t::file_id, Cgns_t::base_id, Cgns_t::zone_id, isol, datatype, "q", q.data(), &ifld);
        std::cout << "ifld=" << ifld << "\n";
        std::vector<float> u( ni, 1.0 );
        cg_field_write( Cgns_t::file_id, Cgns_t::base_id, Cgns_t::zone_id, isol, datatype, "u", u.data(), &ifld);
        std::cout << "ifld=" << ifld << "\n";
        std::vector<float> v( ni, 2.0 );
        cg_field_write( Cgns_t::file_id, Cgns_t::base_id, Cgns_t::zone_id, isol, datatype, "v", v.data(), &ifld);
        std::cout << "ifld=" << ifld << "\n";
    }
}


void WriteZoneInfo( const std::string & zoneName, ZoneType_t zoneType, cgsize_t * isize )
{
    cg_zone_write( Cgns_t::file_id, Cgns_t::base_id, zoneName.c_str(), isize, zoneType, & Cgns_t::zone_id );
}


void TestCgnsLink()
{
    if ( cg_open( "testlink.cgns", CG_MODE_WRITE, &Cgns_t::file_id ) )
    {
        cg_error_exit();
    }

    Cgns_t::cell_dim = 1;
    Cgns_t::phys_dim = 3;

    Cgns_t::WriteBase( "Base");

    std::vector<int> dim_array;
    cgsize_t isize[ 9 ];
    Cgns_t::SetDimArray( dim_array, 10 );
    Cgns_t::SetISize( isize, dim_array );

    ::WriteZoneInfo( "Zone 1", CGNS_ENUMV( Structured ), isize );
    Cgns_t::SetDimArray( dim_array, 20 );
    Cgns_t::SetISize( isize, dim_array );
    ::WriteZoneInfo( "Zone 2", CGNS_ENUMV( Structured ), isize );

    cg_close( Cgns_t::file_id );
}

void cgns_dump_grid()
{
    if ( cg_open( "oneflow-2d.cgns", CG_MODE_WRITE, &Cgns_t::file_id ) )
    {
        cg_error_exit();
    }

    Cgns_t::cell_dim = 2;
    Cgns_t::phys_dim = 3;
    cg_base_write( Cgns_t::file_id, "Base", Cgns_t::cell_dim, Cgns_t::phys_dim, &Cgns_t::base_id );
    std::vector<int> dim_array;
    cgsize_t isize[ 9 ];
    int ni = 11;
    int nj = 2;
    Cgns_t::SetDimArray( dim_array, ni, nj );
    Cgns_t::SetISize( isize, dim_array );
    cg_zone_write( Cgns_t::file_id, Cgns_t::base_id, "Zone 1", isize, CGNS_ENUMV(Structured), &Cgns_t::zone_id );
    cgsize_t ipnts[ 6 ];
    {
        Cgns_t::SetStructuredIpnts( ipnts, 1, 1, 1, nj );
        cg_boco_write( Cgns_t::file_id, Cgns_t::base_id, Cgns_t::zone_id, "left", CGNS_ENUMV( BCInflow ),
            CGNS_ENUMV( PointRange ), 2, ipnts, &Cgns_t::bc_id );
    }
    {
        Cgns_t::SetStructuredIpnts( ipnts, ni, 1, ni, nj );
        cg_boco_write( Cgns_t::file_id, Cgns_t::base_id, Cgns_t::zone_id, "right", CGNS_ENUMV( BCOutflow ),
            CGNS_ENUMV( PointRange ), 2, ipnts, &Cgns_t::bc_id );
    }
    {
        Cgns_t::SetStructuredIpnts( ipnts, 1, 1, ni, 1 );
        cg_boco_write( Cgns_t::file_id, Cgns_t::base_id, Cgns_t::zone_id, "bottom", CGNS_ENUMV( BCSymmetryPlane ),
            CGNS_ENUMV( PointRange ), 2, ipnts, &Cgns_t::bc_id );
    }
    {
        Cgns_t::SetStructuredIpnts( ipnts, 1, nj, ni, nj );
        cg_boco_write( Cgns_t::file_id, Cgns_t::base_id, Cgns_t::zone_id, "top", CGNS_ENUMV( BCSymmetryPlane ),
            CGNS_ENUMV( PointRange ), 2, ipnts, &Cgns_t::bc_id );
    }


    cg_close( Cgns_t::file_id );

}

void cgns_dump_grid_1d()
{
    if ( cg_open( "oneflow-1d.cgns", CG_MODE_WRITE, &Cgns_t::file_id ) )
    {
        cg_error_exit();
    }

    Cgns_t::cell_dim = 1;
    Cgns_t::phys_dim = 3;
    cg_base_write( Cgns_t::file_id, "Base", Cgns_t::cell_dim, Cgns_t::phys_dim, &Cgns_t::base_id );
    std::vector<int> dim_array;
    cgsize_t isize[ 9 ];
    int ni = 11;
    Cgns_t::SetDimArray( dim_array, ni );
    Cgns_t::SetISize( isize, dim_array );
    cg_zone_write( Cgns_t::file_id, Cgns_t::base_id, "Zone 1", isize, CGNS_ENUMV(Structured), &Cgns_t::zone_id );
    cgsize_t ipnts[ 6 ];
    {
        Cgns_t::SetStructuredIpnts( ipnts, 1 );
        cg_boco_write( Cgns_t::file_id, Cgns_t::base_id, Cgns_t::zone_id, "left", CGNS_ENUMV( BCInflow ),
            CGNS_ENUMV( PointRange ), 2, ipnts, &Cgns_t::bc_id );
    }
    {
        Cgns_t::SetStructuredIpnts( ipnts, ni );
        cg_boco_write( Cgns_t::file_id, Cgns_t::base_id, Cgns_t::zone_id, "right", CGNS_ENUMV( BCOutflow ),
            CGNS_ENUMV( PointRange ), 2, ipnts, &Cgns_t::bc_id );
    }

    cg_close( Cgns_t::file_id );

}

#endif