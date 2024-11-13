#include "CgnsGrid.h"
#include "cgnslib.h"
#include "global.h"
#include <iostream>
#include <vector>

Coor::Coor()
{
    ;
}

Coor::~Coor()
{
}

void Coor::DumpCoor()
{
    double * xd = reinterpret_cast<double *>( const_cast<char *>( coord.data() ) );
    for ( int i = 0; i < this->nNodes; ++ i )
    {
        //std::cout << coord[i] << " ";
        std::cout << xd[ i ] << " ";
        if ( ( i + 1 ) % 5 == 0 ) std::cout << "\n";
    }
    std::cout << "\n";
}

void Coor::DumpCoorX( std::vector<double> &x )
{
    double * xd = reinterpret_cast<double *>( const_cast<char *>( coord.data() ) );
    for ( int i = 0; i < this->nNodes; ++ i )
    {
        x[ i ] = xd[ i ];
    }
}

ZoneBcPatch::ZoneBcPatch()
{
    ;
}

ZoneBcPatch::~ZoneBcPatch()
{
}

ZoneBc::ZoneBc()
{
    ;
}

ZoneBc::~ZoneBc()
{
    delete this->patch;
}

ZoneBc1To1::ZoneBc1To1()
{
    ;
}

ZoneBc1To1::~ZoneBc1To1()
{
}

Zone::Zone()
{
    ;
}

Zone::~Zone()
{
    for ( int i = 0; i < bccos.size(); ++ i )
    {
        delete bccos[ i ];
    }

    for ( int i = 0; i < coors.size(); ++ i )
    {
        delete coors[ i ];
    }
}

void ReadCgnsGridBaseZone( const std::string & filename )
{
    int fileId = -1;
    cg_open( filename.c_str(), CG_MODE_READ, &fileId);
    std::cout << "fileId = " << fileId << "\n";
    int nbases = -1;
    cg_nbases( fileId, &nbases );
    std::cout << "nbases = " << nbases << "\n";

    for ( int iBase = 0; iBase < nbases; ++ iBase )
    {
        char basename[ 33 ];
        int baseId = iBase + 1;
        int icelldim = -1;
        int iphysdim = -1;
        cg_base_read( fileId, baseId, basename, &icelldim, &iphysdim );
        std::cout << "icelldim = " << icelldim << " iphysdim = " << iphysdim << "\n";
        int nzones = -1;
        cg_nzones( fileId, baseId, &nzones );
        std::cout << "nzones = " << nzones << "\n";
        for ( int iZone = 0; iZone < nzones; ++ iZone )
        {
            int zoneId = iZone + 1;
            int index_dim = -1;
            cg_index_dim( fileId, baseId, zoneId, &index_dim );
            std::cout << "index_dim = " << index_dim << "\n";

            std::vector<cgsize_t> isize( index_dim * 3 );

            char zonename[33];
            cg_zone_read( fileId, baseId, zoneId, zonename, isize.data() );
            for ( int i = 0; i < isize.size(); ++ i )
            {
                std::cout << "i = " << i << " isize = " << isize[i] << "\n";
            }

            BaseZone baseZone;
            baseZone.zone_name = zonename;

            Global::zone_names.AddBaseZone( baseZone );
        }
    }
    cg_close( fileId );
}

void ReadCgnsGrid( const std::string & filename )
{
    int fileId = -1;
    cg_open( filename.c_str(), CG_MODE_READ, &fileId );
    std::cout << "fileId = " << fileId << "\n";
    int nbases = -1;
    cg_nbases( fileId, &nbases );
    std::cout << "nbases = " << nbases << "\n";

    for ( int iBase = 0; iBase < nbases; ++ iBase )
    {
        char basename[ 33 ];
        int baseId = iBase + 1;
        int icelldim = -1;
        int iphysdim = -1;
        cg_base_read( fileId, baseId, basename, &icelldim, &iphysdim );
        std::cout << "icelldim = " << icelldim << " iphysdim = " << iphysdim << "\n";
        Global::cell_dim = icelldim;
        Global::phys_dim = iphysdim;
        int nzones = -1;
        cg_nzones( fileId, baseId, &nzones );
        std::cout << "nzones = " << nzones << "\n";
        for ( int iZone = 0; iZone < nzones; ++ iZone )
        {
            int zoneId = iZone + 1;
            int index_dim = -1;
            cg_index_dim( fileId, baseId, zoneId, &index_dim );
            std::cout << "index_dim = " << index_dim << "\n";

            std::vector<cgsize_t> isize( index_dim * 3 );

            char zonename[33];
            cg_zone_read( fileId, baseId, zoneId, zonename, isize.data() );
            for ( int i = 0; i < isize.size(); ++ i )
            {
                std::cout << "i = " << i << " isize = " << isize[i] << "\n";
            }

            Zone * zone = new Zone();
            Global::zones.push_back( zone );

            std::vector<cgsize_t> irmin(index_dim);
            std::vector<cgsize_t> irmax(index_dim);
            int nNodes = 1;
            for ( int m = 0; m < index_dim; ++ m )
            {
                /* lower range index */
                irmin[ m ] = 1;
                /* upper range index of vertices */
                irmax[ m ] = isize[ m ];
                nNodes *= irmax[ m ];
                zone->nijk.push_back( isize[ m ] );
            }
            std::cout << "nNodes = " << nNodes << "\n";

            ZoneType_t zoneType;
            cg_zone_type( fileId, baseId, zoneId, &zoneType );
            std::cout << "zoneType = " << zoneType << " ZoneTypeName = " << ZoneTypeName[ zoneType ] << "\n";
            int ncoords = -1;
            cg_ncoords( fileId, baseId, zoneId, &ncoords );
            std::cout << "ncoords = " << ncoords << "\n";

            Grid * grid = new Grid();
            Global::grids.push_back( grid );

            BC * bc = new BC();
            Global::bcs.push_back( bc );

            InterFaceZone * interfacezone = new InterFaceZone();
            Global::interfacezones.push_back( interfacezone );

            BaseZone baseZone;
            baseZone.zone_name = zonename;

            int gZoneId = Global::zone_names.FindBaseZone( baseZone );

            for ( int icoord = 0; icoord < ncoords; ++ icoord )
            {
                int coorId = icoord + 1;
                DataType_t dataType;
                char coordname[33];
                cg_coord_info( fileId, baseId, zoneId, coorId, &dataType, coordname );
                std::cout << "coordname = " << coordname << "\n";
                std::cout << "dataType = " << dataType << " DataTypeName = " << DataTypeName[ dataType ] << "\n";
                //std::vector<char> coord( nNodes * sizeof(double) );
                grid->x.resize( nNodes );

                Coor * coor = new Coor();
                zone->coors.push_back( coor );
                coor->coorname = coordname;
                coor->nNodes = nNodes;
                coor->coord.resize( nNodes * sizeof( double ) );

                cg_coord_read( fileId, baseId, zoneId, coordname, dataType, irmin.data(), irmax.data(), coor->coord.data() );

                coor->DumpCoor();
                if ( icoord == 0 )
                {
                    coor->DumpCoorX( grid->x );
                }
            }

            int nbocos = -1;

            cg_nbocos( fileId, baseId, zoneId, &nbocos );
            std::cout << "nbocos = " << nbocos << "\n";
            for ( int iboco = 0; iboco < nbocos; ++ iboco )
            {
                int bccoId = iboco + 1;
                GridLocation_t location;
                cg_boco_gridlocation_read( fileId, baseId, zoneId, bccoId, &location );
                std::cout << "iboco = " << iboco <<  " location = " << location << " GridLocationName = " << GridLocationName[location] << "\n";

                ZoneBc * zonebc = new ZoneBc();
                zone->bccos.push_back( zonebc );

                char boconame[ 33 ];
                BCType_t bocotype;
                PointSetType_t ptset_type;
                cgsize_t npnts = 0;
                std::vector<int> normalIndex( index_dim, -1 );
                cgsize_t normalListSize = 0;
                DataType_t normalDataType;
                int ndataset = -1;

                cg_boco_info( fileId, baseId, zoneId, bccoId, boconame, &bocotype, &ptset_type,
                    &npnts, normalIndex.data(), &normalListSize, &normalDataType, &ndataset );
                zonebc->bcType = bocotype;
                std::cout << "boconame = " << boconame << " bocotype = " << bocotype << " BCTypeName = " << BCTypeName[ bocotype ] << "\n";
                std::cout << "ptset_type = " << ptset_type <<  " PointSetTypeName = " << PointSetTypeName[ptset_type] << "\n";
                std::cout << "npnts = " << npnts << "\n";
                std::cout << "normalIndex = ";
                for ( int i = 0; i < index_dim; ++ i )
                {
                    std::cout << normalIndex[ i ] << " ";
                }
                std::cout << "\n";
                std::cout << "normalListSize = " << normalListSize << "\n";
                std::cout << "normalDataType = " << normalDataType << " DataTypeName = " << DataTypeName[ normalDataType ] << "\n";
                std::cout << "ndataset = " << ndataset << "\n";

                std::vector<char> normalList( nNodes * iphysdim * sizeof( double ) );

                std::vector<cgsize_t> pnts( npnts * index_dim );
                cg_boco_read( fileId, baseId, zoneId, bccoId, pnts.data(), normalList.data() );
                std::cout << "pnts = ";
                for ( int i = 0; i < pnts.size(); ++ i )
                {
                    std::cout << pnts[ i ] << " ";
                    zonebc->pnts.push_back( pnts[ i ] );
                }
                std::cout << "\n";
                for ( int i = 0; i < pnts.size(); ++ i )
                {
                    Global::boundary_points.AddPoint( grid->x[ pnts[ i ] - 1 ] );
                }

                if ( index_dim == 1 )
                {
                    int p1 = pnts[ 0 ];
                    double x1 = grid->x[ p1 - 1 ];
                    int id = Global::boundary_points.AddPoint( x1 );
                    Face face;
                    std::vector<int> nodes;
                    nodes.push_back( id );
                    face.AddFace( nodes );
                    Global::boundary_faces.AddFace( face );

                    bc->bctypes.push_back( bocotype );
                    bc->faceids.push_back( p1 - 1 );
                }

                double * normal_d = reinterpret_cast<double *>( const_cast<char *>( normalList.data() ) );
            }
            int n1to1 = -1;
            cg_n1to1( fileId, baseId, zoneId, &n1to1 );
            std::cout << "n1to1 = " << n1to1 << "\n";
            for ( int i1to1 = 0; i1to1 < n1to1; ++ i1to1 )
            {
                int i1to1Id = i1to1 + 1;

                ZoneBc1To1 * zonebc_1to1 = new ZoneBc1To1();
                zone->bc1to1s.push_back( zonebc_1to1 );

                char connectname[ 33 ];
                char donorname[ 33 ];
                cgsize_t npnts = 2;
                std::vector<cgsize_t> range( npnts * index_dim );
                std::vector<cgsize_t> donor_range( npnts * index_dim );
                std::vector<int> transform( index_dim );
                cg_1to1_read( fileId, baseId, zoneId, i1to1Id, connectname, donorname, range.data(), donor_range.data(), transform.data() );
                std::cout << "connectname = " << connectname << "\n";
                std::cout << "donorname = " << donorname << "\n";

                std::cout << "range = \n";
                for ( int i = 0; i < range.size(); ++ i )
                {
                    std::cout << range[ i ] << " ";
                }
                std::cout << "\n";
                std::cout << "donor_range = \n";
                for ( int i = 0; i < donor_range.size(); ++ i )
                {
                    std::cout << donor_range[ i ] << " ";
                }
                std::cout << "\n";

                for ( int i = 0; i < range.size(); ++ i )
                {
                    Global::boundary_points.AddPoint( grid->x[ range[ i ] - 1 ] );
                }

                BaseZone baseZone;
                baseZone.zone_name = donorname;

                int gDonorZoneId = Global::zone_names.FindBaseZone( baseZone );

                zonebc_1to1->zoneid = gZoneId;
                zonebc_1to1->transform = transform;
                zonebc_1to1->donor_zoneid = gDonorZoneId;

                for ( int i = 0; i < range.size(); ++ i )
                {
                    zonebc_1to1->pnts.push_back( range[ i ] );
                    zonebc_1to1->donor_pnts.push_back( donor_range[ i ] );
                }

                if ( index_dim == 1 )
                {
                    int p1 = range[ 0 ];
                    double x1 = grid->x[ p1 - 1 ];
                    int id = Global::boundary_points.AddPoint( x1 );
                    Face face;
                    std::vector<int> nodes;
                    nodes.push_back( id );
                    face.AddFace( nodes );
                    Global::boundary_faces.AddFace( face );
                    Global::inter_faces.AddFace( face );

                    bc->bctypes.push_back( BCTypeUserDefined );
                    bc->faceids.push_back( p1 - 1 );

                    interfacezone->face_ids.push_back( p1 );

                    interfacezone->left_zones.push_back( gZoneId );
                    interfacezone->left_cells.push_back( p1 - 1 );

                    int q = donor_range[ 0 ];
                    interfacezone->right_zones.push_back( gDonorZoneId );
                    interfacezone->right_cells.push_back( q - 1 );
                }
            }
        }
    }
    cg_close( fileId );
}