#include "CgnsGrid.h"
#include "Parallel.h"
#include "cgnslib.h"
#include "global.h"
#include "ZoneState.h"
#include <iostream>
#include <vector>

Region::Region()
{
}

Region::~Region()
{
}

Region & Region::operator = ( const Region & rhs )
{
    if ( this == & rhs ) return * this;

    this->start = rhs.start;
    this->end = rhs.end;

    return * this;
}

void Region::SetRegion( std::vector<int> & pnts )
{
    int index_dim = pnts.size() / 2;
    this->start.resize( index_dim );
    this->end.resize( index_dim );
    for ( int m = 0; m < index_dim; ++ m )
    {
        this->start[ m ] = pnts[ m ];
        this->end[ m ] = pnts[ index_dim + m ];
    }
}

void Region::Print()
{
    int nSize = this->start.size();
    std::cout << "start:(";
    for ( int m = 0; m < nSize; ++ m )
    {
        std::cout << this->start[ m ];
        if ( m != nSize - 1 )
        {
            std::cout << ",";
        }
    }
    std::cout << ")\n";
    std::cout << "end  :(";
    for ( int m = 0; m < nSize; ++ m )
    {
        std::cout << this->end[ m ];
        if ( m != nSize - 1 )
        {
            std::cout << ",";
        }
    }
    std::cout << ")\n";
}

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

ZoneBc::ZoneBc()
{
    ;
}

ZoneBc::~ZoneBc()
{
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

int Trans::M[ 3 ][ 3 ];
std::vector<int> Trans::transform;

int Trans::sgn( int x )
{
    if ( x >= 0 )
    {
        return 1;
    }
    else
    {
        return -1;
    }
}

int Trans::del( int x, int y )
{
    if ( std::abs( x ) == std::abs( y ) )
    {
        return 1;
    }
    return 0;
}

void Trans::ZeroMatrix()
{
    int dim = 3;
    for ( int j = 0; j < dim; ++ j )
    {
        for ( int i = 0; i < dim; ++ i )
        {
            Trans::M[ i ][ j ] = 0;
        }
    }
}

void Trans::CalcTransformMatrix()
{
    int dim = Trans::transform.size();
    if ( dim == 1 )
    {
        int a = Trans::transform[ 0 ];
        int sgna = Trans::sgn( a );
        int a1 = Trans::del( a, 1 );
        Trans::M[ 0 ][ 0 ] = sgna * a1;
    }
    else if ( dim == 2 )
    {
        int a = Trans::transform[ 0 ];
        int b = Trans::transform[ 1 ];
        int sgna = Trans::sgn( a );
        int sgnb = Trans::sgn( b );
        int a1 = Trans::del( a, 1 );
        int a2 = Trans::del( a, 2 );
        int b1 = Trans::del( b, 1 );
        int b2 = Trans::del( b, 2 );
        Trans::M[ 0 ][ 0 ] = sgna * a1;
        Trans::M[ 1 ][ 0 ] = sgna * a2;
        Trans::M[ 0 ][ 1 ] = sgnb * b1;
        Trans::M[ 1 ][ 1 ] = sgnb * b2;
    }
    else if ( dim == 3 )
    {
        int a = Trans::transform[ 0 ];
        int b = Trans::transform[ 1 ];
        int c = Trans::transform[ 2 ];
        int sgna = Trans::sgn( a );
        int sgnb = Trans::sgn( b );
        int sgnc = Trans::sgn( c );
        int a1 = Trans::del( a, 1 );
        int a2 = Trans::del( a, 2 );
        int a3 = Trans::del( a, 3 );
        int b1 = Trans::del( b, 1 );
        int b2 = Trans::del( b, 2 );
        int b3 = Trans::del( b, 3 );
        int c1 = Trans::del( c, 1 );
        int c2 = Trans::del( c, 2 );
        int c3 = Trans::del( c, 3 );
        Trans::M[ 0 ][ 0 ] = sgna * a1;
        Trans::M[ 1 ][ 0 ] = sgna * a2;
        Trans::M[ 2 ][ 0 ] = sgna * a3;
        Trans::M[ 0 ][ 1 ] = sgnb * b1;
        Trans::M[ 1 ][ 1 ] = sgnb * b2;
        Trans::M[ 2 ][ 1 ] = sgnb * b3;
        Trans::M[ 0 ][ 2 ] = sgnc * c1;
        Trans::M[ 1 ][ 2 ] = sgnc * c2;
        Trans::M[ 2 ][ 2 ] = sgnc * c3;
    }
}

Transform::Transform()
{
    //int dim = Dim::dim;
    int dim = 1;
    this->diff.resize( dim );
    this->mul.resize( dim );
}

Transform::~Transform()
{
    ;
}

void Transform::Init()
{
    Trans::ZeroMatrix();
    Trans::transform = this->transform;
    Trans::CalcTransformMatrix();

    int dim = 3;
    for ( int j = 0; j < dim; ++ j )
    {
        for ( int i = 0; i < dim; ++ i )
        {
            this->Mt[ i ][ j ] = Trans::M[ i ][ j ];
        }
    }
}

void Transform::MapIndex( std::vector<int> & index1, std::vector<int> & index2 )
{
    int dim = index1.size();
    for ( int m = 0; m < dim; ++ m )
    {
        this->diff[ m ] = index1[ m ] - this->begin1[ m ];
    }

    this->Multiply( diff, this->mul );

    for ( int m = 0; m < dim; ++ m )
    {
        index2[ m ] = this->mul[ m ] + this->begin2[ m ];
    }

}

void Transform::Multiply( std::vector<int> & a, std::vector<int> & b )
{
    int dim = a.size();
    for ( int i = 0; i < dim; ++ i )
    {
        b[ i ] = 0;
        for ( int j = 0; j < dim; ++ j )
        {
            b[ i ] += this->Mt[ i ][ j ] * a[ j ];
        }
    }
}

void ReadCgnsGridBaseZone( const std::string & filename )
{
    int fileId = -1;
    if ( Parallel::pid == Parallel::serverid )
    {
        cg_open( filename.c_str(), CG_MODE_READ, &fileId );
        std::cout << "Parallel::pid = " << Parallel::pid << " ";
        std::cout << "fileId = " << fileId << "\n";
    }
    
    int nbases = -1;
    if ( Parallel::pid == Parallel::serverid )
    {
        cg_nbases( fileId, &nbases );
    }
    MPI_Bcast( &nbases, 1, MPI_INT, Parallel::serverid, MPI_COMM_WORLD );
    std::cout << "Parallel::pid = " << Parallel::pid << " ";
    std::cout << "nbases = " << nbases << "\n";
    
    for ( int iBase = 0; iBase < nbases; ++ iBase )
    {
        char basename[ 33 ];
        int baseId = iBase + 1;
        int icelldim = -1;
        int iphysdim = -1;
        if ( Parallel::pid == Parallel::serverid )
        {
            cg_base_read( fileId, baseId, basename, &icelldim, &iphysdim );
        }

        MPI_Bcast( &icelldim, 1, MPI_INT, Parallel::serverid, MPI_COMM_WORLD );
        MPI_Bcast( &iphysdim, 1, MPI_INT, Parallel::serverid, MPI_COMM_WORLD );

        std::cout << "Parallel::pid = " << Parallel::pid << " ";
        std::cout << "icelldim = " << icelldim << " iphysdim = " << iphysdim << "\n";

        int nzones = -1;
        if ( Parallel::pid == Parallel::serverid )
        {
            cg_nzones( fileId, baseId, &nzones );
        }

        MPI_Bcast( &nzones, 1, MPI_INT, Parallel::serverid, MPI_COMM_WORLD );

        ZoneState::nZones = nzones;
        ZoneState::pids.resize( nzones );

        if ( Parallel::nProc == 1 )
        {
            for ( int iZone = 0; iZone < ZoneState::nZones; ++ iZone )
            {
                ZoneState::pids[ iZone ] = 0;
            }
        }
        else
        {
            for ( int iZone = 0; iZone < ZoneState::nZones; ++ iZone )
            {
                ZoneState::pids[ iZone ] = iZone;
            }
        }

        std::cout << "Parallel::pid = " << Parallel::pid << " ";
        std::cout << "nzones = " << nzones << "\n";

        for ( int iZone = 0; iZone < nzones; ++ iZone )
        {
            int zoneId = iZone + 1;
            std::cout << "Parallel::pid = " << Parallel::pid << " ";
            std::cout << "iZone = " << iZone << "\n";

            int index_dim = -1;
            if ( Parallel::pid == Parallel::serverid )
            {
                cg_index_dim( fileId, baseId, zoneId, &index_dim );
            }
            MPI_Bcast( &index_dim, 1, MPI_INT, Parallel::serverid, MPI_COMM_WORLD );

            std::cout << "Parallel::pid = " << Parallel::pid << " ";
            std::cout << "index_dim = " << index_dim << "\n";

            std::vector<cgsize_t> isize( index_dim * 3 );

            char zonename[ 33 ];
            if ( Parallel::pid == Parallel::serverid )
            {
                cg_zone_read( fileId, baseId, zoneId, zonename, isize.data() );
            }

            MPI_Bcast( zonename, 33, MPI_CHAR, Parallel::serverid, MPI_COMM_WORLD );
            MPI_Bcast( isize.data(), index_dim * 3, MPI_LONG_LONG, Parallel::serverid, MPI_COMM_WORLD );

            std::cout << "Parallel::pid = " << Parallel::pid << " ";
            std::cout << "zonename = " << zonename << "\n";

            for ( int i = 0; i < isize.size(); ++ i )
            {
                std::cout << "Parallel::pid = " << Parallel::pid << " ";
                std::cout << "i = " << i << " isize = " << isize[ i ] << "\n";
            }

            BaseZone baseZone;
            baseZone.zone_name = zonename;

            Global::zone_names.AddBaseZone( baseZone );
        }
    }

    if ( Parallel::pid == Parallel::serverid )
    {
        cg_close( fileId );
    }
}

void ReadCgnsGrid( const std::string & filename )
{
    int fileId = -1;
    if ( Parallel::pid == Parallel::serverid )
    {
        cg_open( filename.c_str(), CG_MODE_READ, &fileId );
        std::cout << "fileId = " << fileId << "\n";
    }
    
    int nbases = -1;
    if ( Parallel::pid == Parallel::serverid )
    {
        cg_nbases( fileId, &nbases );
    }
    MPI_Bcast( &nbases, 1, MPI_INT, Parallel::serverid, MPI_COMM_WORLD );
    std::cout << "Parallel::pid = " << Parallel::pid << " ";
    std::cout << "nbases = " << nbases << "\n";

    for ( int iBase = 0; iBase < nbases; ++ iBase )
    {
        char basename[ 33 ];
        int baseId = iBase + 1;
        int icelldim = -1;
        int iphysdim = -1;
        if ( Parallel::pid == Parallel::serverid )
        {
            cg_base_read( fileId, baseId, basename, &icelldim, &iphysdim );
        }
        MPI_Bcast( &icelldim, 1, MPI_INT, Parallel::serverid, MPI_COMM_WORLD );
        MPI_Bcast( &iphysdim, 1, MPI_INT, Parallel::serverid, MPI_COMM_WORLD );
        std::cout << "Parallel::pid = " << Parallel::pid << " ";
        std::cout << "icelldim = " << icelldim << " iphysdim = " << iphysdim << "\n";

        Global::cell_dim = icelldim;
        Global::phys_dim = iphysdim;

        int nzones = -1;
        if ( Parallel::pid == Parallel::serverid )
        {
            cg_nzones( fileId, baseId, &nzones );
        }
        MPI_Bcast( &nzones, 1, MPI_INT, Parallel::serverid, MPI_COMM_WORLD );
        std::cout << "Parallel::pid = " << Parallel::pid << " ";
        std::cout << "nzones = " << nzones << "\n";

        for ( int iZone = 0; iZone < nzones; ++ iZone )
        {
            int zoneId = iZone + 1;
            int index_dim = -1;
            if ( Parallel::pid == Parallel::serverid )
            {
                cg_index_dim( fileId, baseId, zoneId, &index_dim );
            }
            HXSendRecvData<int>( &index_dim, 1, Parallel::serverid, ZoneState::GetProcID( iZone ) );

            if ( ZoneState::IsValid( iZone ) )
            {
                std::cout << "Parallel::pid = " << Parallel::pid << " ";
                std::cout << "index_dim = " << index_dim << "\n";
            }

            //std::vector<cgsize_t> isize( index_dim * 3 );
            std::vector<cgsize_t> isize;
            char zonename[ 33 ];
            if ( Parallel::pid == Parallel::serverid )
            {
                isize.resize( index_dim * 3 );
                cg_zone_read( fileId, baseId, zoneId, zonename, isize.data() );
            }

            if ( ZoneState::IsValid( iZone ) )
            {
                isize.resize( index_dim * 3 );
            }

            //MPI_Bcast( zonename, 33, MPI_CHAR, Parallel::serverid, MPI_COMM_WORLD );
            //MPI_Bcast( isize.data(), index_dim * 3, MPI_LONG_LONG, Parallel::serverid, MPI_COMM_WORLD );

            HXSendRecvData( zonename, 33, Parallel::serverid, ZoneState::GetProcID( iZone ) );
            HXSendRecvData( isize.data(), index_dim * 3, Parallel::serverid, ZoneState::GetProcID( iZone ) );

            if ( ZoneState::IsValid( iZone ) )
            {
                std::cout << "Parallel::pid = " << Parallel::pid << " ";
                std::cout << "zonename = " << zonename << "\n";

                for ( int i = 0; i < isize.size(); ++ i )
                {
                    std::cout << "Parallel::pid = " << Parallel::pid << " ";
                    std::cout << "i = " << i << " isize = " << isize[ i ] << "\n";
                }
            }

            std::vector<cgsize_t> irmin;
            std::vector<cgsize_t> irmax;
            int nNodes = 1;
            if ( ZoneState::IsValid( iZone ) || Parallel::IsServer() )
            {
                irmin.resize( index_dim );
                irmax.resize( index_dim );
                for ( int m = 0; m < index_dim; ++ m )
                {
                    /* lower range index */
                    irmin[ m ] = 1;
                    /* upper range index of vertices */
                    irmax[ m ] = isize[ m ];
                    nNodes *= irmax[ m ];
                }
                std::cout << "Parallel::pid = " << Parallel::pid << " ";
                std::cout << "nNodes = " << nNodes << "\n";
            }

            ZoneType_t zoneType;
            if ( Parallel::pid == Parallel::serverid )
            {
                cg_zone_type( fileId, baseId, zoneId, &zoneType );
            }
            HXSendRecvData( &zoneType, 1, Parallel::serverid, ZoneState::GetProcID( iZone ) );

            if ( ZoneState::IsValid( iZone ) || Parallel::IsServer() )
            {
                std::cout << "Parallel::pid = " << Parallel::pid << " ";
                std::cout << "zoneType = " << zoneType << " ZoneTypeName = " << ZoneTypeName[ zoneType ] << "\n";
            }

            Zone * zone = nullptr;
            if ( ZoneState::IsValid( iZone ) )
            {
                zone = new Zone();
                Global::zones.push_back( zone );
                for ( int m = 0; m < index_dim; ++ m )
                {
                    zone->nijk.push_back( isize[ m ] );
                }
            }

            Grid * grid = nullptr;
            if ( ZoneState::IsValid( iZone ) )
            {
                grid = new Grid();
                Global::grids.push_back( grid );
            }

            BaseZone baseZone;

            int gZoneId = -1;
            if ( Parallel::IsServer() )
            {
                baseZone.zone_name = zonename;
                gZoneId = Global::zone_names.FindBaseZone( baseZone );
            }

            HXSendRecvData( &gZoneId, 1, Parallel::serverid, ZoneState::GetProcID( iZone ) );

            if ( ZoneState::IsValid( iZone ) || Parallel::IsServer() )
            {
                std::cout << "Parallel::pid = " << Parallel::pid << " ";
                std::cout << "gZoneId = " << gZoneId << "\n";
            }

            int ncoords = -1;
            if ( Parallel::IsServer() )
            {
                cg_ncoords( fileId, baseId, zoneId, &ncoords );
            }
            HXSendRecvData( &ncoords, 1, Parallel::serverid, ZoneState::GetProcID( iZone ) );
            if ( ZoneState::IsValid( iZone ) || Parallel::IsServer() )
            {
                std::cout << "Parallel::pid = " << Parallel::pid << " ";
                std::cout << "ncoords = " << ncoords << "\n";
            }

            for ( int icoord = 0; icoord < ncoords; ++ icoord )
            {
                int coorId = icoord + 1;
                DataType_t dataType;
                char coordname[ 33 ];
                if ( Parallel::IsServer() )
                {
                    cg_coord_info( fileId, baseId, zoneId, coorId, &dataType, coordname );
                }
                HXSendRecvData( coordname, 33, Parallel::serverid, ZoneState::GetProcID( iZone ) );
                HXSendRecvData( &dataType, 1, Parallel::serverid, ZoneState::GetProcID( iZone ) );
                if ( ZoneState::IsValid( iZone ) || Parallel::IsServer() )
                {
                    std::cout << "Parallel::pid = " << Parallel::pid << " ";
                    std::cout << "coordname = " << coordname << "\n";
                    std::cout << "Parallel::pid = " << Parallel::pid << " ";
                    std::cout << "dataType = " << dataType << " DataTypeName = " << DataTypeName[ dataType ] << "\n";
                }

                std::vector<char> coord;

                if ( ZoneState::IsValid( iZone ) || Parallel::IsServer() )
                {
                    coord.resize( nNodes * sizeof( double ) );
                }

                if ( Parallel::IsServer() )
                {
                    cg_coord_read( fileId, baseId, zoneId, coordname, dataType, irmin.data(), irmax.data(), coord.data() );
                }

                HXSendRecvData( coord.data(), coord.size(), Parallel::serverid, ZoneState::GetProcID( iZone ) );

                Coor * coor = nullptr;
                if ( ZoneState::IsValid( iZone ) )
                {
                    coor = new Coor();
                    zone->coors.push_back( coor );
                    coor->coorname = coordname;
                    coor->nNodes = nNodes;
                    coor->coord = coord;

                    std::cout << "Parallel::pid = " << Parallel::pid << " ";
                    coor->DumpCoor();

                    if ( icoord == 0 )
                    {
                        grid->x.resize( nNodes );
                        coor->DumpCoorX( grid->x );
                    }
                }
            }

            int nbocos = -1;
            if ( Parallel::IsServer() )
            {
                cg_nbocos( fileId, baseId, zoneId, &nbocos );
            }
            HXSendRecvData( &nbocos, 1, Parallel::serverid, ZoneState::GetProcID( iZone ) );

            if ( ZoneState::IsValid( iZone ) || Parallel::IsServer() )
            {
                std::cout << "Parallel::pid = " << Parallel::pid << " ";
                std::cout << "nbocos = " << nbocos << "\n";
            }

            for ( int iboco = 0; iboco < nbocos; ++ iboco )
            {
                int bccoId = iboco + 1;
                GridLocation_t location;
                if ( Parallel::IsServer() )
                {
                    cg_boco_gridlocation_read( fileId, baseId, zoneId, bccoId, &location );
                }
                HXSendRecvData( &location, 1, Parallel::serverid, ZoneState::GetProcID( iZone ) );

                if ( ZoneState::IsValid( iZone ) || Parallel::IsServer() )
                {
                    std::cout << "Parallel::pid = " << Parallel::pid << " ";
                    std::cout << "iboco = " << iboco <<  " location = " << location << " GridLocationName = " << GridLocationName[location] << "\n";
                }

                char boconame[ 33 ];
                BCType_t bocotype;
                PointSetType_t ptset_type;
                cgsize_t npnts = 0;
                //std::vector<int> normalIndex( index_dim, -1 );
                std::vector<int> normalIndex;
                cgsize_t normalListSize = 0;
                DataType_t normalDataType;
                int ndataset = -1;
                if ( Parallel::IsServer() )
                {
                    normalIndex.resize( index_dim );
                    cg_boco_info( fileId, baseId, zoneId, bccoId, boconame, &bocotype, &ptset_type,
                        &npnts, normalIndex.data(), &normalListSize, &normalDataType, &ndataset );
                }

                if ( ZoneState::IsValid( iZone ) )
                {
                    normalIndex.resize( index_dim );
                }
                
                HXSendRecvData( boconame, 33, Parallel::serverid, ZoneState::GetProcID( iZone ) );
                HXSendRecvData( &bocotype, 1, Parallel::serverid, ZoneState::GetProcID( iZone ) );
                HXSendRecvData( &ptset_type, 1, Parallel::serverid, ZoneState::GetProcID( iZone ) );
                HXSendRecvData( &npnts, 1, Parallel::serverid, ZoneState::GetProcID( iZone ) );
                HXSendRecvData( &normalDataType, 1, Parallel::serverid, ZoneState::GetProcID( iZone ) );
                HXSendRecvData( &normalListSize, 1, Parallel::serverid, ZoneState::GetProcID( iZone ) );
                HXSendRecvData( normalIndex.data(), index_dim, Parallel::serverid, ZoneState::GetProcID( iZone ) );
                HXSendRecvData( &ndataset, 1, Parallel::serverid, ZoneState::GetProcID( iZone ) );

                if ( ZoneState::IsValid( iZone ) || Parallel::IsServer() )
                {
                    std::cout << "Parallel::pid = " << Parallel::pid << " ";
                    std::cout << "boconame = " << boconame << " bocotype = " << bocotype << " BCTypeName = " << BCTypeName[ bocotype ] << "\n";
                    std::cout << "Parallel::pid = " << Parallel::pid << " ";
                    std::cout << "ptset_type = " << ptset_type <<  " PointSetTypeName = " << PointSetTypeName[ptset_type] << "\n";
                    std::cout << "Parallel::pid = " << Parallel::pid << " ";
                    std::cout << "npnts = " << npnts << "\n";
                    std::cout << "Parallel::pid = " << Parallel::pid << " ";
                    std::cout << "normalListSize = " << normalListSize << "\n";
                    std::cout << "Parallel::pid = " << Parallel::pid << " ";
                    std::cout << "normalDataType = " << normalDataType << " DataTypeName = " << DataTypeName[ normalDataType ] << "\n";
                    std::cout << "Parallel::pid = " << Parallel::pid << " ";
                    std::cout << "normalIndex = ";
                    for ( int i = 0; i < index_dim; ++ i )
                    {
                        std::cout << normalIndex[ i ] << " ";
                    }
                    std::cout << "\n";
                    std::cout << "Parallel::pid = " << Parallel::pid << " ";
                    std::cout << "ndataset = " << ndataset << "\n";
                }

                std::vector<char> normalList;
                std::vector<cgsize_t> pnts;

                if ( ZoneState::IsValid( iZone ) || Parallel::IsServer() )
                {
                    if ( normalDataType == DataTypeNull )
                    {
                        normalList.resize( sizeof( int ) );
                    }
                    else
                    {
                        int nSize = nNodes * index_dim * sizeof( double );
                        normalList.resize( nSize );
                    }
                    pnts.resize( npnts * index_dim );
                }

                if ( Parallel::IsServer() )
                {
                    cg_boco_read( fileId, baseId, zoneId, bccoId, pnts.data(), normalList.data() );
                }

                HXSendRecvData( pnts.data(), npnts * index_dim, Parallel::serverid, ZoneState::GetProcID( iZone ) );

                if ( ZoneState::IsValid( iZone ) || Parallel::IsServer() )
                {
                    std::cout << "Parallel::pid = " << Parallel::pid << " ";
                    std::cout << "pnts = ";
                    for ( int i = 0; i < pnts.size(); ++ i )
                    {
                        std::cout << pnts[ i ] << " ";
                    }
                    std::cout << "\n";
                }
                HXSendRecvData( normalList.data(), normalList.size(), Parallel::serverid, ZoneState::GetProcID(iZone));

                ZoneBc * zonebc = nullptr;
                if ( ZoneState::IsValid( iZone ) )
                {
                    zonebc = new ZoneBc();
                    zone->bccos.push_back( zonebc );
                    zonebc->bcType = bocotype;
                    for ( int i = 0; i < pnts.size(); ++ i )
                    {
                        zonebc->pnts.push_back( pnts[ i ] );
                    }
                }
            }
            int n1to1 = -1;
            if ( Parallel::IsServer() )
            {
                cg_n1to1( fileId, baseId, zoneId, &n1to1 );
            }

            HXSendRecvData( &n1to1, 1, Parallel::serverid, ZoneState::GetProcID( iZone ) );

            if ( ZoneState::IsValid( iZone ) || Parallel::IsServer() )
            {
                std::cout << "Parallel::pid = " << Parallel::pid << " ";
                std::cout << "n1to1 = " << n1to1 << "\n";
            }

            for ( int i1to1 = 0; i1to1 < n1to1; ++ i1to1 )
            {
                int i1to1Id = i1to1 + 1;

                //ZoneBc1To1 * zonebc_1to1 = new ZoneBc1To1();
                //zone->bc1to1s.push_back( zonebc_1to1 );

                char connectname[ 33 ];
                char donorname[ 33 ];
                cgsize_t npnts = 2;
                std::vector<cgsize_t> range;
                std::vector<cgsize_t> donor_range;
                std::vector<int> transform;
                if ( ZoneState::IsValid( iZone ) || Parallel::IsServer() )
                {
                    range.resize( npnts * index_dim );
                    donor_range.resize( npnts * index_dim );
                    transform.resize( index_dim );
                }

                if ( Parallel::IsServer() )
                {
                    cg_1to1_read( fileId, baseId, zoneId, i1to1Id, connectname, donorname, range.data(), donor_range.data(), transform.data() );
                }

                HXSendRecvData( connectname, 33, Parallel::serverid, ZoneState::GetProcID( iZone ) );
                HXSendRecvData( donorname, 33, Parallel::serverid, ZoneState::GetProcID( iZone ) );
                HXSendRecvData( range.data(), range.size(), Parallel::serverid, ZoneState::GetProcID( iZone ) );
                HXSendRecvData( donor_range.data(), donor_range.size(), Parallel::serverid, ZoneState::GetProcID( iZone ) );
                HXSendRecvData( transform.data(), transform.size(), Parallel::serverid, ZoneState::GetProcID( iZone ) );

                if ( ZoneState::IsValid( iZone ) || Parallel::IsServer() )
                {
                    std::cout << "Parallel::pid = " << Parallel::pid << " ";
                    std::cout << "connectname = " << connectname << "\n";
                    std::cout << "Parallel::pid = " << Parallel::pid << " ";
                    std::cout << "donorname = " << donorname << "\n";

                    std::cout << "Parallel::pid = " << Parallel::pid << " ";
                    std::cout << "range = ";
                    for ( int i = 0; i < range.size(); ++ i )
                    {
                        std::cout << range[ i ] << " ";
                    }
                    std::cout << "\n";
                    std::cout << "Parallel::pid = " << Parallel::pid << " ";
                    std::cout << "donor_range = ";
                    for ( int i = 0; i < donor_range.size(); ++ i )
                    {
                        std::cout << donor_range[ i ] << " ";
                    }
                    std::cout << "\n";

                    std::cout << "Parallel::pid = " << Parallel::pid << " ";
                    std::cout << "transform = ";
                    for ( int i = 0; i < transform.size(); ++ i )
                    {
                        std::cout << transform[ i ] << " ";
                    }
                    std::cout << "\n";
                }

                int gDonorZoneId = -1;
                if ( Parallel::IsServer() )
                {
                    BaseZone baseZoneDonor;
                    baseZoneDonor.zone_name = donorname;

                    gDonorZoneId = Global::zone_names.FindBaseZone( baseZoneDonor );
                }

                HXSendRecvData( &gDonorZoneId, 1, Parallel::serverid, ZoneState::GetProcID( iZone ) );

                ZoneBc1To1 * zonebc_1to1 = nullptr;
                if ( ZoneState::IsValid( iZone ) )
                {
                    zonebc_1to1 = new ZoneBc1To1();
                    zone->bc1to1s.push_back( zonebc_1to1 );
                    zonebc_1to1->zoneid = gZoneId;
                    zonebc_1to1->transform = transform;
                    zonebc_1to1->donor_zoneid = gDonorZoneId;

                    for ( int i = 0; i < range.size(); ++ i )
                    {
                        zonebc_1to1->pnts.push_back( range[ i ] );
                        zonebc_1to1->donor_pnts.push_back( donor_range[ i ] );
                    }
                }
            }
        }
    }

    {
        std::cout << "Parallel::pid = " << Parallel::pid << " ";
        std::cout << "Global::zones.size() = " << Global::zones.size() << "\n";
    }

    cg_close( fileId );
}