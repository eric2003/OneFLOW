#include "Global.h"
#include "Grid.h"
#include "ZoneState.h"
#include "Parallel.h"
#include <iostream>

std::vector<Grid *> Global::grids;
std::vector<Field *> Global::fields;

int Global::nt = -1;
int Global::iter = -1;
int Global::cell_dim = -1;
int Global::phys_dim = -1;
int Global::nghost = -1;

std::vector<Zone *> Global::zones;
std::vector<Interface *> Global::interfaces;

std::map<Face, int> Global::faceMap;
std::map<FacePair, int> Global::facePairMap;
std::vector<FacePair> Global::facePairList;
std::vector<FacePair> Global::mpi_facePairList;

std::vector<std::set<int>> Global::donor_zone_sets;
std::vector<std::vector<int>> Global::donor_zones;

InterfaceTopo Global::interfaceTopo;

bool Face::operator < ( const Face & rhs ) const
{
    if ( this->zone != rhs.zone )
    {
        return this->zone < rhs.zone;
    }

    if ( this->i != rhs.i )
    {
        return this->i < rhs.i;
    }

    if ( this->j != rhs.j )
    {
        return this->j < rhs.j;
    }

    return this->k < rhs.k;
}

bool Face::operator == ( const Face & rhs ) const
{
    if ( this->zone != rhs.zone )
    {
        return false;
    }

    if ( this->i != rhs.i )
    {
        return false;
    }

    if ( this->j != rhs.j )
    {
        return false;
    }

    return this->k == rhs.k;
}

void Face::Print()
{
    std::cout << "(" << this->zone << "," << this->i << ")";
}

void FacePair::AddPair( const Face & face1, const Face & face2 )
{
    if ( face1 < face2 )
    {
        this->left = face1;
        this->right = face2;
    }
    else
    {
        this->left = face2;
        this->right = face1;
    }
}

bool FacePair::operator < ( const FacePair & rhs ) const
{
    if ( this->left == rhs.left || this->left == rhs.right )
    {
        return false;
    }

    return this->left < rhs.left;
}

void FacePair::Print()
{
    this->left.Print();
    std::cout << " ";
    this->right.Print();
    std::cout << "\n";
}

void InterfaceTopo::InitNeighborInfo()
{
    this->linkmap.resize( ZoneState::nZones );

    for ( int iZone = 0; iZone < ZoneState::nZones; ++ iZone )
    {
        if ( ! ZoneState::IsValid( iZone ) ) continue;

        int local_zoneid = ZoneState::g2lzoneids[ iZone ];

        Interface * interface = Global::interfaces[ local_zoneid ];

        std::vector<int> & t = this->linkmap[ iZone ];
        t = interface->neighbor_donor_zones;
    }
}

void InterfaceTopo::SwapNeighborInfo()
{
    for ( int iZone = 0; iZone < ZoneState::nZones; ++ iZone )
    {
        int pid = ZoneState::pids[ iZone ];

        std::vector<int> & donor_zones = this->linkmap[ iZone ];
        int nNeighbor = donor_zones.size();

        HXBcastData( &nNeighbor, 1, pid );

        donor_zones.resize( nNeighbor );

        HXBcastData( donor_zones.data(), donor_zones.size(), pid );
    }

    this->SwapNeighborDonorfaces();
}

void InterfaceTopo::SwapNeighborDonorfaces()
{
    int gl = 0;
    for ( int iZone = 0; iZone < ZoneState::nZones; ++ iZone )
    {
        int send_pid = ZoneState::pids[ iZone ];

        std::vector<int> & donor_zones = this->linkmap[ iZone ];
        int ndonor_zones = donor_zones.size();

        for ( int iNei = 0; iNei < ndonor_zones; ++ iNei )
        {
            int donor_zone = donor_zones[ iNei ];
            int recv_pid = ZoneState::pids[ donor_zone ];
            int nInterFaces = 0;
            std::vector<int> donorfaces;

            if ( Parallel::pid == send_pid )
            {
                int local_zoneid = ZoneState::g2lzoneids[ iZone ];

                Interface * interface = Global::interfaces[ local_zoneid ];

                std::vector<std::vector<int>> & neighbor_donorfaces = interface->neighbor_donorfaces;

                std::vector<int> & neighbor_donorface = neighbor_donorfaces[ iNei ];

                nInterFaces = neighbor_donorface.size();

                donorfaces = neighbor_donorface;
            }

            HXSendRecvData( &nInterFaces, 1, send_pid, recv_pid );

            if ( Parallel::pid == recv_pid && send_pid != recv_pid )
            {
                donorfaces.resize( nInterFaces );
            }

            HXSendRecvData( donorfaces.data(), donorfaces.size(), send_pid, recv_pid );

            if ( Parallel::pid == recv_pid )
            {
                int local_donor_zoneid = ZoneState::g2lzoneids[ donor_zone ];
                Interface * interface_recv = Global::interfaces[ local_donor_zoneid ];
                interface_recv->SendGeom( iZone, donorfaces );
            }
        }
    }
}


void Interface::CalcInterface( Transform * transform, std::vector<int> & start, std::vector<int> & end, int donor_zoneid, int nghost )
{
    int ist = start[ 0 ];
    int ied = end[ 0 ];
    int dim = start.size();
    std::vector<int> index1( dim );
    std::vector<int> index2( dim );

    int icount = this->zoneList.size();
    for ( int i = ist; i <= ied; ++ i )
    {
        int faceid = icount;
        this->zoneList.push_back( donor_zoneid );
        this->local_faceids.push_back( faceid );
        index1[ 0 ] = i;
        transform->MapIndex( index1, index2 );
        Face face;
        face.zone = zoneid;
        face.i = i;

        int i_donor = index2[ 0 ];

        Face face_donor;
        face_donor.zone = donor_zoneid;
        face_donor.i = i_donor;

        FacePair facePair;
        facePair.AddPair( face, face_donor );

        Global::facePairList.push_back( facePair );
        int nSize = Global::facePairList.size();
        this->proc_global_faceids.push_back( nSize - 1 );

        if ( i == 1 ) {
            for ( int ig = 1; ig <= nghost; ++ ig )
            {
                ijk_ghosts.push_back( i - ig );
                ijk_donors.push_back( i + ig );
            }
        }
        else {
            for ( int ig = 1; ig <= nghost; ++ ig )
            {
                ijk_ghosts.push_back( i + ig );
                ijk_donors.push_back( i - ig );
            }
        }

        icount ++;
    }
}

void Interface::SendGeom( int zone, std::vector<int> & donorfaces )
{
    Interface * interface = this;

    std::vector<int> & send_to_zones = interface->send_to_zones;

    send_to_zones.push_back( zone );
    interface->donorfaces_for_send.push_back( donorfaces );

    int nface = donorfaces.size();
    std::vector<int> sub_donorijk;
    int index_dim = 1;
    for ( int i = 0; i < nface; ++ i )
    {
        int global_faceid = donorfaces[ i ];
        int local_faceid = interface->global_local_face_map[ global_faceid ];
        int ijkpos = index_dim * local_faceid * Global::nghost;

        for ( int ig = 0; ig < Global::nghost; ++ ig )
        {
            int i_donor_cell = interface->ijk_donors[ ijkpos + ig ];
            sub_donorijk.push_back( i_donor_cell );
        }
        int kkk = 1;
    }

    interface->donorijk_for_send.push_back( sub_donorijk );

    std::vector<double> sub_donordata( sub_donorijk.size() );
    interface->donordata_for_send.push_back( sub_donordata );
}

void Global::InsertFaceMap( const Face & face )
{
    std::map<Face, int>::iterator iter;
    iter = Global::faceMap.find( face );
    if ( iter == Global::faceMap.end() )
    {
        int faceid = Global::faceMap.size();
        Global::faceMap.insert( std::make_pair( face, faceid ) );
    }
}

int Global::InsertFacePairMap( const FacePair & facePair )
{
    std::map<FacePair, int>::iterator iter;
    iter = Global::facePairMap.find( facePair );
    if ( iter == Global::facePairMap.end() )
    {
        int facePairId = Global::facePairMap.size();
        Global::facePairMap.insert( std::make_pair( facePair, facePairId ) );
        return facePairId;
    }
    return iter->second;
}

void Global::AddFacePairList( std::vector<FacePair> & a, std::vector<FacePair> & b )
{
    for ( int i = 0; i < b.size(); ++ i )
    {
        a.push_back( b[ i ] );
    }
}
