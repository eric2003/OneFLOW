#include "global.h"
#include "CgnsGrid.h"

std::vector<Grid *> Global::grids;
std::vector<Field *> Global::fields;

int Global::nt = -1;
int Global::cell_dim = -1;
int Global::phys_dim = -1;

std::vector<Zone *> Global::zones;
BaseZoneList Global::zone_names;
std::vector<Interface *> Global::interfaces;

std::map<Face, int> Global::faceMap;
std::map<FacePair, int> Global::facePairMap;
std::vector<FacePair> Global::facePairList;

std::vector<std::set<int>> Global::donor_zone_sets;
std::vector<std::vector<int>> Global::donor_zones;

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

void Interface::CalcInterface( Transform * transform, std::vector<int> & start, std::vector<int> & end, int donor_zoneid )
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

        Global::InsertFaceMap( face );

        int i_donor = index2[ 0 ];

        Face face_donor;
        face_donor.zone = donor_zoneid;
        face_donor.i = i_donor;

        Global::InsertFaceMap( face_donor );

        FacePair facePair;
        facePair.AddPair( face, face_donor );

        Global::facePairList.push_back( facePair );
        int global_faceid = Global::InsertFacePairMap( facePair );
        this->global_faceids.push_back( global_faceid );

        this->global_local_face_map.insert( std::make_pair( global_faceid, faceid ) );

        int i_ghost_cell = i + 1;
        int i_local_donor_cell = i - 1;

        if ( i == 1 )
        {
            i_ghost_cell = i - 1;
            i_local_donor_cell = i + 1;
        }
        ijk_ghosts.push_back( i_ghost_cell );
        ijk_donors.push_back( i_local_donor_cell );

        icount ++;
    }
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

void Global::SendGeom( int zone, int donorzone, std::vector<int> & donorfaces )
{
    Interface * interface = Global::interfaces[ donorzone ];

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
        int ijkpos = index_dim * local_faceid;
        int i_donor_cell = interface->ijk_donors[ ijkpos + 0 ];
        sub_donorijk.push_back( i_donor_cell );
        int kkk = 1;
    }

    interface->donorijk_for_send.push_back( sub_donorijk );

    std::vector<double> sub_donordata( sub_donorijk.size() );
    interface->donordata_for_send.push_back( sub_donordata );
}

void Global::ReSendGeomTest( int zone, int zone_to_send, std::vector<int> & donorfaces )
{
    Interface * interface = Global::interfaces[ zone_to_send ];
    interface->tmp_recv_donor_zones.push_back( zone );
    interface->tmp_recv_donorfaces.push_back( donorfaces );
}
