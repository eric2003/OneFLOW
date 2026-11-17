#include "global.h"
#include "CgnsGrid.h"

std::vector<Grid *> Global::grids;
std::vector<Field *> Global::fields;

int Global::nt = -1;
int Global::cell_dim = -1;
int Global::phys_dim = -1;

std::vector<Zone *> Global::zones;
BaseZoneList Global::zone_names;
std::vector<Interface *>  Global::interfaces;

void Interface::CalcInterface( Transform * transform, std::vector<int> &start,  std::vector<int> &end, int donor_zoneid )
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
        this->faceidList.push_back( faceid );
        index1[ 0 ] = i;
        transform->MapIndex( index1, index2 );
        int i_donor = index2[ 0 ];

        int i_ghost_cell = i + 1;

        if ( i == 1 )
        {
            i_ghost_cell = i - 1;
        }
        ijk_ghosts.push_back( i_ghost_cell );

        int i_donor_cell = i_donor - 1;

        if ( i_donor == 1 )
        {
            i_donor_cell = i_donor + 1;
        }

        this->ijk_donors.push_back( i_donor_cell );

        icount ++;
    }
}