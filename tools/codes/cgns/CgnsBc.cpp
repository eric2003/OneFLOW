#include <iostream>
#include <string>
#include <vector>
#include "cgnslib.h"
#include "CgnsHeader.h"
#include "CgnsBc.h"
#include "CgnsBase.h"

Boco::Boco()
{
    this->bocoId = -1;
}

Boco::~Boco()
{
    ;
}

void Boco::ReadGridLocation( int fileId, int baseId, int zoneId )
{
    GridLocation_t bcGridLocation;
    cg_boco_gridlocation_read( fileId, baseId, zoneId, this->bocoId, &bcGridLocation );

    this->gridLocation = bcGridLocation;
    this->modifiedLocation = bcGridLocation;

    //if ( cgnsZone->ExistSection( this->name ) )
    //{
    //    this->modifiedLocation = CGNS_ENUMV( FaceCenter );
    //}

    std::cout << "   CGNS Grid Location Name        = " << GridLocationName[ bcGridLocation ] << "\n";
}

void Boco::Read( int fileId, int baseId, int zoneId, int bocoId )
{
    this->bocoId = bocoId;
    // Read the info for this boundary condition.

    this->gridConnType = GridConnectivityTypeNull;

    cg_goto( fileId, baseId, "Zone_t", 1, "ZoneBC_t", 1, "BC_t", this->bocoId, "end" );

    this->ReadGridLocation( fileId, baseId, zoneId );

    cg_boco_id( fileId, baseId, zoneId, this->bocoId, & this->bc_double_id );

    char33 bcRegionName;
    cg_boco_info( fileId, baseId, zoneId, this->bocoId,
        bcRegionName, & this->bcType, & this->pointSetType, & this->nElements,
        this->normalIndex,  & this->normalListSize, & this->normalDataType, & this->nDataSets );

    this->name = bcRegionName;

    std::cout << "   CGNS Boundary Name             = " << bcRegionName << "\n";
    std::cout << "   CGNS Boundary Condition Name   = " << BCTypeName[ this->bcType ] << "\n";

    ZoneType_t zoneType = ::GetZone( fileId, baseId, zoneId )->zoneType;
    if ( zoneType == CGNS_ENUMV( Unstructured ) )
    {
        this->conn.resize( nElements );
    }
    else
    {
        if ( this->nElements == 2 )
        {
            this->nElements = 6;
        }
        this->conn.resize( this->nElements );
    }

    std::cout << "   CGNS PointSet Type Name        = " << PointSetTypeName[ this->pointSetType ] << "\n";

    int cgnsNormalList;

    // Read the element IDs.
    cg_boco_read( fileId, baseId, zoneId, this->bocoId, this->conn.data(), & cgnsNormalList);

    if ( zoneType == CGNS_ENUMV( Unstructured ) )
    {
        if ( this->modifiedLocation == CGNS_ENUMV( Vertex ) )
        {
            std::cout << "   CGNS Boundary Point's Number   = ";
        }
        else
        {
            std::cout << "   CGNS Boundary Element's Number = ";
        }

        if ( this->pointSetType == CGNS_ENUMV( ElementRange ) ||
            this->pointSetType == CGNS_ENUMV( PointRange   ) )
        {
            std::cout << this->nElements << "( " << this->conn[ 1 ] - this->conn[ 0 ] + 1 << " )" << "\n";
        }
        else
        {
            std::cout << this->nElements << "\n";
        }

        if ( this->nElements == 2 )
        {
            std::cout << "   " << this->conn[ 0 ] << " " << this->conn[ 1 ] << "\n";
        }
        else
        {
            std::cout << "   ";
            cgsize_t num = 5;
            for ( int i = 0; i < std::min(num, this->nElements ); ++ i )
            {
                std::cout << this->conn[ i ] << " ";
            }
            std::cout << "\n";
        }

    }
    int kkk = 1;
}

Bc::Bc()
{
    ;
}

Bc::~Bc()
{
    for ( int iBocos = 0; iBocos < this->bocos.size(); ++ iBocos )
    {
        delete this->bocos[ iBocos ];
    }
}

void Bc::AllocateBocos()
{
    for ( int iBoco = 0; iBoco < this->nBocos; ++ iBoco )
    {
        Boco * boco = new Boco();
        this->bocos.push_back( boco );
    }
}

void Bc::ReadBoundaries( int fileId, int baseId, int zoneId )
{
    // Determine the number of boundary conditions for this zone.
    cg_nbocos( fileId, baseId, zoneId, & this->nBocos );

    std::cout << "   nBocos        = " << this->nBocos << std::endl;

    this->AllocateBocos();

    for ( int iBoco = 0; iBoco < this->nBocos; ++ iBoco )
    {
        std::cout << "\n";
        std::cout << "-->iBoco  = " << iBoco << " nBocos = " << nBocos << "\n";
        int bocoId = iBoco + 1;
        Boco * boco = this->bocos[ iBoco ];
        boco->Read( fileId, baseId, zoneId, bocoId );
    }
}
