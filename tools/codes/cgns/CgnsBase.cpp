#include "CgnsBase.h"
#include "CgnsBc.h"
#include "Plot.h"
#include <fstream>
#include <iostream>
#include <iomanip>
#include <algorithm>

std::string GetCgnsFileTypeName( int file_type )
{
    std::string fileTypeName;
    if ( file_type == CG_FILE_ADF )
    {
        fileTypeName = "CG_FILE_ADF";
    }
    else if ( file_type == CG_FILE_HDF5 )
    {
        fileTypeName = "CG_FILE_HDF5";
    }
    else if ( file_type == CG_FILE_ADF2 )
    {
        fileTypeName = "CG_FILE_ADF2";
    }
    else
    {
        fileTypeName = "CG_FILE_NONE";
    }
    return fileTypeName;
}

BBase bbase;

BBase * GetBBase()
{
    return &bbase;
}

Zone * GetZone( int fileId, int BaseId, int zoneId )
{
    return bbase.GetBase( BaseId )->zones[ zoneId - 1];
}

void ReadCgnsFile( const std::string &fileName )
{
    bbase.ReadFile( fileName );
}

BBase::BBase()
{
    this->nBases = -1;
    this->fileId = -1;
    this->file_type = -1;
    this->fileVersion = -1.0;
    this->precision = -1;
}

BBase::~BBase()
{
    for ( int iBase = 0; iBase < this->bases.size(); ++ iBase )
    {
        delete this->bases[ iBase ];
    }
}

void BBase::AllocateBases( int nBases )
{
    this->nBases = nBases;
    for ( int iBase = 0; iBase < this->nBases; ++ iBase )
    {
        Base * base = new Base();
        this->bases.push_back( base );
    }
}


void BBase::ReadFile(const std::string &fileName)
{
    //namespace fs = std::filesystem;
    //std::cout << "Current path is " << fs::current_path() << '\n';

    int openStatus = cg_open( fileName.c_str(), CG_MODE_READ, & this->fileId );
    std::string stars( "**************************************************************" );
    std::cout << stars << "\n";
    std::cout << "   CGNS File Index = " << this->fileId << "\n";

    if ( openStatus != CG_OK )
    {
        cg_error_exit();
    }

    cg_version( fileId, & this->fileVersion );

    std::cout << "   CGNS File Version = " << std::setiosflags( std::ios::fixed ) << std::setprecision( 4 ) << this->fileVersion << "\n";

    cg_precision( fileId, & precision );

    std::cout << "   CGNS Precision = " << precision << "\n";

    cg_get_file_type( fileId, & this->file_type );

    std::cout << "   CGNS File Type = " << this->file_type << " FileTypeName = " << GetCgnsFileTypeName( this->file_type ) << "\n";
    std::cout << stars << "\n";

    cg_nbases( fileId, & this->nBases );
    std::cout << "   Total number of CGNS Base = " << this->nBases << "\n";

    this->AllocateBases( nBases );

    for ( int iBase = 0; iBase < nBases; ++ iBase )
    {
        int baseId = iBase + 1;
        std::vector<Base *> bases;
        Base * base = this->bases[ iBase ];
        base->ReadBase(fileId,  baseId);
    }

    cg_close( fileId );
}

Zone * BBase::GetZone()
{
    return GetBase()->zones[0];
}

Zone * BBase::GetZone( int zoneId )
{
    return GetBase()->zones[zoneId-1];
}

Base::Base()
{
    this->nZones = -1;
    this->double_base_id = -1.0;
    this->celldim = -1;
    this->phydim = -1;
}

Base::~Base()
{
    for ( int iZone = 0; iZone < this->zones.size(); ++ iZone )
    {
        delete this->zones[ iZone ];
    }
}

void Base::AllocateZones( int nZones )
{
    this->nZones = nZones;
    for ( int iZone = 0; iZone < nZones; ++ iZone )
    {
        Zone * zone = new Zone();
        this->zones.push_back( zone );
    }
}

void Base::ReadBaseUnits( int fileId, int baseId )
{
    DataClass_t dclass;
    cg_goto(fileId, 1, "end");
    cg_dataclass_read(&dclass);

    std::cout << "dclass = " << dclass << " dclass name = " << DataClassName[dclass] << "\n";

    CGNS_ENUMT(MassUnits_t) mass;
    CGNS_ENUMT(LengthUnits_t) length;
    CGNS_ENUMT(TimeUnits_t) time;
    CGNS_ENUMT(TemperatureUnits_t) temp;
    CGNS_ENUMT(AngleUnits_t) angle;
    CGNS_ENUMT(ElectricCurrentUnits_t) current;
    CGNS_ENUMT(SubstanceAmountUnits_t) amount;
    CGNS_ENUMT(LuminousIntensityUnits_t) intensity;

    CGNS_ENUMT(MassUnits_t) im;
    CGNS_ENUMT(LengthUnits_t) il;
    CGNS_ENUMT(TimeUnits_t) it;
    CGNS_ENUMT(TemperatureUnits_t) ix;
    CGNS_ENUMT(AngleUnits_t) ia;

    cg_units_read(&im,&il,&it,&ix,&ia);

    std::printf("\nUnits=\n    %s\n    %s\n    %s\n    %s\n    %s\n",
                MassUnitsName[im],LengthUnitsName[il],TimeUnitsName[it],
                TemperatureUnitsName[ix],AngleUnitsName[ia]);

    int nunits = -1;
    cg_nunits(&nunits);

    std::cout << "   CGNS nunits = " << nunits << "\n";

    cg_unitsfull_read(&mass, &length, &time, &temp, &angle, &current, &amount, &intensity);
    std::cout << "mass = " << mass << " MassUnitsName[mass] = " << MassUnitsName[mass]<<"\n";
    std::cout << "length = " << length << " LengthUnitsName[length] = " << LengthUnitsName[length]<<"\n";
    std::cout << "time = " << time << " TimeUnitsName[time] = " << TimeUnitsName[time]<<"\n";
    std::cout << "temp = " << temp << " TemperatureUnitsName[temp] = " << TemperatureUnitsName[temp]<<"\n";
    std::cout << "angle = " << angle<< " AngleUnitsName[angle] = " << AngleUnitsName[angle]<<"\n";
    std::cout << "current = " << current<< " ElectricCurrentUnitsName[current] = " << ElectricCurrentUnitsName[current]<<"\n";
    std::cout << "amount = " << amount<< " SubstanceAmountUnitsName[amount] = " << SubstanceAmountUnitsName[amount]<<"\n";
    std::cout << "intensity = " << intensity << " LuminousIntensityUnitsName[intensity] = " << LuminousIntensityUnitsName[intensity]<<"\n";

    int nexps = -1;
    cg_nexponents (&nexps);
    std::cout << "nexps = " << nexps << "\n";
}

void Base::ReadBase( int fileId, int baseId )
{
    this->baseId = baseId;
    cg_base_id( fileId, baseId, & this->double_base_id );
    std::cout << "   double_base_id = " << this->double_base_id << "\n";
    //Check the cell and physical dimensions of the bases.
    cg_base_read( fileId, baseId, this->name, & celldim, & phydim );
    std::cout << "   baseId = " << baseId << " baseName = " << this->name << "\n";
    std::cout << "   cell dim = " << this->celldim << " physical dim = " << this->phydim << "\n";

    //Read the number of zones in the grid.
    int nZones = -1;
    cg_nzones( fileId, baseId, & nZones );

    std::cout << "** Reading CGNS Grid In Base " << baseId << "\n";
    std::cout << "   Reading CGNS Family Specified BC \n";
    //ReadFamilySpecifiedBc();
    std::cout << "   numberOfCgnsZones       = " << nZones << "\n\n";

    this->ReadBaseUnits( fileId, baseId );

    this->AllocateZones( nZones );

    for ( int iZone = 0; iZone < nZones; ++ iZone )
    {
        std::cout << "==>iZone = " << iZone << " numberOfCgnsZones = " << nZones << "\n";
        int zoneId = iZone + 1;
        Zone * zone = this->zones[iZone];
        zone->ReadZone(fileId, baseId, zoneId);
        zone->DumpZone();
    }
}

Zone::Zone()
{
    this->zoneType = ZoneTypeNull;
    this->nNodes = -1;
    this->nCells = -1;
    this->nFaces = -1;
    this->nCoords = -1;
    this->nSections = -1;
    this->bc = new Bc();
}

Zone::~Zone()
{
    this->DeAllocateCoors();
    this->DeAllocateSections();
    delete this->bc;
}

void Zone::ReadZone( int fileId, int baseId, int zoneId )
{
    this->zoneId = zoneId;
    //Check the zone type
    cg_zone_type( fileId, baseId, zoneId, & this->zoneType );

    std::cout << "   The Zone Type is " << ZoneTypeName[ this->zoneType ] << " Zone" << "\n";

    //Determine the number of vertices and cellVolume elements in this zone
    cg_zone_read( fileId, baseId, zoneId, this->zoneName, this->isize );

    std::cout << "   CGNS Zone Name = " << this->zoneName << "\n";

    this->SetDimensions();

    this->ReadCoordinates( fileId, baseId, zoneId );
    this->ReadElements( fileId, baseId, zoneId );
    this->ReadBoundaries( fileId, baseId, zoneId );
    this->ReadFlowSolution( fileId, baseId, zoneId );
}

void Zone::AllocateSolutions( int nSolutions )
{
    for ( int iSolution = 0; iSolution < nSolutions; ++ iSolution )
    {
        Solution * solution = new Solution();
        this->solutions.push_back( solution );
    }
}

void Zone::ReadFlowSolution( int fileId, int baseId, int zoneId )
{
    int result = cg_nsols(fileId, baseId, zoneId, &this->nSolutions);
    if ( result != CG_OK )
    {
        std::cout << "Get Num of Solution Error." << cg_get_error() << std::endl;
    }
    else
    {
        std::cout << "   nSolutions = " << this->nSolutions << std::endl;
    }

    this->AllocateSolutions( this->nSolutions );

    for ( int iSolution = 0; iSolution < this->nSolutions; ++ iSolution )
    {
        int solutionId = iSolution + 1;
        std::cout << "Zone::ReadFlowSolution solutionId = " << solutionId << std::endl;
        Solution * solution = this->solutions[iSolution];
        solution->ReadField( fileId, baseId, zoneId, solutionId, this );
    }
}

void Zone::DumpZone()
{
    this->nFaces = 0;
    this->totalNumFaceNodes = 0;

    for ( int iSection = 0; iSection < this->nSections; ++ iSection )
    {
        Section * section = this->sections[iSection];
        int sectionId = iSection + 1;
        if ( section->elementType == NGON_n )
        {
            this->nFaces += section->nElements;
            totalNumFaceNodes += section->elementDataSize;
        }
    }
    std::cout << "total faces of zone : " <<  this->nFaces << "\n";
    std::cout << "total number face nodes of zone : " <<  this->totalNumFaceNodes << "\n";
    this->DrawZone();
}

void Zone::DrawZone()
{
    std::cout << "   DrawZone Zone Name  = " << this->zoneName << "\n";
    std::string  fileName = this->zoneName;
    std::replace( fileName.begin(), fileName.end(), ':', '-');
    std::replace( fileName.begin(), fileName.end(), ' ', '-');
    fileName += ".plt";
    std::cout << "   fileName  = " << fileName << "\n";

    std::fstream file;
    file.open(fileName.c_str(), std::ios::out);

    Plot plot;
    std::ostringstream oss;
    plot.PlotZoneMesh( oss, this );

    this->DumpCoor( oss );
    this->DumpFaceNodeNumber( oss );
    this->DumpFaceNodeLink( oss );

    this->CalcFaceElementLink();
    this->DumpFaceElementLink( oss, this->left_elements );
    this->DumpFaceElementLink( oss, this->right_elements );

    this->DumpSectionMesh( oss );

    //std::cout << oss.str() << "\n";

    file << oss.str() << "\n";
    file.close();
    file.clear();
}

void Zone::DumpSectionMesh( std::ostringstream & oss )
{
    for ( int iSection = 0; iSection < this->nSections; ++ iSection )
    {
        Section * section = this->sections[iSection];
        if ( section->elementType == NGON_n )
        {
            section->zone = this;
            section->DumpSectionMesh( oss );
        }
    }
}

void Zone::DumpCoor( std::ostringstream & oss )
{
    int nWords = 5;
    for ( int iCoord = 0; iCoord < this->nCoords; ++ iCoord )
    {
        Coor * coor = this->coors[iCoord];
        if ( coor->dataType == RealSingle )
        {
            float * field = static_cast<float *>( coor->data );
            ::DumpCoor( oss, field, this->nNodes, nWords );
        }
        else
        {
            double * field = static_cast<double *>( coor->data );
            ::DumpCoor( oss, field, this->nNodes, nWords );
        }
    }
}

void Zone::DumpFaceNodeNumber( std::ostringstream & oss )
{
    int icount = 0;
    int nWords = 10;
    for ( int iSection = 0; iSection < this->nSections; ++ iSection )
    {
        Section * section = this->sections[iSection];
        if ( section->elementType == NGON_n )
        {
            for ( int iElem = 0; iElem < section->nElements; ++ iElem )
            {
                int st = section->conn_offsets[ iElem ];
                int ed = section->conn_offsets[ iElem + 1 ];
                int nNode = ed - st;
                oss <<  std::setw(8) << nNode;
                icount ++;
                if ( icount % nWords == 0 ) oss << std::endl;
            }
        }
    }
    if ( icount % nWords != 0 ) oss << std::endl;
}

void Zone::DumpFaceNodeLink( std::ostringstream & oss )
{
    int icount = 0;
    int nWords = 10;
    for ( int iSection = 0; iSection < this->nSections; ++ iSection )
    {
        Section * section = this->sections[iSection];
        if ( section->elementType == NGON_n )
        {
            for ( int i = 0; i < section->conn.size(); ++ i )
            {
                oss <<  std::setw(8) << section->conn[i];
                icount ++;
                if ( icount % nWords == 0 ) oss << std::endl;
            }
        }
    }
    if ( icount % nWords != 0 ) oss << std::endl;
}

void Zone::DumpFaceElementLink( std::ostringstream & oss, std::vector<cgsize_t> &elementId )
{
    int icount = 0;
    int nWords = 10;
    for ( int iFace = 0; iFace < this->nFaces; ++ iFace )
    {
        int eId = elementId[ iFace ];
        oss << eId << " ";
        icount ++;
        if ( icount % nWords == 0 ) oss << std::endl;
    }
    if ( icount % nWords != 0 ) oss << std::endl;
}

void Zone::CalcFaceElementLink()
{
    this->left_elements.resize(this->nFaces,0);
    this->right_elements.resize(this->nFaces,0);

    for ( int iSection = 0; iSection < this->nSections; ++ iSection )
    {
        Section * section = this->sections[iSection];
        if ( section->elementType != NFACE_n ) continue;
        for ( int iElem = 0; iElem < section->nElements; ++ iElem )
        {
            int st = section->conn_offsets[ iElem ];
            int ed = section->conn_offsets[ iElem + 1 ];
            int elemId = iElem + 1;
            for ( int i = st; i < ed; ++ i )
            {
                int faceIndex = section->conn[ i ];
                int fIndex = std::abs( faceIndex ) - 1;
                if ( faceIndex > 0 )
                {
                    left_elements[ fIndex ] = elemId;
                }
                else
                {
                    right_elements[ fIndex ] = elemId;
                }
            }
        }
    }
}


void Zone::SetDimensions()
{
    if ( this->zoneType == CGNS_ENUMV( Unstructured ) )
    {
        this->irmin[ 0 ] = 1;
        this->irmin[ 1 ] = 0;
        this->irmin[ 2 ] = 0;

        this->irmax[ 0 ] = this->isize[ 0 ];
        this->irmax[ 1 ] = 0;
        this->irmax[ 2 ] = 0;

        this->cellSize[ 0 ] = this->isize[ 1 ];

        this->nNodes = this->irmax[ 0 ];
        this->nCells = this->cellSize[ 0 ];
    }
    std::cout << "   numberOfNodes = " << this->nNodes << " numberOfCells = " << this->nCells << "\n";
}


void Zone::ReadCoordinates( int fileId, int baseId, int zoneId )
{
    cg_ncoords( fileId, baseId, zoneId, & this->nCoords );
    std::cout << "   nCoords = " << this->nCoords << "\n";
    this->AllocateCoors( this->nCoords );

    for ( int iCoord = 0; iCoord < this->nCoords; ++ iCoord )
    {
        int coordId = iCoord + 1;
        Coor * coor = this->coors[iCoord];
        coor->ReadCoor( fileId, baseId, zoneId, coordId, this );
    }
}

void Zone::ReadElements( int fileId, int baseId, int zoneId )
{
    // Determine the number of sections for this zone. Note that
    // surface elements can be stored in a cellVolume zone, but they
    // are NOT taken into account in the number obtained from
    // cg_zone_read.

    cg_nsections( fileId, baseId, zoneId, & this->nSections );

    std::cout << "   numberOfCgnsSections = " << this->nSections << "\n";
    this->AllocateSections(this->nSections);

    std::cout << "   Reading Cgns Section Data......\n";
    std::cout << "\n";

    for ( int iSection = 0; iSection < this->nSections; ++ iSection )
    {
        std::cout << "-->iSection     = " << iSection << " numberOfCgnsSections = " << this->nSections << "\n";
        Section * section = this->sections[iSection];
        int sectionId = iSection + 1;
        section->Read( fileId, baseId, zoneId, sectionId );
    }
}

void Zone::ReadBoundaries( int fileId, int baseId, int zoneId )
{
    this->bc->ReadBoundaries( fileId, baseId, zoneId );
}

void Zone::AllocateCoors( int nCoords )
{
    for ( int iCoord = 0; iCoord < nCoords; ++ iCoord )
    {
        Coor * coor = new Coor();
        this->coors.push_back( coor );
    }
}

void Zone::DeAllocateCoors()
{
    for ( int iCoor = 0; iCoor < this->coors.size(); ++ iCoor )
    {
        delete this->coors[ iCoor ];
    }
}

void Zone::AllocateSections( int nSections )
{
    for ( int iSection = 0; iSection < this->nSections; ++ iSection )
    {
        Section * section = new Section();
        this->sections.push_back( section );
    }
}

void Zone::DeAllocateSections()
{
    for ( int iSection = 0; iSection < this->sections.size(); ++ iSection )
    {
        delete this->sections[ iSection ];
    }
}

Coor::Coor()
{
    this->data = 0;
}

Coor::~Coor()
{
    this->DeAllocateData();
}

void Coor::ReadCoor( int fileId, int baseId, int zoneId, int coordId, Zone * zone )
{
    cg_coord_info( fileId, baseId, zoneId, coordId, & this->dataType, this->coorName );
    std::cout << "   coordId = " << coordId << " coorName = " << this->coorName << " dataType = " << dataType << " dataTypeName = " << DataTypeName[ dataType ] << "\n";
    this->AllocateData( zone->nNodes );
    //Read coordinates.
    cg_coord_read( fileId, baseId, zoneId, this->coorName, this->dataType, zone->irmin, zone->irmax, this->data );
    //this->DumpCoor();
}

void Coor::AllocateData( int nNodes )
{
    this->nNodes = nNodes;
    if ( this->dataType == RealSingle )
    {
        this->data = new float [ nNodes ];
    }
    else
    {
        this->data = new double [ nNodes ];
    }
}

void Coor::DeAllocateData()
{
    if ( this->dataType == RealSingle )
    {
        float * f_data  = static_cast< float * >( this->data );
        delete [] f_data;
    }
    else
    {
        double * d_data  = static_cast< double * >( this->data );
        delete [] d_data;
    }
}

void Coor::DumpCoor()
{
    if ( this->dataType == RealSingle )
    {
        int icount = 0;
        for ( int iNode = 0; iNode < this->nNodes; ++ iNode )
        {
            icount ++;
            std::cout << (static_cast<float *>(this->data))[ iNode ] << " ";
            if ( icount % 5 == 0 )  std::cout << "\n";
        }
    }
    else
    {
        int icount = 0;
        for ( int iNode = 0; iNode < this->nNodes; ++ iNode )
        {
            icount ++;
            std::cout << (static_cast<double *>(this->data))[ iNode ] << " ";
            if ( icount % 5 == 0 )  std::cout << "\n";
        }
    }
}

Field::Field()
{
    this->data = 0;
}

Field::~Field()
{
    this->DeAllocateData();
}

void Field::ReadField( int fileId, int baseId, int zoneId, int solutionId, int fieldId, Zone * zone )
{
    /* read DimensionalExponents */
    //cg_goto(fileId, baseId, "Zone_t", 1, "FlowSolution_t", 1, "DataArray_t", fieldId, "end");
    //cg_goto(fileId, baseId, zoneId, "FlowSolution_t", 1, "DataArray_t", fieldId, "end");

    std::cout << "Field::ReadField" << " solutionId = " << solutionId << " fieldId = " << fieldId << std::endl;
    cg_field_info(fileId, baseId, zoneId, solutionId, fieldId, &this->dataType, this->fieldName);
    std::cout << "dataType Value = " << this->dataType << " dataType Name = " << DataTypeName[dataType] << " ";
    std::cout << "fieldName = " << this->fieldName << std::endl;

    cg_field_id(fileId, baseId, zoneId, solutionId, fieldId, &this->doubleFieldId);
    std::cout << "doubleFieldId = " << this->doubleFieldId << std::endl;
    std::cout << "zone->nNodes = " << zone->nNodes << std::endl;

    std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << std::endl;

    this->AllocateData( zone->nNodes );

    std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << std::endl;
    int ierr = cg_field_read( fileId, baseId, zoneId, solutionId, this->fieldName, this->dataType, zone->irmin, zone->irmax, this->data);
    if ( ierr == CG_OK )
    {
        std::cout << "cg_field_read ierr = " << ierr << "\n";
    }

    cg_goto(fileId, baseId, "Zone_t", zoneId, "FlowSolution_t", solutionId, "DataArray_t", fieldId, "end");

    this->nExponents = -1;
    ierr = cg_nexponents (&this->nExponents);

    std::cout << "fieldId = " << fieldId << "  nExponents = " << nExponents << "\n";

    if ( ierr == CG_OK && this->nExponents > 0 )
    {
        this->dimflag = 1;
        this->exponents.resize(this->nExponents);
        std::cout << "cg_nexponents ierr = " << ierr << "\n";
        DataType_t expDataType;
        ierr = cg_exponents_info(&expDataType);
        std::cout << "cg_exponents_info ierr = " << ierr << "\n";
        if ( ierr == CG_OK )
        {
            std::cout << "ierr = " << ierr << "\n";
            std::cout << "expDataType value = " << expDataType << " Name = " << DataTypeName[expDataType] << "\n";
        }
        cg_exponents_read(exponents.data());
        std::cout << "exponents = \n";
        for ( int iexp = 0; iexp < this->nExponents; ++ iexp )
        {
            std::cout << exponents[ iexp ] << " ";
        }
        std::cout << std::endl;
    }
    else
    {
        //DataClass_t dataclass;
        ierr = cg_dataclass_read(&this->dataclass);
        if ( ierr == CG_OK )
        {
            this->dimflag = 0;
            std::cout << "dataclass value = " << dataclass << " Name = " << DataTypeName[this->dataclass] << "\n";
        }
        else
        {
            std::cout << "-----------------------------------\n";
            std::cout << "fieldName = " << this->fieldName << std::endl;
            std::cout << "-----------------------------------\n";
        }
    }

}

void Field::AllocateData( int nNodes )
{
    this->nNodes = nNodes;
    if ( this->dataType == RealSingle )
    {
        this->data = new float [ nNodes ];
    }
    else
    {
        this->data = new double [ nNodes ];
    }
}

void Field::DeAllocateData()
{
    if ( this->dataType == RealSingle )
    {
        float * f_data  = static_cast< float * >( this->data );
        delete [] f_data;
    }
    else
    {
        double * d_data  = static_cast< double * >( this->data );
        delete [] d_data;
    }
}

Solution::Solution()
{
    ;
}

Solution::~Solution()
{
    ;
}

void Solution::AllocateFields( int nFields )
{
    for ( int iField = 0; iField < nFields; ++ iField )
    {
        Field * field = new Field();
        this->fields.push_back( field );
    }
}

void Solution::ReadField( int fileId, int baseId, int zoneId, int solutionId, Zone * zone )
{
    cg_sol_info(fileId, baseId, zoneId, solutionId, this->solutionName, &this->location);
    std::cout << "solutionName = " << this->solutionName << std::endl;

    cg_sol_id(fileId, baseId, zoneId, solutionId, &this->solver_double_id);
    std::cout << "solver_double_id = " << solver_double_id << std::endl;
    int data_dim = -1;
    cgsize_t dim_vals[3];
    cg_sol_size(fileId, baseId, zoneId, solutionId, &data_dim, dim_vals);
    std::cout << "data_dim = " << data_dim << std::endl;

    int result = cg_nfields(fileId, baseId, zoneId, solutionId, &this->nFields);
    if ( result != CG_OK )
    {
        std::cout << "cg_nfields Error:" << cg_get_error() << std::endl;
    }
    else
    {
        std::cout << "nFields = " << this->nFields << std::endl;
    }

    for ( int iField = 0; iField < this->nFields; ++ iField )
    {
        Field * field = new Field();
        this->fields.push_back( field );

        int fieldId = iField + 1;

        field->ReadField( fileId, baseId, zoneId, solutionId, fieldId, zone );
    }
}


Section::Section()
{
    this->startId = -1;
    this->endId = -1;
    this->nElements = -1;
    this->nbndry = -1;
    this->iparentflag = -1;
    this->elementDataSize = -1;
    this->pos_shift = -1;
    this->zone = 0;
}

Section::~Section()
{
    ;
}

void Section::Read( int fileId, int baseId, int zoneId, int sectionId )
{
    cg_section_read( fileId, baseId, zoneId, sectionId, this->sectionName, & this->elementType, & this->startId, & this->endId, & this->nbndry, & this->iparentflag );
    this->nElements = this->endId - this->startId + 1;

    std::cout << "   Section Name  = " << this->sectionName << "\n";
    std::cout << "   Section Type  = " << ElementTypeName[ this->elementType ] << "\n";
    std::cout << "   startId, endId = " << this->startId << " " << this->endId << "\n";
    std::cout << "   nElements      = " << this->nElements << "\n";

    cg_ElementDataSize( fileId, baseId, zoneId, sectionId, & this->elementDataSize );
    std::cout << "   elementDataSize = " << this->elementDataSize << "\n";

    this->conn.resize( this->elementDataSize );

    if ( this->IsMixedSection() )
    {
        this->pos_shift = 1;
        this->conn_offsets.resize(this->nElements+1);
        cg_poly_elements_read ( fileId, baseId, zoneId, sectionId, this->conn.data(), this->conn_offsets.data(), 0 );
        //this->PrintElementInfo();
    }
}

bool Section::IsMixedSection()
{
    bool flag = ( this->elementType == MIXED ||
                 this->elementType == NGON_n ||
                 this->elementType == NFACE_n );
    return flag;
}

void Section::PrintElementInfo()
{
    if ( this->elementType == NGON_n )
    {
        this->PrintPolygonFaceInfo();
    }
    else if ( this->elementType == NFACE_n )
    {
        this->PrintPolyhedronElementInfo();
    }
}

void Section::PrintPolygonFaceInfo()
{
    std::cout << "  Polygon Face Info:\n";
    for ( int iElem = 0; iElem < this->nElements; ++ iElem )
    {
        std::cout << "   iFace = " << iElem + 1 << " , ";
        int iFace = iElem;
        int st = this->conn_offsets[ iFace ];
        int ed = this->conn_offsets[ iFace + 1 ];
        int nNode = ed - st;
        std::cout << " nNode = " << nNode << " Index = ";
        for ( int i = st; i < ed; ++ i )
        {
            std::cout << this->conn[ i ] << " ";
        }
        std::cout << "\n";
    }
}

void Section::PrintPolyhedronElementInfo()
{
    std::cout << "   Polyhedron Element Info:\n";
    int icount = 0;
    for ( int iElem = 0; iElem < this->nElements; ++ iElem )
    {
        std::cout << "   iElem = " << iElem + 1 << " , ";
        int st = this->conn_offsets[ iElem ];
        int ed = this->conn_offsets[ iElem + 1 ];
        int nFace = ed - st;
        std::cout << "nFace = " << nFace << " Index = ";
        for ( int i = st; i < ed; ++ i )
        {
            std::cout << this->conn[ i ] << " ";
        }
        std::cout << "\n";
    }
}

int Section::GetNumFaceNodes()
{
    int numFaceNodes = 0;
    for ( int iElem = 0; iElem < this->nElements; ++ iElem )
    {
        int st = this->conn_offsets[ iElem ];
        int ed = this->conn_offsets[ iElem + 1 ];
        int nNode = ed - st;
        numFaceNodes += nNode;
    }
    std::cout << "numFaceNodes = " << numFaceNodes << " elementDataSize = " << elementDataSize << "\n";
    return numFaceNodes;
}

void Section::DumpSectionMesh( std::ostringstream & oss )
{
    std::cout << "   DumpSectionMesh Section Name  = " << this->sectionName << "\n";
    std::string  fileName = this->sectionName;
    if ( fileName.substr( 0, 9 ) == "Polyfaces" ) return;
    std::replace( fileName.begin(), fileName.end(), ':', '-');
    std::replace( fileName.begin(), fileName.end(), ' ', '-');
    fileName += "-section.plt";
    std::cout << "   DumpSectionMesh fileName  = " << fileName << "\n";

    std::fstream file;
    file.open(fileName.c_str(), std::ios::out);

    std::ostringstream hoss;
    std::ostringstream myoss;

    DumpTecplotHeader( hoss );

    Uns2D uns2d;
    uns2d.CalcUnsFaceInfo( this );
    uns2d.PlotSectionMesh( myoss, this );
    uns2d.DumpFaceNodeLink( myoss, this );
    uns2d.DumpFaceElementLink( myoss, this );

    oss << myoss.str();
    file << hoss.str();
    file << myoss.str();
    file.close();
}
