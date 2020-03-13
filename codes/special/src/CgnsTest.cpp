/*-----------------------this->----------------------------------------------------*\
    OneFLOW - LargeScale Multiphysics Scientific Simulation Environment
    Copyright (C) 2017-2020 He Xin and the OneFLOW contributors.
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

#include "CgnsTest.h"
#include "CgnsFactory.h"
#include "Prj.h"
#include "CgnsZone.h"
#include <iostream>
#include <iomanip>
using namespace std;

BeginNameSpace( ONEFLOW )

CgnsTest::CgnsTest()
{
    ;
}

CgnsTest::~CgnsTest()
{
    ;
}

void CgnsTest::Init()
{
    this->SetDefaultGridName();
}

void CgnsTest::Run()
{
    this->Init();

    //this->WriteSimpleMultiBaseTest();
    //this->ReadSimpleMultiBaseTest();
    //this->WriteEmptyCgnsFile();
    //this->ReadEmptyCgnsFile();
    //this->WriteDescriptor();
    this->ReadDescriptor();
    //this->WriteNondimensionalParameter();
    //this->ReadNondimensionalParameter();
 }

void CgnsTest::Test()
{
}

void CgnsTest::SetDefaultGridName()
{
    string gridName = "/grid/oneflow.cgns";
    string prjFileName = ONEFLOW::GetPrjFileName( gridName );
    cout << " CGNS File Name = " << prjFileName << "\n";

    this->fileName = prjFileName;
}

void CgnsTest::WriteBase( const string & baseName )
{
    int celldim = 3;
    int physdim = 3;
    this->WriteBase( baseName, celldim, physdim );
}

void CgnsTest::WriteBase( const string & baseName, int celldim, int physdim )
{
    int index_base = -1;
    cg_base_write( index_file, baseName.c_str(), celldim, physdim, & index_base );
    cout << " CGNS Base index = " << index_base << "\n";
    curr_base_id = index_base;
}

void CgnsTest::WriteSimpleMultiBaseTest()
{
    this->OpenCgnsFile( CG_MODE_WRITE );

    //this->WriteBase( "base1" );
    //this->WriteBase( "base2" );
    //this->WriteBase( "base3" );
    //this->WriteBase( "base4" );

    //this->WriteBase( "base1" );
    //this->WriteBase( "base2", 2, 3 );
    //this->WriteBase( "base3", 3, 2 );
    //this->WriteBase( "base4", 1, 3 );

    this->WriteBase( "OneFLOW1" );
    this->WriteBase( "OneFLOW 2" );
    this->WriteBase( "CGNS base 3" );
    this->WriteBase( "Fluid" );
    this->WriteBase( "CAE library" );

    this->CloseCgnsFile();
}

void CgnsTest::ReadSimpleMultiBaseTest()
{
    this->OpenCgnsFile( CG_MODE_READ );

    //Determine the of bases in the grid
    cg_nbases( this->index_file, & this->nBases );
    cout << " Total number of CGNS Base = " << this->nBases << "\n";

    for ( int iBase = 0; iBase < this->nBases; ++ iBase )
    {
        double double_base_id;
        int baseId = iBase + 1;
        char basename[ 33 ];
        int cell_dim, phys_dim;
        cg_base_id( this->index_file, baseId, & double_base_id );
        cg_base_read( this->index_file, baseId, basename, & cell_dim, & phys_dim );
        cout << " baseId = " << baseId << " double_base_id = " << double_base_id << " basename = " << basename << "\n";
    }

    this->CloseCgnsFile();
}

void CgnsTest::WriteDescriptor()
{
    this->OpenCgnsFile( CG_MODE_WRITE );

    cout << "Program write_descriptor\n";
    this->WriteBase( "base1" );
    cout << " curr_base_id = " << curr_base_id << "\n";

    //cg_goto must be called or an error will occur
    cg_goto( index_file, curr_base_id, "end" );
    cg_descriptor_write("Information","info1");
    cg_descriptor_write("hello world","haha ! hello world!");

    this->WriteBase( "base2" );
    cout << " curr_base_id = " << curr_base_id << "\n";
    cg_goto( index_file, curr_base_id, "end" );
    cg_descriptor_write("descript1","des1");
    cg_descriptor_write("descript2","des1");
    cg_descriptor_write("descript3","des1");

    this->WriteBase( "base3" );
    cout << " curr_base_id = " << curr_base_id << "\n";
    cg_goto( index_file, curr_base_id, "end" );
    cg_descriptor_write("mydes","mydes1");
    cg_descriptor_write("mydes","mydes2");
    cg_descriptor_write("mydes","mydes3");

    this->CloseCgnsFile();
}

void CgnsTest::ReadBaseDescriptor( int baseIndex )
{
    //go to base node
    cg_goto( this->index_file, baseIndex, "end" );

    //find out how many descriptors are here:
    int ndescriptors = -1;
    cg_ndescriptors( & ndescriptors );
    cout << " ndescriptors = " << ndescriptors << "\n";
    for ( int n = 1; n <= ndescriptors; ++ n )
    {
        //read descriptor
        char *text, name[33];
        cg_descriptor_read( n, name, &text );
        cout << "The descriptor is : " << name << "," << text << "\n";
    }
}

void CgnsTest::ReadDescriptor()
{
    this->OpenCgnsFile( CG_MODE_READ );

    //Determine the of bases in the grid
    cg_nbases( this->index_file, & this->nBases );
    cout << " Total number of CGNS Base = " << this->nBases << "\n";

    for ( int iBase = 0; iBase < this->nBases; ++ iBase )
    {
        int baseId = iBase + 1;
        ReadBaseDescriptor( baseId );
    }

    this->CloseCgnsFile();
}

string CgnsTest::GetCgnsFileTypeName( int file_type )
{
    string fileTypeName;
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


void CgnsTest::OpenCgnsFile( int cgnsOpenMode )
{
    this->OpenCgnsFile( this->fileName, cgnsOpenMode );
}

void CgnsTest::OpenCgnsFile( const string & fileName, int cgnsOpenMode )
{
    int result = cg_open( fileName.c_str(), cgnsOpenMode, & index_file );
    cout << " CGNS File Index = " << index_file << "\n";
    if ( result != CG_OK )
    {
        cg_error_exit();
    }
}

void CgnsTest::CloseCgnsFile()
{
    cg_close( index_file );
}

void CgnsTest::WriteEmptyCgnsFile()
{
    this->OpenCgnsFile( CG_MODE_WRITE );
    this->CloseCgnsFile();
}

void CgnsTest::ReadEmptyCgnsFile()
{
    this->OpenCgnsFile( CG_MODE_READ );

    float fileVersion = -1;
    cg_version( index_file, & fileVersion );

    cout << " CGNS File Version = " << setiosflags( ios::fixed ) << setprecision( 4 ) << fileVersion << "\n";

    int precision = -1;
    cg_precision( index_file, & precision );

    cout << " CGNS Precision = " << precision << "\n";

    int file_type = -1;
    cg_get_file_type( index_file, & file_type );

    cout << " CGNS file_type = " << file_type << " file_type name = " << GetCgnsFileTypeName( file_type ) << "\n";

    this->CloseCgnsFile();
}

void CgnsTest::WriteDouble( const string & varName, const double & varValue )
{
    int nDim = 1;
    cgsize_t ndims[ 1 ] = { 1 };
    cg_array_write( varName.c_str(), CGNS_ENUMV(RealDouble), nDim ,ndims, &varValue );
}

void CgnsTest::WriteNondimensionalParameter()
{
    double xmach,reue,xmv,xmc,rev,rel,renu,rho0;
    double p0,c0,vm0,xlength0,vx,vy,vz;
    double gamma;
    //int index_file,index_base;
    CGNS_ENUMT(DataClass_t) idata;

    cout << "\n";
    cout << "Program write_nondimensional\n";

    //define nondimensional parameters
    xmach    = 4.6;
    reue     = 6000000.;
    xmv      = xmach;
    xmc      = 1.;
    rev      = xmach;
    rel      = 1.;
    renu     = xmach / reue;
    rho0     = 1.;
    gamma    = 1.4;
    p0       = 1./gamma;
    c0       = 1.;
    vm0      = xmach/reue;
    xlength0 = 1.;
    vx       = xmach;
    vy       = 0.0;
    vz       = 0.0;

    this->OpenCgnsFile( CG_MODE_WRITE );

    this->WriteBase( "base1" );

    int index_base = this->curr_base_id;

    //put DataClass under Base
    cg_goto(index_file,index_base,"end");
    cg_dataclass_write( CGNS_ENUMV( NormalizedByUnknownDimensional ) );

    //put ReferenceState under Base
    cg_state_write("ReferenceQuantities");

    //Go to ReferenceState node, write Mach array and its dataclass
    cg_goto(index_file,index_base,"ReferenceState_t",1,"end");

    WriteDouble("Mach", xmach );

    cg_goto(index_file,index_base,"ReferenceState_t",1,"DataArray_t",1,"end");

    cg_dataclass_write( CGNS_ENUMV(NondimensionalParameter) );

    //Go to ReferenceState node, write Reynolds array and its dataclass
    cg_goto(index_file,index_base,"ReferenceState_t",1,"end");
    WriteDouble("Reynolds", reue );
    cg_goto(index_file,index_base,"ReferenceState_t",1,"DataArray_t",2,"end");
    cg_dataclass_write( CGNS_ENUMV(NondimensionalParameter) );
    //Go to ReferenceState node to write reference quantities
    cg_goto(index_file,index_base,"ReferenceState_t",1,"end");

    // First, write reference quantities that make up Mach and Reynolds:
    WriteDouble("Mach_Velocity", xmv );
    WriteDouble("Mach_VelocitySound", xmc );
    WriteDouble("Reynolds_Velocity", rev );
    WriteDouble("Reynolds_Length", rel );
    WriteDouble("Reynolds_ViscosityKinematic", renu );

    //Next, write flow field reference quantities:
    WriteDouble("Density", rho0 );
    WriteDouble("Pressure", p0 );
    WriteDouble("VelocitySound", c0 );
    WriteDouble("ViscosityMolecular", vm0 );
    WriteDouble("LengthReference", xlength0 );

    WriteDouble("VelocityX", vx );
    WriteDouble("VelocityY", vy );
    WriteDouble("VelocityZ", vz );

    WriteDouble("aoa", 5.0 );

    this->CloseCgnsFile();
}

void CgnsTest::GotoBaseBegin( int baseIndex )
{
    cg_goto(index_file,baseIndex,"end");
}

void CgnsTest::ReadNondimensionalParameter()
{
    int narrays, idim;
    char *state, arrayname[33];
    
    this->OpenCgnsFile( CG_MODE_READ );

    //Determine the of bases in the grid
    cg_nbases( this->index_file, & this->nBases );

    cout << " Total number of CGNS Base = " << this->nBases << "\n";

    for ( int iBase = 0; iBase < this->nBases; ++ iBase )
    {
        int baseId = iBase + 1;
        int index_base = 1;
        this->GotoBaseBegin( index_base );

        CGNS_ENUMT(DataClass_t) id;
        cg_dataclass_read( & id );
        cout << "DataClass = " << DataClassName[ id ] << "\n";

        cg_state_read( & state );
        cout << "ReferenceState = " << state << "\n";

        //Go to ReferenceState node, read Mach array and its dataclass
        cg_goto( index_file, index_base, "ReferenceState_t", 1, "end");
        //find out how many data arrays 
        cg_narrays( & narrays );
        cout << " narrays = " << narrays << "\n";
        for ( int n = 1; n <= narrays; ++ n )
        {
            CGNS_ENUMT(DataType_t) idata;
            cgsize_t idimvec;
            cg_array_info( n, arrayname, & idata, & idim, & idimvec );
            //cout << " DataTypeName = " << DataTypeName[ idata ] << "\n";
            double data;
            cg_array_read_as( n, CGNS_ENUMV(RealDouble), & data );
            cout << "Variable = " << arrayname << "\n";
            cout << "   data = " << data << "\n";
        }
    }

    this->CloseCgnsFile();
}



EndNameSpace