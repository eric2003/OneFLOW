/*---------------------------------------------------------------------------*\
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

void CgnsTest::Run()
{
    //CgnsFactory * cgnsFactory = new CgnsFactory();

    //string filename = "test";
    //delete cgnsFactory;

    //this->Test();

    this->ReadNondimensionalParameter();
    //this->WriteDescriptor();
 }

void CgnsTest::Test()
{
}

void CgnsTest::ReadNondimensionalParameter()
{
    string gridName = "/grid/grid_c.cgns";
    string prjFileName = ONEFLOW::GetPrjFileName( gridName );
    cout << " prjFileName = " << prjFileName << "\n";

    int index_file, index_base;

    double data;
    int narrays, n, idim;
    char *state,arrayname[33];
    CGNS_ENUMT(DataClass_t) id;
    CGNS_ENUMT(DataType_t) idata;
    cgsize_t idimvec;

    if ( cg_open( prjFileName.c_str(), CG_MODE_READ, &index_file ) != CG_OK )
    {
        cg_error_exit();
    }

    float fileVersion = -1;

    cg_version( index_file, & fileVersion );

    int precision = -1;
    cg_precision( index_file, & precision );

    int file_type = -1;
    cg_get_file_type( index_file, & file_type );

    int ndescriptors = -1;
    cg_ndescriptors( & ndescriptors );

    //#define CG_FILE_NONE  0
    //#define CG_FILE_ADF   1
    //#define CG_FILE_HDF5  2
    //#define CG_FILE_ADF2  3

    index_base=1;
    /* read DataClass under Base */
    cg_goto(index_file,index_base,"end");
    cg_dataclass_read( & id );

    cout << "DataClass = " << DataClassName[ id ] << "\n";
    if ( id != CGNS_ENUMV( NormalizedByUnknownDimensional ) )
    {
        cout << "Error!  Expecting NormalizedByUnknownDimensional\n";
        return;
    }
    /* read ReferenceState under Base */
    cg_state_read( & state );
    cout << "ReferenceState = " << state << "\n";
    /* Go to ReferenceState node, read Mach array and its dataclass */
    cg_goto(index_file,index_base,"ReferenceState_t",1,"end");
    /* find out how many data arrays */
    cg_narrays( &narrays );
    for ( n=1; n <= narrays; n ++ )
    {
        cg_array_info(n,arrayname,&idata,&idim,&idimvec);
        if (idim != 1 || idimvec != 1)
        {
            cout << "Error! expecting idim,idimvec=1,1\n";
            cout << "   they are idim,idimvec= " << idim << "," << (int)(idimvec) << "\n";
            return;
        }
        cg_array_read_as(n,CGNS_ENUMV(RealDouble),&data);
        cout << "Variable = " << arrayname << "\n";
        cout << "   data = " << data << "\n";
    }

    cg_close(index_file);
}

void CgnsTest::WriteNondimensionalParameter()
{
    double xmach,reue,xmv,xmc,rev,rel,renu,rho0;
    double p0,c0,vm0,xlength0,vx,vy,vz;
    double gamma;
    int index_file,index_base;
    CGNS_ENUMT(DataClass_t) idata;
    cgsize_t nuse;

    cout << "\n";
    cout << "Program write_nondimensional\n";

    //define nondimensional parameters
    xmach=4.6;
    reue=6000000.;
    xmv=xmach;
    xmc=1.;
    rev=xmach;
    rel=1.;
    renu=xmach/reue;
    rho0=1.;
    gamma=1.4;
    p0=1./gamma;
    c0=1.;
    vm0=xmach/reue;
    xlength0=1.;
    vx=xmach;
    vy=0.;
    vz=0.;
    nuse=1;
    // WRITE NONDIMENSIONAL INFO
    /* open CGNS file for modify */
    string gridName = "/grid/grid_c.cgns";
    string prjFileName = ONEFLOW::GetPrjFileName( gridName );
    cout << " prjFileName = " << prjFileName << "\n";

    if ( cg_open( prjFileName.c_str(), CG_MODE_WRITE, &index_file ) != CG_OK )
    {
        cg_error_exit();
    }
    string baseName = "base";
    int celldim = 3;
    int phydim = 3;
    index_base = -1;
    cg_base_write( index_file, baseName.c_str(), celldim, phydim, & index_base );
    //put DataClass under Base
    cg_goto(index_file,index_base,"end");
    cg_dataclass_write( CGNS_ENUMV( NormalizedByUnknownDimensional ) );

    //put ReferenceState under Base
    cg_state_write("ReferenceQuantities");
    //Go to ReferenceState node, write Mach array and its dataclass
    cg_goto(index_file,index_base,"ReferenceState_t",1,"end");
    cg_array_write("Mach",CGNS_ENUMV(RealDouble),1,&nuse,&xmach);
    cg_goto(index_file,index_base,"ReferenceState_t",1,"DataArray_t",1,"end");
    cg_dataclass_write(CGNS_ENUMV(NondimensionalParameter));
    //Go to ReferenceState node, write Reynolds array and its dataclass
    cg_goto(index_file,index_base,"ReferenceState_t",1,"end");
    cg_array_write("Reynolds",CGNS_ENUMV(RealDouble),1,&nuse,&reue);
    cg_goto(index_file,index_base,"ReferenceState_t",1,"DataArray_t",2,"end");
    cg_dataclass_write(CGNS_ENUMV(NondimensionalParameter));
    //Go to ReferenceState node to write reference quantities
    cg_goto(index_file,index_base,"ReferenceState_t",1,"end");
    // First, write reference quantities that make up Mach and Reynolds:
    //Mach_Velocity
    cg_array_write("Mach_Velocity",CGNS_ENUMV(RealDouble),1,&nuse,&xmv);
    //Mach_VelocitySound
    cg_array_write("Mach_VelocitySound",CGNS_ENUMV(RealDouble),1,&nuse,&xmc);
    //Reynolds_Velocity
    cg_array_write("Reynolds_Velocity",CGNS_ENUMV(RealDouble),1,&nuse,&rev);
    //Reynolds_Length
    cg_array_write("Reynolds_Length",CGNS_ENUMV(RealDouble),1,&nuse,&rel);
    //Reynolds_ViscosityKinematic
    cg_array_write("Reynolds_ViscosityKinematic",CGNS_ENUMV(RealDouble),1,&nuse,&renu);

    /* Next, write flow field reference quantities: */
    /* Density */
    cg_array_write("Density",CGNS_ENUMV(RealDouble),1,&nuse,&rho0);
    /* Pressure */
    cg_array_write("Pressure",CGNS_ENUMV(RealDouble),1,&nuse,&p0);
    /* VelocitySound */
    cg_array_write("VelocitySound",CGNS_ENUMV(RealDouble),1,&nuse,&c0);
    /* ViscosityMolecular */
    cg_array_write("ViscosityMolecular",CGNS_ENUMV(RealDouble),1,&nuse,&vm0);
    /* LengthReference */
    cg_array_write("LengthReference",CGNS_ENUMV(RealDouble),1,&nuse,&xlength0);
    /* VelocityX */
    cg_array_write("VelocityX",CGNS_ENUMV(RealDouble),1,&nuse,&vx);
    /* VelocityY */
    //cg_array_write("VelocityY",CGNS_ENUMV(RealDouble),1,&nuse,&vy);
    /* VelocityZ */
    cg_array_write("VelocityZ",CGNS_ENUMV(RealDouble),1,&nuse,&vz);
    /* close CGNS file */
    cg_close(index_file);
    cout << "\n";
    cout <<"Successfully wrote nondimensional info to file grid_c.cgns\n";

}

void CgnsTest::WriteDescriptor()
{
    string gridName = "/grid/grid_c.cgns";
    string prjFileName = ONEFLOW::GetPrjFileName( gridName );
    cout << " prjFileName = " << prjFileName << "\n";

    this->fileName = prjFileName;

    int index_file,index_base;
    char textstring[74];

    cout << "Program write_descriptor\n";

    //if ( cg_open( this->fileName.c_str(), CG_MODE_MODIFY, &index_file ) )
    if ( cg_open( this->fileName.c_str(), CG_MODE_WRITE, &index_file ) )
    {
        cg_error_exit();
    }
    index_base = 1;
    cg_goto( index_file, index_base, "end" );

    /* write descriptor node (user can give any name) */
    string a = "Supersonic vehicle with landing gear\n";
    string b = "M=4.6, Re=6 million";
    string c = a + b;

    //cg_descriptor_write("Information",textstring);
    cg_descriptor_write("Information",c.c_str());
    cg_close(index_file);

    cout << "Successfully wrote descriptor node to file " << this->fileName << "\n";
}


EndNameSpace