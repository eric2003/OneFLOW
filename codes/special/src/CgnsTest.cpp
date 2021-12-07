/*---------------------------------------------------------------------------*\
    OneFLOW - LargeScale Multiphysics Scientific Simulation Environment
    Copyright (C) 2017-2021 He Xin and the OneFLOW contributors.
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
#include "CgnsTestTmp.h"
#include "CgnsFile.h"
#include "CgnsBase.h"
#include "CgnsFactory.h"
#include "Prj.h"
#include "StrUtil.h"
#include "CgnsZone.h"
#include "CgnsZbc.h"
#include "CgnsZbcBoco.h"
#include "CgnsBcBoco.h"
#include <cstring>
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
    this->mytest_read();
    //this->mytest_write();
    //this->read_grid_unst();
    //this->read_bcpnts_unst();
    //this->write_bcpnts_unst();
    //this->write_grid_unst();
    //this->Init();
    //this->WriteTest();
    //this->WriteSimpleMultiBaseTest();
    //this->ReadSimpleMultiBaseTest();
    //this->WriteEmptyCgnsFile();
    //this->WriteBc();
    //this->ReadBc();
    //this->ReadEmptyCgnsFile();
    //this->WriteDescriptor();
    //this->ReadDescriptor();
    //this->Test();
    //this->TestCgnsLink();
    //this->WriteArray();
    //this->ReadArray();
    //this->WriteReferenceState();
    //this->ReadReferenceState();
    //this->WriteConvergence();
    //this->ReadConvergence();
    //this->WriteFlowEqn();
    //this->ReadFlowEqn();
    
}

void CgnsTest::Test()
{
   this->TestCgnsLink();
}

void CgnsTest::SetDefaultGridName()
{
    std::string gridName = "/grid/oneflow.cgns";
    std::string prjFileName = ONEFLOW::GetPrjFileName( gridName );
    std::cout << " CGNS File Name = " << prjFileName << "\n";

    this->fileName = prjFileName;
}

void CgnsTest::WriteSimpleMultiBaseTest()
{
    CgnsFile * cgnsFile = new CgnsFile( "smiplebase.cgns", CG_MODE_WRITE );
    cgnsFile->WriteBase( "OneFLOW1" );
    cgnsFile->WriteBase( "OneFLOW 2" );
    cgnsFile->WriteBase( "CGNS base 3" );
    cgnsFile->WriteBase( "Fluid" );
    cgnsFile->WriteBase( "CAE library" );
    delete cgnsFile;
}

void CgnsTest::ReadSimpleMultiBaseTest()
{
    CgnsFile * cgnsFile = new CgnsFile( "smiplebase.cgns", CG_MODE_READ );
    cgnsFile->ReadBases();
    delete cgnsFile;
}

void CgnsTest::WriteDescriptor()
{
    CgnsFile * cgnsFile = new CgnsFile( "descript.cgns", CG_MODE_WRITE );
    cgnsFile->WriteBaseDescriptor();
    delete cgnsFile;
}

void CgnsTest::ReadDescriptor()
{
    CgnsFile * cgnsFile = new CgnsFile( "descript.cgns", CG_MODE_READ );
    cgnsFile->ReadBaseDescriptor();
    delete cgnsFile;
}

void CgnsTest::WriteEmptyCgnsFile()
{
    //cout << " CgnsTest::WriteEmptyCgnsFile() " << "\n";
    CgnsFile * cgnsFile = new CgnsFile( "empty.cgns", CG_MODE_WRITE );
    //cout << " 111 " << "\n";
    delete cgnsFile;
    //cout << " 222 " << "\n";
}

void CgnsTest::ReadEmptyCgnsFile()
{
    CgnsFile * cgnsFile = new CgnsFile( "empty.cgns", CG_MODE_READ );
    delete cgnsFile;
}

void CgnsTest::WriteDouble( const std::string & varName, const double & varValue )
{
    int nDim = 1;
    cgsize_t ndims[ 1 ] = { 1 };
    cg_array_write( varName.c_str(), CGNS_ENUMV(RealDouble), nDim ,ndims, &varValue );
}

void CgnsTest::SetISize( cgsize_t * isize )
{
    int nijk = 5;
    for ( int n = 0; n < 3; n ++ )
    {
        isize[ n     ] = nijk;
        isize[ n + 3 ] = nijk - 1;
        isize[ n + 6 ] = 0;
    }
}


void CgnsTest::TestCgnsLink()
{
    std::string fname    = "zones.cgns";
    std::string linkname = "zones_link.cgns";

    cgsize_t isize[ 9 ];
    this->SetISize( isize );
    int nZones = 5;

    CgnsFile * fileZone = new CgnsFile( fname, CG_MODE_WRITE );
    CgnsBase * cgnsBase = fileZone->WriteBase( "Base" );

    for ( int iZone = 0; iZone < nZones; ++ iZone )
    {
        std::string name = AddString( "Zone", iZone + 1 );
        cgnsBase->WriteZoneInfo( name, CGNS_ENUMV( Structured ), isize );
    }
    delete fileZone;

    CgnsFile * fileZoneM = new CgnsFile( fname, CG_MODE_MODIFY );
    CgnsBase * cgnsBaseM = fileZoneM->WriteBase( "Base" );

    for ( int iZone = 0; iZone < nZones; ++ iZone )
    {
        std::string name = AddString( "Zone", iZone + 1 );
        cgnsBaseM->WriteZoneInfo( name, CGNS_ENUMV( Structured ), isize );
    }

    delete fileZoneM;

    CgnsFile * fileLink = new CgnsFile( linkname, CG_MODE_WRITE );
    CgnsBase * cgnsBaseLink = fileLink->WriteBase( "Base2" );
    cgnsBaseLink->GoToBase();

    for ( int iZone = 0; iZone < nZones; ++ iZone )
    {
        std::string name     = AddString( "Link to Zone", iZone + 1 );
        std::string linkpath = AddString( "/Base/Zone", iZone + 1 );

        cg_link_write( name.c_str(), fname.c_str(), linkpath.c_str() );
    }

    delete fileLink;
}

void CgnsTest::GetArray( std::vector< std::vector< float > > & myfloat2d )
{
    std::vector< float > a1, a2, a3;
    a1.push_back( 1 );
    a1.push_back( 2 );
    a1.push_back( 3 );

    a2.push_back( 10 );
    a2.push_back( 20 );
    a2.push_back( 30 );
    a2.push_back( 40 );

    a3.push_back( 100 );
    a3.push_back( 200 );
    a3.push_back( 300 );
    a3.push_back( 400 );
    a3.push_back( 500 );

    myfloat2d.push_back( a1 );
    myfloat2d.push_back( a2 );
    myfloat2d.push_back( a3 );
}

void CgnsTest::WriteArray()
{
    std::vector< std::vector< float > > myarray;
    this->GetArray( myarray );

    CgnsFile * cgnsFile = new CgnsFile( "array.cgns", CG_MODE_WRITE );
    CgnsBase * cgnsBase = cgnsFile->WriteBase( "BaseXXX" );
    this->WriteArray( cgnsFile, cgnsBase );
    cgnsBase = cgnsFile->WriteBase( "BaseYYY" );
    this->WriteArray( cgnsFile, cgnsBase );

    delete cgnsFile;
}

void CgnsTest::WriteArray( CgnsFile * cgnsFile, CgnsBase * cgnsBase )
{
    std::vector< std::vector< float > > myarray;
    this->GetArray( myarray );

    cgnsBase->GoToBase();
    cg_user_data_write( "DataYYY" );
    cgnsBase->GoToNode( "UserDefinedData_t", 1 );

    for ( int i = 0; i < myarray.size(); ++ i )
    {
        std::string name = AddString( "MyArray", i + 1 );
        cgsize_t arraysize = myarray[ i ].size();
        cg_array_write( name.c_str(), CGNS_ENUMV( RealSingle ), 1, &arraysize, &myarray[ i ][ 0 ] );
        cgnsFile->GoPath( name );
        cgnsFile->GoPath( ".." );
    }

    cgnsBase->GoToBase();
    cg_user_data_write( "DataZZZ" );
    cgnsBase->GoToNode( "UserDefinedData_t", 2 );

    for ( int i = 0; i < myarray.size(); ++ i )
    {
        std::string name = AddString( "MyArray", i + 1 );
        cgsize_t arraysize = myarray[ i ].size();
        cg_array_write( name.c_str(), CGNS_ENUMV( RealSingle ), 1, &arraysize, &myarray[ i ][ 0 ] );
        cgnsFile->GoPath( name );
        cgnsFile->GoPath( ".." );
    }
}

void CgnsTest::ReadArray()
{
    CgnsFile * cgnsFile = new CgnsFile( "array.cgns", CG_MODE_READ );
    cgnsFile->ReadArray();
    delete cgnsFile;
}

void CgnsTest::WriteReferenceState()
{
    //define nondimensional parameters
    double xmach    = 4.6;
    double reue     = 6000000.;
    double xmv      = xmach;
    double xmc      = 1.;
    double rev      = xmach;
    double rel      = 1.;
    double renu     = xmach / reue;
    double rho0     = 1.;
    double gamma    = 1.4;
    double p0       = 1./gamma;
    double c0       = 1.;
    double vm0      = xmach/reue;
    double xlength0 = 1.;
    double vx       = xmach;
    double vy       = 0.0;
    double vz       = 0.0;

    CgnsFile * cgnsFile = new CgnsFile( "refstate.cgns", CG_MODE_WRITE );
    CgnsBase * cgnsBase1 = cgnsFile->WriteBase( "Base1" );

    cgnsBase1->GoToBase();
    cg_dataclass_write(CGNS_ENUMV(NormalizedByUnknownDimensional));
    cg_state_write("ReferenceQuantities");
    cgnsBase1->GoToNode( "ReferenceState_t", 1 );

    WriteDouble("Mach", xmach );
    WriteDouble("Reynolds", reue );

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

    CgnsBase * cgnsBase2 = cgnsFile->WriteBase( "Base2" );
    cgnsBase2->GoToBase();
    cg_dataclass_write(CGNS_ENUMV(NormalizedByUnknownDimensional));
    cg_state_write("Test1");
    cgnsBase1->GoToNode( "ReferenceState_t", 1 );
    WriteDouble("Mach", xmach );
    WriteDouble("Reynolds", reue );

    CgnsBase * cgnsBase3 = cgnsFile->WriteBase( "Base3" );
    cgnsBase3->GoToBase();
    cg_state_write("Test2");

    delete cgnsFile; 
}

void CgnsTest::ReadReferenceState()
{
    CgnsFile * cgnsFile = new CgnsFile( "refstate.cgns", CG_MODE_READ );
    cgnsFile->ReadReferenceState();
    delete cgnsFile; 
}

void CgnsTest::WriteConvergence()
{
    CgnsFile * cgnsFile = new CgnsFile( "convergence.cgns", CG_MODE_WRITE );
    CgnsBase * cgnsBase = cgnsFile->WriteBase( "Base" );
    cgnsBase->GoToBase();
    const int nIterations = 20;
    std::vector< double > cl( nIterations ), dl( 2 * nIterations );
    /* create history array simple example: */
    for ( int n = 0; n < nIterations; ++ n )
    {
        cl[ n ] = static_cast< float >( n + 1.0 );
    }

    for ( int n = 0; n < 2 * nIterations; ++ n )
    {
        dl[ n ] = - static_cast< float >( n + 1.0 );
    }

    /* create history node (SIDS names it GlobalConvergenceHistory at base level) */
    //cg_convergence_write( nIterations, "" );
    cg_convergence_write( nIterations, "haha" );
    /* go to new history node */
    cgnsBase->GoToNode( "ConvergenceHistory_t", 1 );
    /* write lift coefficient array (user must use SIDS-standard name here) */

    cgsize_t nuse = nIterations;
    cgsize_t muse = 2 * nIterations;
    cg_array_write("CoefLift",CGNS_ENUMV(RealDouble), 1, &nuse, &cl[ 0 ] );
    cg_array_write("DoefLift",CGNS_ENUMV(RealDouble), 1, &muse, &dl[ 0 ] );
    delete cgnsFile; 
}

void CgnsTest::ReadConvergence()
{
    CgnsFile * cgnsFile = new CgnsFile( "convergence.cgns", CG_MODE_READ );
    cgnsFile->ReadConvergence();
    delete cgnsFile; 
}


void CgnsTest::WriteFlowEqn()
{
    float gamma   = 1.4;
    float prandtl = 0.90;

    int idata[6];
    idata[0]=0;
    idata[1]=1;
    idata[2]=0;
    idata[3]=0;
    idata[4]=0;
    idata[5]=0;

    CgnsFile * cgnsFile = new CgnsFile( "floweqn.cgns", CG_MODE_WRITE );
    CgnsBase * cgnsBase = cgnsFile->WriteBase( "Base1" );
    CgnsZone * cgnsZone = cgnsBase->WriteZone( "Zone1" );
    cgnsZone->GoToZone();

    //Create 'FlowEquationSet' node under 'Zone_t'
    // equation dimension = 3
    int ieq_dim = 3;
    cg_equationset_write( ieq_dim );

    //Create 'GoverningEquations' node under 'FlowEquationSet'
    cgnsZone->GoToNode( "FlowEquationSet_t", 1 );
    cg_governing_write( CGNS_ENUMV(NSTurbulent) );

    //Create 'DiffusionModel' node under 'GoverningEquations'
    cgnsZone->GoToNode( "FlowEquationSet_t", 1, "GoverningEquations_t",1 );
    cg_diffusion_write(idata);

    cgsize_t nuse = 1;
    //Create 'GasModel' under 'FlowEquationSet'
    cgnsZone->GoToNode( "FlowEquationSet_t", 1 );
    cg_model_write("GasModel_t",CGNS_ENUMV(Ideal));

    // Create 'SpecificHeatRatio' under GasModel
    cgnsZone->GoToNode( "FlowEquationSet_t", 1, "GasModel_t",1 );
    cg_array_write("SpecificHeatRatio",CGNS_ENUMV(RealSingle), 1, &nuse, &gamma);
    // Create 'DataClass' under 'SpecificHeatRatio'

    cgnsZone->GoToNode( "FlowEquationSet_t", 1, "GasModel_t",1, "DataArray_t",1 );
    cg_dataclass_write(CGNS_ENUMV(NondimensionalParameter));

    //Create 'TurbulenceClosure' under 'FlowEquationSet'
    cgnsZone->GoToNode( "FlowEquationSet_t", 1 );
    cg_model_write("TurbulenceClosure_t",CGNS_ENUMV(EddyViscosity));

    //Create 'PrandtlTurbulent' under 'TurbulenceClosure'
    cgnsZone->GoToNode( "FlowEquationSet_t", 1, "TurbulenceClosure_t",1 );
    cg_array_write("PrandtlTurbulent",CGNS_ENUMV(RealSingle),1,&nuse,&prandtl);
    //Create 'DataClass' under 'PrandtlTurbulent'
    cgnsZone->GoToNode( "FlowEquationSet_t", 1, "TurbulenceClosure_t", 1, "DataArray_t", 1 );
    cg_dataclass_write(CGNS_ENUMV(NondimensionalParameter));

    //Create 'TurbulenceModel' under 'FlowEquationSet'
    cgnsZone->GoToNode( "FlowEquationSet_t", 1 );
    cg_model_write("TurbulenceModel_t",CGNS_ENUMV(OneEquation_SpalartAllmaras));
    delete cgnsFile;
}

void CgnsTest::ReadFlowEqn()
{
    CgnsFile * cgnsFile = new CgnsFile( "floweqn.cgns", CG_MODE_READ );
    cgnsFile->ReadFlowEqn();
    delete cgnsFile;
}

void CgnsTest::WriteTest()
{
    //init_data();
    //int cgfile;

    //cg_open( this->fileName.c_str(), CG_MODE_WRITE, &cgfile );
    //SetCgFile( cgfile );

    //write_structured();
    //write_unstructured();
    //write_mixed();
    //write_mismatched();

    //cg_close( GetCgFile() );
}

int CgnsTest::read_bcpnts_unst()
{
    int index_file,index_base,index_zone,nbocos,ib;
    int normalindex[3],ndataset;
    int i,normallist;
    char boconame[33];
    CGNS_ENUMT(BCType_t) ibocotype;
    CGNS_ENUMT(PointSetType_t) iptset;
    CGNS_ENUMT(DataType_t) normaldatatype;
    CGNS_ENUMT(GridLocation_t) igr;
    const int maxpnts = 960;
    cgsize_t npts,normallistflag;
    cgsize_t ipnts[maxpnts];

    /* READ BOUNDARY CONDITIONS FROM EXISTING CGNS FILE */
    /* open CGNS file for read-only */
    if (cg_open("grid_c.cgns",CG_MODE_READ,&index_file)) cg_error_exit();
    /* we know there is only one base (real working code would check!) */
    index_base=1;
    /* we know there is only one zone (real working code would check!) */
    index_zone=1;
    /* find out number of BCs that exist under this zone */
    cg_nbocos(index_file,index_base,index_zone,&nbocos);
    /* do loop over the total number of BCs */
    for (ib=1; ib <= nbocos; ib++)
    {
        /* find out what BC grid location is (expecting FaceCenter) */
        cg_goto(index_file,index_base,"Zone_t",1,"ZoneBC_t",1,"BC_t",ib,"end");
        cg_gridlocation_read(&igr);
        if (igr == CGNS_ENUMV(FaceCenter))
        {
            printf("\nGridLocation=FaceCenter means BC data refers to elements, not nodes\n");
        }
        /* get BC info */
        cg_boco_info(index_file,index_base,index_zone,ib,boconame,&ibocotype,
            &iptset,&npts,normalindex,&normallistflag,&normaldatatype,&ndataset);
        if (iptset != CGNS_ENUMV(PointList))
        {
            printf("\nError.  For this program, BCs must be std::set up as PointList type %s\n",
                PointSetTypeName[iptset]);
            return 1;
        }
        printf("\nBC number: %i\n",ib);
        printf("   name= %s\n",boconame);
        printf("   type= %s\n",BCTypeName[ibocotype]);
        printf("   no of elements= %i\n",(int)npts);
        if (npts > maxpnts)
        {
            printf("\nError.  Must increase maxpnts to at least %i\n",(int)npts);
            return 1;
        }
        /* read point list in here (nothing done with them in this program) */
        cg_boco_read(index_file,index_base,index_zone,ib,ipnts,&normallist);
        printf("      (these elements read here, but only some printed out:)\n");
        for (i=0; i < 10; i++)
        {
            printf("    ipnts[%i]=%i\n",i,(int)ipnts[i]);
        }
    }
    /* close CGNS file */
    cg_close(index_file);
    printf("\nSuccessfully read BCs (PointList format) from file grid_c.cgns\n");
    return 0;
}

int CgnsTest::write_bcpnts_unst()
{
    int index_file,index_base,index_zone;
    int nelem_start,nelem_end,icount,index_bc,ibc,n;
    const int maxcount = 960;
    cgsize_t ipnts[maxcount],icounts;

    printf("\nProgram write_bcpnts_unst\n");

    /* WRITE BOUNDARY CONDITIONS TO EXISTING CGNS FILE */
    /* open CGNS file for modify */
    if (cg_open("grid_c.cgns",CG_MODE_MODIFY,&index_file)) cg_error_exit();
    /* we know there is only one base (real working code would check!) */
    index_base=1;
    /* we know there is only one zone (real working code would check!) */
    index_zone=1;
    /* we know that for the unstructured zone, the following face elements */
    /* have been defined as inflow (real working code would check!): */
    nelem_start=2561;
    nelem_end=2688;
    icount=0;
    for (n=nelem_start; n <= nelem_end; n++)
    {
        ipnts[icount]=n;
        icount=icount+1;
    }
    if (icount > maxcount)
    {
        printf("\nError. Need to increase maxcount to at least %i\n",icount);
        return 1;
    }
    /* write boundary conditions for ilo face */
    icounts=icount;
    cg_boco_write(index_file,index_base,index_zone,"Ilo",CGNS_ENUMV(BCTunnelInflow),
        CGNS_ENUMV(PointList),icounts,ipnts,&index_bc);
    /* we know that for the unstructured zone, the following face elements */
    /* have been defined as outflow (real working code would check!): */
    nelem_start=2689;
    nelem_end=2816;
    icount=0;
    for (n=nelem_start; n <= nelem_end; n++)
    {
        ipnts[icount]=n;
        icount=icount+1;
    }
    if (icount > maxcount)
    {
        printf("\nError. Need to increase maxcount to at least %i\n",icount);
        return 1;
    }
    /* write boundary conditions for ihi face */
    icounts=icount;
    cg_boco_write(index_file,index_base,index_zone,"Ihi",CGNS_ENUMV(BCExtrapolate),
        CGNS_ENUMV(PointList),icounts,ipnts,&index_bc);
    /* we know that for the unstructured zone, the following face elements */
    /* have been defined as walls (real working code would check!): */
    nelem_start=2817;
    nelem_end=3776;
    icount=0;
    for (n=nelem_start; n <= nelem_end; n++)
    {
        ipnts[icount]=n;
        icount=icount+1;
    }
    if (icount > maxcount)
    {
        printf("\nError. Need to increase maxcount to at least %i\n",icount);
        return 1;
    }
    /* write boundary conditions for wall faces */
    icounts=icount;
    cg_boco_write(index_file,index_base,index_zone,"Walls",CGNS_ENUMV(BCWallInviscid),
        CGNS_ENUMV(PointList),icounts,ipnts,&index_bc);

    /* the above are all face-center locations for the BCs - must indicate this, */
    /* otherwise Vertices will be assumed! */
    for (ibc=1; ibc <= index_bc; ibc++)
    {
        /*    (the following call positions you in BC_t - it assumes there */
        /*    is only one Zone_t and one ZoneBC_t - real working code would check!) */
        cg_goto(index_file,index_base,"Zone_t",1,"ZoneBC_t",1,"BC_t",ibc,"end");
        cg_gridlocation_write(CGNS_ENUMV(FaceCenter));
    }
    /* close CGNS file */
    cg_close(index_file);
    printf("\nSuccessfully added FaceCenter BCs (PointList) to unstructured grid file grid_c.cgns\n");
    return 0;
}

int CgnsTest::write_grid_unst()
{
    const int maxelemi = 20 * 16 * 8;
    const int maxelemj = 1216;

    double x[21*17*9],y[21*17*9],z[21*17*9];
    cgsize_t isize[3][1],ielem[maxelemi][8],jelem[maxelemj][4];
    cgsize_t nelem_start,nelem_end;
    int ni,nj,nk,iset,i,j,k,index_file,icelldim,iphysdim;
    int index_base,index_zone,index_coord,ielem_no;
    int ifirstnode,nbdyelem,index_section;
    char basename[33],zonename[33];

    printf("\nProgram write_grid_unst\n");

    /* create gridpoints for simple example: */
    ni=21;
    nj=17;
    nk=9;
    iset=0;
    for (k=1; k <= nk; k++)
    {
        for (j=1; j <=nj; j++)
        {
            for (i=1; i <= ni; i++)
            {
                x[iset]=(float)i-1.;
                y[iset]=(float)j-1.;
                z[iset]=(float)k-1.;
                iset=iset+1;
            }
        }
    }
    printf("\ncreated simple 3-D grid points\n");

    /* WRITE X, Y, Z GRID POINTS TO CGNS FILE */
    /* open CGNS file for write */
    if (cg_open("grid_c.cgns",CG_MODE_WRITE,&index_file)) cg_error_exit();
    /* create base (user can give any name) */
    strcpy(basename,"Base");
    icelldim=3;
    iphysdim=3;
    cg_base_write(index_file,basename,icelldim,iphysdim,&index_base);
    /* define zone name (user can give any name) */
    strcpy(zonename,"Zone  1");
    /* vertex size */
    isize[0][0]=ni*nj*nk;
    /* cell size */
    isize[1][0]=(ni-1)*(nj-1)*(nk-1);
    /* boundary vertex size (zero if elements not sorted) */
    isize[2][0]=0;
    /* create zone */
    cg_zone_write(index_file,index_base,zonename,isize[0],CGNS_ENUMV(Unstructured),&index_zone);
    /* write grid coordinates (user must use SIDS-standard names here) */
    cg_coord_write(index_file,index_base,index_zone,CGNS_ENUMV(RealDouble),"CoordinateX",
        x,&index_coord);
    cg_coord_write(index_file,index_base,index_zone,CGNS_ENUMV(RealDouble),"CoordinateY",
        y,&index_coord);
    cg_coord_write(index_file,index_base,index_zone,CGNS_ENUMV(RealDouble),"CoordinateZ",
        z,&index_coord);
    /* std::set element connectivity: */
    /* ---------------------------------------------------------- */
    /* do all the HEXA_8 elements (this part is mandatory): */
    /* maintain SIDS-standard ordering */
    ielem_no=0;
    /* index no of first element */
    nelem_start=1;
    for (k=1; k < nk; k++)
    {
        for (j=1; j < nj; j++)
        {
            for (i=1; i < ni; i++)
            {
                /*
                in this example, due to the order in the node numbering, the
                hexahedral elements can be reconstructed using the following
                relationships:
                */
                ifirstnode=i+(j-1)*ni+(k-1)*ni*nj;
                ielem[ielem_no][0]=ifirstnode;
                ielem[ielem_no][1]=ifirstnode+1;
                ielem[ielem_no][2]=ifirstnode+1+ni;
                ielem[ielem_no][3]=ifirstnode+ni;
                ielem[ielem_no][4]=ifirstnode+ni*nj;
                ielem[ielem_no][5]=ifirstnode+ni*nj+1;
                ielem[ielem_no][6]=ifirstnode+ni*nj+1+ni;
                ielem[ielem_no][7]=ifirstnode+ni*nj+ni;
                ielem_no=ielem_no+1;
            }
        }
    }
    /* index no of last element (=2560) */
    nelem_end=ielem_no;
    if (nelem_end > maxelemi)
    {
        printf("\nError, must increase maxelemi to at least %d\n",nelem_end);
        return 1;
    }
    /* unsorted boundary elements */
    nbdyelem=0;
    /* write CGNS_ENUMV(HEXA_8) element connectivity (user can give any name) */
    cg_section_write(index_file,index_base,index_zone,"Elem",CGNS_ENUMV(HEXA_8),nelem_start,
        nelem_end,nbdyelem,ielem[0],&index_section);
    /* ---------------------------------------------------------- */
    /*
    do boundary (QUAD) elements (this part is optional,
    but you must do it if you eventually want to define BCs
    at element faces rather than at nodes):
    maintain SIDS-standard ordering
    */
    /* INFLOW: */
    ielem_no=0;
    /* index no of first element */
    nelem_start=nelem_end+1;
    i=1;
    for (k=1; k < nk; k++)
    {
        for (j=1; j < nj; j++)
        {
            ifirstnode=i+(j-1)*ni+(k-1)*ni*nj;
            jelem[ielem_no][0]=ifirstnode;
            jelem[ielem_no][1]=ifirstnode+ni*nj;
            jelem[ielem_no][2]=ifirstnode+ni*nj+ni;
            jelem[ielem_no][3]=ifirstnode+ni;
            ielem_no=ielem_no+1;
        }
    }
    /* index no of last element */
    nelem_end=nelem_start+ielem_no-1;
    if (ielem_no > maxelemj)
    {
        printf("\nError, must increase maxelemj to at least %d\n",ielem_no);
        return 1;
    }
    /* write QUAD element connectivity for inflow face (user can give any name) */
    cg_section_write(index_file,index_base,index_zone,"InflowElem",CGNS_ENUMV(QUAD_4),nelem_start,
        nelem_end,nbdyelem,jelem[0],&index_section);
    /* OUTFLOW: */
    ielem_no=0;
    /* index no of first element */
    nelem_start=nelem_end+1;
    i=ni-1;
    for (k=1; k < nk; k++)
    {
        for (j=1; j < nj; j++)
        {
            ifirstnode=i+(j-1)*ni+(k-1)*ni*nj;
            jelem[ielem_no][0]=ifirstnode+1;
            jelem[ielem_no][1]=ifirstnode+1+ni;
            jelem[ielem_no][2]=ifirstnode+ni*nj+1+ni;
            jelem[ielem_no][3]=ifirstnode+ni*nj+1;
            ielem_no=ielem_no+1;
        }
    }
    /* index no of last element */
    nelem_end=nelem_start+ielem_no-1;
    if (ielem_no > maxelemj)
    {
        printf("\nError, must increase maxelemj to at least %d\n",ielem_no);
        return 1;
    }
    /* write QUAD element connectivity for outflow face (user can give any name) */
    cg_section_write(index_file,index_base,index_zone,"OutflowElem",CGNS_ENUMV(QUAD_4),nelem_start,
        nelem_end,nbdyelem,jelem[0],&index_section);
    /* SIDEWALLS: */
    ielem_no=0;
    /* index no of first element */
    nelem_start=nelem_end+1;
    j=1;
    for (k=1; k < nk; k++)
    {
        for (i=1; i < ni; i++)
        {
            ifirstnode=i+(j-1)*ni+(k-1)*ni*nj;
            jelem[ielem_no][0]=ifirstnode;
            jelem[ielem_no][1]=ifirstnode+ni*nj;
            jelem[ielem_no][2]=ifirstnode+ni*nj+1;
            jelem[ielem_no][3]=ifirstnode+1;
            ielem_no=ielem_no+1;
        }
    }
    j=nj-1;
    for (k=1; k < nk; k++)
    {
        for (i=1; i < ni; i++)
        {
            ifirstnode=i+(j-1)*ni+(k-1)*ni*nj;
            jelem[ielem_no][0]=ifirstnode+1+ni;
            jelem[ielem_no][1]=ifirstnode+ni;
            jelem[ielem_no][2]=ifirstnode+ni*nj+ni;
            jelem[ielem_no][3]=ifirstnode+ni*nj+1+ni;
            ielem_no=ielem_no+1;
        }
    }
    k=1;
    for (j=1; j < nj; j++)
    {
        for (i=1; i < ni; i++)
        {
            ifirstnode=i+(j-1)*ni+(k-1)*ni*nj;
            jelem[ielem_no][0]=ifirstnode;
            jelem[ielem_no][1]=ifirstnode+1;
            jelem[ielem_no][2]=ifirstnode+1+ni;
            jelem[ielem_no][3]=ifirstnode+ni;
            ielem_no=ielem_no+1;
        }
    }
    k=nk-1;
    for (j=1; j < nj; j++)
    {
        for (i=1; i < ni; i++)
        {
            ifirstnode=i+(j-1)*ni+(k-1)*ni*nj;
            jelem[ielem_no][0]=ifirstnode+ni*nj;
            jelem[ielem_no][1]=ifirstnode+ni*nj+ni;
            jelem[ielem_no][2]=ifirstnode+ni*nj+1+ni;
            jelem[ielem_no][3]=ifirstnode+ni*nj+1;
            ielem_no=ielem_no+1;
        }
    }
    /* index no of last element */
    nelem_end=nelem_start+ielem_no-1;
    if (ielem_no > maxelemj)
    {
        printf("\nError, must increase maxelemj to at least %d\n",ielem_no);
        return 1;
    }
    /* write QUAD element connectivity for sidewall face (user can give any name) */
    cg_section_write(index_file,index_base,index_zone,"SidewallElem",CGNS_ENUMV(QUAD_4),nelem_start,
        nelem_end,nbdyelem,jelem[0],&index_section);
    /* ---------------------------------------------------------- */
    /* close CGNS file */
    cg_close(index_file);
    printf("\nSuccessfully wrote unstructured grid to file grid_c.cgns\n");
    return 0;
}

int CgnsTest::read_grid_unst()
{
    float x[21*17*9],y[21*17*9],z[21*17*9];
    cgsize_t isize[3][1],ielem[20*16*8][8];
    int index_file,index_base,index_zone;
    cgsize_t irmin,irmax,istart,iend;
    int nsections,index_sect,nbndry,iparent_flag;
    cgsize_t iparentdata;
    char zonename[33],sectionname[33];
    CGNS_ENUMT(ElementType_t) itype;

    /* READ X, Y, Z GRID POINTS FROM CGNS FILE */
    /* open CGNS file for read-only */
    if (cg_open("grid_c.cgns",CG_MODE_READ,&index_file)) cg_error_exit();
    /* we know there is only one base (real working code would check!) */
    index_base=1;
    /* we know there is only one zone (real working code would check!) */
    index_zone=1;
    /* get zone size (and name - although not needed here) */
    cg_zone_read(index_file,index_base,index_zone,zonename,isize[0]);
    /* lower range index */
    irmin=1;
    /* upper range index of vertices */
    irmax=isize[0][0];
    /* read grid coordinates */
    cg_coord_read(index_file,index_base,index_zone,"CoordinateX",
        CGNS_ENUMV(RealSingle),&irmin,&irmax,x);
    cg_coord_read(index_file,index_base,index_zone,"CoordinateY",
        CGNS_ENUMV(RealSingle),&irmin,&irmax,y);
    cg_coord_read(index_file,index_base,index_zone,"CoordinateZ",
        CGNS_ENUMV(RealSingle),&irmin,&irmax,z);
    /* find out how many sections */
    cg_nsections(index_file,index_base,index_zone,&nsections);
    printf("\nnumber of sections=%i\n",nsections);
    /* read element connectivity */
    for (index_sect=1; index_sect <= nsections; index_sect++)
    {
        cg_section_read(index_file,index_base,index_zone,index_sect,sectionname,
            &itype,&istart,&iend,&nbndry,&iparent_flag);
        printf("\nReading section data...\n");
        printf("   section name=%s\n",sectionname);
        printf("   section type=%s\n",ElementTypeName[itype]);
        printf("   istart,iend=%i, %i\n",(int)istart,(int)iend);
        if (itype == CGNS_ENUMV(HEXA_8))
        {
            printf("   reading element data for this element\n");
            cg_elements_read(index_file,index_base,index_zone,index_sect,ielem[0], \
                &iparentdata);
        }
        else
        {
            printf("   not reading element data for this element\n");
        }
    }
    /* close CGNS file */
    cg_close(index_file);
    printf("\nSuccessfully read unstructured grid from file grid_c.cgns\n");
    printf("   for example, element 1 is made up of nodes: %i, %i, %i, %i, %i, %i, %i, %i\n",
        (int)ielem[0][0],(int)ielem[0][1],(int)ielem[0][2],(int)ielem[0][3],
        (int)ielem[0][4],(int)ielem[0][5],(int)ielem[0][6],(int)ielem[0][7]);
    printf("   x,y,z of node 357 are: %f, %f, %f\n",x[357],y[357],z[357]);
    printf("   x,y,z of node 1357 are: %f, %f, %f\n",x[1357],y[1357],z[1357]);
    return 0;
}

void CgnsTest::mytest_read()
{
    CgnsFile * cgnsFile = new CgnsFile( "mytest.cgns", CG_MODE_READ );
    int index_base = -1;
    int icelldim = -1;
    int iphysdim = -1;

    cgnsFile->ReadNumberOfBases();
    std::cout << " cgnsFile->nBases = " << cgnsFile->nBases << "\n";

    for ( int iBase = 0; iBase < cgnsFile->nBases; ++ iBase )
    {
        CgnsBase * cgnsBase = cgnsFile->CreateCgnsBase();
        std::cout << " cgnsBase->baseId = " << cgnsBase->baseId << "\n";
        cgnsBase->ReadCgnsBaseBasicInfo();
        cgnsBase->ReadNumberOfCgnsZones();
        std::cout << " cgnsBase->nZones = " << cgnsBase->nZones << "\n";
        for ( int iZone = 0; iZone < cgnsBase->nZones; ++ iZone )
        {
            CgnsZone * cgnsZone = cgnsBase->CreateCgnsZone();
            std::cout << " cgnsZone->zId = " << cgnsZone->zId << "\n";

            cgnsZone->ReadCgnsZoneAttribute();
            cgnsZone->ReadCgnsGridBoundary();
        }
    }

    delete cgnsFile;
}

void CgnsTest::mytest_write()
{
    CgnsFile * cgnsFile = new CgnsFile( "mytest.cgns", CG_MODE_WRITE );
    int icelldim = 3;
    int iphysdim = 3;
    CgnsBase * cgnsBase = cgnsFile->WriteBase( "Base", icelldim, iphysdim );

    cgsize_t isize[ 3 ][ 1 ];

    int nNodes = 1;
    int nCells = 1;

    /* vertex size */
    isize[0][0] = nNodes;
    /* cell size */
    isize[1][0] = nCells;
    /* boundary vertex size (zero if elements not sorted) */
    isize[2][0] = 0;

    cgsize_t ipnts[ 1 ],icounts;
    ipnts[ 0 ] = 0;
    icounts = 1;

    std::string zoneName = "Zone1";
    CgnsZone * cgnsZone = cgnsBase->WriteZoneInfo( zoneName, CGNS_ENUMV(Unstructured), isize[ 0 ] );

    CgnsZbcBoco * cgnsZbcBoco = cgnsZone->cgnsZbc->cgnsZbcBoco;
    CgnsBcBoco * cgnsBcBoco = 0;
    cgnsBcBoco = cgnsZbcBoco->WriteCgnsBoco( "Bc1", CGNS_ENUMV(BCTunnelInflow), CGNS_ENUMV(PointList), icounts, ipnts );
    cgnsBcBoco->WriteGridLocation( CGNS_ENUMV(FaceCenter) );
    cgnsBcBoco = cgnsZbcBoco->WriteCgnsBoco( "Bc2", CGNS_ENUMV(BCTunnelInflow), CGNS_ENUMV(PointList), icounts, ipnts );
    cgnsBcBoco->WriteGridLocation( CGNS_ENUMV(Vertex) );
    cgnsBcBoco = cgnsZbcBoco->WriteCgnsBoco( "Bc3", CGNS_ENUMV(BCTunnelInflow), CGNS_ENUMV(PointList), icounts, ipnts );
    cgnsBcBoco->WriteGridLocation( CGNS_ENUMV(CellCenter) );

    zoneName = "Zone2";
    cgnsZone = cgnsBase->WriteZoneInfo( zoneName, CGNS_ENUMV(Unstructured), isize[ 0 ] );

    cgnsZbcBoco = cgnsZone->cgnsZbc->cgnsZbcBoco;
    cgnsBcBoco = cgnsZbcBoco->WriteCgnsBoco( "Bc_1", CGNS_ENUMV(BCTunnelInflow), CGNS_ENUMV(PointList), icounts, ipnts );
    cgnsBcBoco->WriteGridLocation( CGNS_ENUMV(Vertex) );
    cgnsBcBoco = cgnsZbcBoco->WriteCgnsBoco( "Bc_2", CGNS_ENUMV(BCTunnelInflow), CGNS_ENUMV(PointList), icounts, ipnts );
    cgnsBcBoco->WriteGridLocation( CGNS_ENUMV(Vertex) );
    cgnsBcBoco = cgnsZbcBoco->WriteCgnsBoco( "Bc_3", CGNS_ENUMV(BCTunnelInflow), CGNS_ENUMV(PointList), icounts, ipnts );
    cgnsBcBoco->WriteGridLocation( CGNS_ENUMV(CellCenter) );

    delete cgnsFile;
}


EndNameSpace
