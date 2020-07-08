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
#include "CgnsTestTmp.h"
#include "CgnsFile.h"
#include "CgnsBase.h"
#include "CgnsFactory.h"
#include "Prj.h"
#include "StrUtil.h"
#include "CgnsZone.h"
#include <iostream>
#include <iomanip>
#include <cmath>
using namespace std;

BeginNameSpace( ONEFLOW )

#define TWOPI 6.2831853

/* set STRUCTURED_FACES to write the structured zone bocos for the
mixed case as face range and face list instead of as point range/list
#define STRUCTURED_FACES
*/

/* set UNSTRUCTURED_1TO1 to write the unstructured case
zone connectivities as 1to1 instead of abutting1to1
#define UNSTRUCTURED_1TO1
*/

/* set ABUTTING1TO1_FACES to write the unstructured case
zone connectivites as abutting1to1 with faces(elements)
instead of points. This also writes the mixed case zone 1
to zone 2 connectivity with a face range
#define ABUTTING1TO1_FACES
*/

int CellDim = 3, PhyDim = 3;

int cgfile, cgbase, cgzone;
cgsize_t size[9];

#define NUM_SIDE 5
#define NODE_INDEX(I,J,K) ((I)+NUM_SIDE*(((J)-1)+NUM_SIDE*((K)-1)))
#define CELL_INDEX(I,J,K) ((I)+(NUM_SIDE-1)*(((J)-1)+(NUM_SIDE-1)*((K)-1)))

int num_coord;
float *xcoord, *ycoord, *zcoord;
int num_element, num_face;
cgsize_t *elements, *faces, *parent;

int max_sol;
float *solution;

int npts;
cgsize_t *pts, *d_pts;
float *interp;

char errmsg[128];

//void init_data();
//void write_structured(), write_unstructured();
//void write_mixed(), write_mismatched();

void SetCgFile( int cgfileIn )
{
    cgfile = cgfileIn;
}

int GetCgFile()
{
    return cgfile;
}

void error_exit (char *where)
{
    fprintf (stderr, "ERROR:%s:%s\n", where, cg_get_error());
    exit (1);
}

void init_data()
{
    int n, i, j, k, nn, nf, np;

    /* compute coordinates - make it twice as big for use with cylindrical */

    num_coord = NUM_SIDE * NUM_SIDE * NUM_SIDE;
    xcoord = (float *) malloc (6 * num_coord * sizeof(float));
    if (NULL == xcoord) {
        fprintf(stderr, "malloc failed for coordinates\n");
        exit(1);
    }
    ycoord = xcoord + 2 * num_coord;
    zcoord = ycoord + 2 * num_coord;
    for (n = 0, k = 0; k < NUM_SIDE; k++) {
        for (j = 0; j < NUM_SIDE; j++) {
            for (i = 0; i < NUM_SIDE; i++, n++) {
                xcoord[n] = (float)i;
                ycoord[n] = (float)j;
                zcoord[n] = (float)k;
            }
        }
    }

    /* create solution vector large enough for grid + 1 rind plane */

    max_sol = (NUM_SIDE + 2) * (NUM_SIDE + 2) * (NUM_SIDE + 2);
    solution = (float *) malloc (max_sol * sizeof(float));
    if (NULL == solution) {
        fprintf(stderr, "malloc failed for solution\n");
        exit(1);
    }
    for (n = 0; n < max_sol; n++)
        solution[n] = (float)(n + 1);

    /* compute elements */

    num_element = (NUM_SIDE - 1) * (NUM_SIDE - 1) * (NUM_SIDE - 1);
    elements = (cgsize_t *) malloc (8 * num_element * sizeof(cgsize_t));
    if (NULL == elements) {
        fprintf(stderr, "malloc failed for elements");
        exit(1);
    }
    for (n = 0, k = 1; k < NUM_SIDE; k++) {
        for (j = 1; j < NUM_SIDE; j++) {
            for (i = 1; i < NUM_SIDE; i++) {
                nn = NODE_INDEX(i, j, k);
                elements[n++] = nn;
                elements[n++] = nn + 1;
                elements[n++] = nn + 1 + NUM_SIDE;
                elements[n++] = nn + NUM_SIDE;
                nn += NUM_SIDE * NUM_SIDE;
                elements[n++] = nn;
                elements[n++] = nn + 1;
                elements[n++] = nn + 1 + NUM_SIDE;
                elements[n++] = nn + NUM_SIDE;
            }
        }
    }

    /* compute outside face elements */

    num_face = 6 * (NUM_SIDE - 1) * (NUM_SIDE - 1);
    faces = (cgsize_t *) malloc (4 * num_face * sizeof(cgsize_t));
    parent = (cgsize_t *) malloc (4 * num_face * sizeof(cgsize_t));
    if (NULL == faces || NULL == parent) {
        fprintf(stderr, "malloc failed for elements");
        exit(1);
    }
    for (n = 0; n < 4*num_face; n++)
        parent[n] = 0;
    nf = np = 0;
    n = 2 * num_face;
    i = 1;
    for (k = 1; k < NUM_SIDE; k++) {
        for (j = 1; j < NUM_SIDE; j++) {
            nn = NODE_INDEX(i, j, k);
            faces[nf++]  = nn;
            faces[nf++]  = nn + NUM_SIDE * NUM_SIDE;
            faces[nf++]  = nn + NUM_SIDE * (NUM_SIDE + 1);
            faces[nf++]  = nn + NUM_SIDE;
            parent[np]   = CELL_INDEX(i, j, k);
            parent[np+n] = 5;
            np++;
        }
    }
    i = NUM_SIDE;
    for (k = 1; k < NUM_SIDE; k++) {
        for (j = 1; j < NUM_SIDE; j++) {
            nn = NODE_INDEX(i, j, k);
            faces[nf++]  = nn;
            faces[nf++]  = nn + NUM_SIDE;
            faces[nf++]  = nn + NUM_SIDE * (NUM_SIDE + 1);
            faces[nf++]  = nn + NUM_SIDE * NUM_SIDE;
            parent[np]   = CELL_INDEX(i-1, j, k);
            parent[np+n] = 3;
            np++;
        }
    }
    j = 1;
    for (k = 1; k < NUM_SIDE; k++) {
        for (i = 1; i < NUM_SIDE; i++) {
            nn = NODE_INDEX(i, j, k);
            faces[nf++]  = nn;
            faces[nf++]  = nn + 1;
            faces[nf++]  = nn + 1 + NUM_SIDE * NUM_SIDE;
            faces[nf++]  = nn + NUM_SIDE * NUM_SIDE;
            parent[np]   = CELL_INDEX(i, j, k);
            parent[np+n] = 2;
            np++;
        }
    }
    j = NUM_SIDE;
    for (k = 1; k < NUM_SIDE; k++) {
        for (i = 1; i < NUM_SIDE; i++) {
            nn = NODE_INDEX(i, j, k);
            faces[nf++]  = nn;
            faces[nf++]  = nn + NUM_SIDE * NUM_SIDE;
            faces[nf++]  = nn + 1 + NUM_SIDE * NUM_SIDE;
            faces[nf++]  = nn + 1;
            parent[np]   = CELL_INDEX(i, j-1, k);
            parent[np+n] = 4;
            np++;
        }
    }
    k = 1;
    for (j = 1; j < NUM_SIDE; j++) {
        for (i = 1; i < NUM_SIDE; i++) {
            nn = NODE_INDEX(i, j, k);
            faces[nf++]  = nn;
            faces[nf++]  = nn + NUM_SIDE;
            faces[nf++]  = nn + NUM_SIDE + 1;
            faces[nf++]  = nn + 1;
            parent[np]   = CELL_INDEX(i, j, k);
            parent[np+n] = 1;
            np++;
        }
    }
    k = NUM_SIDE;
    for (j = 1; j < NUM_SIDE; j++) {
        for (i = 1; i < NUM_SIDE; i++) {
            nn = NODE_INDEX(i, j, k);
            faces[nf++]  = nn;
            faces[nf++]  = nn + 1;
            faces[nf++]  = nn + NUM_SIDE + 1;
            faces[nf++]  = nn + NUM_SIDE;
            parent[np]   = CELL_INDEX(i, j, k-1);
            parent[np+n] = 6;
            np++;
        }
    }

    /* connectivity points - make it big enough to hold 4 surfaces */

    npts = NUM_SIDE * NUM_SIDE;
    pts = (cgsize_t *) malloc (12 * npts * sizeof(cgsize_t));
    if (NULL == pts) {
        fprintf(stderr, "malloc failed for connectivity points");
        exit(1);
    }
    d_pts = pts + 6 * npts;

    /* create interpolate data array */

    interp = (float *) malloc (6 * npts * sizeof(float));
    if (NULL == interp) {
        fprintf(stderr, "malloc failed for interpolate array");
        exit(1);
    }
}

void write_reference ()
{
    int n, i, ierr;
    cgsize_t dim = 1;
    float exps[5];

    static struct {
        char *name;
        int dim;
        int exps[5];
        float val;
    } state[] = {
        {"Mach", 0, {0, 0, 0, 0, 0}, (float)0.2},
        {"VelocitySound", 1, {0, 1, -1, 0, 0}, (float)330.0},
        {"VelocityMagnitude", 1, {0, 1, -1, 0, 0}, (float)66.0},
        {"VelocityUnitVectorX", 0, {0, 0, 0, 0, 0}, (float)1.0},
        {"VelocityUnitVectorY", 0, {0, 0, 0, 0, 0}, (float)0.0},
        {"VelocityUnitVectorZ", 0, {0, 0, 0, 0, 0}, (float)0.0},
        {"Reynolds", 0, {0, 0, 0, 0, 0}, (float)3.0e6},
        {"Temperature", 1, {0, 0, 0, 1, 0}, (float)300.0},
        {"Pressure", 1, {1, -1, -2, 0, 0}, (float)1.0e5},
        {"LengthReference", 1, {0, 1, 0, 0, 0}, (float)10.0}
    };

    if (cg_goto(cgfile, cgbase, "end") ||
        cg_state_write("reference state quantities") ||
        cg_goto(cgfile, cgbase, "ReferenceState_t", 1, "end") ||
        cg_dataclass_write(CGNS_ENUMV(Dimensional)) ||
        cg_units_write(CGNS_ENUMV(Kilogram), CGNS_ENUMV(Meter),
            CGNS_ENUMV(Second), CGNS_ENUMV(Kelvin), CGNS_ENUMV(Radian)))
        error_exit("reference state");

    for (n = 0; n < 10; n++) {
        if (cg_goto(cgfile, cgbase, "ReferenceState_t", 1, "end") ||
            cg_array_write(state[n].name, CGNS_ENUMV(RealSingle), 1,
                &dim, &state[n].val) ||
            cg_goto(cgfile, cgbase, "ReferenceState_t", 1,
                "DataArray_t", n+1, "end")) {
            sprintf (errmsg, "reference state data %d", n+1);
            error_exit(errmsg);
        }
        if (state[n].dim) {
            for (i = 0; i < 5; i++)
                exps[i] = (float)state[n].exps[i];
            ierr = cg_exponents_write(CGNS_ENUMV(RealSingle), exps);
        }
        else
            ierr = cg_dataclass_write(CGNS_ENUMV(NondimensionalParameter));
        if (ierr) {
            sprintf (errmsg, "reference state data %d dataclass", n+1);
            error_exit(errmsg);
        }
    }
}

void write_equationset ()
{
    int n, diff[6];
    cgsize_t dim = 1;
    float g = (float)1.4;
    float R = (float)53.352;
    float ts = (float)110.6;
    float mu = (float)1.716e-5;
    float exp = (float)0.666;
    float pt = (float)0.9;
    float exps[5];

    for (n = 0; n < 6; n++)
        diff[n] = 0;
    diff[2] = 1;

    /* flow equation set */

    if (cg_goto(cgfile, cgbase, "end") ||
        cg_equationset_write (3) ||
        cg_goto(cgfile, cgbase, "FlowEquationSet_t", 1, "end") ||
        cg_governing_write(CGNS_ENUMV(NSTurbulent)) ||
        cg_model_write("GasModel_t", CGNS_ENUMV(Ideal)) ||
        cg_model_write("ViscosityModel_t", CGNS_ENUMV(SutherlandLaw)) ||
        cg_model_write("ThermalConductivityModel_t", CGNS_ENUMV(PowerLaw)) ||
        cg_model_write("TurbulenceClosure_t", CGNS_ENUMV(EddyViscosity)) ||
        cg_model_write("TurbulenceModel_t", CGNS_ENUMV(Algebraic_BaldwinLomax)))
        error_exit("flow equation set");

    if (cg_dataclass_write(CGNS_ENUMV(Dimensional)) ||
        cg_units_write(CGNS_ENUMV(Kilogram), CGNS_ENUMV(Meter),
            CGNS_ENUMV(Second), CGNS_ENUMV(Kelvin), CGNS_ENUMV(Radian)))
        error_exit("flow equation set dataclass");

    /* diffusion model */

    if (cg_goto(cgfile, cgbase, "FlowEquationSet_t", 1,
        "GoverningEquations_t", 1, "end") ||
        cg_diffusion_write(diff) ||
        cg_goto(cgfile, cgbase, "FlowEquationSet_t", 1,
            "TurbulenceModel_t", 1, "end") ||
        cg_diffusion_write(diff))
        error_exit("diffusion model");

    /* gas model */

    exps[0] = (float)0.0;
    exps[1] = (float)2.0;
    exps[2] = (float)-2.0;
    exps[3] = (float)-1.0;
    exps[4] = (float)0.0;
    if (cg_goto(cgfile, cgbase, "FlowEquationSet_t", 1,
        "GasModel_t", 1, "end") ||
        cg_dataclass_write(CGNS_ENUMV(DimensionlessConstant)) ||
        cg_array_write("SpecificHeatRatio", CGNS_ENUMV(RealSingle),
            1, &dim, &g) ||
        cg_array_write("IdealGasConstant", CGNS_ENUMV(RealSingle),
            1, &dim, &R) ||
        cg_goto(cgfile, cgbase, "FlowEquationSet_t", 1,
            "GasModel_t", 1, "DataArray_t", 2, "end") ||
        cg_dataclass_write(CGNS_ENUMV(Dimensional)) ||
        cg_units_write(CGNS_ENUMV(Slug), CGNS_ENUMV(Foot),
            CGNS_ENUMV(Second), CGNS_ENUMV(Rankine), CGNS_ENUMV(Radian)) ||
        cg_exponents_write(CGNS_ENUMV(RealSingle), exps))
        error_exit("gas model");

    /* viscosity model */

    exps[1] = (float)0.0;
    exps[2] = (float)0.0;
    exps[3] = (float)1.0;
    if (cg_goto(cgfile, cgbase, "FlowEquationSet_t", 1,
        "ViscosityModel_t", 1, "end") ||
        cg_array_write("SutherlandLawConstant", CGNS_ENUMV(RealSingle),
            1, &dim, &ts) ||
        cg_gopath(cgfile, "SutherlandLawConstant") ||
        cg_exponents_write(CGNS_ENUMV(RealSingle), exps))
        error_exit("viscosity model");

    exps[0] = (float)1.0;
    exps[1] = (float)-1.0;
    exps[2] = (float)-1.0;
    exps[3] = (float)0.0;
    if (cg_gopath(cgfile, "..") ||
        cg_array_write("ViscosityMolecularReference", CGNS_ENUMV(RealSingle),
            1, &dim, &mu) ||
        cg_gopath(cgfile, "ViscosityMolecularReference") ||
        cg_exponents_write(CGNS_ENUMV(RealSingle), exps))
        error_exit("viscosity model");

    /* thermal conductivity model */

    if (cg_goto(cgfile, cgbase, "FlowEquationSet_t", 1,
        "ThermalConductivityModel_t", 1, "end") ||
        cg_array_write("PowerLawExponent", CGNS_ENUMV(RealSingle),
            1, &dim, &exp) ||
        cg_goto(cgfile, cgbase, "FlowEquationSet_t", 1,
            "ThermalConductivityModel_t", 1, "DataArray_t", 1, "end") ||
        cg_dataclass_write(CGNS_ENUMV(DimensionlessConstant)))
        error_exit("thermal conductivity model");

    /* turbulence closure model */

    if (cg_goto(cgfile, cgbase, "FlowEquationSet_t", 1,
        "TurbulenceClosure_t", 1, "end") ||
        cg_array_write("PrandtlTurbulent", CGNS_ENUMV(RealSingle),
            1, &dim, &pt) ||
        cg_goto(cgfile, cgbase, "FlowEquationSet_t", 1,
            "TurbulenceClosure_t", 1, "DataArray_t", 1, "end") ||
        cg_dataclass_write(CGNS_ENUMV(DimensionlessConstant)))
        error_exit("turbulence closure model");
}

void write_coords(int nz)
{
    int k, nn, n, nij, koff, cgcoord;
    float exps[5];

    koff = nz == 1 ? 1 - NUM_SIDE : 0;
    nij = NUM_SIDE * NUM_SIDE;
    for (n = 0, k = 0; k < NUM_SIDE; k++) {
        for (nn = 0; nn < nij; nn++)
            zcoord[n++] = (float)(k + koff);
    }
    for (n = 0; n < 5; n++)
        exps[n] = (float)0.0;
    exps[1] = (float)1.0;

    if (cg_coord_write(cgfile, cgbase, nz, CGNS_ENUMV(RealSingle),
        "CoordinateX", xcoord, &cgcoord) ||
        cg_goto(cgfile, cgbase, "Zone_t", nz, "GridCoordinates_t", 1,
            "CoordinateX", 0, NULL) ||
        cg_exponents_write(CGNS_ENUMV(RealSingle), exps) ||
        cg_coord_write(cgfile, cgbase, nz, CGNS_ENUMV(RealSingle),
            "CoordinateY", ycoord, &cgcoord) ||
        cg_goto(cgfile, cgbase, "Zone_t", nz, "GridCoordinates_t", 1,
            "CoordinateY", 0, NULL) ||
        cg_exponents_write(CGNS_ENUMV(RealSingle), exps) ||
        cg_coord_write(cgfile, cgbase, nz, CGNS_ENUMV(RealSingle),
            "CoordinateZ", zcoord, &cgcoord) ||
        cg_goto(cgfile, cgbase, "Zone_t", nz, "GridCoordinates_t", 1,
            "CoordinateZ", 0, NULL) ||
        cg_exponents_write(CGNS_ENUMV(RealSingle), exps)) {
        sprintf (errmsg, "zone %d coordinates", nz);
        error_exit(errmsg);
    }
}

void write_elements(int nz)
{
    int cgsect;

    if (cg_section_write(cgfile, cgbase, nz, "Elements", CGNS_ENUMV(HEXA_8),
        1, num_element, 0, elements, &cgsect) ||
        cg_section_write(cgfile, cgbase, nz, "Faces", CGNS_ENUMV(QUAD_4),
            num_element+1, num_element+num_face, 0, faces, &cgsect) ||
        cg_parent_data_write(cgfile, cgbase, nz, cgsect, parent)) {
        sprintf (errmsg, "zone %d elements", nz);
        error_exit(errmsg);
    }
}

void write_zone_link(int nz, char *basename, char *nodename)
{
    char pathname[128];

    sprintf(pathname, "/%s/Zone%d/%s", basename, nz, nodename);
    if (cg_goto(cgfile, cgbase, "Zone_t", nz, "end") ||
        cg_link_write(nodename, "", pathname)) {
        sprintf (errmsg, "zone %d link", nz);
        error_exit(errmsg);
    }
}

/*------------------------------------------------------------------------
* structured grid
*------------------------------------------------------------------------*/

void write_structured()
{
    int i, j, k, n, cgconn, cgbc, cgfam, cggeo, cgpart;
    cgsize_t range[6], d_range[6];
    int transform[3];
    int cgsol, rind[6], cgfld;
    char name[33];
    float exps[5];

    printf ("writing structured base\n");
    fflush (stdout);

    if (cg_base_write(cgfile, "Structured", CellDim, PhyDim, &cgbase) ||
        cg_goto(cgfile, cgbase, "end") ||
        cg_descriptor_write("Descriptor", "Multi-block Structured Grid") ||
        cg_dataclass_write(CGNS_ENUMV(Dimensional)) ||
        cg_units_write(CGNS_ENUMV(Kilogram), CGNS_ENUMV(Meter),
            CGNS_ENUMV(Second), CGNS_ENUMV(Kelvin), CGNS_ENUMV(Radian)))
        error_exit("structured base");
    if (cg_simulation_type_write(cgfile, cgbase, CGNS_ENUMV(NonTimeAccurate)))
        error_exit("simulation type");

    write_reference();
    write_equationset();

    /* write zones */

    for (n = 0; n < 3; n++) {
        size[n]   = NUM_SIDE;
        size[n+3] = NUM_SIDE - 1;
        size[n+6] = 0;
    }
    for (n = 1; n <= 2; n++) {
        sprintf(name, "Zone%d", n);
        if (cg_zone_write(cgfile, cgbase, name, size,
            CGNS_ENUMV(Structured), &cgzone)) {
            sprintf (errmsg, "structured zone %d", n);
            error_exit(errmsg);
        }
        write_coords(n);
    }

    /* write zone 1 to zone 2 connectivity as 1to1 */

    for (n = 0; n < 3; n++) {
        range[n] = d_range[n] = 1;
        range[n+3] = d_range[n+3] = NUM_SIDE;
        transform[n] = n + 1;
    }
    range[2] = NUM_SIDE;
    d_range[5] = 1;
    if (cg_1to1_write(cgfile, cgbase, 1, "1to1 -> Zone2", "Zone2",
        range, d_range, transform, &cgconn))
        error_exit("1to1->zone2");

    /* write zone 2 to zone 1 connectivity as Abutting1to1 */

    for (n = 0; n < 3; n++) {
        range[n] = 1;
        range[n+3] = NUM_SIDE;
    }
    range[5] = 1;
    for (n = 0, j = 1; j <= NUM_SIDE; j++) {
        for (i = 1; i <= NUM_SIDE; i++) {
            d_pts[n++] = i;
            d_pts[n++] = j;
            d_pts[n++] = 1;
        }
    }
    if (cg_conn_write(cgfile, cgbase, 2, "Abutting1to1 -> Zone1",
        CGNS_ENUMV(Vertex), CGNS_ENUMV(Abutting1to1),
        CGNS_ENUMV(PointRange), 2, range, "Zone1",
        CGNS_ENUMV(Structured), CGNS_ENUMV(PointListDonor),
        CGNS_ENUMV(Integer), npts, d_pts, &cgconn))
        error_exit("abutting1to1->zone1");

    /* write inlet BC (zone 1) as point range */

    for (n = 0; n < 3; n++) {
        range[n] = 1;
        range[n+3] = NUM_SIDE;
    }
    range[5] = 1;
    if (cg_boco_write(cgfile, cgbase, 1, "Inlet", CGNS_ENUMV(BCInflow),
        CGNS_ENUMV(PointRange), 2, range, &cgbc))
        error_exit("inlet boco");

    /* write outlet BC (zone 2) as point list */

    for (n = 0, j = 1; j <= NUM_SIDE; j++) {
        for (i = 1; i <= NUM_SIDE; i++) {
            pts[n++] = i;
            pts[n++] = j;
            pts[n++] = NUM_SIDE;
        }
    }
    if (cg_boco_write(cgfile, cgbase, 2, "Outlet", CGNS_ENUMV(BCOutflow),
        CGNS_ENUMV(PointList), npts, pts, &cgbc))
        error_exit("outlet boco");

    /* write zone 1 wall BC as point ranges using a family to group them */

    if (cg_family_write(cgfile, cgbase, "WallFamily", &cgfam) ||
        cg_fambc_write(cgfile, cgbase, cgfam, "WallBC",
            CGNS_ENUMV(BCWall), &cgbc) ||
        cg_goto(cgfile, cgbase, "Family_t", 1, "FamilyBC_t", 1, "end") ||
        cg_bcdataset_write("WallBCDataSet", CGNS_ENUMV(BCWall),
            CGNS_ENUMV(Neumann)))
        error_exit("wall family bc");

    /* write out some bogus geometry info for the family */

    if (cg_geo_write(cgfile, cgbase, cgfam, "Geometry",
        "geometry.file", "CADsystem", &cggeo) ||
        cg_part_write(cgfile, cgbase, cgfam, cggeo, "imin part", &cgpart) ||
        cg_part_write(cgfile, cgbase, cgfam, cggeo, "imax part", &cgpart) ||
        cg_part_write(cgfile, cgbase, cgfam, cggeo, "jmin part", &cgpart) ||
        cg_part_write(cgfile, cgbase, cgfam, cggeo, "jmax part", &cgpart))
        error_exit("wall family parts");

    for (n = 0; n < 3; n++) {
        range[n] = 1;
        range[n+3] = NUM_SIDE;
    }

    range[3] = 1;
    if (cg_boco_write(cgfile, cgbase, 1, "imin", CGNS_ENUMV(FamilySpecified),
        CGNS_ENUMV(PointRange), 2, range, &cgbc) ||
        cg_goto(cgfile, cgbase, "Zone_t", 1, "ZoneBC_t", 1,
            "BC_t", cgbc, "end") ||
        cg_multifam_write ("MultiFam Arg 1", "MultiFam Arg 2") ||
        cg_famname_write("WallFamily"))
        error_exit("imin boco");

    range[0] = range[3] = NUM_SIDE;
    if (cg_boco_write(cgfile, cgbase, 1, "imax", CGNS_ENUMV(FamilySpecified),
        CGNS_ENUMV(PointRange), 2, range, &cgbc) ||
        cg_goto(cgfile, cgbase, "Zone_t", 1, "ZoneBC_t", 1,
            "BC_t", cgbc, "end") ||
        //cg_famname_write("WallFamily"))
        cg_famname_write("WallFamily111"))
        error_exit("imax boco");

    range[0] = range[4] = 1;
    if (cg_boco_write(cgfile, cgbase, 1, "jmin", CGNS_ENUMV(FamilySpecified),
        CGNS_ENUMV(PointRange), 2, range, &cgbc) ||
        cg_goto(cgfile, cgbase, "Zone_t", 1, "ZoneBC_t", 1,
            "BC_t", cgbc, "end") ||
        cg_famname_write("WallFamily"))
        error_exit("jmin boco");

    range[1] = range[4] = NUM_SIDE;
    if (cg_boco_write(cgfile, cgbase, 1, "jmax", CGNS_ENUMV(FamilySpecified),
        CGNS_ENUMV(PointRange), 2, range, &cgbc) ||
        cg_goto(cgfile, cgbase, "Zone_t", 1, "ZoneBC_t", 1,
            "BC_t", cgbc, "end") ||
        cg_famname_write("WallFamily"))
        error_exit("jmax boco");

    /* write zone 2 wall BC as point list */

    for (n = 0, k = 1; k <= NUM_SIDE; k++) {
        for (i = 1; i < NUM_SIDE; i++) {
            pts[n++] = i + 1;
            pts[n++] = 1;
            pts[n++] = k;
            pts[n++] = i;
            pts[n++] = NUM_SIDE;
            pts[n++] = k;
        }
        for (j = 1; j < NUM_SIDE; j++) {
            pts[n++] = 1;
            pts[n++] = j;
            pts[n++] = k;
            pts[n++] = NUM_SIDE;
            pts[n++] = j+1;
            pts[n++] = k;
        }
    }
    if (cg_boco_write(cgfile, cgbase, 2, "Walls", CGNS_ENUMV(BCWall),
        CGNS_ENUMV(PointList), n / 3, pts, &cgbc))
        error_exit("zone 2 walls boco");

    /* write solution for zone 1 as vertex with rind points */

    exps[0] = (float)1.0;
    exps[1] = (float)-3.0;
    for (n = 2; n < 5; n++)
        exps[n] = (float)0.0;
    for (n = 0; n < 6; n++)
        rind[n] = 1;
    if (cg_sol_write(cgfile, cgbase, 1, "VertexSolution",
        CGNS_ENUMV(Vertex), &cgsol) ||
        cg_goto(cgfile, cgbase, "Zone_t", 1, "FlowSolution_t", cgsol, "end") ||
        cg_rind_write(rind) ||
        cg_field_write(cgfile, cgbase, 1, cgsol, CGNS_ENUMV(RealSingle),
            "Density", solution, &cgfld) ||
        cg_gorel(cgfile, "Density", 0, NULL) ||
        cg_exponents_write(CGNS_ENUMV(RealSingle), exps))
        error_exit("zone 1 solution");

    /* write solution for zone 2 as cell center with rind points */

    rind[0] = rind[1] = 0;
    if (cg_sol_write(cgfile, cgbase, 2, "CellCenterSolution",
        CGNS_ENUMV(CellCenter), &cgsol) ||
        cg_goto(cgfile, cgbase, "Zone_t", 2, "FlowSolution_t", cgsol, "end") ||
        cg_rind_write(rind) ||
        cg_field_write(cgfile, cgbase, 2, cgsol, CGNS_ENUMV(RealSingle),
            "Density", solution, &cgfld) ||
        cg_gorel(cgfile, "Density", 0, NULL) ||
        cg_exponents_write(CGNS_ENUMV(RealSingle), exps))
        error_exit("zone 2 solution");
}

/*------------------------------------------------------------------------
* unstructured grid
*------------------------------------------------------------------------*/

void write_unstructured()
{
    int n, nelem, cgconn, cgbc;
    cgsize_t range[2];
    int cgsol, cgfld;
    char name[33];
    float exps[5];
#ifdef UNSTRUCTURED_1TO1
    cgsize_t d_range[2];
    int transform;
#else
# ifdef ABUTTING1TO1_FACES
    GridLocation_t location;
# else
    int i, j;
# endif
#endif

    printf ("writing unstructured base\n");
    fflush (stdout);

    if (cg_base_write(cgfile, "Unstructured", CellDim, PhyDim, &cgbase) ||
        cg_goto(cgfile, cgbase, "end") ||
        cg_descriptor_write("Descriptor", "Multi-block Unstructured Grid") ||
        cg_dataclass_write(CGNS_ENUMV(NormalizedByDimensional)) ||
        cg_units_write(CGNS_ENUMV(Kilogram), CGNS_ENUMV(Meter),
            CGNS_ENUMV(Second), CGNS_ENUMV(Kelvin), CGNS_ENUMV(Radian)))
        error_exit("unstructured base");

    /* write zones */

    for (n = 0; n < 9; n++)
        size[n] = 0;
    size[0] = num_coord;
    size[1] = num_element;
    for (n = 1; n <= 2; n++) {
        sprintf(name, "Zone%d", n);
        if (cg_zone_write(cgfile, cgbase, name, size,
            CGNS_ENUMV(Unstructured), &cgzone)) {
            sprintf (errmsg, "unstructured zone %d", n);
            error_exit(errmsg);
        }
        write_coords(n);
        write_elements(n);
    }
    nelem = (NUM_SIDE - 1) * (NUM_SIDE - 1);

#ifdef UNSTRUCTURED_1TO1

    /* write connectivities as 1to1 */

    range[0] = NODE_INDEX(1, 1, NUM_SIDE);
    range[1] = NODE_INDEX(NUM_SIDE, NUM_SIDE, NUM_SIDE);
    d_range[0] = NODE_INDEX(1, 1, 1);
    d_range[1] = NODE_INDEX(NUM_SIDE, NUM_SIDE, 1);
    transform = 1;
    if (cg_1to1_write(cgfile, cgbase, 1, "1to1 -> Zone2", "Zone2",
        range, d_range, &transform, &cgconn))
        error_exit("1to1->zone2");

    if (cg_1to1_write(cgfile, cgbase, 2, "1to1 -> Zone1", "Zone1",
        d_range, range, &transform, &cgconn))
        error_exit("1to1->zone1");

#else
# ifdef ABUTTING1TO1_FACES

    /* zone 1 to zone 2 connectivity as Abutting1to1 with element range */

    range[0] = num_element + num_face - nelem + 1;
    range[1] = num_element + num_face;
    for (n = 0; n < nelem; n++)
        d_pts[n] = range[0] + n - nelem;

    location = FaceCenter;

    /* this fail for version prior to 2.2 - see below */
    if (cg_conn_write(cgfile, cgbase, 1, "Abutting1to1 -> Zone2",
        location, CGNS_ENUMV(Abutting1to1),
        CGNS_ENUMV(PointRange), 2, range, "Zone2",
        CGNS_ENUMV(Unstructured), CGNS_ENUMV(PointListDonor),
        CGNS_ENUMV(Integer), nelem, d_pts, &cgconn))
        error_exit("face center abutting1to1->zone2");

    /* zone 2 to zone 1 connectivity as Abutting1to1 with element list */

    for (n = 0; n < nelem; n++) {
        pts[n]   = num_element + num_face - 2 * nelem + 1 + n;
        d_pts[n] = pts[n] + nelem;
    }
    if (cg_conn_write(cgfile, cgbase, 2, "Abutting1to1 -> Zone1",
        location, CGNS_ENUMV(Abutting1to1),
        CGNS_ENUMV(PointList), nelem, pts, "Zone1",
        CGNS_ENUMV(Unstructured), CGNS_ENUMV(PointListDonor),
        CGNS_ENUMV(Integer), nelem, d_pts, &cgconn))
        error_exit("face center abutting1to1->zone1");

# else

    /* zone 1 to zone 2 connectivity as Abutting1to1 with point range */

    range[0] = NODE_INDEX(1, 1, NUM_SIDE);
    range[1] = NODE_INDEX(NUM_SIDE, NUM_SIDE, NUM_SIDE);
    for (n = 0, j = 1; j <= NUM_SIDE; j++) {
        for (i = 1; i <= NUM_SIDE; i++)
            d_pts[n++] = NODE_INDEX(i, j, 1);
    }
    if (cg_conn_write(cgfile, cgbase, 1, "Abutting1to1 -> Zone2",
        CGNS_ENUMV(Vertex), CGNS_ENUMV(Abutting1to1),
        CGNS_ENUMV(PointRange), 2, range, "Zone2",
        CGNS_ENUMV(Unstructured), CGNS_ENUMV(PointListDonor),
        CGNS_ENUMV(Integer), npts, d_pts, &cgconn))
        error_exit("point range abutting1to1->zone2");

    /* zone 2 to zone 1 connectivity as Abutting1to1 with point list */

    for (n = 0, j = 1; j <= NUM_SIDE; j++) {
        for (i = 1; i <= NUM_SIDE; i++) {
            pts[n] = d_pts[n];
            d_pts[n++] = NODE_INDEX(i, j, NUM_SIDE);
        }
    }
    if (cg_conn_write(cgfile, cgbase, 2, "Abutting1to1 -> Zone1",
        CGNS_ENUMV(Vertex), CGNS_ENUMV(Abutting1to1),
        CGNS_ENUMV(PointList), npts, pts, "Zone1",
        CGNS_ENUMV(Unstructured), CGNS_ENUMV(PointListDonor),
        CGNS_ENUMV(Integer), npts, d_pts, &cgconn))
        error_exit("point list abutting1to1->zone1");

# endif
#endif

    /* write inlet BC (zone 1) and outlet BC (zone 2) as element lists
    and zone 1 and zone 2 walls as element range */

    range[0] = num_element + 1;
    range[1] = num_element + 4 * nelem;
    for (n = 0; n < nelem; n++) {
        pts[n]   = num_element + num_face - 2 * nelem + 1 + n;
        d_pts[n] = pts[n] + nelem;
    }

    if (cg_boco_write(cgfile, cgbase, 1, "Inlet", CGNS_ENUMV(BCInflow),
        CGNS_ENUMV(ElementList), nelem, pts, &cgbc))
        error_exit ("elementlist inlet boco");
    if (cg_boco_write(cgfile, cgbase, 2, "Outlet", CGNS_ENUMV(BCOutflow),
        CGNS_ENUMV(ElementList), nelem, d_pts, &cgbc))
        error_exit ("elementlist outlet boco");
    if (cg_boco_write(cgfile, cgbase, 1, "Walls", CGNS_ENUMV(BCWall),
        CGNS_ENUMV(ElementRange), 2, range, &cgbc))
        error_exit ("elementrange zone 1 walls boco");
    if (cg_boco_write(cgfile, cgbase, 2, "Walls", CGNS_ENUMV(BCWall),
        CGNS_ENUMV(ElementRange), 2, range, &cgbc))
        error_exit ("elementrange zone 2 walls boco");

    /* write solution for zone 1 as vertex and zone 2 as cell center */

    exps[0] = (float)1.0;
    exps[1] = (float)-3.0;
    for (n = 2; n < 5; n++)
        exps[n] = (float)0.0;
    if (cg_sol_write(cgfile, cgbase, 1, "VertexSolution",
        CGNS_ENUMV(Vertex), &cgsol) ||
        cg_field_write(cgfile, cgbase, 1, cgsol, CGNS_ENUMV(RealSingle),
            "Density", solution, &cgfld) ||
        cg_goto(cgfile, cgbase, "Zone_t", 1, "FlowSolution_t", cgsol,
            "Density", 0, NULL) ||
        cg_exponents_write(CGNS_ENUMV(RealSingle), exps))
        error_exit("zone 1 solution");
    if (cg_sol_write(cgfile, cgbase, 2, "CellCenterSolution",
        CGNS_ENUMV(CellCenter), &cgsol) ||
        cg_field_write(cgfile, cgbase, 2, cgsol, CGNS_ENUMV(RealSingle),
            "Density", solution, &cgfld) ||
        cg_goto(cgfile, cgbase, "Zone_t", 2, "FlowSolution_t", cgsol,
            "Density", 0, NULL) ||
        cg_exponents_write(CGNS_ENUMV(RealSingle), exps))
        error_exit("zone 2 solution");
}

/*------------------------------------------------------------------------
* mixed grid
*------------------------------------------------------------------------*/

void write_mixed()
{
    int i, j, k, n, nelem, cgconn, cgbc;
    cgsize_t range[6];
#ifdef ABUTTING1TO1_FACES
    GridLocation_t location;
#endif

    printf ("writing mixed base\n");
    fflush (stdout);

    if (cg_base_write(cgfile, "Mixed", CellDim, PhyDim, &cgbase) ||
        cg_goto(cgfile, cgbase, "end") ||
        cg_descriptor_write("Descriptor",
            "Mixed Structured and Unstructured Grid") ||
        cg_dataclass_write(CGNS_ENUMV(Dimensional)) ||
        cg_units_write(CGNS_ENUMV(Kilogram), CGNS_ENUMV(Meter),
            CGNS_ENUMV(Second), CGNS_ENUMV(Kelvin), CGNS_ENUMV(Radian)))
        error_exit("mixed base");

    /* zone 1 is structured */

    for (n = 0; n < 3; n++) {
        size[n]   = NUM_SIDE;
        size[n+3] = NUM_SIDE - 1;
        size[n+6] = 0;
    }
    if (cg_zone_write(cgfile, cgbase, "StructuredZone", size,
        CGNS_ENUMV(Structured), &cgzone))
        error_exit("structured zone");
    write_zone_link(1, "Structured", "GridCoordinates");

    /* zone 2 is unstructured */

    for (n = 0; n < 9; n++)
        size[n] = 0;
    size[0] = num_coord;
    size[1] = num_element;
    if (cg_zone_write(cgfile, cgbase, "UnstructuredZone", size,
        CGNS_ENUMV(Unstructured), &cgzone))
        error_exit("unstructured zone");
    write_zone_link(2, "Unstructured", "GridCoordinates");
    write_zone_link(2, "Unstructured", "Elements");
    write_zone_link(2, "Unstructured", "Faces");
    nelem = (NUM_SIDE - 1) * (NUM_SIDE - 1);

#ifdef ABUTTING1TO1_FACES

    /* zone 1 -> zone 2 connectivity as face range */

    for (n = 0; n < 3; n++) {
        range[n] = 1;
        range[n+3] = NUM_SIDE - 1;
    }
    range[2] = range[5] = NUM_SIDE;
    for (n = 0; n < nelem; n++)
        d_pts[n] = num_element + num_face - 2 * nelem + 1 + n;

    location = KFaceCenter;
    if (cg_conn_write(cgfile, cgbase, 1, "Structured -> Unstructured",
        location, CGNS_ENUMV(Abutting1to1),
        CGNS_ENUMV(PointRange), 2, range, "UnstructuredZone",
        CGNS_ENUMV(Unstructured), CGNS_ENUMV(PointListDonor),
        CGNS_ENUMV(Integer), nelem, d_pts, &cgconn))
        error_exit("abutting1to1->unstructured face range");

#else

    /* zone 1 -> zone 2 connectivity as point range */

    for (n = 0; n < 3; n++) {
        range[n] = 1;
        range[n+3] = NUM_SIDE;
    }
    range[2] = NUM_SIDE;
    for (n = 0; n < npts; n++)
        d_pts[n] = n + 1;
    if (cg_conn_write(cgfile, cgbase, 1, "Structured -> Unstructured",
        CGNS_ENUMV(Vertex), CGNS_ENUMV(Abutting1to1),
        CGNS_ENUMV(PointRange), 2, range, "UnstructuredZone",
        CGNS_ENUMV(Unstructured), CGNS_ENUMV(PointListDonor),
        CGNS_ENUMV(Integer), npts, d_pts, &cgconn))
        error_exit("abutting1to1->unstructured point range");

#endif

    /* zone 2 -> zone 1 connectivity as point range (k=1 surface) */

    range[0] = 1;
    range[1] = npts;
    for (n = 0, j = 1; j <= NUM_SIDE; j++) {
        for (i = 1; i <= NUM_SIDE; i++) {
            d_pts[n++] = i;
            d_pts[n++] = j;
            d_pts[n++] = NUM_SIDE;
        }
    }
    if (cg_conn_write(cgfile, cgbase, 2, "Unstructured -> Structured",
        CGNS_ENUMV(Vertex), CGNS_ENUMV(Abutting1to1),
        CGNS_ENUMV(PointRange), 2, range, "StructuredZone",
        CGNS_ENUMV(Structured), CGNS_ENUMV(PointListDonor),
        CGNS_ENUMV(Integer), npts, d_pts, &cgconn))
        error_exit("abutting1to1->structured point range");

#ifdef STRUCTURED_FACES

    /* write inlet BC (zone 1) as face range */

    for (n = 0; n < 3; n++) {
        range[n] = 1;
        range[n+3] = NUM_SIDE - 1;
    }
    range[5] = 1;
    if (cg_boco_write(cgfile, cgbase, 1, "Inlet", CGNS_ENUMV( BCInflow ),
        CGNS_ENUMV(PointRange), 2, range, &cgbc) ||
        cg_goto(cgfile, cgbase, "Zone_t", 1, "ZoneBC_t", 1,
            "BC_t", cgbc, "end") ||
        cg_gridlocation_write(CGNS_ENUMV(KFaceCenter))
        error_exit("kfacecenter inlet boco");

#else

    /* write inlet BC (zone 1) as point range */

    for (n = 0; n < 3; n++) {
        range[n] = 1;
        range[n+3] = NUM_SIDE;
    }
    range[5] = 1;
    if (cg_boco_write(cgfile, cgbase, 1, "Inlet", CGNS_ENUMV(BCInflow),
        CGNS_ENUMV(PointRange), 2, range, &cgbc))
        error_exit("point range inlet boco");

#endif

    /* write outlet BC (zone 2) as point range */

    range[0] = num_coord - npts + 1;
    range[1] = num_coord;
    if (cg_boco_write(cgfile, cgbase, 2, "Outlet", CGNS_ENUMV(BCOutflow),
        CGNS_ENUMV(PointRange), 2, range, &cgbc))
        error_exit("point range outlet boco");

#ifdef STRUCTURED_FACES

    /* write zone 1 walls as face list */

    for (n = 0, k = 1; k < NUM_SIDE; k++) {
        for (i = 1; i < NUM_SIDE; i++) {
            pts[n++] = i;
            pts[n++] = 1;
            pts[n++] = k;
            pts[n++] = i;
            pts[n++] = NUM_SIDE;
            pts[n++] = k;
        }
        for (j = 1; j < NUM_SIDE; j++) {
            pts[n++] = 1;
            pts[n++] = j;
            pts[n++] = k;
            pts[n++] = NUM_SIDE;
            pts[n++] = j;
            pts[n++] = k;
        }
    }
    if (cg_boco_write(cgfile, cgbase, 1, "Walls", CGNS_ENUMV(BCWall),
        CGNS_ENUMV(PointList), 4 * nelem, pts, &cgbc) ||
        cg_goto(cgfile, cgbase, "Zone_t", 1, "ZoneBC_t", 1,
            "BC_t", cgbc, "end") ||
        cg_gridlocation_write(CGNS_ENUMV(FaceCenter)))
        error_exit("face list zone 1 walls boco");

#else

    /* write zone 1 wall BC as point list */

    for (n = 0, k = 1; k <= NUM_SIDE; k++) {
        for (i = 1; i < NUM_SIDE; i++) {
            pts[n++] = i + 1;
            pts[n++] = 1;
            pts[n++] = k;
            pts[n++] = i;
            pts[n++] = NUM_SIDE;
            pts[n++] = k;
        }
        for (j = 1; j < NUM_SIDE; j++) {
            pts[n++] = 1;
            pts[n++] = j;
            pts[n++] = k;
            pts[n++] = NUM_SIDE;
            pts[n++] = j+1;
            pts[n++] = k;
        }
    }
    if (cg_boco_write(cgfile, cgbase, 1, "Walls", CGNS_ENUMV(BCWall),
        CGNS_ENUMV(PointList), n / 3, pts, &cgbc))
        error_exit("point list zone 1 walls boco");

#endif

    /* write zone 2 walls as point list */

    for (n = 0, k = 1; k <= NUM_SIDE; k++) {
        for (i = 1; i < NUM_SIDE; i++) {
            pts[n++] = NODE_INDEX(i + 1, 1, k);
            pts[n++] = NODE_INDEX(i, NUM_SIDE, k);
        }
        for (j = 1; j < NUM_SIDE; j++) {
            pts[n++] = NODE_INDEX(1, j, k);
            pts[n++] = NODE_INDEX(NUM_SIDE, j + 1, k);
        }
    }
    if (cg_boco_write(cgfile, cgbase, 2, "Walls", CGNS_ENUMV(BCWall),
        CGNS_ENUMV(PointList), n, pts, &cgbc))
        error_exit("point list zone 2 walls boco");
}

/*------------------------------------------------------------------------
* mismatched zones with cylindrical coordinates
*------------------------------------------------------------------------*/

void write_mismatched()
{
    int i, j, k, n, nj, cgcoord, cgconn, cgbc;
    int ic, jc;
    cgsize_t dims[2];
    float dx, dy, dtheta, r, t, x, y;
    cgsize_t range[6], d_range[6];
    int transform[3];
    CGNS_ENUMT(PointSetType_t) d_type;
    float exps[5];

    printf ("writing mismatched base\n");
    fflush (stdout);

    if (cg_base_write(cgfile, "Mismatched", CellDim, PhyDim, &cgbase) ||
        cg_goto(cgfile, cgbase, "end") ||
        cg_descriptor_write("Descriptor", "Mismatched Grid") ||
        cg_dataclass_write(CGNS_ENUMV(Dimensional)) ||
        cg_units_write(CGNS_ENUMV(Kilogram), CGNS_ENUMV(Meter),
            CGNS_ENUMV(Second), CGNS_ENUMV(Kelvin), CGNS_ENUMV(Radian)))
        error_exit("mismatched base");

    /* zone 1 is cartesian */

    for (n = 0; n < 3; n++) {
        size[n]   = NUM_SIDE;
        size[n+3] = NUM_SIDE - 1;
        size[n+6] = 0;
    }

    dx = 0.5 * (float)(NUM_SIDE - 1);
    dy = 0.5 * (float)(NUM_SIDE - 1);
    for (n = 0, k = 0; k < NUM_SIDE; k++) {
        for (j = 0; j < NUM_SIDE; j++) {
            for (i = 0; i < NUM_SIDE; i++, n++) {
                xcoord[n] -= dx;
                ycoord[n] -= dy;
            }
        }
    }
    if (cg_zone_write(cgfile, cgbase, "CartesianZone", size,
        CGNS_ENUMV(Structured), &cgzone))
        error_exit("cartesion zone");
    write_coords(1);

    /* zone 2 is cylindrical */

    nj = 2 * NUM_SIDE;
    dtheta = (float)TWOPI / (float)(nj - 1);
    for (n = 0, k = 0; k < NUM_SIDE; k++) {
        for (j = 0; j < nj; j++) {
            for (i = 0; i < NUM_SIDE; i++, n++) {
                xcoord[n] = (float)i;
                ycoord[n] = (float)j * dtheta;
                zcoord[n] = (float)k;
            }
        }
    }
    size[1] = nj;
    size[4] = nj - 1;
    for (n = 0; n < 5; n++)
        exps[n] = (float)0.0;
    exps[1] = (float)1.0;
    if (cg_zone_write(cgfile, cgbase, "CylindricalZone", size,
        CGNS_ENUMV(Structured), &cgzone) ||
        cg_coord_write(cgfile, cgbase, cgzone, CGNS_ENUMV(RealSingle),
            "CoordinateR", xcoord, &cgcoord) ||
        cg_goto(cgfile, cgbase, "Zone_t", cgzone, "GridCoordinates", 0,
            "CoordinateR", 0, NULL) ||
        cg_exponents_write(CGNS_ENUMV(RealSingle), exps) ||
        cg_coord_write(cgfile, cgbase, cgzone, CGNS_ENUMV(RealSingle),
            "CoordinateZ", zcoord, &cgcoord) ||
        cg_goto(cgfile, cgbase, "Zone_t", cgzone, "GridCoordinates", 0,
            "CoordinateZ", 0, NULL) ||
        cg_exponents_write(CGNS_ENUMV(RealSingle), exps))
        error_exit("cylindrical zone");
    exps[1] = (float)0.0;
    exps[4] = (float)1.0;
    if (cg_coord_write(cgfile, cgbase, cgzone, CGNS_ENUMV(RealSingle),
        "CoordinateTheta", ycoord, &cgcoord) ||
        cg_goto(cgfile, cgbase, "Zone_t", cgzone, "GridCoordinates", 0,
            "CoordinateTheta", 0, NULL) ||
        cg_exponents_write(CGNS_ENUMV(RealSingle), exps))
        error_exit("cylindrical zone");

    /* zone 1 -> zone 2 connectivity */

    for (n = 0, j = 0; j < NUM_SIDE; j++) {
        for (i = 0; i < NUM_SIDE; i++) {
            x = (float)i - dx;
            y = (float)j - dy;
            r = (float)sqrt (x * x + y * y);
            if (r > (float)(NUM_SIDE - 1)) continue;
            ic = (int)r;
            if (ic >= NUM_SIDE - 1) ic = NUM_SIDE - 2;
            t = (float)atan2 (y, x);
            if (t < 0.0) t += (float)TWOPI;
            jc = (int)(t / dtheta);
            if (jc >= nj - 1) jc = nj - 2;
            pts[n]      = i + 1;
            pts[n+1]    = j + 1;
            pts[n+2]    = NUM_SIDE;
            d_pts[n]    = ic + 1;
            d_pts[n+1]  = jc + 1;
            d_pts[n+2]  = 1;
            interp[n++] = r - (float)ic;
            interp[n++] = t / dtheta - (float)jc;
            interp[n++] = 0.0;
        }
    }

    dims[0] = 3;
    dims[1] = n / 3;
    d_type = CGNS_ENUMV(CellListDonor);
    if (cg_conn_write(cgfile, cgbase, 1, "Cartesian -> Cylindrical",
        CGNS_ENUMV(Vertex), CGNS_ENUMV(Abutting),
        CGNS_ENUMV(PointList), n/3, pts, "CylindricalZone",
        CGNS_ENUMV(Structured), d_type,
        CGNS_ENUMV(Integer), n/3, d_pts, &cgconn) ||
        cg_goto(cgfile, cgbase, "Zone_t", 1, "ZoneGridConnectivity_t", 1,
            "GridConnectivity_t", cgconn, "end") ||
        cg_array_write("InterpolantsDonor", CGNS_ENUMV(RealSingle),
            2, dims, interp))
        error_exit("cartesian->cylindrical connectivity");

    /* zone 2 -> zone 1 connectivity */

    for (n = 0, j = 0; j < nj; j++) {
        for (i = 0; i < NUM_SIDE; i++) {
            r = (float)i;
            t = (float)j * dtheta;
            x = r * (float)cos (t) + dx;
            y = r * (float)sin (t) + dy;
            if (x < 0.0 || x > (float)(NUM_SIDE - 1) ||
                y < 0.0 || y > (float)(NUM_SIDE - 1)) continue;
            ic = (int)x;
            if (ic >= NUM_SIDE - 1) ic = NUM_SIDE - 2;
            jc = (int)y;
            if (jc >= NUM_SIDE - 1) jc = NUM_SIDE - 2;
            pts[n]      = i + 1;
            pts[n+1]    = j + 1;
            pts[n+2]    = 1;
            d_pts[n]    = ic + 1;
            d_pts[n+1]  = jc + 1;
            d_pts[n+2]  = NUM_SIDE;
            interp[n++] = x - (float)ic;
            interp[n++] = y - (float)jc;
            interp[n++] = 0.0;
        }
    }

    dims[0] = 3;
    dims[1] = n / 3;
    d_type = CGNS_ENUMV( CellListDonor );
    if (cg_conn_write(cgfile, cgbase, 2, "Cylindrical -> Cartesian",
        CGNS_ENUMV(Vertex), CGNS_ENUMV(Abutting),
        CGNS_ENUMV(PointList), n/3, pts, "CartesianZone",
        CGNS_ENUMV(Structured) , d_type,
        CGNS_ENUMV(Integer), n/3, d_pts, &cgconn) ||
        cg_goto(cgfile, cgbase, "Zone_t", 2, "ZoneGridConnectivity_t", 1,
            "GridConnectivity_t", cgconn, "end") ||
        cg_array_write("InterpolantsDonor", CGNS_ENUMV(RealSingle),
            2, dims, interp))
        error_exit("cylindrical->cartesian connectivity");

    /* periodic boundary for zone 2 */

    for (n = 0; n < 3; n++) {
        transform[n] = n + 1;
        range[n] = d_range[n] = 1;
        range[n+3] = d_range[n+3] = NUM_SIDE;
    }
    range[4] = 1;
    d_range[1] = d_range[4] = nj;
    if (cg_1to1_write(cgfile, cgbase, 2, "Periodic", "CylindricalZone",
        range, d_range, transform, &cgconn))
        error_exit("periodic 1to1");

    /* write inlet BC (zone 1) as point range */

    range[4] = NUM_SIDE;
    range[5] = 1;
    if (cg_boco_write(cgfile, cgbase, 1, "Inlet", CGNS_ENUMV(BCInflow),
        CGNS_ENUMV(PointRange), 2, range, &cgbc))
        error_exit("inlet boco");

    /* write outlet BC (zone 2) as point range */

    range[4] = nj;
    range[2] = range[5] = NUM_SIDE;
    if (cg_boco_write(cgfile, cgbase, 2, "Outlet", CGNS_ENUMV(BCOutflow),
        CGNS_ENUMV(PointRange), 2, range, &cgbc))
        error_exit("outlet boco");

    /* write zone 1 wall BC as point list */

    for (n = 0, k = 1; k <= NUM_SIDE; k++) {
        for (i = 1; i < NUM_SIDE; i++) {
            pts[n++] = i + 1;
            pts[n++] = 1;
            pts[n++] = k;
            pts[n++] = i;
            pts[n++] = NUM_SIDE;
            pts[n++] = k;
        }
        for (j = 1; j < NUM_SIDE; j++) {
            pts[n++] = 1;
            pts[n++] = j;
            pts[n++] = k;
            pts[n++] = NUM_SIDE;
            pts[n++] = j+1;
            pts[n++] = k;
        }
    }
    if (cg_boco_write(cgfile, cgbase, 1, "Walls", CGNS_ENUMV(BCWall),
        CGNS_ENUMV(PointList), n/3, pts, &cgbc))
        error_exit("zone 1 walls boco");

    /* write zone 2 wall BC as point list */

    for (n = 0, k = 1; k <= NUM_SIDE; k++) {
        for (j = 1; j <= nj; j++) {
            pts[n++] = NUM_SIDE;
            pts[n++] = j;
            pts[n++] = k;
        }
    }
    if (cg_boco_write(cgfile, cgbase, 2, "Walls", CGNS_ENUMV(BCWall),
        CGNS_ENUMV(PointList), n/3, pts, &cgbc))
        error_exit("zone 2 walls boco");
}

EndNameSpace