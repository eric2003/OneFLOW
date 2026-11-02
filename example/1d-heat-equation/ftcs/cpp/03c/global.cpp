#include "global.h"

std::vector<Grid *> Global::grids;
std::vector<Field *> Global::fields;
std::vector<BC *> Global::bcs;
std::vector<InterFaceZone *> Global::interfacezones;
PointFactory Global::boundary_points;
FaceList Global::boundary_faces;
FaceList Global::inter_faces;

int Global::nt = -1;
int Global::cell_dim = -1;
int Global::phys_dim = -1;

std::vector<Zone *> Global::zones;
BaseZoneList Global::zone_names;