#include "global.h"

std::vector<Grid *> Global::grids;
std::vector<Field *> Global::fields;
std::vector<BC *> Global::bcs;
std::vector<InterFaceZone *> Global::interfacezones;
PointFactory Global::boundary_points;
FaceList Global::boundary_faces;
FaceList Global::inter_faces;
BaseZoneList Global::zone_names;
int Global::nt = -1;