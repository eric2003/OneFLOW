#include "cgnslib.h"
#include "CgnsGrid.h"
#include "Partition.h"
#include <iostream>
#include <vector>
#include <numbers>
#include <cmath>
#include <fstream>
#include <iomanip>

int main(int argc, char **argv)
{
    Partition * partition = new Partition();
    partition->Solve();
    delete partition;

    return 0;
}
