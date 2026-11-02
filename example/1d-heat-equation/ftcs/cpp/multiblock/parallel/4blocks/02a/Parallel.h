#pragma once
#ifdef HX_PARALLEL
#include "mpi.h"
#endif
#include <vector>

class Parallel
{
public:
    static int pid;
    static int nProc;
    static int serverid;
    static int tag;
public:
    static void Init();
    static void Finalize();
};