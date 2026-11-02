#pragma once
#include "HXDefine.h"

BeginNameSpace( ONEFLOW )

class BlockMachine
{
public:
    BlockMachine() = default;
    ~BlockMachine() = default;
public:
    void AddLineToFace( int faceId, int position, int lineId );
    void AddFaceToBlock( int blockId, int position, int faceId );
    void GenerateGrid();
    void Reset();
};

extern BlockMachine block_Machine;

EndNameSpace
