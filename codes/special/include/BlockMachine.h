#pragma once
#include "HXDefine.h"

BeginNameSpace( ONEFLOW )

class TextFileParser;

class BlockMachine
{
public:
    BlockMachine() = default;
    ~BlockMachine() = default;
public:
    void ApplyRelation( TextFileParser & textFileParser );
    void AddLineToFace( int faceId, int position, int lineId );
    void AddFaceToBlock( int blockId, int position, int faceId );
    void GenerateGrid();
};

extern BlockMachine block_Machine;

EndNameSpace
