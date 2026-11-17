#pragma once
#include "HXDefine.h"

BeginNameSpace( ONEFLOW )

class TextFileParser;

class DomainMachine
{
public:
    DomainMachine();
    ~DomainMachine();
public:
    IntField bctypeList, bcLineList;
public:
    void Reset();
    void AddBcType( TextFileParser & textFileParser );
    void SetBcType( int id, int bctype );
    int GetBcType( int id );
};

extern DomainMachine domain_Machine;

EndNameSpace
