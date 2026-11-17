#pragma once

#include "HXDefine.h"

#include <string>

BeginNameSpace( ONEFLOW )

class GridLayout;

class GridLayoutParser
{
public:
    GridLayout Parse( const std::string & fileName ) const;
};

EndNameSpace
