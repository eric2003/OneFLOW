/*---------------------------------------------------------------------------*\
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

#include "Category.h"

BeginNameSpace( ONEFLOW )

map< int, int > * Category::data = 0;
Category::Category()
{
    ;
}

Category::~Category()
{
    ;
}

void Category::Init()
{
    if ( ! Category::data )
    {
        Category::data = new map< int, int >();
    }
}

void Category::Free()
{
    if ( ! Category::data ) return;
    delete Category::data;
    Category::data = 0;
}

void Category::AddCategory( int sTid, int category )
{
    Category::Init();
    map< int, int >::iterator iter = Category::data->find( sTid );
    if ( iter == Category::data->end() )
    {
        ( * Category::data )[ sTid ] = category;
    }
}

int Category::GetCategory( int sTid )
{
    map< int, int >::iterator iter = Category::data->find( sTid );
    return iter->second;
}

EndNameSpace