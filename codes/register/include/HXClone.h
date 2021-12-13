/*---------------------------------------------------------------------------*\
    OneFLOW - LargeScale Multiphysics Scientific Simulation Environment
    Copyright (C) 2017-2021 He Xin and the OneFLOW contributors.
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


#pragma once
#include "Configure.h"
#include "HXDefine.h"
#include <map>
#include <string>


BeginNameSpace( ONEFLOW )
class HXClone;

#define IMPLEMENT_CLASS_CLONE( TYPE ) \
HXClone * Clone() const { return new TYPE( * this ); }

#define REGISTER_CLASS( TYPE ) \
    HXClone * TYPE ## _myClass = \
        HXClone::Register( #TYPE, new TYPE() );

#define C_CLASS( TYPE ) C##TYPE

#define REGISTER_DATA_CLASS( TYPE ) \
     REGISTER_CLASS( C##TYPE )

class HXClone
{
public:
    virtual ~HXClone() {}
public:
    virtual HXClone * Clone() const = 0;
public:
    static HXClone * SafeClone( const std::string & type );
    static HXClone * Register( const std::string & type, HXClone * clone );
    static std::map < std::string, HXClone * > * classMap;
    StringField data;
public:
    virtual void Solve(){};
};

#define DEFINE_CLASS( TYPE ) \
void HX##TYPE(); \
class TYPE : public HXClone \
{ \
public: \
    IMPLEMENT_CLASS_CLONE( TYPE ) \
public: \
    void Solve() override \
    { \
        HX##TYPE();\
    } \
};


#define DEFINE_DATA_CLASS( TYPE ) \
void TYPE( StringField & data ); \
class C##TYPE : public HXClone \
{ \
public: \
    IMPLEMENT_CLASS_CLONE( C##TYPE ) \
public: \
    void Solve() override \
    { \
        TYPE( this->data );\
    } \
};


EndNameSpace
