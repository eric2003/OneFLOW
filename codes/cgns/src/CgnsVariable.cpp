/*-----------------------this->----------------------------------------------------*\
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

#include "CgnsVariable.h"
#include "CgnsBase.h"
#include <iostream>
#include <iomanip>
using namespace std;

BeginNameSpace( ONEFLOW )

CgnsVector::CgnsVector()
{
    data = 0;
}

CgnsVector::~CgnsVector()
{
    delete[ ] this->data;
}

void CgnsVector::Create()
{
    cout << " CgnsVector::Create()\n";
    int nSize = dims[ 0 ];
    for ( int i = 1; i < ndim; ++ i )
    {
        nSize *= dims[ i ];
    }

    this->data = new VEC_DATA[ nSize ];
    cout << " nSize = " << nSize << "\n";
}

void CgnsVector::ReadArray()
{
    this->ReadArrayInfo( this->arrayId );
    this->Create();
    this->ReadArrayContent();
}

void CgnsVector::ReadArrayInfo( int arrayId )
{
    this->arrayId = arrayId;
    cg_array_info( arrayId, name, &dataType, & this->ndim, this->dims );

    cout << " iArray = " << arrayId << " arrayName = " << name << "\n";
    cout << " ndim = " << ndim << " dims = ";
    for ( int i = 0; i < ndim; ++ i )
    {
        cout << dims[ i ] << " ";
    }
    cout << "\n";
}

void CgnsVector::ReadArrayContent()
{
    cg_array_read_as( this->arrayId, CGNS_ENUMV( RealDouble ), data );
    this->PrintData();
}

void CgnsVector::PrintData()
{
    int icount = 0;
    for ( int i = 0; i < ndim; ++ i )
    {
        int n = dims[ i ];
        for ( int j = 0; j < n; ++ j )
        {
            cout << this->data[ icount ++ ] << " ";
        }
        cout << "\n";
    }
}

CgnsZVector::CgnsZVector()
{
}

CgnsZVector::~CgnsZVector()
{
    int nArrays = cgnsVectorList.size();
    for ( int iArray = 0; iArray < nArrays; ++ iArray )
    {
        delete cgnsVectorList[ iArray ];
    }
}

CgnsUserData::CgnsUserData( CgnsBase * cgnsBase )
{
    this->cgnsBase = cgnsBase;
}

CgnsUserData::~CgnsUserData()
{
    int nSize = cgnsZVectorList.size();
    for ( int i = 0; i < nSize; ++ i )
    {
        delete cgnsZVectorList[ i ];
    }
}

void CgnsUserData::ReadUserData()
{
    cgnsBase->GoToBase();
    int nUserData = -1;
    cg_nuser_data( & nUserData );

    cout << " nUserData = " << nUserData << "\n";

    for ( int iData = 0; iData < nUserData; ++ iData )
    {
        int iDataId = iData + 1;
        cgnsBase->GoToBase();

        char name[ 33 ];
        cg_user_data_read( iDataId, name );

        cout << " iData = " << iData << " user data name = " << name << "\n";

        cgnsBase->GoToNode( "UserDefinedData_t", iDataId );

        CgnsZVector * cgnsZVector = new CgnsZVector();
        cgnsZVector->ReadArray();
        cgnsZVectorList.push_back( cgnsZVector );
    }
}

void CgnsZVector::ReadArray()
{
    int nArrays = -1;
    cg_narrays( & nArrays );
    cout << " narrays = " << nArrays << "\n";
    this->ReadArray( nArrays );
}

void CgnsZVector::ReadArray( int nArrays )
{
    for ( int iArray = 0; iArray < nArrays; ++ iArray )
    {
        int A = iArray + 1;
        CgnsVector * cgnsVar = new CgnsVector();
        cgnsVar->arrayId = A;
        cgnsVar->ReadArray();

        cgnsVectorList.push_back( cgnsVar );
    }
}

//CgnsVariable::CgnsVariable()
//{
//}
//
//CgnsVariable::~CgnsVariable()
//{
//}

EndNameSpace