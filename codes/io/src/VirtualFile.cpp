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

#include "VirtualFile.h"
#include "DataBaseIO.h"
#include "DataBook.h"
#include <iostream>
using namespace std;

BeginNameSpace( ONEFLOW )

VirtualFile::VirtualFile( fstream * file )
{
    this->file = file;
    this->type = 0;
}

VirtualFile::VirtualFile( DataBook * databook )
{
    this->databook = databook;
    this->type = 1;
}

void VirtualFile::Read( void * data, streamsize size )
{
    char * charData = reinterpret_cast< char * > ( data );
    if ( this->type == 0 )
    {
        ONEFLOW::HXRead( this->file, charData, size );
    }
    else
    {
        ONEFLOW::HXRead( this->databook, charData, size );
    }
}

void VirtualFile::Write( void * data, streamsize size )
{
    char * charData = reinterpret_cast< char * > ( data );

    if ( this->type == 0 )
    {
        ONEFLOW::HXWrite( this->file, charData, size );
    }
    else
    {
        ONEFLOW::HXWrite( this->databook, charData, size );
    }
}

void VirtualFile::MarkSectionBegin()
{
    this->sectionBegin = this->GetCurrentPosition();
}

void VirtualFile::MarkSectionEnd()
{
    this->sectionEnd = this->GetCurrentPosition();
}

void VirtualFile::MoveToPosition( streamsize sectionPosition )
{
    file->seekp( sectionPosition );
}

void VirtualFile::MoveToSectionBegin()
{
    this->MoveToPosition( this->GetSectionBegin() );
}

void VirtualFile::MoveToSectionEnd()
{
    this->MoveToPosition( this->GetSectionEnd() );
}

void VirtualFile::BeginWriteWork()
{
    if ( this->type == 0 )
    {
        this->MarkSectionBegin();
        this->ReservePlaceholder();
    }
}

void VirtualFile::EndWriteWork()
{
    if ( this->type == 0 )
    {
        this->MarkSectionEnd();
        this->MoveToSectionBegin();
        this->WriteDataLength();
        this->MoveToSectionEnd();
    }
}

void VirtualFile::BeginReadWork()
{
    if ( this->type == 0 )
    {
        this->MarkSectionBegin();
        this->ReadDataLength();
    }
    else
    {
        this->databook->MoveToBegin();
    }
}

void VirtualFile::EndReadWork()
{
    if ( this->type == 0 )
    {
        this->MarkSectionEnd();
    }
}

void VirtualFile::ReadDataLength()
{
    this->Read( & dataLength, sizeof( streamsize ) );
}

void VirtualFile::WriteDataLength()
{
    dataLength = GetSectionEnd() - GetSectionBegin() - sizeof( streamsize );

    this->Write( & dataLength, sizeof( streamsize ) );
}

void VirtualFile::ReservePlaceholder()
{
    this->Write( & sectionBegin, sizeof( streamsize ) );
}

streamsize VirtualFile::GetCurrentPosition()
{
    if ( this->type == 0 )
    {
        return file->tellp();
    }
    else
    {
        return 0;
    }
}


EndNameSpace
