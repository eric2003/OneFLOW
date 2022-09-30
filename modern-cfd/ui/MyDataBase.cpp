/*---------------------------------------------------------------------------*\
OneFLOW - LargeScale Multiphysics Scientific Simulation Environment
Copyright (C) 2017-2022 He Xin and the OneFLOW contributors.
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
#include "MyDataBase.h"
#include <QtSql>
#include <QDebug>
#include <QString>

MyDataBase::MyDataBase()
{
    this->connName = "first";
}

MyDataBase::~MyDataBase()
{
    ;
}

bool MyDataBase::CreateConnDataBase( const QString &connName )
{
    qDebug() << QSqlDatabase::drivers();
    this->db = QSqlDatabase::addDatabase( "QMYSQL", connName );
    this->db.setHostName("127.0.0.1"); // Host name
    this->db.setUserName("root");     // User name
    this->db.setPassword("123456");   // Password
    this->db.setPort(3306);           // Port to connect

    if ( ! this->db.open() )
    {
        qDebug() << "Error: " << this->db.lastError().text();
        return false;
    }
    else
    {
        qDebug() << connName << " Database Opened!";
        return false;
    }
}

bool MyDataBase::CreateNewDataBase( const QString & dataBaseName )
{
    if ( ! this->db.isValid() )
    {
        this->CreateConnDataBase( this->connName );
    }
    QString cs;
    cs.append( "create database if not exists " );
    cs.append( dataBaseName );

    this->db.exec( cs );
    this->db.setDatabaseName( dataBaseName );
    return true;
}

void MyDataBase::Run()
{
    this->CreateNewDataBase( "cfd_db" );

    if ( ! this->db.open() )
    {
        qDebug() << "Error: " << db.lastError().text();
        return;
    }

    QSqlQuery query( this->db );
    QString tableName;
    tableName = "cfdprj";
    QString queryString;
    queryString.append( "create table " );
    queryString.append( tableName );
    queryString.append( " (id int primary key, prjname varchar(255) )" );

    bool flag = false;

    flag = query.exec( queryString );

    qDebug() << "---- insert operation start----";
    this->InsertPrjName( tableName, 1, "Sod's Shock Tube");
    this->InsertPrjName( tableName, 2, "Incompressible Driven Cavity");
    this->InsertPrjName( tableName, 3, "Blasius Incompressible Laminar Flat Plate");
    this->InsertPrjName( tableName, 4, "Driver-Seegmiller Incompressible Backward-Facing Step");
    this->InsertPrjName( tableName, 5, "Fraser Subsonic Conical Diffuser");
    this->InsertPrjName( tableName, 6, "Incompressible, Buice Axisymmetric Diffuser");
    this->InsertPrjName( tableName, 7, "NLR Airfoil with Flap");
    this->InsertPrjName( tableName, 8, "Incompressible, Turbulent Flat Plate");
    this->QueryTable( tableName );
    qDebug() << "---- insert operation end----- " << "\n ";

    qDebug() << "---- query operation start----";
    this->QueryPjrRecordByName( tableName, "Sod''s Shock Tube");
    this->QueryPjrRecordByName( tableName, "Blasius Incompressible Laminar Flat Plate");
    this->QueryPjrRecordByName( tableName, "Incompressible, Turbulent Flat Plate");
    qDebug() << "---- query operation end----- " << "\n ";

    qDebug() << "---- modify operation start----- ";
    this->ModifyPjrRecord( tableName, 1, "Shock Tube" );
    this->ModifyPjrRecord( tableName, 2, "Driven Cavity" );
    this->QueryTable( tableName );
    qDebug() << "---- modify operation end----- " << "\n ";

    qDebug() << "---- delete operation start----";
    this->DeletePjrRecord( tableName, "Fraser Subsonic Conical Diffuser");
    this->DeletePjrRecord( tableName, "NLR Airfoil with Flap");
    this->QueryTable( tableName );
    qDebug() << "---- delete operation end----" << "\n ";

}

void MyDataBase::InsertPrjName( const QString & tableName, const int &prjId, const QString &prjName )
{
    QSqlQuery query( this->db );
    QString queryString;
    queryString.append( "INSERT INTO " );
    queryString.append( tableName );
    queryString.append( "(id, prjname ) VALUES (:id, :prjname)" );

    query.prepare( queryString );
    query.bindValue( ":id", prjId );
    query.bindValue( ":prjname", prjName );
    query.exec();
}

void MyDataBase::DeletePjrRecord( const QString & tableName, const QString &prjName )
{
    QSqlQuery query( this->db );
    QString queryString;
    queryString.append( "DELETE FROM " );
    queryString.append( tableName );
    queryString.append( " WHERE prjname=:prjname" );

    query.prepare( queryString );
    query.bindValue(":prjname", prjName);
    query.exec();
}

void MyDataBase::ModifyPjrRecord( const QString & tableName, const int &prjId, const QString &prjName )
{
    QSqlQuery query( this->db );

    QString queryString;
    queryString.append( "update " );
    queryString.append( tableName );
    queryString.append( " set prjname=:prjname WHERE id=:id" );

    query.prepare( queryString );
    query.bindValue(":id", prjId);
    query.bindValue(":prjname", prjName);
    query.exec();
}

void MyDataBase::QueryPjrRecordByName( const QString & tableName, const QString &prjName )
{
    QSqlQuery query( this->db );

    QString queryString;
    queryString.append( "SELECT * FROM " );
    queryString.append( tableName );
    queryString.append( " WHERE prjname='" );
    queryString.append( prjName );
    queryString.append( "'" );

    query.exec( queryString );

    while ( query.next() )
    {
        qDebug() << QString("Id: %1, Prjname: %2")
            .arg(query.value("id").toInt())
            .arg(query.value("prjname").toString());
    }
}

void MyDataBase::QueryTable( const QString & tableName )
{
    QString sql;
    sql.append( "SELECT id, prjname FROM " );
    sql.append( tableName );

    QSqlQuery query( db );
    query.exec( sql );
    while ( query.next() )
    {
        qDebug() << QString("Id: %1, Prjname: %2")
            .arg( query.value("id").toInt() )
            .arg( query.value("prjname").toString() );
    }
}
