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
#include "JsonTest.h"
#include "json/json.h"
#include <iostream>
#include <fstream>
using namespace std;

BeginNameSpace( ONEFLOW )

void demo_write_array();

void demo_write_array()
{
    Json::Value root;
    Json::StreamWriterBuilder writer;
    {
        Json::Value person; 
        person["name"] = "Tocy";
        person["salary"] = 100;
        root[0] = person;
    }

    {
        Json::Value person; 
        person["name"] = "Kit";
        person["salary"] = 89;
        root[1] = person;
    }

    root[2] = "a json note";   
    //string json_file = writer.write(root);
    //string json_file = writer.write(root, &std::cout);
    //cout << "demo write json ==============\n";
    //cout << json_file << endl;
}

void test_demo_write_array()
{
    Json::Value root;
    Json::Value data;
    constexpr bool shouldUseOldWay = false;
    root[ "action" ] = "run";
    data[ "number" ] = 1;
    root[ "data" ] = data;

    Json::StreamWriterBuilder builder;
    const std::string json_file = Json::writeString( builder, root );
    std::cout << json_file << std::endl;
}

void test_demo_write();
void test_demo_write()
{
    Json::Value root;
    Json::Value data;
    root[ "action" ] = "run";
    data[ "number" ] = 1;
    root[ "data" ] = data;

    Json::StreamWriterBuilder writer;
    {
        Json::Value person; 
        person["name"] = "Tocy";
        person["salary"] = 100;
        root[0] = person;
    }

    {
        Json::Value person; 
        person["name"] = "Kit";
        person["salary"] = 89;
        root[1] = person;
    }


    Json::StreamWriterBuilder builder;
    const std::string json_file = Json::writeString( builder, root );
    std::cout << json_file << std::endl;

}

void readFileJson()
{
    Json::CharReaderBuilder reader;
    Json::Value root;

    ifstream f;
    f.open( "test.json", ios::in );
    if ( ! f.is_open() )
    {
        cout << "Open json file error!" << endl;
    }

    JSONCPP_STRING errs;

    bool parse_ok = Json::parseFromStream(reader, f, &root, &errs);

    cout << root.size() << endl;
    string a1 = root[ "name" ].asString();
    string a2 = root[ "age" ].asString();
    string a3 = root[ "sex_is_male" ].asString();
    Json::Value & v = root[ "partner" ];
    string b1 = v[ "partner_name" ].asString();
    string b2 = v[ "partner_age" ].asString();
    string b3 = v[ "partner_sex_is_male" ].asString();

    Json::Value & w = root[ "achievement" ];
    int s = w.size();
    for ( int i = 0; i < s; ++ i )
    {
        string ss = w[ i ].asString();
        cout << " ss = " << ss << "\n";
    }

    ofstream os;
    os.open("1.json", ios::out);
    Json::StreamWriterBuilder builder;
    std::unique_ptr<Json::StreamWriter> writer( builder.newStreamWriter() );
    writer->write(root, &os);

    int kkk = 1;
}

JsonTest::JsonTest()
{
    ;
}

JsonTest::~JsonTest()
{
}

void JsonTest::Run()
{
    readFileJson();
}

EndNameSpace
