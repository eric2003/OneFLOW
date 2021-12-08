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
#include "JsonTest.h"
#include "json/json.h"
#include <iostream>
#include <fstream>


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
    //std::string json_file = writer.write(root);
    //string json_file = writer.write(root, &std::cout);
    //std::cout << "demo write json ==============\n";
    //std::cout << json_file << std::endl;
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

    std::ifstream f;
    f.open( "test.json", std::ios::in );
    if ( ! f.is_open() )
    {
        std::cout << "Open json file error!" << std::endl;
    }

    JSONCPP_STRING errs;

    bool parse_ok = Json::parseFromStream(reader, f, &root, &errs);

    std::cout << root.size() << std::endl;
    std::string a1 = root[ "name" ].asString();
    std::string a2 = root[ "age" ].asString();
    std::string a3 = root[ "sex_is_male" ].asString();
    Json::Value & v = root[ "partner" ];
    std::string b1 = v[ "partner_name" ].asString();
    std::string b2 = v[ "partner_age" ].asString();
    std::string b3 = v[ "partner_sex_is_male" ].asString();

    Json::Value & w = root[ "achievement" ];
    int s = w.size();
    for ( int i = 0; i < s; ++ i )
    {
        std::string ss = w[ i ].asString();
        std::cout << " ss = " << ss << "\n";
    }

    std::ofstream os;
    os.open("1.json", std::ios::out);
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
