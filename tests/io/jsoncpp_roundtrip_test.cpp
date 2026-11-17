#include <gtest/gtest.h>

#include "json/json.h"

#include <sstream>
#include <string>

TEST( JsonCppTest, ParsesNestedDocument )
{
    const std::string input =
        R"({"name":"OneFLOW","version":2,"mesh":{"dimension":2},"steps":[1,2,3]})";
    std::istringstream stream( input );
    Json::CharReaderBuilder reader;
    Json::Value root;
    JSONCPP_STRING errors;

    ASSERT_TRUE( Json::parseFromStream( reader, stream, &root, &errors ) )
        << errors;
    EXPECT_EQ( root[ "name" ].asString(), "OneFLOW" );
    EXPECT_EQ( root[ "version" ].asInt(), 2 );
    EXPECT_EQ( root[ "mesh" ][ "dimension" ].asInt(), 2 );
    ASSERT_EQ( root[ "steps" ].size(), 3u );
    EXPECT_EQ( root[ "steps" ][ 2 ].asInt(), 3 );
}

TEST( JsonCppTest, ReportsMalformedDocument )
{
    const std::string input = R"({"name":"OneFLOW", "steps":[1,2,})";
    std::istringstream stream( input );
    Json::CharReaderBuilder reader;
    Json::Value root;
    JSONCPP_STRING errors;

    EXPECT_FALSE( Json::parseFromStream( reader, stream, &root, &errors ) );
    EXPECT_FALSE( errors.empty() );
}

TEST( JsonCppTest, SerializedDocumentCanBeParsedBack )
{
    Json::Value expected;
    expected[ "task" ] = "Solve";
    expected[ "dimension" ] = 2;
    expected[ "enabled" ] = true;
    expected[ "values" ].append( 0.25 );
    expected[ "values" ].append( 0.5 );

    Json::StreamWriterBuilder writer;
    writer[ "indentation" ] = "";
    const std::string serialized = Json::writeString( writer, expected );

    std::istringstream stream( serialized );
    Json::CharReaderBuilder reader;
    Json::Value actual;
    JSONCPP_STRING errors;

    ASSERT_TRUE( Json::parseFromStream( reader, stream, &actual, &errors ) )
        << errors;
    EXPECT_EQ( actual, expected );
}
