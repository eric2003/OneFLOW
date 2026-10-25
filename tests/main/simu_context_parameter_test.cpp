#include <gtest/gtest.h>

#include "DataBaseType.h"
#include "DataObject.h"
#include "DataPara.h"
#include "SimuContextParameter.h"

#include <memory>
#include <stdexcept>
#include <string>

using namespace ONEFLOW;

namespace
{

void AddStringParameter(
    DataPara & parameters, const std::string & name, const std::string & value )
{
    auto dataObject = std::make_unique< TDataObject< std::string > >( 1 );
    dataObject->CopyValue( &value, 1 );
    parameters.SetDataEntry( std::make_unique< DataEntry >(
        name, HX_STRING, 1, std::move( dataObject ) ) );
}

void AddIntegerParameter(
    DataPara & parameters, const std::string & name, int value )
{
    auto dataObject = std::make_unique< TDataObject< int > >( 1 );
    dataObject->CopyValue( &value, 1 );
    parameters.SetDataEntry( std::make_unique< DataEntry >(
        name, HX_INT, 1, std::move( dataObject ) ) );
}

}

TEST( SimuContextParameterTest, MissingOptionalStringReturnsEmpty )
{
    DataPara parameters;

    EXPECT_EQ(
        SimuContextDetail::ReadOptionalStringParameter(
            parameters, "inputRoot" ),
        "" );
}

TEST( SimuContextParameterTest, ReadsStringParameter )
{
    DataPara parameters;
    AddStringParameter( parameters, "inputRoot", "mesh/input" );

    EXPECT_EQ(
        SimuContextDetail::ReadOptionalStringParameter(
            parameters, "inputRoot" ),
        "mesh/input" );
}

TEST( SimuContextParameterTest, RejectsNonStringAndNamesParameter )
{
    DataPara parameters;
    AddIntegerParameter( parameters, "inputRoot", 12 );

    try
    {
        SimuContextDetail::ReadOptionalStringParameter(
            parameters, "inputRoot" );
        FAIL() << "Expected a type validation exception";
    }
    catch ( const std::invalid_argument & error )
    {
        EXPECT_NE(
            std::string( error.what() ).find( "inputRoot" ),
            std::string::npos );
        EXPECT_NE(
            std::string( error.what() ).find( "must be a string" ),
            std::string::npos );
    }
}
