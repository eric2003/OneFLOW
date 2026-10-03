#include "GridLayoutParser.h"

#include "GridLayout.h"
#include "TextFileParser.h"

#include <stdexcept>
#include <string>

BeginNameSpace( ONEFLOW )

namespace
{

GridDistributionType ParseDistributionType( const std::string & type )
{
    if ( type == "r" ) return GridDistributionType::Ratio;
    if ( type == "d" ) return GridDistributionType::Distance;
    if ( type == "tanh" ) return GridDistributionType::Tanh;
    if ( ! type.empty() && type.front() == 'c' ) return GridDistributionType::Copy;
    if ( ! type.empty() && type.front() == 'e' ) return GridDistributionType::Exponential;

    throw std::runtime_error( "Unsupported grid distribution: " + type );
}

}

GridLayout GridLayoutParser::Parse( const std::string & fileName ) const
{
    GridLayout layout;

    const std::string separator = " =\r\n\t#$,;\"(){}";
    TextFileParser parser;
    parser.OpenPrjFile( fileName, std::ios_base::in );
    parser.SetDefaultSeparator( separator );

    while ( ! parser.ReachTheEndOfFile() )
    {
        if ( ! parser.ReadNextMeaningfulLine() ) break;

        const std::string keyword = parser.ReadNextWord();

        if ( keyword == "Point" )
        {
            GridPointDefinition point;
            point.id = parser.ReadNextDigit< int >();
            point.x = parser.ReadNextDigit< Real >();
            point.y = parser.ReadNextDigit< Real >();
            point.z = parser.ReadNextDigit< Real >();
            layout.points.push_back( point );
        }
        else if ( keyword == "Line" )
        {
            GridLineDefinition line;
            line.id = parser.ReadNextDigit< int >();
            line.p1 = parser.ReadNextDigit< int >();
            line.p2 = parser.ReadNextDigit< int >();
            layout.lines.push_back( line );
        }
        else if ( keyword == "Circle" )
        {
            GridCircleDefinition circle;
            circle.id = parser.ReadNextDigit< int >();
            circle.p1 = parser.ReadNextDigit< int >();
            circle.pc = parser.ReadNextDigit< int >();
            circle.p2 = parser.ReadNextDigit< int >();
            layout.circles.push_back( circle );
        }
        else if ( keyword == "Dim" )
        {
            GridDimensionDefinition dimension;
            dimension.id = parser.ReadNextDigit< int >();
            dimension.pointCount = parser.ReadNextDigit< int >();
            layout.dimensions.push_back( dimension );
        }
        else if ( keyword == "Ds" )
        {
            GridDistributionDefinition distribution;
            distribution.lineId = parser.ReadNextDigit< int >();

            const std::string type = parser.ReadNextWord();
            distribution.type = ParseDistributionType( type );

            if ( distribution.type == GridDistributionType::Copy )
            {
                while ( true )
                {
                    const std::string word = parser.ReadNextWord();
                    if ( word.empty() ) break;
                    distribution.copyLineIds.push_back( StringToDigit< int >( word ) );
                }
            }
            else if ( distribution.type == GridDistributionType::Ratio ||
                      distribution.type == GridDistributionType::Distance ||
                      distribution.type == GridDistributionType::Tanh )
            {
                distribution.startValue = parser.ReadNextDigit< Real >();
                distribution.endValue = parser.ReadNextDigit< Real >();
            }

            layout.distributions.push_back( distribution );
        }
        else if ( keyword == "Boundary" )
        {
            GridBoundaryDefinition boundary;
            boundary.id = parser.ReadNextDigit< int >();
            boundary.boundaryType = parser.ReadNextDigit< int >();
            layout.boundaries.push_back( boundary );
        }
        else if ( keyword == "Add" )
        {
            const std::string relationType = parser.ReadNextWord();
            if ( relationType == "L2F" )
            {
                GridLineToFaceDefinition relation;
                relation.faceId = parser.ReadNextDigit< int >();
                relation.position = parser.ReadNextDigit< int >();
                relation.lineId = parser.ReadNextDigit< int >();
                layout.lineToFaces.push_back( relation );
            }
            else if ( relationType == "F2B" )
            {
                GridFaceToBlockDefinition relation;
                relation.blockId = parser.ReadNextDigit< int >();
                relation.position = parser.ReadNextDigit< int >();
                relation.faceId = parser.ReadNextDigit< int >();
                layout.faceToBlocks.push_back( relation );
            }
            else
            {
                throw std::runtime_error( "Unsupported grid relation: " + relationType );
            }
        }
    }

    parser.CloseFile();
    layout.Validate();
    return layout;
}

EndNameSpace
