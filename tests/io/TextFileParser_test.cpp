#include "gtest/gtest.h"
#include "TextFileParser.h"
#include <fstream>
#include <cstdio>

using namespace ONEFLOW;

namespace
{
    // Writes a temp file with the given content and returns its path.
    std::string WriteTempFile( const std::string & name, const std::string & content )
    {
        std::ofstream out( name );
        out << content;
        out.close();
        return name;
    }
}

// --- Reproduces the exact read pattern MRegister::Register() uses on
// system/*/funcMap.txt and taskMap.txt: "NAME, ClassName, nParams[, extraWord]*"
// per non-empty/non-comment line. ---
TEST( SystemConfigRead, ParsesFuncMapStyleFile )
{
    WriteTempFile( "funcMap_test.txt",
        "\n"
        "FILL_WALL_STRUCT  , CFillWallStruct, 0\n"
        "CALC_WALL_DIST    , CCalcWallDist, 0\n"
        "READ_WALL_DIST    , CReadWallDist, 0\n" );

    TextFileParser parser;
    parser.OpenFile( "funcMap_test.txt", std::ios_base::in );
    parser.SetDefaultSeparator( " =\r\n\t#$,;\"()" );

    std::vector< std::pair< std::string, std::string > > entries;

    while ( ! parser.ReachTheEndOfFile() )
    {
        if ( ! parser.ReadNextNonEmptyLine() ) break;

        std::string actionName = parser.ReadNextWord();
        std::string className  = parser.ReadNextWord();
        int nParameters = parser.ReadNextDigit< int >();

        EXPECT_EQ( nParameters, 0 );
        entries.push_back( { actionName, className } );
    }

    parser.CloseFile();
    std::remove( "funcMap_test.txt" );

    ASSERT_EQ( entries.size(), static_cast< size_t >( 3 ) );
    EXPECT_EQ( entries[ 0 ].first,  "FILL_WALL_STRUCT" );
    EXPECT_EQ( entries[ 0 ].second, "CFillWallStruct" );
    EXPECT_EQ( entries[ 2 ].first,  "READ_WALL_DIST" );
    EXPECT_EQ( entries[ 2 ].second, "CReadWallDist" );
}

// --- fileMap.txt has a nonzero parameter count followed by that many
// trailing words -- exercises the multi-word-per-entry path. ---
TEST( SystemConfigRead, ParsesFileMapStyleFileWithExtraParameters )
{
    WriteTempFile( "fileMap_test.txt",
        "READ_WALL_DIST , CSetFile, 3, walldist_file, ( in , binary )\n"
        "WRITE_WALL_DIST, CSetFile, 4, walldist_file, ( out, binary, trunc )\n" );

    TextFileParser parser;
    parser.OpenFile( "fileMap_test.txt", std::ios_base::in );
    parser.SetDefaultSeparator( " =\r\n\t#$,;\"()" );

    struct Entry { std::string action, cls; std::vector< std::string > params; };
    std::vector< Entry > entries;

    while ( ! parser.ReachTheEndOfFile() )
    {
        if ( ! parser.ReadNextNonEmptyLine() ) break;

        Entry e;
        e.action = parser.ReadNextWord();
        e.cls    = parser.ReadNextWord();
        int nParameters = parser.ReadNextDigit< int >();
        for ( int i = 0; i < nParameters; ++ i )
        {
            e.params.push_back( parser.ReadNextWord() );
        }
        entries.push_back( e );
    }

    parser.CloseFile();
    std::remove( "fileMap_test.txt" );

    ASSERT_EQ( entries.size(), static_cast< size_t >( 2 ) );
    EXPECT_EQ( entries[ 0 ].action, "READ_WALL_DIST" );
    ASSERT_EQ( entries[ 0 ].params.size(), static_cast< size_t >( 3 ) );
    EXPECT_EQ( entries[ 0 ].params[ 0 ], "walldist_file" );
    EXPECT_EQ( entries[ 0 ].params[ 1 ], "in" );
    EXPECT_EQ( entries[ 0 ].params[ 2 ], "binary" );

    EXPECT_EQ( entries[ 1 ].action, "WRITE_WALL_DIST" );
    ASSERT_EQ( entries[ 1 ].params.size(), static_cast< size_t >( 4 ) );
    EXPECT_EQ( entries[ 1 ].params[ 3 ], "trunc" );
}

// --- One name per line, blank lines interleaved (actionName.txt shape). ---
TEST( SystemConfigRead, ParsesOneNamePerLineWithBlankLines )
{
    WriteTempFile( "actionName_test.txt",
        "NO_TASK\n"
        "READ_INSRESTART\n"
        "\n"
        "READ_RESTART\n"
        "DUMP_RESTART\n" );

    TextFileParser parser;
    parser.OpenFile( "actionName_test.txt", std::ios_base::in );
    parser.SetDefaultSeparator( " =\r\n\t#$,;\"" );

    std::vector< std::string > names;
    while ( parser.ReadNextNonEmptyLine() )
    {
        names.push_back( parser.ReadNextWord() );
        if ( parser.ReachTheEndOfFile() ) break;
    }

    parser.CloseFile();
    std::remove( "actionName_test.txt" );

    ASSERT_EQ( names.size(), static_cast< size_t >( 4 ) );
    EXPECT_EQ( names[ 0 ], "NO_TASK" );
    EXPECT_EQ( names[ 3 ], "DUMP_RESTART" );
}

// --- The "N*value" shorthand used by numeric config blocks elsewhere in
// system/*.txt (e.g. "5*0.0"). ---
TEST( SystemConfigRead, ReadsRepeatCountShorthand )
{
    WriteTempFile( "repeat_test.txt", "5*0.0 1.5\n" );

    TextFileParser parser;
    parser.OpenFile( "repeat_test.txt", std::ios_base::in );
    parser.SetDefaultSeparator( " =\r\n\t#$,;\"" );
    parser.ReadNextNonEmptyLine();

    int repeatCount = 0;
    double value = parser.ReadNextDigit< double >( repeatCount );
    EXPECT_EQ( repeatCount, 5 );
    EXPECT_EQ( value, 0.0 );

    double second = parser.ReadNextDigit< double >();
    EXPECT_EQ( second, 1.5 );

    parser.CloseFile();
    std::remove( "repeat_test.txt" );
}

// --- Comment lines ('#' and '//') must be skipped by ReadNextMeaningfulLine,
// used by SkipReadSymbol/SkipReadWholeBlock. ---
TEST( SystemConfigRead, SkipsCommentLinesAndFindsBlock )
{
    WriteTempFile( "block_test.txt",
        "# this whole file is a comment-laden block\n"
        "// another comment style\n"
        "SECTION\n"
        "{\n"
        "  KEY VALUE\n"
        "}\n" );

    TextFileParser parser;
    parser.OpenFile( "block_test.txt", std::ios_base::in );
    parser.SetDefaultSeparator( " =\r\n\t#$,;\"{}" );

    parser.SkipReadSymbol( "SECTION" );
    EXPECT_FALSE( parser.ReachTheEndOfFile() );

    parser.CloseFile();
    std::remove( "block_test.txt" );
}

