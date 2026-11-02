#include "post.h"
#include "Vec1d.h"
#include <fstream>
#include <format>

void DumpCsvFile( const std::string &filename, std::vector<double> &x, std::vector<std::vector<double>> &u )
{
    std::fstream file;
    file.open( filename.c_str(), std::fstream::out );
    for ( int i = 0; i < u.size(); ++ i )
    {
        std::string str = {};
        str += std::format( "{:.16f} ", x[ i ] );
        int nj = u[ i ].size();
        for ( int j = 0; j < nj; ++ j )
        {
            double value = u[ i ][ j ];
            str += std::format( "{:.16f}", value );
            if ( j != nj - 1 )
            {
                str += " ";
            }
        }
        std::format_to(std::ostream_iterator<char>(file), "{}\n", str );
    }
    file.close();
}

void DumpField( int iter, Vec1d & x, Vec1d & u )
{
    std::string total_string = {};
    std::fstream file;
    {
        std::string filename = std::format( "field_final{}.csv", iter );
        file.open( filename.c_str(), std::fstream::out );
    }
    std::string file_string = {};

    for ( int i = 0; i < x.size(); ++ i )
    {
        file_string += std::format( "{:.16f} {:.16f}\n", x[ i ], u[ i ] );
    }

    total_string += file_string;

    std::format_to( std::ostream_iterator<char>( file ), "{}", total_string );
    file.close();
}

