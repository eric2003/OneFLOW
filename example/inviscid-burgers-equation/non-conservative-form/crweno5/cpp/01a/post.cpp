#include "post.h"
#include <fstream>

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
