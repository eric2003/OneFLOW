#include "Solver.h"
#include <fstream>

double compute_l2norm( int ni, std::vector<double> & r )
{
    double rms = 0.0;
    for ( int i = 1; i < ni - 1; ++ i )
    {
        rms += r[ i ] * r[ i ];
    }
    rms = std::sqrt( rms / ( ni - 2 ) );
    return rms;
}

double compute_max_error( int ni, std::vector<double> & u_error )
{
    double val_max = -1;
    int ipos = -1;
    for ( int i = 1; i < ni - 1; ++ i )
    {
        //val_max = std::max( val_max, std::abs( u_error[ i ] ) );
        if ( val_max < std::abs( u_error[ i ] ) )
        {
            ipos = i;
            val_max = std::abs( u_error[ i ] );
        }
    }
    std::cout << " ipos = " << ipos << "\n";
    return val_max;
}

void DumpErrorDetails( std::vector<double> &u_error )
{
    int ni = u_error.size();
    double rms_error = compute_l2norm( ni, u_error );
    double max_error = compute_max_error( ni, u_error );
    std::cout << "max_error = " << std::setprecision(15) << max_error << "\n";
    //create output file for L2-norm
    std::fstream file;
    file.open("output.txt", std::fstream::out);
    std::format_to(std::ostream_iterator<char>(file), "Error details: \n");
    std::format_to(std::ostream_iterator<char>(file), "L-2 Norm = {0}\n", rms_error);
    std::format_to(std::ostream_iterator<char>(file), "Maximum Norm = {0}\n", max_error);
    file.close();
}

void DumpCsvFile( const std::string &filename, std::vector<double> &x, std::vector<double> &ue, std::vector<double> &un, std::vector<double> &uerror )
{
    std::fstream file;
    file.open(filename.c_str(), std::fstream::out);
    std::format_to(std::ostream_iterator<char>(file), "x ue un uerror\n");
    for ( int i = 0; i < x.size(); ++ i )
    {
        std::format_to(std::ostream_iterator<char>(file), "{:.16f} {:.16f} {:.16f} {:.16f}\n", x[i], ue[i], un[i], uerror[i] );
    }
    file.close();
}