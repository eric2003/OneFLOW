#include <iostream>
#include <fstream>
#include <format>
#include <chrono>
#include <Eigen/Dense>

using Eigen::MatrixXd;
using Eigen::VectorXd;

double compute_l2norm( int nx, int ny, MatrixXd &r )
{
    double rms = 0.0;
    for ( int j = 0; j < ny + 1; ++ j )
	{
		for ( int i = 0; i < nx + 1; ++ i )
		{
			rms += r(i, j) * r(i, j);
		}
	}
    rms = sqrt( rms / ( ( nx + 1 ) * ( ny + 1 ) ) );
    return rms;
}

double get_max_value( int nx, int ny, MatrixXd & field )
{
	double max_value =  abs( field( 0, 0 ) );
	for ( int j = 0; j < ny + 1; ++ j )
	{
		for ( int i = 0; i < nx + 1; ++ i )
		{
			max_value = std::max( max_value, abs( field(i, j) ) );
		}
	}
	return max_value;
}

void compute_residual( int nx, int ny, double dx, double dy, MatrixXd &f, MatrixXd &un, MatrixXd &r )
{
	double dx2 = dx * dx;
	double dy2 = dy * dy;

	for ( int j = 1; j < ny; ++ j )
	{
		for ( int i = 1; i < nx; ++ i )
		{
			double d2udx2 = ( un( i - 1, j ) - 2 * un( i, j ) + un( i + 1, j ) ) / ( dx2 );
			double d2udy2 = ( un( i, j - 1 ) - 2 * un( i, j ) + un( i, j + 1 ) ) / ( dy2 );
			r(i, j) = f(i, j) - d2udx2 - d2udy2;
		}
	}
}

void gauss_seidel( double dx, double dy, int nx, int ny, MatrixXd &r, MatrixXd &f, MatrixXd &un, double rms, double init_rms, int max_iter, double tolerance, std::fstream &outfile )
{
	// create text file for writing residual history
	std::fstream residual_plot;
	residual_plot.open( "gs_residual.txt", std::ios_base::out );

	double count = 0.0;

	compute_residual( nx, ny, dx, dy, f, un, r );

	rms = compute_l2norm( nx, ny, r );

	init_rms = rms;
	int iter_count = 0;
	std::cout << iter_count << " " << init_rms << " " << rms / init_rms << "\n";

	double dx2 = dx * dx;
	double dy2 = dy * dy;

	double den = -2.0 / dx2 - 2.0 / dy2;
	for ( iter_count = 1; iter_count < 5 * max_iter; ++ iter_count )
	{
		// compute solution at next time step
		// residual
		for ( int j = 1; j < ny; ++ j )
		{
			for ( int i = 1; i < nx; ++ i )
			{
				double d2udx2 = ( un( i - 1, j ) - 2 * un( i, j ) + un( i + 1, j ) ) / ( dx2 );
				double d2udy2 = ( un( i, j - 1 ) - 2 * un( i, j ) + un( i, j + 1 ) ) / ( dy2 );
				r( i, j ) = f( i, j ) - d2udx2 - d2udy2;

				un( i, j ) = un( i, j ) + r( i, j ) / den;
			}
		}
		compute_residual( nx, ny, dx, dy, f, un, r );

		//compute the l2norm of residual;
		rms = compute_l2norm( nx, ny, r );

		std::format_to(std::ostream_iterator<char>(residual_plot), "{0} {1} {2}\n", iter_count, rms, rms / init_rms );

		count = iter_count;

		std::cout << iter_count <<  " " << rms <<  " "<<  rms / init_rms << "\n";
		if ( ( rms / init_rms ) <= tolerance )
		{
			break;
		}
	}

	double max_r = get_max_value( nx, ny, r );

	std::format_to(std::ostream_iterator<char>(outfile), "L-2 Norm = {0}\n", rms);
	std::format_to(std::ostream_iterator<char>(outfile), "Maximum Norm = {0}\n", max_r);
	std::format_to(std::ostream_iterator<char>(outfile), "Iterations =  {0}\n", count);

	residual_plot.close();
}

void dumpfield( VectorXd&x, VectorXd&y,  MatrixXd &f,  MatrixXd &un,  MatrixXd &ue, int nx, int ny, std::fstream &outfile )
{
	std::format_to( std::ostream_iterator<char>( outfile ), "{0} {1}\n", nx, ny );
	for ( int j = 0; j < ny + 1; ++ j )
	{
		for ( int i = 0; i < nx + 1; ++ i )
		{
			std::format_to(std::ostream_iterator<char>(outfile), "{0} {1} {2} {3} {4}\n",
				x(i), y(j), f(i,j), un(i,j), ue(i,j) );
		}
	}
}
 
int main( int argc, char **argv )
{
	int nx = 512;
	int ny = 512;

	double tolerance = 1.0e-10;
	int	max_iter = 20 * 100000;

	std::cout << "max_iter = " << max_iter << "\n";

	double x_l = 0.0;
	double x_r = 1.0;
	double y_b = 0.0;
	double y_t = 1.0;

	double	dx = ( x_r - x_l ) / nx;
	double	dy = ( y_t - y_b ) / ny;

	VectorXd x( nx + 1 );
	VectorXd y( nx + 1 );
	MatrixXd ue = MatrixXd::Zero( nx + 1, ny + 1 );
	MatrixXd f = MatrixXd::Zero( nx + 1, ny + 1 );
	MatrixXd un = MatrixXd::Zero( nx + 1, ny + 1 );

	for ( int i = 0; i < nx + 1; ++ i )
	{
		x(i) = x_l + dx * ( i );
	}

	for ( int j = 0; j < ny + 1; ++ j )
	{
		y(j) = y_b + dy * ( j );
	}

	//given exact solution
	for ( int j = 0; j < ny + 1; ++ j )
	{
		for ( int i = 0; i < nx + 1; ++ i )
		{
			ue(i, j) = ( x(i) * x(i) - 1.0 ) * ( y(j) * y(j) - 1.0 );
			f(i, j) = -2.0 * ( 2.0 - x(i) * x(i) - y(j) * y(j) );
			un(i, j) = 0.0;
		}
	}

	for ( int i = 0; i < nx + 1; ++ i )
	{
		un( i, 0 ) = ue( i, 0 );
		un( i, ny ) = ue( i, ny );
	}

	for ( int j = 0; j < ny + 1; ++ j )
	{
		un( 0, j ) = ue( 0, j );
		un( nx, j ) = ue( nx, j );
	}

	MatrixXd r = MatrixXd::Zero( nx + 1, ny + 1 );

	double init_rms = 0.0;
	double rms = 0.0;

	//create text file for initial field
	std::fstream field_initial;
	field_initial.open( "field_initial.txt", std::ios_base::out );
	dumpfield( x, y, f, un, ue, nx, ny, field_initial );
	field_initial.close();

	//create output file for L2-norm
	std::fstream outfile;
	outfile.open( "output.txt", std::ios_base::out );
	outfile << "Residual details: \n";

	using clock_type = std::chrono::time_point<std::chrono::system_clock>;
	clock_type time_now, time_old;

	time_old = std::chrono::system_clock::now();

	gauss_seidel( dx, dy, nx, ny, r, f, un, rms, init_rms, max_iter, tolerance, outfile );

	MatrixXd uerror( nx + 1, ny + 1 );

	double rms_error = 0.0;

	uerror = un - ue;

	rms_error = compute_l2norm( nx, ny, uerror );
	double max_error = get_max_value( nx, ny, uerror );

	time_now = std::chrono::system_clock::now();
	double elapsed_ms = std::chrono::duration_cast<std::chrono::milliseconds>
		(time_now-time_old).count();


	double t = elapsed_ms/1000.0;

    std::cout << "Error details:\n";
	std::cout << "L-2 Norm = " << rms_error << "\n";
	std::cout << "Maximum Norm = " << max_error << "\n";
	std::cout << "CPU Time = " << t << " s\n";

	outfile << "Error details: \n";
	std::format_to(std::ostream_iterator<char>(outfile), "L-2 Norm = {0}\n", rms_error);
	std::format_to(std::ostream_iterator<char>(outfile), "Maximum Norm = {0}\n", max_error);
	std::format_to(std::ostream_iterator<char>(outfile), "CPU Time = {0} s\n", t);
	outfile.close();

	//create text file for final field
	std::fstream field_final;
	field_final.open( "field_final.txt", std::ios_base::out );
	dumpfield( x, y, f, un, ue, nx, ny, field_final );
	field_final.close();

    return 0;
}
