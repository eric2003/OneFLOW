#include <iostream>
#include <fstream>
#include <format>
#include <chrono>
#include <numbers>
#include <cmath>
#include <Eigen/Dense>

using Eigen::MatrixXd;
using Eigen::VectorXd;

void boundary( int nx, int ny, MatrixXd & f );
void bc( int nx, int ny, double dx, double dy, MatrixXd & w, MatrixXd & s );
void bc2( int nx, int ny, double dx, double dy, MatrixXd & w, MatrixXd & s );
void cal_stream_function( int nx, int ny, double dx, double dy, MatrixXd & sf, MatrixXd & vt );
void dumpfield( VectorXd & x, VectorXd & y, MatrixXd & field, int nx, int ny, double re, const std::string & filename );
void rhs( int nx, int ny, double dx, double dy, double re, MatrixXd & w, MatrixXd & s, MatrixXd & r );
void rhsnew( int nx, int ny, double dx, double dy, double re, MatrixXd & w, MatrixXd & r );
void numerical( int nx, int ny, int nt, int ns, VectorXd & x, VectorXd & y, double dx, double dy, double dt, double re, MatrixXd & wn );
void numericalnew( int nx, int ny, int nt, int ns, VectorXd & x, VectorXd & y, double dx, double dy, double dt, double re, MatrixXd & wn );
void vm_ic( int nx, int ny, VectorXd & x, VectorXd & y, MatrixXd & w );

void bc(int nx, int ny, double dx, double dy, MatrixXd &w, MatrixXd &s)
{
    //first order approximation
    //boundary condition for vorticity (Hoffmann) left and right
    for ( int j = 0; j < ny + 1; ++ j )
    {
        w(0 ,j) = -2.0*s(1   ,j)/(dx*dx);
        w(nx,j) = -2.0*s(nx-1,j)/(dx*dx);
    }

    //boundary condition for vorticity (Hoffmann) bottom and top
    for ( int i = 0; i < nx + 1; ++ i )
    {
        w(i,0 ) = -2.0*s(i,1   )/(dy*dy);
        w(i,ny) = -2.0*s(i,ny-1)/(dy*dy) - 2.0/dy;
    }
}

void bc2(int nx, int ny, double dx, double dy, MatrixXd &w, MatrixXd &s)
{
    //second order approximation
    //boundary condition for vorticity (Jensen) left and right
    for ( int j = 0; j < ny + 1; ++ j )
	{
        w(0 ,j) = (-4.0*s(1   ,j)+0.5*s(2   ,j))/(dx*dx);
        w(nx,j) = (-4.0*s(nx-1,j)+0.5*s(nx-2,j))/(dx*dx);
	}

    //boundary condition for vorticity (Jensen) bottom and top
    for ( int i = 0; i < nx + 1; ++ i )
	{
        w(i,0 ) = (-4.0*s(i,1   )+0.5*s(i,2   ))/(dy*dy);
        w(i,ny) = (-4.0*s(i,ny-1)+0.5*s(i,ny-2))/(dy*dy) - 3.0/dy;
	}
}

void cal_stream_functionBAK(int nx, int ny, double dx, double dy, MatrixXd &sf, MatrixXd &vt)
{
    //int    maxIt  = 100;
    int    maxIt  = 300;
    double beta   = 1.5;
    double maxErr = 1.0e-6;

    double b = dx / dy;
    double b2 = b * b;
    double dx2 = dx * dx;
    double cc = 1.0 /( 2.0 * ( 1.0 + b2 ) );
    double err0 = 1.0;
	MatrixXd w = MatrixXd::Zero( nx + 1, ny + 1 );
	for ( int iter = 0; iter < maxIt; ++ iter )
	{
		//MatrixXd w = sf; //by SOR iteration
        for ( int j = 0; j < ny + 1; ++ j )
        {
            for ( int i = 0; i < nx + 1; ++ i )
            {
                w(i,j) = sf(i,j);
            }
        }

		for ( int j = 1; j < ny; ++ j )
		{
			for ( int i = 1; i < nx; ++ i )
			{
                double sfnew = cc * ( dx2 * vt(i,j) + sf(i+1,j) + sf(i-1,j) + b2 * ( sf(i,j+1) + sf(i,j-1) ) );
                sf(i,j) = beta * sfnew + ( 1.0 - beta ) * sf(i,j);
			}
		}
        double err=0.0;
		for ( int j = 0; j < ny + 1; ++ j )
		{
			for ( int i = 0; i < nx + 1; ++ i )
			{
                err = err + abs( w(i,j) - sf(i,j) );
			}
		}
        if ( iter == 0 )
		{
           err0 = err;
           if ( err0 == 0 )
		   {
              err0 = 1.0;
		   }
		}
        double rerr = err / err0;
        //stop if iteration has converged
        if ( rerr <= 0.001 || err < maxErr )
		{
            break;
		}
	}
}

void cal_stream_functionBBB(int nx, int ny, double dx, double dy, MatrixXd &sf, MatrixXd &vt)
{
    //int    maxIt  = 200;
    int    maxIt  = 1;
    double beta   = 1.5;
    //double maxErr = 1.0e-6;
    double maxErr = 1.0e-8;

    double b = dx / dy;
    double b2 = b * b;
    double dx2 = dx * dx;
    double cc = 1.0 /( 2.0 * ( 1.0 + b2 ) );
    double err0 = 1.0;
    MatrixXd w = MatrixXd::Zero( nx + 1, ny + 1 );
    for ( int iter = 0; iter < maxIt; ++ iter )
    {
        //MatrixXd w = sf; //by SOR iteration
        for ( int j = 0; j < ny + 1; ++ j )
        {
            for ( int i = 0; i < nx + 1; ++ i )
            {
                w(i,j) = sf(i,j);
            }
        }

        for ( int j = 1; j < ny; ++ j )
        {
            for ( int i = 1; i < nx; ++ i )
            {
                double sfnew = cc * ( dx2 * vt(i,j) + sf(i+1,j) + sf(i-1,j) + b2 * ( sf(i,j+1) + sf(i,j-1) ) );
                sf(i,j) = beta * sfnew + ( 1.0 - beta ) * sf(i,j);
            }
        }
        double err=0.0;
        for ( int j = 0; j < ny + 1; ++ j )
        {
            for ( int i = 0; i < nx + 1; ++ i )
            {
                err = err + abs( w(i,j) - sf(i,j) );
            }
        }
        if ( iter == 0 )
        {
            err0 = err;
            if ( err0 == 0 )
            {
                err0 = 1.0;
            }
        }
        double rerr = err / err0;
        //stop if iteration has converged
        //if ( rerr <= 0.001 || err < maxErr )
        if ( err < maxErr )
        {
            break;
        }
    }
}

void cal_stream_function( int nx, int ny, double dx, double dy, MatrixXd & sf, MatrixXd & vt )
{
    int    maxIt  = 3200;
    double beta   = 1.5;
    double maxErr = 1.0e-6;

    double b = dx / dy;
    double b2 = b * b;
    double dx2 = dx * dx;
    double cc = 1.0 /( 2.0 * ( 1.0 + b2 ) );
    double err0 = 1.0;
    MatrixXd w = MatrixXd::Zero( nx + 1, ny + 1 );
    for ( int iter = 0; iter < maxIt; ++ iter )
    {
        for ( int j = 0; j < ny + 1; ++ j )
        {
            for ( int i = 0; i < nx + 1; ++ i )
            {
                w(i,j) = sf(i,j);
            }
        }

        for ( int j = 1; j < ny; ++ j )
        {
            for ( int i = 1; i < nx; ++ i )
            {
                double sfnew = cc * ( dx2 * vt(i,j) + sf(i+1,j) + sf(i-1,j) + b2 * ( sf(i,j+1) + sf(i,j-1) ) );
                sf(i,j) = beta * sfnew + ( 1.0 - beta ) * sf(i,j);
            }
        }
        double err=0.0;
        for ( int j = 0; j < ny + 1; ++ j )
        {
            for ( int i = 0; i < nx + 1; ++ i )
            {
                err = err + abs( w(i,j) - sf(i,j) );
            }
        }
        if ( iter == 0 )
        {
            err0 = err;
            if ( err0 == 0 )
            {
                err0 = 1.0;
            }
        }
        double rerr = err / err0;
        //stop if iteration has converged
        //if ( rerr <= 0.001 || err < maxErr )
        if ( err < maxErr )
        {
            break;
        }
    }
}


void dumpfield( VectorXd&x, VectorXd&y,  MatrixXd &field, int nx, int ny, double re, const std::string &filename )
{
    std::fstream file;
    file.open( filename.c_str(), std::ios_base::out);

    std::format_to( std::ostream_iterator<char>( file ), "{0} {1} {2}\n", nx, ny, re );
    for ( int j = 0; j < ny + 1; ++ j )
    {
        for ( int i = 0; i < nx + 1; ++ i )
        {
            std::format_to(std::ostream_iterator<char>(file), "{0} {1} {2:e}\n",
                x(i), y(j), field(i,j) );
        }
    }

    file.close();
}

//-----------------------------------------------------------------------------
// Calculate right hand term of the inviscid Burgers equation
// r = -J(w,ψ) + ν ∇^2(w)
//-----------------------------------------------------------------------------
void rhs(int nx, int ny, double dx, double dy, double re, MatrixXd &w, MatrixXd &s, MatrixXd &r )
{
    //Arakawa numerical scheme for Jacobian
    double aa = 1.0/(re*dx*dx);
    double bb = 1.0/(re*dy*dy);
    double gg = 1.0/(4.0*dx*dy);
    double hh = 1.0/3.0;
    
	for ( int j = 1; j < ny; ++ j )
	{
		for ( int i = 1; i < nx; ++ i )
		{
            double j1 = gg *( ( w(i+1,j) - w(i-1,j) ) * ( s(i,j+1) - s(i,j-1) ) -
                              ( w(i,j+1) - w(i,j-1) ) * ( s(i+1,j) - s(i-1,j) ) );
            
            double j2 = gg * ( w(i+1,j) * ( s(i+1,j+1) - s(i+1,j-1) ) -
                               w(i-1,j) * ( s(i-1,j+1) - s(i-1,j-1) ) -
                               w(i,j+1) * ( s(i+1,j+1) - s(i-1,j+1) ) +
                               w(i,j-1) * ( s(i+1,j-1) - s(i-1,j-1) ) );
            
            double j3 = gg * ( w(i+1,j+1) * ( s(i,j+1) - s(i+1,j) ) -
                               w(i-1,j-1) * ( s(i-1,j) - s(i,j-1) ) -
                               w(i-1,j+1) * ( s(i,j+1) - s(i-1,j) ) +
                               w(i+1,j-1) * ( s(i+1,j) - s(i,j-1) ) );
            
            double jac = ( j1 + j2 + j3 ) * hh;
            
            //Central difference for Laplacian
            r(i,j) = -jac + aa * ( w(i+1,j) - 2.0 * w(i,j) + w(i-1,j) ) +
                            bb * ( w(i,j+1) - 2.0 * w(i,j) + w(i,j-1) );
		}
	}
}

void rhsnew(int nx, int ny, double dx, double dy, double re, MatrixXd &w, MatrixXd &r )
{
    MatrixXd s = MatrixXd::Zero( nx + 2, ny + 2 );
    cal_stream_function( nx, ny, dx, dy, s, w );

    //Arakawa numerical scheme for Jacobian
    double aa = 1.0/(re*dx*dx);
    double bb = 1.0/(re*dy*dy);
    double gg = 1.0/(4.0*dx*dy);
    double hh = 1.0/3.0;

    for ( int j = 1; j < ny; ++ j )
    {
        for ( int i = 1; i < nx; ++ i )
        {
            double j1 = gg *( ( w(i+1,j) - w(i-1,j) ) * ( s(i,j+1) - s(i,j-1) ) -
                ( w(i,j+1) - w(i,j-1) ) * ( s(i+1,j) - s(i-1,j) ) );

            double j2 = gg * ( w(i+1,j) * ( s(i+1,j+1) - s(i+1,j-1) ) -
                w(i-1,j) * ( s(i-1,j+1) - s(i-1,j-1) ) -
                w(i,j+1) * ( s(i+1,j+1) - s(i-1,j+1) ) +
                w(i,j-1) * ( s(i+1,j-1) - s(i-1,j-1) ) );

            double j3 = gg * ( w(i+1,j+1) * ( s(i,j+1) - s(i+1,j) ) -
                w(i-1,j-1) * ( s(i-1,j) - s(i,j-1) ) -
                w(i-1,j+1) * ( s(i,j+1) - s(i-1,j) ) +
                w(i+1,j-1) * ( s(i+1,j) - s(i,j-1) ) );

            double jac = ( j1 + j2 + j3 ) * hh;

            //Central difference for Laplacian
            r(i,j) = -jac + aa * ( w(i+1,j) - 2.0 * w(i,j) + w(i-1,j) ) +
                bb * ( w(i,j+1) - 2.0 * w(i,j) + w(i,j-1) );
        }
    }
}

//-----------------------------------------------------------------------------
// Compute numerical solution
//   - Time integration using Runge-Kutta third order
//   - 2nd-order finite difference discretization
//-----------------------------------------------------------------------------
void numerical( int nx, int ny, int nt, int ns, VectorXd & x, VectorXd & y, double dx, double dy, double dt, double re, MatrixXd & wn )
{
    MatrixXd wt = MatrixXd::Zero( nx + 2, ny + 2 ); // temporary array during RK3 integration
    MatrixXd r  = MatrixXd::Zero( nx + 1, ny + 1 ); // right hand side
    MatrixXd sn = MatrixXd::Zero( nx + 1, ny + 1 );
    MatrixXd ut = MatrixXd::Zero( nx + 1, ny + 1 );

    int m = 1; // record index
    int freq = int( nt / ns );
    std::cout << "freq = " << freq << "\n";

    for ( int k = 0; k < nt; ++ k )
    {
        //Compute right-hand-side from vorticity
        rhs( nx, ny, dx, dy, re, wn, sn, r );

        for ( int j = 1; j < ny; ++ j )
        {
            for ( int i = 1; i < nx; ++ i )
            {
                wt( i, j ) = wn( i, j ) + dt * r( i, j );
            }
        }

        boundary( nx, ny, wt );

        //compute streamfunction from vorticity
        cal_stream_function( nx, ny, dx, dy, sn, wt );

        //Compute right-hand-side from vorticity
        rhs( nx, ny, dx, dy, re, wt, sn, r );

        for ( int j = 1; j < ny; ++ j )
        {
            for ( int i = 1; i < nx; ++ i )
            {
                wt( i, j ) = 0.75 * wn( i, j ) + 0.25 * wt( i, j ) + 0.25 * dt * r( i, j );
            }
        }


        boundary( nx, ny, wt );
        //compute streamfunction from vorticity
        cal_stream_function( nx, ny, dx, dy, sn, wt );
        //Compute right-hand-side from vorticity
        rhs( nx, ny, dx, dy, re, wt, sn, r );

        for ( int j = 1; j < ny; ++ j )
        {
            for ( int i = 1; i < nx; ++ i )
            {
                wn( i, j ) = ( 1.0 / 3.0 ) * wn( i, j ) + ( 2.0 / 3.0 ) * wt( i, j ) + ( 2.0 / 3.0 ) * dt * r( i, j );
            }
        }

        boundary( nx, ny, wn );

        //compute streamfunction from vorticity
        cal_stream_function( nx, ny, dx, dy, sn, wn );

        if ( (k+1) % freq == 0 )
        {
            std::cout << k+1 << "\n";
            std::string filename = "vm" + std::to_string( m ) + ".txt";
            std::cout << "filename = " << filename << "\n";

            for ( int j = 0; j < ny + 1; ++ j )
            {
                for ( int i = 0; i < nx + 1; ++ i )
                {
                    ut( i, j ) = wn( i + 1, j + 1 );
                }
            }
            dumpfield( x, y, ut, nx, ny, re, filename );
            m = m + 1;
        }
    }
}

void numericalnew( int nx, int ny, int nt, int ns, VectorXd & x, VectorXd & y, double dx, double dy, double dt, double re, MatrixXd & wn )
{
    MatrixXd wt = MatrixXd::Zero( nx + 2, ny + 2 ); // temporary array during RK3 integration
    MatrixXd r  = MatrixXd::Zero( nx + 1, ny + 1 ); // right hand side
    MatrixXd sn = MatrixXd::Zero( nx + 1, ny + 1 );
    MatrixXd ut = MatrixXd::Zero( nx + 1, ny + 1 );

    int m = 1; // record index
    int freq = int( nt / ns );
    std::cout << "freq = " << freq << "\n";

    for ( int k = 0; k < nt; ++ k )
    {
        //Compute right-hand-side from vorticity
        rhsnew( nx, ny, dx, dy, re, wn, r );

        for ( int j = 1; j < ny; ++ j )
        {
            for ( int i = 1; i < nx; ++ i )
            {
                wt( i, j ) = wn( i, j ) + dt * r( i, j );
            }
        }

        boundary( nx, ny, wt );

        //compute streamfunction from vorticity
        //cal_stream_function( nx, ny, dx, dy, sn, wt );

        //Compute right-hand-side from vorticity
        rhsnew( nx, ny, dx, dy, re, wt, r );

        for ( int j = 1; j < ny; ++ j )
        {
            for ( int i = 1; i < nx; ++ i )
            {
                wt( i, j ) = 0.75 * wn( i, j ) + 0.25 * wt( i, j ) + 0.25 * dt * r( i, j );
            }
        }


        boundary( nx, ny, wt );
        //compute streamfunction from vorticity
        //cal_stream_function( nx, ny, dx, dy, sn, wt );
        //Compute right-hand-side from vorticity
        //rhs( nx, ny, dx, dy, re, wt, sn, r );
        rhsnew( nx, ny, dx, dy, re, wt, r );

        for ( int j = 1; j < ny; ++ j )
        {
            for ( int i = 1; i < nx; ++ i )
            {
                wn( i, j ) = ( 1.0 / 3.0 ) * wn( i, j ) + ( 2.0 / 3.0 ) * wt( i, j ) + ( 2.0 / 3.0 ) * dt * r( i, j );
            }
        }

        boundary( nx, ny, wn );

        //compute streamfunction from vorticity
        //cal_stream_function( nx, ny, dx, dy, sn, wn );

        if ( (k+1) % freq == 0 )
        {
            std::cout << k+1 << "\n";
            std::string filename = "vm" + std::to_string( m ) + ".txt";
            std::cout << "filename = " << filename << "\n";

            for ( int j = 0; j < ny + 1; ++ j )
            {
                for ( int i = 0; i < nx + 1; ++ i )
                {
                    ut( i, j ) = wn( i + 1, j + 1 );
                }
            }
            dumpfield( x, y, ut, nx, ny, re, filename );
            m = m + 1;
        }
    }
}

//initial condition for vortex merger problem
void vm_ic( int nx, int ny, VectorXd & x, VectorXd & y, MatrixXd &w )
{
    using namespace std::numbers;
    double sigma = pi;
    double xc1 = pi - pi / 4.0;
    double yc1 = pi;
    double xc2 = pi + pi / 4.0;
    double yc2 = pi;
    //std::cout << "xc1 = " << xc1 << " yc1 = " << yc1 << "\n";
    //std::cout << "xc2 = " << xc2 << " yc2 = " << yc2 << "\n";

    for ( int j = 1; j < ny + 2; ++ j )
    {
        for ( int i = 1; i < nx + 2; ++ i )
        {
            double dx1 = x( i - 1 ) - xc1;
            double dy1 = y( j - 1 ) - yc1;
            double c1 = exp( - sigma * ( dx1 * dx1 + dy1 * dy1) );
            double dx2 = x( i - 1 ) - xc2;
            double dy2 = y( j - 1 ) - yc2;
            double c2 = exp( - sigma * ( dx2 * dx2 + dy2 * dy2) );

            w( i, j ) = c1 + c2;
            //std::cout << "i,j=" << i << " " << j <<" c1 = " << c1 << " c2 = " << c2 << "\n";
        }
    }
}

void boundary( int nx, int ny, MatrixXd &f )
{
    //ghost points
    for ( int j = 0; j < ny + 2; ++ j )
    {
        f( 0   , j ) = f( nx, j );
        f( nx+1, j ) = f( 1 , j );
    }

    for ( int i = 0; i < nx + 2; ++ i )
    {
        f( i, 0    ) = f( i, ny );
        f( i, ny+1 ) = f( i, 1  );
    }

}

int main( int argc, char **argv )
{
    int nx = 128;
    int ny = 128;

    double dt = 0.01;
    double tf = 20.0;
    int nt = int(tf/dt);
    double re = 1000.0;
    int ns = 10;

    MatrixXd wn  = MatrixXd::Zero( nx + 2, ny + 2 );
    MatrixXd un  = MatrixXd::Zero( nx + 1, ny + 1 );
    MatrixXd un0 = MatrixXd::Zero( nx + 1, ny + 1 );
    MatrixXd ue  = MatrixXd::Zero( nx + 1, ny + 1 );
    MatrixXd uerror  = MatrixXd::Zero( nx + 1, ny + 1 );

    double xmin = 0.0;
    double xmax = 2.0*std::numbers::pi;
    double ymin = 0.0;
    double ymax = 2.0*std::numbers::pi;
    double dx = ( xmax - xmin ) / nx;
    double dy = ( ymax - ymin ) / ny;

    VectorXd x( nx + 1 );
    VectorXd y( nx + 1 );

    for ( int i = 0; i < nx + 1; ++ i )
    {
        x(i) = xmin + dx * ( i );
    }

    for ( int j = 0; j < ny + 1; ++ j )
    {
        y(j) = ymin + dy * ( j );
    }


    using clock_type = std::chrono::time_point<std::chrono::system_clock>;
    clock_type time_now, time_old;

    time_old = std::chrono::system_clock::now();

    vm_ic( nx, ny, x, y, wn );
    boundary( nx, ny, wn );

    for ( int j = 0; j < ny + 1; ++ j )
    {
        for ( int i = 0; i < nx + 1; ++ i )
        {
            un0(i,j) = wn(i+1,j+1);
        }
    }

    dumpfield( x, y, un0, nx, ny, re, "vm0.txt" );

    //numerical( nx, ny, nt, ns, x, y, dx, dy, dt, re, wn );
    numericalnew( nx, ny, nt, ns, x, y, dx, dy, dt, re, wn );

    for ( int j = 0; j < ny + 1; ++ j )
    {
        for ( int i = 0; i < nx + 1; ++ i )
        {
            un(i,j) = wn(i+1,j+1);
        }
    }

    dumpfield( x, y, un, nx, ny, re, "field_final.txt" );

    time_now = std::chrono::system_clock::now();
    double elapsed_ms = std::chrono::duration_cast<std::chrono::milliseconds>
                        ( time_now - time_old ).count();

    double t = elapsed_ms / 1000.0;

    std::cout << "CPU Time = " << t << " s\n";

    return 0;
}
