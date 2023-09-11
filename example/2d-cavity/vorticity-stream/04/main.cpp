#include <iostream>
#include <fstream>
#include <format>
#include <chrono>
#include <Eigen/Dense>

using Eigen::MatrixXd;
using Eigen::VectorXd;

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

void cal_stream_function(int nx, int ny, double dx, double dy, MatrixXd &sf, MatrixXd &vt, int MaxIt, double MaxErr, double beta)
{
    double b = dx / dy;
    double b2 = b * b;
    double dx2 = dx * dx;
    double cc = 1.0 /( 2.0 * ( 1.0 + b2 ) );
    double err0 = 1.0;
	MatrixXd w = MatrixXd::Zero( nx + 1, ny + 1 );
	for ( int iter = 0; iter < MaxIt; ++ iter )
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
        if ( rerr <= 0.001 || err < MaxErr )
		{
            break;
		}
	}
}

void dumpfield( VectorXd&x, VectorXd&y,  MatrixXd &f1,  MatrixXd &f2, int nx, int ny, std::fstream &outfile )
{
    std::format_to( std::ostream_iterator<char>( outfile ), "{0} {1}\n", nx, ny );
    for ( int j = 0; j < ny + 1; ++ j )
    {
        for ( int i = 0; i < nx + 1; ++ i )
        {
            std::format_to(std::ostream_iterator<char>(outfile), "{0} {1} {2} {3}\n",
                x(i), y(j), f1(i,j), f2(i,j) );
        }
    }
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

//-----------------------------------------------------------------------------
// Compute numerical solution
//   - Time integration using Runge-Kutta third order
//   - 2nd-order finite difference discretization
//-----------------------------------------------------------------------------
void numerical(int nx, int ny, int nt, double dx, double dy, double dt, double re,
    MatrixXd &wn, MatrixXd &sn, std::vector<double> &rms, int MaxIt, double MaxErr,double beta)
{
    MatrixXd wt = MatrixXd::Zero( nx + 1, ny + 1 ); // temporary array during RK3 integration
    MatrixXd r  = MatrixXd::Zero( nx + 1, ny + 1 ); // right hand side
    MatrixXd sp = MatrixXd::Zero( nx + 1, ny + 1 ); // old streamfunction

    for ( int k = 0; k < nt; ++ k )
    {
        for ( int j = 0; j < ny + 1; ++ j )
        {
            for ( int i = 0; i < nx + 1; ++ i )
            {
                sp(i,j) = sn(i,j);
            }
        }

        //Compute right-hand-side from vorticity
        rhs(nx,ny,dx,dy,re,wn,sn,r);

        for ( int j = 1; j < ny; ++ j )
        {
            for ( int i = 1; i < nx; ++ i )
            {
                wt(i,j) = wn(i,j) + dt * r(i,j);
            }
        }				

        bc2(nx,ny,dx,dy,wt,sn);

        //compute streamfunction from vorticity
        cal_stream_function( nx, ny, dx, dy, sn, wt, MaxIt, MaxErr, beta );

        //Compute right-hand-side from vorticity
        rhs(nx,ny,dx,dy,re,wt,sn,r);

        for ( int j = 1; j < ny; ++ j )
        {
            for ( int i = 1; i < nx; ++ i )
            {
                wt(i,j) = 0.75*wn(i,j) + 0.25*wt(i,j) + 0.25*dt*r(i,j);
            }
        }


        bc2(nx,ny,dx,dy,wt,sn);
        //compute streamfunction from vorticity
        cal_stream_function( nx, ny, dx, dy, sn, wt, MaxIt, MaxErr, beta );
        //Compute right-hand-side from vorticity
        rhs(nx,ny,dx,dy,re,wt,sn,r);

        for ( int j = 1; j < ny; ++ j )
        {
            for ( int i = 1; i < nx; ++ i )
            {
                wn(i,j) = (1.0/3.0)*wn(i,j) + (2.0/3.0)*wt(i,j) + (2.0/3.0)*dt*r(i,j);
            }
        }


        bc2( nx, ny, dx, dy, wn, sn );

        //compute streamfunction from vorticity
        cal_stream_function(nx,ny,dx,dy,sn,wt,MaxIt,MaxErr,beta);

        rms[k] = 0.0;

        for ( int j = 0; j < ny + 1; ++ j )
        {
            for ( int i = 0; i < nx + 1; ++ i )
            {
                double ds = sn(i,j) - sp(i,j);

                rms[k] = rms[k] + ds * ds;
            }
        }

        rms[ k ] = sqrt( rms[ k ] / ( ( nx + 1 ) * ( ny + 1 ) ) );
        std::cout << std::setprecision( 20 );
        std::cout << k << " " << rms[k] << "\n";
    }
}
 
int main( int argc, char **argv )
{
    int nx = 64;
    int ny = 64;

    double dt = 0.001;
    double tf = 10.0;
    int nt = int(tf/dt);
    double re = 100.0;

    //parameters for SOR iteration
    int MaxIt=100;
    double Beta=1.5;
    double MaxErr=1.0e-6;

    MatrixXd wn = MatrixXd::Zero( nx + 1, ny + 1 );
    MatrixXd sn = MatrixXd::Zero( nx + 1, ny + 1 );
	
	std::vector<double> rms(nt,0);

    double xmin = 0.0;
    double xmax = 1.0;
    double ymin = 0.0;
    double ymax = 1.0;
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

    numerical(nx,ny,nt,dx,dy,dt,re,wn,sn,rms,MaxIt,MaxErr,Beta);

    time_now = std::chrono::system_clock::now();
    double elapsed_ms = std::chrono::duration_cast<std::chrono::milliseconds>
                        ( time_now - time_old ).count();

    double t = elapsed_ms / 1000.0;

    std::fstream residual_plot;
    residual_plot.open( "residual_plot.txt", std::ios_base::out );
	for ( int n = 0; n < nt; ++ n )
	{
        std::format_to(std::ostream_iterator<char>(residual_plot),"{0} {1}\n", n+1,rms[n] );
	}
	residual_plot.close();

    std::cout << "CPU Time = " << t << " s\n";

    //create text file for final field
    std::fstream field_final;
    field_final.open( "field_final.txt", std::ios_base::out );
    dumpfield( x, y, wn, sn, nx, ny, field_final );
    field_final.close();

    return 0;
}
