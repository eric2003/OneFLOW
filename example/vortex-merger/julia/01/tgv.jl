using CPUTime
using Printf
using Plots
using FFTW

font = Plots.font("Times New Roman", 18)
plot(guidefont=font, xtickfont=font, ytickfont=font, legendfont=font)

#-----------------------------------------------------------------------------#
# Compute L-2 norm for a vector
#-----------------------------------------------------------------------------#
function compute_l2norm(nx, ny, r)
    rms = 0.0
    # println(residual)
    for j = 1:ny+1 for i = 1:nx+1
        rms = rms + r[i,j]^2
    end end
    # println(rms)
    rms = sqrt(rms/((nx+1)*(ny+1)))
    return rms
end

#-----------------------------------------------------------------------------#
# Fast poisson solver for periodic domain
#-----------------------------------------------------------------------------#
function fps(nx,ny,dx,dy,f)
    eps = 1.0e-6

    kx = Array{Float64}(undef,nx)
    ky = Array{Float64}(undef,ny)

    data = Array{Complex{Float64}}(undef,nx,ny)
    data1 = Array{Complex{Float64}}(undef,nx,ny)
    e = Array{Complex{Float64}}(undef,nx,ny)

    u = Array{Complex{Float64}}(undef,nx,ny)

    aa = -2.0/(dx*dx) - 2.0/(dy*dy)
    bb = 2.0/(dx*dx)
    cc = 2.0/(dy*dy)

    #wave number indexing
    hx = 2.0*pi/nx

    for i = 1:Int64(nx/2)
        kx[i] = hx*(i-1.0)
        kx[i+Int64(nx/2)] = hx*(i-Int64(nx/2)-1)
    end
    kx[1] = eps

    ky = kx

    for i = 1:nx
        for j = 1:ny
            data[i,j] = complex(f[i,j],0.0)
        end
    end

    e = fft(data)
    e[1,1] = 0.0
    for i = 1:nx
        for j = 1:ny
            data1[i,j] = e[i,j]/(aa + bb*cos(kx[i]) + cc*cos(ky[j]))
        end
    end

    u = real(ifft(data1))

    return u
end

#-----------------------------------------------------------------------------#
# Compute numerical solution
#   - Time integration using Runge-Kutta third order
#   - 2nd-order finite difference discretization
#-----------------------------------------------------------------------------#
function numerical(nx,ny,nt,dx,dy,dt,re,wn)

    wt = Array{Float64}(undef, nx+2, ny+2) # temporary array during RK3 integration
    r = Array{Float64}(undef, nx+2, ny+2)

    for k = 1:nt
        println(k)
        # Compute right-hand-side from vorticity
        rhs(nx,ny,dx,dy,re,wn,r)

        for i = 2:nx+1 for j = 2:ny+1
            wt[i,j] = wn[i,j] + dt*r[i,j]
        end end

        # periodic BC
        wt[nx+2,:] = wt[2,:]
        wt[:,ny+2] = wt[:,2]

        # ghost points
        wt[1,:] = wt[nx+1,:]
        wt[:,1] = wt[:,ny+1]

        # Compute right-hand-side from vorticity
        rhs(nx,ny,dx,dy,re,wt,r)

        for i = 2:nx+1 for j = 2:ny+1
            wt[i,j] = 0.75*wn[i,j] + 0.25*wt[i,j] + 0.25*dt*r[i,j]
        end end

        # periodic BC
        wt[nx+2,:] = wt[2,:]
        wt[:,ny+2] = wt[:,2]

        # ghost points
        wt[1,:] = wt[nx+1,:]
        wt[:,1] = wt[:,ny+1]

        # Compute right-hand-side from vorticity
        rhs(nx,ny,dx,dy,re,wt,r)

        for i = 2:nx+1 for j = 2:ny+1
            wn[i,j] = (1.0/3.0)*wn[i,j] + (2.0/3.0)*wt[i,j] + (2.0/3.0)*dt*r[i,j]
        end end

        # periodic BC
        wn[nx+2,:] = wn[2,:]
        wn[:,ny+2] = wn[:,2]

        # ghost points
        wn[1,:] = wn[nx+1,:]
        wn[:,1] = wn[:,ny+1]
    end

    return wn[2:nx+2,2:ny+2]
end

#-----------------------------------------------------------------------------#
# Calculate right hand term of the inviscid Burgers equation
# r = -J(w,ψ) + ν ∇^2(w)
#-----------------------------------------------------------------------------#
function rhs(nx,ny,dx,dy,re,w,r)

    # compute streamfunction from vorticity
    s = Array{Float64}(undef, nx+2, ny+2)
    f = Array{Float64}(undef, nx, ny)
    f = -w[2:nx+1,2:ny+1]

    s[2:nx+1,2:ny+1] = fps(nx,ny,dx,dy,f)

    # periodic BC
    s[nx+2,:] = s[2,:]
    s[:,ny+2] = s[:,2]

    # ghost points
    s[1,:] = s[nx+1,:]
    s[:,1] = s[:,ny+1]

    # Arakawa numerical scheme for Jacobian
    aa = 1.0/(re*dx*dx)
    bb = 1.0/(re*dy*dy)
    gg = 1.0/(4.0*dx*dy)
    hh = 1.0/3.0

    for i = 2:nx+1 for j = 2:ny+1
        j1 = gg*((w[i+1,j]-w[i-1,j])*(s[i,j+1]-s[i,j-1]) -
                 (w[i,j+1]-w[i,j-1])*(s[i+1,j]-s[i-1,j]))

        j2 = gg*(w[i+1,j]*(s[i+1,j+1]-s[i+1,j-1]) -
                 w[i-1,j]*(s[i-1,j+1]-s[i-1,j-1]) -
                 w[i,j+1]*(s[i+1,j+1]-s[i-1,j+1]) +
                 w[i,j-1]*(s[i+1,j-1]-s[i-1,j-1]))

        j3 = gg*(w[i+1,j+1]*(s[i,j+1]-s[i+1,j]) -
                 w[i-1,j-1]*(s[i-1,j]-s[i,j-1]) -
            	 w[i-1,j+1]*(s[i,j+1]-s[i-1,j]) +
            	 w[i+1,j-1]*(s[i+1,j]-s[i,j-1]))

        jac = (j1+j2+j3)*hh

        #Central difference for Laplacian
        r[i,j] = -jac + aa*(w[i+1,j]-2.0*w[i,j]+w[i-1,j]) +
                        bb*(w[i,j+1]-2.0*w[i,j]+w[i,j-1])
        end end
end

# compute exact solution for TGV problem
function exact_tgv(nx,ny,x,y,time,re)
    ue = Array{Float64}(undef, nx+1, ny+1)

    nq = 4.0
    for i = 1:nx+1 for j = 1:ny+1
        ue[i,j] = 2.0*nq*cos(nq*x[i])*cos(nq*y[j])*
                  exp(-2.0*nq*nq*time/re)
    end end
    return ue
end


#---------------------------------------------------------------------------#
# main program
#---------------------------------------------------------------------------#
nx = 64
ny = 64

x_l = 0.0
x_r = 2.0*pi
y_b = 0.0
y_t = 2.0*pi

dx = (x_r-x_l)/nx
dy = (y_t-y_b)/ny

dt = 0.01
tf = 1.0
nt = tf/dt
re = 10.0

x = Array{Float64}(undef, nx+1)
y = Array{Float64}(undef, ny+1)

for i = 1:nx+1
    x[i] = dx*(i-1)
end
for i = 1:ny+1
    y[i] = dy*(i-1)
end

wn = Array{Float64}(undef, nx+2, ny+2)
un = Array{Float64}(undef, nx+1, ny+1)
ue = Array{Float64}(undef, nx+1, ny+1)
uerror = Array{Float64}(undef, nx+1, ny+1)

time = 0.0

wn[2:nx+2,2:ny+2] = exact_tgv(nx,ny,x,y,time,re)
# ghost points
wn[1,:] = wn[nx+1,:]
wn[:,1] = wn[:,ny+1]

val, t, bytes, gctime, memallocs = @timed begin

un = numerical(nx,ny,nt,dx,dy,dt,re,wn)

end

print("CPU Time = ", t);

time = tf
ue = exact_tgv(nx,ny,x,y,time,re)

uerror = un-ue

rms_error = compute_l2norm(nx, ny, uerror)
max_error = maximum(abs.(uerror))

println("Error details:");
println("L-2 Norm = ", rms_error);
println("Maximum Norm = ", max_error);


p1 = contour(x, y, transpose(ue), fill=true,xlabel="\$X\$", ylabel="\$Y\$", title="Exact")
p2 = contour(x, y, transpose(un), fill=true,xlabel="\$X\$", ylabel="\$Y\$", title="Numerical")
p3 = plot(p1,p2, size = (1300, 600))
savefig(p3,"tgv.pdf")
