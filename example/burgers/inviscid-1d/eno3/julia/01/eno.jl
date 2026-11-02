using OffsetArrays

# Global constants and variables
const nx = 40
const ighost = 10
const iorder = 3
const ishift = ighost + 1
const ist = 1 + ishift
const ied = nx + ishift
const ntcell = nx + ishift + ighost
const isize = iorder * (iorder + 1)
const pi = 3.14159265358979323846

# Global arrays with zero-based indexing
il = OffsetArray(zeros(Int, nx + 1), 0:nx)
ir = OffsetArray(zeros(Int, nx + 1), 0:nx)
coef = OffsetArray(zeros(iorder + 1, iorder), 0:iorder, 0:iorder-1)
dd = OffsetArray(zeros(ighost, ntcell + 1), 0:ighost-1, 0:ntcell)
up1_2m = OffsetArray(zeros(nx + 1), 0:nx)
up1_2p = OffsetArray(zeros(nx + 1), 0:nx)
flux = OffsetArray(zeros(nx + 1), 0:nx)
res = OffsetArray(zeros(nx), 0:nx-1)
dt = 0.0

# Mesh module variables
xstart = 0.0
xend = 0.0
dx = 0.0
x = OffsetArray(zeros(ntcell + 2), 0:ntcell+1)
xcc = OffsetArray(zeros(ntcell + 1), 0:ntcell)

# Field module variables
u = OffsetArray(zeros(ntcell + 1), 0:ntcell)
un = OffsetArray(zeros(ntcell + 1), 0:ntcell)

function residual(q)
    reconstruction(q)
    engquist_osher_flux(up1_2m, up1_2p, flux)
    for i in 0:nx-1
        res[i] = -(flux[i + 1] - flux[i]) / dx
    end
end

function reconstruction(q)
    # Choose the stencil by ENO method
    dd[0, 1:ntcell] .= q[1:ntcell]  # Matches Python’s dd[0, 1:ntcell + 1]
    
    for m in 1:iorder-1
        for j in 1:ntcell-1
            dd[m, j] = dd[m-1, j+1] - dd[m-1, j]
        end
    end
    
    for i in 0:nx
        il[i] = i
        ir[i] = i + 1
        for m in 1:iorder-1
            if abs(dd[m, il[i]-1+ishift]) <= abs(dd[m, il[i]+ishift])
                il[i] -= 1
            end
            if abs(dd[m, ir[i]-1+ishift]) <= abs(dd[m, ir[i]+ishift])
                ir[i] -= 1
            end
        end
    end
    
    # Reconstruction u(j+1/2)
    for i in 0:nx
        k1 = il[i]
        k2 = ir[i]
        l1 = i - k1 + 1
        l2 = i - k2 + 1
        up1_2m[i] = 0.0
        up1_2p[i] = 0.0
        for m in 0:iorder-1
            up1_2m[i] += q[k1 + ishift + m] * coef[l1, m]
            up1_2p[i] += q[k2 + ishift + m] * coef[l2, m]
        end
    end
end

function engquist_osher_flux(up1_2m, up1_2p, flux)
    for i in 0:nx
        if up1_2m[i] >= 0
            if up1_2p[i] >= 0
                flux[i] = 0.5 * up1_2m[i] * up1_2m[i]
            else
                flux[i] = 0.5 * (up1_2m[i] * up1_2m[i] + up1_2p[i] * up1_2p[i])
            end
        else
            if up1_2p[i] >= 0
                flux[i] = 0.0
            else
                flux[i] = 0.5 * up1_2p[i] * up1_2p[i]
            end
        end
    end
end

function boundary(u)
    for i in -ighost:0
        u[ishift + i] = u[ied + i]
    end
    for i in 1:ighost
        u[ied + i] = u[ishift + i]
    end
end

function update_oldfield(qn, q)
    qn .= q
end

function init_coef()
    coef[0, :] = [ 11.0/6.0, -7.0/6.0,  1.0/3.0 ]
    coef[1, :] = [  1.0/3.0,  5.0/6.0, -1.0/6.0 ]
    coef[2, :] = [ -1.0/6.0,  5.0/6.0,  1.0/3.0 ]
    coef[3, :] = [  1.0/3.0, -7.0/6.0, 11.0/6.0 ]	
end

function init_mesh()
    global xstart, xend, dx, x, xcc
    xstart = -1.0
    xend = 1.0
    dx = (xend - xstart) / nx
    xstart0 = xstart - ishift * dx
    
    for i in 1:ntcell+1
        x[i] = xstart0 + (i - 1) * dx
    end
    
    for i in 1:ntcell
        xcc[i] = 0.5 * (x[i] + x[i + 1])
    end
end

function init_field()
    for i in ist:ied
        u[i] = 0.25 + 0.5 * sin(pi * xcc[i])
    end
    boundary(u)
    update_oldfield(un, u)
end

function runge_kutta_3()
    global u, un, dt
    residual(u)
    for i in 0:nx-1
        j = i + 1 + ishift
        u[j] = u[j] + dt * res[i]
    end
    boundary(u)
    
    residual(u)
    for i in 0:nx-1
        j = i + 1 + ishift
        u[j] = 0.75 * un[j] + 0.25 * u[j] + 0.25 * dt * res[i]
    end
    boundary(u)
    
    residual(u)
    c1, c2, c3 = 1.0/3.0, 2.0/3.0, 2.0/3.0
    for i in 0:nx-1
        j = i + 1 + ishift
        u[j] = c1 * un[j] + c2 * u[j] + c3 * dt * res[i]
    end
    boundary(u)
    update_oldfield(un, u)
end

function visualize()
    open("solution_total.plt", "w") do f1
        for i in 1:ntcell
            println(f1, "$(xcc[i])\t$(u[i])")
        end
    end
    
    open("solution.plt", "w") do f2
        for i in ist:ied
            println(f2, "$(xcc[i])\t$(u[i])")
        end
    end
end

function main()
    global dt
    init_coef()
    init_mesh()
    init_field()
    
    println("Input T: ")
    simu_time = parse(Float64, readline())
    dt = dx * 0.5
    t = 0.0
    
    while t < simu_time
        runge_kutta_3()
        t += dt
        if t + dt > simu_time
            dt = simu_time - t
        end
    end
    
    println(t)
    visualize()
end

# Run the main function
main()