import numpy as np
import matplotlib.pyplot as plt

def Weighted_Jacobi(v, f, n, niter, omega):
    v2 = v.copy()
    a = 0.5 * omega
    b = 1.0 - omega
    c = a    
    d = a
    for iter in range(0, niter):
        for i in range(1, n):
            v2[i] = a * v[i-1] + b * v[i] + c * v[i+1] + d * f[i]
        v = v2.copy()
    return v
    
def CalcResidual(v, f, n):
    v2 = v.copy()
    a = -1.0
    b = 2.0
    c = -1.0
    for i in range(1, n):
        av = a * v[i-1] +  b * v[i] + c * v[i+1]
        v2[i] = f[i] - av
    return v2
    
def Restrict(rh, nc):
    r2h = np.zeros( nc + 1 )
    for i in range(1, nc):
        r2h[i] = rh[2*i]
    return r2h

def Interpolate(e2h, nc):
    eh = np.zeros( 2*nc + 1 )
    for i in range(1, nc):
        eh[2*i] = e2h[i]
        eh[2*i+1] = 0.5 * ( e2h[i] + e2h[i+1] )
    return eh
    
def Correct(vh, eh, n):
    vh_new = np.zeros( n + 1 )
    for i in range(1, n):
        vh_new[i] = vh[i] + eh[i]
    return vh_new
    
def VectorNorm(vh, n): 
    return np.max(np.abs(vh))
    
def IsCoarseGrid(level, nCoarse):
    return level == nCoarse
    
def VCycle(vh, fh, n, level, nCoarse, omega):
    nIter1 = 3
    nIter2 = 3
    vh  = Weighted_Jacobi(vh, fh, n, nIter1, omega)
    if not IsCoarseGrid( level, nCoarse ):
        vh  = Weighted_Jacobi(vh, fh, n, nIter2, omega)
        rh  = CalcResidual(vh, fh, n)
        f2h = Restrict( rh, int(n/2) )
        v2h = np.zeros( int(n/2) + 1 )
        v2h = VCycle( v2h, f2h, int(n/2), level + 1, nCoarse, omega )
        eh  = Interpolate( v2h, int(n/2) )
        vh  = Correct(vh, eh, n)
    
    vh  = Weighted_Jacobi(vh, fh, n, nIter2, omega)
    return vh
    
def FMG_VCycle_Level_1(vh, fh, n, omega):
    nIter = 3
    vIter = 1
    nCoarse = 1
    rh  = CalcResidual(vh, fh, n)
    f2h = Restrict( rh, int(n/2) )
    e2h = np.zeros( int(n/2) + 1 )
    e2h = Weighted_Jacobi(e2h, f2h, int(n/2), nIter, omega)

    eh  = Interpolate( e2h, int(n/2) )
    vh = Correct(vh, eh, n)
    
    level = 0
    for i in range(0, vIter):
        vh = VCycle( vh, fh, n, level, nCoarse, omega )
    return vh
    
def FMG_VCycle_Level_2(vh, fh, n, omega):
    nIter = 3
    vIter = 1
    nCoarse = 2
    level = 0
    
    rh  = CalcResidual(vh, fh, n)
    r2h = Restrict( rh, int(n/2) )
    r4h = Restrict( r2h, int(n/4) )
    e4h = np.zeros( int(n/4) + 1 )
    e4h = Weighted_Jacobi(e4h, r4h, int(n/4), nIter, omega)
    
    e2h = Interpolate( e4h, int(n/4) )

    for i in range(0, vIter):
        e2h = VCycle( e2h, r2h, int(n/2), level+1, nCoarse, omega )
    
    eh = Interpolate( e2h, int(n/2) )
    vh = Correct(vh, eh, n)
    
    for i in range(0, vIter):
        vh = VCycle( vh, fh, n, level, nCoarse, omega )
    return vh
    
def FMG_VCycle_Level_3(vh, fh, n, omega):
    nIter = 3
    vIter = 1
    nCoarse = 3
    level = 0
    
    rh  = CalcResidual(vh, fh, n)
    r2h = Restrict( rh, int(n/2) )
    r4h = Restrict( r2h, int(n/4) )
    r8h = Restrict( r4h, int(n/8) )
    e8h = np.zeros( int(n/8) + 1 )
    e8h = Weighted_Jacobi(e8h, r8h, int(n/8), nIter, omega)
    
    e4h = Interpolate( e8h, int(n/8) )
    for i in range(0, vIter):
        e4h = VCycle( e4h, r4h, int(n/4), level+2, nCoarse, omega )

    e2h = Interpolate( e4h, int(n/4) )

    for i in range(0, vIter):
        e2h = VCycle( e2h, r2h, int(n/2), level+1, nCoarse, omega )
    
    eh = Interpolate( e2h, int(n/2) )
    vh = Correct(vh, eh, n)
    
    for i in range(0, vIter):
        vh = VCycle( vh, fh, n, level, nCoarse, omega )
    return vh
    
def FMG_VCycle_Level_4(vh, fh, n, omega):
    nIter = 3
    vIter = 1
    nCoarse = 4
    level = 0
    
    rh  = CalcResidual(vh, fh, n)
    r2h = Restrict( rh, int(n/2) )
    r4h = Restrict( r2h, int(n/4) )
    r8h = Restrict( r4h, int(n/8) )
    r16h = Restrict( r8h, int(n/16) )

    e16h = np.zeros( int(n/16) + 1 )
    e16h = Weighted_Jacobi(e16h, r16h, int(n/16), nIter, omega)
    
    e8h = Interpolate( e16h, int(n/16) )
    for i in range(0, vIter):
        e8h = VCycle( e8h, r8h, int(n/8), level+3, nCoarse, omega )

    e4h = Interpolate( e8h, int(n/8) )
    for i in range(0, vIter):
        e4h = VCycle( e4h, r4h, int(n/4), level+2, nCoarse, omega )

    e2h = Interpolate( e4h, int(n/4) )

    for i in range(0, vIter):
        e2h = VCycle( e2h, r2h, int(n/2), level+1, nCoarse, omega )
    
    eh = Interpolate( e2h, int(n/2) )
    vh = Correct(vh, eh, n)
    
    for i in range(0, vIter):
        vh = VCycle( vh, fh, n, level, nCoarse, omega )
    return vh
    
def FMG_VCycle_Level_4_NEW(vh, fh, n, omega):
    nIter = 3
    vIter = 1
    nCoarse = 4
    level = 0
    
    rh  = CalcResidual(vh, fh, n)
    r2h = Restrict( rh, int(n/2) )
    r4h = Restrict( r2h, int(n/4) )
    r8h = Restrict( r4h, int(n/8) )
    r16h = Restrict( r8h, int(n/16) )

    e16h = np.zeros( int(n/16) + 1 )
    for i in range(0, vIter):
        e16h = VCycle( e16h, r16h, int(n/16), level+4, nCoarse, omega )  

    e8h = Interpolate( e16h, int(n/16) )
    for i in range(0, vIter):
        e8h = VCycle( e8h, r8h, int(n/8), level+3, nCoarse, omega )

    e4h = Interpolate( e8h, int(n/8) )
    for i in range(0, vIter):
        e4h = VCycle( e4h, r4h, int(n/4), level+2, nCoarse, omega )

    e2h = Interpolate( e4h, int(n/4) )

    for i in range(0, vIter):
        e2h = VCycle( e2h, r2h, int(n/2), level+1, nCoarse, omega )
    
    eh = Interpolate( e2h, int(n/2) )
    vh = Correct(vh, eh, n)
    
    for i in range(0, vIter):
        vh = VCycle( vh, fh, n, level, nCoarse, omega )
    return vh
    
def FMG_VCycle( eh, rh, n, level, nCoarse, omega, vIter ):
    for i in range(0, vIter):
        eh = VCycle( eh, rh, n, level, nCoarse, omega )
    return eh
        
def FMG_VCycle_Recursive_Level_4(vh, fh, n, omega):
    nIter = 3
    vIter = 1
    nCoarse = 4
    level = 0
    
    rh  = CalcResidual(vh, fh, n)
    r2h = Restrict( rh, int(n/2) )
    r4h = Restrict( r2h, int(n/4) )
    r8h = Restrict( r4h, int(n/8) )
    r16h = Restrict( r8h, int(n/16) )

    e16h = np.zeros( int(n/16) + 1 )
    e16h = FMG_VCycle( e16h, r16h, int(n/16), level+4, nCoarse, omega, vIter )

    e8h = Interpolate( e16h, int(n/16) )
    e8h = FMG_VCycle( e8h, r8h, int(n/8), level+3, nCoarse, omega, vIter )
    
    e4h = Interpolate( e8h, int(n/8) )
    e4h = FMG_VCycle( e4h, r4h, int(n/4), level+2, nCoarse, omega, vIter )
    
    e2h = Interpolate( e4h, int(n/4) )
    e2h = FMG_VCycle( e2h, r2h, int(n/2), level+1, nCoarse, omega, vIter )

    eh = Interpolate( e2h, int(n/2) )
    vh = Correct(vh, eh, n)
    
    for i in range(0, vIter):
        vh = VCycle( vh, fh, n, level, nCoarse, omega )
    return vh    
    
def GetEh( n, rr, nLevel, nCoarse, omega, vIter ):
    m1 = int( np.power(2,nLevel) )
    nc = int( n / m1 )
    exh = np.zeros( nc + 1 )
    for i in range(0, nLevel+1):
        kLevel = nLevel - i
        #print("kLevel=",kLevel)
        km = int( np.power(2,kLevel) )
        #print("km=",km)
        kPoint = int(n/km)
        #print("kPoint=",kPoint)
        
        if i != 0:
            exh = Interpolate( exh, int(kPoint/2) )
            
        exh = FMG_VCycle( exh, rr[kLevel], kPoint, kLevel, nCoarse, omega, vIter )

    return exh

def FMG_VCycle_Recursive_Level_4_VERSION(vh, fh, n, omega):
    nIter = 3
    vIter = 1
    nCoarse = 4
    level = 0
    nLevel = 4
    
    rh  = CalcResidual(vh, fh, n)
    r2h = Restrict( rh, int(n/2) )
    r4h = Restrict( r2h, int(n/4) )
    r8h = Restrict( r4h, int(n/8) )
    r16h = Restrict( r8h, int(n/16) )
    
    rr = []
    rr.append( rh )
    rr.append( r2h )
    rr.append( r4h )
    rr.append( r8h )
    rr.append( r16h )
    
    eh = GetEh( n, rr, nLevel, nCoarse, omega, vIter )
    vh = Correct(vh, eh, n)
    return vh

def basic_visual_setup(ax):    
    ax.spines['right'].set_color('none')
    ax.spines['top'].set_color('none')
    ax.set_xlabel("X")
    ax.set_ylabel("Y",rotation='horizontal')
    ax.xaxis.set_label_coords(1.02, 0.02)
    ax.yaxis.set_label_coords(0.0, 1.0)
    ax.set_ylim(-1.0, 1.0)
    return


n = 64
m = int(n/2)
dx = 1.0 / n
h = dx

fh = np.zeros( n + 1 )

x = np.zeros( n + 1 )
y = np.zeros( n + 1 )

for i in range(0, n+1):
    x[i] = i * h
    y16m = np.sin( i * 16 * np.pi / n )
    y40m = np.sin( i * 40 * np.pi / n )
    y[i] = 0.5 * ( y16m + y40m )
    
omega = 2.0/3.0

vh = y.copy()

nCycle = 1000

res = np.zeros( nCycle + 1 )
res[0] = VectorNorm( vh, n )

for iCycle in range(0, nCycle):
    #vh = FMG_VCycle_Level_4_NEW(vh, fh, n, omega)
    #vh = FMG_VCycle_Recursive_Level_4(vh, fh, n, omega)
    vh = FMG_VCycle_Recursive_Level_4_VERSION(vh, fh, n, omega)
    
    res[iCycle+1] = VectorNorm( vh, n )
    
    
print("res=\n",res)

fig = plt.figure('OneFLOW-CFD Solver Full Multigrid V-Cycle Scheme', figsize=(8, 6), dpi=100)

ax = fig.add_subplot(1,3,1)
basic_visual_setup( ax )
mycolor  = 'blue'
mymarker = 'none'
mylabel  = 'initial field'
ax.plot(x,y,color='black',linewidth=1,marker=mymarker,markerfacecolor='none',label=mylabel)

ax = fig.add_subplot(1,3,2)
basic_visual_setup( ax )

mycolor  = 'blue'
mymarker = 'o'
mylabel  = 'coarse grid'
ax.plot(x,y,color='black',linewidth=1,linestyle='dotted',marker='none',markerfacecolor='none',label=mylabel)
ax.plot(x,vh,color=mycolor,linewidth=2,marker='none',markerfacecolor='none',label=mylabel)

ax = fig.add_subplot(1,3,3)
ii = np.linspace(0, nCycle, nCycle + 1)
ax.set_yscale("log", base=10)
ax.plot(ii,res,color='black',linewidth=2,linestyle='solid',marker='none',markerfacecolor='none',label='res')

fig.tight_layout()
plt.show()
