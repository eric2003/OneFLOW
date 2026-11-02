import numpy as np
import matplotlib.pyplot as plt
import csv
import sys

class MyPlot:
    def __init__( self ):
        self.x = []
        self.r = []
        self.u = []
        self.m = []
        self.p = []
        self.ReadData()
    
    def AddData( self, xm, rm, um, mm, pm ):
        self.x.append( xm )
        self.r.append( rm )
        self.u.append( um )
        self.m.append( mm )
        self.p.append( pm )
        
    def ReadData( self ):
        with open('../sod_theory.plt', 'r') as f:
            for index, line in enumerate(f):
                words = line.strip().split()
                self.x.append( float(words[0]) )
                self.r.append( float(words[1]) )
                self.u.append( float(words[2]) )
                self.m.append( float(words[3]) )
                self.p.append( float(words[4]) )
        self.ComputeEnergy()
                
    def ComputeEnergy( self ):
        num = len(self.x)
        self.e = np.zeros( num )
        print("self.e.len=", len(self.e))
        gama = 1.4
        for i in range(0, num ):
            um = self.u[i]
            rm = self.r[i]
            pm = self.p[i]
            self.e[i] = (1.0/(gama-1.0))* pm/rm + 0.5 * ( um * um )
    
    def PlotTheory( self ):
        plt.figure("Exact solution for the Sod's shock-tube problem", figsize=(10, 8), dpi=100)
        plt.subplot(2, 2, 1)
        plt.plot(self.x, self.r, linewidth=1.0, label="density")
        plt.xlabel("$x$")
        plt.ylabel(r"$\rho$")
        plt.legend()
        
        plt.subplot(2, 2, 2)
        plt.plot(self.x, self.u, linewidth=1.0, label="velocity")
        plt.xlabel("$x$")
        plt.ylabel(r"$u$")
        plt.legend()
        
        plt.subplot(2, 2, 3)
        plt.plot(self.x, self.m, linewidth=1.0, label="mach number")
        plt.xlabel("$x$")
        plt.ylabel(r"$m$")
        plt.legend()
        
        plt.subplot(2, 2, 4)
        plt.plot(self.x, self.p, linewidth=1.0, label="pressure")
        plt.xlabel("$x$")
        plt.ylabel(r"$p$")
        plt.legend()
        
        plt.tight_layout()
        
        plt.show()
        
    def PlotCompare( self, x, q, title ):
        numPoints = len( q[:,0] )
        print("numPoints=",numPoints)
        
        rr = np.zeros( numPoints )
        uu = np.zeros( numPoints )
        pp = np.zeros( numPoints )
        ee = np.zeros( numPoints )
        
        gama = 1.4
        for i in range( 0, numPoints ):
            rho  = q[ i, 0 ]
            rhou = q[ i, 1 ]
            rhoe = q[ i, 2 ]
            rr[i] = rho
            uu[i] = rhou / rho
            ee[i] = rhoe / rho
            pp[i] = ( gama - 1.0 ) * ( rhoe - 0.5 * rho * ( uu[i] * uu[i] ) )
          
        sizes = 4
        #plt.figure("Sod's shock-tube problem+Rusanov Scheme+WENO-5 reconstruction", figsize=(10, 8), dpi=100)
        plt.figure(title, figsize=(10, 8), dpi=100)
        plt.subplot(2, 2, 1)
        plt.plot(self.x, self.r, color='black', linewidth=1.0, label="theory")
        plt.scatter(x, rr,  marker= "o", s=sizes, facecolors='none', edgecolors='blue', label="OneFLOW-CFD" )
            
        plt.xlabel("$x$")
        plt.ylabel(r"$\rho$")
        plt.legend()
        
        plt.subplot(2, 2, 2)
        plt.plot(self.x, self.u, color='black', linewidth=1.0, label="theory")
        plt.scatter(x, uu,  marker= "o", s=sizes, facecolors='none', edgecolors='red', label="OneFLOW-CFD" )
        plt.xlabel("$x$")
        plt.ylabel(r"$u$")
        plt.legend()
        
        plt.subplot(2, 2, 3)
        plt.plot(self.x, self.e, color='black', linewidth=1.0, label="theory")
        plt.scatter(x, ee,  marker= "o", s=sizes, facecolors='none', edgecolors='green', label="OneFLOW-CFD" )
        plt.xlabel("$x$")
        plt.ylabel(r"$E$")
        plt.legend()
        
        plt.subplot(2, 2, 4)
        plt.plot(self.x, self.p, color='black', linewidth=1.0, label="theory")
        plt.scatter(x, pp, marker= "o", s=sizes, facecolors='none', edgecolors='orange', label="OneFLOW-CFD" )
        plt.xlabel("$x$")
        plt.ylabel(r"$p$")
        plt.legend()
        
        plt.tight_layout()
        
        plt.show()

nvar = len(sys.argv)
print('nvar=',nvar)
print('sys.argv=',sys.argv)

scheme = 'Rusanov'
nt = 2000
if nvar >= 2:
    scheme = sys.argv[1]
    print('scheme=',scheme)
    
if nvar >= 3:
    mt = sys.argv[2]
    nt = int(mt)
    print('nt=',nt)

with open('field_final0.csv', newline='') as csvfile:
    readCSV = csv.reader(csvfile, delimiter= ' ')
    icount = 0
    for row in readCSV:
        icount += 1

ni = icount
print("ni=",ni)

nm = 3

q = np.zeros( (ni, nm ) )
x = np.zeros( (ni) )


filename = 'field_final'+str(nt)+'.csv'

with open(filename, newline='') as csvfile:
    readCSV = csv.reader(csvfile, delimiter= ' ')
    i = 0
    for row in readCSV:
        x[i]    = float(row[0])
        q[i][0] = float(row[1])
        q[i][1] = float(row[2])
        q[i][2] = float(row[3])
        i += 1

title = "Sod's shock-tube problem+ " + scheme + " Scheme+WENO-5 reconstruction"
mplot = MyPlot()
mplot.PlotCompare(x, q, title)
