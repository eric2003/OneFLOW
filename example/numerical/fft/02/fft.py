import numpy as np
import matplotlib.pyplot as plt
    
def FFT(P):
    #P - [p0,p1,...,pn-1] coeff representation
    n = len(P) # n is a power of 2
    #print("len(P)=",n)
    if n == 1:
        return P

    nh = int( n / 2 )
    #print("nh=",nh)
    Pe = np.zeros( nh, dtype=complex )
    Po = np.zeros( nh, dtype=complex )
    for j in range( nh ):
        #print("j=",j)
        Pe[ j ] = P[ 2 * j ]
        Po[ j ] = P[ 2 * j + 1 ]
    #print("Pe=",Pe)
    #print("Po=",Po)
    ye = FFT(Pe)
    yo = FFT(Po)
    y = np.zeros( n, dtype=complex )
    
    wn = np.exp(-1j * 2 * np.pi / n)    
    w  = complex(1, 0) 
    for j in range( nh ):
        #print( "j=", j, "nh=", nh )
        y[j   ] = ye[j] + w * yo[j]
        y[j+nh] = ye[j] - w * yo[j]
        w *= wn
    return y
    
def IFFT(P):
    #P - [p0,p1,...,pn-1] coeff representation
    n = len(P) # n is a power of 2
    #print("len(P)=",n)
    if n == 1:
        return P

    nh = int( n / 2 )

    Pe = np.zeros( nh, dtype=complex )
    Po = np.zeros( nh, dtype=complex )
    for j in range( nh ):
        Pe[ j ] = P[ 2 * j ]
        Po[ j ] = P[ 2 * j + 1 ]
    ye = FFT(Pe)
    yo = FFT(Po)
    y = np.zeros( n, dtype=complex )
    
    wn = np.exp(1j * 2 * np.pi / n)
    w  = complex(1, 0) 
    for j in range( nh ):
        y[j   ] = ye[j] + w * yo[j]
        y[j+nh] = ye[j] - w * yo[j]
        w *= wn
    coef = 1.0 / n    
    for i in range( n ):
        y[i] *= coef
    return y    
  
p=[0,1,2,3]
print("p=",p)
q=FFT(p)
print("q=",q)

pp = IFFT(q)
print("pp=",pp)

# using np.fft() method
a = np.fft.fft(p)
print("a=",a)
b = np.fft.ifft(a)
print("b=",b)
    
