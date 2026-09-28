import numpy as np

n = 5
A = np.zeros( ( n, n ) )
U = np.zeros( ( n, n ) )

for j in range(0, n):
    for i in range(0, n):
      if ( i== j ):
          A[i,j] = 2
          if ( i != 0 ):
              A[i,j-1] = -1
          if ( i != n-1 ):
              U[i,j+1] = 1

print("A=\n",A)
print("U=\n",U)

AInv = np.linalg.inv(A)

print("AInv=\n",AInv)

RG= AInv@U

B=A@AInv

print("B=\n",B)
print("RG=\n",RG)
