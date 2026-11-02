import numpy as np

n = 8
A = np.zeros( ( n, n ) )

for j in range(0, n):
    for i in range(0, n):
      if ( i== j ):
          A[i,j] = 2
          if ( i != 0 ):
              A[i,j-1] = -1

print("A=\n",A)

AInv = np.linalg.inv(A)

print("AInv=\n",AInv)

B=A@AInv

print("B=\n",B)
