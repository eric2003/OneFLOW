import numpy as np

A1 = [[1, -1], [21, -20]]
print("A1=\n",A1)

b1 = [-1,-19]
print("b1=\n",b1)

x1 = np.linalg.solve(A1, b1)
d1 = A1@x1
print("d1=\n",d1)


print("x1=\n",x1)

A2 = [[1, -1], [3, -1]]
print("A2=\n",A2)

b2 = [-1,1]

print("b2=\n",b2)

x2 = np.linalg.solve(A2, b2)
print("x2=\n",x2)

d2 = A2@x2
print("d2=\n",d2)


