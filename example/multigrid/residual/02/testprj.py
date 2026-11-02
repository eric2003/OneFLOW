import numpy as np

A1=np.array([[1, -1], [21, -20]])
print("A1=\n",A1)

b1 = np.array([-1,-19])
print("b1=\n",b1)

x1 = np.linalg.solve(A1, b1)
d1 = A1@x1
print("d1=\n",d1)


print("x1=\n",x1)

A2 = np.array([[1, -1], [3, -1]])
print("A2=\n",A2)

b2 = np.array([-1,1])

print("b2=\n",b2)

x2 = np.linalg.solve(A2, b2)
print("x2=\n",x2)

d2 = A2@x2
print("d2=\n",d2)

v1=np.array([1.95,3])
v2=v1

r1 = b1-A1@v1
r2 = b2-A2@v2

r1m = np.linalg.norm(r1)
r2m = np.linalg.norm(r2)

print("norm(r1)=\n",r1m)
print("norm(r2)=\n",r2m)
