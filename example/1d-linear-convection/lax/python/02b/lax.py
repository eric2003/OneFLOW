import numpy as np
import matplotlib.pyplot as plt

N = 41
t_min = 0
t_max = 1
k = (t_max-t_min)/N
x_min = 0
x_max = 1
a = 1
h = (x_max-x_min)/N
alpha = a*k/(2*h)
x = np.arange(x_min, x_max+2*h, h)
print(f'{x=}')
u_0 = np.zeros_like(x)
u_0[np.where((x>=0.6) & (x <= 0.8))] = 1.0
u_lf = u_0.copy()
u_cs = u_0.copy()
u_num1 = u_0.copy()
u_num2 = u_0.copy()
t = 10000

plt.clf()

for i in range(10):
    for j in range(N):
        u_num1[j] = (u_lf[j+1]+u_lf[j-1])/2 - alpha*(u_lf[j+1]-u_lf[j-1])
        u_num2[j] = u_cs[j] - alpha*(u_cs[j+1]-u_cs[j-1])
    u_lf = u_num1.copy()
    u_cs = u_num2.copy()

u_ex = np.zeros_like(x)
u_ex[np.where((x>=0.6) & (x <= 0.8))] = 1.0

u_ex=u_ex%len(u_ex)
u_lf=u_lf%len(u_lf)
u_cs=u_cs%len(u_cs)

plt.plot(x, u_ex, 'r-', label="Exact solution", fillstyle='none')
plt.plot(x, u_lf, 'o', label="Lax-Friedrichs", fillstyle='none')
plt.plot(x, u_cs, '^', label="Central scheme", fillstyle='none')
plt.axis((0, 1, -.5, 1.5))
plt.legend(loc='lower right')
plt.suptitle("t = %1d" % (t))

plt.show()