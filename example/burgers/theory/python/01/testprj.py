import numpy as np
import matplotlib.pyplot as plt

# 参数设置
nu = 0.1  # 粘性系数
x = np.linspace(-1, 1, 400)  # x 范围从 -1 到 1
t = 0.1  # 时间 t

# 计算解析解
u = 2 * np.pi * nu * np.tan(np.pi * x)

# 绘图
plt.figure(figsize=(10, 6))
plt.plot(x, u, label=f't = {t}', linewidth=2)
plt.title('解析解 u(x, t) = 2 * pi * nu * tan(pi * x)', fontsize=14)
plt.xlabel('x', fontsize=12)
plt.ylabel('u(x, t)', fontsize=12)
plt.legend(fontsize=12)
plt.grid(True)
plt.show()