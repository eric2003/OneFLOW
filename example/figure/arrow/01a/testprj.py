import matplotlib.pyplot as plt
import numpy as np

# 定义函数值
x = np.array([0, 1, 2, 3])
y = np.array([1, 4, 9, 16])

# 绘制原始点
plt.scatter(x, y, color='blue', label='Data Points')

# 绘制0阶差商
for i in range(len(x)):
    plt.plot([x[i], x[i]], [y[i], y[i]], 'k--')  # 绘制水平线

# 绘制1阶差商
for i in range(len(x) - 1):
    plt.plot([x[i], x[i+1]], [y[i], y[i+1]], 'r--')  # 绘制斜线

# 绘制2阶差商
for i in range(len(x) - 2):
    plt.plot([x[i], x[i+2]], [y[i], y[i+2]], 'g--')  # 绘制更长的斜线

# 添加图例
plt.legend(['Data Points', '0th Order Differences', '1st Order Differences', '2nd Order Differences'])

# 添加标题和标签
plt.title('Difference Quotient Visualization')
plt.xlabel('x')
plt.ylabel('f(x)')

# 显示图形
plt.show()