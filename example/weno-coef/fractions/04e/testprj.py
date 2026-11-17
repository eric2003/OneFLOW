import numpy as np

# 定义三个二维数组
a1 = np.array([[1, 1, 1]])
a2 = np.array([[2, 2, 2]])
a3 = np.array([[3, 3, 3]])

# 使用 vstack 函数将这三个数组垂直堆叠成一个矩阵
matrix = np.vstack((a1, a2, a3))

# 打印结果
print("Assembled Matrix:")
print(matrix)