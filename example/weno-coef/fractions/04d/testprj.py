import numpy as np

# 定义列表 ii
ii = [-2, -1, 0]

# 将列表 ii 转化为一个 1x3 的 NumPy 二维数组 matrix
matrix = np.array(ii).reshape(1, -1)

# 打印结果
print("Matrix:")
print(matrix)