import numpy as np
from fractions import Fraction

# 定义原始矩阵 original_matrix
original_matrix = np.array([[1, 2], [1, 4]])

inverse = np.linalg.inv(original_matrix)

print(f'{original_matrix=}')
print(f'{inverse=}')

# 计算两个矩阵的乘积
product = np.dot(original_matrix, inverse)
print(f'{product=}')

m = np.array([1, 1.0/2])

m1 = np.dot(m, inverse)

print(f'{m1=}')