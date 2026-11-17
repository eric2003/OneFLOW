import numpy as np
from fractions import Fraction

# 定义原始矩阵
original_matrix = np.array([[1, -1, 13/12],
                           [1, 0, 1/12],
                           [1, 1, 13/12]])

# 定义计算出的逆矩阵
inverse_matrix = np.array([[Fraction(-1, 24), Fraction(13, 12), Fraction(-1, 24)],
                          [Fraction(-1, 2), Fraction(0, 1), Fraction(1, 2)],
                          [Fraction(1, 2), Fraction(-1, 1), Fraction(1, 2)]])

# 将 Fraction 转换为浮点数进行矩阵乘法
original_matrix = original_matrix.astype(float)
inverse_matrix = inverse_matrix.astype(float)

# 计算两个矩阵的乘积
product = np.dot(original_matrix, inverse_matrix)

# 检查乘积是否接近单位矩阵
identity_matrix = np.eye(3)

# 打印结果
print("Product of the matrices:")
print(product)
print("\nIdentity matrix:")
print(identity_matrix)
print("\nDifference:")
print(np.abs(product - identity_matrix) < 1e-9)