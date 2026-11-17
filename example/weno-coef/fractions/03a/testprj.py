import numpy as np
from fractions import Fraction

# 定义矩阵 A 的元素
A = np.array([
    [Fraction(1, 1), Fraction(-2, 1), Fraction(49, 12)],
    [Fraction(1, 1), Fraction(-1, 1), Fraction(13, 12)],
    [Fraction(1, 1), Fraction(0, 1), Fraction(1, 12)]
], dtype=object)

# 打印矩阵 A
print("Matrix A:")
for row in A:
    print("[" + ", ".join(str(frac) for frac in row) + "]")
    
    
# 计算 A * A
product_AA = np.dot(A, A)

# 打印矩阵 A * A
print("Matrix A * A:")
for row in product_AA:
    print("[" + ", ".join(str(frac) for frac in row) + "]")    
    
# 计算矩阵 A 的逆矩阵
A_inv = np.linalg.inv(A).astype(object)

# 打印矩阵 A 的逆矩阵
print("Inverse of Matrix A:")
for row in A_inv:
    print("[" + ", ".join(str(Fraction(int(x.real), int(y))) if y != 0 else str(int(x.real)) for x, y in zip(row, [1, 1, 1])) + "]")    