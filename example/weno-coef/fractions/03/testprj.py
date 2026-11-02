import numpy as np
from fractions import Fraction

def inverse_matrix(matrix):
    inverse = np.linalg.inv(matrix)
    inverse_fraction = [[Fraction(inverse[i, j]).limit_denominator() for j in range(len(inverse))] for i in range(len(inverse))]
    return inverse_fraction

matrix = np.array([[1, -2, 49/12], [1, -1, 13/12], [1, 0, 1/12]])
inverse = inverse_matrix(matrix)
for row in inverse:
    print(row)