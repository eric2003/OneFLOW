import numpy as np
from fractions import Fraction

def calculate_determinant(matrix):
    det = np.linalg.det(matrix)
    return Fraction(det).limit_denominator()

def calculate_adjugate(matrix):
    cofactor_matrix = np.linalg.inv(matrix).T
    return cofactor_matrix.tolist()

def inverse_matrix(matrix):
    det = calculate_determinant(matrix)
    if det == 0:
        raise ValueError("Matrix is singular and does not have an inverse.")
    adjugate = calculate_adjugate(matrix)
    inverse = [[Fraction(adjugate[i][j]) / det for j in range(len(adjugate))] for i in range(len(adjugate))]
    return inverse

matrix = np.array([[1, -1, 13/12], [1, 0, 1/12], [1, 1, 13/12]])
inverse = inverse_matrix(matrix)
for row in inverse:
    print(row)