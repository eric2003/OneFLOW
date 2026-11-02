import sympy as sp
import numpy as np

# 定义变量
v = sp.symbols('v0:6')  # v0, v1, v2, v3, v4, v5
# 为了对应原表达式，设定索引映射：
# v[i]   -> v2
# v[i+1] -> v3
# v[i+2] -> v4
# v[i-1] -> v1
# v[i-2] -> v0
# 我们用通用符号，之后替换

i = 2  # 中间索引为2，则 v[i]=v2, v[i+1]=v3, v[i+2]=v4, v[i-1]=v1, v[i-2]=v0

# 定义原表达式
beta0_expr = 10*v[2]**2/3 - 31*v[2]*v[3]/3 + 11*v[2]*v[4]/3 + 25*v[3]**2/3 - 19*v[3]*v[4]/3 + 4*v[4]**2/3
beta1_expr = 13*v[2]**2/3 - 13*v[2]*v[1]/3 - 13*v[2]*v[3]/3 + 4*v[1]**2/3 + 5*v[1]*v[3]/3 + 4*v[3]**2/3
beta2_expr = 10*v[2]**2/3 + 11*v[2]*v[0]/3 - 31*v[2]*v[1]/3 + 4*v[0]**2/3 - 19*v[0]*v[1]/3 + 25*v[1]**2/3

# 将二次型转换为矩阵形式
def quad_form_to_matrix(expr, var_list):
    """返回二次型 expr 关于变量 var_list 的对称矩阵 A"""
    n = len(var_list)
    A = sp.zeros(n, n)
    for i in range(n):
        for j in range(i, n):
            coeff = expr.coeff(var_list[i]*var_list[j])
            if i == j:
                A[i,j] = coeff
            else:
                # 交叉项系数，二次型中 x_i x_j 系数对应矩阵 (i,j) 和 (j,i) 各一半
                A[i,j] = coeff/2
                A[j,i] = coeff/2
    return A

# 对 beta0
var0 = [v[2], v[3], v[4]]
A0 = quad_form_to_matrix(beta0_expr, var0)
print("A0 =", A0)

# 对 beta1
var1 = [v[1], v[2], v[3]]
A1 = quad_form_to_matrix(beta1_expr, var1)
print("A1 =", A1)

# 对 beta2
var2 = [v[0], v[1], v[2]]
A2 = quad_form_to_matrix(beta2_expr, var2)
print("A2 =", A2)

# 特征值分解找平方和
def sos_from_matrix(A, vars):
    """将半正定矩阵A分解为 sum_i (linear_form)^2"""
    # 对称矩阵的特征值分解
    A_np = np.array(A, dtype=float)
    eigvals, eigvecs = np.linalg.eigh(A_np)
    # 只取正特征值
    sos_terms = []
    for val, vec in zip(eigvals, eigvecs.T):
        if abs(val) > 1e-10:
            linear = sum(float(c)*x for c, x in zip(vec, vars))
            sos_terms.append((val, linear))
    return sos_terms

print("\n--- Beta0 SOS ---")
sos0 = sos_from_matrix(A0, var0)
for val, lin in sos0:
    print(f"{val:.6f} * ({lin})^2")

print("\n--- Beta1 SOS ---")
sos1 = sos_from_matrix(A1, var1)
for val, lin in sos1:
    print(f"{val:.6f} * ({lin})^2")

print("\n--- Beta2 SOS ---")
sos2 = sos_from_matrix(A2, var2)
for val, lin in sos2:
    print(f"{val:.6f} * ({lin})^2")

# 也可以尝试用符号 Cholesky 分解（但要求正定，这里可能半正定）
# 用 sympy 的 LDL 分解（对半正定有效）
def sos_ldl(A, vars):
    """LDL^T 分解得到平方和"""
    L, D = A.LDLdecomposition()
    # A = L * D * L.T
    # 那么 x^T A x = (sqrt(D) L^T x)^T (sqrt(D) L^T x)
    n = A.rows
    terms = []
    for j in range(n):
        if D[j, j] != 0:
            linear = sum(L[i, j]*vars[i] for i in range(n))
            terms.append((D[j, j], linear))
    return terms

print("\n--- 使用 LDL 分解 ---")
print("Beta0:")
terms0 = sos_ldl(A0, var0)
for coeff, lin in terms0:
    print(f"{coeff} * ({lin})^2")
print("Beta1:")
terms1 = sos_ldl(A1, var1)
for coeff, lin in terms1:
    print(f"{coeff} * ({lin})^2")
print("Beta2:")
terms2 = sos_ldl(A2, var2)
for coeff, lin in terms2:
    print(f"{coeff} * ({lin})^2")