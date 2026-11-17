import sympy as sp
import numpy as np

# 定义变量
v = sp.symbols('v0:6')  # v0, v1, v2, v3, v4, v5
# 索引映射：
# v[i]   -> v2
# v[i+1] -> v3
# v[i+2] -> v4
# v[i-1] -> v1
# v[i-2] -> v0

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
print("\nA1 =", A1)

# 对 beta2
var2 = [v[0], v[1], v[2]]
A2 = quad_form_to_matrix(beta2_expr, var2)
print("\nA2 =", A2)

# 方法1：使用正交对角化（特征值分解）得到有理数形式
def sos_rational_from_matrix(A, vars):
    """将半正定矩阵A分解为有理数形式的平方和"""
    # 计算特征值和特征向量（符号计算）
    eig_data = A.eigenvects()
    
    sos_terms = []
    for eigval, multiplicity, eigvecs in eig_data:
        if eigval > 0:  # 只取正特征值
            for vec in eigvecs:
                # 将特征向量化为最简整数形式
                vec = sp.Matrix(vec)
                # 找到最小公倍数使得所有系数为整数
                denoms = [sp.fraction(sp.simplify(c))[1] for c in vec if c != 0]
                lcm_denom = 1
                for d in denoms:
                    if d != 1:
                        lcm_denom = sp.lcm(lcm_denom, d)
                
                # 乘以最小公倍数得到整数系数
                vec_int = vec * lcm_denom
                # 提取系数的最大公约数
                coeffs = [vec_int[i] for i in range(vec_int.rows)]
                gcd_coeff = abs(sp.gcd(coeffs))
                if gcd_coeff > 1:
                    vec_int = vec_int / gcd_coeff
                
                # 构造线性组合
                linear_expr = sum(vec_int[i] * vars[i] for i in range(len(vars)))
                
                # 计算系数
                # 验证：vec^T A vec = eigval * vec^T vec
                # 但我们需要的是 eigval * (vec^T x)^2 / (vec^T vec)
                vec_norm2 = sum(c**2 for c in vec)
                coefficient = eigval / vec_norm2 * (lcm_denom / gcd_coeff)**2
                
                # 简化系数为分数形式
                coefficient = sp.nsimplify(coefficient)
                
                sos_terms.append((coefficient, linear_expr))
    
    return sos_terms

# 方法2：使用完成平方的方法
def complete_square(expr, vars):
    """使用完成平方的方法得到SOS表示"""
    # 这是一个启发式方法，尝试常见的线性组合模式
    # 对于WENO格式，通常的形式是 a*(p*v[i] + q*v[i+1] + r*v[i+2])^2 + b*(s*v[i] + t*v[i+1] + u*v[i+2])^2
    
    # 对于beta0，我们知道形式应该是：
    # c1*(v[i] - 2*v[i+1] + v[i+2])^2 + c2*(3*v[i] - 4*v[i+1] + v[i+2])^2
    # 让我们用符号求解系数
    
    if vars == var0:  # beta0
        # 尝试形式: c1*(a*v2 + b*v3 + c*v4)^2 + c2*(d*v2 + e*v3 + f*v4)^2
        c1, c2, a, b, c, d, e, f = sp.symbols('c1 c2 a b c d e f')
        target = c1*(a*v[2] + b*v[3] + c*v[4])**2 + c2*(d*v[2] + e*v[3] + f*v[4])**2
        
        # 展开并比较系数
        target_expanded = sp.expand(target)
        
        # 收集系数方程
        coeff_eqs = []
        # v2^2 系数
        coeff_eqs.append(sp.Eq(target_expanded.coeff(v[2]**2), beta0_expr.coeff(v[2]**2)))
        # v3^2 系数
        coeff_eqs.append(sp.Eq(target_expanded.coeff(v[3]**2), beta0_expr.coeff(v[3]**2)))
        # v4^2 系数
        coeff_eqs.append(sp.Eq(target_expanded.coeff(v[4]**2), beta0_expr.coeff(v[4]**2)))
        # v2*v3 系数
        coeff_eqs.append(sp.Eq(target_expanded.coeff(v[2]*v[3]), beta0_expr.coeff(v[2]*v[3])))
        # v2*v4 系数
        coeff_eqs.append(sp.Eq(target_expanded.coeff(v[2]*v[4]), beta0_expr.coeff(v[2]*v[4])))
        # v3*v4 系数
        coeff_eqs.append(sp.Eq(target_expanded.coeff(v[3]*v[4]), beta0_expr.coeff(v[3]*v[4])))
        
        # 添加归一化约束以减少自由度
        # 让第一个线性组合的系数为简单整数，比如 (1, -2, 1)
        coeff_eqs.append(sp.Eq(a, 1))
        coeff_eqs.append(sp.Eq(b, -2))
        coeff_eqs.append(sp.Eq(c, 1))
        
        # 求解
        sol = sp.nsolve(coeff_eqs, [c1, c2, a, b, c, d, e, f], [13/12, 1/4, 1, -2, 1, 3, -4, 1])
        
        return sol
    
    elif vars == var1:  # beta1
        # 形式: c1*(v1 - 2*v2 + v3)^2 + c2*(v1 - v3)^2
        c1, c2 = sp.symbols('c1 c2')
        target = c1*(v[1] - 2*v[2] + v[3])**2 + c2*(v[1] - v[3])**2
        target_expanded = sp.expand(target)
        
        # 收集系数方程
        coeff_eqs = []
        coeff_eqs.append(sp.Eq(target_expanded.coeff(v[1]**2), beta1_expr.coeff(v[1]**2)))
        coeff_eqs.append(sp.Eq(target_expanded.coeff(v[2]**2), beta1_expr.coeff(v[2]**2)))
        coeff_eqs.append(sp.Eq(target_expanded.coeff(v[3]**2), beta1_expr.coeff(v[3]**2)))
        coeff_eqs.append(sp.Eq(target_expanded.coeff(v[1]*v[2]), beta1_expr.coeff(v[1]*v[2])))
        coeff_eqs.append(sp.Eq(target_expanded.coeff(v[1]*v[3]), beta1_expr.coeff(v[1]*v[3])))
        coeff_eqs.append(sp.Eq(target_expanded.coeff(v[2]*v[3]), beta1_expr.coeff(v[2]*v[3])))
        
        # 求解
        sol = sp.solve(coeff_eqs, [c1, c2])
        return sol
    
    elif vars == var2:  # beta2
        # 形式: c1*(v0 - 2*v1 + v2)^2 + c2*(v0 - 4*v1 + 3*v2)^2
        c1, c2 = sp.symbols('c1 c2')
        target = c1*(v[0] - 2*v[1] + v[2])**2 + c2*(v[0] - 4*v[1] + 3*v[2])**2
        target_expanded = sp.expand(target)
        
        # 收集系数方程
        coeff_eqs = []
        coeff_eqs.append(sp.Eq(target_expanded.coeff(v[0]**2), beta2_expr.coeff(v[0]**2)))
        coeff_eqs.append(sp.Eq(target_expanded.coeff(v[1]**2), beta2_expr.coeff(v[1]**2)))
        coeff_eqs.append(sp.Eq(target_expanded.coeff(v[2]**2), beta2_expr.coeff(v[2]**2)))
        coeff_eqs.append(sp.Eq(target_expanded.coeff(v[0]*v[1]), beta2_expr.coeff(v[0]*v[1])))
        coeff_eqs.append(sp.Eq(target_expanded.coeff(v[0]*v[2]), beta2_expr.coeff(v[0]*v[2])))
        coeff_eqs.append(sp.Eq(target_expanded.coeff(v[1]*v[2]), beta2_expr.coeff(v[1]*v[2])))
        
        # 求解
        sol = sp.solve(coeff_eqs, [c1, c2])
        return sol

print("\n=== 使用完成平方的方法 ===")
print("\n--- Beta0 ---")
sol0 = complete_square(beta0_expr, var0)
if sol0:
    print("系数解:", sol0)
    # 重构表达式
    if isinstance(sol0, list):
        c1_val, c2_val = sol0[0][0], sol0[0][1]
        beta0_reconstructed = c1_val*(v[2] - 2*v[3] + v[4])**2 + c2_val*(3*v[2] - 4*v[3] + v[4])**2
        print("重构的Beta0:", beta0_reconstructed)
        print("是否等于原式:", sp.simplify(beta0_expr - beta0_reconstructed) == 0)

print("\n--- Beta1 ---")
sol1 = complete_square(beta1_expr, var1)
if sol1:
    print("系数解:", sol1)
    # 重构表达式
    c1_val, c2_val = list(sol1.values())[0]
    beta1_reconstructed = c1_val*(v[1] - 2*v[2] + v[3])**2 + c2_val*(v[1] - v[3])**2
    print("重构的Beta1:", beta1_reconstructed)
    print("是否等于原式:", sp.simplify(beta1_expr - beta1_reconstructed) == 0)

print("\n--- Beta2 ---")
sol2 = complete_square(beta2_expr, var2)
if sol2:
    print("系数解:", sol2)
    # 重构表达式
    c1_val, c2_val = list(sol2.values())[0]
    beta2_reconstructed = c1_val*(v[0] - 2*v[1] + v[2])**2 + c2_val*(v[0] - 4*v[1] + 3*v[2])**2
    print("重构的Beta2:", beta2_reconstructed)
    print("是否等于原式:", sp.simplify(beta2_expr - beta2_reconstructed) == 0)

print("\n=== 总结 ===")
print("Beta0 = 13/12 * (v[i] - 2v[i+1] + v[i+2])^2 + 1/4 * (3v[i] - 4v[i+1] + v[i+2])^2")
print("Beta1 = 13/12 * (v[i-1] - 2v[i] + v[i+1])^2 + 1/4 * (v[i-1] - v[i+1])^2")
print("Beta2 = 13/12 * (v[i-2] - 2v[i-1] + v[i])^2 + 1/4 * (v[i-2] - 4v[i-1] + 3v[i])^2")

# 验证
print("\n=== 验证 ===")
beta0_target = sp.Rational(13,12)*(v[2] - 2*v[3] + v[4])**2 + sp.Rational(1,4)*(3*v[2] - 4*v[3] + v[4])**2
print("Beta0 匹配:", sp.simplify(beta0_expr - beta0_target) == 0)

beta1_target = sp.Rational(13,12)*(v[1] - 2*v[2] + v[3])**2 + sp.Rational(1,4)*(v[1] - v[3])**2
print("Beta1 匹配:", sp.simplify(beta1_expr - beta1_target) == 0)

beta2_target = sp.Rational(13,12)*(v[0] - 2*v[1] + v[2])**2 + sp.Rational(1,4)*(v[0] - 4*v[1] + 3*v[2])**2
print("Beta2 匹配:", sp.simplify(beta2_expr - beta2_target) == 0)