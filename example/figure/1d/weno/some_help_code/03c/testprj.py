from sympy import symbols, Rational, expand, simplify, nsimplify
from sympy.core.numbers import Float

# 1. 定义符号变量
v2, v3, v4 = symbols('v2 v3 v4', real=True)

# 2. 配置参数（控制有理数简洁性）
DECIMAL_PRECISION = 8  # 保留8位小数（可调整，平衡精度和简洁性）
#MAX_DENOMINATOR = 1000000  # 有理数分母最大值（避免超大数字）
MAX_DENOMINATOR = 1000

# 3. 原始系数
# 第一个表达式：0.255002 * (0.64295988*v2 + 0.11435482*v3 - 0.75731471*v4)^2（截断后）
coeff1_out = 0.255002
coeff1_in = [0.642959882291173, 0.114354823786141, -0.757314706077309]

# 第二个表达式：12.744998 * (-0.50325864*v2 + 0.80844891*v3 - 0.30519027*v4)^2（截断后）
coeff2_out = 12.744998
coeff2_in = [-0.503258637711058, 0.808448910533936, -0.305190272822878]

def float_to_simple_rational(num):
    """
    将浮点数转换为简洁的有理数：
    1. 截断浮点精度，消除噪声
    2. 限制分母最大值，避免虚假大数字
    """
    # 步骤1：截断浮点数精度（保留DECIMAL_PRECISION位小数）
    #truncated = round(num, DECIMAL_PRECISION)
    truncated = num
    # 步骤2：转换为有理数，限制分母最大值
    return Rational(truncated).limit_denominator(MAX_DENOMINATOR)

def normalize_linear_coeffs_simple(coeffs):
    """
    优化版：将线性系数转为最小整数（基于简洁有理数）
    """
    # 步骤1：转简洁有理数
    rat_coeffs = [float_to_simple_rational(c) for c in coeffs]
    
    # 步骤2：提取公分母，转为整数
    denoms = [c.denominator for c in rat_coeffs]
    lcm_denom = 1
    for d in denoms:
        lcm_denom = lcm_denom * d // gcd(lcm_denom, d)
    int_coeffs = [c * lcm_denom for c in rat_coeffs]
    
    # 步骤3：提取最大公约数，化为最小整数
    abs_ints = [abs(int(c)) for c in int_coeffs if c != 0]
    gcd_int = abs_ints[0] if abs_ints else 1
    for num in abs_ints[1:]:
        gcd_int = gcd(gcd_int, num)
    
    # 最小整数系数 + 提取的公因子
    min_coeffs = [int(c / gcd_int) for c in int_coeffs]  # 确保是整数类型
    factor_out = Rational(gcd_int, lcm_denom)
    
    return min_coeffs, factor_out

# 辅助函数：GCD和LCM
def gcd(a, b):
    while b:
        a, b = b, a % b
    return a

# ------------------------ 处理第一个表达式 ------------------------
min_coeffs1, factor1 = normalize_linear_coeffs_simple(coeff1_in)
rat_out1 = float_to_simple_rational(coeff1_out)
linear1 = min_coeffs1[0]*v2 + min_coeffs1[1]*v3 + min_coeffs1[2]*v4
expr1 = rat_out1 * (factor1 **2) * (linear1** 2)
expr1_simplified = simplify(expr1)

# ------------------------ 处理第二个表达式 ------------------------
min_coeffs2, factor2 = normalize_linear_coeffs_simple(coeff2_in)
rat_out2 = float_to_simple_rational(coeff2_out)
linear2 = min_coeffs2[0]*v2 + min_coeffs2[1]*v3 + min_coeffs2[2]*v4
expr2 = rat_out2 * (factor2 **2) * (linear2** 2)
expr2_simplified = simplify(expr2)

# ------------------------ 输出结果 ------------------------
print(f"配置：保留{DECIMAL_PRECISION}位小数，分母最大为{MAX_DENOMINATOR}")
print("\n=== 第一个表达式处理结果 ===")
print(f"平方内最小整数线性组合：{linear1}")
print(f"平方外最简有理系数：{expr1_simplified.coeff(linear1**2)}")
print(f"完整等价表达式：{expr1_simplified}")

print("\n=== 第二个表达式处理结果 ===")
print(f"平方内最小整数线性组合：{linear2}")
print(f"平方外最简有理系数：{expr2_simplified.coeff(linear2**2)}")
print(f"完整等价表达式：{expr2_simplified}")

# ------------------------ 等价性验证 ------------------------
def verify_equivalence(expr_sym, coeff_out, coeff_in, v_vals=[1,1,1]):
    """验证符号表达式与原始浮点表达式的数值等价性"""
    # 原始浮点值
    linear_float = sum(c*v for c, v in zip(coeff_in, v_vals))
    float_val = coeff_out * (linear_float **2)
    # 符号表达式值
    sym_val = expr_sym.subs({v2:v_vals[0], v3:v_vals[1], v4:v_vals[2]}).evalf()
    # 误差阈值（根据截断精度调整）
    error = abs(float_val - sym_val)
    return error < 10**(-DECIMAL_PRECISION)

print("\n=== 等价性验证 ===")
print(f"第一个表达式（v2=1,v3=1,v4=1）：误差={verify_equivalence(expr1_simplified, coeff1_out, coeff1_in, [1,1,1])}")
print(f"第二个表达式（v2=1,v3=1,v4=1）：误差={verify_equivalence(expr2_simplified, coeff2_out, coeff2_in, [1,1,1])}")