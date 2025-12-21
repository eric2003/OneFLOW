from fractions import Fraction

def calculate_crj(r: int, j: int, k: int) -> Fraction:
    """
    计算系数 c_{r,j} 的值
    
    参数:
        r: 公式中的 r 参数
        j: 求和索引 j
        k: 模板阶数 k
        
    返回:
        Fraction: 系数的有理数值
    """
    result = Fraction(0, 1)
    # 外层求和：m 从 j+1 到 k
    for m in range(j + 1, k + 1):
        # 计算分子部分
        numerator = 0
        for l in range(0, k + 1):
            if l == m:
                continue
            product = 1
            # 计算连乘积
            for q in range(0, k + 1):
                if q == m or q == l:
                    continue
                product *= (r - q + 1)
            numerator += product
        
        # 计算分母部分
        denominator = 1
        for l in range(0, k + 1):
            if l == m:
                continue
            denominator *= (m - l)
        
        # 累加当前项
        result += Fraction(numerator, denominator)
    
    return result

def expand_formula_with_values(k: int, r: int) -> str:
    """
    展开公式 v_{i+1/2}^{(r)} = \sum_{j=0}^{k-1} c_{r,j} \bar{v}_{i-r+j}
    其中 c_{r,j} 使用 calculate_crj(r, j, k) 计算的实际值
    
    参数:
        k: 模板阶数
        r: 公式中的 r 参数
        
    返回:
        str: 展开后的 LaTeX 字符串
    """
    terms = []
    for j in range(k):
        # 计算实际系数值
        coeff_fraction = calculate_crj(r, j, k)
        
        # 如果系数为 0，跳过该项
        if coeff_fraction == 0:
            continue
        
        # 判断符号并取绝对值
        is_negative = (coeff_fraction < 0)
        abs_fraction = abs(coeff_fraction)
        
        # 将 Fraction 转换为 LaTeX 格式
        if abs_fraction.denominator == 1:
            # 整数情况
            if abs_fraction.numerator == 1:
                # 系数为 ±1，省略数字
                coeff_str = ""
            else:
                coeff_str = str(abs_fraction.numerator)
        else:
            # 分数情况
            coeff_str = f"\\frac{{{abs_fraction.numerator}}}{{{abs_fraction.denominator}}}"
        
        # 计算下标偏移量
        offset = j - r
        
        # 生成下标字符串
        if offset == 0:
            subscript = "i"
        elif offset > 0:
            subscript = f"i+{offset}"
        else:
            subscript = f"i{offset}"  # offset 为负时自动带 -
        
        # 变量部分
        v_term = r"\bar{v}_{" + subscript + "}"
        
        # 构建最终项字符串
        if coeff_str == "":
            # 系数为 ±1
            if is_negative:
                term = "-" + v_term
            else:
                term = v_term
        else:
            # 系数不为 ±1
            term = coeff_str + r"\," + v_term
        
        terms.append(term)
    
    # 如果所有项都是 0
    if not terms:
        right_side = "0"
    else:
        # 处理符号连接，确保格式美观
        formatted_terms = []
        for i, term in enumerate(terms):
            if i == 0:
                # 第一项直接添加
                formatted_terms.append(term)
            else:
                # 后续项根据符号添加连接符
                if term.startswith('-'):
                    formatted_terms.append(' - ' + term[1:])  # 去掉负号
                else:
                    formatted_terms.append(' + ' + term)
        
        right_side = "".join(formatted_terms)
    
    # 左侧：v_{i+1/2}^{(r)}
    left = f"v_{{i+1/2}}^{{({r})}}"
    
    # 组合方程
    equation = left + r" = " + right_side
    
    return equation

# 示例
if __name__ == "__main__":
    print("k=3, r=0 的展开结果：")
    print(expand_formula_with_values(3, 0))
    print("\nk=3, r=1 的展开结果：")
    print(expand_formula_with_values(3, 1))
    print("\nk=4, r=2 的展开结果：")
    print(expand_formula_with_values(4, 2))