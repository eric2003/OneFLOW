from fractions import Fraction

def calculate_crj(r: int, j: int, k: int) -> Fraction:
    """
    计算系数 c_{r,j} 的值
    
    参数:
        r: 公式中的r参数
        j: 求和索引j
        k: 模板阶数
        
    返回:
        Fraction: 系数的有理数值
    """
    result = Fraction(0, 1)
    # 外层求和：m从j+1到k
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
    r"""
    展开公式 v_{i+1/2}^{(r)} = \sum_{j=0}^{k-1} c_{r,j} \bar{v}_{i-r+j}
    其中 c_{r,j} 使用 calculate_crj(r, j, k) 计算的实际值
    
    参数:
        k: 模板阶数
        r: 公式中的r参数
        
    返回:
        str: 展开后的LaTeX字符串
    """
    terms = []
    for j in range(k):
        # 计算实际系数值
        coeff_fraction = calculate_crj(r, j, k)
        
        # 如果系数为0，跳过该项
        if coeff_fraction == 0:
            continue
        
        # 判断符号并取绝对值
        is_negative = (coeff_fraction < 0)
        abs_fraction = abs(coeff_fraction)
        
        # 将Fraction转换为LaTeX格式
        if abs_fraction.denominator == 1:
            # 整数情况
            if abs_fraction.numerator == 1:
                # 系数为±1，省略数字
                coeff_str = ""
            else:
                coeff_str = str(abs_fraction.numerator)
        else:
            # 分数情况
            coeff_str = f"\\frac{{{abs_fraction.numerator}}}{{{abs_fraction.denominator}}}"
        
        # 计算下标偏移量
        offset = j - r
        
        # 生成下标字符串（使用\phantom实现宽度对齐）
        # i      -> i\phantom{+0} （占位但不显示，宽度与i+1相同）
        # i+1    -> i+1
        # i-1    -> i-1
        if offset == 0:
            subscript = r"i\hphantom{+0}"  # 关键修改：占位对齐
        elif offset > 0:
            subscript = f"i+{offset}"
        else:
            subscript = f"i{offset}"  # offset为负时自动带-
        
        # 变量部分
        v_term = r"\bar{v}_{" + subscript + "}"
        
        # 构建最终项字符串（修正符号处理）
        sign_prefix = "-" if is_negative else ""
        
        if coeff_str == "":
            term = sign_prefix + v_term
        else:
            term = sign_prefix + coeff_str + r"\," + v_term
        
        terms.append(term)
    
    # 如果所有项都是0
    if not terms:
        right_side = "0"
    else:
        # 处理符号连接，确保格式美观
        formatted_terms = []
        for i, term in enumerate(terms):
            # 后续项根据符号添加连接符
            if term.startswith('-'):
                formatted_terms.append(' - ' + term[1:])  # 去掉负号，添加分隔符
            else:
                formatted_terms.append(' + ' + term)
        
        right_side = "".join(formatted_terms)
    
    # 左侧：v_{i+1/2}^{(r)}
    left = f"v_{{i+1/2}}^{{({r})}}"
    
    return f"{left} = {right_side}"

def generate_latex_array(k: int) -> str:
    """
    生成给定k值时，r从k-1到-1的所有展开公式
    使用LaTeX的array环境格式化输出
    
    参数:
        k: 模板阶数
        
    返回:
        str: 包含所有r值的LaTeX array字符串
    """
    equations = []
    
    # r从k-1到-1（降序）
    for r in range(k-1, -2, -1):
        eq = expand_formula_with_values(k, r)
        equations.append(eq)
    
    # 构建array环境，左对齐
    array_content = " \\\\\n".join([f"  {eq}" for eq in equations])
    return f"\\begin{{array}}{{l}}\n{array_content}\n\\end{{array}}"

# 主函数示例
if __name__ == "__main__":
    # 示例：k=3 时的多行公式输出
    print("="*60)
    print("k=3 时的展开公式组（r从2到-1）：")
    print("="*60)
    print(generate_latex_array(3))
    
    print("\n" + "="*60 + "\n")
    
    # 示例：k=4 时的多行公式输出
    print("="*60)
    print("k=4 时的展开公式组（r从3到-1）：")
    print("="*60)
    print(generate_latex_array(4))
    
    print("\n" + "="*60 + "\n")
    
    # 示例：k=5 时的多行公式输出
    print("="*60)
    print("k=5 时的展开公式组（r从4到-1）：")
    print("="*60)
    print(generate_latex_array(5))