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
    for m in range(j + 1, k + 1):
        numerator = 0
        for l in range(0, k + 1):
            if l == m:
                continue
            product = 1
            for q in range(0, k + 1):
                if q == m or q == l:
                    continue
                product *= (r - q + 1)
            numerator += product
        
        denominator = 1
        for l in range(0, k + 1):
            if l == m:
                continue
            denominator *= (m - l)
        
        result += Fraction(numerator, denominator)
    
    return result

def generate_aligned_fraction(abs_fraction: Fraction, max_num_len: int, max_den_len: int) -> str:
    """
    生成分数字符串，带有对齐填充
    
    参数:
        abs_fraction: 分数的绝对值
        max_num_len: 最大分子字符串长度
        max_den_len: 最大分母字符串长度
        
    返回:
        str: 对齐的LaTeX分数字符串
    """
    if abs_fraction.denominator == 1:
        # 整数情况
        if abs_fraction.numerator == 1:
            return ""  # ±1时省略
        else:
            return str(abs_fraction.numerator)
    
    # 分数情况
    num_str = str(abs_fraction.numerator)
    den_str = str(abs_fraction.denominator)
    
    # 计算需要的填充长度
    num_pad_len = max_num_len - len(num_str)
    den_pad_len = max_den_len - len(den_str)
    
    # 构建填充字符串（使用\hphantom）
    # \hphantom{0} 会创建一个与"0"等宽的空白
    num_padding = r"\hphantom{" + "0" * num_pad_len + "}" if num_pad_len > 0 else ""
    den_padding = r"\hphantom{" + "0" * den_pad_len + "}" if den_pad_len > 0 else ""
    
    # 返回带填充的分数
    return f"\\frac{{{num_padding}{num_str}}}{{{den_padding}{den_str}}}"

def expand_formula_with_values(k: int, r: int, max_widths: dict) -> str:
    r"""
    展开公式 v_{i+1/2}^{(r)} = \sum_{j=0}^{k-1} c_{r,j} \bar{v}_{i-r+j}
    其中 c_{r,j} 使用 calculate_crj(r, j, k) 计算的实际值
    
    参数:
        k: 模板阶数
        r: 公式中的r参数
        max_widths: 包含全局最大宽度的字典 {'num': int, 'den': int}
        
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
        
        # 生成分数字符串（带对齐）
        coeff_str = generate_aligned_fraction(abs_fraction, max_widths['num'], max_widths['den'])
        
        # 计算下标偏移量
        offset = j - r
        
        # 生成下标字符串（使用\hphantom实现宽度对齐）
        if offset == 0:
            subscript = r"i\hphantom{+0}"  # 占位对齐
        elif offset > 0:
            subscript = f"i+{offset}"
        else:
            subscript = f"i{offset}"  # offset为负时自动带-
        
        # 变量部分
        v_term = r"\bar{v}_{" + subscript + "}"
        
        # 构建最终项字符串
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
    
    return f"{left} = {right_side}"

def generate_latex_array(k: int) -> str:
    """
    生成给定k值时，r从k-1到-1的所有展开公式
    使用LaTeX的array环境格式化输出，带全局分数对齐
    
    参数:
        k: 模板阶数
        
    返回:
        str: 包含所有r值的LaTeX array字符串
    """
    # 第一步：计算所有公式的系数，找出全局最大分子分母宽度
    max_num_len = 0
    max_den_len = 0
    
    for r in range(k-1, -2, -1):
        for j in range(k):
            coeff = calculate_crj(r, j, k)
            if coeff == 0:
                continue
            
            if coeff.denominator != 1:
                max_num_len = max(max_num_len, len(str(abs(coeff.numerator))))
                max_den_len = max(max_den_len, len(str(abs(coeff.denominator))))
    
    max_widths = {'num': max_num_len, 'den': max_den_len}
    
    # 第二步：生成所有对齐的公式
    equations = []
    for r in range(k-1, -2, -1):
        eq = expand_formula_with_values(k, r, max_widths)
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