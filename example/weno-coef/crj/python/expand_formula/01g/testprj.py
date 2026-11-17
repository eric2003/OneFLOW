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

def generate_latex_alignat(k: int) -> str:
    """
    生成给定k值时，r从k-1到-1的所有展开公式
    使用LaTeX的alignat环境格式化输出
    
    参数:
        k: 模板阶数
        
    返回:
        str: 包含所有r值的LaTeX alignat字符串
    """
    rows = []
    
    # r从k-1到-1（降序）
    for r in range(k-1, -2, -1):
        terms_info = []
        for j in range(k):
            coeff_fraction = calculate_crj(r, j, k)
            
            if coeff_fraction == 0:
                continue
            
            # 符号和绝对值
            sign = -1 if coeff_fraction < 0 else 1
            abs_coeff = abs(coeff_fraction)
            
            # 生成系数字符串
            if abs_coeff.denominator == 1:
                if abs_coeff.numerator == 1:
                    coeff_str = ""  # 系数为1，省略数字
                else:
                    coeff_str = str(abs_coeff.numerator)
            else:
                coeff_str = f"\\frac{{{abs_coeff.numerator}}}{{{abs_coeff.denominator}}}"
            
            # 生成变量下标
            offset = j - r
            if offset == 0:
                subscript = "i"
            elif offset > 0:
                subscript = f"i+{offset}"
            else:
                subscript = f"i{offset}"
            
            # 变量部分
            v_term = f"\\bar{{v}}_{{{subscript}}}"
            
            terms_info.append({
                'sign': sign,
                'coeff_str': coeff_str,
                'v_term': v_term
            })
        
        # 构建行内容
        if not terms_info:
            row = f"  v_{{i+1/2}}^{{({r})}} &= 0"
        else:
            # 左侧
            left_part = f"  v_{{i+1/2}}^{{({r})}} &="
            
            # 处理每一项
            term_parts = []
            for idx, term_info in enumerate(terms_info):
                # 符号处理
                if idx == 0:
                    # 第一项
                    if term_info['sign'] == 1:
                        sign_str = "\\phantom{+}"
                    else:
                        sign_str = "-"
                else:
                    # 后续项
                    if term_info['sign'] == 1:
                        sign_str = "+"
                    else:
                        sign_str = "-"
                
                # 系数和变量
                coeff_str = term_info['coeff_str']
                v_term = term_info['v_term']
                
                if coeff_str == "":
                    term_str = f"{sign_str} {v_term}"
                else:
                    term_str = f"{sign_str} {coeff_str}{v_term}"
                
                term_parts.append(term_str)
            
            # 用 && 连接各项（第一项前不加&&）
            right_side = term_parts[0]
            for part in term_parts[1:]:
                right_side += f" && {part}"
            
            row = f"{left_part} {right_side}"
        
        rows.append(row)
    
    # 构建alignat环境
    n_columns = k  # 项数
    array_content = " \\\\\n".join(rows)
    return f"\\begin{{alignat}}{{{n_columns}}}\n{array_content}\n\\end{{alignat}}"

# 主函数示例
if __name__ == "__main__":
    # 示例：k=3 时的多行公式输出
    print("="*60)
    print("k=3 时的展开公式组（r从2到-1）：")
    print("="*60)
    print(generate_latex_alignat(3))
    
    print("\n" + "="*60 + "\n")
    
    # 示例：k=4 时的多行公式输出
    print("="*60)
    print("k=4 时的展开公式组（r从3到-1）：")
    print("="*60)
    print(generate_latex_alignat(4))
    
    print("\n" + "="*60 + "\n")
    
    # 示例：k=5 时的多行公式输出
    print("="*60)
    print("k=5 时的展开公式组（r从4到-1）：")
    print("="*60)
    print(generate_latex_alignat(5))