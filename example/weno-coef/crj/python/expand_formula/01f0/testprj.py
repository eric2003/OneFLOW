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

def calculate_max_widths(k: int) -> dict:
    """
    计算所有公式中最大的分子、分母宽度
    
    参数:
        k: 模板阶数
        
    返回:
        dict: {'num': 最大分子长度, 'den': 最大分母长度}
    """
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
    
    return {'num': max_num_len, 'den': max_den_len}

def generate_latex_alignat(k: int) -> str:
    """
    生成给定k值时，r从k-1到-1的所有展开公式
    使用LaTeX的alignat环境格式化输出
    
    参数:
        k: 模板阶数
        
    返回:
        str: 包含所有r值的LaTeX alignat字符串
    """
    max_widths = calculate_max_widths(k)
    
    lines = []
    
    for r in range(k-1, -2, -1):
        # 左侧：v_{i+1/2}^{(r)}
        left = f"v_{{i+1/2}}^{{({r})}}"
        # 右侧开始
        right_parts = []
        
        # 收集所有项数据
        terms_data = []
        for j in range(k):
            coeff = calculate_crj(r, j, k)
            if coeff == 0:
                continue
            
            is_negative = (coeff < 0)
            abs_coeff = abs(coeff)
            coeff_str = abs_coeff
            
            # 下标部分
            offset = j - r
            if offset == 0:
                subscript = r"i\hphantom{+0}"
            elif offset > 0:
                subscript = f"i+{offset}"
            else:
                subscript = f"i{offset}"
            
            v_term = r"\bar{v}_{" + subscript + "}"
            terms_data.append((is_negative, coeff_str, v_term))
        
        # 组装右侧表达式
        if not terms_data:
            right_parts.append("0")
        else:
            for idx, (is_negative, coeff_str, v_term) in enumerate(terms_data):
                # 符号处理
                if idx == 0:
                    # 第一项：符号与项在同一单元格
                    sign = "-" if is_negative else r"\phantom{+}"
                    term = f"{sign}{coeff_str}{v_term}" if coeff_str else f"{sign}{v_term}"
                    right_parts.append(term)
                else:
                    # 后续项：符号和项分开对齐
                    sign = "-" if is_negative else "+"
                    term = f"{coeff_str}{v_term}" if coeff_str else v_term
                    right_parts.append(f"&&{sign} {term}")
        
        # 组合成行
        line = f"  {left} &= " + " ".join(right_parts) + r" \\"
        lines.append(line)
    
    # 计算alignat列数：每项需要2列（符号列和项列），但第一项的符号和项在同一列
    max_terms = max(len([j for j in range(k) if calculate_crj(r, j, k) != 0]) for r in range(k-1, -2, -1))
    alignat_cols = 2 * max_terms - 1
    
    content = "\n".join(lines)
    return f"\\begin{{alignat}}{{{alignat_cols}}}\n{content}\n\\end{{alignat}}"

# 主函数示例
if __name__ == "__main__":
    print("="*60)
    print("k=3 时的展开公式组（r从2到-1）：")
    print("="*60)
    print(generate_latex_alignat(3))
    
    print("\n" + "="*60 + "\n")
    
    print("="*60)
    print("k=4 时的展开公式组（r从3到-1）：")
    print("="*60)
    print(generate_latex_alignat(4))