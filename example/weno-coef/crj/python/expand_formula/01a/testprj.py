def expand_formula(k: int, r: int) -> str:
    """
    展开公式 v_{i+1/2} = \sum_{j=0}^{k-1} c_{rj} \bar{v}_{i-r+j}
    返回展开后的 LaTeX 字符串
    """
    terms = []
    for j in range(k):
        # 系数：改为 c_{r,j} 形式（添加逗号）
        coeff = f"c_{{{r},{j}}}"
        
        # 计算下标偏移量
        offset = j - r
        
        # 生成下标字符串（无空格，标准写法 i-2、i、i+3 等）
        if offset == 0:
            subscript = "i"
        elif offset > 0:
            subscript = f"i+{offset}"
        else:
            subscript = f"i{offset}"      # offset 为负时自动带 -
        
        # 变量部分
        v_term = r"\bar{v}_{" + subscript + "}"
        
        # 整项：在系数和变量之间加薄空格 \,
        term = coeff + r"\," + v_term
        terms.append(term)
    
    # 左侧：改为 v_{i+1/2}^{(r)} 形式，其中 r 用真实值替换
    left = f"v_{{i+1/2}}^{{({r})}}"
    
    # 用 " + " 连接所有项（自动处理只有一项或首项为负的情况）
    equation = left + r" = " + " + ".join(terms)
    
    return equation


# 示例
print("k=3, r=0 的展开结果：")
print(expand_formula(3, 0))
print("\nk=3, r=1 的展开结果：")
print(expand_formula(3, 1))
print("\nk=4, r=2 的展开结果：")
print(expand_formula(4, 2))