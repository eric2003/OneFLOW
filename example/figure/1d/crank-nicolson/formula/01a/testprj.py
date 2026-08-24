
from typing import List, Union

def build_tridiag_latex(
    size: int,
    param: str = "r",
    show_first: int = 2,
    show_last: int = 2,
    mid_shown: Union[List[int], int, None] = None,
    mid_index: str = "i",
    matrix_name: str = "B",
    vector_name: str = r"\mathbf{u}^n",
    env: str = "pmatrix"
) -> dict:
    """
    动态生成三对角矩阵 B 及其乘积 B·u^n 的 LaTeX 代码。
    
    参数:
        size:       矩阵维度 (N-1)
        param:      参数符号，如 "r", "\\alpha", "\\lambda" 等
        show_first: 开头显示的具体行数
        show_last:  末尾显示的具体行数
        mid_shown:  中间额外显示的行数（仅控制段落是否存在，不影响内容）
        mid_index:  中间通用行的索引符号，如 "i", "j", "k"
        matrix_name: 矩阵名称
        vector_name: 向量名称
        env:        矩阵环境，如 "pmatrix", "bmatrix"
    """
    p = param
    m = size
    
    def fmt_u(sub, sup="n"):
        """智能格式化 u 的下标：纯数字不加花括号，复合表达式加花括号"""
        sub_str = str(sub)
        if sub_str.isdigit():
            return "u_" + sub_str + "^" + sup
        else:
            return "u_{" + sub_str + "}^" + sup
    
    half_p = "\\dfrac{" + p + "}{2}"
    one_minus_p = "(1-" + p + ")"
    
    # 通用中间行
    def generic_row(idx_name):
        return (
            half_p + fmt_u(idx_name + "-1") + " + " +
            one_minus_p + fmt_u(idx_name) + " + " +
            half_p + fmt_u(idx_name + "+1")
        )
    
    # 开头具体数字行（0-based 行号 i）
    def specific_row(i):
        terms = []
        if i > 0:
            terms.append(half_p + fmt_u(i))
        terms.append(one_minus_p + fmt_u(i + 1))
        if i < m - 1:
            terms.append(half_p + fmt_u(i + 2))
        return " + ".join(terms)
    
    # 末尾 N 表达式行（0-based 行号 i）
    def tail_row(i):
        terms = []
        if i > 0:
            prev_1based = i
            terms.append(half_p + fmt_u("N-" + str((m + 1) - prev_1based)))
        center_1based = i + 1
        terms.append(one_minus_p + fmt_u("N-" + str((m + 1) - center_1based)))
        if i < m - 1:
            next_1based = i + 2
            terms.append(half_p + fmt_u("N-" + str((m + 1) - next_1based)))
        return " + ".join(terms)
    
    # 规范化 mid_shown
    if mid_shown is None:
        mid_count = 0
    elif isinstance(mid_shown, int):
        mid_count = 1
    else:
        mid_count = len(mid_shown)
    
    has_head = show_first > 0
    has_mid = mid_count > 0
    has_tail = show_last > 0
    
    # 构建 product 行（按段落组织，段落间插入 \vdots）
    product_lines = []
    
    # Head 段
    if has_head:
        for i in range(min(show_first, m)):
            product_lines.append(specific_row(i))
    
    # Head → Mid
    if has_head and has_mid:
        product_lines.append(r"\vdots")
    
    # Mid 段
    if has_mid:
        for _ in range(mid_count):
            product_lines.append(generic_row(mid_index))
            # Mid 段内部行之间也插入 \vdots（除了最后一行）
            if _ < mid_count - 1:
                product_lines.append(r"\vdots")
    
    # Mid → Tail 或 Head → Tail
    if has_tail and (has_head or has_mid):
        product_lines.append(r"\vdots")
    
    # Tail 段
    if has_tail:
        tail_start = max(m - show_last, 0)
        for i in range(tail_start, m):
            product_lines.append(tail_row(i))
    
    # 构建 B 矩阵行（类似段落结构）
    def b_matrix_row(i):
        row = [""] * m
        row[i] = "1-" + p
        if i > 0:   row[i-1] = half_p
        if i < m-1: row[i+1] = half_p
        return " & ".join(row)
    
    matrix_lines = []
    
    if has_head:
        for i in range(min(show_first, m)):
            matrix_lines.append(b_matrix_row(i))
    
    if has_head and has_mid:
        matrix_lines.append(r"\ddots")
    
    if has_mid:
        for _ in range(mid_count):
            matrix_lines.append(r"\ddots")
            if _ < mid_count - 1:
                matrix_lines.append(r"\ddots")
    
    if has_tail and (has_head or has_mid):
        matrix_lines.append(r"\ddots")
    
    if has_tail:
        tail_start = max(m - show_last, 0)
        for i in range(tail_start, m):
            matrix_lines.append(b_matrix_row(i))
    
    newline_sep = r"\\" + "\n"
    matrix_body = newline_sep.join(matrix_lines)
    product_body = newline_sep.join(product_lines)
    
    matrix_latex = matrix_name + " = \\begin{" + env + "}\n" + matrix_body + "\n\\end{" + env + "}"
    product_latex = matrix_name + vector_name + " = \\begin{" + env + "}\n" + product_body + "\n\\end{" + env + "}"
    
    return {"matrix": matrix_latex, "product": product_latex}


# ========== 示例 1：完全匹配用户期望的格式 ==========
print("=" * 70)
print("示例 1：完全匹配用户期望 (size=6, N=7)")
print("=" * 70)
result = build_tridiag_latex(
    size=6,
    param="r",
    show_first=2,
    show_last=2,
    mid_shown=[3],
    mid_index="i"
)
print(result["product"])
print()

# ========== 示例 2：参数换成 α ==========
print("=" * 70)
print("示例 2：参数换成 α (size=8, N=9)")
print("=" * 70)
result2 = build_tridiag_latex(
    size=8,
    param=r"\alpha",
    show_first=2,
    show_last=2,
    mid_shown=[4],
    mid_index="j"
)
print(result2["product"])
print()

# ========== 示例 3：多个中间行 ==========
print("=" * 70)
print("示例 3：多个中间行 (size=12, N=13)")
print("=" * 70)
result3 = build_tridiag_latex(
    size=12,
    param="r",
    show_first=2,
    show_last=2,
    mid_shown=[4, 7],  # 两个中间行
    mid_index="k"
)
print(result3["product"])
print()

# ========== 示例 4：B 矩阵 ==========
print("=" * 70)
print("示例 4：B 矩阵 (size=6, N=7)")
print("=" * 70)
print(result["matrix"])
