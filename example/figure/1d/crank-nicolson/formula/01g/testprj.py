
from typing import List, Union, Optional
import re


def make_expr(base: str, offset: int) -> str:
    """
    智能构造表达式 base±offset，自动简化。
    支持 base 为 "N", "N+1", "N+2", "M-1" 等形式。
    
    示例:
        make_expr("N", 0)    → "N"
        make_expr("N", 1)    → "N-1"
        make_expr("N+1", 1)  → "N"
        make_expr("N+1", 2)  → "N-1"
        make_expr("M+2", 1)  → "M+1"
    """
    offset = int(offset)
    if offset == 0:
        return base
    
    match = re.match(r'^([A-Za-z]+)([+-]\d+)?$', base)
    if match:
        var = match.group(1)
        base_offset_str = match.group(2)
        base_val = int(base_offset_str) if base_offset_str else 0
        total = base_val - offset
        
        if total == 0:
            return var
        elif total > 0:
            return f"{var}+{total}"
        else:
            return f"{var}{total}"
    else:
        if offset < 0:
            return f"{base}+{abs(offset)}"
        else:
            return f"{base}-{offset}"


def build_tridiag_latex(
    N: Optional[int] = None,
    size: Optional[int] = None,
    param: str = "r",
    show_first: int = 2,
    show_last: int = 2,
    mid_count: int = 1,
    mid_index: str = "i",
    N_expr: Optional[str] = None,
    u_symbol: str = "u",
    superscript: str = "n",
    start_index: int = 0,
    matrix_name: str = "B",
    vector_name: Optional[str] = None,
    env: str = "pmatrix",
    verbose: bool = True
) -> dict:
    """
    动态生成三对角矩阵 B 及其乘积 B·u^n 的 LaTeX 代码。
    支持自定义网格点下标起始值（0 或 1）。
    区分真实代数维数和显示维数。
    
    参数:
        N:           真实网格规模参数（总网格点数 = N + 1，内部点 = N - 1）。
                     若提供，则真实内部点数量 = N - 1。
                     若 N 为 None，则使用 size 参数。
        size:        显示用的内部点数量（向后兼容）。若 N 为 None，则 size 作为真实内部点数量。
        param:       参数符号，如 "r", "\\alpha", "\\lambda" 等
        show_first:  开头显示的具体行数
        show_last:   末尾显示的具体行数
        mid_count:   中间连续通用行数（下标连续递增，内部无省略号）
        mid_index:   中间第一行的中心索引符号，如 "i", "j", "k"
        N_expr:      右端点的 LaTeX 表示。若 None，自动使用 "N"。
        u_symbol:    向量符号，如 "u", "v", "w", "\\phi"
        superscript: 上标，如 "n", "k", "*", "{n+1}"
        start_index: 网格点起始下标，0 或 1
        matrix_name: 矩阵名称
        vector_name: 向量名称。若 None，自动生成为 \\mathbf{u}^n 形式
        env:         矩阵环境，如 "pmatrix", "bmatrix"
        verbose:     是否打印网格信息
    
    返回:
        dict，包含 "matrix", "product", "info" 三个键
    """
    if start_index not in (0, 1):
        raise ValueError("start_index 必须是 0 或 1")
    
    # 确定真实网格规模
    if N is not None and size is not None:
        # 两者都提供，检查一致性
        if size != N - 1:
            raise ValueError(f"N={N} 与 size={size} 不一致: 期望 size={N-1}")
        m = N - 1
    elif N is not None:
        m = N - 1
    elif size is not None:
        m = size
        N = m + 1  # 推导 N
    else:
        raise ValueError("必须提供 N 或 size 之一")
    
    if N_expr is None:
        N_expr = "N"
    
    p = param
    u = u_symbol
    s = superscript
    
    # 边界点
    left_bound = start_index
    right_bound = start_index + m + 1  # = start_index + N
    
    # 内部点列表（真实）
    interior_points = list(range(start_index + 1, right_bound))
    
    # 显示用的内部点（仅用于 info 展示）
    displayed_interior = interior_points[:show_first]
    if mid_count > 0:
        displayed_interior.append("...")
    displayed_interior.extend(interior_points[-show_last:])
    
    # 网格信息
    info = {
        "N": N,
        "start_index": start_index,
        "left_bound": left_bound,
        "right_bound": right_bound,
        "boundary_points": [left_bound, right_bound],
        "interior_points_full": interior_points,  # 完整内部点列表
        "num_interior": m,  # 真实内部点数量
        "num_boundary": 2,
        "total_points": m + 2,  # 总网格点数
        "matrix_dim": f"{m} × {m}",  # 真实矩阵维度
        "show_first": show_first,
        "show_last": show_last,
        "displayed_rows": show_first + show_last + mid_count,
    }
    
    if vector_name is None:
        if u_symbol.startswith("\\"):
            vector_name = f"\\boldsymbol{{{u_symbol}}}^{superscript}"
        else:
            vector_name = f"\\mathbf{{{u_symbol}}}^{superscript}"
    
    if verbose:
        print("=" * 60)
        print("【网格信息】")
        print(f"  起始下标: {start_index}")
        print(f"  左端点: {left_bound}")
        print(f"  右端点: {N_expr} (={right_bound})")
        print(f"  边界点: [{left_bound}, {right_bound}]")
        if m <= 10:
            print(f"  内部点: {interior_points}")
        else:
            print(f"  内部点: [{interior_points[0]}, {interior_points[1]}, ..., {interior_points[-2]}, {interior_points[-1]}]")
        print(f"  真实内部点数量: N-1 (={m})")
        print(f"  显示内部点数量: {show_first + show_last}")
        print(f"  总网格点数: N+1 (={m + 2})")
        print(f"  真实矩阵维度: (N-1) × (N-1) (={m} × {m})")
        print("=" * 60)
    
    # ========== LaTeX 格式化工具 ==========
    def fmt_u(sub, sup=None):
        """格式化 u 的下标和上标"""
        if sup is None:
            sup = s
        sub_str = str(sub)
        sub_part = f"_{sub_str}" if sub_str.isdigit() else f"_{{{sub_str}}}"
        sup_part = f"^{sup}" if sup else ""
        return f"{u}{sub_part}{sup_part}"
    
    half_p = f"\\dfrac{{{p}}}{{2}}"
    one_minus_p = f"(1-{p})"
    
    # ========== B·u^n 向量 ==========
    
    # 开头具体行（第 k 行，0-based，对应内部点 interior_points[k]）
    def head_row(k):
        j = interior_points[k]
        terms = []
        if k > 0:
            terms.append(half_p + fmt_u(j - 1))
        terms.append(one_minus_p + fmt_u(j))
        if k < m - 1:
            terms.append(half_p + fmt_u(j + 1))
        return " + ".join(terms)
    
    # 中间通用行（第 j 行，中心索引递增）
    def mid_row(j):
        if j == 0:
            center, prev, nxt = mid_index, f"{mid_index}-1", f"{mid_index}+1"
        else:
            center = f"{mid_index}+{j}"
            prev = mid_index if j == 1 else f"{mid_index}+{j-1}"
            nxt = f"{mid_index}+{j+1}"
        return (
            half_p + fmt_u(prev) + " + " +
            one_minus_p + fmt_u(center) + " + " +
            half_p + fmt_u(nxt)
        )
    
    # 末尾行（倒数第 j 行，j 从 show_last 到 1）
    def tail_row(j):
        terms = []
        center = right_bound - j
        prev = center - 1
        nxt = center + 1
        
        # 左邻居
        if prev > left_bound:
            if start_index == 0:
                prev_expr = make_expr(N_expr, j + 1)
            else:
                prev_expr = make_expr(N_expr, j)
            terms.append(half_p + fmt_u(prev_expr))
        
        # 中心
        if start_index == 0:
            center_expr = make_expr(N_expr, j)
        else:
            center_expr = make_expr(N_expr, j - 1)
        terms.append(one_minus_p + fmt_u(center_expr))
        
        # 右邻居
        if nxt < right_bound:
            if start_index == 0:
                nxt_expr = make_expr(N_expr, j - 1)
            else:
                nxt_expr = make_expr(N_expr, j - 2)
            terms.append(half_p + fmt_u(nxt_expr))
        
        return " + ".join(terms)
    
    # 构建 product 行
    product_lines = []
    
    # Head 段
    for k in range(min(show_first, m)):
        product_lines.append(head_row(k))
    
    # Head → Mid 省略号
    if show_first > 0 and mid_count > 0:
        product_lines.append(r"\vdots")
    
    # Mid 段（连续行，内部无省略号）
    for j in range(mid_count):
        product_lines.append(mid_row(j))
    
    # Mid → Tail 省略号
    if show_last > 0 and (show_first > 0 or mid_count > 0):
        product_lines.append(r"\vdots")
    
    # Tail 段
    for j in range(show_last, 0, -1):
        product_lines.append(tail_row(j))
    
    # ========== B 矩阵 ==========
    matrix_lines = []
    
    # Head 段
    for i in range(min(show_first, m)):
        elems = []
        if i > 0:
            elems.append(f"\\dfrac{{{p}}}{{2}}")
        elems.append(f"1-{p}")
        if i < m - 1:
            elems.append(f"\\dfrac{{{p}}}{{2}}")
        matrix_lines.append(" & ".join(elems))
    
    # Head → Mid/Tail 省略号
    if (show_first > 0 or show_last > 0) and m > max(show_first, show_last):
        matrix_lines.append(r"\ddots & \ddots & \ddots")
    
    # Tail 段（带缩进）
    for j in range(show_last):
        actual_i = m - show_last + j
        indent = j + 1
        row = [""] * indent
        if actual_i > 0:
            row.append(f"\\dfrac{{{p}}}{{2}}")
        row.append(f"1-{p}")
        if actual_i < m - 1:
            row.append(f"\\dfrac{{{p}}}{{2}}")
        matrix_lines.append(" & ".join(row))
    
    # 组装 LaTeX
    newline = r"\\" + "\n"
    product_body = newline.join(product_lines)
    matrix_body = newline.join(matrix_lines)
    
    product_latex = f"{matrix_name}{vector_name} = \\begin{{{env}}}\n{product_body}\n\\end{{{env}}}"
    matrix_latex = f"{matrix_name} = \\begin{{{env}}}\n{matrix_body}\n\\end{{{env}}}"
    
    return {
        "matrix": matrix_latex,
        "product": product_latex,
        "info": info
    }


# ========== 示例 1：N=7（小例子，start_index=0）==========
print("=" * 70)
print("示例 1：N=7, start_index=0（内部点 1~6）")
print("=" * 70)
result1 = build_tridiag_latex(
    N=7,
    param="r",
    show_first=2,
    show_last=2,
    mid_count=1,
    mid_index="i",
    N_expr="N",
    start_index=0,
    u_symbol="u",
    superscript="n"
)
print(result1["product"])
print()

# ========== 示例 2：N=10000（大例子，start_index=0）==========
print("=" * 70)
print("示例 2：N=10000, start_index=0（大矩阵）")
print("=" * 70)
result2 = build_tridiag_latex(
    N=10000,
    param="r",
    show_first=2,
    show_last=2,
    mid_count=1,
    mid_index="i",
    N_expr="N",
    start_index=0,
    u_symbol="u",
    superscript="n"
)
print(result2["product"])
print()

# ========== 示例 3：N=7, start_index=1（内部点 2~7）==========
print("=" * 70)
print("示例 3：N=7, start_index=1（内部点 2~7）")
print("=" * 70)
result3 = build_tridiag_latex(
    N=7,
    param="r",
    show_first=2,
    show_last=2,
    mid_count=1,
    mid_index="i",
    N_expr="N",
    start_index=1,
    u_symbol="u",
    superscript="n"
)
print(result3["product"])
print()

# ========== 示例 4：N=10000, start_index=1（大矩阵）==========
print("=" * 70)
print("示例 4：N=10000, start_index=1（大矩阵）")
print("=" * 70)
result4 = build_tridiag_latex(
    N=10000,
    param="r",
    show_first=2,
    show_last=2,
    mid_count=1,
    mid_index="i",
    N_expr="N",
    start_index=1,
    u_symbol="u",
    superscript="n"
)
print(result4["product"])
print()

# ========== 示例 5：使用 size 参数（向后兼容）==========
print("=" * 70)
print("示例 5：使用 size 参数（向后兼容，N 自动推导为 7）")
print("=" * 70)
result5 = build_tridiag_latex(
    size=6,
    param="r",
    show_first=2,
    show_last=2,
    mid_count=1,
    mid_index="i",
    N_expr="N",
    start_index=0,
    u_symbol="u",
    superscript="n"
)
print(result5["product"])
print()

# ========== 示例 6：B 矩阵 ==========
print("=" * 70)
print("示例 6：B 矩阵（N=10000）")
print("=" * 70)
print(result2["matrix"])
