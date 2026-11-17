import sympy as sp
import re

i_sym, r = sp.symbols('i r', integer=True)

def _extract_index_from_vbar_name(name):
    match = re.search(r'\\overline\{v\}_\{(.+)\}', name)
    if not match:
        return name
    index_latex = match.group(1)
    try:
        # 注意：这里必须传入所有可能符号！
        return sp.simplify(eval(index_latex, {'i': i_sym, 'r': r}))
    except Exception as e:
        print(f"Eval failed for '{index_latex}': {e}")
        return index_latex

def _sort_key_from_index(index_expr):
    if isinstance(index_expr, str):
        return (1, index_expr)  # 字符串放后面
    try:
        offset = sp.simplify(index_expr - i_sym)
        if offset.is_number:
            return (0, float(offset))
        else:
            return (1, str(index_expr))
    except:
        return (1, str(index_expr))

# 测试用例
test_names = [
    r"\overline{v}_{i - 2}",
    r"\overline{v}_{i + 1}",
    r"\overline{v}_{i}",
    r"\overline{v}_{i - 1}",
    r"\overline{v}_{i + 2}",
]

print("Test 1: Pure i + const")
indices = [_extract_index_from_vbar_name(name) for name in test_names]
print("Indices:", indices)
sorted_names = [name for _, name in sorted(zip(indices, test_names), key=lambda x: _sort_key_from_index(x[0]))]
print("Sorted:", sorted_names)

# Test with r substituted (like r=0)
print("\nTest 2: After r=0 substitution (should be i, i+1, i+2)")
names_r0 = [
    r"\overline{v}_{i}",      # j=0
    r"\overline{v}_{i + 1}",  # j=1
    r"\overline{v}_{i + 2}",  # j=2
]
indices2 = [_extract_index_from_vbar_name(name) for name in names_r0]
print("Indices:", indices2)
sorted_names2 = [name for _, name in sorted(zip(indices2, names_r0), key=lambda x: _sort_key_from_index(x[0]))]
print("Sorted:", sorted_names2)