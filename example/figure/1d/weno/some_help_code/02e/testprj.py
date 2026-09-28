import json
from sympy import symbols, Rational, expand, latex, sympify

# 定义符号变量（使用LaTeX下标名称）
v_im2, v_im1, v_i, v_ip1, v_ip2 = symbols('v_{i-2} v_{i-1} v_i v_{i+1} v_{i+2}')

def read_and_expand_formulas(file_path):
    """
    从JSON文件中读取公式表达式，展开并输出指定的LaTeX格式。
    
    输入文件格式：JSON对象，键为"beta0", "beta1", "beta2"，值为SymPy兼容的字符串表达式。
    示例文件内容：
    {
        "beta0": "Rational(13, 12) * (v_i - 2 * v_ip1 + v_ip2)**2 + Rational(1, 4) * (3 * v_i - 4 * v_ip1 + v_ip2)**2",
        "beta1": "Rational(13, 12) * (v_im1 - 2 * v_i + v_ip1)**2 + Rational(1, 4) * (v_im1 - v_ip1)**2",
        "beta2": "Rational(13, 12) * (v_im2 - 2 * v_im1 + v_i)**2 + Rational(1, 4) * (v_im2 - 4 * v_im1 + 3 * v_i)**2"
    }
    """
    with open(file_path, 'r') as f:
        formulas = json.load(f)
    
    # 定义局部变量字典，用于sympify解析表达式
    locals_dict = {
        'v_im2': v_im2,
        'v_im1': v_im1,
        'v_i': v_i,
        'v_ip1': v_ip1,
        'v_ip2': v_ip2
    }
    
    results = {}
    for name, expr_str in formulas.items():
        expr = sympify(expr_str, locals=locals_dict)
        expanded = expand(expr)
        results[name] = latex(expanded)
    
    # 输出指定的LaTeX数组格式
    output = r'\begin{array}{l}' + '\n'
    for i, name in enumerate(['beta0', 'beta1', 'beta2']):
        output += f'β{name[4:]} = {results[name]}'
        if i < 2:
            output += r'\\'
        output += '\n'
    output += r'\end{array}'
    print(output)

# 使用示例：替换为实际文件路径
read_and_expand_formulas('formulas.json')