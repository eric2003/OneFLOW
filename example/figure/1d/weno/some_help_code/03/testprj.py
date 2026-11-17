import sympy as sp

class LatexParser:
    def __init__(self):
        self.symbols = {}
    
    def register_symbols(self, symbol_list):
        """预定义符号"""
        for sym in symbol_list:
            self.symbols[sym] = sp.symbols(sym)
    
    def parse(self, latex_str, use_backup=False):
        """
        解析 LaTeX 字符串
        use_backup: 如果主要方法失败，是否使用备选方案
        """
        # 方法1: 尝试 sympy 的 parse_latex
        try:
            expr = parse_latex(latex_str)
            return expr, "sympy_parse_latex"
        except:
            if not use_backup:
                raise
        
        # 方法2: 转换为 sympy 友好格式
        try:
            sympy_friendly = self._latex_to_sympy_friendly(latex_str)
            expr = sp.sympify(sympy_friendly, locals=self.symbols)
            return expr, "converted_sympify"
        except:
            pass
        
        # 方法3: 尝试 latex2sympy2 (如果安装)
        try:
            from latex2sympy2 import latex2sympy
            expr = latex2sympy(latex_str)
            return expr, "latex2sympy2"
        except ImportError:
            print("警告: 未安装 latex2sympy2")
        except Exception:
            pass
        
        raise ValueError(f"无法解析 LaTeX: {latex_str}")
    
    def _latex_to_sympy_friendly(self, latex_str):
        """内部转换方法"""
        # 简化的转换逻辑
        import re
        
        replacements = [
            (r'\\frac{(.*?)}{(.*?)}', r'(\1)/(\2)'),
            (r'\\sqrt{(.*?)}', r'sqrt(\1)'),
            (r'\^', r'**'),
            (r'\\cdot', r'*'),
            (r'\\times', r'*'),
            (r'\\div', r'/'),
        ]
        
        result = latex_str
        for pattern, replacement in replacements:
            result = re.sub(pattern, replacement, result)
        
        # 移除括号
        result = result.replace('{', '(').replace('}', ')')
        
        return result

# 使用示例
parser = LatexParser()
parser.register_symbols(['x', 'y', 'z', 'alpha', 'beta'])

test_cases = [
    (r"x^2 + y^2", "简单幂运算"),
    (r"\frac{x + 1}{y - 2}", "分数"),
    (r"\sin(\alpha) \cdot \cos(\beta)", "三角函数"),
]

for latex_str, description in test_cases:
    try:
        expr, method = parser.parse(latex_str, use_backup=True)
        print(f"{description}:")
        print(f"  LaTeX: {latex_str}")
        print(f"  方法: {method}")
        print(f"  结果: {expr}")
        print()
    except Exception as e:
        print(f"{description} 解析失败: {e}")