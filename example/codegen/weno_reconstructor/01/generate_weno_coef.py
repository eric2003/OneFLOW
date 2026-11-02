#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import yaml
from pathlib import Path
import sys

def fraction_str_to_code(frac) -> str:
    """Convert '3/2' or 1 to '3.0/2.0' or '1.0'."""
    if isinstance(frac, str):
        if '/' in frac:
            a, b = frac.split('/')
            return f"{float(a):g}.0/{float(b):g}.0"
        else:
            return f"{float(frac):g}.0"
    else:
        # Handle int/float from YAML
        return f"{float(frac):g}.0"
        
def build_equation_str(idx, coeffs, fraction_converter):
    """
    构建格式正确的方程字符串（避免多余的+号）
    
    参数：
        idx: 当前索引（用于计算变量名v的编号）
        coeffs: 当前索引对应的系数列表（字符串形式）
        fraction_converter: 分数字符串转代码格式的函数
    
    返回：
        格式化后的方程右侧字符串（如 "1/2*v1-1/3*v2+2/5*v3"）
    """
    # 转换系数格式
    code_coeffs = [fraction_converter(c) for c in coeffs]
    
    # 生成带正确符号的项
    terms = []
    for i in range(len(code_coeffs)):
        coeff = code_coeffs[i]
        var_name = f"v{idx + i + 1}"  # 计算变量名（v1/v2/v3...）
        
        if coeff.startswith('-'):
            # 负数项：直接拼接，保留负号
            terms.append(f"{coeff}*{var_name}")
        else:
            # 正数项：添加+号（第一项后续处理）
            terms.append(f"+{coeff}*{var_name}")
    
    # 处理第一项的多余+号
    if terms and terms[0].startswith('+'):
        terms[0] = terms[0][1:]
    
    # 拼接最终表达式
    return "".join(terms)        
    
def generate_weno_function_code(func_name: str, yaml_path: str, weno_order: int, lr: int) -> list:
    """Generate Python function code from YAML coefficient file."""
    with open(yaml_path, 'r', encoding='utf-8') as f:
        data = yaml.safe_load(f)

    if not data:
        raise ValueError(f"YAML file {yaml_path} is empty")
        

    eno_order = weno_order // 2 + 1
    print(f"weno_order,eno_order={weno_order,eno_order}")
    
    k = eno_order  # spatial order = k
    
    scheme = data["weno_schemes"][k]
    print(f"scheme={scheme}")
    
    inputs = scheme["inputs"]
    input_str = ", ".join(inputs)

    lines = []
    lines.append(f"def {func_name}(self, {input_str}):")

    is_reversed = False    
    if not is_reversed:
        # r = k-1, k-2, ..., 0, -1
        r_vals = list(range(k - 1, -2, -1))
    else:
        # r = -1, 0, 1, ..., k-1
        r_vals = list(range(-1, k))
        
   
    print(f"r_vals={r_vals}")
    lines.append(f"    eps = 1e-6")
    
    
    # Compute betas
    betas = scheme["betas"]
    for name, expr in betas.items():
        lines.append(f"    {name} = {expr}")    
    
    # Linear weights
    weights = scheme["weights"]["d"]
    weights_iter = reversed(weights) if lr == 1 else weights
    
    for i, d in enumerate(weights_iter):
        d_code = fraction_str_to_code(d)
        lines.append(f"    d{i} = {d_code}")
            
    for idx in range( eno_order ):
        lines.append(f"    alpha{idx} = d{idx} / (eps + beta{idx})**2")

    alpha = " + ".join(f"alpha{i}" for i in range(eno_order) )
    
    lines.append(f"    alpha = {alpha}")
    for idx in range( eno_order ):
        lines.append(f"    w{idx} = alpha{idx} / alpha")
        
    # 获取 ENO 系数
    eno_data = data["eno"]  #  ENO
    coeffs_list = eno_data[eno_order]
    
    ishift = 1 if lr == 1 else 0
    
    for idx in range( eno_order ):
        coeffs = coeffs_list[idx+ishift]
        
        code_coeffs = [fraction_str_to_code(c) for c in coeffs]
        #print(f"idx,coeffs={idx,coeffs}")
       
        result = "+".join(f"{code_coeffs[i]}*v{idx+i+1}" for i in range(len(coeffs)))
        
        equation_right = build_equation_str(idx, coeffs, fraction_str_to_code)
        #print(f"q{idx} = {result}")
        print(f"q{idx} = {equation_right}")
        r = r_vals[idx]
        
        lines.append(f"    q{idx} = {equation_right}  # r={r}")
        
    f = " + ".join(f"w{i} * q{i}" for i in range(eno_order) )
    lines.append(f"    return {f}")
        
    return lines

def main():
    input_yaml = "eno_weno_coefficients.yaml"
    output_file = "weno_reconstructors.py"
    
    weno_order = 3

    # Write combined output
    full_lines = [
        "# Auto-generated WENO coefficient initializers",
        "# DO NOT EDIT MANUALLY",
        ""
    ]

    for weno_order in [1,3,5]:
        left_lines = generate_weno_function_code(
            f"_reconstruct_weno{weno_order}_left",
            input_yaml,
            weno_order=weno_order,
            lr=-1
        )
        
        right_lines = generate_weno_function_code(
            f"_reconstruct_weno{weno_order}_right",
            input_yaml,
            weno_order=weno_order,
            lr=1
        )        
    
        full_lines.extend(left_lines)
        full_lines.append("")
        full_lines.extend(right_lines)
        full_lines.append("")    

    with open(output_file, 'w', encoding='utf-8') as f:
        f.write("\n".join(full_lines))

    print(f"Generated {output_file} from:")
    print(f"  - {input_yaml}")

if __name__ == "__main__":
    main()