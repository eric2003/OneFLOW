import math

def divide_as_simplified_fraction(a, b):
    if b == 0:
        raise ValueError("分母不能为零！")
    
    # 计算最大公约数（GCD）
    gcd_value = math.gcd(a, b)
    
    # 化简分数
    simplified_num = a // gcd_value
    simplified_den = b // gcd_value
    
    # 返回最简分数形式
    return f"{simplified_num}/{simplified_den}"

# 示例
a = 10
b = 4
result = divide_as_simplified_fraction(a, b)
print(f"{a} / {b} = {result}")