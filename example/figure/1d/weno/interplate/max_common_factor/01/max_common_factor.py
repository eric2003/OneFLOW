from fractions import Fraction
from math import gcd
from functools import reduce

def extract_max_common_factor(numbers):
    """
    从数值列表中提取最大的公共因子，使得剩余部分为互质整数列表
    
    参数:
        numbers: 数值列表（可包含整数、浮点数、字符串分数等）
        
    返回:
        tuple: (factor, simplified_list)
               factor: Fraction类型，提取的公共因子
               simplified_list: 整数列表，最简形式（gcd=1）
    """
    # 1. 将所有输入转换为Fraction，确保精确的有理数运算
    #    Fraction(0.5) = 1/2, Fraction(-1) = -1/1, Fraction("1/3") = 1/3
    fractions = [Fraction(x) for x in numbers]
    
    # 2. 特殊情况处理
    if not fractions:  # 空列表
        return Fraction(1, 1), []
    
    if all(f == 0 for f in fractions):  # 全零列表
        return Fraction(1, 1), [0] * len(fractions)
    
    # 3. 提取所有分子和分母
    numerators = [f.numerator for f in fractions]
    denominators = [f.denominator for f in fractions]
    
    # 4. 计算分子的最大公约数(GCD)
    #    reduce函数连续应用gcd: gcd(gcd(p1,p2), p3), ...
    numerator_gcd = reduce(gcd, numerators)
    
    # 5. 计算分母的最小公倍数(LCM)
    def lcm(a, b):
        """计算两个数的最小公倍数：lcm(a,b) = |a×b| / gcd(a,b)"""
        if a == 0 or b == 0:
            return 0
        return abs(a * b) // gcd(a, b)
    
    #    连续应用lcm: lcm(lcm(q1,q2), q3), ...
    denominator_lcm = reduce(lcm, denominators)
    
    # 6. 最大公共因子 = 分子GCD / 分母LCM
    common_factor = Fraction(numerator_gcd, denominator_lcm)
    
    # 7. 计算简化后的列表
    simplified_fractions = [f / common_factor for f in fractions]
    
    # 8. 验证并转换为整数列表
    #    理论上所有分母都应为1，因为除以了最大公共因子
    simplified_integers = [sf.numerator for sf in simplified_fractions]
    
    # 9. 额外验证：确保结果整数列表互质（gcd=1）
    verification_gcd = reduce(gcd, simplified_integers)
    
    #    如果verification_gcd≠1，说明还能继续提取，调整结果
    if verification_gcd != 1:
        true_factor = common_factor * verification_gcd
        simplified_integers = [x // verification_gcd for x in simplified_integers]
        return true_factor, simplified_integers
    
    return common_factor, simplified_integers


# 测试函数
def demonstrate():
    """演示各种情况"""
    test_cases = [
        ("用户示例", [Fraction(1,2), -1, Fraction(1,2)]),
        ("整数提取", [2, -4, 2]),
        ("分数列表", [Fraction(1,3), Fraction(2,3), -1]),
        ("浮点数", [0.5, -1.0, 0.5]),
        ("互质整数", [1, 2, 3]),
        ("负分数", [Fraction(-2,3), Fraction(4,3), -2]),
    ]
    
    for name, case in test_cases:
        factor, simplified = extract_max_common_factor(case)
        print(f"【{name}】")
        print(f"  原始列表: {case}")
        print(f"  提取系数: {factor} (小数形式: {float(factor)})")
        print(f"  简化列表: {simplified}")
        print(f"  验证除法: {[Fraction(case[i]) / factor for i in range(len(case))]}")
        print(f"  简化列表gcd: {reduce(gcd, simplified)}\n")

# 运行演示
demonstrate()