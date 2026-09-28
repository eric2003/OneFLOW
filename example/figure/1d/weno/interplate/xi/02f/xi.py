from fractions import Fraction

def integral_xi(x, j):
    return (x ** (j + 1)) / (j + 1)

data = []
for j in range(-2, 3):
    xia = Fraction(j) - Fraction(1, 2)
    xib = Fraction(j) + Fraction(1, 2)
    row = [j, f"[{xia}, {xib}]"]
    for i in range(5):
        val = integral_xi(xib, i) - integral_xi(xia, i)
        row.append(0 if val == 0 else val)   # 0 显示为整数 0，其他保持 Fraction
    data.append(row)

import pandas as pd
df = pd.DataFrame(data, columns=['j', 'interval'] + [f'i={i}' for i in range(5)])

# 关键：设置 pandas 显示选项 + 直接 print
pd.set_option('display.max_columns', None)
pd.set_option('display.width', None)
pd.set_option('display.colheader_justify', 'center')

print(df.to_string(index=False))