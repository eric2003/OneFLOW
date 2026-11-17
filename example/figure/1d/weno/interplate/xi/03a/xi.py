import numpy as np
from fractions import Fraction

def integral_xi(x, j):
    return (x ** (j + 1)) / (j + 1)

vv = [-1,0,1]
isize = len(vv)
print(f"{'j':>3}  {'interval':>14}  " + "  ".join(f"i={i:^8}" for i in range(isize)))
print("-" * 72)

arrays_list = []
for j in vv:
    xia = Fraction(j) - Fraction(1,2)
    xib = Fraction(j) + Fraction(1,2)
    print(f"{j:>3}  [{xia:>4},{xib:>4}]  ", end="")
    a_list = []
    for i in range( isize ):
        val = integral_xi(xib, i) - integral_xi(xia, i)
        a_list.append(val)
        s = "0" if val == 0 else str(val)
        print(f"{s:^11}", end="")
    print()
    arrays_list.append(a_list)
    
# 使用 vstack 函数将列表中的数组堆叠成一个矩阵
matrix = np.vstack(arrays_list)
print(f'matrix={matrix}')