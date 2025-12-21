from fractions import Fraction

def integral_xi(x, j):
    return (x ** (j + 1)) / (j + 1)

def calxi(j):
    half = Fraction(1, 2)
    xia = j - half
    xib = j + half
    return xia, xib

jst = -2
jed = 2

# 先计算所有值，找出每列需要的最大宽度（可选，更极致对齐）
# 这里直接用固定宽度，也已经非常整齐

print(f"{'j':>3}  {'interval':>14}  i=0      i=1      i=2      i=3       i=4")
print("-" * 68)

for j in range(jst, jed + 1):
    xia, xib = calxi(j)
    # 把区间格式化为固定宽度字符串
    interval_str = f"[{xia:>4},{xib:>4}]"
    
    line = f"{j:>2}   {interval_str}  "
    
    for i in range(0, 5):  # i 从 0 到 4
        v1 = integral_xi(xia, i)
        v2 = integral_xi(xib, i)
        diff = v2 - v1
        
        if diff == 0:
            s = "0"
        elif diff == 1:
            s = "1"
        elif diff == -1:
            s = "-1"
        else:
            s = str(diff)
        
        # 每列固定 11 个字符宽度，居中对齐（足够容纳 1441/80 这类最长的）
        line += f"{s:^11}"
    
    print(line)