from fractions import Fraction

def integral_xi(x, j):
    return (x ** (j + 1)) / (j + 1)

def calxi(j):
    half = Fraction(1, 2)
    return j - half, j + half

jst, jed = -2, 2

# Step 1: 生成所有行数据（先单独处理 interval 的三部分：左括号、分数值、右括号）
rows = []
left_brackets  = []   # 全部是 "["
right_brackets = []   # 全部是 "]"
intervals_a    = []   # xia
intervals_b    = []   # xib

for j in range(jst, jed + 1):
    xia, xib = calxi(j)
    row = [str(j)]
    for i in range(5):
        diff = integral_xi(xib, i) - integral_xi(xia, i)
        s = "0" if diff == 0 else ("1" if diff == 1 else ("-1" if diff == -1 else str(diff)))
        row.append(s)
    rows.append(row)
    
    left_brackets.append("[")
    intervals_a.append(str(xia))
    intervals_b.append(str(xib))
    right_brackets.append("]")

# Step 2: 标题（j 单独一列，interval 拆成 [ | 值 | ] 三列）
headers = ["j", "", "interval", "", "", "i=0", "i=1", "i=2", "i=3", "i=4"]
# 对应列： 0:j   1:[   2:xia   3:xib   4:]   5~9:i=0~i=4

# Step 3: 计算每一列最大宽度
col_widths = [0] * len(headers)

# j 列 + 数值列
for i, h in enumerate(headers):
    col_widths[i] = max(col_widths[i], len(h))

for row, a, b in zip(rows, intervals_a, intervals_b):
    col_widths[0] = max(col_widths[0], len(row[0]))      # j
    col_widths[2] = max(col_widths[2], len(a))           # xia
    col_widths[3] = max(col_widths[3], len(b))           # xib
    for k in range(5):                                   # i=0~4
        col_widths[5 + k] = max(col_widths[5 + k], len(row[1 + k]))

# 括号列固定宽度 1 就够了，但我们也算进去
col_widths[1] = max(col_widths[1], len("["))   # 1
col_widths[4] = max(col_widths[4], len("]"))   # 1

# 为了美观，给每列再 +1 空格间隔（除了最右侧可以不加）
widths_with_space = [w + 1 for w in col_widths]

# Step 4: 打印函数
def print_line(parts, aligns=None):
    if aligns is None:
        aligns = ['^'] * len(parts)       # 默认居中
    line = ""
    for p, w, align in zip(parts, widths_with_space, aligns):
        line += f"{p:{align}{w}}"
    print(line.rstrip())   # 去掉行末多余空格

# 表头（interval 跨三列，我们手动合并显示）
print_line(headers[:5] + headers[5:], aligns=['^','^','^','^','^','^','^','^','^','^'])
print_line(["", "interval","","","","i=0","i=1","i=2","i=3","i=4"], 
           aligns=['<','^','^','^','>','^','^','^','^','^'])  # interval 标题居中跨三列
print("-" * (sum(widths_with_space) - 1))

# 数据行
for j_row, lb, a, b, rb, data_row in zip(rows, left_brackets, intervals_a, intervals_b, right_brackets, rows):
    parts = [data_row[0], lb, a, b, rb] + data_row[1:]
    print_line(parts)
