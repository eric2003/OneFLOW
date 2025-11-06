import numpy as np
import matplotlib.pyplot as plt
import os

# 设置目录路径
directory = '../'

# 初始化一个空列表来存储文件名
solution_files = []

# 遍历目录中的所有文件和文件夹
for filename in os.listdir(directory):
    # 检查文件名是否以 "solution" 开头，并且以 ".plt" 结尾
    if filename.startswith("solution") and filename.endswith(".plt"):
        # 将匹配的文件名添加到列表中
        solution_files.append(filename)

# 打印结果
print("Found the following solution files:")
for file in solution_files:
    print(file)
    
    
def read_solution_file(file_name):
    x_list = []
    u_list = []
    
    with open(file_name, newline='') as csvfile:
        icount = 0
        for line in csvfile:
            # 去除首尾空格，按连续空格分割
            row = line.strip().split()    
            x_list.append(row[0])
            u_list.append(row[1])
            icount += 1
    ni = icount
    print("ni=",ni)
    x = np.zeros( ni )
    u = np.zeros( ni )

    for i in range(0, ni):
        x[i] = float(x_list[i])
        u[i] = float(u_list[i])
    return x, u

plt.figure("OneFLOW-CFD Solver", figsize=(10, 6))
plt.xlabel('x')
plt.ylabel('u')

for file in solution_files:
    filename = '../' + file
    x, u = read_solution_file(filename)
    print(f'{x=}')
    print(f'{u=}')

    plt.plot(x, u, marker='o', linestyle='-', color='b', markersize=8, linewidth=2, alpha=0.7, label=f'Numerical (Rusanov)')
 
T  = 0.625
rk = 2
plt.legend()
plt.grid(True, color='gray', linestyle='--', linewidth=0.5, alpha=0.7)    
plt.tight_layout()
plt.show()      
