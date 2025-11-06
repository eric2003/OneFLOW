import numpy as np
import matplotlib.pyplot as plt
import os

# 设置目录路径，'.' 表示当前目录
directory = '.'

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