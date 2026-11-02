import numpy as np

# 创建一个列表来存储所有数组
arrays_list = []

# 动态地添加数组到列表中
arrays_list.append(np.array([[1, 1, 1]]))
arrays_list.append(np.array([[2, 2, 2]]))
arrays_list.append(np.array([[3, 3, 3]]))
# 如果有更多数组，可以继续添加
# arrays_list.append(np.array([[4, 4, 4]]))
# arrays_list.append(np.array([[5, 5, 5]]))

# 使用 vstack 函数将列表中的所有数组垂直堆叠成一个矩阵
matrix = np.vstack(arrays_list)

# 打印结果
print("Assembled Matrix:")
print(matrix)