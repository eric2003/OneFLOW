def read_matrix_from_file(file_path):
    matrix = []
    with open(file_path, 'r') as file:
        for line in file:
            # 分割每行的数据
            elements = line.strip().split(', ')
            row = []
            for element in elements:
                # 分割分子和分母
                numerator, denominator = element.split('/')
                # 确保分子和分母是浮点数形式
                row.append(f"{float(numerator):.1f}/{float(denominator):.1f}")
            matrix.append(row)
    return matrix

# 示例：读取文件
file_path = 'matrix_data.txt'
matrix = read_matrix_from_file(file_path)

# 输出结果
for i in range(len(matrix)):
    for j in range(len(matrix[i])):
        print(f"coef({i}, {j}) = {matrix[i][j]}")