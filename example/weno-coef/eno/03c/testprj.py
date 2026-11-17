def read_matrix_from_file(file_path):
    matrix = []
    with open(file_path, 'r') as file:
        for line in file:
            # 去掉每行首尾的方括号和空格
            line = line.strip().strip('[]')
            # 分割每行的数据
            elements = line.split(', ')
            row = []
            for element in elements:
                # 分割分子和分母
                numerator, denominator = element.split('/')
                # 确保分子和分母是浮点数形式
                row.append(f"{float(numerator):.1f}/{float(denominator):.1f}")
            matrix.append(row)
    return matrix

def write_matrix_to_file(matrix, output_file_path):
    with open(output_file_path, 'w') as file:
        for i in range(len(matrix)):
            for j in range(len(matrix[i])):
                file.write(f"coef({i}, {j}) = {matrix[i][j]}\n")

# 示例：读取文件
file_path = 'matrix_data.txt'
output_file_path = 'matrix_output.txt'
matrix = read_matrix_from_file(file_path)

# 输出结果到文件
write_matrix_to_file(matrix, output_file_path)

# 打印结果到控制台
for i in range(len(matrix)):
    for j in range(len(matrix[i])):
        print(f"coef({i}, {j}) = {matrix[i][j]}")