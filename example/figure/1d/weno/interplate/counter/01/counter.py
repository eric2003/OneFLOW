from collections import Counter

def sort_indices_with_counts(index_list):
    """
    统计下标频次并排序
    
    返回: (排序后的下标列表, 对应的次数列表)
    """
    freq_dict = Counter(index_list)
    sorted_items = sorted(freq_dict.items())
    indices, counts = zip(*sorted_items)  # 解压元组
    return list(indices), list(counts)

# 使用示例
index_list = [0, 1, 3, 2, 5, 1]
indices, counts = sort_indices_with_counts(index_list)

print(f"原始列表: {index_list}")
print(f"排序下标: {indices}")  # [0, 1, 2, 3, 5]
print(f"出现次数: {counts}")   # [1, 2, 1, 1, 1]