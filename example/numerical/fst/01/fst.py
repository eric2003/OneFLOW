import numpy as np

def fft_fst(signal):
    n = len(signal)
    if n % 2 != 0:
        raise ValueError("输入信号的长度必须是偶数")

    if n == 2:
        # 当输入信号长度为2时，直接计算FST结果
        fst_result = np.zeros(2)
        fst_result[0] = signal[0] + signal[1]
        fst_result[1] = signal[0] - signal[1]
        return fst_result

    # 对输入信号进行拆分
    even = signal[0::2]
    odd = signal[1::2]

    # 递归计算FST
    fst_even = fft_fst(even)
    fst_odd = fft_fst(odd)

    # 计算旋转因子
    angles = np.exp(-1j * np.pi * np.arange(n) / n)

    # 奇偶合并
    fst_result = np.zeros(n)
    fst_result[:n//2] = fst_even + angles[:n//2] * fst_odd
    fst_result[n//2:] = fst_even - angles[n//2:] * fst_odd

    return fst_result

# 示例使用
signal = np.array([1, 2, 3, 4, 5, 6, 7, 8])
fst_result = fft_fst(signal)
print("FFT FST结果:", fst_result)