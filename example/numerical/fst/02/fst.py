import numpy as np
from scipy.fft import fft, fftfreq, dst

signal = np.array([1, 2, 3, 4, 5, 6, 7, 8])
spectrum = dst(signal, type=1) / (2 * len(signal))  # 计算FST

print("spectrum:", spectrum)