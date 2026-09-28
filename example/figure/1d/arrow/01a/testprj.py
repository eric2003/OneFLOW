import matplotlib.pyplot as plt

fig, ax = plt.subplots(figsize=(6, 6))

# 使用 annotate 绘制箭头
ax.annotate('', xy=(0.5, 0.8), xytext=(0.5, 0.2),
            arrowprops=dict(arrowstyle='->', 
                           color='red',
                           lw=3,
                           shrinkA=0, shrinkB=0))

# 或者使用 fancy 箭头样式
ax.annotate('向上箭头', xy=(0.5, 0.8), xytext=(0.5, 0.2),
            arrowprops=dict(arrowstyle='fancy',
                           fc='green', ec='darkgreen',
                           connectionstyle="arc3,rad=0"))

ax.set_xlim(0, 1)
ax.set_ylim(0, 1)
plt.title('向上的箭头 - plt.annotate()')
plt.grid(True, alpha=0.3)
plt.show()