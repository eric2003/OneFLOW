import matplotlib.pyplot as plt

fig, ax = plt.subplots(figsize=(10, 8))

# 不同样式的向上箭头示例
arrow_styles = [
    ('简单箭头', (0.2, 0.2), '->', 'blue'),
    ('燕尾箭头', (0.4, 0.2), '<->', 'red'),
    ('填充箭头', (0.6, 0.2), '-|>', 'green'),
    ('宽头箭头', (0.8, 0.2), 'wedge,tail_width=0.7', 'purple')
]

for label, start_pos, style, color in arrow_styles:
    x, y = start_pos
    ax.annotate(label, xy=(x, y+0.6), xytext=(x, y),
                arrowprops=dict(arrowstyle=style,
                               color=color,
                               lw=2,
                               shrinkA=5, shrinkB=5))
    ax.text(x-0.08, y-0.05, f'起点: ({x},{y})', fontsize=8)

ax.set_xlim(0, 1)
ax.set_ylim(0, 1)
plt.title('不同样式的向上箭头')
plt.grid(True, alpha=0.3)
plt.show()