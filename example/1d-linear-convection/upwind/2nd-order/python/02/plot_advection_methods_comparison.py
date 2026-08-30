import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Circle, Arrow, FancyArrowPatch
import matplotlib.patches as mpatches

# 设置中文字体
plt.rcParams['font.sans-serif'] = ['Microsoft YaHei', 'SimHei', 'Arial Unicode MS']
plt.rcParams['axes.unicode_minus'] = False

# 创建图示
fig = plt.figure(figsize=(16, 12))

# ========== 图1: 一阶迎风格式 ==========
ax1 = plt.subplot(3, 3, 1)
ax1.set_xlim(-1, 3)
ax1.set_ylim(-1, 3)
ax1.axis('off')
ax1.set_title('一阶迎风格式\n$u_i^{n+1} = u_i^n - \\nu (u_i^n - u_{i-1}^n)$', 
              fontsize=12, fontweight='bold')

# 绘制网格点
points_x = [0, 1, 2]
points_y = [1, 1, 1]
for i, (px, py) in enumerate(zip(points_x, points_y)):
    ax1.plot(px, py, 'o', markersize=20, color='lightblue', markeredgecolor='blue', markeredgewidth=2)
    ax1.text(px, py-0.3, f'$u_{i-1}$' if i==0 else f'$u_i$' if i==1 else f'$u_{i+1}$', 
             ha='center', fontsize=12, fontweight='bold')

# 绘制箭头
arrow1 = Arrow(1, 1, -0.8, 0, width=0.15, color='red')
ax1.add_patch(arrow1)
ax1.text(0.5, 1.3, '$-\\nu(u_i - u_{i-1})$', ha='center', fontsize=11, color='red')

# 结果点
ax1.plot(1, 0, 's', markersize=20, color='lightgreen', markeredgecolor='green', markeredgewidth=2)
ax1.text(1, -0.3, '$u_i^{n+1}$', ha='center', fontsize=12, fontweight='bold')
arrow2 = Arrow(1, 1, 0, -0.8, width=0.15, color='green')
ax1.add_patch(arrow2)

# ========== 图2: 二阶TVD格式 ==========
ax2 = plt.subplot(3, 3, 2)
ax2.set_xlim(-1, 3)
ax2.set_ylim(-1, 3)
ax2.axis('off')
ax2.set_title('二阶TVD格式(minmod)\n界面重构法', fontsize=12, fontweight='bold')

# 绘制网格点
for i, px in enumerate([0, 1, 2]):
    ax2.plot(px, 1, 'o', markersize=20, color='lightblue', markeredgecolor='blue', markeredgewidth=2)
    ax2.text(px, 1.3, f'$u_{i-1}$' if i==0 else f'$u_i$' if i==1 else f'$u_{i+1}$', 
             ha='center', fontsize=11)

# 绘制斜率
ax2.plot([0, 1], [1, 1.5], '--', color='purple', linewidth=2)
ax2.plot([1, 2], [1, 1.5], '--', color='purple', linewidth=2)
ax2.text(0.5, 1.6, '$\\sigma_{i-1}$', ha='center', fontsize=10, color='purple')
ax2.text(1.5, 1.6, '$\\sigma_i$', ha='center', fontsize=10, color='purple')

# 界面值
ax2.plot(0.5, 1.25, 's', markersize=15, color='orange', markeredgecolor='red')
ax2.plot(1.5, 1.25, 's', markersize=15, color='orange', markeredgecolor='red')
ax2.text(0.5, 0.9, '$u_{i-1/2}$', ha='center', fontsize=10, color='red')
ax2.text(1.5, 0.9, '$u_{i+1/2}$', ha='center', fontsize=10, color='red')

# 通量差
arrow3 = FancyArrowPatch((1.5, 1.25), (1, 0.3), arrowstyle='->', mutation_scale=20, 
                         linewidth=2, color='red')
ax2.add_patch(arrow3)
arrow4 = FancyArrowPatch((0.5, 1.25), (1, 0.3), arrowstyle='->', mutation_scale=20, 
                         linewidth=2, color='blue')
ax2.add_patch(arrow4)
ax2.text(1, 0, '$u_i^{n+1}$', ha='center', fontsize=12, fontweight='bold', 
         bbox=dict(boxstyle='round', facecolor='lightgreen'))

# ========== 图3: Lax-Wendroff格式 ==========
ax3 = plt.subplot(3, 3, 3)
ax3.set_xlim(-1, 3)
ax3.set_ylim(-1, 3)
ax3.axis('off')
ax3.set_title('Lax-Wendroff中心格式\n二阶精度', fontsize=12, fontweight='bold')

# 绘制网格点
for i, px in enumerate([0, 1, 2]):
    ax3.plot(px, 1, 'o', markersize=20, color='lightblue', markeredgecolor='blue', markeredgewidth=2)
    ax3.text(px, 1.3, f'$u_{i-1}$' if i==0 else f'$u_i$' if i==1 else f'$u_{i+1}$', 
             ha='center', fontsize=11)

# 中心差分
ax3.plot([0, 2], [1, 1], '-', color='purple', linewidth=3, alpha=0.5)
ax3.text(1, 1.5, '中心差分\n$O(\\Delta x^2)$', ha='center', fontsize=10, 
         bbox=dict(boxstyle='round', facecolor='yellow', alpha=0.5))

# 结果
ax3.plot(1, 0, 's', markersize=20, color='lightgreen', markeredgecolor='green', markeredgewidth=2)
ax3.text(1, -0.3, '$u_i^{n+1}$', ha='center', fontsize=12, fontweight='bold')

# 收敛箭头
for px in [0, 2]:
    arrow = Arrow(px, 1, (1-px)*0.3, -0.7, width=0.12, color='purple', alpha=0.6)
    ax3.add_patch(arrow)

# ========== 图4: 周期性边界 ==========
ax4 = plt.subplot(3, 3, 4)
ax4.set_xlim(-1, 6)
ax4.set_ylim(-1, 2)
ax4.axis('off')
ax4.set_title('周期性边界条件处理', fontsize=12, fontweight='bold')

# 绘制周期性网格
nx_demo = 5
for i in range(nx_demo):
    x_pos = i
    ax4.plot(x_pos, 1, 'o', markersize=15, color='lightblue', markeredgecolor='blue', markeredgewidth=2)
    ax4.text(x_pos, 1.3, f'$u_{i}$', ha='center', fontsize=10)
    ax4.text(x_pos, 0.7, str(i), ha='center', fontsize=8, color='gray')

# 特殊标记边界
ax4.plot(0, 1, 'o', markersize=20, color='red', markeredgecolor='darkred', markeredgewidth=3)
ax4.plot(4, 1, 'o', markersize=20, color='red', markeredgecolor='darkred', markeredgewidth=3)
ax4.text(0, 1.6, 'i=0 (左边界)', ha='center', fontsize=10, color='darkred', fontweight='bold')
ax4.text(4, 1.6, 'i=4 (右边界)', ha='center', fontsize=10, color='darkred', fontweight='bold')

# 周期性连接
curve_x = np.linspace(0, 4, 100)
curve_y = 1.5 + 0.3 * np.sin(np.pi * curve_x / 4)
ax4.plot(curve_x, curve_y, '--', color='red', linewidth=2, alpha=0.7)
ax4.text(2, 2.0, '周期性连接: u[0]的左邻居是u[4]', ha='center', fontsize=10, 
         bbox=dict(boxstyle='round', facecolor='pink', alpha=0.5))

# 箭头显示
arrow_left = Arrow(-0.5, 1, -0.3, 0, width=0.1, color='red')
ax4.add_patch(arrow_left)
ax4.text(-0.8, 1.3, '左邻居\nu[i-1]=u[4]', ha='center', fontsize=9, color='red')

arrow_right = Arrow(4.5, 1, 0.3, 0, width=0.1, color='red')
ax4.add_patch(arrow_right)
ax4.text(4.8, 1.3, '右邻居\nu[i+1]=u[0]', ha='center', fontsize=9, color='red')

# ========== 图5: 实际计算示例 ==========
ax5 = plt.subplot(3, 3, 5)
ax5.set_xlim(-0.5, 3.5)
ax5.set_ylim(-0.5, 2.5)
ax5.axis('off')
ax5.set_title('实际计算示例 (i=1, n时刻)', fontsize=12, fontweight='bold')

# 假设值
u_im1, u_i, u_ip1 = 0.5, 1.0, 0.8
courant = 0.5

# 绘制
y_pos = [1.5, 1.0, 0.5]
values = [u_im1, u_i, u_ip1]
labels = ['$u_{i-1}^n$', '$u_i^n$', '$u_{i+1}^n$']

for i, (y, val, label) in enumerate(zip(y_pos, values, labels)):
    x = 1 + i
    ax5.plot(x, y, 'o', markersize=25, color='lightblue', markeredgecolor='blue', markeredgewidth=2)
    ax5.text(x, y, f'{val:.1f}', ha='center', va='center', fontsize=12, fontweight='bold')
    ax5.text(x, y-0.4, label, ha='center', fontsize=10)

# 一阶迎风计算
u_new_1st = u_i - courant * (u_i - u_im1)
ax5.plot(2, -0.2, 's', markersize=20, color='lightgreen', markeredgecolor='green')
ax5.text(2, -0.2, f'{u_new_1st:.1f}', ha='center', va='center', fontsize=12, fontweight='bold')
ax5.text(2, -0.6, '一阶迎风结果', ha='center', fontsize=10, color='green', fontweight='bold')

# 公式
ax5.text(2, 2.2, f'计算: $u_i^{{n+1}} = {u_i} - {courant}({u_i} - {u_im1})$', 
         ha='center', fontsize=10, bbox=dict(boxstyle='round', facecolor='yellow', alpha=0.5))

# ========== 图6: minmod函数详解 ==========
ax6 = plt.subplot(3, 3, 6)
ax6.set_xlim(-2, 2)
ax6.set_ylim(-2, 2)
ax6.axis('off')
ax6.set_title('minmod限幅器函数\n$\\text{minmod}(a,b)$', fontsize=12, fontweight='bold')

# 绘制坐标轴
ax6.plot([-2, 2], [0, 0], 'k-', linewidth=1)
ax6.plot([0, 0], [-2, 2], 'k-', linewidth=1)
ax6.text(2.1, 0, 'a', fontsize=12)
ax6.text(0, 2.1, 'b', fontsize=12)

# 绘制minmod区域
# 第一象限
theta = np.linspace(0, np.pi/2, 100)
r = np.minimum(np.abs(np.cos(theta)), np.abs(np.sin(theta)))
x = r * np.cos(theta)
y = r * np.sin(theta)
ax6.fill(x, y, color='lightblue', alpha=0.5)

# 第三象限
x3 = -r * np.cos(theta)
y3 = -r * np.sin(theta)
ax6.fill(x3, y3, color='lightblue', alpha=0.5)

# 标注
ax6.text(1, 1, 'minmod(a,b)\n= min(|a|,|b|)', ha='center', fontsize=10, 
         bbox=dict(boxstyle='round', facecolor='white'))
ax6.text(-1, -1, 'minmod(a,b)\n= -min(|a|,|b|)', ha='center', fontsize=10,
         bbox=dict(boxstyle='round', facecolor='white'))
ax6.text(1.5, -0.5, 'minmod = 0\n(符号相反)', ha='center', fontsize=10, color='red',
         bbox=dict(boxstyle='round', facecolor='pink'))

# ========== 图7: 数值耗散与色散 ==========
ax7 = plt.subplot(3, 3, 7)
ax7.set_xlim(0, 4)
ax7.set_ylim(-0.5, 1.5)
ax7.set_xlabel('x', fontsize=11)
ax7.set_ylabel('u', fontsize=11)
ax7.set_title('数值耗散(一阶) vs 色散(中心格式)', fontsize=11, fontweight='bold')
ax7.grid(True, alpha=0.3)

x_demo = np.linspace(0, 4, 200)
# 精确解
u_exact_demo = np.zeros_like(x_demo)
u_exact_demo[(x_demo >= 1.5) & (x_demo <= 2.5)] = 1.0

# 一阶迎风(耗散)
u_1st_demo = np.zeros_like(x_demo)
u_1st_demo[(x_demo >= 1.4) & (x_demo <= 2.6)] = 1.0
# 添加平滑
from scipy.ndimage import gaussian_filter1d
u_1st_demo = gaussian_filter1d(u_1st_demo, sigma=2)

# 中心格式(色散)
u_central_demo = np.zeros_like(x_demo)
u_central_demo[(x_demo >= 1.5) & (x_demo <= 2.5)] = 1.0
# 添加振荡
u_central_demo += 0.1 * np.sin(20 * x_demo) * np.exp(-((x_demo-2)**2)/0.1)

ax7.plot(x_demo, u_exact_demo, 'k-', linewidth=2, label='精确解')
ax7.plot(x_demo, u_1st_demo, 'b-', linewidth=2, label='一阶迎风(耗散)')
ax7.plot(x_demo, u_central_demo, 'r-', linewidth=2, label='中心格式(色散)')
ax7.legend(fontsize=9)

# ========== 图8: CFL条件 ==========
ax8 = plt.subplot(3, 3, 8)
ax8.set_xlim(0, 1)
ax8.set_ylim(0, 1)
ax8.axis('off')
ax8.set_title('CFL条件\n$\\nu = c\\Delta t / \\Delta x \\leq 1$', fontsize=12, fontweight='bold')

# 绘制特征线
for i in range(5):
    x_start = i * 0.2
    x_end = x_start + 0.3
    ax8.plot([x_start, x_end], [0.9, 0.1], '-', color='red', linewidth=2, alpha=0.7)
    ax8.text(x_start, 0.95, f'{i}', ha='center', fontsize=10)

ax8.text(0.5, 0.5, '数值依赖域\n必须包含\n物理依赖域', ha='center', fontsize=11,
         bbox=dict(boxstyle='round', facecolor='yellow', alpha=0.5))

# 网格
for i in range(6):
    x = i * 0.2
    ax8.plot([x, x], [0, 1], '--', color='gray', alpha=0.3)
    ax8.plot(x, 0.5, 'o', markersize=8, color='blue', alpha=0.5)

# ========== 图9: 总结表格 ==========
ax9 = plt.subplot(3, 3, 9)
ax9.set_xlim(0, 1)
ax9.set_ylim(0, 1)
ax9.axis('off')
ax9.set_title('方法对比总结', fontsize=12, fontweight='bold')

table_data = [
    ['方法', '精度', '耗散', '色散', '稳定性'],
    ['一阶迎风', '1阶', '强', '弱', '无条件稳定(CFL≤1)'],
    ['二阶TVD', '2阶', '中', '弱', '无条件稳定(CFL≤1)'],
    ['Lax-Wendroff', '2阶', '无', '强', 'CFL≤1'],
]

table = ax9.table(cellText=table_data, cellLoc='center', loc='center', fontsize=9)
table.auto_set_font_size(False)
table.set_fontsize(9)
table.scale(1.2, 1.2)

# 高亮表头
for i in range(5):
    table[(0, i)].set_facecolor('lightgray')
    table[(0, i)].set_text_props(weight='bold')

plt.tight_layout()
plt.savefig('numerical_methods_illustration.png', dpi=300, bbox_inches='tight')
plt.show()

print("图示已生成！包含9个子图详细说明各种数值方法")