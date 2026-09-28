"""
文件名: 05_interactive_cfd_plot.py
功能: 创建交互式CFD可视化
包含: 使用Plotly创建可交互的CFD图像
注意: 需要安装plotly库: pip install plotly
"""

import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import plotly.io as pio

def create_interactive_cfd_plot():
    """创建交互式CFD可视化"""
    
    print("正在创建交互式CFD可视化...")
    print("这将在浏览器中打开一个交互式图表。")
    
    # 创建数据
    x = np.linspace(0, 10, 100)
    
    # 精确解和数值解
    u_exact = np.sin(x * 0.8) * np.exp(-0.1*x)
    
    # 添加不同数值格式的"数值解"
    np.random.seed(42)
    u_upwind = u_exact + 0.08*np.random.randn(len(x))
    u_central = u_exact + 0.05*np.sin(5*x)*0.3  # 模拟振荡
    u_quick = u_exact + 0.02*np.random.randn(len(x))
    
    # 创建子图
    fig = make_subplots(
        rows=2, cols=2,
        subplot_titles=('Upwind Scheme (1st Order)', 
                       'Central Difference Scheme (2nd Order)',
                       'QUICK Scheme (3rd Order)',
                       'Grid Point Stencil Dependencies'),
        vertical_spacing=0.12,
        horizontal_spacing=0.1
    )
    
    # 1. 迎风格式
    fig.add_trace(
        go.Scatter(
            x=x, y=u_exact,
            mode='lines',
            name='Exact Solution',
            line=dict(color='black', width=3, dash='solid'),
            showlegend=True,
            legendgroup="exact"
        ),
        row=1, col=1
    )
    
    fig.add_trace(
        go.Scatter(
            x=x, y=u_upwind,
            mode='lines',
            name='Upwind Scheme',
            line=dict(color='red', width=2, dash='dash'),
            showlegend=True,
            legendgroup="upwind"
        ),
        row=1, col=1
    )
    
    # 2. 中心差分格式
    fig.add_trace(
        go.Scatter(
            x=x, y=u_exact,
            mode='lines',
            name='Exact Solution',
            line=dict(color='black', width=3, dash='solid'),
            showlegend=False,
            legendgroup="exact"
        ),
        row=1, col=2
    )
    
    fig.add_trace(
        go.Scatter(
            x=x, y=u_central,
            mode='lines',
            name='Central Difference',
            line=dict(color='blue', width=2, dash='dash'),
            showlegend=True,
            legendgroup="central"
        ),
        row=1, col=2
    )
    
    # 3. QUICK格式
    fig.add_trace(
        go.Scatter(
            x=x, y=u_exact,
            mode='lines',
            name='Exact Solution',
            line=dict(color='black', width=3, dash='solid'),
            showlegend=False,
            legendgroup="exact"
        ),
        row=2, col=1
    )
    
    fig.add_trace(
        go.Scatter(
            x=x, y=u_quick,
            mode='lines',
            name='QUICK Scheme',
            line=dict(color='green', width=2, dash='dash'),
            showlegend=True,
            legendgroup="quick"
        ),
        row=2, col=1
    )
    
    # 4. 网格点模板依赖关系
    # 创建网格点
    grid_points = np.array([0, 1, 2, 3, 4])
    point_names = ['u_{i-2}', 'u_{i-1}', 'u_i', 'u_{i+1}', 'u_{i+2}']
    
    # 添加网格点
    fig.add_trace(
        go.Scatter(
            x=grid_points,
            y=[0, 0, 0, 0, 0],
            mode='markers+text',
            name='Grid Points',
            marker=dict(size=15, color='gray'),
            text=point_names,
            textposition="top center",
            showlegend=False
        ),
        row=2, col=2
    )
    
    # 添加箭头表示依赖关系
    # 迎风格式箭头（从i-1到i）
    fig.add_annotation(
        x=1, y=0,
        ax=2, ay=0,
        xref="x4", yref="y4",
        axref="x4", ayref="y4",
        showarrow=True,
        arrowhead=2,
        arrowsize=1,
        arrowwidth=2,
        arrowcolor="red",
        row=2, col=2
    )
    
    # 中心差分箭头（从i-1和i+1到i）
    fig.add_annotation(
        x=1, y=-0.05,
        ax=2, ay=-0.05,
        xref="x4", yref="y4",
        axref="x4", ayref="y4",
        showarrow=True,
        arrowhead=2,
        arrowsize=1,
        arrowwidth=2,
        arrowcolor="blue",
        row=2, col=2
    )
    
    fig.add_annotation(
        x=3, y=-0.05,
        ax=2, ay=-0.05,
        xref="x4", yref="y4",
        axref="x4", ayref="y4",
        showarrow=True,
        arrowhead=2,
        arrowsize=1,
        arrowwidth=2,
        arrowcolor="blue",
        row=2, col=2
    )
    
    # QUICK格式箭头（从i-2, i-1, i+1到i）
    fig.add_annotation(
        x=0, y=0.05,
        ax=2, ay=0.05,
        xref="x4", yref="y4",
        axref="x4", ayref="y4",
        showarrow=True,
        arrowhead=2,
        arrowsize=1,
        arrowwidth=2,
        arrowcolor="green",
        row=2, col=2
    )
    
    fig.add_annotation(
        x=1, y=0.05,
        ax=2, ay=0.05,
        xref="x4", yref="y4",
        axref="x4", ayref="y4",
        showarrow=True,
        arrowhead=2,
        arrowsize=1,
        arrowwidth=2,
        arrowcolor="green",
        row=2, col=2
    )
    
    fig.add_annotation(
        x=3, y=0.05,
        ax=2, ay=0.05,
        xref="x4", yref="y4",
        axref="x4", ayref="y4",
        showarrow=True,
        arrowhead=2,
        arrowsize=1,
        arrowwidth=2,
        arrowcolor="green",
        row=2, col=2
    )
    
    # 添加文本标签
    fig.add_annotation(
        x=1.5, y=-0.15,
        text="Upwind",
        showarrow=False,
        font=dict(color="red", size=12),
        row=2, col=2
    )
    
    fig.add_annotation(
        x=2, y=-0.25,
        text="Central",
        showarrow=False,
        font=dict(color="blue", size=12),
        row=2, col=2
    )
    
    fig.add_annotation(
        x=2, y=0.2,
        text="QUICK",
        showarrow=False,
        font=dict(color="green", size=12),
        row=2, col=2
    )
    
    # 更新布局
    fig.update_layout(
        title_text="Interactive CFD Visualization: Comparison of Convection Schemes",
        title_font_size=20,
        title_x=0.5,
        height=900,
        width=1200,
        showlegend=True,
        legend=dict(
            yanchor="top",
            y=0.99,
            xanchor="left",
            x=1.02
        )
    )
    
    # 更新轴标签
    fig.update_xaxes(title_text="Position x", row=1, col=1)
    fig.update_yaxes(title_text="Velocity u", row=1, col=1)
    
    fig.update_xaxes(title_text="Position x", row=1, col=2)
    fig.update_yaxes(title_text="Velocity u", row=1, col=2)
    
    fig.update_xaxes(title_text="Position x", row=2, col=1)
    fig.update_yaxes(title_text="Velocity u", row=2, col=1)
    
    fig.update_xaxes(title_text="Grid Point Index", row=2, col=2, range=[-0.5, 4.5])
    fig.update_yaxes(title_text="", row=2, col=2, range=[-0.3, 0.3], showticklabels=False)
    
    # 保存为HTML文件
    html_filename = "05_interactive_cfd.html"
    pio.write_html(fig, html_filename, auto_open=True)
    
    print(f"✓ 交互式图表已保存为 '{html_filename}'")
    print("正在浏览器中打开...")
    
    return fig

def create_interactive_grid_visualization():
    """创建交互式网格可视化"""
    
    print("\n创建交互式网格可视化...")
    
    # 创建网格数据
    n_cells = 5
    dx = 1.0
    x_vertices = np.linspace(0, n_cells*dx, n_cells + 1)
    x_centers = (x_vertices[:-1] + x_vertices[1:]) / 2
    
    # 创建子图
    fig = make_subplots(
        rows=1, cols=2,
        subplot_titles=('Vertex-centered Storage', 'Cell-centered Storage'),
        horizontal_spacing=0.15
    )
    
    # 1. 顶点中心存储
    # 添加顶点
    fig.add_trace(
        go.Scatter(
            x=x_vertices,
            y=[0]*len(x_vertices),
            mode='markers+text',
            name='Vertices',
            marker=dict(size=15, color='red', symbol='circle'),
            text=[f'u{i}' for i in range(len(x_vertices))],
            textposition="top center",
            showlegend=False
        ),
        row=1, col=1
    )
    
    # 添加控制体矩形（使用形状）
    for i in range(n_cells):
        fig.add_shape(
            type="rect",
            x0=x_vertices[i],
            y0=-0.1,
            x1=x_vertices[i+1],
            y1=0.1,
            line=dict(color="blue", width=2),
            fillcolor="rgba(0, 0, 255, 0.1)",
            row=1, col=1
        )
        
        # 添加控制体标签
        fig.add_annotation(
            x=(x_vertices[i] + x_vertices[i+1])/2,
            y=0.15,
            text=f"Cell {i}",
            showarrow=False,
            font=dict(size=10),
            row=1, col=1
        )
    
    # 2. 单元中心存储
    # 添加单元中心
    fig.add_trace(
        go.Scatter(
            x=x_centers,
            y=[0]*len(x_centers),
            mode='markers+text',
            name='Cell Centers',
            marker=dict(size=15, color='green', symbol='circle'),
            text=[f'u{i}' for i in range(len(x_centers))],
            textposition="top center",
            showlegend=False
        ),
        row=1, col=2
    )
    
    # 添加控制体矩形
    for i in range(n_cells):
        fig.add_shape(
            type="rect",
            x0=x_vertices[i],
            y0=-0.1,
            x1=x_vertices[i+1],
            y1=0.1,
            line=dict(color="orange", width=2),
            fillcolor="rgba(255, 165, 0, 0.1)",
            row=1, col=2
        )
        
        # 添加控制体标签
        fig.add_annotation(
            x=(x_vertices[i] + x_vertices[i+1])/2,
            y=0.15,
            text=f"Cell {i}",
            showarrow=False,
            font=dict(size=10),
            row=1, col=2
        )
    
    # 标记边界
    fig.add_annotation(
        x=x_vertices[0],
        y=0.25,
        text="Boundary",
        showarrow=True,
        arrowhead=2,
        font=dict(color="red", size=12),
        row=1, col=1
    )
    
    fig.add_annotation(
        x=x_vertices[-1],
        y=0.25,
        text="Boundary",
        showarrow=True,
        arrowhead=2,
        font=dict(color="red", size=12),
        row=1, col=1
    )
    
    # 更新布局
    fig.update_layout(
        title_text="Interactive CFD Grid Visualization",
        title_font_size=20,
        title_x=0.5,
        height=500,
        width=1000,
        showlegend=False
    )
    
    # 更新轴
    for col in [1, 2]:
        fig.update_xaxes(title_text="Position x", row=1, col=col, range=[-0.5, n_cells*dx+0.5])
        fig.update_yaxes(title_text="", row=1, col=col, range=[-0.3, 0.3], showticklabels=False)
    
    # 保存为HTML
    html_filename = "05_interactive_grid.html"
    pio.write_html(fig, html_filename, auto_open=True)
    
    print(f"✓ 网格可视化已保存为 '{html_filename}'")
    
    return fig

def create_interactive_boundary_conditions():
    """创建交互式边界条件可视化"""
    
    print("\n创建交互式边界条件可视化...")
    
    n_cells = 5
    dx = 1.0
    
    # 创建图形
    fig = go.Figure()
    
    # 真实计算点
    x_real = np.linspace(0, n_cells*dx, n_cells + 1)
    
    # 虚拟点
    x_ghost_left = [-dx]
    x_ghost_right = [n_cells*dx + dx]
    
    # 添加点
    fig.add_trace(go.Scatter(
        x=x_real,
        y=[0]*len(x_real),
        mode='markers+text',
        name='Real Points',
        marker=dict(size=15, color='blue', symbol='circle'),
        text=[f'u{i}' for i in range(len(x_real))],
        textposition="top center"
    ))
    
    fig.add_trace(go.Scatter(
        x=x_ghost_left + x_ghost_right,
        y=[0, 0],
        mode='markers+text',
        name='Ghost Cells',
        marker=dict(size=15, color='red', symbol='square'),
        text=['Ghost L', 'Ghost R'],
        textposition="top center"
    ))
    
    # 添加区域
    fig.add_vrect(
        x0=-dx, x1=0,
        fillcolor="red", opacity=0.1,
        line_width=0,
        annotation_text="Ghost Cell Region",
        annotation_position="top left"
    )
    
    fig.add_vrect(
        x0=n_cells*dx, x1=n_cells*dx+dx,
        fillcolor="red", opacity=0.1,
        line_width=0,
        annotation_text="Ghost Cell Region"
    )
    
    fig.add_vrect(
        x0=0, x1=n_cells*dx,
        fillcolor="green", opacity=0.1,
        line_width=0,
        annotation_text="Computational Domain",
        annotation_position="top"
    )
    
    # 添加边界线
    fig.add_vline(x=0, line_width=3, line_dash="dash", line_color="red")
    fig.add_vline(x=n_cells*dx, line_width=3, line_dash="dash", line_color="red")
    
    # 添加边界条件说明
    fig.add_annotation(
        x=-dx/2, y=0.2,
        text="Dirichlet BC:<br>u = u_boundary",
        showarrow=False,
        font=dict(size=11),
        align="center"
    )
    
    fig.add_annotation(
        x=n_cells*dx + dx/2, y=0.2,
        text="Neumann BC:<br>∂u/∂x = 0",
        showarrow=False,
        font=dict(size=11),
        align="center"
    )
    
    # 更新布局
    fig.update_layout(
        title="Boundary Conditions and Ghost Cells",
        xaxis_title="Position x",
        yaxis_title="",
        height=500,
        width=800,
        showlegend=True,
        yaxis=dict(showticklabels=False, range=[-0.3, 0.4])
    )
    
    # 保存为HTML
    html_filename = "05_interactive_boundary.html"
    pio.write_html(fig, html_filename, auto_open=True)
    
    print(f"✓ 边界条件可视化已保存为 '{html_filename}'")
    
    return fig

def main():
    """主函数"""
    print("=" * 60)
    print("Interactive CFD Visualization with Plotly")
    print("=" * 60)
    print("\n选项:")
    print("1. 对流格式比较")
    print("2. 网格存储方式")
    print("3. 边界条件处理")
    print("4. 全部创建")
    
    try:
        choice = input("\n请选择 (1-4, 默认=1): ").strip()
        if choice == "":
            choice = "1"
        choice = int(choice)
    except:
        choice = 1
    
    if choice == 1:
        fig = create_interactive_cfd_plot()
    elif choice == 2:
        fig = create_interactive_grid_visualization()
    elif choice == 3:
        fig = create_interactive_boundary_conditions()
    elif choice == 4:
        print("\n创建所有交互式可视化...")
        fig1 = create_interactive_cfd_plot()
        fig2 = create_interactive_grid_visualization()
        fig3 = create_interactive_boundary_conditions()
        print("\n✓ 所有可视化已创建完成！")
    else:
        fig = create_interactive_cfd_plot()
    
    print("\n" + "=" * 60)
    print("说明:")
    print("- 交互式图表已保存为HTML文件")
    print("- 可以在任何浏览器中打开查看")
    print("- 支持缩放、平移、悬停查看数据点")
    print("=" * 60)

if __name__ == "__main__":
    # 检查是否安装了plotly
    try:
        import plotly
        main()
    except ImportError:
        print("错误: 需要安装plotly库")
        print("请运行: pip install plotly")
        print("或: pip install plotly numpy")