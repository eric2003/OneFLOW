"""
文件名: 05_interactive_cfd_plot.py
功能: 创建交互式CFD可视化
包含: 使用Plotly创建可交互的CFD图像
注意: 需要安装plotly库: pip install plotly
"""

import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots

def create_interactive_plot():
    """创建交互式CFD可视化"""
    
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
        subplot_titles=('Upwind Scheme', 
                       'Central Difference Scheme',
                       'QUICK Scheme',
                       'Grid Point Stencil'),
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
            showlegend=True
        ),
        row=1, col=1
    )
    
    fig.add_trace(
        go.Scatter(
            x=x, y=u_upwind,
            mode='lines+markers',
            name='Upwind',
            line=dict(color='red', width=2, dash='dash'),
            marker=dict(size=4, symbol='circle'),
            showlegend=True
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
            showlegend=False
        ),
        row=1, col=2
    )
    
    fig.add_trace(
        go.Scatter(
            x=x, y=u_central,
            mode='lines+markers',
            name='Central',
            line=dict(color='blue', width=2, dash='dash'),
            marker=dict(size=4, symbol='square'),
            showlegend=True
        ),
        row=1