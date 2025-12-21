import matplotlib.pyplot as plt
import numpy as np
import inflect

class CFDPlotter:
    """CFD可视化工具类：解耦绘图逻辑"""
    def __init__(self):
        # 预设样式（统一管理）
        self.default_styles = {
            "numerical": {"color": "blue", "linestyle": "-", "marker": "o", "markerfacecolor": "none"},
            "analytical": {"color": "red", "linestyle": "--", "marker": "", "linewidth": 1.5},
            "comparison": [
                {"color": "black", "linestyle": "-", "marker": "o", "markerfacecolor": "none"},
                {"color": "blue", "linestyle": "--", "marker": "s", "markerfacecolor": "none"},
                {"color": "green", "linestyle": ":", "marker": "^", "markerfacecolor": "none"},
            ]
        }
        self.p = inflect.engine()

    def plot_quick(self, cfd_result, title=None, show=True, save_path=None):
        """轻量即时绘图（快速验证结果）"""
        plt.figure(figsize=(10, 6))

        # 自动生成标题
        if title is None:
            rk_str = self.p.ordinal(cfd_result["config"]["rk_order"])
            title = (f'1D Convection (t={cfd_result["config"]["final_time"]:.3f})\n'
                     f'{cfd_result["config"]["order"]}th-order {cfd_result["config"]["scheme"].upper()} + {rk_str}-order RK')

        # 绘制数值解
        plt.plot(
            cfd_result["x"], cfd_result["numerical"],
            label=f'Numerical ({cfd_result["config"]["scheme"].upper()})',
            **self.default_styles["numerical"],
            markersize=5, linewidth=0.5
        )

        # 绘制解析解
        plt.plot(
            cfd_result["x"], cfd_result["analytical"],
            label='Analytical',
            **self.default_styles["analytical"]
        )

        # 通用样式
        self._set_common_style(title)

        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
        if show:
            plt.show()
        plt.close()

    def plot_comparison(self, result_list, title=None, show=True, save_path=None):
        """多格式/多精度对比绘图"""
        plt.figure(figsize=(10, 6))

        # 自动生成标题
        if title is None:
            schemes = [f'{r["config"]["scheme"].upper()}{r["config"]["order"]}' for r in result_list]
            rk_str = self.p.ordinal(result_list[0]["config"]["rk_order"])
            title = (f'1D Convection Comparison (t={result_list[0]["config"]["final_time"]:.3f})\n'
                     f'{", ".join(schemes)} + {rk_str}-order RK')

        # 绘制多个数值解
        for i, res in enumerate(result_list):
            style = self.default_styles["comparison"][i % len(self.default_styles["comparison"])]
            label = f'Numerical ({res["config"]["scheme"].upper()}{res["config"]["order"]})'
            plt.plot(
                res["x"], res["numerical"],
                label=label,
                **style,
                markersize=5, linewidth=0.5
            )

        # 绘制解析解
        plt.plot(
            result_list[0]["x"], result_list[0]["analytical"],
            label='Analytical',
            **self.default_styles["analytical"]
        )

        # 通用样式
        self._set_common_style(title)

        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
        if show:
            plt.show()
        plt.close()

    def _set_common_style(self, title):
        """统一设置图表样式"""
        plt.title(title, fontsize=12)
        plt.xlabel('x', fontsize=10)
        plt.ylabel('u', fontsize=10)
        plt.legend(fontsize=9)
        plt.grid(True, color='gray', linestyle='--', linewidth=0.5, alpha=0.7)
        plt.tight_layout()

# 快捷函数：ENO/WENO对比绘图
def plot_eno_weno_comparison(eno_result, weno_result, save_path=None):
    plotter = CFDPlotter()
    plotter.plot_comparison(
        result_list=[eno_result, weno_result],
        save_path=save_path
    )