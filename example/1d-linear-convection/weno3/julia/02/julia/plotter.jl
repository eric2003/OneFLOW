# julia/plotter.jl
"""
CFDPlotter 的 Julia 实现（通过 PythonCall.jl 调用 Matplotlib）
确保与 Python plotter.py 行为完全一致
"""

using PythonCall

# 初始化 Python 环境（加载 matplotlib, inflect）
const plt = pyimport("matplotlib.pyplot")
const inflect = pyimport("inflect")

mutable struct CFDPlotter
    default_styles::Dict{String, Any}
    p::Py
end

function CFDPlotter()
	default_styles = Dict{String, Any}(
		"numerical" => Dict(
			:color => "blue",
			:linestyle => "-",
			:marker => "o",
			:markerfacecolor => "none"
		),
		"analytical" => Dict(
			:color => "red",
			:linestyle => "--",
			:marker => "",
			:linewidth => 1.5
		),
		"comparison" => [
			Dict(:color => "black", :linestyle => "-", :marker => "o", :markerfacecolor => "none"),
			Dict(:color => "blue", :linestyle => "--", :marker => "s", :markerfacecolor => "none"),
			Dict(:color => "green", :linestyle => ":", :marker => "^", :markerfacecolor => "none")
		]
	)	
    p = inflect.engine()
    CFDPlotter(default_styles, p)
end

"""
轻量即时绘图（快速验证结果）
"""
function plot_quick(plotter::CFDPlotter, cfd_result::Dict; title=nothing, show=true, save_path=nothing)
    plt.figure("OneFLOW-CFD Solver", figsize=(10, 6))

    # 自动生成标题
    actual_title = title
    if actual_title === nothing
        rk_order = cfd_result["config"]["rk_order"]
        rk_str = plotter.p.ordinal(rk_order)
        final_time = cfd_result["config"]["final_time"]
        order = cfd_result["config"]["order"]
        scheme = uppercase(cfd_result["config"]["scheme"])
        actual_title = "1D Convection (t=$(final_time))\n$(order)th-order $(scheme) + $(rk_str)-order RK"
    end

    # 绘制数值解
    plt.plot(
        cfd_result["x"], cfd_result["numerical"];
        label="Numerical ($(uppercase(cfd_result["config"]["scheme"])))",
        plotter.default_styles["numerical"]...,
        markersize=5, linewidth=0.5
    )

    # 绘制解析解
    plt.plot(
        cfd_result["x"], cfd_result["analytical"];
        label="Analytical",
        plotter.default_styles["analytical"]...
    )

    # 通用样式
    _set_common_style(plotter, actual_title)

    if save_path !== nothing
        plt.savefig(save_path, dpi=300, bbox_inches="tight")
    end
    if show
        plt.show()
    end
    plt.close()
end

"""
多格式/多精度对比绘图
"""
function plot_comparison(plotter::CFDPlotter, result_list::Vector{Dict{String, Any}}; title=nothing, show=true, save_path=nothing)
    plt.figure("OneFLOW-CFD Solver", figsize=(10, 6))

    # 自动生成标题
    actual_title = title
    if actual_title === nothing
        schemes = [uppercase(r["config"]["scheme"]) * string(r["config"]["order"]) for r in result_list]
        rk_order = result_list[1]["config"]["rk_order"]
        rk_str = plotter.p.ordinal(rk_order)
        final_time = result_list[1]["config"]["final_time"]
        actual_title = "1D Convection Comparison (t=$(final_time))\n$(join(schemes, ", ")) + $(rk_str)-order RK"
    end

    # 绘制多个数值解
    for (i, res) in enumerate(result_list)
        style = plotter.default_styles["comparison"][mod1(i, length(plotter.default_styles["comparison"]))]
        label = "Numerical ($(uppercase(res["config"]["scheme"]))$(res["config"]["order"]))"
        plt.plot(
            res["x"], res["numerical"];
            label=label,
            style...,
            markersize=5, linewidth=0.5
        )
    end

    # 绘制解析解
    plt.plot(
        result_list[1]["x"], result_list[1]["analytical"];
        label="Analytical",
        plotter.default_styles["analytical"]...
    )

    _set_common_style(plotter, actual_title)

    if save_path !== nothing
        plt.savefig(save_path, dpi=300, bbox_inches="tight")
    end
    if show
        plt.show()
    end
    plt.close()
end

function _set_common_style(plotter::CFDPlotter, title::String)
    plt.title(title, fontsize=12)
    plt.xlabel("x", fontsize=10)
    plt.ylabel("u", fontsize=10)
    plt.legend(fontsize=9)
    plt.grid(true, color="gray", linestyle="--", linewidth=0.5, alpha=0.7)
    plt.tight_layout()
end

"""
快捷函数：ENO/WENO对比绘图
"""
function plot_eno_weno_comparison(eno_result::Dict, weno_result::Dict; save_path=nothing)
    plotter = CFDPlotter()
    plot_comparison(plotter, [eno_result, weno_result]; save_path=save_path)
end