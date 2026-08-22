"""CR3BP 场景、轨迹、事件、能量和局部放大图的 Makie 可视化工具。"""
module Visualization

# Export functionality
export plot_traj_rot!, plot_traj_inertial!, plot_traj_ε2!
export cr3bp_set_scn!


using ..Dynamics

include("magInsert.jl")

# External dependencies
using DifferentialEquations, LinearAlgebra, CairoMakie

const ps = Bcr4bp_Aux()
const μ = ps.EMRot.μ
const x_L₁  = solve_L1_L2_x(μ , which = :L1)
const x_L₂  = solve_L1_L2_x(μ , which = :L2)
const scatter_args = (;markersize = 10  , strokewidth = 0.5, strokecolor = :black)

const Marker = Dict(
    :collision => :star8,
    :p1collision => :star8,
    :p2collision => :star8,
    :enter => :utriangle,
    :escape => :hexagon,
    :energy_zero => :rect,
    :origin => :circle,
    :intersect => :star4,
    :perigee => :dtriangle,
    :apogee => :utriangle,
    :perilune => :dtriangle,
    :apolune => :utriangle,
    :y_zero_cross => :cross,
)

@inline _event_marker(code::Symbol) = get(Marker, code, :circle)

"""
    cr3bp_set_scn!(ax::Axis; lims=(0.5, 1.5, -0.5, 0.5)) -> Axis
    cr3bp_set_scn!(ax::Axis3; lims=(0.5, 1.5, -0.5, 0.5, -0.5, 0.5)) -> Axis3

在 Makie 二维或三维轴上绘制地球、月球以及 L1–L5 标记，并设置显示范围。

# Arguments
- `ax`: 要原位修改的 `Axis` 或 `Axis3`。

# Keywords
- `lims`: 二维轴使用 `(xmin, xmax, ymin, ymax)`；三维轴使用
  `(xmin, xmax, ymin, ymax, zmin, zmax)`。

# Returns
返回输入轴 `ax`，便于继续链式配置。维数、数值类型或上下限顺序不合法时抛出
`ArgumentError`。
"""
function cr3bp_set_scn!(ax::Axis3; lims = (0.5, 1.5, -0.5, 0.5, -0.5, 0.5))
    _validate_scene_lims(lims, 6)

    mass_parameter = ps.EMRot.μ
    x_l1 = solve_L1_L2_x(mass_parameter; which = :L1)
    x_l2 = solve_L1_L2_x(mass_parameter; which = :L2)
    x_l4_l5 = 0.5 - mass_parameter
    y_l4_l5 = sqrt(3) / 2
    earth_radius = 6378.0 / ps.dim.EMRot_l
    moon_radius = 1737.4 / ps.dim.EMRot_l

    scatter!(
        ax,
        [x_l4_l5, x_l4_l5, x_l1, x_l2],
        [y_l4_l5, -y_l4_l5, 0, 0],
        zeros(4);
        marker = '+',
        markersize = 10,
        color = :red,
    )
    mesh!(ax, Sphere(Point3f(-mass_parameter, 0, 0), earth_radius); color = :dodgerblue)
    mesh!(ax, Sphere(Point3f(1 - mass_parameter, 0, 0), moon_radius); color = :grey70)

    xlims!(ax, lims[1], lims[2])
    ylims!(ax, lims[3], lims[4])
    zlims!(ax, lims[5], lims[6])
    return ax
end

"""验证二维/三维场景 limits 的长度、数值类型与上下界顺序。"""
function _validate_scene_lims(lims, expected_length)
    length(lims) == expected_length || throw(ArgumentError(
        "lims for this axis must contain $expected_length values; got $(length(lims)).",
    ))
    all(isreal, lims) || throw(ArgumentError("lims must contain only real values."))

    for index in 1:2:expected_length
        lims[index] < lims[index + 1] || throw(ArgumentError(
            "Each lower limit in lims must be smaller than its upper limit.",
        ))
    end
    return nothing
end

function cr3bp_set_scn!(ax::Axis; lims = (0.5, 1.5 , -0.5 , 0.5))
    _validate_scene_lims(lims, 4)
    x_L₄ , y_L₄ =  0.5-μ ,  √3 /2
    x_L₅ , y_L₅ =  0.5-μ , -√3 /2
    scatter!(ax , [x_L₄, x_L₅] , [y_L₄, y_L₅] ; marker = '+' , markersize = 10 , color =:red )

    
    R_earth = 6378.0  / ps.dim.EMRot_l   
    R_moon  = 1737.4  / ps.dim.EMRot_l   

    # 地球（蓝色），位于 (-μ, 0)
    poly!(ax, Circle(Point2f(-μ, 0), R_earth); color = :dodgerblue, strokecolor = :black, strokewidth = 1)

    # 月球（灰色），位于 (1-μ, 0)
    poly!(ax, Circle(Point2f(1 - μ, 0), R_moon);  color = :grey70,      strokecolor = :black, strokewidth = 1)

    scatter!(ax , [x_L₁, x_L₂] , [0, 0] ; marker = '+' , markersize = 10 , color =:red )
    # scatter!(ax , [-μ, 1-μ] , [0, 0] ; marker = Marker[:origin] ,color = :grey80,  scatter_args... )

    xlims!(ax,  lims[1] , lims[2])
    ylims!(ax , lims[3] , lims[4])
    return ax
end

# 旋转系

"""
    plot_traj_rot!(ax, sol, events; kwargs...) -> Lines

在旋转坐标系中绘制 ODE 轨迹、起点、平动点以及已记录事件。

# Arguments
- `ax`: Makie 二维轴。
- `sol`: 支持稠密插值的 ODE 解。
- `events`: `Event` 集合；未知事件代码使用圆形标记。

# Keywords
额外关键字原样传给 `lines!`。

# Returns
返回轨迹对应的 Makie `Lines` 图元。
"""
function plot_traj_rot!(ax , sol, events; args...)
    t_range = range(sol.t[1], sol.t[end], length=1500)
    u_cart = sol(t_range).u
    xs = getindex.(u_cart, 1)
    ys = getindex.(u_cart, 2)
    trajectory_plot = lines!(ax, xs, ys; args...)

    scatter!(ax , xs[1] , ys[1] ; marker = Marker[:origin] , scatter_args...)
    scatter!(ax , [x_L₁, x_L₂] , [0, 0] ; marker = '+' , markersize = 6 , color =:red )

    for event in events
        u_cart_event = event.state
        scatter!(ax, u_cart_event[1], u_cart_event[2];
                 marker=_event_marker(event.code), scatter_args...)
    end
    return trajectory_plot
end

"""
    plot_traj_inertial!(ax, sol, events; kwargs...) -> Lines

将平面旋转系轨迹转换到第二主天体中心惯性系后绘制，并在相同坐标系标记事件。

# Returns
返回惯性轨迹的 Makie `Lines` 图元；额外关键字传给 `lines!`。
"""
function plot_traj_inertial!(ax , sol, events; args...)
    t_range = range(sol.t[1], sol.t[end], length=1500)
    u_cart = sol(t_range).u
    u_inertial = cr3bp_rotating_to_inertial.(μ , t_range , u_cart , center = :p2)

    # 惯性系
    xs = getindex.(u_inertial, 1)
    ys = getindex.(u_inertial, 2)

    trajectory_plot = lines!(ax, xs, ys; args...)

    scatter_args = (;markersize = 10 , color = :transparent , strokewidth = 0.5, strokecolor = :black)

    scatter!(ax , xs[1] , ys[1] ; marker = Marker[:origin] , scatter_args...)

    for event in events
        u_cart_event = event.state
        u_cart_event_inertial = cr3bp_rotating_to_inertial(μ , event.time , u_cart_event , center = :p2)
        scatter!(ax, u_cart_event_inertial[1], u_cart_event_inertial[2];
                 marker=_event_marker(event.code), scatter_args...)
    end
    return trajectory_plot
end

"""
    plot_traj_ε2!(ax, sol, events; kwargs...) -> (energy_plot, rate_plot)

绘制相对第二主天体的两体比机械能 `ε₂` 及其时间导数，并在曲线上标记事件。

# Returns
返回 `(energy_plot, rate_plot)` 两个 Makie `Lines` 图元。额外关键字传给两条曲线的
`lines!` 调用，导数曲线额外使用虚线。
"""
function plot_traj_ε2!(ax , sol, events; args...)
    t_range = range(sol.t[1], sol.t[end], length=1500)
    u_cart = sol(t_range).u

    ε2 = compute_ε2.(u_cart ,t_range , μ)
    ε2_dot = compute_ε2_dot.(u_cart , t_range , μ, ps)

    l1 = lines!(ax , t_range , ε2; args...)
    l2 = lines!(ax , t_range , ε2_dot ; linestyle = :dash, args...)
    scatter!(ax, 0, 0; marker=Marker[:energy_zero], scatter_args...)

    for event in events
        marker = _event_marker(event.code)
        scatter!(ax, event.time, compute_ε2(event.state, event.time, μ);
                 marker, scatter_args...)
        scatter!(ax, event.time, compute_ε2_dot(event.state, event.time, μ, ps);
                 marker, scatter_args...)
    end

    axislegend(ax, [l1, l2], ["ε2 ", "ε2_dot"], position = :rt)
    return l1, l2
end

end
