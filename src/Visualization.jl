"""
Visualization模块

提供轨迹可视化功能，包括旋转系、惯性系和能量参数的可视化。
包含plot_traj_rot!、plot_traj_inertial!和plot_traj_ε2!等绘图函数。
"""
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

const Marker = Dict(:collision => :star8,
            :escape => :hexagon,
            :energy_zero => :rect,
            :origin => :circle,
            :intersect => :star4)

"""
    cr3bp_set_scn!(ax; lims=(0.5, 1.5, -0.5, 0.5))

在二维 Makie 轴上绘制地球、月球及 L1–L5 标记，并设置显示范围。
"""
function cr3bp_set_scn!(ax; lims = (0.5, 1.5 , -0.5 , 0.5))
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

end

# 旋转系

"画旋转系轨迹"
# 在旋转坐标系中绘制轨迹、起点和已记录事件。
function plot_traj_rot!(ax , sol, events; args...)
    t_range = range(sol.t[1], sol.t[end], length=1500)
    u_cart = sol(t_range).u
    xs = [s[1]  for s in u_cart]
    ys = [s[2]  for s in u_cart]
    lines!(ax , xs , ys ; args...)

    scatter!(ax , xs[1] , ys[1] ; marker = Marker[:origin] , scatter_args...)
    scatter!(ax , [x_L₁, x_L₂] , [0, 0] ; marker = '+' , markersize = 6 , color =:red )

    for event in events
        u_cart_event = event.state
        scatter!(ax , u_cart_event[1] ,u_cart_event[2]  ; marker = Marker[event.code] ,scatter_args... )
    end
end

"画p2惯性系轨迹"
# 将轨迹转换到以第二主天体为中心的惯性系后绘制。
function plot_traj_inertial!(ax , sol, events; args...)
    t_range = range(sol.t[1], sol.t[end], length=1500)
    u_cart = sol(t_range).u
    u_inertial = cr3bp_rotating_to_inertial.(μ , t_range , u_cart , center = :p2)

    # 惯性系
    xs = [s[1]  for s in u_inertial]
    ys = [s[2]  for s in u_inertial]

    lines!(ax , xs , ys; args...)

    scatter_args = (;markersize = 10 , color = :transparent , strokewidth = 0.5, strokecolor = :black)

    scatter!(ax , xs[1] , ys[1] ; marker = Marker[:origin] , scatter_args...)

    for event in events
        u_cart_event = event.state
        u_cart_event_inertial = cr3bp_rotating_to_inertial(μ , event.time , u_cart_event , center = :p2)
        scatter!(ax , u_cart_event_inertial[1] ,u_cart_event_inertial[2]  ; marker = Marker[event.code] ,scatter_args... )
    end
end

"画ε2及其导数"
# 绘制相对第二主天体的能量 `ε₂` 及其时间导数，并标记事件。
function plot_traj_ε2!(ax , sol, events; args...)
    t_range = range(sol.t[1], sol.t[end], length=1500)
    u_cart = sol(t_range).u

    ε2 = compute_ε2.(u_cart ,t_range , μ)
    ε2_dot = compute_ε2_dot.(u_cart , t_range , μ, ps)

    l1 = lines!(ax , t_range , ε2; args...)
    l2 = lines!(ax , t_range , ε2_dot ; linestyle = :dash, args...)
    scatter!(ax , 0 ,0  ; marker = Marker[:energy_zero],scatter_args... )

    for event in events
        scatter!(ax , event.time ,compute_ε2(event.state ,event.time, μ)  ; marker = Marker[event.code],scatter_args... )
        scatter!(ax , event.time ,compute_ε2_dot(event.state ,event.time , μ , ps)  ; marker = Marker[event.code] ,scatter_args... )
    end

    axislegend(ax, [l1, l2], ["ε2 ", "ε2_dot"], position = :rt)
end

end
