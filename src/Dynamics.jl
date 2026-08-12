"""
Dynamics module - Core dynamics equations for the circular restricted three-body problem.

This module provides equations of motion, state transition matrix propagation,
and energy computations for CR3BP and BCR4BP trajectory analysis.
"""
module Dynamics

# Export functionality

using DifferentialEquations
using LinearAlgebra
using StaticArrays
using ForwardDiff
using CairoMakie

export cr3bp_eqm!
export bcr4bp_eqm!, bcr4bp_eqmEMRot2D!

export ode_params , integration

include("dynamic_parameters.jl") 

"""
    ode_params

封装一次数值积分所需的动力学函数、辅助参数、求解器关键字和可选绘图配置。
"""
struct ode_params 
    dynamics    :: Function
    dynamicsargs :: Union{Nothing , NamedTuple}
    aux         :: NamedTuple
end

"""使用 CR3BP 默认辅助参数和 `1e-12` 误差容限构造 `ode_params`。"""
ode_params(dynamics) = ode_params(dynamics ,(;abstol = 1e-12 , reltol = 1e-12) , Bcr4bp_Aux())

"""
    integration(x0, tspan, parameters; cb=nothing, odeargs=parameters.dynamicsargs, interp_num=200)

积分 `parameters.dynamics`，返回终端状态、等间隔采样时间和状态。若配置了绘图轴，
会同时在该轴绘制二维轨迹；`cb` 可传入 DifferentialEquations 回调。
"""

@inline function integration(x₀ ,tspan, p::ode_params ; cb = nothing , odeargs = p.dynamicsargs, interp_num = 200 )

    prob = ODEProblem(p.dynamics , x₀ ,   tspan  , p.aux )
    sol = solve(prob , Vern7() ,callback = cb ; odeargs...)
    if sol.retcode ≠ ReturnCode.Success
        return nothing, nothing , nothing
    end

    tt = range(sol.t[1] , sol.t[end] , interp_num)
    uu = sol(tt).u

    return sol.u[end] , tt , uu
end


"""
    cr3bp_eqm!(du, u, p, t)

CR3BP 运动方程的维度分发入口，支持 4 维平面状态、6 维空间状态及其 STM 扩展状态。
"""
function cr3bp_eqm! end

"""根据状态长度调用对应的二维、三维或状态转移矩阵 CR3BP 方程。"""
function cr3bp_eqm!(du::AbstractArray, u::AbstractArray, p, t)
    n = length(u)
    if n == 4
        return cr3bp_eqm2D!(du, u, p, t)
    elseif n == 6
        return cr3bp_eqm3D!(du, u, p, t)
    elseif n == 20 
        return cr3bp_stm2D!(du, u, p, t)
    elseif n == 42 
        return cr3bp_stm3D!(du, u, p, t)
    else
        error("不支持的状态维度: $n")
    end
end

"""计算旋转坐标系中四维平面 CR3BP 状态 `[x, y, vx, vy]` 的导数。"""
function cr3bp_eqm2D!(du, u, p, t)  
    μ = p.EMRot.μ

    x , y , vx , vy  = u 

    r1_3pow = ((x + μ)^2 + y^2)^(1.5)  
    r2_3pow = ((x + μ - 1)^2 + y^2 )^(1.5)  
    
    du[1] =  vx  
    du[2] =  vy  
    du[3] =  x + 2*vy - (1 - μ) * (x + μ) / r1_3pow - μ*(x + μ - 1) / r2_3pow  
    du[4] =  y - 2*vx - (1 - μ) * y / r1_3pow - μ * y / r2_3pow  
    
    return nothing  
end 

"""计算旋转坐标系中六维空间 CR3BP 状态的导数。"""
function cr3bp_eqm3D!(du, u, p, t)  
    μ = p.EMRot.μ

    x , y , z , vx , vy , vz  = u 

    r1_3pow = ((x + μ)^2 + y^2 + z^2)^(1.5)  
    r2_3pow = ((x + μ - 1)^2 + y^2 + z^2)^(1.5)  
    
    du[1] =  vx  
    du[2] =  vy  
    du[3] =  vz  
    du[4] =  x + 2*vy - (1 - μ) * (x + μ) / r1_3pow - μ*(x + μ - 1) / r2_3pow  
    du[5] =  y - 2*vx - (1 - μ) * y / r1_3pow - μ * y / r2_3pow  
    du[6] =  -(1 - μ) * z / r1_3pow - μ * z / r2_3pow  
    
    return nothing  
end  

"""
    cr3bp_stm2D!(du, u, p, t)

计算平面 CR3BP 状态及 4×4 状态转移矩阵（STM）的联合导数；输入状态长度为 20。
"""
function cr3bp_stm2D!(du, u, p, t)
    # 提取质量参数
    μ = p.EMRot.μ
    
    # 提取状态变量
    x, y, vx, vy = u[1], u[2], u[3], u[4]
    
    # 计算到两个主天体的距离
    r1 = hypot(x + μ, y)
    r2 = hypot(x + μ - 1.0, y)
    
    # 计算距离的幂次（避免重复计算）
    r1³ = r1^3
    r2³ = r2^3
    r1⁵ = r1^5
    r2⁵ = r2^5
    
    # 预计算常用项
    dx1 = x + μ
    dx2 = x + μ - 1.0
    μ_comp = 1.0 - μ  # μ的补数
    
    term1 = 3.0 * μ_comp / r1⁵
    term2 = 3.0 * μ / r2⁵
    
    # 计算雅可比矩阵的非零元素
    G11 = 1.0 - μ_comp / r1³ - μ / r2³ + term1 * dx1^2 + term2 * dx2^2
    G12 = term1 * dx1 * y + term2 * dx2 * y
    G22 = 1.0 - μ_comp / r1³ - μ / r2³ + term1 * y^2 + term2 * y^2
    
    # 构建雅可比矩阵
    Df = @SMatrix [
        0.0  0.0  1.0   0.0;
        0.0  0.0  0.0   1.0;
        G11  G12  0.0   2.0;
        G12  G22 -2.0   0.0
    ]
    
    # 提取状态转移矩阵并计算其导数
    Φ = reshape(view(u, 5:20), 4, 4)
    dΦ = Df * Φ
    
    # 计算状态导数
    du[1] = vx
    du[2] = vy
    du[3] = x + 2.0 * vy - μ_comp * dx1 / r1³ - μ * dx2 / r2³
    du[4] = y - 2.0 * vx - μ_comp * y / r1³ - μ * y / r2³
    
    # 存储STM导数
    
    du[5:20] .=  reduce(vcat ,dΦ )
    
    return nothing
end

"""
    cr3bp_stm3D!(du, u, p, t)

计算空间 CR3BP 状态及 6×6 状态转移矩阵（STM）的联合导数；输入状态长度为 42。
"""
function cr3bp_stm3D!(du, u, p, t)
    μ = p.EMRot.μ
    x, y, z, vx, vy, vz = view(u, 1:6)
    dx1, dx2 = x + μ, x + μ - 1.0
    r1sq, r2sq = dx1^2 + y^2 + z^2, dx2^2 + y^2 + z^2
    r1_3, r2_3 = r1sq^(3 / 2), r2sq^(3 / 2)
    r1_5, r2_5 = r1sq^(5 / 2), r2sq^(5 / 2)
    μ1 = 1 - μ

    Uxx = 1 - μ1 / r1_3 - μ / r2_3 + 3μ1 * dx1^2 / r1_5 + 3μ * dx2^2 / r2_5
    Uyy = 1 - μ1 / r1_3 - μ / r2_3 + 3μ1 * y^2 / r1_5 + 3μ * y^2 / r2_5
    Uzz = -μ1 / r1_3 - μ / r2_3 + 3μ1 * z^2 / r1_5 + 3μ * z^2 / r2_5
    Uxy = 3μ1 * dx1 * y / r1_5 + 3μ * dx2 * y / r2_5
    Uxz = 3μ1 * dx1 * z / r1_5 + 3μ * dx2 * z / r2_5
    Uyz = 3μ1 * y * z / r1_5 + 3μ * y * z / r2_5

    Df = [
        0 0 0 1 0 0;
        0 0 0 0 1 0;
        0 0 0 0 0 1;
        Uxx Uxy Uxz 0 2 0;
        Uxy Uyy Uyz -2 0 0;
        Uxz Uyz Uzz 0 0 0;
    ]
    Φ = reshape(view(u, 7:42), 6, 6)
    dΦ = Df * Φ

    du[1:6] .= (vx, vy, vz,
        x + 2vy - μ1 * dx1 / r1_3 - μ * dx2 / r2_3,
        y - 2vx - μ1 * y / r1_3 - μ * y / r2_3,
        -μ1 * z / r1_3 - μ * z / r2_3)
    du[7:42] .= vec(dΦ)
    return nothing
end

"""
    bcr4bp_eqm!(du, u, p, t)

地月旋转系 BCR4BP 的维度分发入口，支持 4 维平面、6 维空间、20 维平面+STM 与
42 维空间+STM 状态。
"""
function bcr4bp_eqm!(du, u, p, t)
    n = length(u)
    if n == 4
        return bcr4bp_eqm2D!(du, u, p, t)
    elseif n == 6
        return bcr4bp_eqm3D!(du, u, p, t)
    elseif n == 20
        return bcr4bp_stm2D!(du, u, p, t)
    elseif n == 42
        return bcr4bp_stm3D!(du, u, p, t)
    end
    throw(ArgumentError("BCR4BP 仅支持长度为 4、6、20 或 42 的状态，得到 $n"))
end

"""
    _bcr4bp_terms(x, y, z, p, t)

返回 BCR4BP 的位置相关加速度项及势函数 Hessian。太阳的间接项与位置无关，
因此只进入 `gx`、`gy`，不进入 Hessian。
"""
function _bcr4bp_terms(x, y, z, p, t)
    μ, μs, a, ωs = p.EMRot.μ, p.EMRot.mus, p.EMRot.as, p.EMRot.ws
    dx1, dy1, dz1 = x + μ, y, z
    dx2, dy2, dz2 = x + μ - 1, y, z
    θs = ωs * t
    xs, ys = a * cos(θs), a * sin(θs)
    dxs, dys, dzs = x - xs, y - ys, z

    r1sq = dx1^2 + dy1^2 + dz1^2
    r2sq = dx2^2 + dy2^2 + dz2^2
    rssq = dxs^2 + dys^2 + dzs^2
    r1_3, r2_3, rs_3 = r1sq^(3 / 2), r2sq^(3 / 2), rssq^(3 / 2)
    r1_5, r2_5, rs_5 = r1sq^(5 / 2), r2sq^(5 / 2), rssq^(5 / 2)
    μ1 = 1 - μ

    gx = x - μ1 * dx1 / r1_3 - μ * dx2 / r2_3 - μs * (dxs / rs_3 + xs / a^3)
    gy = y - μ1 * dy1 / r1_3 - μ * dy2 / r2_3 - μs * (dys / rs_3 + ys / a^3)
    gz = -μ1 * dz1 / r1_3 - μ * dz2 / r2_3 - μs * dzs / rs_3

    Uxx = 1 - μ1 / r1_3 - μ / r2_3 - μs / rs_3 +
          3μ1 * dx1^2 / r1_5 + 3μ * dx2^2 / r2_5 + 3μs * dxs^2 / rs_5
    Uyy = 1 - μ1 / r1_3 - μ / r2_3 - μs / rs_3 +
          3μ1 * dy1^2 / r1_5 + 3μ * dy2^2 / r2_5 + 3μs * dys^2 / rs_5
    Uzz = -μ1 / r1_3 - μ / r2_3 - μs / rs_3 +
          3μ1 * dz1^2 / r1_5 + 3μ * dz2^2 / r2_5 + 3μs * dzs^2 / rs_5
    Uxy = 3μ1 * dx1 * dy1 / r1_5 + 3μ * dx2 * dy2 / r2_5 + 3μs * dxs * dys / rs_5
    Uxz = 3μ1 * dx1 * dz1 / r1_5 + 3μ * dx2 * dz2 / r2_5 + 3μs * dxs * dzs / rs_5
    Uyz = 3μ1 * dy1 * dz1 / r1_5 + 3μ * dy2 * dz2 / r2_5 + 3μs * dys * dzs / rs_5
    return (; gx, gy, gz, Uxx, Uyy, Uzz, Uxy, Uxz, Uyz)
end

"""计算四维平面 BCR4BP 状态 `[x, y, vx, vy]` 的导数。"""
function bcr4bp_eqm2D!(du, u, p, t)
    x, y, vx, vy = view(u, 1:4)
    terms = _bcr4bp_terms(x, y, zero(x), p, t)
    du .= (vx, vy, terms.gx + 2vy, terms.gy - 2vx)
    return nothing
end

"""计算六维空间 BCR4BP 状态 `[x, y, z, vx, vy, vz]` 的导数。"""
function bcr4bp_eqm3D!(du, u, p, t)
    x, y, z, vx, vy, vz = view(u, 1:6)
    terms = _bcr4bp_terms(x, y, z, p, t)
    du .= (vx, vy, vz, terms.gx + 2vy, terms.gy - 2vx, terms.gz)
    return nothing
end

"""计算平面 BCR4BP 状态与 4×4 STM 的联合导数。"""
function bcr4bp_stm2D!(du, u, p, t)
    x, y, vx, vy = view(u, 1:4)
    terms = _bcr4bp_terms(x, y, zero(x), p, t)
    Df = [0 0 1 0; 0 0 0 1; terms.Uxx terms.Uxy 0 2; terms.Uxy terms.Uyy -2 0]
    Φ = reshape(view(u, 5:20), 4, 4)
    du[1:4] .= (vx, vy, terms.gx + 2vy, terms.gy - 2vx)
    du[5:20] .= vec(Df * Φ)
    return nothing
end

"""计算空间 BCR4BP 状态与 6×6 STM 的联合导数。"""
function bcr4bp_stm3D!(du, u, p, t)
    x, y, z, vx, vy, vz = view(u, 1:6)
    terms = _bcr4bp_terms(x, y, z, p, t)
    Df = [
        0 0 0 1 0 0;
        0 0 0 0 1 0;
        0 0 0 0 0 1;
        terms.Uxx terms.Uxy terms.Uxz 0 2 0;
        terms.Uxy terms.Uyy terms.Uyz -2 0 0;
        terms.Uxz terms.Uyz terms.Uzz 0 0 0;
    ]
    Φ = reshape(view(u, 7:42), 6, 6)
    du[1:6] .= (vx, vy, vz, terms.gx + 2vy, terms.gy - 2vx, terms.gz)
    du[7:42] .= vec(Df * Φ)
    return nothing
end


# 事件函数集成
include("dynamic_energy.jl")
include("dynamic_events.jl")
include("rotation.jl")

end
