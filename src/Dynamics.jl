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

export cr3bp_eqm!
export bcr4bp_eqm!

export ode_params , integration

include("dynamic_parameters.jl") 

"""
    ode_params(dynamics, dynamicsargs, aux)

封装一次数值积分所需的动力学函数、求解器关键字和辅助参数。

# Arguments
- `dynamics`: 符合 `f(du, u, p, t)` 约定的原位动力学函数。
- `dynamicsargs`: 传递给 `DifferentialEquations.solve` 的关键字 `NamedTuple`，或
  `nothing`。
- `aux`: 作为 ODE 参数 `p` 传递给动力学函数的辅助参数。

# Fields
- `dynamics`: 动力学函数。
- `dynamicsargs`: 默认求解器关键字。
- `aux`: 模型辅助参数。
"""
struct ode_params 
    dynamics    :: Function
    dynamicsargs :: Union{Nothing , NamedTuple}
    aux         :: NamedTuple
end

"""
    ode_params(dynamics) -> ode_params

使用 `Bcr4bp_Aux()` 和 `abstol = reltol = 1e-12` 构造积分参数。
"""
ode_params(dynamics) = ode_params(dynamics ,(;abstol = 1e-12 , reltol = 1e-12) , Bcr4bp_Aux())

"""
    integration(x0, tspan, parameters; cb=nothing,
                odeargs=parameters.dynamicsargs, interp_num=200)
        -> (final_state, times, states)

使用 `Vern7` 积分 `parameters.dynamics`，并对稠密解进行等时间间隔采样。

# Arguments
- `x0`: 初始状态向量。
- `tspan`: 二元积分区间 `(t0, tf)`。
- `parameters::ode_params`: 动力学函数、模型参数与默认求解器选项。

# Keywords
- `cb=nothing`: DifferentialEquations 回调或 `CallbackSet`。
- `odeargs=parameters.dynamicsargs`: 覆盖默认的 `solve` 关键字。
- `interp_num=200`: 输出采样点数量，必须为正整数。

# Returns
成功时返回 `(final_state, times, states)`：终端状态、采样时间范围以及对应状态
向量。求解器失败时返回 `(nothing, nothing, nothing)`。
"""
function integration(
    x₀,
    tspan,
    p::ode_params;
    cb=nothing,
    odeargs=p.dynamicsargs,
    interp_num::Integer=200,
)
    interp_num > 0 || throw(ArgumentError("interp_num must be positive."))
    solve_kwargs = isnothing(odeargs) ? NamedTuple() : odeargs
    prob = ODEProblem(p.dynamics , x₀ ,   tspan  , p.aux )
    sol = solve(prob, Vern7(); callback=cb, solve_kwargs...)
    if !SciMLBase.successful_retcode(sol)
        return nothing, nothing , nothing
    end

    tt = range(sol.t[1], sol.t[end]; length=interp_num)
    uu = sol(tt).u

    return sol.u[end] , tt , uu
end


"""
    cr3bp_eqm!(du, u, p, t) -> nothing

CR3BP 运动方程的维度分发入口，支持 4 维平面状态、6 维空间状态及其 STM 扩展状态。

# Arguments
- `du`: 写入状态导数的预分配数组。
- `u`: 4、6、20 或 42 元状态；后两种在物理状态后附带列优先 STM。
- `p`: 含 `p.EMRot.μ` 的模型参数，通常由 `Bcr4bp_Aux()` 创建。
- `t`: 当前无量纲时间；CR3BP 自治方程中不显式使用。

# Returns
返回 `nothing`，结果原位写入 `du`。不支持的状态长度会抛出异常。
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

    r1sq = (x + μ)^2 + y^2
    r2sq = (x + μ - 1)^2 + y^2
    r1_3pow = r1sq * sqrt(r1sq)
    r2_3pow = r2sq * sqrt(r2sq)
    
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

    r1sq = (x + μ)^2 + y^2 + z^2
    r2sq = (x + μ - 1)^2 + y^2 + z^2
    r1_3pow = r1sq * sqrt(r1sq)
    r2_3pow = r2sq * sqrt(r2sq)
    
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
    
    # 计算状态导数
    du[1] = vx
    du[2] = vy
    du[3] = x + 2.0 * vy - μ_comp * dx1 / r1³ - μ * dx2 / r2³
    du[4] = y - 2.0 * vx - μ_comp * y / r1³ - μ * y / r2³
    
    # 存储STM导数
    
    dΦ = reshape(view(du, 5:20), 4, 4)
    mul!(dΦ, Df, Φ)
    
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
    r1, r2 = sqrt(r1sq), sqrt(r2sq)
    r1_3, r2_3 = r1sq * r1, r2sq * r2
    r1_5, r2_5 = r1_3 * r1sq, r2_3 * r2sq
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

    du[1:6] .= (vx, vy, vz,
        x + 2vy - μ1 * dx1 / r1_3 - μ * dx2 / r2_3,
        y - 2vx - μ1 * y / r1_3 - μ * y / r2_3,
        -μ1 * z / r1_3 - μ * z / r2_3)
    dΦ = reshape(view(du, 7:42), 6, 6)
    mul!(dΦ, Df, Φ)
    return nothing
end

"""
    bcr4bp_eqm!(du, u, p, t) -> nothing

地月旋转系 BCR4BP 的维度分发入口，支持 4 维平面、6 维空间、20 维平面+STM 与
42 维空间+STM 状态。

# Arguments
- `du`: 写入状态导数的预分配数组。
- `u`: 4、6、20 或 42 元状态向量。
- `p`: `Bcr4bp_Aux()` 返回的地月—太阳模型参数。
- `t`: 无量纲时间，用于计算太阳相位。

# Returns
返回 `nothing`，结果原位写入 `du`；不支持的状态长度抛出 `ArgumentError`。
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
    r1, r2, rs = sqrt(r1sq), sqrt(r2sq), sqrt(rssq)
    r1_3, r2_3, rs_3 = r1sq * r1, r2sq * r2, rssq * rs
    r1_5, r2_5, rs_5 = r1_3 * r1sq, r2_3 * r2sq, rs_3 * rssq
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
    dΦ = reshape(view(du, 5:20), 4, 4)
    mul!(dΦ, Df, Φ)
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
    dΦ = reshape(view(du, 7:42), 6, 6)
    mul!(dΦ, Df, Φ)
    return nothing
end


# 事件函数集成
include("dynamic_energy.jl")
include("dynamic_events.jl")
include("rotation.jl")

end
