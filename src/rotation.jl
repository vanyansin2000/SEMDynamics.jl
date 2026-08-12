"""
    旋转系柱坐标 / 月心惯性系坐标的相互转换。
    作为Utils.jl的一部分
"""

export cr3bp_inertial_to_rotating , cr3bp_rotating_to_inertial  
export polar_angle_p2 , unwrap_phase

"""
    cr3bp_inertial_to_rotating(μ, tt, state; center=:p1)
    cr3bp_inertial_to_rotating(μ, tt, state; center=:p2)

圆形限制性三体问题中从惯性系到旋转系的坐标转换。

将惯性坐标系中的状态向量转换到以指定主天体为中心的旋转坐标系中。

# 参数
- `μ::Real`: 系统质量比，μ = m₂/(m₁ + m₂)
- `tt::Real`: 时间或旋转角度（弧度）
- `state::SVector{4, T}`: 4维状态向量 [x, y, vx, vy]
- `states::AbstractMatrix`: N×4 状态矩阵，用于批量转换
- `center::Symbol`: 中心天体，可选 `:p1` 或 `:p2`，默认为 `:p1`

# 返回值
- 单个状态转换：返回 `SVector{4, T}`
- 批量状态转换：返回 `Matrix`

支持以P1或P2为中心的坐标系转换，自动处理坐标偏移。
"""
function cr3bp_inertial_to_rotating end

cr3bp_inertial_to_rotating(μ::Real , tt::Real, state::AbstractVector;  args...) = 
            cr3bp_inertial_to_rotating(μ::Real , tt::Real, SVector{4}(state);  args...)

# 单一状态向量版本
function cr3bp_inertial_to_rotating(μ::Real , tt::Real, state::SVector{4, T}; center = :p1) where T <: Real
    s, c = sincos(tt)
    x, y, vx, vy = state

    if center == :p1
        [
        x*c - μ + y*s,       # x_rot
        y*c - x*s,           # y_rot
        vx*c + y*c + vy*s - x*s,  # vx_rot
        vy*c - x*c - vx*s - y*s   # vy_rot
        ]
    elseif center == :p2
        [
        x*c - μ + y*s + 1.0,  # x_rot (P2中心调整)
        y*c - x*s,           # y_rot
        vx*c + y*c + vy*s - x*s,  # vx_rot
        vy*c - x*c - vx*s - y*s   # vy_rot
        ]
    else
        try
          error("the center can be :p1 or :p2 as Symbol Type")
        catch err
          showerror(stdout, err)
        end
        nothing 
    end
end


"""
    cr3bp_rotating_to_inertial(μ, tt, state; center=:p1)
    cr3bp_rotating_to_inertial(μ, tt, states; center=:p1)

圆形限制性三体问题中从旋转系到惯性系的坐标转换。

将旋转坐标系中的状态向量转换到惯性坐标系中。

# 参数
- `μ::Real`: 系统质量比，μ = m₂/(m₁ + m₂)
- `tt::Real`: 时间或旋转角度（弧度）
- `state::SVector{4, T}`: 4维状态向量 [x, y, vx, vy]
- `states::AbstractMatrix`: N×4 状态矩阵，用于批量转换
- `center::Symbol`: 中心天体，可选 `:p1` 或 `:p2`，默认为 `:p1`

# 返回值
- 单个状态转换：返回 `SVector{4, T}`
- 批量状态转换：返回 `Matrix`

此为 `cr3bp_inertial_to_rotating` 的逆变换，支持相同的参数选项。
"""
function cr3bp_rotating_to_inertial end

cr3bp_rotating_to_inertial(μ::Real , tt::Real, state::AbstractVector; args...) = 
            cr3bp_rotating_to_inertial(μ::Real , tt::Real, SVector{4}(state); args...)

# 单一状态向量版本
function cr3bp_rotating_to_inertial(μ::Real , tt::Real, state::SVector{4, T}; center = :p1) where T <: Real
    s, c = sincos(tt)
    x, y, vx, vy = state
    
       if center == :p1
    [
        c*(μ + x) - y*s,          # x_inertial
        s*(μ + x) + y*c,           # y_inertial  
        vx*c - s*(μ + x) - y*c - vy*s,  # vx_inertial
        c*(μ + x) + vy*c + vx*s - y*s   # vy_inertial
    ]
    elseif center == :p2
    [
        c*(μ + x - 1.0) - y*s,          # x_inertial (P2中心调整)
        s*(μ + x - 1.0) + y*c,           # y_inertial
        vx*c - s*(μ + x - 1.0) - y*c - vy*s,  # vx_inertial
        c*(μ + x - 1.0) + vy*c + vx*s - y*s   # vy_inertial
    ]
    else
        try
          error("the center can be :p1 or :p2 as Symbol Type")
        catch err
          showerror(stdout, err)
        end
        nothing 
    end
end

"计算旋转系下相位 及其导数"
@inline function polar_angle_p2(μ, u)
    ρ, θ, ρdot, θdot = u
    r = 1 + ρ
    s, c = sincos(θ)
    
    x = r * c - μ
    y = r * s
    Ψ = atan(y, x - (1 - μ)) 

    num = -s*ρdot + (r^2 - r*c)*θdot
    den = r^2 - 2r*c + 1
    dΨ = num/den
    return Ψ , dΨ
end

"将sol中的相位角展开"
function unwrap_phase(μ , sol::ODESolution )
us = sol.u   # 每个元素是状态向量 u = (ρ, θ, ρ̇, θ̇, ...)

# 2) 旋转系相位角（包裹角）
Ψ = map((u)->polar_angle_p2(μ, u)[1] , us)

  # φ_R(t_i)
# 3) unwrap 得到累计相位 Ψ_unwrapped
Ψ_unwrapped = similar(Ψ)
Ψ_unwrapped[1] = Ψ[1]
for j in Iterators.drop(eachindex(Ψ), 1)
    i = prevind(Ψ, j)
    d = mod(Ψ[j] - Ψ[i] + π, 2π) - π
    Ψ_unwrapped[j] = Ψ_unwrapped[i] + d
end
    return Ψ_unwrapped
end
