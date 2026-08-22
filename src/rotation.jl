"""
CR3BP 旋转系、主天体中心惯性系与相位角工具。
"""

export cr3bp_inertial_to_rotating , cr3bp_rotating_to_inertial  
export polar_angle_p2 , unwrap_phase

"""
    cr3bp_inertial_to_rotating(μ, t, state; center=:p1) -> SVector{4}

将以指定主天体为原点的平面惯性状态转换到 CR3BP 旋转系。

# Arguments
- `μ`: 系统质量比 ``m₂/(m₁+m₂)``。
- `t`: 无量纲时间，同时也是旋转角（弧度）。
- `state`: 惯性状态 `[x, y, vx, vy]`。

# Keywords
- `center=:p1`: 惯性状态原点，可取 `:p1` 或 `:p2`。

# Returns
返回旋转系静态向量 `[x, y, vx, vy]`。输入中心非法时抛出 `ArgumentError`。
"""
function cr3bp_inertial_to_rotating end

function cr3bp_inertial_to_rotating(μ::Real, t::Real, state::AbstractVector; kwargs...)
    length(state) == 4 || throw(ArgumentError("state must contain [x, y, vx, vy]."))
    return cr3bp_inertial_to_rotating(μ, t, SVector{4}(state); kwargs...)
end

# 单一状态向量版本
function cr3bp_inertial_to_rotating(
    μ::Real,
    t::Real,
    state::SVector{4,T};
    center::Symbol=:p1,
) where {T<:Real}
    center in (:p1, :p2) || throw(ArgumentError("center must be :p1 or :p2."))
    s, c = sincos(t)
    x, y, vx, vy = state
    x_offset = center === :p1 ? -μ : 1 - μ
    return SVector(
        x * c + y * s + x_offset,
        y * c - x * s,
        vx * c + y * c + vy * s - x * s,
        vy * c - x * c - vx * s - y * s,
    )
end


"""
    cr3bp_rotating_to_inertial(μ, t, state; center=:p1) -> SVector{4}

将平面 CR3BP 旋转状态转换到以指定主天体为原点的惯性系。

# Arguments
- `μ`: 系统质量比。
- `t`: 无量纲时间/旋转角。
- `state`: 旋转系状态 `[x, y, vx, vy]`。

# Keywords
- `center=:p1`: 输出惯性坐标系的原点，可取 `:p1` 或 `:p2`。

# Returns
返回惯性系静态向量 `[x, y, vx, vy]`；该函数是
`cr3bp_inertial_to_rotating` 的逆变换。
"""
function cr3bp_rotating_to_inertial end

function cr3bp_rotating_to_inertial(μ::Real, t::Real, state::AbstractVector; kwargs...)
    length(state) == 4 || throw(ArgumentError("state must contain [x, y, vx, vy]."))
    return cr3bp_rotating_to_inertial(μ, t, SVector{4}(state); kwargs...)
end

# 单一状态向量版本
function cr3bp_rotating_to_inertial(
    μ::Real,
    t::Real,
    state::SVector{4,T};
    center::Symbol=:p1,
) where {T<:Real}
    center in (:p1, :p2) || throw(ArgumentError("center must be :p1 or :p2."))
    s, c = sincos(t)
    x, y, vx, vy = state
    x_relative = center === :p1 ? μ + x : μ + x - 1
    return SVector(
        c * x_relative - y * s,
        s * x_relative + y * c,
        vx * c - s * x_relative - y * c - vy * s,
        c * x_relative + vy * c + vx * s - y * s,
    )
end

"""
    polar_angle_p2(μ, u) -> (angle, angular_rate)

把极坐标形式 `u = [ρ, θ, ρ̇, θ̇]` 转换为相对第二主天体的方位角及其导数。

# Returns
返回 `(Ψ, dΨ)`，分别为包裹到 `[-π, π]` 的角度和无量纲角速度。
"""
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

"""
    unwrap_phase(μ, sol) -> Vector

计算 ODE 解中每个状态相对第二主天体的方位角，并消除相邻样本之间的 ``2π``
跳变。

# Arguments
- `μ`: CR3BP 质量比。
- `sol`: 状态按 `[ρ, θ, ρ̇, θ̇, ...]` 排列的 `ODESolution`。

# Returns
返回与 `sol.u` 等长的累计连续相位向量。
"""
function unwrap_phase(μ, sol::ODESolution)
    Ψ = map(u -> first(polar_angle_p2(μ, u)), sol.u)
    isempty(Ψ) && return Ψ
    Ψ_unwrapped = similar(Ψ)
    Ψ_unwrapped[1] = Ψ[1]
    for j in Iterators.drop(eachindex(Ψ), 1)
        i = prevind(Ψ, j)
        d = mod(Ψ[j] - Ψ[i] + π, 2π) - π
        Ψ_unwrapped[j] = Ψ_unwrapped[i] + d
    end
    return Ψ_unwrapped
end
