

export compute_ε2 , compute_ε2_dot , compute_jacobi
export solve_L1_L2_x


"""
    compute_ε2(u, t, μ) -> Real

计算旋转系状态相对第二主天体的惯性两体比机械能。

# Arguments
- `u`: 长度为 4、6、20 或 42 的状态；STM 部分被忽略。
- `t`: 无量纲时间/旋转角。
- `μ`: 第二主天体的 CR3BP 质量比。

# Returns
返回 `ε₂ = ‖v₂‖²/2 - μ/‖r₂‖`。负值表示相对于该两体近似为束缚状态。
"""
@inline function compute_ε2(u , t , μ)

    dimension = _state_dimension(u)
    state = view(u, 1:2dimension)
    u_inertial = cr3bp_rotating_to_inertial(μ, t, state; center=:p2)
    r2 = norm(view(u_inertial, 1:dimension))
    v2_pow_2 = sum(abs2, view(u_inertial, dimension+1:2dimension))

    return  0.5*v2_pow_2 - μ/r2
end

"""
    compute_ε2_dot(u, t, μ, p; dynamics=cr3bp_eqm!) -> Real

利用状态导数与自动微分计算 `compute_ε2` 沿轨迹的时间导数。
默认采用 CR3BP；四体轨迹应传入 `dynamics=bcr4bp_eqm!`。

# Arguments
- `u`: 长度为 4、6、20 或 42 的状态；STM 部分被忽略。
- `t`: 当前无量纲时间。
- `μ`: CR3BP 质量比。
- `p`: 动力学参数，通常为 `Bcr4bp_Aux()`。

# Returns
返回标量 `dε₂/dt`。
"""
@inline function compute_ε2_dot(u, t, μ, p; dynamics=cr3bp_eqm!)
    dimension = _state_dimension(u)
    state = collect(view(u, 1:2dimension))
    du = similar(state)
    dynamics(du, state, p, t)
    # 旋转变换不改变 r、v 的范数，因此 ε₂ 没有显式时间项。
    ∂ε_∂u = ForwardDiff.gradient(state -> compute_ε2(state, t, μ), state)
    return dot(∂ε_∂u,  du)
end


"""
    compute_jacobi(u, μ) -> Real

计算平面或空间 CR3BP 状态的 Jacobi 常数。

# Arguments
- `u`: `[x, y, vx, vy]` 或 `[x, y, z, vx, vy, vz]`，可为普通向量或
  `StaticArrays.SVector`。
- `μ`: CR3BP 质量比。

# Returns
返回与输入数值类型兼容的 Jacobi 常数标量。

# Throws
普通向量长度不是 4 或 6 时抛出 `ArgumentError`。
"""
@inline function compute_jacobi(u::SVector{4} , μ)
    x, y, vx, vy = u
    r1 = sqrt((x + μ)^2    + y^2 )
    r2 = sqrt((x + μ - 1)^2 + y^2 )
    JC = -(vx^2 + vy^2) + (x^2 + y^2) + 2 * (1 - μ) / r1 + 2 * μ / r2
    return JC
end

@inline function compute_jacobi(u::SVector{6} , μ)
    x, y , z, vx , vy , vz = u
    r1 = sqrt((x + μ)^2    + y^2 + z^2)
    r2 = sqrt((x + μ - 1)^2 + y^2 + z^2)
    v2 = vx^2 + vy^2 + vz^2
    return x^2 + y^2 + 2*(1-μ)/r1 + 2*μ/r2 - v2
end

function compute_jacobi(u::AbstractVector, μ)
    length(u) == 4 && return compute_jacobi(SVector{4}(u), μ)
    length(u) == 6 && return compute_jacobi(SVector{6}(u), μ)
    throw(ArgumentError("Jacobi constant requires a 4- or 6-element state."))
end


"""
    Hamiltonian(u, θ, aux) -> (HEM, HEMdot, Γ)

计算四体模型下瞬时地月 Jacobi/Hamiltonian 型积分量。

状态格式：
- 2D：`u = [x, y, vx, vy]`
- 3D：`u = [x, y, z, vx, vy, vz]`

`θ` 为太阳在地月旋转坐标系中的相位角，单位 rad。
"""
@inline function Hamiltonian_EM(u::SVector{4}, θ::Real, aux)
    x, y, vx, vy = u

    μ  = aux.EMRot.μ
    ms = aux.EMRot.mus
    as = aux.EMRot.as
    ws = aux.EMRot.ws

    cθ, sθ = cos(θ), sin(θ)

    r1 = sqrt((x + μ)^2 + y^2)
    r2 = sqrt((x + μ - 1)^2 + y^2)
    rs = sqrt((x - as*cθ)^2 + (y - as*sθ)^2)

    Γ = (
        (x^2 + y^2) / 2
        + (1 - μ) / r1
        + μ / r2
        + ms / rs
        - ms / as^2 * (x*cθ + y*sθ)
    )

    HEM = -(vx^2 + vy^2) + 2Γ

    HEMdot = (
        -2 * ws * ms * as
        * (x*sθ - y*cθ)
        * (inv(rs^3) - inv(as^3))
    )

    return HEM, HEMdot, Γ
end

@inline function Hamiltonian_EM(u::SVector{6}, θ::Real, aux )
    x, y, z, vx, vy, vz = u

    μ  = aux.EMRot.μ
    ms = aux.EMRot.mus
    as = aux.EMRot.as
    ws = aux.EMRot.ws

    cθ, sθ = cos(θ), sin(θ)

    r1 = sqrt((x + μ)^2 + y^2 + z^2)
    r2 = sqrt((x + μ - 1)^2 + y^2 + z^2)
    rs = sqrt((x - as*cθ)^2 + (y - as*sθ)^2 + z^2)

    # 旋转离心势只包含 x²+y²，不包含 z²
    Γ = (
        (x^2 + y^2) / 2
        + (1 - μ) / r1
        + μ / r2
        + ms / rs
        - ms / as^2 * (x*cθ + y*sθ)
    )

    HEM = -(vx^2 + vy^2 + vz^2) + 2Γ

    # 太阳位于参考 xy 平面时，此表达式与二维形式一致
    HEMdot = (
        -2 * ws * ms * as
        * (x*sθ - y*cθ)
        * (inv(rs^3) - inv(as^3))
    )

    return HEM, HEMdot, Γ
end

function Hamiltonian(u::AbstractVector, θ::Real, aux ; coordinate = :EMRot)
    if coordinate == :EMRot
        if length(u) == 4
            return Hamiltonian_EM(SVector{4}(u), θ, aux)
        elseif length(u) == 6
            return Hamiltonian_EM(SVector{6}(u), θ, aux)
        end
    end

    throw(ArgumentError(
        "State must contain 4 elements [x,y,vx,vy] " *
        "or 6 elements [x,y,z,vx,vy,vz]."
    ))
end


"""
    f_collinear(μ, x) -> Real

返回共线平衡点方程在横坐标 `x` 处的残差；主要供 `solve_L1_L2_x` 使用。
"""
function f_collinear(μ, x)
    r1 = abs(x + μ)
    r2 = abs(x - (1 - μ))
    x - (1 - μ)*(x + μ)/r1^3 - μ*(x - (1 - μ))/r2^3
end
"""
    fprime_collinear(μ, x) -> Real

返回共线平衡点残差对 `x` 的解析导数。
"""
function fprime_collinear(μ, x)
    r1 = abs(x + μ)
    r2 = abs(x - (1 - μ))
    return 1 + 2(1 - μ) / r1^3 + 2μ / r2^3
end

"""
    solve_L1_L2_x(μ; which=:L1, tol=1e-14, maxit=100) -> Real

用 Newton 法求 CR3BP 的 L1 或 L2 共线平衡点横坐标。

# Arguments
- `μ`: CR3BP 质量比。

# Keywords
- `which=:L1`: 选择 `:L1` 或 `:L2`。
- `tol=1e-14`: Newton 步长的收敛阈值。
- `maxit=100`: 最大迭代次数。

# Returns
返回所选平动点在地月旋转系中的无量纲横坐标。

# Throws
`which` 非法或在 `maxit` 次迭代内未收敛时抛出异常。
"""
function solve_L1_L2_x(μ; which=:L1, tol=1e-14, maxit=100)
    which in (:L1, :L2) || throw(ArgumentError("which must be :L1 or :L2."))
    δ = (μ/3)^(1/3)
    x = which === :L1 ? (1 - μ - δ) : (1 - μ + δ)
    for _ in 1:maxit
        fx = f_collinear(μ, x)
        dfx = fprime_collinear(μ, x)
        step = fx/dfx
        x -= step
        if abs(step) < tol
            return x
        end
    end
    error("Newton 未收敛：$which")
end

