

export compute_ε2 , compute_ε2_dot , compute_jacobi
export solve_L1_L2_x


"""计算相对第二主天体的惯性两体比机械能 `ε₂`。"""
@inline function compute_ε2(u , t , μ)

    u_inertial = cr3bp_rotating_to_inertial(μ , t, u; center = :p2) 

    r2 = hypot(u_inertial[1]  , u_inertial[2])
    v2_pow_2 = u_inertial[3]^2+ u_inertial[4]^2

    return  0.5*v2_pow_2 - μ/r2
end

"""用动力学导数和自动微分计算 `ε₂` 的时间导数。"""
@inline function compute_ε2_dot(u , t , μ ,p)
    du = similar(u) 
    cr3bp_eqm2D!(du , u , p , t)
    ∂ε_∂u = ForwardDiff.gradient((u)->compute_ε2(u , t , μ), u)
    return dot(∂ε_∂u,  du)
end


"""计算四维平面 CR3BP 状态的 Jacobi 常数。"""
@inline function compute_jacobi(u::SVector{4} , μ)
    r1 = sqrt((x + μ)^2    + y^2 )
    r2 = sqrt((x + μ - 1)^2 + y^2 )
    JC = -(u[3]^2 + u[4]^2) + (u[1]^2 + u[2]^2) + 2 * (1 - μ) / r1 + 2 * μ / r2
    return JC
end

"""计算六维空间 CR3BP 状态的 Jacobi 常数。"""
@inline function compute_jacobi(u::SVector{6} , μ)
    x, y , z, vx , vy , vz = u
    r1 = sqrt((x + μ)^2    + y^2 + z^2)
    r2 = sqrt((x + μ - 1)^2 + y^2 + z^2)
    v2 = vx^2 + vy^2 + vz^2
    return x^2 + y^2 + 2*(1-μ)/r1 + 2*μ/r2 - v2
end

"""将一般向量转换为静态向量后计算 Jacobi 常数。"""
compute_jacobi(u::AbstractVector , μ) = compute_jacobi(SVector{length(u)}(u), μ)


"""返回共线平衡点方程在坐标 `x` 处的函数值。"""
function f_collinear(μ, x)
    r1 = abs(x + μ)
    r2 = abs(x - (1 - μ))
    x - (1 - μ)*(x + μ)/r1^3 - μ*(x - (1 - μ))/r2^3
end
"""以中心差分近似共线平衡点方程的导数。"""
fprime_collinear(μ, x; h=1e-8) = (f_collinear(μ, x + h) - f_collinear(μ, x - h))/(2h)

"""
    solve_L1_L2_x(μ; which=:L1, tol=1e-14, maxit=100)

用 Newton 法求 CR3BP 的 L1 或 L2 共线平衡点横坐标。
"""
function solve_L1_L2_x(μ; which=:L1, tol=1e-14, maxit=100)
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

