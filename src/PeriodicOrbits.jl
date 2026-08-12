module PeriodicOrbits

using CairoMakie
using DifferentialEquations
using LinearAlgebra
using NonlinearSolve

using ..Dynamics

export AbstractPeriodicOrbit, PeriodicOrbit, generate_DRO, orbit_shooting, plotPO!

const DRO_REFERENCE_PERIOD = pi
const DRO_REFERENCE_GUESS = [1.1754, -0.4943]
const DEFAULT_TOLERANCE = 1e-12

"""
    AbstractPeriodicOrbit

所有周期轨道对象的抽象父类型。新的轨道族可通过定义其子类型并提供相同的
`x0`、`P` 和数值解语义接入绘图与后处理接口。
"""
abstract type AbstractPeriodicOrbit end

"""
    PeriodicOrbit{F,S} <: AbstractPeriodicOrbit

一个已经由打靶修正的周期轨道。

字段 `x0`、`C` 和 `P` 分别为旋转坐标系初值、Jacobi 常数和无量纲周期；
`f` 是动力学函数，`sol` 是在一个完整周期 `[0, P]` 上的 ODE 解。
"""
struct PeriodicOrbit{F,S} <: AbstractPeriodicOrbit
    x0::Vector{Float64}
    C::Float64
    P::Float64
    f::F
    sol::S
end

"""
    integrate_orbit(f, x0, P; solver=Vern7(), kwargs...) -> ODESolution

在 CR3BP 辅助参数下，从 `x0` 将动力学函数 `f` 积分一个正周期 `P`。关键字
参数会原样传给 `solve`，因而可覆盖误差容限或选择不同求解器。
"""
function integrate_orbit(
    f::Function,
    x0::AbstractVector{<:Real},
    P::Real;
    solver=Vern7(),
    kwargs...,
)
    P > 0 || throw(ArgumentError("周期 P 必须为正数，得到 $P"))
    problem = ODEProblem(f, Float64.(x0), (0.0, Float64(P)), Bcr4bp_Aux())
    return solve(problem, solver; dense=true, save_everystep=true, kwargs...)
end

"""
    make_periodic_orbit(f, x0, P; kwargs...) -> PeriodicOrbit

由已修正的初值构造 `PeriodicOrbit`，并重新积分一个周期以得到可用于绘图和
插值的解。此函数不会再次执行打靶修正。
"""
function make_periodic_orbit(
    f::Function,
    x0::AbstractVector{<:Real},
    P::Real;
    abstol::Real=DEFAULT_TOLERANCE,
    reltol::Real=DEFAULT_TOLERANCE,
    kwargs...,
)
    sol = integrate_orbit(f, x0, P; abstol, reltol, kwargs...)
    μ = Bcr4bp_Aux().EMRot.μ
    return PeriodicOrbit(Float64.(x0), compute_jacobi(x0, μ), Float64(P), f, sol)
end

"""
    orbit_shooting(residual!, initial_guess, parameters; tol=1e-12) -> Vector{Float64}

通用的非线性打靶求根接口。`residual!` 应遵循 SciML 的
`residual!(residual, unknowns, parameters)` 约定；本函数只负责构造并求解
`NonlinearProblem`。未来的轨道族可以复用该接口，并替换自己的残差和参数化。

若 Newton 迭代未收敛，将抛出错误而不是返回未经修正的“初值”。
"""
function orbit_shooting(
    residual!::Function,
    initial_guess::AbstractVector{<:Real},
    parameters;
    tol::Real=DEFAULT_TOLERANCE,
)
    problem = NonlinearProblem(residual!, Float64.(initial_guess), parameters)
    solution = solve(problem, NewtonRaphson(); abstol=tol, reltol=tol)
    SciMLBase.successful_retcode(solution) || error(
        "打靶未收敛（残差范数=$(norm(solution.resid))）。请提供更接近目标轨道的 seed，或减小 continuation_step。",
    )
    return Float64.(solution.u)
end

"""
    dro_half_period_residual!(residual, u, parameters)

平面 DRO 的半周期对称打靶残差。未知量 `u = [x, vy]` 表示
`[x, 0, 0, 0, vy, 0]`；在 `P/2` 时要求 `y = vx = 0`。周期固定在
`parameters.P`，因此该残差可用于按指定周期求解轨道。
"""
function dro_half_period_residual!(residual, u, parameters)
    x, vy = u
    P, dynamics, aux = parameters.P, parameters.dynamics, parameters.aux
    state = [x, 0.0, 0.0, 0.0, vy, 0.0]
    problem = ODEProblem(dynamics, state, (0.0, P / 2), aux)
    solution = solve(
        problem,
        Vern7();
        abstol=parameters.tol,
        reltol=parameters.tol,
        save_everystep=false,
    )
    final_state = solution.u[end]
    residual[1] = final_state[2]
    residual[2] = final_state[4]
    return nothing
end

"""
    orbit_shooting(dynamics, x_guess, vy_guess, P; tol=1e-12) -> Vector{Float64}

使用平面 DRO 的半周期对称条件，在给定周期 `P` 下修正 `x_guess` 和 `vy_guess`。
返回六维旋转坐标系初值 `[x, 0, 0, 0, vy, 0]`。
"""
function orbit_shooting(
    dynamics::Function,
    x_guess::Real,
    vy_guess::Real,
    P::Real;
    tol::Real=DEFAULT_TOLERANCE,
)
    P > 0 || throw(ArgumentError("周期 P 必须为正数，得到 $P"))
    parameters = (P=Float64(P), dynamics=dynamics, aux=Bcr4bp_Aux(), tol=Float64(tol))
    corrected = orbit_shooting(dro_half_period_residual!, [x_guess, vy_guess], parameters; tol)
    return [corrected[1], 0.0, 0.0, 0.0, corrected[2], 0.0]
end

"""
    _continue_dro_to_period(dynamics, P, seed, seed_period; continuation_step, tol)

从一个已知 DRO seed 沿周期参数走到 `P`。该实现仅在内存中保留最近两个修正点，
利用割线预测下一个初值后再打靶；它避免了 BifurcationKit 的全分支延拓以及任何
磁盘缓存，同时比每步都重复使用原始 seed 更稳定。
"""
function _continue_dro_to_period(
    dynamics::Function,
    P::Real,
    seed::AbstractVector{<:Real},
    seed_period::Real;
    continuation_step::Real,
    tol::Real,
)
    continuation_step > 0 || throw(ArgumentError("continuation_step 必须为正数"))
    length(seed) == 2 || throw(ArgumentError("DRO seed 必须是 [x, vy] 两个元素"))
    P > 0 || throw(ArgumentError("周期 P 必须为正数，得到 $P"))
    seed_period > 0 || throw(ArgumentError("seed_period 必须为正数"))

    nsteps = ceil(Int, abs(P - seed_period) / continuation_step)
    periods = range(Float64(seed_period), Float64(P); length=nsteps + 1)
    previous_u = Float64.(seed)
    previous_P = first(periods)
    penultimate_u = nothing
    penultimate_P = nothing

    # 先修正 seed 本身：用户提供的近似初值也会成为严格周期解。
    corrected_x0 = orbit_shooting(dynamics, previous_u[1], previous_u[2], previous_P; tol)
    previous_u = corrected_x0[[1, 5]]

    for current_P in Iterators.drop(periods, 1)
        guess = if isnothing(penultimate_u)
            previous_u
        else
            previous_u .+ (current_P - previous_P) / (previous_P - penultimate_P) .* (previous_u .- penultimate_u)
        end
        corrected_x0 = orbit_shooting(dynamics, guess[1], guess[2], current_P; tol)
        corrected_u = corrected_x0[[1, 5]]
        penultimate_u, penultimate_P = previous_u, previous_P
        previous_u, previous_P = corrected_u, current_P
    end
    return [previous_u[1], 0.0, 0.0, 0.0, previous_u[2], 0.0]
end

"""
    generate_DRO(; P=pi, seed=DRO_REFERENCE_GUESS, seed_period=pi,
                   continuation_step=0.05, tol=1e-12) -> PeriodicOrbit

生成具有目标无量纲周期 `P` 的平面远距逆行轨道（DRO）。默认 seed 是 2:1 DRO，
因此 `generate_DRO(P=pi)` 直接给出其打靶修正后的轨道；初值可由 `.x0` 取得。

当 `P` 不等于 `seed_period` 时，函数采用短步长、割线预测的周期参数延拓，并在
每一步执行固定周期打靶。整个过程不读写轨道数据文件，也不会建立完整轨道族。
若已有其他轨道族的近似解，可将其二维 `[x, vy]` seed 和对应 `seed_period` 传入。
"""
function generate_DRO(
    ;
    P::Real=DRO_REFERENCE_PERIOD,
    seed::AbstractVector{<:Real}=DRO_REFERENCE_GUESS,
    seed_period::Real=DRO_REFERENCE_PERIOD,
    continuation_step::Real=0.05,
    tol::Real=DEFAULT_TOLERANCE,
)
    x0 = _continue_dro_to_period(
        cr3bp_eqm!,
        P,
        seed,
        seed_period;
        continuation_step,
        tol,
    )
    return make_periodic_orbit(cr3bp_eqm!, x0, P; abstol=tol, reltol=tol)
end

"""
    plotPO!(ax, orbit; kwargs...)

在 Makie 轴 `ax` 上绘制单条周期轨道的三维坐标曲线。额外关键字参数将传给
`lines!`，例如 `color` 或 `linewidth`。
"""
function plotPO!(ax, orbit::AbstractPeriodicOrbit; kwargs...)
    times = range(orbit.sol.t[1], orbit.sol.t[end]; length=2000)
    states = orbit.sol(times).u
    orb = lines!(ax, getindex.(states, 1), getindex.(states, 2), getindex.(states, 3); kwargs...)
    return orb
end

"""
    plotPO!(ax, orbits::AbstractVector{<:AbstractPeriodicOrbit}; kwargs...)

在同一 Makie 轴上批量绘制周期轨道。各轨道以 `NaN` 分隔，避免 Makie 将不同
轨道误连为一条线。
"""
function plotPO!(ax, orbits::AbstractVector{<:AbstractPeriodicOrbit}; kwargs...)
    xs, ys, zs = Float64[], Float64[], Float64[]
    for orbit in orbits
        times = range(orbit.sol.t[1], orbit.sol.t[end]; length=1000)
        states = orbit.sol(times).u
        append!(xs, getindex.(states, 1)); push!(xs, NaN)
        append!(ys, getindex.(states, 2)); push!(ys, NaN)
        append!(zs, getindex.(states, 3)); push!(zs, NaN)
    end
    orb = lines!(ax, xs, ys, zs; kwargs...)
    return orb
end

end # module PeriodicOrbits
