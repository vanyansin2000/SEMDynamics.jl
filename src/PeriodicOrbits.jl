"""Differential correction, continuation, storage, and plotting of CR3BP periodic orbits."""
module PeriodicOrbits

using CairoMakie
using DifferentialEquations
using LinearAlgebra
using NonlinearSolve

using ..Dynamics

export AbstractPeriodicOrbit, PeriodicOrbit
export generate_DRO, generate_halo, generate_nrho_9_2
export orbit_shooting, plotPO!

const DEFAULT_INTEGRATION_TOLERANCE = 1e-12
const DEFAULT_SHOOTING_TOLERANCE = 1e-10
const DEFAULT_HALO_CONTINUATION_STEP = 0.005
const DEFAULT_MINIMUM_HALO_Z = 1e-8
const LUNAR_SYNODIC_PERIOD = 2pi / abs(Bcr4bp_Aux().EMRot.ws)
const NRHO_9_2_PERIOD = 2LUNAR_SYNODIC_PERIOD / 9

struct OrbitReferenceData
    C::Union{Nothing,Float64}
    P::Float64
    seed::Vector{Float64}
end

# Every reference is stored as a complete six-state symmetry-plane crossing
# [x, 0, z, 0, vy, 0]. PeriodicOrbit.C is recomputed after correction.
const ORBIT_REFERENCE_DATA = Dict{NTuple{3,Symbol},OrbitReferenceData}(
    (:dro, :default, :default) => OrbitReferenceData(
        nothing, pi, [1.1754, 0.0, 0.0, 0.0, -0.4943, 0.0],
    ),
    (:halo, :southern, :L1) => OrbitReferenceData(
        nothing,
        2.74494694510099,
        [0.823379815960398, 0.0, -0.016954328, 0.0, 0.130975375154282, 0.0],
    ),
    (:halo, :northern, :L1) => OrbitReferenceData(
        nothing,
        2.744967911,
        [0.823379821, 0.0, 0.017046171, 0.0, 0.131024353, 0.0],
    ),
    (:halo, :southern, :L2) => OrbitReferenceData(
        nothing,
        1.767711412,
        [0.987068513, 0.0, 0.017543723, 0.0, 1.14211561, 0.0],
    ),
    (:halo, :northern, :L2) => OrbitReferenceData(
        nothing,
        1.75900893,
        [0.987068596, 0.0, -0.017201977, 0.0, 1.153880419, 0.0],
    ),
)
const DRO_REFERENCE_DATA = ORBIT_REFERENCE_DATA[(:dro, :default, :default)]
const DRO_REFERENCE_PERIOD = DRO_REFERENCE_DATA.P
const DRO_REFERENCE_GUESS = DRO_REFERENCE_DATA.seed

"""
    AbstractPeriodicOrbit

所有周期轨道对象的抽象父类型。新的轨道族可通过定义其子类型并提供相同的
`x0`、`P` 和数值解语义接入绘图与后处理接口。
"""
abstract type AbstractPeriodicOrbit end

"""
    PeriodicOrbit{F,S} <: AbstractPeriodicOrbit

一个已经由打靶修正的周期轨道。

# Fields
- `x0::Vector{Float64}`: 旋转坐标系中的六维对称面初值。
- `C::Float64`: 初值对应的 Jacobi 常数。
- `P::Float64`: 无量纲完整周期。
- `f`: 用于传播轨道的原位动力学函数。
- `sol`: 在完整周期 `[0, P]` 上、带稠密插值的 ODE 解。
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

# Returns
返回带稠密插值的 `ODESolution`；`P ≤ 0` 时抛出 `ArgumentError`。
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

# Keywords
- `abstol`, `reltol`: 完整周期积分容差，默认均为 `1e-12`。
- 其余关键字传给 `integrate_orbit`。

# Returns
返回包含初值、Jacobi 常数、周期、动力学函数和 ODE 解的 `PeriodicOrbit`。
"""
function make_periodic_orbit(
    f::Function,
    x0::AbstractVector{<:Real},
    P::Real;
    abstol::Real=DEFAULT_INTEGRATION_TOLERANCE,
    reltol::Real=DEFAULT_INTEGRATION_TOLERANCE,
    kwargs...,
)
    sol = integrate_orbit(f, x0, P; abstol, reltol, kwargs...)
    μ = Bcr4bp_Aux().EMRot.μ
    return PeriodicOrbit(Float64.(x0), compute_jacobi(x0, μ), Float64(P), f, sol)
end

"""
    orbit_shooting(residual!, initial_guess, parameters; tol=1e-10)
        -> Vector{Float64}

通用的非线性打靶求根接口。`residual!` 应遵循 SciML 的
`residual!(residual, unknowns, parameters)` 约定；本函数只负责构造并求解
`NonlinearProblem`。未来的轨道族可以复用该接口，并替换自己的残差和参数化。

若 Newton 迭代未收敛，将抛出错误而不是返回未经修正的“初值”。

# Arguments
- `residual!`: SciML 原位残差 `residual!(residual, unknowns, parameters)`。
- `initial_guess`: 非线性未知量初值。
- `parameters`: 原样传给残差的参数对象。

# Returns
返回 Newton 修正后的未知量向量。
"""
function orbit_shooting(
    residual!::Function,
    initial_guess::AbstractVector{<:Real},
    parameters;
    tol::Real=DEFAULT_SHOOTING_TOLERANCE,
)
    problem = NonlinearProblem(residual!, Float64.(initial_guess), parameters)
    solution = solve(problem, NewtonRaphson(); abstol=tol, reltol=tol)
    SciMLBase.successful_retcode(solution) || error(
        "打靶未收敛（残差范数=$(norm(solution.resid))）。请提供更接近目标轨道的 seed，或减小 continuation_step。",
    )
    return Float64.(solution.u)
end

const SYMMETRY_RESIDUAL_INDICES = (2, 4, 6)

"""
    _state_from_unknowns(u) -> Vector

把统一打靶未知量 `[x, z, vy]` 展开为对称面六维状态 `[x, 0, z, 0, vy, 0]`。
"""
function _state_from_unknowns(u)
    length(u) == 3 || throw(ArgumentError("Shooting unknowns must be [x, z, vy]."))
    x, z, vy = u
    zero_x = zero(x)
    return [x, zero_x, z, zero_x, vy, zero_x]
end

"""
    _unknowns_from_state(state) -> Vector{Float64}

从完整六维对称面状态提取统一打靶未知量 `[x, z, vy]`。
"""
function _unknowns_from_state(state)
    length(state) == 6 || throw(ArgumentError(
        "Periodic-orbit seeds must contain the full six-state [x, y, z, vx, vy, vz].",
    ))
    return Float64[state[1], state[3], state[5]]
end

"""
    _symmetric_half_period_residual!(residual, u, parameters) -> nothing

传播半周期并原位写入 `[y(P/2), vx(P/2), vz(P/2)]` 对称残差。
"""
function _symmetric_half_period_residual!(residual, u, parameters)
    state = _state_from_unknowns(u)
    problem = ODEProblem(
        parameters.dynamics,
        state,
        (0.0, parameters.P / 2),
        parameters.aux,
    )
    integration_tol = hasproperty(parameters, :integration_tol) ?
                      parameters.integration_tol : parameters.tol
    solution = solve(
        problem,
        Vern7();
        abstol=integration_tol,
        reltol=integration_tol,
        save_everystep=false,
    )
    final_state = solution.u[end]
    for (residual_index, state_index) in enumerate(SYMMETRY_RESIDUAL_INDICES)
        residual[residual_index] = final_state[state_index]
    end
    return nothing
end

"""
    _shoot_symmetric_orbit(dynamics, seed_state, P; tol, z_sign=0, min_abs_z=0)
        -> Vector{Float64}

在固定周期下修正一个六维对称轨道初值，并可检查非平面分支的符号和最小振幅。
"""
function _shoot_symmetric_orbit(
    dynamics::Function,
    seed_state::AbstractVector{<:Real},
    P::Real;
    tol::Real,
    z_sign::Int=0,
    min_abs_z::Real=0.0,
)
    length(seed_state) == 6 || throw(ArgumentError(
        "Periodic-orbit seeds must contain the full six-state [x, y, z, vx, vy, vz].",
    ))
    P > 0 || throw(ArgumentError("Orbit period P must be positive; got $P."))
    min_abs_z >= 0 || throw(ArgumentError("min_abs_z must be nonnegative."))
    z_sign in (-1, 0, 1) || throw(ArgumentError("z_sign must be -1, 0, or 1."))

    parameters = (
        P=Float64(P),
        dynamics=dynamics,
        aux=Bcr4bp_Aux(),
        tol=Float64(tol),
        integration_tol=Float64(DEFAULT_INTEGRATION_TOLERANCE),
    )
    corrected = orbit_shooting(
        _symmetric_half_period_residual!,
        _unknowns_from_state(seed_state),
        parameters;
        tol,
    )
    if z_sign != 0 && (sign(corrected[2]) != z_sign || abs(corrected[2]) < min_abs_z)
        error(
            "The shooting corrector converged to the planar z = 0 branch. " *
            "Provide a non-planar seed closer to the requested orbit or reduce continuation_step.",
        )
    end
    return _state_from_unknowns(corrected)
end

"""
    _continue_symmetric_orbit(dynamics, P, seed, seed_period;
                              continuation_step, tol, z_sign=0, min_abs_z=0)

Shared period-continuation engine for all symmetric periodic orbits. Every seed
is a full six-state, and every correction uses the same unknowns `[x, z, vy]`
and residuals `[y(P/2), vx(P/2), vz(P/2)]`. The optional `z_sign` and
`min_abs_z` guard a requested non-planar branch against convergence to the
coexisting planar solution.
"""
function _continue_symmetric_orbit(
    dynamics::Function,
    P::Real,
    seed::AbstractVector{<:Real},
    seed_period::Real;
    continuation_step::Real,
    tol::Real,
    z_sign::Int=0,
    min_abs_z::Real=0.0,
)
    continuation_step > 0 || throw(ArgumentError("continuation_step must be positive."))
    P > 0 || throw(ArgumentError("Orbit period P must be positive; got $P."))
    seed_period > 0 || throw(ArgumentError("seed_period must be positive."))
    length(seed) == 6 || throw(ArgumentError(
        "Periodic-orbit seeds must contain the full six-state [x, y, z, vx, vy, vz].",
    ))

    nsteps = ceil(Int, abs(P - seed_period) / continuation_step)
    periods = range(Float64(seed_period), Float64(P); length=nsteps + 1)
    previous_state = Float64.(seed)
    previous_P = first(periods)
    penultimate_u = nothing
    penultimate_P = nothing

    previous_state = _shoot_symmetric_orbit(
        dynamics, previous_state, previous_P; tol, z_sign, min_abs_z,
    )
    previous_u = _unknowns_from_state(previous_state)

    for current_P in Iterators.drop(periods, 1)
        guess_u = if isnothing(penultimate_u)
            copy(previous_u)
        else
            previous_u .+ (current_P - previous_P) / (previous_P - penultimate_P) .* (previous_u .- penultimate_u)
        end
        guess_state = _state_from_unknowns(guess_u)
        corrected_state = _shoot_symmetric_orbit(
            dynamics, guess_state, current_P; tol, z_sign, min_abs_z,
        )
        corrected_u = _unknowns_from_state(corrected_state)
        penultimate_u, penultimate_P = previous_u, previous_P
        previous_u, previous_P = corrected_u, current_P
        previous_state = corrected_state
    end
    return previous_state
end

"""
    dro_half_period_residual!(residual, u, parameters)

兼容旧名称的对称周期轨道残差。未知量统一为 `u = [x, z, vy]`，并在
`P/2` 时要求 `y = vx = vz = 0`。即使默认 DRO 是平面的，`z` 仍参与
打靶和延拓。
"""
function dro_half_period_residual!(residual, u, parameters)
    return _symmetric_half_period_residual!(residual, u, parameters)
end

"""
    orbit_shooting(dynamics, state_guess, P; tol=1e-10,
                   z_sign=0, min_abs_z=0) -> Vector{Float64}

使用统一的半周期对称条件，在给定周期 `P` 下修正六维初值。输入必须是完整的
`[x, y, z, vx, vy, vz]`，其中对称面上的 `y`、`vx`、`vz` 应为零；返回修正后的
六维初值 `[x, 0, z, 0, vy, 0]`。

# Keywords
- `tol=1e-10`: 非线性残差收敛容差。
- `z_sign=0`: 设为 `±1` 时要求修正结果保持对应的 `z` 符号。
- `min_abs_z=0`: 非平面分支允许的最小 `|z|`。

# Returns
返回修正后的六维 `Vector{Float64}`。
"""
function orbit_shooting(
    dynamics::Function,
    state_guess::AbstractVector{<:Real},
    P::Real;
    tol::Real=DEFAULT_SHOOTING_TOLERANCE,
    z_sign::Int=0,
    min_abs_z::Real=0.0,
)
    return _shoot_symmetric_orbit(
        dynamics, state_guess, P; tol, z_sign, min_abs_z,
    )
end

"""
    halo_half_period_residual!(residual, u, parameters)

Half-period symmetry residual for a spatial CR3BP halo orbit. The unknowns
`u = [x, z, vy]` define the symmetry-plane crossing
`[x, 0, z, 0, vy, 0]`. At `P/2`, the corrected orbit must satisfy
`y = vx = vz = 0`.
"""
function halo_half_period_residual!(residual, u, parameters)
    return _symmetric_half_period_residual!(residual, u, parameters)
end

"""返回指定南/北分支和 L1/L2 平动点对应的内置 Halo 参考数据。"""
function _halo_reference(branch::Symbol, lp::Symbol)
    branch in (:northern, :southern) || throw(ArgumentError(
        "branch must be :northern or :southern; got $branch.",
    ))
    lp in (:L1, :L2) || throw(ArgumentError("lp must be :L1 or :L2; got $lp."))
    return ORBIT_REFERENCE_DATA[(:halo, branch, lp)]
end

"""
    generate_halo(; branch=:northern, lp=:L1, P=nothing, seed=nothing,
                    seed_period=nothing, continuation_step=0.005,
                    tol=1e-10, min_abs_z=1e-8) -> PeriodicOrbit

Generate a spatial CR3BP halo orbit around `lp`, using either the
`:northern` or `:southern` family branch. The reference initial condition is
an `x-z` symmetry-plane crossing of the form `[x, 0, z, 0, vy, 0]`.

If `P` is omitted, the period associated with the selected reference seed is
used. A different period is reached by short-step continuation with a secant
predictor and half-period differential correction. A custom seed is always a
full six-state `[x, y, z, vx, vy, vz]`; unless `seed_period` is supplied, it is
assumed to correspond to the target period. `min_abs_z` prevents a non-planar
branch from being silently accepted as the coexisting planar solution.

# Keywords
- `branch`: `:northern` or `:southern`.
- `lp`: `:L1` or `:L2`.
- `P=nothing`: target full period; the reference period is used when omitted.
- `seed`, `seed_period`: optional full six-state seed and its period.
- `continuation_step=0.005`: maximum period spacing between corrections.
- `tol=1e-10`: shooting tolerance.
- `min_abs_z=1e-8`: minimum accepted out-of-plane amplitude.

# Returns
Returns a corrected `PeriodicOrbit` propagated over one full period.
"""
function generate_halo(
    ;
    branch::Symbol=:northern,
    lp::Symbol=:L1,
    P::Union{Nothing,Real}=nothing,
    seed::Union{Nothing,AbstractVector{<:Real}}=nothing,
    seed_period::Union{Nothing,Real}=nothing,
    continuation_step::Real=DEFAULT_HALO_CONTINUATION_STEP,
    tol::Real=DEFAULT_SHOOTING_TOLERANCE,
    min_abs_z::Real=DEFAULT_MINIMUM_HALO_Z,
)
    reference = _halo_reference(branch, lp)
    target_period = isnothing(P) ? reference.P : Float64(P)
    initial_seed = isnothing(seed) ? reference.seed : seed
    initial_period = if isnothing(seed_period)
        isnothing(seed) ? reference.P : target_period
    else
        Float64(seed_period)
    end

    length(initial_seed) == 6 || throw(ArgumentError(
        "Periodic-orbit seeds must contain the full six-state [x, y, z, vx, vy, vz].",
    ))
    z_sign = Int(sign(initial_seed[3]))
    z_sign != 0 || throw(ArgumentError(
        "A halo seed must have nonzero z so its northern or southern branch can be preserved.",
    ))

    x0 = _continue_symmetric_orbit(
        cr3bp_eqm!,
        target_period,
        initial_seed,
        initial_period;
        continuation_step,
        tol,
        z_sign,
        min_abs_z,
    )
    return make_periodic_orbit(cr3bp_eqm!, x0, target_period; abstol=tol, reltol=tol)
end

"""
    generate_nrho_9_2(; branch=:southern, continuation_step=0.005,
                        tol=1e-10, min_abs_z=1e-8) -> PeriodicOrbit

Generate the selected northern or southern 9:2 near-rectilinear halo orbit
(NRHO) by continuing the corresponding L2 halo branch to
`P = 4pi / (9abs(ws))`, where `ws` is the Sun's angular rate in the Earth--Moon
rotating frame. The 9:2 designation refers to nine spacecraft revolutions in
two lunar synodic periods, not two sidereal Earth--Moon periods.

# Returns
Returns the selected L2 `PeriodicOrbit` at the 9:2 lunar-synodic resonant
period. Keywords control the branch, continuation step, tolerance, and minimum
out-of-plane amplitude as in `generate_halo`.
"""
function generate_nrho_9_2(
    ;
    branch::Symbol=:southern,
    continuation_step::Real=DEFAULT_HALO_CONTINUATION_STEP,
    tol::Real=DEFAULT_SHOOTING_TOLERANCE,
    min_abs_z::Real=DEFAULT_MINIMUM_HALO_Z,
)
    return generate_halo(
        ;
        branch,
        lp=:L2,
        P=NRHO_9_2_PERIOD,
        continuation_step,
        tol,
        min_abs_z,
    )
end

"""
    generate_DRO(; P=pi, seed=DRO_REFERENCE_GUESS, seed_period=pi,
                   continuation_step=0.05, tol=1e-10) -> PeriodicOrbit

生成具有目标无量纲周期 `P` 的平面远距逆行轨道（DRO）。默认 seed 是 2:1 DRO，
因此 `generate_DRO(P=pi)` 直接给出其打靶修正后的轨道；初值可由 `.x0` 取得。

当 `P` 不等于 `seed_period` 时，函数采用短步长、割线预测的周期参数延拓，并在
每一步执行固定周期打靶。整个过程不读写轨道数据文件，也不会建立完整轨道族。
若已有其他轨道族的近似解，可将其完整六维 seed 和对应 `seed_period` 传入；`z`
与 `x`、`vy` 一样参与每一步打靶及延拓。

# Keywords
- `P=pi`: 目标完整周期。
- `seed`, `seed_period`: 完整六维参考初值及其周期。
- `continuation_step=0.05`: 相邻周期修正的最大间距。
- `tol=1e-10`: 打靶容差。

# Returns
返回在一个完整周期上积分的 `PeriodicOrbit`。
"""
function generate_DRO(
    ;
    P::Real=DRO_REFERENCE_PERIOD,
    seed::AbstractVector{<:Real}=DRO_REFERENCE_GUESS,
    seed_period::Real=DRO_REFERENCE_PERIOD,
    continuation_step::Real=0.05,
    tol::Real=DEFAULT_SHOOTING_TOLERANCE,
)
    x0 = _continue_symmetric_orbit(
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
    plotPO!(ax, orbit; kwargs...) -> Lines

在 Makie 轴 `ax` 上绘制单条周期轨道的三维坐标曲线。额外关键字参数将传给
`lines!`，例如 `color` 或 `linewidth`。

# Returns
返回 Makie `Lines` 图元。
"""
function plotPO!(ax, orbit::AbstractPeriodicOrbit; kwargs...)
    times = range(orbit.sol.t[1], orbit.sol.t[end]; length=2000)
    states = orbit.sol(times).u
    orb = lines!(ax, getindex.(states, 1), getindex.(states, 2), getindex.(states, 3); kwargs...)
    return orb
end

"""
    plotPO!(ax, orbits::AbstractVector{<:AbstractPeriodicOrbit}; kwargs...) -> Lines

在同一 Makie 轴上批量绘制周期轨道。各轨道以 `NaN` 分隔，避免 Makie 将不同
轨道误连为一条线。

# Returns
返回包含全部轨道的单个 Makie `Lines` 图元。
"""
function plotPO!(ax, orbits::AbstractVector{<:AbstractPeriodicOrbit}; kwargs...)
    xs, ys, zs = Float64[], Float64[], Float64[]
    sizehint!(xs, 1001 * length(orbits))
    sizehint!(ys, 1001 * length(orbits))
    sizehint!(zs, 1001 * length(orbits))
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
