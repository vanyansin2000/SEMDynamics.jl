
export cb_p2collision, cb_p1collision, cb_escape, cb_enter
export cb_energy_zero, cb_xy_zero_cross, cb_y_zero_cross
export cb_perilune, cb_apolune, cb_apse_p1, cb_apse_p2
export Event , dynamic_events
export ODESolveFailedError

# 捕获条件
# 初始化记录存储（可在积分参数中存储）
"""
    Event(code, time, state)

记录一个检测到的动力学事件。

# Fields
- `code::Symbol`: 事件类别，例如 `:enter`、`:perilune` 或 `:p2collision`。
- `time::Number`: 事件发生的无量纲积分时间。
- `state::AbstractVector`: 事件根处的状态快照。
"""
struct Event
    code ::Symbol
    time  ::Number
    state ::AbstractVector{<:Number}
end

"""
    dynamic_events() -> Vector{Event}

创建空事件容器，供本模块的 callback 在积分期间写入 `Event`。
"""
dynamic_events() = Vector{Event}(undef , 0);

# 积分失败错误
"""
    ODESolveFailedError(msg)

表示 ODE 积分未能成功结束的异常；`msg` 保存面向用户的错误说明。
"""
struct ODESolveFailedError <: Exception
    msg ::String
end

Base.showerror(io::IO, error::ODESolveFailedError) = print(io, error.msg)

"""返回状态的物理空间维数，并兼容包含 STM 的 20/42 维扩展状态。"""
@inline function _state_dimension(u::AbstractVector)
    length(u) in (4, 20) && return 2
    length(u) in (6, 42) && return 3
    throw(ArgumentError("事件回调仅支持长度为 4、6、20 或 42 的状态，得到 $(length(u))"))
end

"""
    state_from_r1_r2(u, μ) -> (r1, r2)

计算状态到两个主天体的距离。支持平面状态（长度 4 或 20）和空间状态（长度 6 或
42）；STM 扩展部分会被自动忽略。
"""
@inline function state_from_r1_r2(u::AbstractVector{<:Real} , μ)
    dimension = _state_dimension(u)
    x, y = u[1], u[2]
    z = dimension == 3 ? u[3] : zero(x)
    return sqrt((x + μ)^2 + y^2 + z^2), sqrt((x + μ - 1)^2 + y^2 + z^2)
end


"""返回相对第二主天体的径向速度指标、二阶导数和距离。"""
@inline function apse_p2(integrator)
    return _apse(integrator, :p2)
end

"""返回相对第一主天体的径向速度指标、二阶导数和距离。"""
@inline function apse_p1(integrator)
    return _apse(integrator, :p1)
end

"""
    _apse(integrator, primary) -> (r_dot_times_r, radial_second_derivative, distance)

计算相对指定主天体的拱点判据。状态和加速度索引由 `_state_dimension` 自动选择，
故同一近/远拱点回调可用于平面、空间及其 STM 扩展积分。
"""
@inline function _apse(integrator, primary::Symbol)
    u = integrator.u
    t = integrator.t
    μ = integrator.p.EMRot.μ
    du = similar(u)
    integrator.f(du, u, integrator.p, t)
    if _state_dimension(u) == 2
        rx = primary === :p1 ? u[1] + μ : u[1] - (1 - μ)
        ry, vx, vy, ax, ay = u[2], u[3], u[4], du[3], du[4]
        distance = hypot(rx, ry)
        radial_velocity = rx * vx + ry * vy
        radial_acceleration = vx^2 + vy^2 + rx * ax + ry * ay
        return radial_velocity, radial_acceleration, distance
    end

    rx = primary === :p1 ? u[1] + μ : u[1] - (1 - μ)
    ry, rz = u[2], u[3]
    vx, vy, vz = u[4], u[5], u[6]
    ax, ay, az = du[4], du[5], du[6]
    distance = sqrt(rx^2 + ry^2 + rz^2)
    radial_velocity = rx * vx + ry * vy + rz * vz
    radial_acceleration = vx^2 + vy^2 + vz^2 + rx * ax + ry * ay + rz * az
    return radial_velocity, radial_acceleration, distance
end


##处理积分事件
##————————————————————————————————————————————————————————————————————————————————————
#停止积分事件 1碰撞
# 撞击事件条件函数
"""第二主天体表面的有符号距离；负值表示位于碰撞区域内。"""
function p2_collision_condition(u, t, integrator , scale)
    p = integrator.p  # 获取系统参数
    μ = p.EMRot.μ

    _ , r2 = state_from_r1_r2(u , μ)
    return r2 - scale * p.EMRot.r_p2
end

"""第一主天体表面的有符号距离；负值表示位于碰撞区域内。"""
function p1_collision_condition(u, t, integrator , scale )
    p = integrator.p  # 获取系统参数
    μ = p.EMRot.μ
    r1 , _ = state_from_r1_r2(u, μ)
    return r1 - scale * p.EMRot.r_p1
end

# 撞击事件影响函数：终止积分
"""记录第二主天体碰撞事件，并按需终止积分。"""
function affect_p2collision!(integrator , event_storage::Vector{Event}, terminate) 
    code = :p2collision
    time = integrator.t
    state = integrator.u 
    event = Event(code, time, copy(state))
    push!(event_storage , event)

    if terminate
        terminate!(integrator)
    end
    return nothing
end

"""记录第一主天体碰撞事件，并按需终止积分。"""
function affect_p1collision!(integrator , event_storage::Vector{Event}, terminate) 
    code = :p1collision
    time = integrator.t
    state = integrator.u 
    event = Event(code, time, copy(state))
    push!(event_storage , event)

    if terminate
        terminate!(integrator)
    end
    return nothing
end


"""
    cb_p2collision(event_storage; terminate=true, scale=1.0) -> ContinuousCallback

检测从外向内穿越第二主天体表面的事件，并记录 `:p2collision`。

# Keywords
- `terminate=true`: 记录后是否终止积分。
- `scale=1.0`: 相对默认无量纲半径 `p.EMRot.r_p2` 的缩放系数。

# Returns
返回仅响应负向穿越的 `ContinuousCallback`。
"""
function cb_p2collision(event_storage::Vector{Event} ;terminate = true , scale = 1.0)
    condition_wrapper(u, t, integrator) = p2_collision_condition(u, t, integrator, scale)
    affect_wrapper = (integrator) -> affect_p2collision!(integrator, event_storage, terminate)
    # Collision entry is a positive-to-negative crossing.  Use the negative
    # crossing callback and locate the surface before the singular interior.
    return ContinuousCallback(condition_wrapper, nothing, affect_wrapper; rootfind=true)
end

"""
    cb_p1collision(event_storage; terminate=true, scale=1.0) -> ContinuousCallback

检测从外向内穿越第一主天体表面的事件并记录 `:p1collision`。关键字含义与
`cb_p2collision` 相同。
"""
function cb_p1collision(event_storage::Vector{Event} ;terminate = true, scale = 1.0)
    condition_wrapper(u, t, integrator) = p1_collision_condition(u, t, integrator, scale)
    affect_wrapper = (integrator) -> affect_p1collision!(integrator, event_storage, terminate)
    return ContinuousCallback(condition_wrapper, nothing, affect_wrapper; rootfind=true)
end



# 
##————————————————————————————————————————————————————————————————————————————————————
#停止积分事件 2 进入或逃逸第二主天体作用球
"""第二主天体作用球的有符号距离；球内为负，球外为正。"""
function p2_sphere_condition(u, t, integrator, scale)
    p = integrator.p
    μ = p.EMRot.μ
    _ , r2 = state_from_r1_r2(u , μ)
    return r2 - scale * p.EMRot.d2lim
end

"""记录第二主天体作用球的进入或逃逸事件，并按需终止积分。"""
function affect_p2_sphere!(
    integrator,
    event_storage::Vector{Event},
    code::Symbol,
    terminate::Bool,
)
    push!(event_storage, Event(code, integrator.t, copy(integrator.u)))
    terminate && terminate!(integrator)
    return nothing
end

"""
    cb_escape(event_storage; scale=1.0, terminate=true, root_find=true)

检测轨迹从内向外穿越第二主天体作用球的事件。球半径为
`scale * integrator.p.EMRot.d2lim`，事件代码为 `:escape`。

# Returns
返回仅响应正向穿越的 `ContinuousCallback`；`terminate` 控制是否终止积分，
`root_find` 控制是否定位精确事件根。
"""
function cb_escape(
    event_storage::Vector{Event};
    scale::Real=1.0,
    terminate::Bool=true,
    root_find::Bool=true,
)
    scale > 0 || throw(ArgumentError("scale must be positive."))
    condition_wrapper(u, t, integrator) = p2_sphere_condition(u, t, integrator, scale)
    affect_wrapper = integrator -> affect_p2_sphere!(
        integrator, event_storage, :escape, terminate,
    )
    return ContinuousCallback(
        condition_wrapper,
        affect_wrapper,
        nothing;
        rootfind=root_find,
        save_positions=(false, true),
    )
end

"""
    cb_enter(event_storage; scale=1.0, terminate=true, root_find=true)

检测轨迹从外向内穿越第二主天体作用球的事件，是 `cb_escape` 的反向事件。
球半径为 `scale * integrator.p.EMRot.d2lim`，事件代码为 `:enter`。

# Returns
返回仅响应负向穿越的 `ContinuousCallback`。
"""
function cb_enter(
    event_storage::Vector{Event};
    scale::Real=1.0,
    terminate::Bool=true,
    root_find::Bool=true,
)
    scale > 0 || throw(ArgumentError("scale must be positive."))
    condition_wrapper(u, t, integrator) = p2_sphere_condition(u, t, integrator, scale)
    affect_wrapper = integrator -> affect_p2_sphere!(
        integrator, event_storage, :enter, terminate,
    )
    return ContinuousCallback(
        condition_wrapper,
        nothing,
        affect_wrapper;
        rootfind=root_find,
        save_positions=(false, true),
    )
end
##————————————————————————————————————————————————————————————————————————————————————
# ε2 零交叉事件条件函数
"""返回相对第二主天体能量 `ε₂`，用于检测其零交叉。"""
function energy_zero_condition(u, t, integrator)
    μ = integrator.p.EMRot.μ
    ε2 = compute_ε2(u , t , μ)
    return ε2  # 当 ε2 接近零时触发事件
end

# ε2 零交叉事件影响函数：记录时间和状态
"""记录能量零交叉事件，并按需终止积分。"""
function affect_energy_zero!(integrator, dynamic_events::Vector{Event}, terminate)
    
    code = :energy_zero
    time = integrator.t
    state = integrator.u 
    event = Event(code, time, copy(state))
    push!(dynamic_events , event)
    
    if terminate
        terminate!(integrator)
    end
    return nothing
end

# 创建 ε2 零交叉事件回调（不终止积分）
# cb_energy_zero = ContinuousCallback(energy_zero_condition, affect_energy_zero!)

"""
    cb_energy_zero(event_storage; terminate=false) -> ContinuousCallback

检测相对第二主天体的两体能量 `ε₂=0` 穿越，记录 `:energy_zero`，并按
`terminate` 决定是否结束积分。
"""
function cb_energy_zero(event_storage::Vector{Event}; terminate = false)
    affect_wrapper = (integrator) -> affect_energy_zero!(integrator, event_storage, terminate)
    return ContinuousCallback(energy_zero_condition, affect_wrapper;save_positions = (true, false))
end


# x，y零交叉事件条件函数
"""返回 `(x + 1 - μ)y`，用于检测坐标轴组合零交叉。"""
function xy_zero_cross_condition(u, t, integrator)
    μ = integrator.p.EMRot.μ
    return (u[1] + 1-μ) * u[2]  # 当[x-(1-mu)]*y=0触发事件
end

# x，y零交叉事件影响函数
"""`xy` 零交叉的无状态回调动作。"""
affect_xy_zero_cross!(integrator) = nothing


# 创建  x，y零交叉事件回调（不终止积分）
"""
    cb_xy_zero_cross() -> ContinuousCallback

创建 `(x + 1 - μ)y = 0` 的无状态连续回调。该回调不执行根定位，也不记录事件。
"""
cb_xy_zero_cross() = ContinuousCallback(xy_zero_cross_condition,
 affect_xy_zero_cross!; rootfind=false, save_positions = (false, true))


"""返回 `y`，用于检测 xz 对称平面的穿越。"""
function y_zero_cross_condition(u, t, integrator)
    return u[2]  # 当[x-(1-mu)]*y=0触发事件
end

# x，y零交叉事件影响函数
"""记录 y=0 截面穿越事件，并按需终止积分。"""
function affect_y_zero_cross!(integrator, dynamic_events::Vector{Event} , terminate)
    
    code = :y_zero_cross
    time = integrator.t
    state = integrator.u 
    event = Event(code, time, copy(state))
    push!(dynamic_events , event)
    
    if terminate
        terminate!(integrator)
    end
    return nothing
end


"""
    cb_y_zero_cross(event_storage; terminate=false) -> ContinuousCallback

检测 `y=0`（xz 对称面）穿越，记录 `:y_zero_cross`，并按需终止积分。
"""
function cb_y_zero_cross(event_storage::Vector{Event}; terminate = false)
    affect_wrapper = (integrator) -> affect_y_zero_cross!(integrator, event_storage, terminate)
    return ContinuousCallback(y_zero_cross_condition, affect_wrapper;  rootfind = true, save_positions = (false, true))
end

##————————————————————————————————————————————————————————————————————————————————————
# condition 只负责定位 r·v = 0 的过零点（与方向无关）
"""返回相对第二主天体的径向速度内积 `r⋅v`。"""
@inline function rv_p2_condition(u, t, integrator)
    μ = integrator.p.EMRot.μ
    rx = u[1] - (1 - μ)
    if _state_dimension(u) == 2
        return rx * u[3] + u[2] * u[4]
    end
    return rx * u[4] + u[2] * u[5] + u[3] * u[6]
end

"""返回相对第一主天体的径向速度内积 `r⋅v`。"""
@inline function rv_p1_condition(u, t, integrator)
    μ = integrator.p.EMRot.μ
    rx = u[1] + μ
    if _state_dimension(u) == 2
        return rx * u[3] + u[2] * u[4]
    end
    return rx * u[4] + u[2] * u[5] + u[3] * u[6]
end

##————————————————————————————————————————————————————————————————————————————————————
"""按径向二阶导数符号分类并记录指定主天体的近/远拱点。"""
function _affect_apse!(
    integrator,
    event_storage::Vector{Event};
    primary::Symbol,
    near_code::Symbol,
    far_code::Symbol,
    terminate_near::Bool=false,
    terminate_far::Bool=false,
    near_scale=(0.0, Inf),
    far_scale=(0.0, Inf),
)
    _, d2, d = _apse(integrator, primary)

    if d2 > 0
        lb, ub = near_scale
        if lb < d < ub
            push!(event_storage, Event(near_code, integrator.t, copy(integrator.u)))
            terminate_near && terminate!(integrator)
        end
    elseif d2 < 0
        lb, ub = far_scale
        if lb < d < ub
            push!(event_storage, Event(far_code, integrator.t, copy(integrator.u)))
            terminate_far && terminate!(integrator)
        end
    end

    return nothing
end

"""构造由 `r⋅v=0` 定位、同时检测近点和远点的通用连续回调。"""
function _cb_apse(
    event_storage::Vector{Event};
    primary::Symbol,
    near_code::Symbol,
    far_code::Symbol,
    terminate_near::Bool=false,
    terminate_far::Bool=false,
    near_scale=(0.0, Inf),
    far_scale=(0.0, Inf),
    root_find::Bool=true,
)
    affect_wrapper = integrator -> _affect_apse!(
        integrator,
        event_storage;
        primary,
        near_code,
        far_code,
        terminate_near,
        terminate_far,
        near_scale,
        far_scale,
    )
    condition = primary === :p1 ? rv_p1_condition : rv_p2_condition

    return ContinuousCallback(
        condition,
        affect_wrapper,
        affect_wrapper;
        rootfind = root_find,
        save_positions = (false, true),
    )
end

"""
    cb_apse_p1(event_storage; terminate_perigee=false, terminate_apogee=false,
               perigee_scale=(0, Inf), apogee_scale=(0, Inf), root_find=true)

检测相对第一主天体的近地点与远地点，分别记录 `:perigee` 和 `:apogee`。

# Keywords
- `terminate_perigee`, `terminate_apogee`: 是否在对应事件处终止。
- `perigee_scale`, `apogee_scale`: 允许记录的距离开区间 `(lower, upper)`。
- `root_find=true`: 是否由回调定位 `r⋅v=0` 的根。

# Returns
返回同时响应正、负方向积分的 `ContinuousCallback`。
"""
function cb_apse_p1(
    event_storage::Vector{Event};
    terminate_perigee::Bool=false,
    terminate_apogee::Bool=false,
    perigee_scale=(0.0, Inf),
    apogee_scale=(0.0, Inf),
    root_find::Bool=true,
)
    return _cb_apse(
        event_storage;
        primary=:p1,
        near_code=:perigee,
        far_code=:apogee,
        terminate_near=terminate_perigee,
        terminate_far=terminate_apogee,
        near_scale=perigee_scale,
        far_scale=apogee_scale,
        root_find,
    )
end

"""
    cb_apse_p2(event_storage; terminate_perilune=false, terminate_apolune=false,
               perilune_scale=(0, Inf), apolune_scale=(0, Inf), root_find=true)

检测相对第二主天体的近月点与远月点，分别记录 `:perilune` 和 `:apolune`。
接口与 `cb_apse_p1` 对称。

# Returns
返回 `ContinuousCallback`；各距离区间和终止关键字分别控制近月点与远月点。
"""
function cb_apse_p2(
    event_storage::Vector{Event};
    terminate_perilune::Bool=false,
    terminate_apolune::Bool=false,
    perilune_scale=(0.0, Inf),
    apolune_scale=(0.0, Inf),
    root_find::Bool=true,
)
    return _cb_apse(
        event_storage;
        primary=:p2,
        near_code=:perilune,
        far_code=:apolune,
        terminate_near=terminate_perilune,
        terminate_far=terminate_apolune,
        near_scale=perilune_scale,
        far_scale=apolune_scale,
        root_find,
    )
end

"""
    cb_perilune(event_storage; terminate=false, scale=(0, Inf), root_find=true)
        -> ContinuousCallback

创建仅记录 `:perilune` 的兼容回调；新代码可优先使用 `cb_apse_p2`。
"""
function cb_perilune(
    event_storage::Vector{Event};
    terminate::Bool=false,
    scale=(0.0, Inf),
    root_find::Bool=true,
)
    return cb_apse_p2(
        event_storage;
        terminate_perilune=terminate,
        perilune_scale=scale,
        apolune_scale=(Inf, Inf),
        root_find,
    )
end

"""
    cb_apolune(event_storage; terminate=false, scale=(0, Inf), root_find=true)
        -> ContinuousCallback

创建仅记录 `:apolune` 的回调；也可使用完整的 `cb_apse_p2` 接口。
"""
function cb_apolune(
    event_storage::Vector{Event};
    terminate::Bool=false,
    scale=(0.0, Inf),
    root_find::Bool=true,
)
    return cb_apse_p2(
        event_storage;
        terminate_apolune=terminate,
        perilune_scale=(Inf, Inf),
        apolune_scale=scale,
        root_find,
    )
end

##————————————————————————————————————————————————————————————————————————————————————
# # ε2 旋转系中，一圈事件条件函数
# function revolution_condition(u, t, integrator)
#     μ = integrator.p.μ
#     Ψ = polar_angle_p2(μ , u)
#     Ψ - 
#     return ε2  # 当 ε2 接近零时触发事件
# end

# # ε2 零交叉事件影响函数：记录时间和状态
# function affect_energy_zero!(integrator, dynamic_events::Vector{Event} ; terminate = false)
    
#     code = :energy_zero
#     time = integrator.t
#     state = integrator.u 
#     event = Event(code, time, copy(state))
#     push!(dynamic_events , event)
    
#     if terminate
#         terminate!(integrator)
#     end
#     return nothing
# end

# # 创建 ε2 零交叉事件回调（不终止积分）
# # cb_energy_zero = ContinuousCallback(energy_zero_condition, affect_energy_zero!)

# "ε2 零交叉事件回调"
# function cb_energy_zero(event_storage::Vector{Event})
#     affect_wrapper = (integrator) -> affect_energy_zero!(integrator, event_storage)
#     return ContinuousCallback(revolution_condition, affect_wrapper;save_positions = (true, false))
# end
