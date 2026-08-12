
export cb_p2collision, cb_p1collision, cb_escape , cb_energy_zero , cb_xy_zero_cross , cb_y_zero_cross
export cb_perilune , cb_apse_p1
export Event , dynamic_events
export ODESolveFailedError

# 捕获条件
# 初始化记录存储（可在积分参数中存储）
"""记录一个检测到的动力学事件：事件代码、发生时间和状态快照。"""
struct Event
    code ::Symbol
    time  ::Number
    state ::AbstractVector{<:Number}
end

"""创建用于收集积分事件的空 `Vector{Event}`。"""
dynamic_events() = Vector{Event}(undef , 0);

# 积分失败错误
"""表示 ODE 积分未能成功结束的异常类型。"""
struct ODESolveFailedError <: Exception
    msg ::String
end

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
"""判断状态是否进入第二主天体半径乘以 `scale` 的碰撞区域。"""
function p2_collision_condition(u, t, integrator , scale)
    p = integrator.p  # 获取系统参数
    μ = p.EMRot.μ

    _ , r2 = state_from_r1_r2(u , μ)
    return r2 -  scale*p.EMRot.r_p2 < 0  # 当返回true时触发事件
end

"""判断状态是否进入第一主天体半径乘以 `scale` 的碰撞区域。"""
function p1_collision_condition(u, t, integrator , scale )
    p = integrator.p  # 获取系统参数
    μ = p.EMRot.μ
    r1 , _ = state_from_r1_r2(u, μ)
    return r1 -  scale*p.EMRot.r_p1 < 0  # 当返回true时触发事件
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


"创建撞击事件回调"
# cb_collision = DiscreteCallback(collision_condition, affect_collision!)
# 创建第二主天体碰撞的连续回调，并将事件写入 `event_storage`。
function cb_p2collision(event_storage::Vector{Event} ;terminate = true , scale = 1.0)
    condition_wrapper(u, t, integrator) = p2_collision_condition(u, t, integrator, scale)
    affect_wrapper = (integrator) -> affect_p2collision!(integrator, event_storage, terminate)
    return ContinuousCallback(condition_wrapper, affect_wrapper , rootfind = false)
end

"""创建第一主天体碰撞的连续回调，并将事件写入 `event_storage`。"""
function cb_p1collision(event_storage::Vector{Event} ;terminate = true, scale = 1.0)
    condition_wrapper(u, t, integrator) = p1_collision_condition(u, t, integrator, scale)
    affect_wrapper = (integrator) -> affect_p1collision!(integrator, event_storage, terminate)
    return ContinuousCallback(condition_wrapper, affect_wrapper , rootfind = false)
end



# 
##————————————————————————————————————————————————————————————————————————————————————
#停止积分事件 2 逃逸地月空间
# 逃逸事件条件函数
"""判断轨迹相对第二主天体是否超出逃逸距离阈值。"""
function escape_condition(u, t, integrator, scale)
    p = integrator.p
    μ = p.EMRot.μ
    _ , r2 = state_from_r1_r2(u , μ)
    return r2  - scale * p.EMRot.d2lim > 0   # 当返回true时触发事件
end

# 撞击事件影响函数：终止积分
"""记录逃逸事件并立即终止积分。"""
function affect_escape!(integrator , event_storage::Vector{Event}) 
    code = :escape
    time = integrator.t
    state = integrator.u 
    event = Event(code, time, copy(state))
    push!(event_storage , event)
    terminate!(integrator)
    return nothing
end

"逃逸事件回调"
# cb_escape = DiscreteCallback(escape_condition, affect_escape!)
# 创建基于离散条件的逃逸事件回调。
function cb_escape(event_storage::Vector{Event}; scale = 1.0 )
    condition_wrapper(u, t, integrator) = escape_condition(u, t, integrator, scale)
    affect_wrapper = (integrator) -> affect_escape!(integrator, event_storage)
    return DiscreteCallback(condition_wrapper, affect_wrapper)
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

"ε2 零交叉事件回调"
# 创建相对第二主天体能量零交叉的连续回调。
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
"""创建不记录状态的 `xy` 零交叉连续回调。"""
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


"""创建记录 y=0 截面穿越的连续回调。"""
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
# 1. 近月点（近 P2 点）：r2·v2 = 0，d̈ > 0（极小值），lb < r2 < ub 才记录
"""在第二主天体近拱点且距离位于 `scale` 区间内时记录事件。"""
function affect_perilune!(integrator, event_storage::Vector{Event}, terminate, scale)
    rv, d2, d = apse_p2(integrator)
    lb, ub = scale
    if d2 > 0 && lb < d < ub        # 二阶导 > 0 ⇒ 近点（与积分方向无关）
        push!(event_storage, Event(:perilune, integrator.t, copy(integrator.u)))
        terminate && terminate!(integrator)
    end
    return nothing
end

"近月点（近 P2 点）：r·v=0 且 d̈>0，仅当 lb < r2 < ub 记录"
# 创建第二主天体近拱点（perilune）检测回调。
function cb_perilune(event_storage::Vector{Event}; terminate = false, scale = (0.0, Inf) , root_find = true)
    affect_wrapper = (integrator) -> affect_perilune!(integrator, event_storage, terminate, scale)
    # 两个方向都接同一个回调，靠二阶导筛选近点，避免正/逆向积分误判
    return ContinuousCallback(rv_p2_condition, affect_wrapper, affect_wrapper;
                              rootfind = root_find, save_positions = (false, true))
end

##————————————————————————————————————————————————————————————————————————————————————
function affect_apse_p1!(
    integrator,
    event_storage::Vector{Event};
    terminate_perigee = false,
    terminate_apogee = false,
    perigee_scale = (0.0, Inf),
    apogee_scale = (0.0, Inf),
)
    _, d2, d = apse_p1(integrator)

    if d2 > 0
        lb, ub = perigee_scale
        if lb < d < ub
            push!(event_storage, Event(:perigee, integrator.t, copy(integrator.u)))
            terminate_perigee && terminate!(integrator)
        end
    elseif d2 < 0
        lb, ub = apogee_scale
        if lb < d < ub
            push!(event_storage, Event(:apogee, integrator.t, copy(integrator.u)))
            terminate_apogee && terminate!(integrator)
        end
    end

    return nothing
end

function cb_apse_p1(
    event_storage::Vector{Event};
    terminate_perigee = false,
    terminate_apogee = false,
    perigee_scale = (0.0, Inf),
    apogee_scale = (0.0, Inf),
    root_find = true,
)
    affect_wrapper = integrator -> affect_apse_p1!(
        integrator,
        event_storage;
        terminate_perigee = terminate_perigee,
        terminate_apogee = terminate_apogee,
        perigee_scale = perigee_scale,
        apogee_scale = apogee_scale,
    )

    return ContinuousCallback(
        rv_p1_condition,
        affect_wrapper,
        affect_wrapper;
        rootfind = root_find,
        save_positions = (false, true),
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
