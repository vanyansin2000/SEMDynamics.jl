"""Numerical utilities for trajectory intersections and curve sampling."""
module Utils

export curve_intersect , uniform_sample_curve
using DifferentialEquations
using LinearAlgebra , Statistics

"""
    curve_intersect(trajectory1, trajectory2; samples=800) -> (states1, states2)

在两条二维 ODE 轨迹的线性采样折线上寻找交点，按第一条轨迹的时间排序后返回两侧
对应的插值状态；若没有交点则返回 `(nothing, nothing)`。

# Arguments
- `trajectory1`, `trajectory2`: 具有稠密插值的 `ODESolution`；每个状态的前两个
  分量作为平面坐标。

# Keywords
- `samples=800`: 每条轨迹用于折线近似的采样点数，至少为 2。

# Returns
存在交点时返回 `(states1, states2)`，两者按 `trajectory1` 的交点时间排序；无交点
时返回 `(nothing, nothing)`。该算法检测的是采样折线交点，精度受 `samples` 控制。
"""
function curve_intersect(
    Traj1::ODESolution,
    Traj2::ODESolution;
    samples::Integer=800,
)
    samples >= 2 || throw(ArgumentError("samples must be at least 2."))
    t1 = range(Traj1.t[1], Traj1.t[end]; length=samples)
    t2 = range(Traj2.t[1], Traj2.t[end]; length=samples)

    u1 = Traj1(t1).u
    u2 = Traj2(t2).u

    x1 = getindex.(u1, 1); y1 = getindex.(u1, 2)
    x2 = getindex.(u2, 1); y2 = getindex.(u2, 2)

    intersections = find_intersections(x1, y1, x2, y2, t1, t2, Traj1, Traj2)

    isempty(intersections) && return nothing, nothing

    # 按 Traj1(正向轨迹)时间升序排序，确保交点按发生先后排列 
    sort!(intersections, by = p -> p[1])

    t1_points = [p[1] for p in intersections]
    t2_points = [p[3] for p in intersections]

    return Traj1(t1_points), Traj2(t2_points)
end

"""
    find_intersections(x1, y1, x2, y2, t1, t2, trajectory1, trajectory2)
        -> Vector{Tuple}

枚举两组二维折线段的相交对，返回两条轨迹的交点时间和插值状态。
"""
function find_intersections(x1, y1, x2, y2, t1, t2, Traj1, Traj2)
    time_type = promote_type(eltype(t1), eltype(t2))
    state1_type = typeof(Traj1(first(t1)))
    state2_type = typeof(Traj2(first(t2)))
    intersections = Tuple{time_type,state1_type,time_type,state2_type}[]
    n1, n2 = length(x1)-1, length(x2)-1
    
    for i in 1:n1, j in 1:n2
        # 边界框预检查
        (min_x1, max_x1) = extrema((x1[i], x1[i+1]))
        (min_y1, max_y1) = extrema((y1[i], y1[i+1]))
        (min_x2, max_x2) = extrema((x2[j], x2[j+1]))
        (min_y2, max_y2) = extrema((y2[j], y2[j+1]))
        
        # 如果边界框不重叠，跳过详细计算[1](@ref)
        max_x1 < min_x2 && continue
        min_x1 > max_x2 && continue
        max_y1 < min_y2 && continue
        min_y1 > max_y2 && continue
        
        # 计算线段交点
        intersect_point, t, u = segment_intersection(
            x1[i], y1[i], x1[i+1], y1[i+1],
            x2[j], y2[j], x2[j+1], y2[j+1])
        
        if intersect_point !== nothing && 0 <= t <= 1 && 0 <= u <= 1
            # 计算交点对应的时间
            t1_intersect = t1[i] + t * (t1[i+1] - t1[i])
            t2_intersect = t2[j] + u * (t2[j+1] - t2[j])
            
            # 获取交点状态
            state1 = Traj1(t1_intersect)
            state2 = Traj2(t2_intersect)
            
            push!(intersections, (t1_intersect, state1, t2_intersect, state2))
        end
    end
    
    return intersections
end

"""
    segment_intersection(x1, y1, x2, y2, x3, y3, x4, y4)
        -> (point, t, u)

返回两条二维线段的交点及各自线性参数；平行或不相交时 `point` 为 `nothing`。
"""
function segment_intersection(x1, y1, x2, y2, x3, y3, x4, y4)
    # 向量叉积法计算线段交点[1](@ref)
    denom = (x1 - x2) * (y3 - y4) - (y1 - y2) * (x3 - x4)
    
    abs(denom) < 1e-12 && return nothing, 0.0, 0.0
    
    t = ((x1 - x3) * (y3 - y4) - (y1 - y3) * (x3 - x4)) / denom
    u = -((x1 - x2) * (y1 - y3) - (y1 - y2) * (x1 - x3)) / denom
    
    if 0 <= t <= 1 && 0 <= u <= 1
        x = x1 + t * (x2 - x1)
        y = y1 + t * (y2 - y1)
        return (x, y), t, u
    end
    
    return nothing, 0.0, 0.0
end


"""
    uniform_sample_curve(points, n_samples; start_idx=nothing) -> Vector{Int}

先用最近邻路径排序点集，再按该路径的累计弧长近似均匀选点。

# Arguments
- `points`: 每行一个点的 `N × D` 实数矩阵。
- `n_samples`: 目标采样数量，必须位于 `1:N`。

# Keywords
- `start_idx=nothing`: 路径起点的原始行索引；省略时选择离点集质心最远的点。

# Returns
返回采样点在原矩阵中的行索引。相邻目标弧长映射到同一点时会去重，因此极端情况下
返回数量可能少于 `n_samples`。
"""
function uniform_sample_curve(points::AbstractMatrix, n_samples::Int; start_idx=nothing)
    N = size(points, 1)
    1 <= n_samples <= N || throw(ArgumentError("n_samples must be between 1 and $N."))
    size(points, 2) > 0 || throw(ArgumentError("points must have at least one column."))
    if !isnothing(start_idx)
        1 <= start_idx <= N || throw(BoundsError(points, (start_idx, :)))
    end
    
    # === 1. 点集排序 (Nearest Neighbor Path) ===
    # 注意：对于复杂曲线（如螺旋或U型），简单的最近邻可能会出错。
    # 如果数据很乱，可能需要更高级的流形学习算法。
    
    indices_ordered = Int[]
    visited = falses(N)
    
    # 自动寻找端点：通常是距离几何中心最远的点之一，或者是主要方向上的极值点
    if isnothing(start_idx)
        # 简单策略：取x最小的点作为起点（假设曲线大概是横向的）
        # 或者：取距离重心最远的点
        center = vec(mean(points, dims=1))
        dists_to_center = Vector{Float64}(undef, N)
        @inbounds for i in 1:N
            distance² = 0.0
            for column in axes(points, 2)
                distance² += abs2(points[i, column] - center[column])
            end
            dists_to_center[i] = distance²
        end
        current_idx = argmax(dists_to_center) # 选一个“端点”
    else
        current_idx = start_idx
    end
    
    push!(indices_ordered, current_idx)
    visited[current_idx] = true
    
    for _ in 1:N-1
        # 寻找最近的未访问点
        min_dist = Inf
        next_idx = -1
        
        # 暴力搜索 (O(N^2))，对于 N < 5000 足够快
        # 如果 N 很大，建议使用 KDTree
        @inbounds for i in 1:N
            if !visited[i]
                distance² = 0.0
                for column in axes(points, 2)
                    distance² += abs2(points[i, column] - points[current_idx, column])
                end
                if distance² < min_dist
                    min_dist = distance²
                    next_idx = i
                end
            end
        end
        
        if next_idx != -1
            push!(indices_ordered, next_idx)
            visited[next_idx] = true
            current_idx = next_idx
        else
            break
        end
    end
    
    # === 2. 计算累积弧长 ===
    cum_dist = zeros(Float64, length(indices_ordered))
    @inbounds for path_index in 2:length(indices_ordered)
        previous_index = indices_ordered[path_index - 1]
        current_index = indices_ordered[path_index]
        distance² = 0.0
        for column in axes(points, 2)
            distance² += abs2(points[current_index, column] - points[previous_index, column])
        end
        cum_dist[path_index] = cum_dist[path_index - 1] + sqrt(distance²)
    end
    total_length = cum_dist[end]
    
    # === 3. 均匀采样 ===
    # 目标弧长位置
    target_dists = range(0, total_length, length=n_samples)
    
    selected_indices_in_ordered = Int[]
    
    # 对每个目标位置，寻找原数据中弧长最接近的那个点
    # 这里使用简单的搜索，因为 cum_dist 是单调的
    for target in target_dists
        upper = searchsortedfirst(cum_dist, target)
        selected = if upper <= 1
            1
        elseif upper > length(cum_dist)
            length(cum_dist)
        elseif target - cum_dist[upper - 1] <= cum_dist[upper] - target
            upper - 1
        else
            upper
        end
        push!(selected_indices_in_ordered, selected)
    end
    # 映射回原始索引
    final_indices = indices_ordered[selected_indices_in_ordered]
    
    # 去重（可选）：如果采样过密，同一个点可能被选中两次
    unique!(final_indices)
    
    return final_indices
end

end

