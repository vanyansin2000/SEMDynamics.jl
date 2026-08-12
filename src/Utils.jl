module Utils

export curve_intersect , uniform_sample_curve
using DifferentialEquations
using LinearAlgebra , Statistics

"""
    curve_intersect(trajectory1, trajectory2)

在两条二维 ODE 轨迹的线性采样折线上寻找交点，按第一条轨迹的时间排序后返回两侧
对应的插值状态；若没有交点则返回 `(nothing, nothing)`。
"""
function curve_intersect(Traj1::ODESolution, Traj2::ODESolution)
    t1 = range(Traj1.t[1], Traj1.t[end], 800)
    t2 = range(Traj2.t[1], Traj2.t[end], 800)

    u1 = Traj1(t1).u
    u2 = Traj2(t2).u

    x1 = [u[1] for u in u1]; y1 = [u[2] for u in u1]
    x2 = [u[1] for u in u2]; y2 = [u[2] for u in u2]

    intersections = find_intersections(x1, y1, x2, y2, t1, t2, Traj1, Traj2)

    isempty(intersections) && return nothing, nothing

    # 按 Traj1(正向轨迹)时间升序排序，确保交点按发生先后排列 
    sort!(intersections, by = p -> p[1])

    t1_points = [p[1] for p in intersections]
    t2_points = [p[3] for p in intersections]

    return Traj1(t1_points), Traj2(t2_points)
end

"""枚举两组二维折线段的相交对，并计算每个交点在两条 ODE 解上的插值状态。"""
function find_intersections(x1, y1, x2, y2, t1, t2, Traj1, Traj2)
    intersections = []
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

"""返回两条二维线段的交点和各自线性参数；平行或不相交时返回 `nothing`。"""
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
    uniform_sample_curve(points, n_samples; start_idx=nothing)

在分布呈曲线的点集中均匀采样。
- points: N x 2 的矩阵或 Point2f 数组
- n_samples: 需要采样的点数
- start_idx: (可选) 指定曲线起始点的索引

返回: 采样出的点的原始索引列表 Vector{Int}
"""
function uniform_sample_curve(points::AbstractMatrix, n_samples::Int; start_idx=nothing)
    N = size(points, 1)
    @assert n_samples <= N "采样数不能大于点总数"
    
    # === 1. 点集排序 (Nearest Neighbor Path) ===
    # 注意：对于复杂曲线（如螺旋或U型），简单的最近邻可能会出错。
    # 如果数据很乱，可能需要更高级的流形学习算法。
    
    indices_ordered = Int[]
    visited = falses(N)
    
    # 自动寻找端点：通常是距离几何中心最远的点之一，或者是主要方向上的极值点
    if isnothing(start_idx)
        # 简单策略：取x最小的点作为起点（假设曲线大概是横向的）
        # 或者：取距离重心最远的点
        center = mean(points, dims=1)
        dists_to_center = vec(sum(abs2, points .- center, dims=2))
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
        p_curr = points[current_idx, :]
        
        for i in 1:N
            if !visited[i]
                d = norm(points[i, :] - p_curr)
                if d < min_dist
                    min_dist = d
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
    ordered_points = points[indices_ordered, :]
    # 计算相邻点距离
    diffs = ordered_points[2:end, :] .- ordered_points[1:end-1, :]
    dists = sqrt.(sum(abs2, diffs, dims=2))
    
    # 累积弧长 (从0开始)
    cum_dist = [0.0; cumsum(dists , dims =1)]
    total_length = cum_dist[end]
    
    # === 3. 均匀采样 ===
    # 目标弧长位置
    target_dists = range(0, total_length, length=n_samples)
    
    selected_indices_in_ordered = Int64[]
    
    # 对每个目标位置，寻找原数据中弧长最接近的那个点
    # 这里使用简单的搜索，因为 cum_dist 是单调的
    for target in target_dists
        # 找到 cum_dist 中与 target 差值最小的索引
        idx_in_ordered = argmin(abs.(cum_dist .- target))
        push!(selected_indices_in_ordered, idx_in_ordered[1])
    end
    
    println(selected_indices_in_ordered )
    # 映射回原始索引
    final_indices = indices_ordered[selected_indices_in_ordered]
    
    # 去重（可选）：如果采样过密，同一个点可能被选中两次
    unique!(final_indices)
    
    return final_indices
end

end

