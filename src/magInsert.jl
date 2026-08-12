

using CairoMakie

export MagInset

# ─────────────────────────────────────────────────────────────────────────────
# 数据坐标 → Figure 归一化坐标
# ─────────────────────────────────────────────────────────────────────────────
"""将 `ax` 内的数据坐标转换为所属 Figure 的归一化坐标。"""
function data2norm(ax::Axis, x::Real, y::Real)
    fig_vp = ax.parent.scene.viewport[]
    ax_vp  = ax.scene.viewport[]
    xlims  = ax.finallimits[].origin[1], ax.finallimits[].origin[1] + ax.finallimits[].widths[1]
    ylims  = ax.finallimits[].origin[2], ax.finallimits[].origin[2] + ax.finallimits[].widths[2]
    px = ax_vp.origin[1] + ax_vp.widths[1] * (x - xlims[1]) / (xlims[2] - xlims[1])
    py = ax_vp.origin[2] + ax_vp.widths[2] * (y - ylims[1]) / (ylims[2] - ylims[1])
    nx = (px - fig_vp.origin[1]) / fig_vp.widths[1]
    ny = (py - fig_vp.origin[2]) / fig_vp.widths[2]
    return nx, ny
end

"""返回归一化边界框中 `"NW"`、`"SE"` 等角标对应的位置。"""
function corner_pos(bbox_norm::NTuple{4,Float64}, corner::String)
    xmin, xmax, ymin, ymax = bbox_norm
    cx = corner[2] == 'W' ? xmin : xmax
    cy = corner[1] == 'S' ? ymin : ymax
    return cx, cy
end

# ─────────────────────────────────────────────────────────────────────────────
# MagInset（不做任何改动）
#
# 新增说明：
#   ★ inset_ax 支持两种模式：
#     1. 默认：内部自动用 bbox 创建浮动插图（原有行为）
#     2. 传入 inset_ax 参数：使用外部已创建的 Axis，只向其中填充数据
#        用法：MagInset(fig, ax, zoom_area, inset_pos, lines; inset_ax=my_ax)
# ─────────────────────────────────────────────────────────────────────────────
"""
    MagInset(fig, ax, zoom_area, inset_pos, conn_lines=[]; kwargs...) -> Axis

在主轴 `ax` 的 `zoom_area` 加边框，并创建或填充一个显示该区域的放大轴。可用
`conn_lines` 连接两处角点，或通过 `inset_ax` 使用调用方预先创建的插图轴。
"""
function MagInset(
    fig::Figure,
    ax::Axis,
    zoom_area::NTuple{4,<:Real},
    inset_pos::NTuple{4,<:Real},
    conn_lines::Vector{Tuple{String,String}} = Tuple{String,String}[];
    inset_xscale::Function  = identity,
    inset_yscale::Function  = identity,
    zoom_rect_color         = :black,
    zoom_rect_lw            = 1.5,
    connect_line_lw         = 1.0,
    connect_line_color      = :black,
    connect_line_style      = :dash,
    inset_ax::Union{Axis,Nothing} = nothing,
    inset_kwargs...
)
    xmin_z, xmax_z, ymin_z, ymax_z = Float64.(zoom_area)

    # ★★ 关键修复 0：先让主轴完成一次正常的 limits 计算，然后锁死它。
    #    锁死后 update_state_before_display! 不会再对主轴跑 autolimits，
    #    也就不会对任何 plot 的 boundingbox 做 log10 而触发 DomainError。
    CairoMakie.Makie.update_state_before_display!(fig)  # 此刻主轴还没有方框，正常计算
    cur = ax.finallimits[]                               # 当前主轴真实显示范围
    ax.limits[] = (cur.origin[1], cur.origin[1] + cur.widths[1],
                   cur.origin[2], cur.origin[2] + cur.widths[2])  # ← 显式锁死主轴 limits

    # ── 1. 复制前快照（不含方框）──
    plots_to_copy = copy(ax.scene.plots)

    # ── 2. 画方框（此时主轴 limits 已锁，方框不会触发主轴 autolimits）──
    rect_x = [xmin_z, xmax_z, xmax_z, xmin_z, xmin_z]
    rect_y = [ymin_z, ymin_z, ymax_z, ymax_z, ymin_z]
    rect_plot = lines!(ax, rect_x, rect_y; color = zoom_rect_color, linewidth = zoom_rect_lw)
    translate!(rect_plot, 0, 0, 100)

    # ── 3. 确定 inset_ax ──
    if isnothing(inset_ax)
        nx1, ny1 = data2norm(ax, inset_pos[1], inset_pos[3])
        nx2, ny2 = data2norm(ax, inset_pos[2], inset_pos[4])
        fig_px   = fig.scene.viewport[]
        W, H     = fig_px.widths
        bbox     = BBox(nx1 * W, nx2 * W, ny1 * H, ny2 * H)
        inset_ax = Axis(fig;
            bbox         = bbox,
            xscale       = inset_xscale,
            yscale       = inset_yscale,
            xgridvisible = false,
            ygridvisible = false,
            inset_kwargs...
        )
        _inset_pos_for_lines = inset_pos
    else
        inset_ax.xscale[] = inset_xscale
        inset_ax.yscale[] = inset_yscale
        _inset_pos_for_lines = nothing
    end

    # ── 4. 复制数据（快照，跳过方框）──
    for plt in plots_to_copy
        plt === rect_plot && continue
        _copy_plot_to_axis!(inset_ax, plt)
    end

    # ── 5. 设置 inset 范围（显式锁，避免 inset 也跑 autolimits）──
    #     inset y 是 log，xmin_z/xmax_z 是负的→x 必须 identity；ymin_z/ymax_z 正的→y log 安全
    xlims!(inset_ax, xmin_z, xmax_z)
    ylims!(inset_ax, ymin_z, ymax_z)

    # ── 6. 连接线 ──
    if !isempty(conn_lines) && !isnothing(_inset_pos_for_lines)
        CairoMakie.Makie.update_state_before_display!(fig)
        fig_px = fig.scene.viewport[]
        W, H   = fig_px.widths
        zx1, zy1 = data2norm(ax, xmin_z, ymin_z)
        zx2, zy2 = data2norm(ax, xmax_z, ymax_z)
        zoom_bbox_norm = (zx1, zx2, zy1, zy2)
        ix1, iy1 = data2norm(ax, _inset_pos_for_lines[1], _inset_pos_for_lines[3])
        ix2, iy2 = data2norm(ax, _inset_pos_for_lines[2], _inset_pos_for_lines[4])
        inset_bbox_norm = (ix1, ix2, iy1, iy2)
        for (z_corner, i_corner) in conn_lines
            lx1, ly1 = corner_pos(zoom_bbox_norm,  z_corner)
            lx2, ly2 = corner_pos(inset_bbox_norm, i_corner)
            linesegments!(fig.scene,
                [Point2f(lx1 * W, ly1 * H), Point2f(lx2 * W, ly2 * H)];
                color = connect_line_color, linewidth = connect_line_lw,
                linestyle = connect_line_style,
            )
        end
    end

    return inset_ax
end

"""按已支持的 Makie 图元类型复制 `plt` 到插图轴，避免在两个 scene 间共享图元。"""
function _copy_plot_to_axis!(ax::Axis, plt)
    T = typeof(plt)

    # ── Lines / LineSegments ───────────────────────────────────────────────
    if T <: Makie.Lines || T <: Makie.LineSegments
        lines!(ax, plt[1][];
            color     = plt.color[],
            linewidth = plt.linewidth[],
            linestyle = plt.linestyle[],   # ★ 原版漏掉了 linestyle
        )

    # ── Scatter ────────────────────────────────────────────────────────────
    elseif T <: Makie.Scatter
        scatter!(ax, plt[1][];
            color       = plt.color[],
            markersize  = plt.markersize[],
            marker      = plt.marker[],
            strokecolor = plt.strokecolor[],
            strokewidth = plt.strokewidth[],
        )

    # ── Poly（Circle / Polygon 等几何体）── ★ 新增 ─────────────────────────
    elseif T <: Makie.Poly
        poly!(ax, plt[1][];
            color       = plt.color[],
            strokecolor = plt.strokecolor[],
            strokewidth = plt.strokewidth[],
        )

    # ── Band ───────────────────────────────────────────────────────────────
    elseif T <: Makie.Band
        band!(ax, plt[1][], plt[2][], plt[3][];
            color = plt.color[],
        )

    # ── Text 标注 ── ★ 新增 ───────────────────────────────────────────────
    elseif T <: Makie.Text
        text!(ax, plt[1][];
            position  = plt.position[],
            color     = plt.color[],
            fontsize  = plt.fontsize[],
            align     = plt.align[],
        )

    # ── Heatmap ────────────────────────────────────────────────────────────
    elseif T <: Makie.Heatmap
        heatmap!(ax, plt[1][], plt[2][], plt[3][];
            colormap  = plt.colormap[],
            colorrange = plt.colorrange[],
        )

    # ── Mesh / Surface ─────────────────────────────────────────────────────
    # ★ 原版 plot!(ax.scene, plt) 会将同一个 plot 对象挂载到两个 scene，
    #   导致渲染时 transformation 错乱或双重绘制。
    #   改为按类型重新绘制，避免共享同一 Observable。
    elseif T <: Makie.Surface
        surface!(ax, plt[1][], plt[2][], plt[3][];
            colormap = plt.colormap[],
        )

    elseif T <: Makie.Mesh
        mesh!(ax, plt[1][];
            color    = plt.color[],
        )

    # ── 兜底：跳过未知类型（如内部辅助 plot）─────────────────────────────
    # else
    #     @warn "MagInset: 未处理的 plot 类型 $(T)，已跳过"
    end
end


# # ═════════════════════════════════════════════════════════════════════════════
# #  演示示例
# # ═════════════════════════════════════════════════════════════════════════════
# function demo()
#     # 生成演示数据
#     x  = range(0, 4π, length = 400)
#     y1 = sin.(x)
#     y2 = sin.(x) .* exp.(-0.15 .* x)

#     fig = Figure(size = (800, 500))
#     ax  = Axis(fig[1, 1];
#         xlabel = "x",
#         ylabel = "y",
#         title  = "CairoMakie")

#     lines!(ax, x, y1; color = :steelblue,  linewidth = 2, label = "sin(x)")
#     lines!(ax, x, y2; color = :tomato,     linewidth = 2, label = "sin(x)·e^(-0.15x)")
#     axislegend(ax; position = :rt)

#     # 放大区域（数据坐标）：x ∈ [5.5, 7.5], y ∈ [-1.05, 0.3]
#     zoom_area  = (5., 10., -1., 0.)
#     # 插图位置（数据坐标）：右上角区域
#     inset_pos  = (0.0, 2.5, -1.0, -0.5)

#     # 连接线：zoom 区域的 NW 角 → inset 的 SW 角
#     #          zoom 区域的 SW 角 → inset 的 SE 角（可按需调整）
#     conn_lines = [("NW", "SW"), ("NE", "SE")]

#     inset_ax = MagInset(fig, ax, zoom_area, inset_pos, conn_lines;
#         zoom_rect_color    = :black,
#         zoom_rect_lw       = 1.5,
#         connect_line_lw    = 1.0,
#         connect_line_color = :gray40,
#         connect_line_style = :dash,
#         # 以下为 inset Axis 的额外参数
#         xticklabelsize = 9,
#         yticklabelsize = 9,
#         spinewidth     = 1.2,
#     )

#     save("mag_inset_demo.pdf", fig)
#     save("mag_inset_demo.png", fig, px_per_unit = 2)
#     display(fig)
#     return fig, inset_ax
# end

# demo()
