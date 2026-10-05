# 使用指南

## 动力学方程

`cr3bp_eqm!` 和 `bcr4bp_eqm!` 根据状态长度分派到平面、空间或带状态转移矩阵的方程。
平面状态采用 `[x, y, vx, vy]`，空间状态采用 `[x, y, z, vx, vy, vz]`。

```@example guide
using SEMDynamics

aux = Bcr4bp_Aux()
u = [0.8, 0.1, 0.02, -0.03]
du = similar(u)
cr3bp_eqm!(du, u, aux, 0.0)
du
```

## 积分返回值和终止

`integration` 默认保持三元返回值 `(final_state, times, states)`。求解失败或 callback
主动终止均返回三个 `nothing`，适合打靶时拒绝碰撞轨迹。需要保留事件终止结果时：

```julia
u0 = [0.8, 0.1, 0.02, -0.03]
events = dynamic_events()
sol = integration(u0, (0.0, 10.0), ode_params(cr3bp_eqm!);
                  cb=cb_p2collision(events), return_solution=true)
sol.retcode  # Success、Terminated 或数值失败的返回码
sol.t[end]   # 实际终端时间；事件类别、时间和状态见 events
```

`return_solution=true` 原样返回 SciML 的 `ODESolution`，包括失败时的部分解，
由调用者判断是否可用；不会执行等间隔采样。

## 坐标变换

```@example guide
μ = aux.EMRot.μ
inertial = [0.2, -0.1, 0.03, 0.04]
rotating = cr3bp_inertial_to_rotating(μ, 0.3, inertial; center=:p2)
cr3bp_rotating_to_inertial(μ, 0.3, rotating; center=:p2)
```

坐标转换支持四维平面和六维空间物理状态，返回同维数的 `SVector`；不转换 STM。
`compute_ε2` 和 `compute_ε2_dot` 支持 4/6/20/42 维输入并忽略 STM。
能量导数默认采用 CR3BP；四体轨迹需使用
`compute_ε2_dot(u, t, μ, aux; dynamics=bcr4bp_eqm!)`。

## 周期轨道

`generate_DRO` 以给定周期生成平面 DRO。默认参考周期为 `pi`：

```julia
orbit = generate_DRO(P=pi)
orbit.x0
orbit.C
```

`generate_halo` 使用 `branch=:northern` 或 `:southern` 以及 `lp=:L1` 或
`:L2` 选择空间 Halo 轨道族：

```julia
halo = generate_halo(branch=:northern, lp=:L1)
nrho = generate_nrho_9_2()
```

Halo 初值位于 x-z 对称面，形式为 `[x, 0, z, 0, vy, 0]`。传入 `P` 可沿所选
分支延拓到指定的无量纲周期。`generate_nrho_9_2()` 使用月球会合周期定义的
`P=4pi/(9abs(ws))`；Halo 打靶默认容差为 `1e-10`，默认周期步长为 `0.005`。
所有 Halo 与 DRO seed 均须是完整六维状态，且使用同一组 `[x, z, vy]` 打靶
未知量与 `[y, vx, vz] = 0` 半周期条件。

`orbit.sol` 保存 `[-orbit.P, 2orbit.P]` 的稠密解，供相邻周期查询以避免外插。
从 `t=0` 的初值反向传播到 `-P` 后再构造完整缓存，保持初值的零时刻相位；其误差受
打靶和积分误差影响。此约定面向已修正的 CR3BP 周期轨道。打靶容差 `tol` 与最终缓存传播分开，后者使用 `abstol=reltol=1e-12`。

## 动力学事件

`cb_p1collision` 和 `cb_p2collision` 使用 `DiscreteCallback`，在积分步末检查
是否进入 `scale * aux.EMRot.r_p1/r_p2` 所定义的区域，不定位边界根。
`terminate=false` 时，区域内的每个积分步都会记录事件。默认地月参数中的
`aux.dim.r_p2` 为 1937.4 km；碰撞判据使用这个配置值。

连续事件回调支持 `rootfind`，也接受旧拼写 `root_find`；同时提供两者时以
`root_find` 为准。

`cb_enter` 与 `cb_escape` 分别检测向内、向外穿越第二主天体作用球；`scale`
用于缩放默认球半径，`terminate` 控制是否在事件处终止积分：

```julia
events = dynamic_events()
sphere_callbacks = CallbackSet(
    cb_enter(events; terminate=false),
    cb_escape(events; terminate=false),
)
```

第一、第二主天体的近远拱点采用对称接口：

```julia
cb_apse_p1(events; terminate_perigee=false, terminate_apogee=false)
cb_apse_p2(events; terminate_perilune=false, terminate_apolune=false)
```

第二个回调记录 `:perilune` 和 `:apolune`。兼容接口 `cb_perilune` 仍然可用，
并新增了仅检测远月点的 `cb_apolune`。

## 本地构建文档

首次构建时，在仓库根目录执行：

```julia
julia --project=docs -e 'using Pkg; Pkg.develop(PackageSpec(path=pwd())); Pkg.instantiate()'
julia --project=docs docs/make.jl
```

生成页面位于 `docs/build/`。该目录是构建产物，不应提交到 Git。
