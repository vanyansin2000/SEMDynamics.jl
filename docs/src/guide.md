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

## 坐标变换

```@example guide
μ = aux.EMRot.μ
inertial = [0.2, -0.1, 0.03, 0.04]
rotating = cr3bp_inertial_to_rotating(μ, 0.3, inertial; center=:p2)
cr3bp_rotating_to_inertial(μ, 0.3, rotating; center=:p2)
```

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

## 动力学事件

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
