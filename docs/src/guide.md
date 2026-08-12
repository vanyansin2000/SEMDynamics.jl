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

## 本地构建文档

首次构建时，在仓库根目录执行：

```julia
julia --project=docs -e 'using Pkg; Pkg.develop(PackageSpec(path=pwd())); Pkg.instantiate()'
julia --project=docs docs/make.jl
```

生成页面位于 `docs/build/`。该目录是构建产物，不应提交到 Git。
