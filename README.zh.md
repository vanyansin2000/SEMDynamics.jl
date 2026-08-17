# SEMDynamics.jl

[English](README.md) | 中文

| 文档 | 构建 | 许可证 |
|:---:|:---:|:---:|
| [![稳定版文档](https://img.shields.io/badge/docs-stable-blue.svg)](https://vanyansin2000.github.io/SEMDynamics.jl/stable/) [![开发版文档](https://img.shields.io/badge/docs-dev-blue.svg)](https://vanyansin2000.github.io/SEMDynamics.jl/dev/) | [![CI](https://github.com/vanyansin2000/SEMDynamics.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/vanyansin2000/SEMDynamics.jl/actions/workflows/CI.yml) | [![许可证：MPL-2.0](https://img.shields.io/badge/license-MPL--2.0-green.svg)](LICENSE) |

`SEMDynamics.jl` 是一个用于地月受限多体问题数值轨迹分析的 Julia 包。它提供圆型受限三体问题
（CR3BP）和地月双圆受限四体问题（BCR4BP）的运动方程，以及状态转移矩阵传播、事件回调、
坐标变换、轨迹绘图和平面远距逆行轨道（DRO）的微分修正功能。

除非另有说明，状态均以无量纲旋转坐标系表示，时间和周期也均为无量纲量。

## 安装

本包已经注册到 Julia General Registry，可直接安装：

```julia
using Pkg
Pkg.add("SEMDynamics")
```

SEMDynamics 要求 Julia 1.12 或更高版本。

## 快速开始

积分一条平面 CR3BP 轨迹：

```julia
using SEMDynamics

u0 = [0.8, 0.1, 0.02, -0.03] # [x, y, vx, vy]
ode_args = (; reltol = 1e-12 , abstol = 1e-12) # 默认值
uf , ts , us = integration(u0, (0.0, 10.0), ode_params(cr3bp_eqm! , ode_args , Bcr4bp_Aux()))
# 或者简写为
uf , ts , us = integration(u0, (0.0, 10.0), ode_params(cr3bp_eqm!))
```

### 生成平面 DRO

`generate_DRO` 使用微分修正生成具有指定完整周期的平面 CR3BP 远距逆行轨道。默认种子是周期为
`pi` 的 2:1 DRO。

```julia
using SEMDynamics

orbit = generate_DRO(P=pi)
orbit.x0 # 旋转坐标系中的六维初始状态
orbit.P  # 无量纲完整周期
orbit.C  # Jacobi 常数
```

当目标周期偏离参考种子时，求解器会沿周期进行短步长延拓，每一步使用割线预测和 Newton 修正：

```julia
orbit = generate_DRO(P=3.4, continuation_step=0.05)
```

较小的 `continuation_step` 通常更稳健，但需要更多次修正。也可以显式提供更接近目标轨道的二元种子
`[x, vy]` 及其对应周期：

```julia
orbit = generate_DRO(
    P=3.4,
    seed=[1.1754, -0.4943],
    seed_period=pi,
)
```

目前高层生成器仅支持平面 CR3BP DRO，其初始状态为 `[x, 0, 0, 0, vy, 0]`。可以使用底层接口
`orbit_shooting(residual!, initial_guess, parameters)` 为其他轨道族定义微分修正。

### 检测事件

本包提供与 DifferentialEquations.jl 兼容的回调函数。

```julia
using SEMDynamics
using DifferentialEquations

aux = Bcr4bp_Aux()
events = dynamic_events()
callback = CallbackSet(cb_p2collision(events), cb_escape(events))

u0 = [1.05, 0.0, 0.0, 0.2]
problem = ODEProblem(bcr4bp_eqm!, u0, (0.0, 10.0), aux)
solution = solve(problem, Vern7(); callback, reltol=1e-12, abstol=1e-12)
```

### 绘制轨道

```julia
using CairoMakie
using SEMDynamics

orbit = generate_DRO(P=pi)
figure = Figure()
axis = Axis(figure[1, 1], aspect=AxisAspect(1))
cr3bp_set_scn!(axis)
plotPO!(axis, orbit; color=:crimson, linewidth=2)
display(figure)
```

## 状态约定

| 模型 | 状态 | 带 STM 的状态 |
| --- | --- | --- |
| 平面 CR3BP/BCR4BP | `[x, y, vx, vy]` | 20 个元素：状态之后为 4x4 STM |
| 空间 CR3BP/BCR4BP | `[x, y, z, vx, vy, vz]` | 42 个元素：状态之后为 6x6 STM |

`cr3bp_eqm!` 和 `bcr4bp_eqm!` 根据状态长度自动分派。对于 CR3BP，`compute_jacobi` 用于计算
Jacobi 常数。

## 主要接口

| 接口 | 说明 |
| --- | --- |
| `Bcr4bp_Aux()` | 创建地月模型的物理常数和无量纲参数。 |
| `cr3bp_eqm!` | 平面、空间及带 STM 状态的 CR3BP 运动方程。 |
| `bcr4bp_eqm!` | 平面、空间及带 STM 状态的地月 BCR4BP 运动方程。 |
| `generate_DRO(P=...)` | 生成经过微分修正的平面 CR3BP DRO。 |
| `orbit_shooting(...)` | 通用非线性打靶接口。 |
| `compute_jacobi` | 计算 CR3BP Jacobi 常数。 |
| `cb_*` | 碰撞、逃逸、能量穿越和拱点事件回调。 |
| `cr3bp_inertial_to_rotating` / `cr3bp_rotating_to_inertial` | 在惯性系和旋转系之间转换 CR3BP 状态。 |
| `plotPO!`、`plot_traj_rot!`、`plot_traj_inertial!` | 用于周期轨道或积分轨迹的 Makie 绘图工具。 |

可以在 Julia 帮助模式中查看导出方法及其参数：

```julia
help?> generate_DRO
```

## 文档

在线文档位于
[vanyansin2000.github.io/SEMDynamics.jl](https://vanyansin2000.github.io/SEMDynamics.jl/stable/)。

在仓库根目录中执行以下命令可在本地构建文档：

```sh
julia --project=docs -e 'using Pkg; Pkg.develop(PackageSpec(path=pwd())); Pkg.instantiate()'
julia --project=docs docs/make.jl
```

生成的网站写入 `docs/build/`。

## 测试

在仓库根目录运行测试套件：

```sh
julia --project=. -e 'using Pkg; Pkg.test()'
```

## 贡献

欢迎提交 issue 和 pull request。行为变更请同时提供测试，并使用 Julia docstring 记录公开 API。

## 许可证

SEMDynamics.jl 按照 [Mozilla Public License 2.0](LICENSE) 发布。
