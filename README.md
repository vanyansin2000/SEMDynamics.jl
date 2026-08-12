# SEMDynamics.jl

| **Documentation** | **Tests** | **License** |
|:---:|:---:|:---:|
| [![Stable documentation](https://img.shields.io/badge/docs-stable-blue.svg)](https://vanyansin2000.github.io/SEMDynamics.jl/stable/) [![Development documentation](https://img.shields.io/badge/docs-dev-blue.svg)](https://vanyansin2000.github.io/SEMDynamics.jl/dev/) | [![CI](https://github.com/vanyansin2000/SEMDynamics.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/vanyansin2000/SEMDynamics.jl/actions/workflows/CI.yml) | [![License: MPL-2.0](https://img.shields.io/badge/license-MPL--2.0-green.svg)](LICENSE) |

`SEMDynamics.jl` 是一个用于地月受限多体动力学数值计算的 Julia 包。当前提供圆型受限三体问题（CR3BP）、地月旋转系双圆四体问题（BCR4BP）的动力学方程、事件检测、坐标转换、轨迹可视化，以及按指定周期生成平面远距逆行轨道（DRO）的工具。

> 本项目使用无量纲的旋转坐标系；除非另有说明，时间和周期均为无量纲量。

## 安装

在 Julia REPL 中切换到包模式（按 `]`）后执行：

```julia
pkg> develop path/to/semdynamics
pkg> instantiate
```

或在项目根目录中：

```julia
julia --project=.
```

然后运行：

```julia
using Pkg
Pkg.instantiate()
using SEMDynamics
```

## 快速开始

### 按周期生成 DRO

`generate_DRO` 以 2:1 平面 DRO 的近似初值为默认 seed。传入目标周期 `P` 后，函数对每个固定周期执行半周期对称打靶，返回一个 `PeriodicOrbit`：

```julia
using SEMDynamics

dro_2_1 = generate_DRO(P=pi)
x0 = dro_2_1.x0       # 六维旋转系初值
period = dro_2_1.P    # pi
jacobi = dro_2_1.C
```

当目标周期偏离 `pi` 时，函数在内存中按周期做短步长延拓：每一步都用前两步的割线预测初值，再执行 Newton 打靶修正。因此不需要预先生成、保存或读取完整的轨道族。

```julia
dro = generate_DRO(P=3.4, continuation_step=0.05)
```

`continuation_step` 越小通常越稳健，但会增加打靶次数。若已经有更接近目标族的二维 seed `[x, vy]`，可显式指定其周期：

```julia
dro = generate_DRO(
    P=3.4,
    seed=[1.1754, -0.4943],
    seed_period=pi,
)
```

目前该接口只实现平面 CR3BP DRO，初值形式为 `[x, 0, 0, 0, vy, 0]`。通用的
`orbit_shooting(residual!, initial_guess, parameters)` 已保留，未来轨道族可通过定义自己的残差、参数化和 seed 接入，而无需重新设计求根层。

### 轨迹积分和事件检测

```julia
using SEMDynamics
using DifferentialEquations

aux = Bcr4bp_Aux()
x0 = [1.05, 0.0, 0.0, 0.2]
events = dynamic_events()
callback = CallbackSet(cb_p2collision(events), cb_escape(events))

problem = ODEProblem(cr3bp_eqm!, x0, (0.0, 10.0), aux)
solution = solve(problem, Vern7(); callback, reltol=1e-12, abstol=1e-12)
```

### 可视化

```julia
using CairoMakie

figure = Figure()
axis = Axis(figure[1, 1], aspect=AxisAspect(1))
cr3bp_set_scn!(axis)
plotPO!(axis, dro_2_1; color=:crimson, linewidth=2)
display(figure)
```

对于二维积分结果，也可使用 `plot_traj_rot!`、`plot_traj_inertial!` 和 `plot_traj_ε2!`。

## 主要接口

| 接口 | 说明 |
| --- | --- |
| `cr3bp_eqm!` | CR3BP 平面、空间和 STM 状态方程的分发入口。 |
| `bcr4bp_eqm!` | 地月旋转系 BCR4BP 方程；按状态长度自动处理 2D、3D、2D+STM 和 3D+STM。 |
| `bcr4bp_eqmEMRot2D!` | 旧版二维 BCR4BP 接口，保留以兼容既有代码。 |
| `Bcr4bp_Aux()` | 创建物理常数与无量纲归一化参数。 |
| `generate_DRO(P=...)` | 按目标周期打靶生成平面 DRO。 |
| `orbit_shooting(...)` | 可扩展的通用非线性打靶求解入口。 |
| `compute_jacobi` | 计算 CR3BP Jacobi 常数。 |
| `cb_*` | 碰撞、逃逸、能量零交叉和近/远拱点事件回调。 |
| `cr3bp_inertial_to_rotating` / `cr3bp_rotating_to_inertial` | CR3BP 惯性系与旋转系的状态转换。 |

可通过 Julia REPL 的帮助模式查看每个导出的函数签名和参数说明：

```julia
help?> generate_DRO
```

## 项目结构

```text
src/
├── SEMDynamics.jl       # 顶层模块与重导出
├── Dynamics.jl          # 动力学方程和积分参数
├── dynamic_parameters.jl# 物理参数与无量纲尺度
├── dynamic_energy.jl    # 能量、Jacobi 常数与平衡点
├── dynamic_events.jl    # DifferentialEquations 事件回调
├── rotation.jl          # 坐标和相位转换
├── PeriodicOrbits.jl    # 平面 DRO 打靶与绘图
├── Visualization.jl     # 轨迹与场景绘图
└── Utils.jl             # 曲线交点与采样工具
test/
└── runtests.jl          # 标准测试入口
docs/
├── make.jl              # Documenter 构建入口
└── src/                 # 文档 Markdown 源文件
examples/
└── periodic_orbits.jl   # 周期轨道与绘图示例
```

## 数值约定与限制

- CR3BP/BCR4BP 状态按 `[x, y, z, vx, vy, vz]`（空间）或 `[x, y, vx, vy]`（平面）排列；加入 STM 后总长度分别为 42 和 20。
- `generate_DRO` 的 `P` 是完整周期，不是半周期。
- DRO 打靶使用 `y(P/2)=0` 与 `vx(P/2)=0` 的平面对称条件。
- 默认误差容限为 `1e-12`；对较远目标周期，建议减小 `continuation_step` 并检查最终轨道闭合误差。
- 目前未提供其他周期轨道族的专用生成器，也不再维护轨道数据文件缓存。

## 开发与验证

在项目根目录运行标准测试：

```sh
julia --project=. -e 'using Pkg; Pkg.test()'
```

本地构建 Documenter 文档：

```sh
julia --project=docs -e 'using Pkg; Pkg.develop(PackageSpec(path=pwd())); Pkg.instantiate()'
julia --project=docs docs/make.jl
```

更完整的发布检查见 [`RELEASE_CHECKLIST.md`](RELEASE_CHECKLIST.md)。

## 许可证

本项目采用仓库中的 [MPL-2.0 许可证](LICENSE)。
