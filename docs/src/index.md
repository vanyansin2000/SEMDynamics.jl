# SEMDynamics.jl

`SEMDynamics.jl` 提供地月受限多体动力学数值计算工具，包括 CR3BP、BCR4BP、事件检测、
坐标变换、轨迹可视化和周期轨道计算。

```@contents
Pages = ["tutorial.md", "guide.md", "api.md"]
Depth = 2
```

## 安装

在包注册前，可直接从仓库安装：

```julia
using Pkg
Pkg.add(url="https://github.com/vanyansin2000/SEMDynamics.jl")
```

注册后，可使用标准包安装命令：

```julia
using Pkg
Pkg.add("SEMDynamics")
```

## 最小示例

```@example quickstart
using SEMDynamics

aux = Bcr4bp_Aux()
state = [0.8, 0.1, 0.02, -0.03]
derivative = similar(state)
cr3bp_eqm!(derivative, state, aux, 0.0)
derivative
```

## 文档约定

函数说明直接取自源码中的 Julia docstring。新增公开接口时，请在函数或类型定义前编写 docstring，
并将其加入 [API 参考](@ref)。

第一次使用本包时，建议从[入门教程](@ref)开始；教程覆盖积分、事件、坐标转换、
周期轨道和二维/三维绘图的完整基础流程。
