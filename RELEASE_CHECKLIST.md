# SEMDynamics.jl 注册与发布清单

## 当前仓库已具备

- 标准包入口：`src/SEMDynamics.jl`
- 标准测试入口：`test/runtests.jl`
- 独立 Documenter 环境：`docs/Project.toml`、`docs/make.jl`、`docs/src/`
- GitHub Actions：最低支持 Julia 1.10 与最新稳定版的测试，以及文档构建/部署
- OSI 批准的 MPL-2.0 `LICENSE`
- 所有直接依赖和 Julia 的有限上界 `[compat]`

## 首次注册前

1. 确认 GitHub 仓库公开，默认分支为 `main`，URL 为
   `https://github.com/vanyansin2000/SEMDynamics.jl.git`。
2. 检查 `Project.toml` 中的 `name`、`uuid`、`version`、`authors` 和依赖。
   首次注册版本建议使用当前语义化版本，不要在注册提交后覆盖同一版本标签。
3. 推送本仓库改动，确认 `CI` 和 `Documentation` 工作流通过。
4. 在仓库 Settings → Pages 中允许 GitHub Actions/`gh-pages` 发布页面；首次成功部署后检查
   `https://vanyansin2000.github.io/SEMDynamics.jl/dev/`。
5. 在本地或 CI 中运行：

   ```julia
   using Pkg
   Pkg.test()
   ```

6. 在 Julia REPL 中临时安装仓库，确认不依赖开发目录中的 `Manifest.toml`：

   ```julia
   using Pkg
   Pkg.activate(; temp=true)
   Pkg.add(url="https://github.com/vanyansin2000/SEMDynamics.jl")
   using SEMDynamics
   ```

7. 安装 GitHub Apps：
   [JuliaRegistrator](https://github.com/JuliaRegistries/Registrator.jl) 和
   [Julia TagBot](https://github.com/JuliaRegistries/TagBot)。
8. 在目标提交或关联 issue 下评论：

   ```text
   @JuliaRegistrator register
   ```

9. 跟踪 Registrator 创建的 General PR，处理 RegistryCI 的 AutoMerge 检查。
   新包存在等待期属于正常流程。

## 每次发布前

1. 更新 `Project.toml` 的版本号，遵循语义化版本。
2. 更新 README、文档、测试和变更记录（后续建议新增 `CHANGELOG.md`）。
3. 运行 `Pkg.test()` 和 `julia --project=docs docs/make.jl`。
4. 推送后等待 CI 通过，再触发 Registrator。
5. 注册合并后由 TagBot 创建版本标签和 GitHub Release；确认稳定版文档已部署。

## 后续可选增强

- 增加 Aqua.jl 的包质量测试，以及 JET.jl 的静态分析。
- 增加 Codecov 覆盖率上传。
- 为绘图 API 增加无头环境测试或参考图测试。
- 将较重的可视化依赖拆分为 Julia package extension，以缩短基础包加载与 CI 时间。
- 为主要公开接口补充英文或中英双语 docstring，以方便 General 用户检索。
