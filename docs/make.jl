using Documenter
using SEMDynamics

DocMeta.setdocmeta!(SEMDynamics, :DocTestSetup, :(using SEMDynamics); recursive=true)

makedocs(
    sitename="SEMDynamics.jl",
    modules=[SEMDynamics],
    format=Documenter.HTML(
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://vanyansin2000.github.io/SEMDynamics.jl/stable/",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "首页" => "index.md",
        "使用指南" => "guide.md",
        "API 参考" => "api.md",
    ],
    checkdocs=:exports,
)

deploydocs(
    repo="github.com/vanyansin2000/SEMDynamics.jl.git",
    devbranch="main",
    push_preview=true,
)
