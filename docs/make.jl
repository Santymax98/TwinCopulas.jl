using TwinCopulas
using Documenter

DocMeta.setdocmeta!(TwinCopulas, :DocTestSetup, :(using TwinCopulas); recursive=true)

makedocs(
    modules = [TwinCopulas],
    authors = "Santiago Jimenez Ramos",
    sitename = "TwinCopulas.jl",
    format = Documenter.HTML(
        canonical = "https://Santymax98.github.io/TwinCopulas.jl",
        edit_link = "main",
        assets = String[],
    ),
    pages = [
        "Home" => ["index.md", "starting.md"],
        "Guide" => [
            "Getting started" => "Copula_Sklar.md",
            "Archimedean Copulas" => "Archimedean/Archimedean_theory.md",
            "Elliptical Copulas" => "Elliptical/Elliptical_theory.md",
            "Extreme Value Copulas" => "Extreme/Extreme_Value_theory.md",
        ],
        "Available Models" => [
            "Archimedean Copulas" => "Archimedean/Available_Archimedean_models.md",
            "Elliptical Copulas" => "Elliptical/Available_Elliptical_models.md",
            "Extreme Value Copulas" => "Extreme/Available_Extreme_models.md",
        ]
    ],
)

deploydocs(
    repo = "github.com/Santymax98/TwinCopulas.jl",
    branch = "main",
    devbranch = "main",
    target = "site",
    forcepush = true
)

deploydocs(;
    repo="github.com/Santymax98/TwinCopulas.jl",
    devbranch="main",
)
