using Documenter, DocumenterMarkdown
using ReactiveDynamics

makedocs(;
    format = Documenter.HTML(; prettyurls = false, edit_link = "main"),
    sitename = "ReactiveDynamics.jl API Documentation",
    pages = ["index.md"],
)

deploydocs(; repo = "github.com/Merck/ReactiveDynamics.jl.git")
