using Documenter, DocumenterMarkdown
using ReactiveDynamics

makedocs(format = Documenter.HTML(prettyurls = false, edit_link = nothing),
         sitename = "ReactiveDynamics.jl API Documentation",
         build = "build_html", pages = ["index.md"])

makedocs(format = Markdown(), sitename = "ReactiveDynamics.jl API Documentation",
         build = "build_md", pages = ["index.md"])
