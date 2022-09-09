using Documenter, DocumenterMarkdown
using DyVE

# makedocs(sitename="DyVE Interface Documentation") # HTML output
makedocs(format = Markdown(), sitename="DyVE Documentation", pages = [
    "index.md"]) # MD output
# makedocs(format = Documenter.LaTeX(platform="none"), sitename="DyVE Documentation") # PDF output