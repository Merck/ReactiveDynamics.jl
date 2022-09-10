using Documenter, DocumenterMarkdown
using ReactiveDynamics

# makedocs(sitename="ReactiveDynamics Interface Documentation") # HTML output
makedocs(format = Markdown(), sitename="ReactiveDynamics Documentation", pages = [
    "index.md"]) # MD output
# makedocs(format = Documenter.LaTeX(platform="none"), sitename="ReactiveDynamics Documentation") # PDF output