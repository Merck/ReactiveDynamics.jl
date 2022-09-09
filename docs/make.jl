using Documenter, DocumenterMarkdown
using ReactionDynamics

# makedocs(sitename="ReactionDynamics Interface Documentation") # HTML output
makedocs(format = Markdown(), sitename="ReactionDynamics Documentation", pages = [
    "index.md"]) # MD output
# makedocs(format = Documenter.LaTeX(platform="none"), sitename="ReactionDynamics Documentation") # PDF output