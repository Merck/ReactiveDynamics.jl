# ════════════════════════════════════════════════════════════════════════════════════════
# make.jl — build the ReactiveDynamics.jl documentation site.
#
#   docs/literate/**/*.jl  (Literate sources, migrated from demo/)  →  docs/src/**/*.md
#   docs/src/**/*.md        (+ autodocs/explanation)                 →  docs/build/  (HTML)
#
# This replaces the old DocumenterMarkdown backend. The Literate pre-pass runs FIRST, with
# `execute = false` here (Documenter re-executes the generated `@example` blocks at HTML-build
# time, so output is captured into the page and a code error fails the build). To preview a
# single tutorial as a self-contained HTML page WITHOUT the whole site — the fast authoring
# loop and the exemplar's acceptance check — use build_literate.jl instead.
#
#   First time:  julia --project=docs -e 'using Pkg; Pkg.instantiate()'
#   Build:       julia --project=docs docs/make.jl
#
# NB: this is the site scaffold the DOCS_CHARTER.md tracker (facet E1) builds toward. As
# tutorials/case-studies/reference/explanation pages land, extend LITERATE_TUTORIALS and the
# `pages` tree below; the charter's §10 tracker is the source of truth for what is wired yet.
# ════════════════════════════════════════════════════════════════════════════════════════

using Documenter, Literate
using ReactiveDynamics

const HERE = @__DIR__
const LITERATE_DIR = joinpath(HERE, "literate")
const GEN_DIR = joinpath(HERE, "src")

# ── Literate pre-pass ──────────────────────────────────────────────────────────────────
# Convert each Literate `.jl` to a Documenter-flavored markdown page under docs/src/**, so
# the code blocks become `@example` blocks Documenter executes during the HTML build.
# `documenter = true` emits the Documenter `@example`/`@meta` flavor. Add entries here as the
# charter's tutorial/case-study facets land.
const LITERATE_TUTORIALS = [
    # (source relative to docs/literate,           output subdir under docs/src)
    ("tutorials/introductory.jl", "tutorials"),
]

for (src, outsub) in LITERATE_TUTORIALS
    Literate.markdown(
        joinpath(LITERATE_DIR, src), joinpath(GEN_DIR, outsub);
        documenter = true, execute = false, credit = false,
    )
end

# ── Site ────────────────────────────────────────────────────────────────────────────────
makedocs(;
    sitename = "ReactiveDynamics.jl",
    format = Documenter.HTML(; prettyurls = get(ENV, "CI", "false") == "true", edit_link = "main"),
    modules = [ReactiveDynamics],
    warnonly = true,   # scaffold stage: don't fail on cross-references to pages not yet authored
    pages = [
        "Home" => "index.md",
        "Tutorials" => [
            "Introductory" => "tutorials/introductory.md",
            # "Advanced"  => "tutorials/advanced.md",     # charter A2
            # "Expert"    => "tutorials/expert.md",       # charter A3
        ],
        # "Case studies" => [...],                         # charter B1–B3
        # "Reference"    => [...],                         # charter C1–C2
        # "Explanation"  => [...],                         # charter D1–D2
    ],
)

deploydocs(; repo = "github.com/Merck/ReactiveDynamics.jl.git")
