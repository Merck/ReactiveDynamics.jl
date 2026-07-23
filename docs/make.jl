# ════════════════════════════════════════════════════════════════════════════════════════
# make.jl — build the ReactiveDynamics.jl documentation site.
#
#   docs/literate/**/*.jl  (Literate sources, migrated from demo/)  →  docs/src/**/*.md
#   docs/src/**/*.md        (+ autodocs)                             →  docs/build/  (HTML)
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
# tutorials/case-studies/reference pages land, extend LITERATE_TUTORIALS and the `pages` tree
# below; the charter's §10 tracker is the source of truth for what is wired yet. There is no
# on-site explanation quadrant: the "why" links out to the normative spec/ and the two papers
# (charter D1/D2 dropped 2026-07-21; see §8).
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
    ("tutorials/introductory.jl", "tutorials"),         # charter A1
    ("tutorials/advanced.jl", "tutorials"),             # charter A2
    ("tutorials/expert.jl", "tutorials"),               # charter A3
    ("deep_dives/serialization.jl", "deep_dives"),      # charter A4
    ("deep_dives/composition.jl", "deep_dives"),        # charter A5
    ("case_studies/marginal_scientist.jl", "case_studies"),  # charter B1 (flagship)
    ("case_studies/inlicensing_value.jl", "case_studies"),   # charter B2
    ("case_studies/kill_a_program.jl", "case_studies"),      # charter B3
]

for (src, outsub) in LITERATE_TUTORIALS
    Literate.markdown(
        joinpath(LITERATE_DIR, src), joinpath(GEN_DIR, outsub);
        documenter = true, execute = false, credit = false,
    )
end

# ── Pin & order the offered themes ──────────────────────────────────────────────────────────
# Documenter ships six themes and offers all of them (the theme picker + the copied CSS both
# read the hardcoded `HTMLWriter.THEMES` vector — there is no `HTML(; themes=…)` kwarg in 1.x).
# Our brand layer (assets/rd-theme.css) repaints all three themes we keep; the darker catppuccin
# flavours are dropped. ORDER MATTERS: Documenter marks `THEMES[1]` as the primary (default,
# light-preference) theme and `THEMES[2]` as the primary-dark (dark-OS-preference) fallback
# (HTMLWriter.jl ~L1137). We keep the BRANDED warm `documenter-light` as the default — it coheres
# with the whole warm-neutral identity system (spec/design_system.html; the teal/amber/rose
# figure hues are tuned against warm paper, not catppuccin's cool blue-grey) — with
# `documenter-dark` as the dark fallback and `catppuccin-latte` offered as an alternative in the
# picker. We OVERWRITE the vector in that exact order (a plain `filter!` preserves the stock
# light→dark→latte order, which happens to give the same default, but being explicit documents
# the intent). THEMES is a mutable Vector shared by every theme code path, so assigning into it
# in place also drops the unused catppuccin CSS from the build.
let want = ["documenter-light", "documenter-dark", "catppuccin-latte"]
    empty!(Documenter.HTMLWriter.THEMES)
    append!(Documenter.HTMLWriter.THEMES, want)
end

# ── Site ────────────────────────────────────────────────────────────────────────────────
makedocs(;
    sitename = "ReactiveDynamics.jl",
    format = Documenter.HTML(;
        prettyurls = get(ENV, "CI", "false") == "true", edit_link = "main",
        # House brand overrides (International Typographic Style — teal accent, off-white
        # nav, semantic figure hues), a thin layer over the two stock themes. See
        # docs/src/assets/rd-theme.css for what it repaints and why. The sidebar logo flips
        # by theme — `assets/logo.svg` (ink mark, light off-white nav) and
        # `assets/logo-dark.svg` (reversed-out, dark nav) — both the "firing glyph" from
        # the Claude Design identity board (spec/design_system.html).
        assets = [
            "assets/rd-theme.css",
            # Tints the ".jl" extension of the sidebar wordmark teal — Documenter
            # emits the sitename as a bare text node, so a tiny script splits off
            # the ".jl" suffix into a `.rd-jl` span that rd-theme.css paints.
            "assets/rd-logo.js",
            # SVG favicon — `assets` infers class from extension and only knows css/js,
            # so an .svg icon must be passed as an explicit :ico-class HTMLAsset.
            asset("assets/favicon.svg"; class = :ico, islocal = true),
        ],
    ),
    modules = [ReactiveDynamics],
    warnonly = true,   # scaffold stage: don't fail on cross-references to pages not yet authored
    pages = [
        "Home" => "index.md",
        "Tutorials" => [
            "Introductory: your first model" => "tutorials/introductory.md",        # charter A1
            "Advanced: structured-token portfolios" => "tutorials/advanced.md",      # charter A2
            "Expert: coupled heterogeneous systems" => "tutorials/expert.md",        # charter A3
            "Deep dive: serialization" => "deep_dives/serialization.md",             # charter A4
            "Deep dive: composition" => "deep_dives/composition.md",                 # charter A5
        ],
        "Case studies" => [                                      # charter B1–B3
            "Marginal eNPV of the Nth scientist" => "case_studies/marginal_scientist.md",
            "In-licensing asset value" => "case_studies/inlicensing_value.md",
            "When to kill a program" => "case_studies/kill_a_program.md",
        ],
        "Reference" => [                                         # charter C1–C2
            "Authoring" => "reference/authoring.md",
            "Structured tokens" => "reference/structured_tokens.md",
            "Rules & actions" => "reference/rules_actions.md",
            "Construction & simulation" => "reference/construction_simulation.md",
            "Composition" => "reference/composition.md",
            "Serialization" => "reference/serialization.md",
            "JSON model schema" => "reference/json_schema.md",
            "Analysis & visualization" => "reference/analysis_viz.md",
            "AlgebraicAgents coupling" => "reference/aa_coupling.md",
        ],
        # No "Explanation" nav: charter D1/D2 dropped 2026-07-21 — the "why" links out to the
        # normative spec/CONTRACT_DRAFT.md + spec/adr/ and the two companion papers (charter §8).
    ],
)

deploydocs(; repo = "github.com/Merck/ReactiveDynamics.jl.git", push_preview = true)
