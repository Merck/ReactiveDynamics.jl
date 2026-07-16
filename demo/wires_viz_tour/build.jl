# ════════════════════════════════════════════════════════════════════════════════════════
# build.jl — regenerate the wires_viz_tour HTML from its Literate source.
#
#   demo/wires_viz_tour/wires_literate.jl  (Literate script) → wires_literate.html
#
# The script shares wires_model.jl (the wired-hierarchy build + the two Graphviz renders + the
# data). Run:
#
#     julia --project=demo/wires_viz_tour demo/wires_viz_tour/build.jl
#
# First time (resolve the demo env — a few minutes):
#     julia --project=demo/wires_viz_tour -e 'using Pkg; Pkg.instantiate()'
#
# Long compiles buffer output; when driving this from a watchdog'd shell, launch in the background
# and poll the log. The render executes the tutorial end-to-end (real Graphviz SVGs via `dot`).
# ════════════════════════════════════════════════════════════════════════════════════════

using Pkg
Pkg.activate(@__DIR__)

import Literate
import Markdown

const HERE = @__DIR__
const LITERATE_SRC = joinpath(HERE, "wires_literate.jl")
const LITERATE_HTML = joinpath(HERE, "wires_literate.html")

banner(t) = (println(); println("="^80); println(t); println("="^80))

# ────────────────────────────────────────────────────────────────────────────────────────
# Literate script → markdown → self-contained HTML
#
# Literate.jl's native targets are markdown / notebook / script (no HTML target in this version), so
# we render Literate → markdown WITH executed outputs (`execute = true`), which writes the Plots
# chart and each diagram as sidecar `.svg` files referenced by `![](name-N.svg)`. We then convert the
# markdown to HTML with the stdlib `Markdown` (CommonMark is not in the env) and INLINE every sidecar
# SVG in place of its `<img>` tag, yielding ONE self-contained HTML file with the diagrams embedded.
# ────────────────────────────────────────────────────────────────────────────────────────
function build_literate()
    banner("Literate script → $(basename(LITERATE_HTML))")

    mddir = mktempdir()
    # Literate `cd`s into the output dir while executing, so `@__DIR__` inside the script would
    # resolve to `mddir`. Pin the real demo dir via ENV so the script's `include(wires_model.jl)`
    # finds it (the script reads `WIRES_VIZ_DIR`, falling back to `@__DIR__` when run directly).
    ENV["WIRES_VIZ_DIR"] = HERE
    Literate.markdown(
        LITERATE_SRC, mddir; execute = true, flavor = Literate.CommonMarkFlavor(),
        credit = true, name = "wires_literate"
    )
    mdfile = joinpath(mddir, "wires_literate.md")
    md = read(mdfile, String)

    body = Markdown.html(Markdown.parse(md))

    # Inline each sidecar SVG: replace <img src="name-N.svg" ...> with the SVG file's contents.
    for f in readdir(mddir)
        endswith(f, ".svg") || continue
        svg = read(joinpath(mddir, f), String)
        i = findfirst("<svg", svg)
        svg = i === nothing ? svg : svg[first(i):end]
        # Markdown.html emits <img src="f" alt="" /> ; splice the raw SVG in its place.
        body = replace(
            body, Regex("<img src=\"" * escape_regex(f) * "\"[^>]*/>") =>
                "<div class=\"diagram\">" * svg * "</div>"
        )
    end

    html = wrap_html("ReactiveDynamics — wires & ports (Literate)", body)
    write(LITERATE_HTML, html)
    return println("wrote $(LITERATE_HTML)  ($(filesize(LITERATE_HTML)) bytes)")
end

escape_regex(s) = replace(s, r"([.\-\[\]()+*?^$\\])" => s"\\\1")

# A minimal self-contained HTML wrapper (readable typography + light styling), no external assets.
function wrap_html(title, body)
    css = """
    body { max-width: 860px; margin: 2rem auto; padding: 0 1.2rem;
           font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Helvetica, Arial, sans-serif;
           line-height: 1.55; color: #1a1a1a; }
    h1 { font-size: 1.9rem; border-bottom: 2px solid #eee; padding-bottom: .4rem; }
    h2 { font-size: 1.35rem; margin-top: 2.2rem; color: #16324f; }
    pre { background: #f6f8fa; padding: .8rem 1rem; border-radius: 6px; overflow-x: auto;
          font-size: .85rem; line-height: 1.4; }
    code { font-family: "SF Mono", Menlo, Consolas, monospace; }
    table { border-collapse: collapse; margin: 1rem 0; font-size: .9rem; }
    th, td { border: 1px solid #ddd; padding: .35rem .7rem; text-align: right; }
    th { background: #f0f4f8; }
    .diagram { text-align: center; margin: 1.4rem 0; padding: 1rem; background: #fbfcfd;
               border: 1px solid #eef1f4; border-radius: 8px; }
    .diagram svg { max-width: 100%; height: auto; }
    img { max-width: 100%; }
    hr { border: none; border-top: 1px solid #eee; margin: 2.5rem 0 1rem; }
    """
    return """<!DOCTYPE html>
    <html lang="en"><head><meta charset="utf-8">
    <meta name="viewport" content="width=device-width, initial-scale=1">
    <title>$title</title>
    <style>$css</style>
    </head><body>
    $body
    </body></html>
    """
end

build_literate()

banner("DONE — HTML regenerated")
println("  $(basename(LITERATE_HTML))     $(filesize(LITERATE_HTML)) bytes")
