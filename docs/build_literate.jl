# ════════════════════════════════════════════════════════════════════════════════════════
# build_literate.jl — render ONE Literate tutorial to a self-contained HTML page, executed.
#
#   docs/literate/<path>.jl  →  docs/preview/<name>.html   (executed end-to-end)
#
# This is the fast authoring loop and the exemplar's acceptance check: it runs the tutorial
# top to bottom (so a code error is a hard failure) and captures the printed output into a
# single self-contained HTML file — no site build, no server. It is the same executed-Literate
# → markdown → HTML recipe as demo/wires_viz_tour/build.jl, generalized.
#
#   julia --project=docs docs/build_literate.jl [literate/tutorials/introductory.jl]
#
# The full site (all pages, cross-links, autodocs) is built by make.jl instead.
# ════════════════════════════════════════════════════════════════════════════════════════

using Literate
import Markdown
import Base64

const HERE = @__DIR__
const DEFAULT_SRC = joinpath(HERE, "literate", "tutorials", "introductory.jl")
const OUT_DIR = joinpath(HERE, "preview")

src = isempty(ARGS) ? DEFAULT_SRC : (isabspath(ARGS[1]) ? ARGS[1] : joinpath(HERE, ARGS[1]))
name = first(splitext(basename(src)))
isdir(OUT_DIR) || mkpath(OUT_DIR)

banner(t) = (println(); println("="^80); println(t); println("="^80))

# Literate → markdown WITH executed outputs. `execute = true` runs every code block; a Julia
# error aborts the render, so a green build IS the "every block executes clean" acceptance
# criterion. CommonMarkFlavor keeps the markdown backend-agnostic for the stdlib renderer below.
banner("Executing $(basename(src)) → $(name).html")
mddir = mktempdir()
Literate.markdown(
    src, mddir; execute = true, flavor = Literate.CommonMarkFlavor(),
    credit = false, name = name,
)
escape_regex(s) = replace(s, r"([.\-\[\]()+*?^$\\|{}])" => s"\\\1")

# Inline each sidecar figure Literate wrote for a plot so the HTML is self-contained: SVGs
# are spliced in raw; PNG/JPEG are base64 data-URI'd. (Same trick as demo/wires_viz_tour/build.jl.)
# Wrapped in a function so the accumulating `html` reassignment has a clear local scope.
function inline_figures(html, dir)
    for f in readdir(dir)
        ext = lowercase(last(splitext(f)))
        ext in (".svg", ".png", ".jpg", ".jpeg") || continue
        tag = if ext == ".svg"
            svg = read(joinpath(dir, f), String)
            i = findfirst("<svg", svg)
            "<div class=\"figure\">" * (i === nothing ? svg : svg[first(i):end]) * "</div>"
        else
            mime = ext == ".png" ? "image/png" : "image/jpeg"
            data = Base64.base64encode(read(joinpath(dir, f)))
            "<div class=\"figure\"><img src=\"data:$(mime);base64,$(data)\" /></div>"
        end
        html = replace(html, Regex("<img src=\"" * escape_regex(f) * "\"[^>]*/>") => tag)
    end
    return html
end

md = read(joinpath(mddir, "$(name).md"), String)
body = inline_figures(Markdown.html(Markdown.parse(md)), mddir)

# A minimal self-contained HTML wrapper — readable typography, no external assets.
function wrap_html(title, body)
    css = """
    body { max-width: 860px; margin: 2rem auto; padding: 0 1.2rem;
           font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Helvetica, Arial, sans-serif;
           line-height: 1.55; color: #1a1a1a; }
    h1 { font-size: 2rem; border-bottom: 2px solid #eee; padding-bottom: .4rem; }
    h2 { font-size: 1.4rem; margin-top: 2.4rem; color: #16324f; }
    h3 { font-size: 1.1rem; margin-top: 1.6rem; color: #16324f; }
    pre { background: #f6f8fa; padding: .8rem 1rem; border-radius: 6px; overflow-x: auto;
          font-size: .85rem; line-height: 1.4; }
    code { font-family: "SF Mono", Menlo, Consolas, monospace; }
    p code, li code { background: #f0f4f8; padding: .05rem .3rem; border-radius: 4px; }
    blockquote { border-left: 3px solid #cbd5e0; margin: 1rem 0; padding: .2rem 1rem; color: #444; }
    hr { border: none; border-top: 1px solid #eee; margin: 2.5rem 0 1rem; }
    .figure { text-align: center; margin: 1.4rem 0; }
    .figure svg, .figure img { max-width: 100%; height: auto; }
    """
    return """<!DOCTYPE html>
    <html lang="en"><head><meta charset="utf-8">
    <meta name="viewport" content="width=device-width, initial-scale=1">
    <title>$title</title><style>$css</style></head><body>
    $body
    </body></html>
    """
end

out = joinpath(OUT_DIR, "$(name).html")
write(out, wrap_html("ReactiveDynamics.jl — $(name)", body))
banner("DONE")
println("  wrote $(out)  ($(filesize(out)) bytes)")
