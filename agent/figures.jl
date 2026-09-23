#!/usr/bin/env julia
#
# figures.jl — the two structural diagrams on the overview page.
#
#   julia --project=. agent/figures.jl
#
# A study's figures carry measurements and belong to `scripts/chartkit.jl`. These two carry structure
# instead — the phase sequence, and the surface a user reaches it through — so there is no number in
# them to drift. What they CAN drift against is the charter's tracker, which is why they are generated
# from the tables below rather than drawn by hand: every shipped/planned mark sits once, next to the
# W-id it cites, and a status change is a one-line edit followed by a re-run.
#
# Palette, type stacks and SVG primitives come from `scripts/chartkit.jl`, so a diagram is
# indistinguishable from a chart on the same page and neither pulls in a plotting stack.

const HERE = @__DIR__
const FIGDIR = joinpath(HERE, "figures")

include(joinpath(HERE, "scripts", "chartkit.jl"))
const CK = ChartKit

const W = CK.WIDTH
const PAD = CK.PAD

mkpath(FIGDIR)

# ── shared pieces ─────────────────────────────────────────────────────────────────────────────────

# A right-pointing arrow on the rail between two cards, and a down-pointing one into a lower band.
arrow_r(x1, x2, y) = CK.line(x1, y, x2 - 5, y; stroke = CK.HAIR, sw = 1.4) *
    "<path d=\"M $(x2 - 6) $(y - 4) L $x2 $y L $(x2 - 6) $(y + 4) Z\" fill=\"$(CK.HAIR)\"/>"

arrow_d(x, y1, y2; stroke = CK.HAIR) = CK.line(x, y1, x, y2 - 5; stroke = stroke, sw = 1.4) *
    "<path d=\"M $(x - 4) $(y2 - 6) L $x $y2 L $(x + 4) $(y2 - 6) Z\" fill=\"$stroke\"/>"

# A status chip. Shipped is filled because it is the claim a reader should be able to check today;
# planned is outlined and carries its tracker id, so the diagram never implies more than exists.
function chip(x, y, label; shipped = false)
    w = CK.tw(label, 8.4; mono = true) + 12
    box = shipped ?
        CK.rect(x, y, w, 15; fill = "#E8F3F2", stroke = "#BFDEDB", sw = 1.0) :
        CK.rect(x, y, w, 15; fill = CK.PAPER, stroke = CK.HAIR, sw = 1.0)
    return box * CK.txt(
        x + 6, y + 11, label;
        size = 8.4, family = CK.MONO, fill = shipped ? CK.ACCENT_D : CK.GREY
    )
end

# ── figure 1: the phase sequence ──────────────────────────────────────────────────────────────────

# phase · the skill that runs it · the artifact it leaves on disk · the gate, in one clause
const PHASES = [
    ("FRAME", "rd-frame", "BRIEF.md", "the qualifying test: does something bind"),
    ("AUTHOR", "rd-author", "model.rdj.json", "validate returns an empty vector"),
    ("DIAGNOSE", "rd-doctor", "DOCTOR.md", "nine plausibility checks, or each failure stated"),
    ("REPORT", "rd-report", "FINDINGS.md", "every number carries a spread and a seed"),
]

function pipeline_svg()
    title = "One question, four artifacts, and a gate that can stop the sequence"
    subtitle = "Each phase ends in a file on disk. A phase may not be skipped, and a returning session re-enters at whichever phase the study directory already implies."
    head = CK.head_height(title, subtitle)

    boxw, boxh, gap = 164.0, 78.0, 18.0
    x0 = PAD + (W - 2PAD - (4boxw + 3gap)) / 2
    ytop = head + 20.0
    parts = [CK.head_block(title, subtitle)]

    for (i, (phase, skill, artifact, gate)) in enumerate(PHASES)
        x = x0 + (i - 1) * (boxw + gap)
        last = i == length(PHASES)
        push!(parts, CK.rect(x, ytop, boxw, boxh; fill = CK.PAPER, stroke = CK.HAIR, sw = 1.0))
        push!(parts, CK.rect(x, ytop, boxw, 2.5; fill = last ? CK.INK : CK.ACCENT))
        push!(parts, CK.txt(x + 12, ytop + 23, phase; size = 12.5, weight = 600))
        push!(
            parts, CK.txt(
                x + boxw - 12, ytop + 22, skill;
                size = 8.6, family = CK.MONO, fill = CK.ACCENT_D, anchor = "end"
            )
        )
        push!(
            parts, CK.txt(
                x + 12, ytop + 39, artifact;
                size = 9.6, family = CK.MONO, fill = CK.AMBER
            )
        )
        for (j, l) in enumerate(CK.wrap_lines(gate, boxw - 24, 8.8; maxlines = 3))
            push!(parts, CK.txt(x + 12, ytop + 55 + 10.5 * (j - 1), l; size = 8.8, fill = CK.GREY))
        end
        last || push!(parts, arrow_r(x + boxw + 3, x + boxw + gap - 3, ytop + boxh / 2))
    end

    # The branch is the whole point of drawing this rather than tabulating it: DIAGNOSE is allowed to
    # refuse. A checklist that can only pass is decoration.
    xd = x0 + 2 * (boxw + gap) + boxw / 2
    ybr = ytop + boxh + 30
    push!(parts, arrow_d(xd, ytop + boxh + 3, ybr, stroke = CK.ROSE))
    bx, bw = x0 + boxw + gap, 3boxw + 2gap
    push!(parts, CK.rect(bx, ybr, bw, 44; fill = "#FDF3F4", stroke = CK.ROSE, sw = 1.0))
    push!(
        parts, CK.txt(
            bx + 14, ybr + 19, "A red check blocks the headline number — it is not a footnote.";
            size = 10.2, weight = 600
        )
    )
    push!(
        parts, CK.txt(
            bx + 14, ybr + 34,
            "Nothing contended ⇒ recalibrate, or say a spreadsheet is the right tool and stop.";
            size = 9.2, fill = CK.GREY
        )
    )

    return CK.shell(join(parts), ybr + 44 + 8; title = title, subtitle = subtitle)
end

# ── figure 2: the integration surface ─────────────────────────────────────────────────────────────

# how a user reaches the layer · what that path is made of · status, with the tracker id
const REACH = [
    (
        "Claude Code", ["/rd:study, /rd:frame, …", "PostToolUse validate hook", "plugin + marketplace"],
        "PLANNED · W10, W12", false,
    ),
    (
        "Any MCP host", ["rd — local stdio server", "6 tools + spec resources", "persistent Julia session"],
        "PLANNED · W8", false,
    ),
    (
        "Any harness, or by hand", ["five skills, progressive", "Bash over the CLI below", "portable Markdown + Julia"],
        "SHIPPED", true,
    ),
]

function surface_svg()
    title = "The integration surface: one CLI, three ways in, nothing remote"
    subtitle = "The valuable layer is the portable one. Skills and scripts run in any harness; MCP and the plugin add reach, not capability."
    head = CK.head_height(title, subtitle)

    parts = [CK.head_block(title, subtitle)]
    y = head + 22.0
    push!(
        parts, CK.txt(
            PAD, y - 7, "HOW A USER REACHES IT";
            size = 8.4, family = CK.MONO, fill = CK.GREY_LT
        )
    )

    colw = (W - 2PAD) / 3
    bandh = 104.0
    push!(parts, CK.rect(PAD, y, W - 2PAD, bandh; fill = CK.PAPER, stroke = CK.HAIR, sw = 1.0))
    for (i, (name, lines, status, shipped)) in enumerate(REACH)
        x = PAD + (i - 1) * colw
        i == 1 || push!(parts, CK.line(x, y + 10, x, y + bandh - 10; stroke = CK.HAIR))
        push!(parts, CK.txt(x + 16, y + 24, name; size = 11.5, weight = 600))
        for (j, l) in enumerate(lines)
            push!(
                parts, CK.txt(
                    x + 16, y + 42 + 13 * (j - 1), l;
                    size = 8.9, fill = CK.GREY, family = startswith(l, "/rd") || startswith(l, "rd ") ? CK.MONO : CK.SANS
                )
            )
        end
        push!(parts, chip(x + 16, y + bandh - 24, status; shipped = shipped))
    end

    # Everything above converges on one process boundary, and that is the design: the MCP server and
    # the plugin are transports over the same commands a shell already runs.
    ycli = y + bandh + 26
    for i in 1:3
        push!(parts, arrow_d(PAD + (i - 0.5) * colw, y + bandh + 3, ycli))
    end
    clix, cliw = PAD + 28, W - 2PAD - 56
    push!(parts, CK.rect(clix, ycli, cliw, 66; fill = "#F7FBFB", stroke = CK.ACCENT, sw = 1.4))
    push!(
        parts, CK.txt(
            clix + 16, ycli + 22, "agent/scripts/rdx.jl";
            size = 11.6, family = CK.MONO, weight = 600, fill = CK.ACCENT_D
        )
    )
    push!(parts, chip(clix + cliw - CK.tw("SHIPPED", 8.4; mono = true) - 28, ycli + 11, "SHIPPED"; shipped = true))
    push!(
        parts, CK.txt(
            clix + 16, ycli + 40, "validate · run · diagnose · ensemble · effect · draw";
            size = 9.6, family = CK.MONO, fill = CK.INK
        )
    )
    push!(
        parts, CK.txt(
            clix + 16, ycli + 55, "JSON on stdout · exit 0 ok / 1 findings / 2 error · --host=FILE is the one code seam";
            size = 8.9, fill = CK.GREY
        )
    )

    yeng = ycli + 66 + 26
    push!(parts, arrow_d(W / 2, ycli + 69, yeng))
    push!(parts, CK.rect(PAD, yeng, W - 2PAD, 58; fill = CK.INK))
    push!(parts, CK.txt(PAD + 16, yeng + 24, "ReactiveDynamics.jl"; size = 11.6, weight = 600, fill = "#F4F4F2"))
    push!(
        parts, CK.txt(
            PAD + 16, yeng + 42,
            "pure validate → Vector{Diagnostic} · eval-free JSON IR, closed whitelist · (model, seed) determines a run";
            size = 9.0, fill = "#C7C7C4"
        )
    )
    push!(
        parts, CK.txt(
            W - PAD - 16, yeng + 24, "the trust boundary";
            size = 8.6, family = CK.MONO, fill = CK.ACCENT_LT, anchor = "end"
        )
    )

    return CK.shell(join(parts), yeng + 58 + 8; title = title, subtitle = subtitle)
end

# ── build ─────────────────────────────────────────────────────────────────────────────────────────

for (name, svg) in ("pipeline" => pipeline_svg(), "surface" => surface_svg())
    path = joinpath(FIGDIR, name * ".svg")
    write(path, svg)
    println("figures: wrote $(relpath(path, dirname(HERE))) — $(round(sizeof(svg) / 1024, digits = 1)) KiB")
end
