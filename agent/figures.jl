#!/usr/bin/env julia
#
# figures.jl — the structural diagrams on the overview and copilot pages.
#
#   julia --project=. agent/figures.jl
#
# A study's figures carry measurements and belong to `scripts/chartkit.jl`. These carry structure
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

# ── figure 2: the same sequence, for the copilot page ─────────────────────────────────────────────

# `pipeline.svg` names the skill and the file each phase leaves behind, which is what a reader
# evaluating the layer needs and exactly what a reader evaluating the OFFER does not. Same four steps,
# same refusal branch, stated in the language of the decision instead of the language of the tool.
const PROMISE = [
    ("Describe the decision", "A conversation: the choice, what drives it, how the parts interact. If arithmetic answers it, we stop here."),
    ("We build your operation", "A running model of how it behaves — shared capacity, delays, feedback, decisions taken mid-course, cash."),
    ("We try to break it", "Nine checks ask whether the model does what you said the business does. A failure blocks the answer."),
    ("A number you can defend", "A central estimate with the evidence behind it — and a seed that reproduces it exactly."),
]

function promise_svg()
    title = "From a question to a number you can take to a board"
    subtitle = "Four steps, and the third one is allowed to stop the other three."
    head = CK.head_height(title, subtitle)

    gap = 16.0
    boxw = (W - 2PAD - 3gap) / 4
    boxh = 120.0
    ytop = head + 20.0
    parts = [CK.head_block(title, subtitle)]

    for (i, (name, detail)) in enumerate(PROMISE)
        x = PAD + (i - 1) * (boxw + gap)
        last = i == length(PROMISE)
        push!(parts, CK.rect(x, ytop, boxw, boxh; fill = CK.PAPER, stroke = CK.HAIR, sw = 1.0))
        push!(parts, CK.rect(x, ytop, boxw, 2.5; fill = last ? CK.INK : CK.ACCENT))
        push!(
            parts, CK.txt(
                x + 13, ytop + 25, "0" * string(i);
                size = 9.2, family = CK.MONO, fill = last ? CK.INK : CK.ACCENT_D, weight = 600
            )
        )
        for (j, l) in enumerate(CK.wrap_lines(name, boxw - 26, 11.4; maxlines = 2))
            push!(parts, CK.txt(x + 13, ytop + 44 + 14 * (j - 1), l; size = 11.4, weight = 600))
        end
        for (j, l) in enumerate(CK.wrap_lines(detail, boxw - 26, 8.8; maxlines = 4))
            push!(parts, CK.txt(x + 13, ytop + 74 + 10.5 * (j - 1), l; size = 8.8, fill = CK.GREY))
        end
        last || push!(parts, arrow_r(x + boxw + 2, x + boxw + gap - 2, ytop + boxh / 2))
    end

    xd = PAD + 2 * (boxw + gap) + boxw / 2
    ybr = ytop + boxh + 30
    push!(parts, arrow_d(xd, ytop + boxh + 3, ybr, stroke = CK.ROSE))
    bx, bw = PAD + boxw + gap, 3boxw + 2gap
    push!(parts, CK.rect(bx, ybr, bw, 44; fill = "#FDF3F4", stroke = CK.ROSE, sw = 1.0))
    push!(
        parts, CK.txt(
            bx + 14, ybr + 19, "If the model does not behave like your business, you get the diagnosis, not a number.";
            size = 10.2, weight = 600
        )
    )
    push!(
        parts, CK.txt(
            bx + 14, ybr + 34, "It names the check that failed and what it means: fix the model, or revisit the brief.";
            size = 9.2, fill = CK.GREY
        )
    )

    return CK.shell(join(parts), ybr + 44 + 8; title = title, subtitle = subtitle)
end

# ── figure 3: the integration surface ─────────────────────────────────────────────────────────────

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

# ── figure 4: what the copilot is made of ─────────────────────────────────────────────────────────

# component family · what it is for, in one clause · its members · status, with the tracker id. The
# division is the charter's §5: skills carry procedure, subagents buy context isolation, MCP buys reach
# and a persistent session, hooks buy enforcement a prompt cannot. Only the first ships today, and the
# figure has to say so in the same breath as it draws the rest.
const COMPONENTS = [
    (
        "Skills", "the procedure: what to ask, what to check",
        ["rd-study  (orchestrator)", "rd-frame · rd-author", "rd-doctor · rd-report"], "SHIPPED", true,
    ),
    (
        "Subagents", "isolated context for the heavy lifting",
        ["rd-modeler  authors the model", "rd-critic  audits the prose", "rd-calibrator  (phase 2)"],
        "PLANNED · W9", false,
    ),
    (
        "MCP server rd", "reach beyond Bash, a warm Julia session",
        ["rd_validate · rd_run", "rd_ensemble · rd_diagnose", "rd_effect · rd_draw", "+ the spec as resources"],
        "PLANNED · W8", false,
    ),
    (
        "Hook + commands", "enforcement no prompt can promise",
        ["validate on each *.rdj.json", "/rd:study · /rd:frame · …", "plugin + marketplace"],
        "PLANNED · W10, W12", false,
    ),
]

function architecture_svg()
    title = "ReactiveDynamics.jl Modeling Copilot — what it is made of"
    subtitle = "A thin layer of procedure over an engine that checks everything. The shipped path is the portable one; the planned pieces add trust (a critic, a validate hook) and reach (MCP, plugin)."
    head = CK.head_height(title, subtitle)
    parts = [CK.head_block(title, subtitle)]

    # the person and the host they talk through
    y = head + 16.0
    push!(parts, CK.rect(PAD, y, W - 2PAD, 40; fill = "#FAF8F3", stroke = CK.HAIR, sw = 1.0))
    push!(parts, CK.txt(PAD + 16, y + 17, "You — an analyst, a strategy lead, a VP"; size = 11.2, weight = 600))
    push!(
        parts, CK.txt(
            PAD + 16, y + 31, "in plain business language: the decision, what drives it, how things interact over time";
            size = 8.9, fill = CK.GREY
        )
    )
    push!(
        parts, CK.txt(
            W - PAD - 16, y + 17, "Claude Code today · Desktop / any MCP host next";
            size = 8.6, family = CK.MONO, fill = CK.ACCENT_D, anchor = "end"
        )
    )

    # the plugin: four component families side by side
    yp = y + 40 + 24
    push!(parts, arrow_d(W / 2, y + 43, yp))
    gap = 10.0
    colw = (W - 2PAD - 3gap) / 4
    ph = 132.0
    push!(parts, CK.txt(PAD, yp - 6, "THE COPILOT"; size = 8.4, family = CK.MONO, fill = CK.GREY_LT))
    for (i, (name, role, members, status, shipped)) in enumerate(COMPONENTS)
        x = PAD + (i - 1) * (colw + gap)
        push!(
            parts, CK.rect(
                x, yp, colw, ph; fill = shipped ? "#F7FBFB" : CK.PAPER,
                stroke = shipped ? CK.ACCENT : CK.HAIR, sw = shipped ? 1.4 : 1.0
            )
        )
        push!(parts, CK.txt(x + 11, yp + 20, name; size = 11.0, weight = 600))
        for (j, l) in enumerate(CK.wrap_lines(role, colw - 22, 8.6; maxlines = 2))
            push!(parts, CK.txt(x + 11, yp + 35 + 10.5 * (j - 1), l; size = 8.6, fill = CK.GREY))
        end
        for (j, l) in enumerate(members)
            push!(parts, CK.txt(x + 11, yp + 68 + 12.5 * (j - 1), l; size = 8.0, family = CK.MONO, fill = CK.INK))
        end
        push!(parts, chip(x + 11, yp + ph - 22, status; shipped = shipped))
    end

    # knowledge the components read at the moment of need
    yk = yp + ph + 10
    push!(parts, CK.rect(PAD, yk, W - 2PAD, 30; fill = CK.INSET))
    push!(parts, CK.txt(PAD + 12, yk + 19, "Knowledge"; size = 9.4, weight = 600))
    push!(
        parts, CK.txt(
            PAD + 82, yk + 19,
            "normative spec + ADRs, cited by section · schema reference GENERATED from the engine, drift-gated in CI";
            size = 8.7, fill = CK.GREY
        )
    )

    # one process boundary, then the engine
    ycli = yk + 30 + 24
    push!(parts, arrow_d(W / 2, yk + 33, ycli))
    clix, cliw = PAD + 28, W - 2PAD - 56
    push!(parts, CK.rect(clix, ycli, cliw, 46; fill = "#F7FBFB", stroke = CK.ACCENT, sw = 1.4))
    push!(parts, CK.txt(clix + 14, ycli + 19, "rdx"; size = 11.4, family = CK.MONO, weight = 600, fill = CK.ACCENT_D))
    push!(
        parts, CK.txt(
            clix + 48, ycli + 19, "validate · run · diagnose · ensemble · effect · draw";
            size = 9.4, family = CK.MONO
        )
    )
    push!(
        parts, CK.txt(
            clix + 14, ycli + 35, "one CLI under every path above · JSON out · exit 0 ok / 1 findings / 2 error";
            size = 8.8, fill = CK.GREY
        )
    )
    push!(parts, chip(clix + cliw - CK.tw("SHIPPED", 8.4; mono = true) - 26, ycli + 10, "SHIPPED"; shipped = true))

    yeng = ycli + 46 + 22
    push!(parts, arrow_d(W / 2, ycli + 49, yeng))
    push!(parts, CK.rect(PAD, yeng, W - 2PAD, 58; fill = CK.INK))
    push!(parts, CK.txt(PAD + 16, yeng + 23, "ReactiveDynamics.jl"; size = 11.6, weight = 600, fill = "#F4F4F2"))
    push!(
        parts, CK.txt(
            PAD + 16, yeng + 41,
            "model = eval-free JSON document · pure validate → diagnostics with paths · (model, seed) fixes a run";
            size = 9.0, fill = "#C7C7C4"
        )
    )
    push!(
        parts, CK.txt(
            W - PAD - 16, yeng + 23, "on AlgebraicAgents.jl";
            size = 8.6, family = CK.MONO, fill = CK.ACCENT_LT, anchor = "end"
        )
    )

    return CK.shell(join(parts), yeng + 58 + 8; title = title, subtitle = subtitle)
end

# ── figure 5: who does what, and what passes between them ─────────────────────────────────────────

# phase · skill · its job, in one clause · the file it leaves · the rdx verbs it calls · its gate. The
# architecture figure says what the copilot is made of; this one says how the parts hand work to each
# other, which is the question an engineer asks next.
const ROLES = [
    (
        "FRAME", "rd-frame", "interviews the owner; does the qualifying arithmetic",
        "BRIEF.md", "no rdx: plain arithmetic", "something binds, or stop here",
    ),
    (
        "AUTHOR", "rd-author", "writes one JSON model per option, calibrates it",
        "model.rdj.json + arm_*", "validate · run", "zero validate errors; each arm differs only in its lever",
    ),
    (
        "DIAGNOSE", "rd-doctor", "tries to break the model before anyone trusts it",
        "DOCTOR.md", "diagnose", "nine checks pass, or each failure explained",
    ),
    (
        "REPORT", "rd-report", "measures the options against each other",
        "FINDINGS.md → REPORT.html", "effect · ensemble · draw", "every number carries a spread and a seed",
    ),
]

# the planned pieces, each attached to the phase it serves: (phase, name, what it does, W-id)
const ATTACH = [
    ("AUTHOR", "rd-modeler", "a subagent: runs write → validate → repair in its own context, hands back only the model", "W9"),
    ("AUTHOR", "validate hook", "re-checks every *.rdj.json write; no prompt can skip it", "W10"),
    ("REPORT", "rd-critic", "a subagent with fresh context: the model against the brief, the prose against the numbers", "W9"),
]

function roles_svg()
    title = "Inside the copilot: who does what, and what passes between them"
    subtitle = "One orchestrator, four phase skills that each leave a file, and one CLI underneath. Work passes between the phases only as files on disk, which is what makes each one re-runnable and auditable."
    head = CK.head_height(title, subtitle)
    parts = [CK.head_block(title, subtitle)]

    # the orchestrator, and the door a user comes in through
    y = head + 16.0
    push!(parts, CK.rect(PAD, y, W - 2PAD, 44; fill = "#F7FBFB", stroke = CK.ACCENT, sw = 1.4))
    push!(parts, CK.txt(PAD + 14, y + 18, "rd-study"; size = 11.4, family = CK.MONO, weight = 600, fill = CK.ACCENT_D))
    push!(parts, CK.txt(PAD + 86, y + 18, "the orchestrator"; size = 11.0, weight = 600))
    push!(
        parts, CK.txt(
            PAD + 14, y + 34,
            "reads the study directory, resumes at the phase it implies, enforces every gate, never skips a phase";
            size = 8.8, fill = CK.GREY
        )
    )
    cw = CK.tw("/rd:study command · W12", 8.4; mono = true) + 12
    push!(parts, chip(W - PAD - 12 - cw, y + 8, "/rd:study command · W12"))

    # the four phase skills
    gap = 12.0
    colw = (W - 2PAD - 3gap) / 4
    yc = y + 44 + 26
    ch = 152.0
    for i in 1:4
        push!(parts, arrow_d(PAD + (i - 0.5) * (colw + gap) - gap / 2, y + 47, yc))
    end
    for (i, (phase, skill, job, writes, calls, gate)) in enumerate(ROLES)
        x = PAD + (i - 1) * (colw + gap)
        push!(parts, CK.rect(x, yc, colw, ch; fill = CK.PAPER, stroke = CK.HAIR, sw = 1.0))
        push!(parts, CK.rect(x, yc, colw, 2.5; fill = i == 4 ? CK.INK : CK.ACCENT))
        push!(parts, CK.txt(x + 11, yc + 20, phase; size = 11.6, weight = 600))
        push!(parts, CK.txt(x + colw - 11, yc + 19, skill; size = 8.6, family = CK.MONO, fill = CK.ACCENT_D, anchor = "end"))
        for (j, l) in enumerate(CK.wrap_lines(job, colw - 22, 8.8; maxlines = 2))
            push!(parts, CK.txt(x + 11, yc + 36 + 10.5 * (j - 1), l; size = 8.8, fill = CK.INK))
        end
        rows = (("writes", writes, CK.AMBER, true), ("calls", calls, CK.ACCENT_D, true), ("gate", gate, CK.GREY, false))
        yy = yc + 68
        for (lab, val, col, mono) in rows
            push!(parts, CK.txt(x + 11, yy, uppercase(lab); size = 7.4, family = CK.MONO, fill = CK.GREY_LT))
            for (j, l) in enumerate(CK.wrap_lines(val, colw - 22, mono ? 9.6 : 8.4; maxlines = 2))
                push!(parts, CK.txt(x + 11, yy + 11 + 10.2 * (j - 1), l; size = 8.4, family = mono ? CK.MONO : CK.SANS, fill = col))
            end
            yy += lab == "gate" ? 0 : 27
        end
        i == 4 || push!(parts, arrow_r(x + colw + 1, x + colw + gap - 1, yc + 26))
    end

    # the one backward edge: a failed check returns the work to AUTHOR, never forward to REPORT
    xa = PAD + (colw + gap) + colw / 2
    xd = PAD + 2 * (colw + gap) + colw / 2
    yl = yc + ch + 16
    push!(
        parts, "<path d=\"M $xd $(yc + ch) L $xd $yl L $xa $yl L $xa $(yc + ch + 6)\" fill=\"none\" stroke=\"$(CK.ROSE)\" stroke-width=\"1.4\"/>" *
            "<path d=\"M $(xa - 4) $(yc + ch + 7) L $xa $(yc + ch + 1) L $(xa + 4) $(yc + ch + 7) Z\" fill=\"$(CK.ROSE)\"/>"
    )
    push!(
        parts, CK.txt(
            (xa + xd) / 2, yl + 12, "a failed check sends the model back to AUTHOR, with the check that failed";
            size = 8.6, fill = CK.ROSE, anchor = "middle"
        )
    )

    # the planned attachments, each under the phase it serves
    ya = yl + 26
    push!(parts, CK.txt(PAD, ya - 5, "PLANNED, ATTACHED TO THE PHASE THEY SERVE"; size = 8.0, family = CK.MONO, fill = CK.GREY_LT))
    ah = 80.0
    # AUTHOR carries two attachments, so they span its column and the one to its left
    slots = [(PAD, 2colw + gap), (PAD + 2 * (colw + gap), colw), (PAD + 3 * (colw + gap), colw)]
    for (k, (phase, name, what, wid)) in enumerate(ATTACH)
        x, w = slots[k]
        push!(parts, "<rect x=\"$x\" y=\"$ya\" width=\"$w\" height=\"$ah\" fill=\"$(CK.PAPER)\" stroke=\"$(CK.GREY_LT)\" stroke-width=\"1\" stroke-dasharray=\"4 3\"/>")
        push!(parts, CK.txt(x + 11, ya + 17, name; size = 9.6, family = CK.MONO, weight = 600))
        push!(parts, CK.txt(x + 11, ya + 30, "SERVES " * phase; size = 7.4, family = CK.MONO, fill = CK.GREY_LT))
        push!(parts, chip(x + w - 11 - CK.tw(wid, 8.4; mono = true) - 12, ya + 6, wid))
        for (j, l) in enumerate(CK.wrap_lines(what, w - 22, 8.5; maxlines = 3))
            push!(parts, CK.txt(x + 11, ya + 45 + 10.5 * (j - 1), l; size = 8.5, fill = CK.GREY))
        end
    end

    # two doors onto the same verbs
    yr = ya + ah + 24
    push!(parts, arrow_d(W / 2, ya + ah + 3, yr))
    rw = W - 2PAD
    push!(parts, CK.rect(PAD, yr, rw, 62; fill = "#F7FBFB", stroke = CK.ACCENT, sw = 1.4))
    push!(parts, CK.txt(PAD + 14, yr + 19, "rdx"; size = 11.4, family = CK.MONO, weight = 600, fill = CK.ACCENT_D))
    push!(parts, CK.txt(PAD + 48, yr + 19, "validate · run · diagnose · ensemble · effect · draw"; size = 9.4, family = CK.MONO))
    push!(parts, chip(W - PAD - 14 - CK.tw("SHIPPED", 8.4; mono = true) - 12, yr + 8, "SHIPPED"; shipped = true))
    push!(
        parts, CK.txt(
            PAD + 14, yr + 36, "Today every skill reaches it through Bash: JSON on stdout, exit 0 ok / 1 findings / 2 error.";
            size = 8.8, fill = CK.GREY
        )
    )
    push!(
        parts, CK.txt(
            PAD + 14, yr + 51, "Next, the same six verbs as MCP tools (rd_validate, rd_run, …) over a warm Julia session, for hosts with no terminal.";
            size = 8.8, fill = CK.GREY
        )
    )
    push!(parts, chip(W - PAD - 14 - CK.tw("MCP · W8", 8.4; mono = true) - 12, yr + 40, "MCP · W8"))

    yeng = yr + 62 + 22
    push!(parts, arrow_d(W / 2, yr + 65, yeng))
    push!(parts, CK.rect(PAD, yeng, W - 2PAD, 44; fill = CK.INK))
    push!(parts, CK.txt(PAD + 16, yeng + 18, "ReactiveDynamics.jl"; size = 11.2, weight = 600, fill = "#F4F4F2"))
    push!(
        parts, CK.txt(
            PAD + 16, yeng + 34, "validate is pure and returns every problem with its path · (model, seed) fixes a run, byte for byte";
            size = 8.8, fill = "#C7C7C4"
        )
    )
    push!(
        parts, CK.txt(
            W - PAD - 16, yeng + 18, "on AlgebraicAgents.jl";
            size = 8.4, family = CK.MONO, fill = CK.ACCENT_LT, anchor = "end"
        )
    )

    return CK.shell(join(parts), yeng + 44 + 8; title = title, subtitle = subtitle)
end

# ── build ─────────────────────────────────────────────────────────────────────────────────────────

for (name, svg) in (
        "pipeline" => pipeline_svg(), "promise" => promise_svg(), "surface" => surface_svg(),
        "architecture" => architecture_svg(), "roles" => roles_svg(),
    )
    path = joinpath(FIGDIR, name * ".svg")
    write(path, svg)
    println("figures: wrote $(relpath(path, dirname(HERE))) — $(round(sizeof(svg) / 1024, digits = 1)) KiB")
end
