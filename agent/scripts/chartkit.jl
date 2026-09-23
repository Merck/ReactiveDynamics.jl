#!/usr/bin/env julia
#
# chartkit.jl — data-driven SVG figures for the self-contained session page.
#
#   julia --project=. agent/scripts/chartkit.jl SPEC.json [--out=fig.svg]   # --out=- writes stdout
#
# render_session.jl inlines its figures, so a figure has to be SVG text with no external anything. The
# alternatives were tried and rejected: a base64 PNG from a plotting stack bloats the page, pixelates
# on the retina screens this gets read on and prints badly, while hand-authored SVG cannot be
# regenerated when the numbers move — which would quietly void the artifact's one real claim, that
# every figure traces to a command that was actually run. Hence: vectors in, SVG out, over the same
# palette and font stacks as the page's CSS, so an embedded chart is indistinguishable from page
# furniture and carries no dependency on it.
#
# Four forms, because a study needs four: `bars_error` (the workhorse — an effect per lever, with the
# null ones drawn AS nulls), `waterfall` (a value bridge), `paired_dots` (two arms over one seed
# stream, the figure that earns the error bars) and `trace_ceiling` (a resource against its cap, with
# saturation shaded).
#
# Two ways in, one code path: call a builder for an SVG `String`, or hand the CLI a spec
# `{"chart": "bars_error", …}`. Field access is deliberately duck-typed (`fld`) so a JSON object and a
# NamedTuple travel the same road and the CLI is not a second implementation. A missing required key
# or an unknown chart type raises `ChartError` and exits 2 — an empty figure in an exec deck is worse
# than a crash, because nobody notices it.
#
# No randomness anywhere: dot jitter is derived from the sample index, so a spec renders byte-identical
# every time and the page can be diffed. No dependency beyond JSON and Printf.
#
# Spec keys per chart are listed above each builder and repeated in USAGE (`--help`).

module ChartKit

using JSON
using Printf

# ── palette and type, lifted verbatim from render_session.jl's CSS custom properties ───────────────

const INK = "#0A0A0A"
const PAPER = "#FFFFFF"
const INSET = "#F4F4F2"
const ACCENT = "#0A8A84"
const ACCENT_D = "#076A65"
const ACCENT_LT = "#3FB5AF"
const AMBER = "#B5740A"
const ROSE = "#E5243B"
const GREY = "#6A6A68"
const GREY_LT = "#9A9A98"
const HAIR = "#C7C7C4"

const SANS = "Inter, 'Helvetica Neue', Helvetica, Arial, system-ui, sans-serif"
const MONO = "'JetBrains Mono', 'SF Mono', ui-monospace, monospace"

const WIDTH = 760      # the page column is 50rem; wider gets scaled down in flow and blurs the type
const PAD = 18

# ── errors ────────────────────────────────────────────────────────────────────────────────────────

struct ChartError <: Exception
    msg::String
end

Base.showerror(io::IO, e::ChartError) = print(io, "chartkit: ", e.msg)

fail(msg) = throw(ChartError(msg))

# ── field access: JSON objects and NamedTuples read the same ──────────────────────────────────────

_raw(d::AbstractDict, k::Symbol) = haskey(d, String(k)) ? d[String(k)] : get(d, k, nothing)
_raw(d, k::Symbol) = hasproperty(d, k) ? getproperty(d, k) : nothing

# JSON `null` arrives as `nothing`, so absent and null are the same thing to every caller below.
function fld(d, k::Symbol, default = nothing)
    v = _raw(d, k)
    return v === nothing ? default : v
end

function fnum(d, k::Symbol, ctx; default = nothing)
    v = fld(d, k, default)
    v === nothing && fail("$ctx is missing the required numeric key \"$k\"")
    v isa Real && return float(v)
    return fail("$ctx has a non-numeric \"$k\" ($(repr(v))) — chartkit does not parse numbers out of strings")
end

function fstr(d, k::Symbol, ctx; default = nothing)
    v = fld(d, k, default)
    v === nothing && fail("$ctx is missing the required key \"$k\"")
    return string(v)
end

function fvec(d, k::Symbol, ctx)
    v = fld(d, k)
    v isa AbstractVector || fail("$ctx needs \"$k\" as an array of numbers")
    out = Float64[]
    for (i, e) in enumerate(v)
        e isa Real || fail("$ctx has a non-numeric \"$k\"[$i] ($(repr(e)))")
        push!(out, float(e))
    end
    return out
end

function frows(spec, k::Symbol, ctx)
    v = fld(spec, k)
    v isa AbstractVector || fail("$ctx needs \"$k\" as a non-empty array of objects")
    isempty(v) && fail("$ctx needs at least one entry in \"$k\"")
    return v
end

# ── numbers an exec artifact can quote ────────────────────────────────────────────────────────────

const CURRENCY = ("\$", "€", "£", "¥")

function split_unit(unit)
    u = string(unit)
    for c in CURRENCY
        startswith(u, c) && return (c, u[(ncodeunits(c) + 1):end])
    end
    return ("", u)
end

# Digits from magnitude: money reads at 0-1 decimals, a rate at 2, and nothing ever prints as
# 18.991000000000003. Pass `digits` to pin a published figure to its published precision.
function auto_digits(v)
    a = abs(float(v))
    (isfinite(a) && a > 0) || return 0
    a >= 100 && return 0
    a >= 10 && return 1
    a >= 1 && return 2
    return clamp(1 - floor(Int, log10(a)), 2, 6)
end

function group_thousands(digits::AbstractString)
    n = length(digits)
    io = IOBuffer()
    for (i, c) in enumerate(digits)
        i > 1 && (n - i + 1) % 3 == 0 && print(io, ',')
        print(io, c)
    end
    return String(take!(io))
end

# The one formatter every chart prints through. The unit wears itself: "\$M" → \$41M, "%" → 41%,
# "FTE" → 41 FTE, "/yr" → 41/yr. `signed` is for deltas, where the sign IS the finding.
function fmt(v; unit = "", signed = false, digits = nothing)
    isfinite(v) || return "n/a"
    d = digits === nothing ? auto_digits(v) : max(0, Int(digits))
    r = round(float(v); digits = d)
    r == 0 && (r = 0.0)                     # a rounded-away negative must not print as "-0"
    body = Printf.format(Printf.Format("%.$(d)f"), abs(r))
    parts = split(body, '.')
    body = group_thousands(parts[1]) * (length(parts) > 1 ? "." * parts[2] : "")
    sgn = r < 0 ? "-" : (signed ? "+" : "")
    pre, suf = split_unit(unit)
    gap = (isempty(pre) && !isempty(suf) && isletter(first(suf))) ? " " : ""
    return sgn * pre * body * gap * suf
end

# ── geometry, kept safe against degenerate input ──────────────────────────────────────────────────

# Coordinates are rounded and non-finites are pinned to 0: a NaN in a path silently kills the whole
# element in some renderers, so it must never reach the output.
function num(x)
    v = float(x)
    isfinite(v) || (v = 0.0)
    r = round(v; digits = 2)
    abs(r) < 1.0e9 && r == round(r) && return string(Int(r))
    return string(r)
end

esc_text(s) = replace(
    string(s),
    "&" => "&amp;", "<" => "&lt;", ">" => "&gt;", "\"" => "&quot;"
)

# A single value, an all-zero series or a flat trace must still produce a finite axis.
function axis_range(vals; include_zero = true, pad_frac = 0.06, clamp_zero = true)
    vs = Float64[float(v) for v in vals if isfinite(v)]
    isempty(vs) && return (0.0, 1.0)
    lo, hi = minimum(vs), maximum(vs)
    if include_zero
        lo = min(lo, 0.0)
        hi = max(hi, 0.0)
    end
    if hi <= lo
        p = max(abs(hi), 1.0) * 0.5
        return (lo - p, hi + p)
    end
    p = (hi - lo) * pad_frac
    # With one-signed data, leave zero flush at the plot edge instead of wasting a third of the width
    # on empty space the reader has to scan past.
    return (clamp_zero && lo >= 0 ? 0.0 : lo - p, clamp_zero && hi <= 0 ? 0.0 : hi + p)
end

function scaler(lo, hi, p0, p1)
    d = float(hi) - float(lo)
    abs(d) < 1.0e-12 && (d = 1.0)
    return v -> p0 + (float(v) - float(lo)) / d * (p1 - p0)
end

function nice_step(raw)
    (isfinite(raw) && raw > 0) || return 1.0
    mag = 10.0^floor(log10(raw))
    f = raw / mag
    m = f <= 1 ? 1.0 : f <= 2 ? 2.0 : f <= 2.5 ? 2.5 : f <= 5 ? 5.0 : 10.0
    return m * mag
end

function nice_ticks(lo, hi; target = 6)
    (isfinite(lo) && isfinite(hi) && hi > lo) || return (Float64[float(lo)], 0)
    step = nice_step((hi - lo) / max(target, 2))
    t = ceil(lo / step - 1.0e-9) * step
    out = Float64[]
    while t <= hi + step * 1.0e-6
        push!(out, abs(t) < step * 1.0e-9 ? 0.0 : t)
        t += step
    end
    isempty(out) && push!(out, lo)
    return (out, step >= 1 ? 0 : clamp(1 - floor(Int, log10(step)), 1, 4))
end

# Text metrics without a font: Inter averages ~0.55em of advance, the mono stack ~0.60em. Only used to
# size gutters, so a few percent off costs air, not a collision.
tw(s, size; mono = false) = length(string(s)) * size * (mono ? 0.60 : 0.55)

# ── SVG primitives ────────────────────────────────────────────────────────────────────────────────

function txt(x, y, s; size = 12.0, fill = INK, anchor = "start", family = SANS, weight = nothing)
    a = "<text x=\"$(num(x))\" y=\"$(num(y))\" font-family=\"$family\" font-size=\"$(num(size))\" fill=\"$fill\""
    anchor == "start" || (a *= " text-anchor=\"$anchor\"")
    weight === nothing || (a *= " font-weight=\"$weight\"")
    return a * ">" * esc_text(s) * "</text>"
end

function line(x1, y1, x2, y2; stroke = HAIR, sw = 1.0, dash = "", cap = "butt")
    a = "<line x1=\"$(num(x1))\" y1=\"$(num(y1))\" x2=\"$(num(x2))\" y2=\"$(num(y2))\" stroke=\"$stroke\" stroke-width=\"$(num(sw))\""
    isempty(dash) || (a *= " stroke-dasharray=\"$dash\"")
    cap == "butt" || (a *= " stroke-linecap=\"$cap\"")
    return a * "/>"
end

function rect(x, y, w, h; fill = "none", stroke = "none", sw = 1.0, opacity = 1.0)
    w < 0 && (x += w; w = -w)
    h < 0 && (y += h; h = -h)
    a = "<rect x=\"$(num(x))\" y=\"$(num(y))\" width=\"$(num(w))\" height=\"$(num(h))\" fill=\"$fill\""
    stroke == "none" || (a *= " stroke=\"$stroke\" stroke-width=\"$(num(sw))\"")
    opacity < 1 && (a *= " opacity=\"$(num(opacity))\"")
    return a * "/>"
end

function dot(cx, cy, r; fill = ACCENT, stroke = "none", sw = 1.0, opacity = 1.0)
    a = "<circle cx=\"$(num(cx))\" cy=\"$(num(cy))\" r=\"$(num(r))\" fill=\"$fill\""
    stroke == "none" || (a *= " stroke=\"$stroke\" stroke-width=\"$(num(sw))\"")
    opacity < 1 && (a *= " fill-opacity=\"$(num(opacity))\"")
    return a * "/>"
end

function polyline(pts; stroke = ACCENT, sw = 2.0)
    body = join(("$(num(p[1])),$(num(p[2]))" for p in pts), " ")
    return "<polyline points=\"$body\" fill=\"none\" stroke=\"$stroke\" stroke-width=\"$(num(sw))\" stroke-linejoin=\"round\" stroke-linecap=\"round\"/>"
end

# Bars get a rounded data-end and a square baseline end. A zero-length bar still gets a 0.4px sliver
# so an all-zero series reads as "measured zero" rather than "figure failed to draw".
function hbar(xb, xv, y, h; fill = ACCENT, stroke = "none", sw = 1.0, r = 4.0)
    s = xv < xb ? -1.0 : 1.0
    xe = xb + s * max(abs(xv - xb), 0.4)
    rr = min(r, abs(xe - xb), h / 2)
    xa = xe - s * rr
    d = "M $(num(xb)) $(num(y)) H $(num(xa)) Q $(num(xe)) $(num(y)) $(num(xe)) $(num(y + rr)) V $(num(y + h - rr)) Q $(num(xe)) $(num(y + h)) $(num(xa)) $(num(y + h)) H $(num(xb)) Z"
    a = "<path d=\"$d\" fill=\"$fill\""
    stroke == "none" || (a *= " stroke=\"$stroke\" stroke-width=\"$(num(sw))\"")
    return a * "/>"
end

function vbar(x, w, yb, yv; fill = ACCENT, stroke = "none", sw = 1.0, r = 4.0)
    s = yv < yb ? -1.0 : 1.0
    ye = yb + s * max(abs(yv - yb), 0.4)
    rr = min(r, abs(ye - yb), w / 2)
    ya = ye - s * rr
    d = "M $(num(x)) $(num(yb)) V $(num(ya)) Q $(num(x)) $(num(ye)) $(num(x + rr)) $(num(ye)) H $(num(x + w - rr)) Q $(num(x + w)) $(num(ye)) $(num(x + w)) $(num(ya)) V $(num(yb)) Z"
    a = "<path d=\"$d\" fill=\"$fill\""
    stroke == "none" || (a *= " stroke=\"$stroke\" stroke-width=\"$(num(sw))\"")
    return a * "/>"
end

function whisker_h(xlo, xhi, yc; cap = 7.0, stroke = INK, sw = 1.5)
    return line(xlo, yc, xhi, yc; stroke = stroke, sw = sw) *
        line(xlo, yc - cap / 2, xlo, yc + cap / 2; stroke = stroke, sw = sw) *
        line(xhi, yc - cap / 2, xhi, yc + cap / 2; stroke = stroke, sw = sw)
end

function whisker_v(xc, ylo, yhi; cap = 7.0, stroke = INK, sw = 1.5)
    return line(xc, ylo, xc, yhi; stroke = stroke, sw = sw) *
        line(xc - cap / 2, ylo, xc + cap / 2, ylo; stroke = stroke, sw = sw) *
        line(xc - cap / 2, yhi, xc + cap / 2, yhi; stroke = stroke, sw = sw)
end

# ── the shell every chart shares ──────────────────────────────────────────────────────────────────

const SUB_SIZE = 11.5

# A subtitle carries the provenance — seeds, discount rate, what the whisker is — which is the sentence
# a sceptical reader checks first, so it gets wrapped rather than shortened or silently clipped.
function wrap_lines(s, maxw, size; maxlines = 3)
    str = strip(string(s))
    isempty(str) && return String[]
    lines, cur = String[], ""
    for w in split(str)
        trial = isempty(cur) ? String(w) : cur * " " * w
        if tw(trial, size) <= maxw || isempty(cur)
            cur = trial
        else
            push!(lines, cur)
            cur = String(w)
            length(lines) == maxlines && break
        end
    end
    length(lines) < maxlines && !isempty(cur) && push!(lines, cur)
    # Anything past the last line is elided, so an over-long subtitle degrades visibly instead of
    # running off the viewBox where nobody sees that it went missing.
    consumed = sum(length, lines; init = 0) + length(lines) - 1
    consumed < length(str) && (lines[end] = rstrip(lines[end]) * "…")
    return lines
end

sub_lines(subtitle, width) = wrap_lines(subtitle, width - 2 * PAD, SUB_SIZE)

function head_height(title, subtitle; width = WIDTH)
    lines = sub_lines(subtitle, width)
    isempty(title) && isempty(lines) && return 0
    isempty(lines) && return 30
    return (isempty(title) ? 21 + 9 : 38 + 8) + 15 * (length(lines) - 1)
end

function head_block(title, subtitle; width = WIDTH)
    parts = String[]
    isempty(title) || push!(parts, txt(PAD, 21, title; size = 14.5, weight = 600))
    for (i, l) in enumerate(sub_lines(subtitle, width))
        push!(parts, txt(PAD, (isempty(title) ? 21 : 38) + 15 * (i - 1), l; size = SUB_SIZE, fill = GREY))
    end
    return join(parts)
end

# `<title>` first, so the accessible name is the chart's own label; `width="100%"` with a fixed height
# and a viewBox means the figure scales inside the page column without ever distorting.
function shell(body, height; width = WIDTH, title = "", subtitle = "")
    name = isempty(subtitle) ? title : (isempty(title) ? subtitle : "$title — $subtitle")
    io = IOBuffer()
    print(
        io, "<svg xmlns=\"http://www.w3.org/2000/svg\" viewBox=\"0 0 $(num(width)) $(num(height))\"",
        " width=\"100%\" height=\"$(num(height))\" preserveAspectRatio=\"xMidYMid meet\" role=\"img\">"
    )
    print(io, "<title>", esc_text(isempty(name) ? "chart" : name), "</title>")
    print(io, body)
    print(io, "</svg>\n")
    return String(take!(io))
end

# Wrap a category label onto at most two lines; a label too long for two gets elided rather than
# allowed to collide with its neighbour.
function wrap2(s, maxchars)
    str = string(s)
    length(str) <= maxchars && return (str, "")
    words = split(str)
    a, b = String[], String[]
    for w in words
        if isempty(b) && (isempty(a) || length(join([a; w], " ")) <= maxchars)
            push!(a, w)
        else
            push!(b, w)
        end
    end
    second = join(b, " ")
    length(second) > maxchars && (second = second[1:max(nextind(second, 0, min(maxchars - 1, length(second))) - 1, 1)] * "…")
    return (join(a, " "), second)
end

# ── bars_error ────────────────────────────────────────────────────────────────────────────────────
#
# spec: {"chart": "bars_error", "series": [{"label", "value", "se", "note"}], "unit", "digits",
#        "zero_label", "x_label", "sortby": "desc"|"asc"|"none", "order": [label…], "null_k", "title",
#        "subtitle", "width"}

function bars_error(
        series; unit = "", digits = nothing, title = "", subtitle = "", x_label = "",
        zero_label = "no effect", sortby = "desc", order = [], null_k = 2.0, width = WIDTH
    )
    isempty(series) && fail("bars_error needs at least one entry in \"series\"")
    ctx = "bars_error series"
    labels = [fstr(s, :label, "$ctx[$i]") for (i, s) in enumerate(series)]
    values = [fnum(s, :value, "$ctx[$i]") for (i, s) in enumerate(series)]
    ses = [fnum(s, :se, "$ctx[$i]"; default = 0.0) for (i, s) in enumerate(series)]
    notes = [string(fld(s, :note, "")) for s in series]
    any(<(0), ses) && fail("$ctx has a negative \"se\"; a standard error cannot be negative")

    idx = collect(eachindex(series))
    if !isempty(order)
        rank = Dict(string(order[i]) => i for i in eachindex(order))
        idx = sort(idx; by = i -> (get(rank, labels[i], length(order) + i), i))
    elseif sortby == "desc"
        idx = sort(idx; by = i -> (-values[i], i))
    elseif sortby == "asc"
        idx = sort(idx; by = i -> (values[i], i))
    elseif sortby != "none"
        fail("bars_error \"sortby\" must be \"desc\", \"asc\" or \"none\" (got $(repr(sortby)))")
    end

    # A bar whose ±null_k·se interval covers zero is the method saying "no effect". It has to read as a
    # deliberate null — hollow fill, grey label, an n.s. tag — and never as a small positive.
    isnull = [ses[i] > 0 && abs(values[i]) <= null_k * ses[i] for i in eachindex(series)]

    vlabels = [fmt(values[i]; unit = unit, signed = true, digits = digits) for i in eachindex(series)]
    row, bar_h = 32, 18
    lw = clamp(maximum(tw(l, 12.3) for l in labels) + 12, 90, 0.42 * width)
    vg = clamp(maximum(tw(v, 11.2; mono = true) for v in vlabels) + (any(isnull) ? 40 : 14), 60, 160)
    hh = head_height(title, subtitle; width = width)
    x0, x1 = PAD + lw, width - PAD - vg
    x1 <= x0 + 60 && (x1 = x0 + 60)
    y0 = hh + 12
    plot_h = row * length(series)
    y1 = y0 + plot_h
    height = y1 + 26 + (isempty(x_label) ? 0 : 17) + PAD

    # The axis has to contain the WIDEST thing drawn, which is the outer ±null_k·se whisker — size it
    # to ±1 se and a null's whisker escapes the plot and strikes through the row label beside it.
    wk = max(null_k, 1.0)
    lo, hi = axis_range(vcat(values .+ wk .* ses, values .- wk .* ses, 0.0))
    X = scaler(lo, hi, x0, x1)
    zx = X(0.0)

    P = String[head_block(title, subtitle; width = width)]
    for (k, i) in enumerate(idx)
        yc = y0 + (k - 0.5) * row
        iseven(k) && push!(P, rect(x0, y0 + (k - 1) * row, x1 - x0, row; fill = INSET))
        push!(
            P, txt(
                x0 - 10, yc + (isempty(notes[i]) ? 4.2 : 0.5), labels[i];
                size = 12.3, anchor = "end", fill = isnull[i] ? GREY : INK
            )
        )
        isempty(notes[i]) || push!(P, txt(x0 - 10, yc + 12, notes[i]; size = 9.6, anchor = "end", fill = GREY_LT))
    end
    push!(P, line(zx, y0, zx, y1; stroke = INK, sw = 1.0))

    for (k, i) in enumerate(idx)
        yc = y0 + (k - 0.5) * row
        bx = X(values[i])
        if isnull[i]
            push!(P, hbar(zx, bx, yc - bar_h / 2, bar_h; fill = INSET, stroke = HAIR, sw = 1.0))
        else
            push!(P, hbar(zx, bx, yc - bar_h / 2, bar_h; fill = values[i] < 0 ? ROSE : ACCENT))
        end
        if ses[i] > 0
            null_k > 1 && push!(
                P, whisker_h(
                    X(values[i] - null_k * ses[i]), X(values[i] + null_k * ses[i]), yc;
                    cap = 5.0, stroke = GREY, sw = 0.9
                )
            )
            push!(P, whisker_h(X(values[i] - ses[i]), X(values[i] + ses[i]), yc; stroke = INK, sw = 1.6))
        end
        # The value clears everything drawn on its row: the bar, the outer whisker, and the zero rule —
        # which is where a negative bar's label has to live, since its own side of the rule is full.
        tx = max(bx, ses[i] > 0 ? X(values[i] + wk * ses[i]) : bx, zx) + 9
        anchor = "start"
        push!(
            P, txt(
                tx, yc + 4, vlabels[i];
                size = 11.2, family = MONO, fill = isnull[i] ? GREY : INK, anchor = anchor,
                weight = isnull[i] ? nothing : 500
            )
        )
        isnull[i] && push!(
            P, txt(
                tx + tw(vlabels[i], 11.2; mono = true) + 7, yc + 4, "n.s.";
                size = 9.8, family = SANS, fill = GREY_LT, anchor = anchor
            )
        )
    end

    push!(P, line(x0, y1 + 0.5, x1, y1 + 0.5; stroke = HAIR, sw = 1.0))
    ticks, td = nice_ticks(lo, hi)
    for t in ticks
        push!(P, line(X(t), y1 + 0.5, X(t), y1 + 4; stroke = HAIR, sw = 1.0))
        push!(P, txt(X(t), y1 + 15, fmt(t; digits = td); size = 9.8, family = MONO, fill = GREY, anchor = "middle"))
    end
    if !isempty(zero_label)
        near_edge = zx - x0 < 46
        push!(
            P, txt(
                near_edge ? zx + 5 : zx, y0 - 3, zero_label;
                size = 9.8, fill = GREY, anchor = near_edge ? "start" : "middle"
            )
        )
    end
    isempty(x_label) || push!(P, txt((x0 + x1) / 2, y1 + 30, x_label; size = 10.6, fill = GREY, anchor = "middle"))

    return shell(join(P), height; width = width, title = title, subtitle = subtitle)
end

# ── waterfall ─────────────────────────────────────────────────────────────────────────────────────
#
# spec: {"chart": "waterfall", "start": {"label", "value", "se"}, "steps": [{"label", "value", "se"}],
#        "total": {"label", "value", "se"}, "unit", "digits", "y_label", "title", "subtitle", "width"}
# `total` may be omitted, in which case it is the arithmetic close of start + Σ steps. Whiskers are
# ±1 se in ink over ±`se_k` se in grey, the same convention `bars_error` draws.

function waterfall(
        start, steps, total = nothing; unit = "", digits = nothing, title = "", subtitle = "",
        y_label = "", se_k = 2.0, width = WIDTH
    )
    s_label = fstr(start, :label, "waterfall \"start\""; default = "Start")
    s_value = fnum(start, :value, "waterfall \"start\"")
    s_se = fnum(start, :se, "waterfall \"start\""; default = 0.0)
    isempty(steps) && fail("waterfall needs at least one entry in \"steps\"")
    st_labels = [fstr(s, :label, "waterfall steps[$i]") for (i, s) in enumerate(steps)]
    st_values = [fnum(s, :value, "waterfall steps[$i]") for (i, s) in enumerate(steps)]
    st_ses = [fnum(s, :se, "waterfall steps[$i]"; default = 0.0) for (i, s) in enumerate(steps)]
    t_label = total === nothing ? "Total" : fstr(total, :label, "waterfall \"total\""; default = "Total")
    t_value = total === nothing ? s_value + sum(st_values) : fnum(total, :value, "waterfall \"total\"")
    t_se = total === nothing ? 0.0 : fnum(total, :se, "waterfall \"total\""; default = 0.0)

    n = length(steps) + 2
    levels = Float64[s_value]                              # the running cumulative after each bar
    for v in st_values
        push!(levels, levels[end] + v)
    end
    all_labels = vcat([s_label], st_labels, [t_label])
    all_ses = vcat([s_se], st_ses, [t_se])

    hh = head_height(title, subtitle; width = width)
    x0 = PAD + (isempty(y_label) ? 44 : 60)
    x1 = width - PAD - 6
    y0 = hh + 18
    plot_h = 232
    y1 = y0 + plot_h
    height = y1 + 44 + (isempty(y_label) ? 0 : 0) + PAD

    # Every bar's whisker is centred on the level it ends at, and the axis has to contain the widest of
    # them — the outer ±se_k·se — or a whisker runs off the top of the plot.
    wk = max(se_k, 1.0)
    tops = vcat(levels, [t_value])
    extremes = vcat(tops, [0.0], tops .+ wk .* all_ses, tops .- wk .* all_ses)
    lo, hi = axis_range(extremes; pad_frac = 0.1)
    Y = scaler(lo, hi, y1, y0)
    zy = Y(0.0)
    slot = (x1 - x0) / n
    cw = min(56.0, slot * 0.62)

    P = String[head_block(title, subtitle; width = width)]
    ticks, td = nice_ticks(lo, hi; target = 5)
    for t in ticks
        push!(P, txt(x0 - 8, Y(t) + 3.4, fmt(t; digits = td); size = 9.8, family = MONO, fill = GREY, anchor = "end"))
    end
    push!(P, line(x0, zy, x1, zy; stroke = INK, sw = 1.0))
    isempty(y_label) || push!(
        P, "<g transform=\"rotate(-90 $(num(PAD + 4)) $(num((y0 + y1) / 2)))\">" *
            txt(PAD + 4, (y0 + y1) / 2, y_label; size = 10.6, fill = GREY, anchor = "middle") * "</g>"
    )

    for k in 1:n
        cx = x0 + (k - 0.5) * slot
        bx = cx - cw / 2
        base, top, fill = if k == 1
            (0.0, s_value, GREY)
        elseif k == n
            (0.0, t_value, ACCENT_D)
        else
            (levels[k - 1], levels[k], st_values[k - 1] < 0 ? AMBER : ACCENT_LT)
        end
        push!(P, vbar(bx, cw, Y(base), Y(top); fill = fill))
        k < n && push!(P, line(bx + cw, Y(levels[k]), x0 + (k + 0.5) * slot - cw / 2, Y(levels[k]); stroke = HAIR, sw = 1.0))

        se = all_ses[k]
        if se > 0
            wk > 1 && push!(P, whisker_v(cx, Y(top - wk * se), Y(top + wk * se); cap = 4.5, stroke = GREY, sw = 0.9))
            push!(P, whisker_v(cx, Y(top - se), Y(top + se); cap = 6.0, stroke = INK, sw = 1.4))
        end

        v = k == 1 ? s_value : (k == n ? t_value : st_values[k - 1])
        vl = fmt(v; unit = unit, signed = 1 < k < n, digits = digits)
        # Clear whichever is visually higher — the bar's top edge, or the whisker, which for a bar that
        # steps DOWN hangs off the bottom and must not push the label halfway up the plot.
        edge = min(Y(base), Y(top), se > 0 ? Y(top + wk * se) : Inf)
        ly = edge - 7
        ly < y0 + 8 && (ly = max(Y(base), Y(top)) + 14)
        push!(P, txt(cx, ly, vl; size = 10.6, family = MONO, fill = INK, weight = 500, anchor = "middle"))

        l1, l2 = wrap2(all_labels[k], max(Int(floor(slot / 5.8)), 8))
        push!(P, txt(cx, y1 + 17, l1; size = 10.2, fill = k == 1 || k == n ? INK : GREY, anchor = "middle"))
        isempty(l2) || push!(P, txt(cx, y1 + 29, l2; size = 10.2, fill = k == 1 || k == n ? INK : GREY, anchor = "middle"))
    end

    return shell(join(P), height; width = width, title = title, subtitle = subtitle)
end

# ── paired_dots ───────────────────────────────────────────────────────────────────────────────────
#
# spec: {"chart": "paired_dots", "arms": [{"label", "samples": [...]}], "delta", "delta_se", "unit",
#        "digits", "x_label", "title", "subtitle", "width"}
# With two arms and no "delta", the difference of means is computed; "delta_se" is never invented.

const PHI = 0.6180339887498949   # jitter from a golden-ratio sequence: spread, but fully determined

function paired_dots(
        arms; delta = nothing, delta_se = nothing, unit = "", digits = nothing, title = "",
        subtitle = "", x_label = "", width = WIDTH
    )
    isempty(arms) && fail("paired_dots needs at least one entry in \"arms\"")
    labels = [fstr(a, :label, "paired_dots arms[$i]") for (i, a) in enumerate(arms)]
    samples = [fvec(a, :samples, "paired_dots arms[$i]") for (i, a) in enumerate(arms)]
    any(isempty, samples) && fail("paired_dots needs a non-empty \"samples\" array for every arm")
    means = [sum(s) / length(s) for s in samples]
    colors = [ACCENT, AMBER, ACCENT_D, GREY]

    d = delta === nothing ? (length(arms) == 2 ? means[2] - means[1] : nothing) : float(delta)
    dse = delta_se === nothing ? nothing : float(delta_se)

    hh = head_height(title, subtitle; width = width)
    lab_h, band_h, gap_h = 20, 56, length(arms) == 2 ? 40 : 14
    x0 = PAD + 8
    x1 = width - PAD - 8
    y0 = hh + 10
    track = lab_h + band_h
    plot_h = track * length(arms) + gap_h * (length(arms) - 1)
    y1 = y0 + plot_h
    height = y1 + 28 + (isempty(x_label) ? 0 : 17) + PAD

    lo, hi = axis_range(vcat(samples...); include_zero = false, pad_frac = 0.07, clamp_zero = false)
    X = scaler(lo, hi, x0, x1)

    P = String[head_block(title, subtitle; width = width)]
    band_top = Float64[]
    for (i, lab) in enumerate(labels)
        ty = y0 + (i - 1) * (track + gap_h)
        push!(band_top, ty + lab_h)
        push!(P, dot(x0 + 4, ty + 9, 4.0; fill = colors[mod1(i, length(colors))]))
        push!(P, txt(x0 + 13, ty + 13, lab; size = 12.2, weight = 550))
        push!(
            P, txt(
                x1, ty + 13,
                "mean " * fmt(means[i]; unit = unit, digits = digits) * "  ·  n=$(length(samples[i]))";
                size = 10.4, family = MONO, fill = GREY, anchor = "end"
            )
        )
        push!(P, rect(x0, ty + lab_h, x1 - x0, band_h; fill = INSET))
    end

    if lo < 0 < hi
        for (i, _) in enumerate(labels)
            push!(P, line(X(0.0), band_top[i], X(0.0), band_top[i] + band_h; stroke = HAIR, sw = 1.0))
        end
    end

    inner = band_h - 16
    for (i, s) in enumerate(samples)
        col = colors[mod1(i, length(colors))]
        yc = band_top[i] + band_h / 2
        for (j, v) in enumerate(s)
            jit = (mod(j * PHI, 1.0) - 0.5) * inner
            push!(P, dot(X(v), yc + jit, 3.0; fill = col, stroke = PAPER, sw = 0.7, opacity = 0.62))
        end
        push!(P, line(X(means[i]), band_top[i] - 2, X(means[i]), band_top[i] + band_h + 2; stroke = INK, sw = 3.0))
    end

    # The bracket is the point of the figure: two means on one seed stream, and the gap between them
    # quoted with the standard error that makes it a claim rather than an anecdote.
    if length(arms) == 2 && d !== nothing
        mx1, mx2 = X(means[1]), X(means[2])
        yb = band_top[1] + band_h
        mid = yb + gap_h / 2
        push!(P, line(mx1, yb + 2, mx1, mid; stroke = GREY, sw = 1.0))
        push!(P, line(mx2, mid, mx2, band_top[2] - 2; stroke = GREY, sw = 1.0))
        push!(P, line(mx1, mid, mx2, mid; stroke = GREY, sw = 1.0))
        lab = "Δ " * fmt(d; unit = unit, signed = true, digits = digits) *
            (dse === nothing ? "" : " ± " * fmt(dse; unit = unit, digits = digits))
        cx = clamp((mx1 + mx2) / 2, x0 + tw(lab, 11.4; mono = true) / 2, x1 - tw(lab, 11.4; mono = true) / 2)
        push!(P, txt(cx, mid - 5, lab; size = 11.4, family = MONO, fill = INK, weight = 500, anchor = "middle"))
    end

    push!(P, line(x0, y1 + 0.5, x1, y1 + 0.5; stroke = HAIR, sw = 1.0))
    ticks, td = nice_ticks(lo, hi)
    for t in ticks
        push!(P, line(X(t), y1 + 0.5, X(t), y1 + 4; stroke = HAIR, sw = 1.0))
        push!(P, txt(X(t), y1 + 15, fmt(t; digits = td); size = 9.8, family = MONO, fill = GREY, anchor = "middle"))
    end
    isempty(x_label) || push!(P, txt((x0 + x1) / 2, y1 + 30, x_label; size = 10.6, fill = GREY, anchor = "middle"))

    return shell(join(P), height; width = width, title = title, subtitle = subtitle)
end

# ── trace_ceiling ─────────────────────────────────────────────────────────────────────────────────
#
# spec: {"chart": "trace_ceiling", "lines": [{"label", "y": [...], "emphasis": bool}], "ceiling",
#        "x": [...], "ceiling_label", "x_label", "y_label", "shade_at", "unit", "digits", "title",
#        "subtitle", "width"}

function trace_ceiling(
        x, lines; ceiling, ceiling_label = "", x_label = "", y_label = "", shade_at = 0.98,
        unit = "", digits = nothing, title = "", subtitle = "", width = WIDTH
    )
    isempty(lines) && fail("trace_ceiling needs at least one entry in \"lines\"")
    labels = [fstr(l, :label, "trace_ceiling lines[$i]") for (i, l) in enumerate(lines)]
    ys = [fvec(l, :y, "trace_ceiling lines[$i]") for (i, l) in enumerate(lines)]
    emph = [fld(l, :emphasis, false) == true for l in lines]
    xs = x === nothing ? collect(1.0:length(ys[1])) : Float64[float(v) for v in x]
    isempty(xs) && fail("trace_ceiling needs a non-empty \"x\"")
    for (i, y) in enumerate(ys)
        length(y) == length(xs) ||
            fail("trace_ceiling lines[$i] (\"$(labels[i])\") has $(length(y)) y values but \"x\" has $(length(xs))")
    end
    cap = float(ceiling)
    any(emph) || (emph = fill(true, length(lines)))

    hh = head_height(title, subtitle; width = width)
    rg = clamp(maximum(tw(l, 10.4) for l in labels) + 22, 70, 0.34 * width)
    x0 = PAD + (isempty(y_label) ? 42 : 58)
    x1 = width - PAD - rg
    x1 <= x0 + 80 && (x1 = x0 + 80)
    y0 = hh + 14
    plot_h = 244
    y1 = y0 + plot_h
    height = y1 + 28 + (isempty(x_label) ? 0 : 17) + PAD

    ylo, yhi = axis_range(vcat(vcat(ys...), [cap]); pad_frac = 0.0)
    yhi += (yhi - ylo) * 0.1
    Y = scaler(ylo, yhi, y1, y0)
    X = scaler(minimum(xs), maximum(xs), x0, x1)
    dx = length(xs) > 1 ? (x1 - x0) / (length(xs) - 1) : 0.0

    P = String[head_block(title, subtitle; width = width)]

    # Saturated runs: a cap of zero cannot be hit, so it is not shaded rather than shading everything.
    if cap > 0
        thr = cap * float(shade_at)
        hot = [any(ys[j][i] >= thr for j in eachindex(ys) if emph[j]) for i in eachindex(xs)]
        runs = Tuple{Int, Int}[]
        i = 1
        while i <= length(hot)
            if hot[i]
                j = i
                while j < length(hot) && hot[j + 1]
                    j += 1
                end
                push!(runs, (i, j))
                i = j + 1
            else
                i += 1
            end
        end
        for (a, b) in runs
            xa, xb = X(xs[a]) - dx / 2, X(xs[b]) + dx / 2
            push!(P, rect(max(xa, x0), y0, min(xb, x1) - max(xa, x0), y1 - y0; fill = ACCENT, opacity = 0.08))
            push!(P, rect(max(xa, x0), y1 - 3, min(xb, x1) - max(xa, x0), 3; fill = AMBER))
        end
        if !isempty(runs)
            a, b = argmax(r -> r[2] - r[1], runs)
            share = sum(hot) / length(hot) * 100
            # Annotate the rug, not the ceiling line: just under the cap is exactly where the trace that
            # is pinned TO the cap lives, so a label there lands on top of the data it describes.
            note = "at ceiling " * fmt(share; unit = "%", digits = 0) * " of the run"
            half = tw(note, 9.8) / 2
            cx = (max(X(xs[a]) - dx / 2, x0) + min(X(xs[b]) + dx / 2, x1)) / 2
            push!(
                P, txt(
                    clamp(cx, x0 + half, max(x1 - half, x0 + half)), y1 - 9, note;
                    size = 9.8, fill = AMBER, anchor = "middle"
                )
            )
        end
    end

    ticks, td = nice_ticks(ylo, yhi; target = 5)
    for t in ticks
        push!(P, txt(x0 - 8, Y(t) + 3.4, fmt(t; digits = td); size = 9.8, family = MONO, fill = GREY, anchor = "end"))
    end
    push!(P, line(x0, y1 + 0.5, x1, y1 + 0.5; stroke = HAIR, sw = 1.0))
    xticks, xtd = nice_ticks(minimum(xs), maximum(xs); target = 6)
    for t in xticks
        push!(P, line(X(t), y1 + 0.5, X(t), y1 + 4; stroke = HAIR, sw = 1.0))
        push!(P, txt(X(t), y1 + 15, fmt(t; digits = xtd); size = 9.8, family = MONO, fill = GREY, anchor = "middle"))
    end

    # The ceiling is a limit, not a gridline: dashed so it never reads as data, and labelled in the
    # gutter where it cannot collide with the trace.
    push!(P, line(x0, Y(cap), x1, Y(cap); stroke = INK, sw = 1.3, dash = "5 3"))

    gutter = Tuple{Float64, String, String, Float64}[]   # (y, text, colour, key-length)
    push!(
        gutter,
        (
            Y(cap), isempty(ceiling_label) ? fmt(cap; unit = unit, digits = digits) * " ceiling" : ceiling_label,
            INK, 0.0,
        )
    )
    for (j, lab) in enumerate(labels)
        pts = [(X(xs[i]), Y(ys[j][i])) for i in eachindex(xs)]
        col = emph[j] ? ACCENT : GREY
        if length(pts) > 1
            push!(P, polyline(pts; stroke = col, sw = emph[j] ? 2.0 : 1.4))
        end
        push!(P, dot(pts[end][1], pts[end][2], 3.4; fill = col, stroke = PAPER, sw = 1.6))
        push!(gutter, (pts[end][2], lab, emph[j] ? INK : GREY, 10.0))
    end
    # End labels detach from their marks if they are stacked, so nudge only as far as legibility needs.
    sort!(gutter; by = g -> g[1])
    for i in 2:length(gutter)
        gy, gt, gc, gk = gutter[i]
        gy - gutter[i - 1][1] < 12 && (gutter[i] = (gutter[i - 1][1] + 12, gt, gc, gk))
    end
    for (gy, gt, gc, gk) in gutter
        gk > 0 && push!(P, line(x1 + 5, gy, x1 + 5 + gk, gy; stroke = gc == INK ? ACCENT : GREY_LT, sw = 2.0))
        push!(P, txt(x1 + 6 + gk + (gk > 0 ? 5 : 0), gy + 3.4, gt; size = 10.4, fill = gc))
    end

    isempty(x_label) || push!(P, txt((x0 + x1) / 2, y1 + 30, x_label; size = 10.6, fill = GREY, anchor = "middle"))
    isempty(y_label) || push!(
        P, "<g transform=\"rotate(-90 $(num(PAD + 4)) $(num((y0 + y1) / 2)))\">" *
            txt(PAD + 4, (y0 + y1) / 2, y_label; size = 10.6, fill = GREY, anchor = "middle") * "</g>"
    )

    return shell(join(P), height; width = width, title = title, subtitle = subtitle)
end

# ── the spec surface ──────────────────────────────────────────────────────────────────────────────

const CHARTS = ("bars_error", "waterfall", "paired_dots", "trace_ceiling")

const USAGE = """
usage: julia --project=. agent/scripts/chartkit.jl SPEC.json [--out=fig.svg]   (--out=- for stdout)

Every spec takes: "chart" (required), "title", "subtitle", "unit" ("\$M", "%", "FTE", "/yr"),
"digits" (pin the printed precision), "width" (default 760).

  bars_error     "series": [{"label", "value", "se", "note"}]  (required)
                 "x_label", "zero_label", "sortby": desc|asc|none, "order": [label…], "null_k" (2.0)
  waterfall      "start": {"label", "value", "se"}, "steps": [{"label", "value", "se"}]  (required)
                 "total": {"label", "value", "se"} (default: start + Σ steps), "y_label", "se_k" (2.0)
  paired_dots    "arms": [{"label", "samples": [...]}]  (required)
                 "delta" (default: difference of means), "delta_se", "x_label"
  trace_ceiling  "lines": [{"label", "y": [...], "emphasis"}], "ceiling"  (required)
                 "x" (default 1:n), "ceiling_label", "x_label", "y_label", "shade_at" (0.98)
"""

function from_spec(spec)
    kind = string(fld(spec, :chart, ""))
    isempty(kind) && fail("spec has no \"chart\" key; expected one of $(join(CHARTS, ", "))")
    kind in CHARTS || fail("unknown chart type $(repr(kind)); expected one of $(join(CHARTS, ", "))")
    title = string(fld(spec, :title, ""))
    subtitle = string(fld(spec, :subtitle, ""))
    unit = string(fld(spec, :unit, ""))
    digits = fld(spec, :digits, nothing)
    width = float(fld(spec, :width, WIDTH))

    if kind == "bars_error"
        return bars_error(
            frows(spec, :series, "bars_error");
            unit = unit, digits = digits, title = title, subtitle = subtitle,
            x_label = string(fld(spec, :x_label, "")),
            zero_label = string(fld(spec, :zero_label, "no effect")),
            sortby = string(fld(spec, :sortby, fld(spec, :sort, "desc"))), order = fld(spec, :order, []),
            null_k = float(fld(spec, :null_k, 2.0)), width = width
        )
    elseif kind == "waterfall"
        start = fld(spec, :start)
        start === nothing && fail("waterfall is missing the required key \"start\"")
        return waterfall(
            start, frows(spec, :steps, "waterfall"), fld(spec, :total);
            unit = unit, digits = digits, title = title, subtitle = subtitle,
            y_label = string(fld(spec, :y_label, "")), se_k = float(fld(spec, :se_k, 2.0)),
            width = width
        )
    elseif kind == "paired_dots"
        return paired_dots(
            frows(spec, :arms, "paired_dots");
            delta = fld(spec, :delta), delta_se = fld(spec, :delta_se),
            unit = unit, digits = digits, title = title, subtitle = subtitle,
            x_label = string(fld(spec, :x_label, "")), width = width
        )
    else
        ceiling = fld(spec, :ceiling)
        ceiling === nothing && fail("trace_ceiling is missing the required key \"ceiling\"")
        ceiling isa Real || fail("trace_ceiling \"ceiling\" must be a number (got $(repr(ceiling)))")
        return trace_ceiling(
            fld(spec, :x), frows(spec, :lines, "trace_ceiling");
            ceiling = float(ceiling), ceiling_label = string(fld(spec, :ceiling_label, "")),
            x_label = string(fld(spec, :x_label, "")), y_label = string(fld(spec, :y_label, "")),
            shade_at = float(fld(spec, :shade_at, 0.98)), unit = unit, digits = digits,
            title = title, subtitle = subtitle, width = width
        )
    end
end

function main(argv)
    (isempty(argv) || argv[1] in ("-h", "--help")) && (print(stderr, USAGE); return isempty(argv) ? 2 : 0)
    src = argv[1]
    out = ""
    for a in argv[2:end]
        startswith(a, "--out=") || (println(stderr, "chartkit: unexpected argument $(repr(a))"); return 2)
        out = a[(length("--out=") + 1):end]
    end
    isfile(src) || (println(stderr, "chartkit: no such file: $src"); return 2)

    spec = try
        JSON.parsefile(src)
    catch e
        println(stderr, "chartkit: $src is not valid JSON — ", sprint(showerror, e))
        return 2
    end
    svg = try
        from_spec(spec)
    catch e
        e isa ChartError || rethrow()
        println(stderr, "chartkit: ", e.msg, "\n  in $src")
        return 2
    end

    isempty(out) && (out = replace(src, r"\.json$" => "") * ".svg")
    if out == "-"
        print(stdout, svg)
    else
        write(out, svg)
        println("chartkit: wrote $out — $(fld(spec, :chart, "")), $(round(sizeof(svg) / 1024, digits = 1)) KiB, self-contained")
    end
    return 0
end

end # module ChartKit

if abspath(PROGRAM_FILE) == @__FILE__
    exit(ChartKit.main(ARGS))
end
