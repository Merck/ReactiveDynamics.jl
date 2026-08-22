# ── ADR 0009 / CONTRACT §11 — hierarchical refinement & open-port composition ─────────────────
#
# The VERTICAL axis the flat §7 composition was silent on: declare a fragment's boundary as open
# PORTS (§A), splice a finer sub-model into a coarse transition plug-compatibly (§B refine/abstract),
# check the granularity ladder with advisory diagnostics (§C), and author compactly (§D @pipeline /
# @process, §E @compose). Everything here is AUTHORING-time and additive — it produces a plain
# ModelSpec that constructs/serializes/simulates exactly as a hand-written flat model. The enabling
# mechanism is the ADR-0003 Phase-2 ArcSpec FK-repoint: place identification is repointing an
# integer `place` FK, not string surgery. All of §11 is FORBIDDEN on a live/stepping model (it
# reindexes) — these operate on a static ReactionNetwork, never a ReactionNetworkProblem.

export refine!, refine, abstract!, abstract_transitions, set_port_role!, @port, @compose,
    @pipeline, @process, compose, refinement_diagnostics, port_role

# ── §A: port-role authoring ───────────────────────────────────────────────────────────────────

"""
    set_port_role!(net, name => role, …)

Set the open-port `role ∈ (:private, :input, :output, :shared)` (CONTRACT §11.1) of one or more
place by name. `:private` (default) auto-namespaces on compose; `:input`/`:output` are open ports
matched by `@compose`; `:shared` is identified by bare name (first-class `@catchall`).
"""
function set_port_role!(net::ReactionNetwork, pairs::Pair{Symbol, Symbol}...)
    for (name, role) in pairs
        role in PORT_ROLES ||
            error("set_port_role!: role must be one of $(PORT_ROLES), got $(repr(role))")
        i = find_index(name, net)
        i === nothing && error("set_port_role!: no place named $(repr(name))")
        net[i, :placeRole] = role
    end
    return net
end

"""
    @port net A => input  B => input  C => output  clock => shared

Declarative sugar for `set_port_role!`: tag place with a port role via `place => role` pairs
(role ∈ input/output/shared; anything unlisted keeps its default :private). Each pair is written with
`=>` (not `=`, which macro-call syntax parses as a keyword argument).
"""
macro port(netex, pairs...)
    call = Expr(:block)
    valid = (:input, :output, :shared, :private)
    for p in pairs
        (Meta.isexpr(p, :call) && p.args[1] === :(=>)) ||
            error("@port: each entry must be `place => role`, got $(p)")
        pl = p.args[2]
        role = p.args[3]
        role in valid || error("@port: role must be one of $(valid), got $(role)")
        push!(call.args, :(set_port_role!($(esc(netex)), $(QuoteNode(pl)) => $(QuoteNode(role)))))
    end
    push!(call.args, :($(esc(netex))))
    return call
end

# ── §E: @compose — port-connected composition (the §7 join, with a declared boundary) ─────────
#
# `compose(f1, f2, …)` is `merge_networks!`/@join PLUS automatic port matching: each fragment's `output`
# ports are identified with same-named `input` ports of the other fragments by the §7.4/J7 FK-repoint
# (via equalize!, which now repoints ArcSpec FKs — ADR 0003 Phase 2), `private` places are
# namespaced (m__X), and `shared` places are identified by bare name (prepend! skips them). Because
# it composes already-parsed ModelSpecs it CLOSES the §7/J4 (:E/:obs dropped — merge_networks! now merges
# them) and J9 (undefined include_model — never taken) bugs en route.

"""
    compose(fragments…; namespace=true)

Compose model fragments by matching open ports. `output` ports are identified with same-named
`input` ports across fragments (FK-repoint), `shared` places by bare name, `private` places are
namespaced per fragment. Returns a new `ReactionNetwork`. `@compose f1 f2 …` is the macro form.
"""
function compose(fragments::ReactionNetwork...)
    isempty(fragments) && return ReactionNetwork()
    # Collect, per fragment, its open-port place names (input/output) and shared names BEFORE any
    # namespacing, so we know which bare names to re-identify after the namespaced union.
    portnames = Set{Symbol}()
    for f in fragments
        for i in row_ids(f, :S)
            r = port_role(f, i)
            (is_open_port(r) || r === :shared) && push!(portnames, f[i, :placeName])
        end
    end

    merged = ReactionNetwork()
    # union each fragment under its own namespace. prepend! leaves `shared` place bare; open ports
    # (input/output) are namespaced here, then re-identified below by matching the ORIGINAL name.
    portmap = Dict{Symbol, Vector{Symbol}}()   # original port name → its namespaced aliases in merged
    for (k, f) in enumerate(fragments)
        name = Symbol("f", k)
        # remember each fragment's open-port original names → their namespaced form
        for i in row_ids(f, :S)
            r = port_role(f, i)
            if is_open_port(r)
                orig = f[i, :placeName]
                push!(get!(portmap, orig, Symbol[]), normalize_name(orig, name))
            end
        end
        merge_networks!(merged, f, name)
    end

    # Identify open ports that appear (as the same original name) in ≥2 fragments: their namespaced
    # aliases collapse to one place via equalize! (FK-repoint). A port present in only one fragment
    # stays a namespaced (dangling) open port — §C validate warns on it.
    eqs = Vector{Any}()
    for (orig, aliases) in portmap
        length(aliases) < 2 && continue
        block = Any[(:alias, orig)]
        for a in unique(aliases)
            push!(block, (:catchall, a))
        end
        push!(eqs, block)
    end
    isempty(eqs) || equalize!(merged, eqs)

    populate_arcs!(merged)
    return merged
end

"""
    @compose f1 f2 …

Macro form of [`compose`](@ref): compose declared-port model fragments (each an expression evaluating
to a `ReactionNetwork`). The explicit-boundary counterpart of `@join`; `@join`/`@equalize`
remain the manual no-declared-ports path.
"""
macro compose(exs...)
    return :(compose($(map(esc, exs)...)))
end

# ── §B: refine / abstract — boundary-matched splice via FK-repoint ────────────────────────────
#
# `refine!(spec, T, sub; ports)` replaces the coarse transition named `T` with the sub-model `sub`,
# plug-compatibly at its boundary, in four authoring-time structural moves (CONTRACT §11.2):
#   1. namespace `sub`'s `private` place (leave input/output/shared un-prefixed for matching);
#   2. identify `sub`'s open ports with the parent's boundary places per `ports` by FK-repoint;
#   3. append `sub`'s transitions + remaining places/params/obs/EVENTS (this also merges :E/:obs);
#   4. remove the coarse transition `T` (and its ArcSpec rows).
# Because the boundary places keep their indices/names/attributes, every transition NOT in {T}∪sub
# is structurally unchanged (Invariant 1, plug-compatibility). Forbidden on a live model (reindexes).

"""
    refine!(spec, transition, submodel; ports = Dict(boundary_places => sub_port, …))

Splice `submodel` into the coarse `transition` (named `Symbol`) of `spec`, identifying each of the
submodel's open ports (`sub_port`) with the parent boundary places (`boundary_places`) given in
`ports`. Mutates and returns `spec`. Authoring-time only.
"""
function refine!(
        spec::ReactionNetwork, transition::Symbol, submodel::ReactionNetwork;
        ports::AbstractDict = Dict{Symbol, Symbol}()
    )
    # locate the coarse transition row by name
    ti = findfirst(i -> spec[i, :transName] === transition, collect(row_ids(spec, :T)))
    ti === nothing && error("refine!: no transition named $(repr(transition)) in the parent spec")

    name = Symbol(transition, :__sub)

    # Validate the port map: each parent boundary place and each sub port must exist.
    for (boundary, subport) in ports
        find_index(boundary, spec) === nothing &&
            error("refine!: boundary place $(repr(boundary)) not found in parent")
        find_index(subport, submodel) === nothing &&
            error("refine!: port place $(repr(subport)) not found in submodel")
    end

    # Moves 1+2 are delegated to merge_networks!'s own namespacing + equation-alias mechanism (§7.4/J7):
    # build one eqs block PER port that aliases the sub's port place to the parent boundary name, so
    # `prepend!`/`normalize_name` rename the port to the boundary name (bare) while every PRIVATE
    # place is namespaced `<name>__X`. merge_networks! then merges the boundary-named port onto the
    # existing parent row (incident by placeName) — the structural FK-repoint — and appends the rest.
    # `shared`-role sub places are left bare by prepend! (§A) and merge onto any same-named parent row.
    eqs = Any[]
    for (boundary, subport) in ports
        push!(eqs, Any[(:alias, boundary), (:catchall, subport)])
    end

    # Move 3: append the sub into the parent by namespaced name-merge with the port aliases. merge_networks!
    # merges :E/:obs uniformly (WS-3), so the sub's events/observables come along.
    merge_networks!(spec, submodel, name, eqs)

    # Move 4: remove the coarse transition T (drop its :T row). Its arc relation lived only in
    # its :trans Expr, so dropping the row removes it; the sub's transitions now carry the dynamics.
    rem_rows!(spec, :T, [ti])

    populate_arcs!(spec)
    return spec
end

"""
    refine(spec, transition, submodel; ports = Dict(boundary_places => sub_port, …)) -> ReactionNetwork

Non-mutating convenience over [`refine!`](@ref): splice `submodel` into the coarse `transition` of a `deepcopy` of `spec`, returning the refined copy and leaving `spec` untouched. Same port-matching semantics and authoring-time-only restriction as `refine!`.
"""
refine(spec::ReactionNetwork, transition::Symbol, submodel::ReactionNetwork; kwargs...) =
    refine!(deepcopy(spec), transition, submodel; kwargs...)

"""
    abstract_transitions(spec, transitions, into; boundary)

Inverse of `refine!`: collapse a connected set of sub-transitions (by name) into a single coarse
transition named `into`, whose boundary reaction line consumes/produces the given `boundary`
place. A structural convenience for round-tripping the granularity ladder; the collapsed coarse
transition's attributes (cycletime/pos/cost) are the caller's to summarize (§C advises on drift).
"""
function abstract_transitions(
        spec::ReactionNetwork, transitions::Vector{Symbol}, into::Symbol;
        lhs::Vector{Symbol} = Symbol[], rhs::Vector{Symbol} = Symbol[],
        attrs::AbstractDict = Dict{Symbol, Any}()
    )
    spec = deepcopy(spec)
    tis = Int[]
    for tn in transitions
        i = findfirst(j -> spec[j, :transName] === tn, collect(row_ids(spec, :T)))
        i === nothing && error("abstract_transitions: no transition named $(repr(tn))")
        push!(tis, i)
    end
    # build the coarse reaction line LHS --> RHS from the boundary places
    lhs_ex = isempty(lhs) ? :∅ : foldl((a, b) -> :($a + $b), lhs)
    rhs_ex = isempty(rhs) ? :∅ : foldl((a, b) -> :($a + $b), rhs)
    line = :($lhs_ex --> $rhs_ex)
    ni = add_row!(spec, :T; trans = line, transName = into)
    for (k, v) in attrs
        spec[ni, k] = v
    end
    assign_defaults!(spec)
    rem_rows!(spec, :T, sort(tis))
    populate_arcs!(spec)
    return spec
end

"""
    abstract!(spec, transitions, into; lhs, rhs, attrs)

Alias for [`abstract_transitions`](@ref): collapse a connected set of sub-transitions into one coarse transition named `into`. The inverse of [`refine!`](@ref).
"""
const abstract! = abstract_transitions

# ── §C: advisory boundary-consistency diagnostics ────────────────────────────────────────────
#
# Refinement does NOT claim the fine model is behaviorally equivalent to the coarse one (that needs a
# bisimulation the framework can't check). These are ADVISORY warnings comparing the coarse
# transition's attributes to aggregates of the refinement where computable (CONTRACT §11.3): a
# port-balance check (every input port consumed by ≥1 sub-transition LHS, every output produced by
# ≥1 sub-transition RHS) and, on a LINEAR chain, coarse.cycletime ≈ Σ sub cycletimes and
# coarse.pos ≈ Π sub PoS. Returns a Vector{String} of warnings (empty = clean); the author overrides
# by ignoring them — refinement legitimately changes dynamics.

"""
    refinement_diagnostics(spec, submodel, coarse_attrs; ports) -> Vector{String}

Advisory §11.3 diagnostics for splicing `submodel` into a coarse transition described by
`coarse_attrs` (a Dict of e.g. `:transCycleTime`, `:transProbOfSuccess`). Warnings only.
"""
function refinement_diagnostics(
        submodel::ReactionNetwork, coarse_attrs::AbstractDict;
        ports::AbstractDict = Dict{Symbol, Symbol}(), tol = 0.25
    )
    warns = String[]

    # port-balance: an `input` port should be consumed by some sub-transition LHS; an `output` port
    # produced by some RHS. We check via the promoted ArcSpec table on a populated copy.
    sub = deepcopy(submodel)
    populate_arcs!(sub)
    lhs_places = Set(r.place for r in arcs(sub) if r.side === :lhs && r.place > 0)
    rhs_places = Set(r.place for r in arcs(sub) if r.side === :rhs && r.place > 0)
    for i in row_ids(sub, :S)
        role = port_role(sub, i)
        nm = sub[i, :placeName]
        if role === :input && !(i in lhs_places)
            push!(warns, "input port $(nm) is not consumed by any sub-transition (dangling input)")
        elseif role === :output && !(i in rhs_places)
            push!(warns, "output port $(nm) is not produced by any sub-transition (dangling output)")
        end
    end

    # linear-chain aggregate checks (best-effort; skipped if any needed attr is non-numeric/absent).
    cts = Float64[]
    poss = Float64[]
    ok = true
    for i in row_ids(sub, :T)
        ct = sub[i, :transCycleTime]
        ps = sub[i, :transProbOfSuccess]
        (ct isa Number && ps isa Number) || (ok = false; break)
        push!(cts, Float64(ct))
        push!(poss, Float64(ps))
    end
    if ok && !isempty(cts)
        if haskey(coarse_attrs, :transCycleTime) && coarse_attrs[:transCycleTime] isa Number
            cc = Float64(coarse_attrs[:transCycleTime])
            sc = sum(cts)
            (cc == 0 || abs(cc - sc) <= tol * max(cc, sc)) ||
                push!(warns, "coarse cycletime $(cc) ≉ Σ sub cycletimes $(sc) (>$(round(Int, tol * 100))% drift)")
        end
        if haskey(coarse_attrs, :transProbOfSuccess) && coarse_attrs[:transProbOfSuccess] isa Number
            cp = Float64(coarse_attrs[:transProbOfSuccess])
            sp = prod(poss)
            abs(cp - sp) <= tol ||
                push!(warns, "coarse prob_of_success $(cp) ≉ Π sub PoS $(round(sp, digits = 3)) (>$(tol) drift)")
        end
    end

    return warns
end

# ── §D: @pipeline sugar and @process modules ──────────────────────────────────────────────────
#
# @pipeline expands a chain of phases into N `flow`-genesis routing transitions (CONTRACT §2.8 flow:
# the upstream phase is an upfront-consumed LHS, so the transition fires only when the prior phase
# produced a token — token-flow, not an independent Poisson clock). Each edge `From => To : (ct, pos,
# res)` becomes `res_expr + From --> To` with the per-edge cycletime/pos. @process is a named
# parameterized ModelSpec fragment, instantiated by eval-free parameter substitution and composed by
# ports (§E).

"""
    @pipeline name begin
        A => B : (ct=1.0, pos=0.4, res = 2*@conserved(scientist))
        B => C : (ct=2.0, pos=0.6)
    end

Expand a phase chain into `flow` routing transitions (§2.8). Each `From => To : (ct, pos, res)` edge
becomes a transition consuming `From` (+ optional `res` resources) and producing `To`, carrying the
per-edge cycletime/prob_of_success. Returns a `ReactionNetwork`.
"""
macro pipeline(nameex, block)
    Meta.isexpr(block, :block) || error("@pipeline: expected a begin…end block of `From => To : opts`")
    # Build a standard @reaction_network authoring block, one reaction line per edge, and reuse
    # the full parse pipeline (rate expansion, place extraction, attr handling) via get_data.
    lines = Expr(:block)
    for stmt in block.args
        stmt isa LineNumberNode && continue
        # `From => To : (opts…)` parses as `Pair(From, (To : opts))` — a `=>` call whose RHS is a
        # `:` (colon) call `(To, opts)`. Unwrap both.
        (Meta.isexpr(stmt, :call) && stmt.args[1] === :(=>)) ||
            error("@pipeline: each edge must be `From => To : (ct=…, pos=…, res=…)`, got $(stmt)")
        from = stmt.args[2]
        colonex = stmt.args[3]
        (Meta.isexpr(colonex, :call) && colonex.args[1] === :(:)) ||
            error("@pipeline: edge must carry `: (ct=…, pos=…)` options, got $(colonex)")
        to = colonex.args[2]
        opts = colonex.args[3]
        optd = Dict{Symbol, Any}()
        if Meta.isexpr(opts, :tuple)
            for kv in opts.args
                Meta.isexpr(kv, :(=)) && (optd[kv.args[1]] = kv.args[2])
            end
        elseif Meta.isexpr(opts, :(=))
            optd[opts.args[1]] = opts.args[2]
        end
        ct = get(optd, :ct, 1.0)
        pos = get(optd, :pos, 1.0)
        tname = Symbol("flow_", from, "_", to)
        # flow genesis (§2.8): the upstream phase `From` is consumed upfront (LHS), so the transition
        # fires only when a token exists there; a high nominal @deterministic rate makes the token
        # gate (not a Poisson draw) bound firing. Optional `res` resources add to the LHS.
        lhs = haskey(optd, :res) ? :($(optd[:res]) + $from) : from
        # Build the reaction-line tuple exactly as the DSL author would WRITE it, so get_data parses
        # it identically to a hand-authored `@reaction_network` line: `@deterministic(1e6), From
        # --> To, name => <tname>, cycletime => ct, probability => pos`. `name`'s value must be a BARE
        # identifier Symbol (get_transitions! reads `exs[ix].args[3]`), NOT a QuoteNode — matching the
        # parser's `transName => :flow_…` output. `-->` is the reaction arrow the DSL normalizes.
        line = Expr(
            :tuple,
            Expr(:macrocall, Symbol("@deterministic"), LineNumberNode(0, :pipeline), 1.0e6),
            Expr(:-->, lhs, to),
            Expr(:call, :(=>), :name, tname),
            Expr(:call, :(=>), :cycletime, ct),
            Expr(:call, :(=>), :probability, pos),
        )
        push!(lines.args, line)
    end
    # Delegate to the create.jl parse pipeline exactly as @reaction_network does. The reaction
    # lines are built from literal phase symbols (no caller-scope variables), so nothing needs esc;
    # reference RD's own ReactionNetwork/get_data by module-qualified name.
    return :(ReactiveDynamics.ReactionNetwork(ReactiveDynamics.get_data($(QuoteNode(lines)))...))
end

"""
    @process name(params…) = begin <reaction lines> end

Define a reusable parameterized model-fragment factory. Expands to a function `name(params…)` that
returns a `ReactionNetwork`. Inside the body, write ordinary reaction lines (as in
`@reaction_network`); each occurrence of a PARAMETER name is substituted by its call-time value
(a place symbol, a number, …) into the reaction-line AST BEFORE parsing — eval-free
(`replace_in_expr`), sidestepping the DSL's lack of `\$`-interpolation. Compose instances by ports
with `@compose` (§E).

```julia
@process phase_gate(inp, outp; ct, pos) = begin
    1.0, inp --> outp, name => g, cycletime => ct, probability => pos
end
gate_a = phase_gate(:Phase1, :Phase2; ct=2.0, pos=0.6)
```
"""
macro process(defex)
    (Meta.isexpr(defex, :(=)) || Meta.isexpr(defex, :function)) ||
        error("@process: expected `name(params…) = begin <reaction lines> end`")
    sig = defex.args[1]
    body = defex.args[2]
    Meta.isexpr(sig, :call) ||
        error("@process: signature must be `name(params…)`, got $(sig)")
    # collect parameter NAMES (positional + keyword) to substitute into the reaction block.
    pnames = Symbol[]
    for a in sig.args[2:end]
        if a isa Symbol
            push!(pnames, a)
        elseif Meta.isexpr(a, :parameters)          # keyword params `; ct, pos`
            for kw in a.args
                kw isa Symbol && push!(pnames, kw)
                Meta.isexpr(kw, :kw) && push!(pnames, kw.args[1])
            end
        elseif Meta.isexpr(a, :kw)
            push!(pnames, a.args[1])
        end
    end
    # Emit: name(params…) = ReactionNetwork(get_data(<block with params replaced>)…). The
    # block is quoted, then each param symbol is replaced by its runtime value via replace_in_expr
    # (create.jl) — a structural substitution, no eval.
    blockq = QuoteNode(body)
    subs = Expr(:vect, (Expr(:call, :(=>), QuoteNode(p), p) for p in pnames)...)
    newfn = Expr(
        :function, esc(sig),
        :(
            ReactiveDynamics.ReactionNetwork(
                ReactiveDynamics.get_data(
                    ReactiveDynamics.replace_in_expr($blockq, $(esc(subs))...)
                )...
            )
        )
    )
    return newfn
end
