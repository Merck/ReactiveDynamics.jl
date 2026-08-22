# Phase-1 semantic tests — Hierarchical refinement & open-port composition (CONTRACT §11 / ADR 0009).
#
# The vertical granularity axis: declare a fragment's boundary as open PORTS (§A), splice a finer
# sub-model into a coarse transition plug-compatibly (§B refine/abstract), advisory boundary checks
# (§C), and compact authoring (§D @pipeline/@process, §E @compose). All AUTHORING-time and additive —
# it produces a plain ReactionNetwork. Built on the ADR-0003 Phase-2 ReactantSpec FK-repoint.

using ReactiveDynamics, Test
using Random, Distributions, DataFrames
# ADR 0015 renamed the store shim verbs to store vocabulary AND unexported them; import the ones
# used bare here (nrows = old nparts, row_ids = old parts).
using ReactiveDynamics: nrows, row_ids

const RD = ReactiveDynamics

@testset "Refinement & open-port composition (§11 / ADR 0009)" begin

    # ── §A: open-port role annotation ───────────────────────────────────────────────────────────
    @testset "§A: port roles default to :private and are set by set_port_role!/@port" begin
        f = @reaction_network begin
            1.0, A --> B, name => t
        end
        # default role is :private for every species
        @test RD.port_role(f, :A) == :private
        @test RD.port_role(f, :B) == :private
        RD.set_port_role!(f, :A => :input, :B => :output)
        @test RD.port_role(f, :A) == :input
        @test RD.port_role(f, :B) == :output
        @test RD.is_open_port(:input) && RD.is_open_port(:output)
        @test !RD.is_open_port(:private) && !RD.is_open_port(:shared)
        # @port sugar tags groups; unlisted species keep :private
        g = @reaction_network begin
            1.0, X --> Y, name => t
            1.0, Y --> Z, name => u
        end
        @port g X => input  Z => output  Y => shared
        @test RD.port_role(g, :X) == :input
        @test RD.port_role(g, :Z) == :output
        @test RD.port_role(g, :Y) == :shared
        # illegal role rejected
        @test_throws Exception RD.set_port_role!(f, :A => :bogus)
    end

    # ── §E: @compose — port-connected composition (also closes §7/J4 :E/:obs, J9 include_model) ──
    @testset "§E: @compose identifies output↔input ports by FK-repoint; private namespaced" begin
        f1 = @reaction_network begin
            1.0, raw --> mid, name => step1
        end
        RD.set_port_role!(f1, :raw => :input, :mid => :output)
        f2 = @reaction_network begin
            1.0, mid --> product, name => step2
        end
        RD.set_port_role!(f2, :mid => :input, :product => :output)
        m = @compose f1 f2
        names = m[:, :placeName]
        # the shared output→input port `mid` collapses to ONE species (FK-repoint)
        @test count(==(:mid), names) == 1
        # private/other species are namespaced per fragment
        @test :f1__raw in names
        @test :f2__product in names
        # both transitions survive (structural append, §7/J2)
        @test nrows(m, :T) == 2
        # the promoted ReactantSpec table is FK-exact: every static FK resolves and the two
        # transitions route through the single shared `mid` index.
        rs = RD.reactant_specs(m)
        @test all(r -> r.species == 0 || 1 <= r.species <= nrows(m, :S), rs)
        midix = RD.find_index(:mid, m)
        @test count(r -> r.species == midix, rs) == 2   # produced by step1, consumed by step2
    end

    @testset "§E/J4: @compose merges events (:E) and observables (:obs) of the fragments" begin
        f1 = @reaction_network begin
            1.0, A --> B, name => t1
            (B > 5) && (B -= 1)
        end
        f2 = @reaction_network begin
            1.0, C --> D, name => t2
        end
        n_ev = nrows(f1, :E)
        @test n_ev >= 1
        m = @compose f1 f2
        @test nrows(m, :E) == n_ev            # the event survived composition (was dropped pre-WS-3)
    end

    # ── §B: refine — boundary-matched splice via FK-repoint ──────────────────────────────────────
    @testset "§B: refine! splices a sub-model into a coarse transition, plug-compatibly" begin
        coarse = @reaction_network begin
            1.0, P1 --> P2, name => phase2
            1.0, P2 --> P3, name => phase3
        end
        sub = @reaction_network begin
            1.0, entry --> work, name => screen
            1.0, work --> exit, name => filing
        end
        RD.set_port_role!(sub, :entry => :input, :exit => :output)
        r = refine(coarse, :phase2, sub; ports = Dict(:P1 => :entry, :P2 => :exit))
        tnames = [r[i, :transName] for i in row_ids(r, :T)]
        # Invariant: the coarse transition T is removed…
        @test !(:phase2 in tnames)
        # …its sub-transitions are spliced in (namespaced)…
        @test any(n -> occursin("screen", string(n)), tnames)
        @test any(n -> occursin("filing", string(n)), tnames)
        # …and every OTHER transition is structurally unchanged (plug-compatibility, Invariant 1).
        @test :phase3 in tnames
        # boundary species keep their names/indices (P1, P2 unchanged; P3 untouched).
        @test :P1 in r[:, :placeName]
        @test :P2 in r[:, :placeName]
        @test :P3 in r[:, :placeName]
        # the sub's PRIVATE species is namespaced (not leaked as a bare name).
        @test any(n -> occursin("work", string(n)), r[:, :placeName])
        @test !(:work in r[:, :placeName])
        # refine is non-mutating on the input (refine = refine! on a deepcopy).
        @test :phase2 in [coarse[i, :transName] for i in row_ids(coarse, :T)]
        # the promoted table is FK-exact after the splice.
        RD.populate_reactant_specs!(r)
        @test all(x -> x.species == 0 || 1 <= x.species <= nrows(r, :S), RD.reactant_specs(r))
    end

    # ── §B round-trip: a refined spec serializes/reloads as a flat model (Invariant 5) ───────────
    @testset "§B/Invariant 5: a refined spec round-trips through JSON as a flat model" begin
        coarse = @reaction_network begin
            1.0, P1 --> P2, name => phase2
        end
        sub = @reaction_network begin
            1.0, entry --> exit, name => step
        end
        RD.set_port_role!(sub, :entry => :input, :exit => :output)
        r = refine(coarse, :phase2, sub; ports = Dict(:P1 => :entry, :P2 => :exit))
        @prob_params r
        json = RD.to_json_model(r; meta = Dict{String, Any}("tspan" => 5.0))
        r2 = RD.build_network_from_dict(RD.JSON.parse(json))
        # the reloaded flat model has the same species and transition counts (refinement left no
        # runtime trace — it is a plain ModelSpec).
        @test nrows(r2, :S) == nrows(r, :S)
        @test nrows(r2, :T) == nrows(r, :T)
        @test Set(r2[:, :placeName]) == Set(r[:, :placeName])
    end

    # ── §D: @pipeline expands to flow-genesis routing transitions ────────────────────────────────
    @testset "§D: @pipeline expands a phase chain to N flow transitions carrying (ct, pos)" begin
        p = @pipeline Project begin
            Discovery => Phase1:(ct = 1.0, pos = 0.4)
            Phase1 => Phase2:(ct = 2.0, pos = 0.6)
            Phase2 => Market:(ct = 1.0, pos = 0.9)
        end
        @test nrows(p, :T) == 3                       # one transition per edge
        @test Set(p[:, :placeName]) == Set([:Discovery, :Phase1, :Phase2, :Market])
        # each edge carries its per-edge cycletime / prob_of_success.
        cts = Dict(p[i, :transName] => p[i, :transCycleTime] for i in row_ids(p, :T))
        poss = Dict(p[i, :transName] => p[i, :transProbOfSuccess] for i in row_ids(p, :T))
        @test cts[:flow_Discovery_Phase1] == 1.0 && poss[:flow_Discovery_Phase1] == 0.4
        @test cts[:flow_Phase1_Phase2] == 2.0 && poss[:flow_Phase1_Phase2] == 0.6
        @test cts[:flow_Phase2_Market] == 1.0 && poss[:flow_Phase2_Market] == 0.9
        # a flow transition consumes its upstream phase (upfront LHS) — the §2.8 flow idiom.
        RD.populate_reactant_specs!(p)
        p1ix = RD.find_index(:Phase1, p)
        @test any(r -> r.species == p1ix && r.side == :lhs, RD.reactant_specs(p))
    end

    # ── §D: @process — reusable parameterized fragment (eval-free param substitution) ────────────
    @testset "§D: @process instantiates a parameterized fragment; instances compose by ports" begin
        @process phase_gate(inp, outp; ct, pos) = begin
            1.0, inp --> outp, name => g, cycletime => ct, probability => pos
        end
        ga = phase_gate(:Phase1, :Phase2; ct = 2.0, pos = 0.6)
        @test Set(ga[:, :placeName]) == Set([:Phase1, :Phase2])
        @test nrows(ga, :T) == 1
        @test ga[1, :transCycleTime] == 2.0
        @test ga[1, :transProbOfSuccess] == 0.6
        # two instances sharing a port compose into a chain.
        gb = phase_gate(:Phase2, :Phase3; ct = 3.0, pos = 0.5)
        RD.set_port_role!(ga, :Phase1 => :input, :Phase2 => :output)
        RD.set_port_role!(gb, :Phase2 => :input, :Phase3 => :output)
        model = @compose ga gb
        @test count(==(:Phase2), model[:, :placeName]) == 1   # the shared port identified
        @test nrows(model, :T) == 2
    end

    # ── §C: advisory boundary-consistency diagnostics (warnings, not equivalence proofs) ─────────
    @testset "§C: refinement_diagnostics flags a dangling port and aggregate drift" begin
        # a sub whose cycletimes sum to 1.0, spliced under a coarse transition claiming ct=5.0 → drift.
        sub = @reaction_network begin
            1.0, entry --> work, name => s1, cycletime => 0.5
            1.0, work --> exit, name => s2, cycletime => 0.5
        end
        RD.set_port_role!(sub, :entry => :input, :exit => :output)
        warns = RD.refinement_diagnostics(sub, Dict(:transCycleTime => 5.0))
        @test any(w -> occursin("cycletime", w), warns)      # Σ sub ct (1.0) ≉ coarse ct (5.0)
        # a well-matched refinement (ct sums ≈ coarse, ports balanced) yields NO warnings.
        clean = RD.refinement_diagnostics(sub, Dict(:transCycleTime => 1.0))
        @test isempty(clean)
        # a dangling input port (declared :input but consumed by no sub-transition LHS) is flagged.
        # `sink` only ever appears on the RHS (produced), so declaring it an :input is a dangling port.
        bad = @reaction_network begin
            1.0, src --> sink, name => s
        end
        RD.set_port_role!(bad, :sink => :input)   # declared input but only ever PRODUCED (RHS)
        @test any(w -> occursin("dangling input", w), RD.refinement_diagnostics(bad, Dict()))
    end

    # ── Invariant 4 (closure): refine/@compose map ModelSpec(s) → a ModelSpec usable downstream ──
    @testset "Invariant 4: refine/compose results are valid inputs to further composition" begin
        base = @reaction_network begin
            1.0, P1 --> P2, name => phase2
        end
        sub = @reaction_network begin
            1.0, entry --> exit, name => step
        end
        RD.set_port_role!(sub, :entry => :input, :exit => :output)
        refined = refine(base, :phase2, sub; ports = Dict(:P1 => :entry, :P2 => :exit))
        # the refined model is a valid ReactionNetwork that can be @join'd again.
        extra = @reaction_network begin
            1.0, P2 --> P3, name => downstream
        end
        combined = @join refined extra
        @test combined isa RD.ReactionNetwork
        @test nrows(combined, :T) == nrows(refined, :T) + 1
    end

end
