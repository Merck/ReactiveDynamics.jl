# Phase-0 semantic test suite — aggregator.
#
# These tests encode the Phase-0 modeling contract (docs/CONTRACT_DRAFT.md) and the ADRs
# (docs/adr/). Unlike the legacy tutorial smoke tests (which only assert "did not throw"),
# these are real assertions over the engine's operational semantics.
#
# Two tiers (see each file's header):
#   T1-characterization — runs against the CURRENT engine; locks in behavior, or pins a known
#                         bug via @test_broken so the suite is green-when-expected.
#   T2-acceptance       — encodes TARGET behavior per the contract; references not-yet-built
#                         APIs (progressive_fill!, seed=, ReactantSpec FK-repoint, …) and is
#                         wrapped with @test_skip / commented blocks so the file still loads.
#
# As Phase 1 implements each piece, flip the corresponding @test_skip / @test_broken to @test.

using SafeTestsets

@safetestset "Allocation / Conservation / Lifecycle" begin
    include("allocation_conservation_lifecycle.jl")
end
@safetestset "Modality / Genesis" begin
    include("modality_genesis.jl")
end
@safetestset "Determinism / Composition / Bug-pins" begin
    include("determinism_composition_bugs.jl")
end
@safetestset "Reference models (SIR / toy-pharma / rNPV)" begin
    include("reference_models.jl")
end
@safetestset "Rules / Decisions (endogenous channel)" begin
    include("rules_decisions.jl")
end
@safetestset "Token filtration (@select / @advance)" begin
    include("token_filtration.jl")
end
@safetestset "Declarative initial state + checkpoint" begin
    include("initial_state.jl")
end
