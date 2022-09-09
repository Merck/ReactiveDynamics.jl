# reaction network DSL: CREATE part; reaction line and event parsing 

export @ReactionNetwork

using MacroTools: prewalk, postwalk, striplines, isexpr
using Symbolics: build_function, get_variables

empty_set = Set{Symbol}([:∅])
fwd_arrows = Set{Symbol}([:>, :→, :↣, :↦, :⇾, :⟶, :⟼, :⥟, :⥟, :⇀, :⇁, :⇒, :⟾])
bwd_arrows = Set{Symbol}([:<, :←, :↢, :↤, :⇽, :⟵, :⟻, :⥚, :⥞, :↼, :↽, :⇐, :⟽, Symbol("<--")])
double_arrows = Set{Symbol}([:↔, :⟷, :⇄, :⇆, :⇌, :⇋, :⇔, :⟺, Symbol("<-->")])

arrows = fwd_arrows ∪ bwd_arrows ∪ double_arrows ∪ [:-->]
ifs = [:&&, :if]

reserved_sampling_macros = [:register, :sample, :take]

struct Reactant
    species::Symbol
    stoich::SampleableValues
    modality::Set{Symbol}
end

struct FoldedReactionStruct
    rate::SampleableValues
    reaction::SampleableValues
end

struct Event
    trigger::SampleableValues
    action::SampleableValues
end

# Declares symbols which may neither be used as parameters not varriables.
forbidden_symbols = [:t, :π, :pi, :ℯ, :im, :nothing, :∅]

"""
Macro that takes an expression corresponding to a reaction network and outputs an instance of `TheoryReactionNetwork` that can be converted to a `DiscreteProblem` or solved directly.

Most arrows accepted (both right, left, and bi-drectional arrows). Use 0 or ∅ for annihilation/creation to/from nothing.

Custom functions and sampleable objects can be used as numeric parameters. Note that these have to be accessible from ReactionDynamics's source code.

# Examples
```julia
acs = @ReactionNetwork begin
    1.0, X ⟶ Y
    1.0, X ⟶ Y, priority=>6., prob=>.7, capacity=>3.
    1.0, ∅ --> (Poisson(.3γ)X, Poisson(.5)Y)
    (XY > 100) && (XY -= 1)
end
@push acs 1.0 X ⟶ Y 
@prob_init acs X=1 Y=2 XY=α
@prob_params acs γ=1 α=4
@solve_and_plot acs
```
"""
macro ReactionNetwork end

macro ReactionNetwork() make_ReactionNetwork(:()) end

macro ReactionNetwork(ex) make_ReactionNetwork(ex; eval_module=__module__) end

macro ReactionNetwork(ex, args...)
    make_ReactionNetwork(generate(Expr(:braces, ex, args...); eval_module=__module__); eval_module=__module__)
end

function make_ReactionNetwork(ex::Expr; eval_module=@__MODULE__)
    blockex = generate(ex; eval_module); blockex = unblock_shallow!(blockex)
    
    :(ReactionNetwork(get_data($(QuoteNode(blockex)))...))
end

### Functions that process the input and rephrase it as a reaction system ###
function esc_dollars!(ex)
    if ex isa Expr
        if ex.head == :$
            return esc(:($(ex.args[1])))
        else
            for i = 1:length(ex.args)
                ex.args[i] = esc_dollars!(ex.args[i])
            end            
        end
    end
    ex
end

symbolize(pairex) = pairex isa Number ? pairex : (pairex.args[2] => pairex.args[3])

function get_data(ex)
    trans = []; evs = []; reactants = []; pcs = []
    if isexpr(ex, :block)
        ex = striplines(ex); esc_dollars!(ex)
        foreach(l -> get_data!(trans, reactants, pcs, evs, isexpr(l, :tuple) ? l.args : [l]), ex.args)
    elseif ex != :() get_data!(trans, reactants, pcs, evs, ex) end

    trans, reactants, pcs, evs
end

get_data!(trans, reactants, pcs, evs, exs) =
    (length(exs) == 0 && return; exs[1] isa Expr && (exs[1].head ∈ ifs) ? get_events!(evs, normalize_pcs!(pcs, exs[1])) : 
        get_transitions!(trans, reactants, pcs, exs))

function get_events!(evs, ex)
    if ex.head == :&& push!(evs, Event(ex.args...))
    else recursively_expand_actions!(evs, Expr(:call, :&), ex) end
end

function recursively_expand_actions!(evs, condex, event)
    if isexpr(event, :if)
        condex_ = deepcopy(condex)
        push!(condex_.args, event.args[1]); push!(evs, Event(condex_, event.args[2]))
        push!(condex.args, Expr(:call, :!, event.args[1]))
        length(event.args) >= 3 && recursively_expand_actions!(evs, condex, event.args[3])
    else push!(evs, Event(condex, event)) end
end

function prune_rate(rate)
    if (isexpr(rate, :macrocall) && (macroname(rate) ∈ prettynames[:transCycleTime]))
        :(1/$(rate.args[3]))
    else rate end
end

function get_transitions!(trans, reactants, pcs, exs)
    args = empty(defargs[:T])
    (rate, r_line) = exs[1:2]
    rxs = prune_reaction_line!(pcs, reactants, r_line); rate = prune_rate(rate)
    rxs = rxs isa Tuple ? tuple.(rate.args, rxs) : ((rate, rxs),)

    exs = exs[3:end]; empty!(args)
    ix = 1; while ix <= length(exs) 
        (!isa(exs[ix], Expr) || (exs[ix].head != :call)) && (ix += 1; continue)
        karg = (xi = findfirst(k -> exs[ix].args[2] ∈ k, prettynames); isnothing(xi) && (ix += 1; continue); xi)
        push!(args, karg => normalize_pcs!(pcs, exs[ix].args[3])); deleteat!(exs, ix)
    end
    args = merge(defargs[:T], args)

    append!(trans, tuple.(rxs, Ref(args))); trans
end

function replace_in_expr(expr, pairs...)
    dict = Dict(pairs...)

    prewalk(ex -> haskey(dict, ex) ? dict[ex] : ex, expr)
end

function normalize_pcs!(pcs, expr)
    postwalk(expr) do ex
        isexpr(ex, :macrocall) && macroname(ex) == :register &&
            (push!(pcs, deepcopy(ex)); ex.args[1] = Symbol("@", :take); ex.args = ex.args[1:3])
        if isexpr(ex, :macrocall) && macroname(ex) == :register
            r_sym = gensym()
            (push!(pcs, (ex_ = deepcopy(ex); insert!(ex_.args, 3, r_sym); ex_)); ex.args[1] = Symbol("@", :take); ex.args = [r_sym; ex.args[1:2]])
        end
        ex
    end
end

function get_reaction_line(expr)
    biarrow = nothing; prewalk(ex --> (ex ∈ double_arrows && (biarrow = ex); ex), expr)
    !isnothing(biarrow) ? [expr] : [replace_in_expr(expr, biarrow => :⟶), replace_in_expr(expr, biarrow => :⟵)]
end

function prune_reaction_line!(pcs, reactants, line)
    line isa Expr && (line.head == :-->) && (line = Expr(:call,:→,line.args[1],line.args[2]))
    if isexpr(line, :macrocall)
        lines = copy(line.args[3:end]); line.args = line.args[1:2]
        for l in lines
            biarrow = nothing; prewalk(ex -> (ex ∈ double_arrows && (biarrow = ex); ex), l)
            append!(line.args, isnothing(biarrow) ? [l] : [replace_in_expr(l, biarrow => :⟶), replace_in_expr(l, biarrow => :⟵)])
        end
        for i in 3:length(line.args)
            line.args[i] = isexpr(line.args[i], :tuple) ? Expr(:tuple, line.args[i].args[1], prune_reaction_line!(pcs, reactants, line.args[i].args[2])) :
                prune_reaction_line!(pcs, reactants, line.args[i])
        end
    elseif line isa Expr && line.args[1] ∈ union(fwd_arrows, bwd_arrows)
        line.args[2:3] = recursively_find_reactants!.(Ref(reactants), Ref(pcs), line.args[2:3])
    elseif line isa Expr && line.args[1] ∈ double_arrows
        biarrow = nothing; prewalk(ex -> (ex ∈ double_arrows && (biarrow = ex); ex), line)
        line = prune_reaction_line!.(Ref(pcs), Ref(reactants), (replace_in_expr(line, biarrow => :⟶), replace_in_expr(line, biarrow => :⟵)))
    end

    line
end

function recursively_find_reactants!(reactants, pcs, ex::SampleableValues)
    if typeof(ex) != Expr || isexpr(ex, :.) || (ex.head == :escape)   
        if (ex == 0 || in(ex, empty_set)) return :∅
        else push!(reactants, recursively_expand_dots(ex)) end
    elseif ex.args[1] == :*
        recursively_find_reactants!(reactants, pcs, ex.args[end])
        foreach(i -> ex.args[i] = normalize_pcs!(pcs, ex.args[i]), 2:(length(ex.args)-1))
    elseif ex.args[1] == :+
        for i = 2:length(ex.args) recursively_find_reactants!(reactants, pcs, ex.args[i]) end
    elseif isexpr(ex, :macrocall) && macroname(ex) == :choose
        for i in 3:length(ex.args)
            recursively_find_reactants!(reactants, pcs, isexpr(ex.args[i], :tuple) ? ex.args[i].args[2] : ex.args[i])
        end
    elseif isexpr(ex, :macrocall)
        recursively_find_reactants!(reactants, pcs, ex.args[3])
    else push!(reactants, underscorize(ex)) end

    ex
end