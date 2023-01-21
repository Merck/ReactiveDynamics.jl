module ReactiveDynamics

using Reexport

# hierarchical co-integration framework
@reexport using AlgebraicAgents
import AlgebraicAgents: @crayon_str
# compact model declarations
@reexport using GeneratedExpressions
# expression manipulation
using MacroTools: prewalk, postwalk, striplines, isexpr

using UUIDs
using Distributions

# methods shared across the code
include("utils.jl")

# agent types
## convenience macro to set up interface
include("interface/params_interface.jl")
## reaction network (umbrella)
include("agents/network.jl")
### actions (stateful transitions, instantious events)
include("agents/actions.jl")
### species agent
include("agents/species.jl")
### sampled process (convenience wrap)
include("agents/sampleable.jl")
## DSL parsing
include("interface/prettynames.jl")
include("interface/parsing_utils.jl")
include("interface/parse_transitions.jl")
include("interface/parse_others.jl")
### turn expressions into functions of `(agent, network)`
include("interface/interpret.jl")
### macro interface
include("interface/macros.jl")
export @ReactionNetwork
export @transitions, @species, @events, @sampleables

end
