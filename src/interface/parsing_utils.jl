# parts of the code were taken from Catalyst.jl and adapted

# reaction arrows
empty_set = Set([:∅])
fwd_arrows = Set([:>, :→, :↣, :↦, :⇾, :⟶, :⟼, :⥟, :⥟, :⇀, :⇁, :⇒, :⟾])
bwd_arrows = Set([:<, :←, :↢, :↤, :⇽, :⟵, :⟻, :⥚, :⥞, :↼, :↽, :⇐, :⟽, Symbol("<--")])
double_arrows = Set([:↔, :⟷, :⇄, :⇆, :⇌, :⇋, :⇔, :⟺, Symbol("<-->")])

arrows = fwd_arrows ∪ bwd_arrows ∪ double_arrows ∪ [:-->]

# Returns the length of a expression tuple, or 1 if it is not an expression tuple (probably a  Symbol/Numerical).
function tup_leng(ex)
    (typeof(ex) == Expr && ex.head == :tuple) && (return length(ex.args))
    return 1
end

# Gets the ith element in a expression tuple, or returns the input itself if it is not an expression tuple (probably a  Symbol/Numerical).
function get_tup_arg(ex, i::Int)
    (tup_leng(ex) == 1) && (return ex)
    return ex.args[i]
end

# Makes dollars be interpolated in the caller's scope.
function esc_dollars!(ex)
    if ex isa Expr
        if ex.head == :$
            return esc(:($(ex.args[1])))
        else
            for i in eachindex(ex.args)
                ex.args[i] = esc_dollars!(ex.args[i])
            end
        end
    end

    ex
end
