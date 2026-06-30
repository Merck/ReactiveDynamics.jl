# RDArrowExt — the byte-faithful Arrow writer for rectangular export artifacts (ADR 0013 §C /
# CONTRACT §14.3).
#
# Loaded automatically when `Arrow` is available (declared in Project.toml `[weakdeps]`/`[extensions]`).
# The export bundle (src/export.jl) writes the always-available JSON + CSV core unconditionally and
# calls `ReactiveDynamics._arrow_write(path, table)` for the faithful columnar `.arrow` siblings of
# the rectangular artifacts (trajectory, ledger, ensemble summary). The core defines `_arrow_write`
# as a no-op that returns `nothing` (so an Arrow-less run silently skips the `.arrow` files); this
# extension overrides it to actually write, so Arrow output appears exactly when the user has Arrow.
module RDArrowExt

using ReactiveDynamics
using Arrow

# Override the core's no-op: write any Tables.jl-compatible table (a DataFrame) to `path` as Arrow.
function ReactiveDynamics._arrow_write(path::AbstractString, table)
    Arrow.write(path, table)
    return path
end

end # module RDArrowExt
