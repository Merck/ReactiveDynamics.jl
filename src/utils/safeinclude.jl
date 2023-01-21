using Test

export @safeinclude

macro safeinclude(args...)
    length(args) != 2 && error("invalid arguments to `@safeinclude`")
    name, ex = args
    quote
        path = pwd()
        td = mktempdir()
        cp(dirname($ex), td; force = true)
        cd(td)
        @testset $name begin include(joinpath(td, basename($ex))) end
        cd(path)
    end
end
