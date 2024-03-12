using Plots

function plot_df(df::DataFrames.DataFrame, t_ix = 1)
    data = Matrix(df)
    t = @view data[:, t_ix]
    data_ = @view data[:, setdiff(1:size(data, 2), (t_ix,))]
    colnames = reshape(DataFrames.names(df)[setdiff(1:size(data, 2), (t_ix,))], 1, :)

    return Plots.plot(t, data_; labels = colnames, xlabel = "t")
end

# plot reduction
function AlgebraicAgents._draw(
    prob::ReactionNetworkProblem,
    vars = string.(prob.acs[:, :specName]);
    kwargs...,
)
    p = plot()
    for var in vars
        p = plot!(
            p,
            prob.sol[!, "t"],
            prob.sol[!, var];
            label = "$var",
            xlabel = "time",
            ylabel = "quantity",
            kwargs...,
        )
    end
    return p
end
