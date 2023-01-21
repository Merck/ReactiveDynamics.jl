using ReactiveDynamics

# import a model, solve the problem
@import_model model.toml sir_acs
@assert @isdefined sir_acs
@assert isdefined(ReactiveDynamics, :tdecay)

prob = @problematize sir_acs
sol = @solve prob trajectories=20

# export, import the solution
@export_solution sol
@import_solution sol.jld2

# export the same model (w/o registered functions)
@export_model sir_acs modell.toml

# load multiple models
@load_models models.txt

# export solution as a DataFrame, .csv
@export_as_table sol
@export_csv sol sol.csv

# another test
@import_model model2.toml sir_acs_2
