using ReactiveDynamics

# import a model, solve the problem
@import_network model.toml sir_acs
@assert @isdefined sir_acs
@assert isdefined(ReactiveDynamics, :tdecay)

prob = @problematize sir_acs
sol = @solve prob trajectories=20

@import_network "csv/model.csv" sir_acs_
@assert @isdefined sir_acs_
@assert isdefined(ReactiveDynamics, :foo)

prob_ = @problematize sir_acs_
sol_ = @solve prob_ trajectories=20

# export, import the solution
@export_solution sol
@import_solution sol.jld2

# export the same model (w/o registered functions)
@export_network sir_acs modell.toml
mkpath("csv_");
@export_network sir_acs "csv_/model.csv";

# load multiple models
@load_models models.txt

# export solution as a DataFrame, .csv
@export_solution_as_table sol
@export_solution_as_csv sol sol.csv

# another test
@import_network model2.toml sir_acs_2
