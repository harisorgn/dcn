using StaticArrays
using DifferentialEquations 
using DiffEqCallbacks
using ODEInterface 	

include("constants.jl")
include("read.jl")
include("auxiliary.jl")
include("strahler.jl")
include("comp.jl")
include("solve.jl")

data = read_genesis("./data")

reduced_comp_v = full_model(data) ;
#reduced_comp_v = strahler(data, strahler_thrs = 3) ;

comp_v = initialise_compartments(reduced_comp_v) ;

sol = model_solve(V0, Ca0, comp_v)
