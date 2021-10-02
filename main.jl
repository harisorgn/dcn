using StaticArrays
using DifferentialEquations 
using DiffEqCallbacks
using ODEInterface 	
using BenchmarkTools

include("constants.jl")
include("auxiliary.jl")
include("strahler.jl")
include("comp.jl")
include("solve.jl")

#reduced_comp_v = strahler(3) ;
reduced_comp_v = full_model() ;

comp_v = init_comp(reduced_comp_v) ;

model_solve(V0, Ca0, comp_v)
