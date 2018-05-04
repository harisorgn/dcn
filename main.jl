
include("constants.jl")
include("strahler.jl")
include("comp.jl")
include("solve.jl")

#reduced_comp_v = strahler(3) ;
reduced_comp_v = full_model() ;

comp_v = init_comp(reduced_comp_v) ;

#=
comp_s = Array{abstract_comp, 1}() ;

area = pi * reduced_comp_v[1].d^2.0 ;
Ca_vol = pi / 3.0 * (3.0 * reduced_comp_v[1].d^2.0 * Ca_thick - 
						6.0 * reduced_comp_v[1].d * Ca_thick^2.0 + 
						4.0 * Ca_thick^3.0) ;

push!(comp_s, soma_t(area, Ca_vol, [], [], n_eq_s, 1.0)) ; 
=#

#Profile.clear_malloc_data()

model_solve(V0, Ca0, comp_v)
