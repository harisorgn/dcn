
include("constants.jl")
include("strahler.jl")
include("comp.jl")
include("solve.jl")


# most distal compartment :
# p1b2b1b1b1b2b2b2b2[1] p1b2b1b1b1b2b2b2b2[0]   -0.62    181.5   127.06     0.58

#=
file = open("ddend.txt", "r")
ddend_dim = readdlm(file) ;
close(file)

file = open("pdend.txt", "r")
pdend_dim = readdlm(file) ;
close(file)

file = open("axis.txt", "r")
axis_dim = readdlm(file) ;
close(file)

ddend_len = 0 ;
for i = 1 : size(ddend_dim,1)-1

	ddend_len += sqrt((ddend_dim[end - i + 1,1] - ddend_dim[end - i,1])^2 + 
					  (ddend_dim[end - i + 1,2] - ddend_dim[end - i,2])^2 +
					  (ddend_dim[end - i + 1,3] - ddend_dim[end - i,3])^2) ;

end
#ddend_len *= 4.0 ;
ddend_diam = mean(ddend_dim[2:end,4]) ;

pdend_len = 0 ;
for i = 1 : size(pdend_dim,1)-1

	pdend_len += sqrt((pdend_dim[end - i + 1,1] - pdend_dim[end - i,1])^2 + 
					  (pdend_dim[end - i + 1,2] - pdend_dim[end - i,2])^2 +
					  (pdend_dim[end - i + 1,3] - pdend_dim[end - i,3])^2) ;

end

pdend_diam = mean(pdend_dim[2:end,4]) ;

axis_len = 0 ;
for i = 1 : size(axis_dim,1)-1

	axis_len += sqrt((axis_dim[end - i + 1,1] - axis_dim[end - i,1])^2 + 
					 (axis_dim[end - i + 1,2] - axis_dim[end - i,2])^2 +
					 (axis_dim[end - i + 1,3] - axis_dim[end - i,3])^2) ;

end

axis_diam = mean(axis_dim[2:end,4]) ;

axhill_len = 5.0 ;
axhill_diam = 4.75 ;

soma_len = 0.0 ;
soma_diam = 21.597 ;

comp_sym = [:s, :axhill, :axis, :pdend, :ddend] ;
#comp_sym = [:s, :pdend] ;
#comp_sym = [:s] ;

diam_v = Array{Float64}(length(comp_sym));
diam_v = [soma_diam * 1e-6, axhill_diam * 1e-6, axis_diam * 1e-6, pdend_diam * 1e-6, ddend_diam * 1e-6] ;
#diam_v = [soma_diam* 1e-6, pdend_diam * 1e-6] ;
#diam_v = [soma_diam* 1e-6] ;

len_v = Array{Float64}(length(comp_sym));
len_v = [soma_len * 1e-6, axhill_len * 1e-6, axis_len * 1e-6, pdend_len * 1e-6, ddend_len * 1e-6] ;
#len_v = [soma_len * 1e-6, pdend_len * 1e-6] ;
#len_v = [soma_len * 1e-6] ;

C = Array{Int64}(length(comp_sym), length(comp_sym))

C = [1 1 0 1 0 ;
	 1 1 1 0 0 ;
	 0 1 1 0 0 ;
	 1 0 0 1 1 ;
	 0 0 0 1 1 ] ;

#C = [1 1 ;
#	 1 1 ] ;
#C = [1] ;

comp_v = initialise_compartments(comp_sym, diam_v, len_v, C) ;

Profile.clear_malloc_data()

model_solve(V0, Ca0, comp_v)
=#



reduced_comp_v = strahler(3) ;

comp_v = init_comp(reduced_comp_v) ;

Profile.clear_malloc_data()

model_solve(V0, Ca0, comp_v)
