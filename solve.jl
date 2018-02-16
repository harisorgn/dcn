

using DifferentialEquations 
using ODEInterface 	
#using Plots
#using PlotlyJS
#using Rsvg

function model_dynamics(du, u, I, t)

	idx_stride = 1 ;
	for i = 1 : length(comp_v)

		du[idx_stride : idx_stride + comp_v[i].n_eq - 1] = compartment_dynamics(comp_v[i], 
																			u[idx_stride : idx_stride + comp_v[i].n_eq - 1],
																			u[comp_v[i].V_connect_idx]) ;
		idx_stride += comp_v[i].n_eq ;

	end

	#println(u)
	return du 

end

function initial_cond(comp::soma_t, V0::Float64, Ca0::Float64)

	u0 = [V0 ;
		m_Naf(V0) ; 
		h_Naf(V0) ; 
		m_Nap(V0) ; 
		h_Nap(V0) ; 
		m_h(V0) ; 
		m_fKdr(V0) ; 
		m_sKdr(V0) ; 
		m_CaHVA(V0) ; 
		m_CaLVA(V0) ;
		h_CaLVA(V0) ;
		z_Sk(Ca0) ;
		Ca0 ]

end

function initial_cond(comp::axon_hill_t, V0::Float64, Ca0::Float64)

	u0 = [V0 ;
		m_Naf(V0) ; 
		h_Naf(V0) ; 
		m_fKdr(V0) ; 
		m_sKdr(V0)] ; 

end

function initial_cond(comp::axon_is_t, V0::Float64, Ca0::Float64)

	u0 = [V0 ;
		m_Naf(V0) ; 
		h_Naf(V0) ; 
		m_fKdr(V0) ; 
		m_sKdr(V0)] ; 

end

function initial_cond(comp::prox_dend_t, V0::Float64, Ca0::Float64)

	u0 = [V0 ;
		m_Naf(V0) ; 
		h_Naf(V0) ; 
		m_h(V0) ; 
		m_fKdr(V0) ; 
		m_sKdr(V0) ; 
		m_CaHVA(V0) ; 
		m_CaLVA(V0) ;
		h_CaLVA(V0) ;
		z_Sk(Ca0) ;
		Ca0 ]

end

function initial_cond(comp::dist_dend_t, V0::Float64, Ca0::Float64)

	u0 = [V0 ;
		m_h(V0) ; 
		m_CaHVA(V0) ; 
		m_CaLVA(V0) ;
		h_CaLVA(V0) ;
		z_Sk(Ca0) ;
		Ca0 ]

end

function model_solve(I::Float64, V0::Float64, Ca0::Float64)

	sum_eq = 0 ;
	for i = 1 : length(comp_v)

		sum_eq += comp_v[i].n_eq ;

	end

	u0 = zeros(Float64, sum_eq) ;
	du = zeros(Float64, sum_eq) ;

	idx_stride = 1 ;
	for i = 1 : length(comp_v)

		u0[idx_stride : idx_stride + comp_v[i].n_eq - 1] = initial_cond(comp_v[i], V0, Ca0) ;

		idx_stride += comp_v[i].n_eq ;

	end

	tspan = (0.0, 0.4) ;

	prob = ODEProblem(model_dynamics, u0, tspan, I) ;
	sol = solve(prob, alg=:radau, dt = 1.0e-6, reltol=1e-6, abstol=1e-6);

	#plotlyjs(size = (1000,1000))

	#plt = PlotlyJS.plot(sol, vars = (0,1), legend = false)
	#PlotlyJS.savefig(plt,"model_volt.png")

	#open("~/Documents/DCN/src/sol.txt") do f
	    
	#end
	file = open("sol.txt", "w+")
	writedlm(file, [sol.t sol[:,:]'])
	close(file)

	return
end

