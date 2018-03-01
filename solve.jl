

using DifferentialEquations 
using ODEInterface 	
#using Plots
#using PlotlyJS
#using Rsvg

function model_dynamics(du, u, p, t)

	idx_stride = 1 ;
	for i = 1 : length(comp_v)

		du[idx_stride : idx_stride + comp_v[i].n_eq - 1] = compartment_dynamics(comp_v[i], 
																			u[idx_stride : idx_stride + comp_v[i].n_eq - 1],
																			u[comp_v[i].V_connect_idx],
																			u[end] ) ;
		idx_stride += comp_v[i].n_eq ;

	end

	du[end] = 0.0 ;
	#println(du)
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

function pulse_initiate_event(u, t, integrator)
	t >= t_pulse_on && t <= t_pulse_on + t_pulse_dur ;
end

function pulse_initiate!(integrator)
	integrator.u[end] = Ie ;
end

function pulse_end_event(u, t, integrator)
	t >= t_pulse_on + t_pulse_dur ;
end

function pulse_end!(integrator)
	integrator.u[end] = 0.0 ;
end

function fig2(comp::soma_t, u::Array{Float64}, V_connect::Array{Float64})

	currents = [ g_fKdr_s * u[7]^4 * (u[1] - E_K) * comp.area + g_sKdr_s * u[8]^4 * (u[1] - E_K) * comp.area,
				 g_Sk_s * u[12] * (u[1] - E_K) * comp.area,
				 p_CaHVA_s * u[9]^3 * z_CaHVA^2 * F^2 * u[1] * (Ca_in - Ca_out * exp(-z_CaHVA * F * u[1] / (R * Temp))) / 
					(R * Temp * (1.0 - exp(-z_CaHVA * F * u[1] / (R * Temp)))) * comp.area,
				g_Naf_s * u[2]^3 * u[3] * (u[1] - E_Na) * comp.area + g_Nap_s * u[4]^3 * u[5] * (u[1] - E_Na) * comp.area,
				-1.0./(comp.Ra_connect_v[2]) * (u[1] - V_connect[2]),
				g_TNC_s * (u[1] - E_TNC) * comp.area,
				g_fKdr_s * u[7]^4 * (u[1] - E_K) * comp.area + g_sKdr_s * u[8]^4 * (u[1] - E_K) * comp.area - 1.0./(comp.Ra_connect_v[2]) * (u[1] - V_connect[2])] ;

	return currents
end

function model_solve(V0::Float64, Ca0::Float64)

	sum_eq = 0 ;
	for i = 1 : length(comp_v)

		sum_eq += comp_v[i].n_eq ;

	end

	u0 = zeros(Float64, sum_eq + 1) ;
	du = zeros(Float64, sum_eq + 1) ;

	idx_stride = 1 ;
	for i = 1 : length(comp_v)

		u0[idx_stride : idx_stride + comp_v[i].n_eq - 1] = initial_cond(comp_v[i], V0, Ca0) ;

		idx_stride += comp_v[i].n_eq ;

	end
	u0[end] = 0.0 ;
	
	tspan = (0.0, 2.5) ;

	prob = ODEProblem(model_dynamics, u0, tspan) ;

	cb_pulse_init = DiscreteCallback(pulse_initiate_event, pulse_initiate!, save_positions = (false, false)) ;
	cb_pulse_end = DiscreteCallback(pulse_end_event, pulse_end!, save_positions = (false, false)) ;
	cb = CallbackSet(cb_pulse_init, cb_pulse_end) ;
	
	sol = solve(prob, alg=:radau, maxiter = typemax(UInt128), callback = cb, dt = 5.0e-6, force_dtmin = true);#, reltol=1e-6, abstol=1e-6);
	
	#plotlyjs(size = (1000,1000))

	#plt = PlotlyJS.plot(sol, vars = (0,1), legend = false)
	#PlotlyJS.savefig(plt,"model_volt.png")

	#open("~/Documents/DCN/src/sol.txt") do f
	    
	#end
	file = open("sol.txt", "w+")
	writedlm(file, [sol.t sol[:,:]'])
	close(file)

	soma_currents = Array{Float64}(length(sol.t), 11) ;
	fig2_currents = Array{Float64}(length(sol.t), 7) ;
	for i = 1 : length(sol.t)

		soma_currents[i, :] = compartment_currents(comp_v[1], sol[1:13, i], sol[comp_v[1].V_connect_idx, i]) ;
		fig2_currents[i, :] = fig2(comp_v[1], sol[1:13, i], sol[comp_v[1].V_connect_idx, i]) ;

	end

	file = open("soma_cur.txt", "w+")
	writedlm(file, [sol.t soma_currents])
	close(file)

	file = open("fig2_cur.txt", "w+")
	writedlm(file, [sol.t fig2_currents])
	close(file)

	return
end

