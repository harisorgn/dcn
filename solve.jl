
function model_dynamics(du, u, p, t)

	idx_stride = 1 ;
	for i = 1 : length(p)
		
		du[idx_stride : idx_stride + p[i].n_eq - 1] = compartment_dynamics(p[i], 
																			u[idx_stride : idx_stride + p[i].n_eq - 1],
																			u[p[i].V_connect_idx],
																			u[end] ) ;
		idx_stride += p[i].n_eq ;

	end

	du[end] = 0.0 ;
	return du 

end

initial_cond(comp::soma_t, V0::Float64, Ca0::Float64) = [
	V0 ;m_Naf(V0) ; h_Naf(V0) ; m_Nap(V0) ; h_Nap(V0) ; m_h(V0) ; m_fKdr(V0) ; m_sKdr(V0) ; 
	m_CaHVA(V0) ; m_CaLVA(V0) ;h_CaLVA(V0) ;z_Sk(Ca0) ;Ca0 ]

initial_cond(comp::axon_hill_t, V0::Float64, Ca0::Float64) = [
	V0 ;m_Naf(V0) ; h_Naf(V0) ; m_fKdr(V0) ; m_sKdr(V0)] 

initial_cond(comp::axon_is_t, V0::Float64, Ca0::Float64) = [
	V0 ;m_Naf(V0) ; h_Naf(V0) ; m_fKdr(V0) ; m_sKdr(V0)]  

initial_cond(comp::axon_in_t, V0::Float64, Ca0::Float64) = [V0]  

initial_cond(comp::prox_dend_t, V0::Float64, Ca0::Float64) = [
	V0 ;m_Naf(V0) ; h_Naf(V0) ; m_h(V0) ; m_fKdr(V0) ; m_sKdr(V0) ; m_CaHVA(V0) ; m_CaLVA(V0) ;
	h_CaLVA(V0) ;z_Sk(Ca0) ;Ca0 ]

initial_cond(comp::dist_dend_t, V0::Float64, Ca0::Float64) = [
	V0 ;m_h(V0) ; m_CaHVA(V0) ; m_CaLVA(V0) ;h_CaLVA(V0) ;z_Sk(Ca0) ;Ca0 ]



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

function model_solve(V0::Float64, Ca0::Float64, comp_v::Array{abstract_comp})

	sum_eq = 0 ;
	for i = 1 : length(comp_v)

		sum_eq += comp_v[i].n_eq ;

	end

	u0 = zeros(Float64, sum_eq + 1) ;

	idx_stride = 1 ;
	for i = 1 : length(comp_v)

		u0[idx_stride : idx_stride + comp_v[i].n_eq - 1] = initial_cond(comp_v[i], V0, Ca0) ;

		idx_stride += comp_v[i].n_eq ;

	end
	u0[end] = 0.0 ;
	
	tspan = (0.0, 2.5) ;

	prob = ODEProblem(model_dynamics, u0, tspan, comp_v) ;

	cb_pulse_init = DiscreteCallback(pulse_initiate_event, pulse_initiate!, save_positions = (false, false)) ;
	cb_pulse_end = DiscreteCallback(pulse_end_event, pulse_end!, save_positions = (false, false)) ;
	cb = CallbackSet(cb_pulse_init, cb_pulse_end) ;
	
	sol = @time solve(prob, alg = :radau, callback = cb, dt = 5.0e-6, 
					save_everystep = true, timeseries_steps = 100 , save_idxs = [1],
					maxiters = typemax(UInt128));

	file = open("sol.bin", "w+")
	#writedlm(file, Float32([sol.t sol[:,:]']))
	writedlm(file, [sol.t sol[1,:]])
	close(file)

	return sol
end

