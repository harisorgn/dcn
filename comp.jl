
#include("constants.jl")

abstract type abstract_comp end

struct soma_t <: abstract_comp
	area ::Float64
	Ca_vol ::Float64
	Ra_connect_v ::Array{Float64}
	V_connect_idx ::Array{Int64}
	n_eq ::Int64
end

struct axon_hill_t <: abstract_comp
	area ::Float64
	Ca_vol ::Float64
	Ra_connect_v ::Array{Float64}
	V_connect_idx ::Array{Int64}
	n_eq ::Int64
end

struct axon_is_t <: abstract_comp
	area ::Float64
	Ca_vol ::Float64
	Ra_connect_v ::Array{Float64}
	V_connect_idx ::Array{Int64}
	n_eq ::Int64
end

struct prox_dend_t <: abstract_comp
	area ::Float64
	Ca_vol ::Float64
	Ra_connect_v ::Array{Float64}
	V_connect_idx ::Array{Int64}
	n_eq ::Int64
end

struct dist_dend_t <: abstract_comp
	area ::Float64
	Ca_vol ::Float64
	Ra_connect_v ::Array{Float64}
	V_connect_idx ::Array{Int64}
	n_eq ::Int64
end


function compartment_initialise(comp_sym::Array{Symbol}, diam_v::Array{Float64}, len_v::Array{Float64}, C::Array{Int64})

	comp_v = Array{abstract_comp}(length(comp_sym)) ;

	area_v = Array{Float64}(length(comp_sym)) ;
	Ca_vol_v = Array{Float64}(length(comp_sym)) ;
	Ra_v = Array{Float64}(length(comp_sym)) ;
	V_idx = Array{Int64}(length(comp_sym)) ;

	idx_stride = 1 ;
	for i = 1 : length(comp_sym)

		if comp_sym[i] == :s

			area_v[i] = pi * diam_v[i]^2 ;
			Ca_vol_v[i] = pi / 3.0 * (3.0 * diam_v[i]^2 * Ca_thick - 6.0 * diam_v[i] * Ca_thick^2 + 4.0 * Ca_thick^3) ;
			Ra_v[i] = 8.0 * RA / (pi * diam_v[i]) ;
			V_idx[i] = idx_stride ;

			idx_stride += n_eq_s ;

		elseif comp_sym[i] == :axhill

			area_v[i] = pi * diam_v[i] * len_v[i] ;
			Ca_vol_v[i] = pi * len_v[i] * (diam_v[i] * Ca_thick - Ca_thick^2) ;
			Ra_v[i] = 4.0 * RA * len_v[i] / (pi * diam_v[i]^2)
			V_idx[i] = idx_stride ;

			idx_stride += n_eq_ax ;

		elseif comp_sym[i] == :axis

			area_v[i] = pi * diam_v[i] * len_v[i] ;
			Ca_vol_v[i] = pi * len_v[i] * (diam_v[i] * Ca_thick - Ca_thick^2) ;
			Ra_v[i] = 4.0 * RA * len_v[i] / (pi * diam_v[i]^2)
			V_idx[i] = idx_stride ;

			idx_stride += n_eq_ax ;

		elseif comp_sym[i] == :pdend

			area_v[i] = pi * diam_v[i] * len_v[i] ;
			Ca_vol_v[i] = pi * len_v[i] * (diam_v[i] * Ca_thick - Ca_thick^2) ;
			Ra_v[i] = 4.0 * RA * len_v[i] / (pi * diam_v[i]^2)
			V_idx[i] = idx_stride ;

			idx_stride += n_eq_pd ;

		elseif comp_sym[i] == :ddend

			area_v[i] = pi * diam_v[i] * len_v[i] ;
			Ca_vol_v[i] = pi * len_v[i] * (diam_v[i] * Ca_thick - Ca_thick^2) ;
			Ra_v[i] = 4.0 * RA * len_v[i] / (pi * diam_v[i]^2)
			V_idx[i] = idx_stride ;

			idx_stride += n_eq_dd ;
			
		else

			println("Invalid compartment type \n")
			return
		end

	end

	for i = 1 : length(comp_sym)

		Ra_connect_idx = find(C[i,:]) ;
		Ra_connect_v = Array{Float64}(length(Ra_connect_idx)) ;

		for j = 1 : length(Ra_connect_idx)

			if Ra_connect_idx[j] <= i 

				Ra_connect_v[j] = Ra_v[i] ;

			else

				Ra_connect_v[j] = Ra_v[Ra_connect_idx[j]] ;

			end

		end

		if comp_sym[i] == :s

			comp_v[i] = soma_t(area_v[i], Ca_vol_v[i], Ra_connect_v, V_idx[find(C[i,:])], n_eq_s) ;

		elseif comp_sym[i] == :axhill

			comp_v[i] = axon_hill_t(area_v[i], Ca_vol_v[i], Ra_connect_v, V_idx[find(C[i,:])], n_eq_ax) ;

		elseif comp_sym[i] == :axis

			comp_v[i] = axon_is_t(area_v[i], Ca_vol_v[i], Ra_connect_v, V_idx[find(C[i,:])], n_eq_ax) ;

		elseif comp_sym[i] == :pdend

			comp_v[i] = prox_dend_t(area_v[i], Ca_vol_v[i], Ra_connect_v, V_idx[find(C[i,:])], n_eq_pd) ;

		elseif comp_sym[i] == :ddend

			comp_v[i] = dist_dend_t(area_v[i], Ca_vol_v[i], Ra_connect_v, V_idx[find(C[i,:])], n_eq_dd) ;

		end

	end
	return comp_v
end


function compartment_dynamics(comp::soma_t, u::Array{Float64}, V_connect::Array{Float64})

	du = Array{Float64}(n_eq_s) ;
	V_rep_v = repmat([u[1]], length(V_connect))

	du[1] = (I - g_Naf_s * u[2]^3 * u[3] * (u[1] - E_Na) * comp.area - 
				 g_Nap_s * u[4]^3 * u[5] * (u[1] - E_Na) * comp.area -
				 g_TNC_s * (u[1] - E_TNC) * comp.area -
				 g_h_s * u[6]^2 * (u[1] - E_h) * comp.area - 
				 g_fKdr_s * u[7]^4 * (u[1] - E_K) * comp.area - 
				 g_sKdr_s * u[8]^4 * (u[1] - E_K) * comp.area -
				 p_CaHVA_s * u[9]^3 * z_CaHVA^2 * F^2 * u[1] * (Ca_in - Ca_out * exp(-z_CaHVA * F * u[1] / (R * Temp))) / 
						(R * Temp * (1.0 - exp(-z_CaHVA * F * u[1] / (R * Temp)))) * comp.area -
				 g_CaLVA_s * u[10]^2 * u[11] * (u[1] - E_Ca) * comp.area - 
				 g_Sk_s * u[12] * (u[1] - E_K) * comp.area -
				 (u[1] - Em) / (RM_s / comp.area) -
				 1.0./(comp.Ra_connect_v') * (V_rep_v - V_connect)) / (CM * comp.area);

	du[2] = (m_Naf(u[1]) - u[2]) / t_m_Naf(u[1]) ;
	du[3] = (h_Naf(u[1]) - u[3]) / t_h_Naf(u[1]) ;
	du[4] = (m_Nap(u[1]) - u[4]) / t_m_Nap ;
	du[5] = (h_Nap(u[1]) - u[5]) / t_h_Nap(u[1]) ;
	du[6] = (m_h(u[1]) - u[6]) / t_m_h ;
	du[7] = (m_fKdr(u[1]) - u[7]) / t_m_fKdr(u[1]) ;
	du[8] = (m_sKdr(u[1]) - u[8]) / t_m_sKdr(u[1]) ;
	du[9] = (m_CaHVA(u[1]) - u[9]) / t_m_CaHVA(u[1]) ;
	du[10] = (m_CaLVA(u[1]) - u[10]) / t_m_CaLVA(u[1]) ;
	du[11] = (h_CaLVA(u[1]) - u[11]) / t_h_CaLVA(u[1]) ;
	du[12] = (z_Sk(u[13]) - u[12]) / t_z_Sk(u[13]) ;
	du[13] = (k_Ca_s / comp.Ca_vol) * p_CaHVA_s * u[9]^3 * z_CaHVA^2 * F^2 * u[1] * (Ca_in - Ca_out * exp(-z_CaHVA * F * u[1] / (R * Temp))) / 
			(R * Temp * (1.0 - exp(-z_CaHVA * F * u[1] / (R * Temp)))) * comp.area - (u[13] - Ca_base) / t_Ca ;

	return du

end

function compartment_dynamics(comp::axon_hill_t, u::Array{Float64}, V_connect::Array{Float64})

	du = Array{Float64}(n_eq_ax) ;
	V_rep_v = repmat([u[1]], length(V_connect))

	du[1] = (  - g_Naf_ax * u[2]^3 * u[3] * (u[1] - E_Na) * comp.area - 
				 g_TNC_ax * (u[1] - E_TNC) * comp.area -
				 g_fKdr_ax * u[4]^4 * (u[1] - E_K) * comp.area - 
				 g_sKdr_ax * u[5]^4 * (u[1] - E_K) * comp.area -
				 (u[1] - Em) / (RM_ax / comp.area) -
				 1.0./(comp.Ra_connect_v') * (V_rep_v - V_connect)) / (CM * comp.area);

	du[2] = (m_Naf(u[1]) - u[2]) / t_m_Naf(u[1]) ;
	du[3] = (h_Naf(u[1]) - u[3]) / t_h_Naf(u[1]) ;
	du[4] = (m_fKdr(u[1]) - u[4]) / t_m_fKdr(u[1]) ;
	du[5] = (m_sKdr(u[1]) - u[5]) / t_m_sKdr(u[1]) ;
	
	return du

end

function compartment_dynamics(comp::axon_is_t, u::Array{Float64}, V_connect::Array{Float64})

	du = Array{Float64}(n_eq_ax) ;
	V_rep_v = repmat([u[1]], length(V_connect))

	du[1] = (  - g_Naf_ax * u[2]^3 * u[3] * (u[1] - E_Na) * comp.area - 
				 g_TNC_ax * (u[1] - E_TNC) * comp.area -
				 g_fKdr_ax * u[4]^4 * (u[1] - E_K) * comp.area - 
				 g_sKdr_ax * u[5]^4 * (u[1] - E_K) * comp.area -
				 (u[1] - Em) / (RM_ax / comp.area) -
				 1.0./(comp.Ra_connect_v') * (V_rep_v - V_connect)) / (CM * comp.area);

	du[2] = (m_Naf(u[1]) - u[2]) / t_m_Naf(u[1]) ;
	du[3] = (h_Naf(u[1]) - u[3]) / t_h_Naf(u[1]) ;
	du[4] = (m_fKdr(u[1]) - u[4]) / t_m_fKdr(u[1]) ;
	du[5] = (m_sKdr(u[1]) - u[5]) / t_m_sKdr(u[1]) ;
	
	return du

end

function compartment_dynamics(comp::prox_dend_t, u::Array{Float64}, V_connect::Array{Float64})

	du = Array{Float64}(n_eq_pd) ;
	V_rep_v = repmat([u[1]], length(V_connect))

	du[1] = (  - g_Naf_pd * u[2]^3 * u[3] * (u[1] - E_Na) * comp.area - 
				 g_TNC_pd * (u[1] - E_TNC) * comp.area -
				 g_h_pd * u[4]^2 * (u[1] - E_h) * comp.area - 
				 g_fKdr_pd * u[5]^4 * (u[1] - E_K) * comp.area - 
				 g_sKdr_pd * u[6]^4 * (u[1] - E_K) * comp.area -
				 p_CaHVA_pd * u[7]^3 * z_CaHVA^2 * F^2 * u[1] * (Ca_in - Ca_out * exp(-z_CaHVA * F * u[1] / (R * Temp))) / 
						(R * Temp * (1.0 - exp(-z_CaHVA * F * u[1] / (R * Temp)))) * comp.area -
				 g_CaLVA_pd * u[8]^2 * u[9] * (u[1] - E_Ca) * comp.area - 
				 g_Sk_pd * u[10] * (u[1] - E_K) * comp.area -
				 (u[1] - Em) / (RM_d / comp.area) -
				 1.0./(comp.Ra_connect_v') * (V_rep_v - V_connect)) / (CM * comp.area);

	du[2] = (m_Naf(u[1]) - u[2]) / t_m_Naf(u[1]) ;
	du[3] = (h_Naf(u[1]) - u[3]) / t_h_Naf(u[1]) ;
	du[4] = (m_h(u[1]) - u[4]) / t_m_h ;
	du[5] = (m_fKdr(u[1]) - u[5]) / t_m_fKdr(u[1]) ;
	du[6] = (m_sKdr(u[1]) - u[6]) / t_m_sKdr(u[1]) ;
	du[7] = (m_CaHVA(u[1]) - u[7]) / t_m_CaHVA(u[1]) ;
	du[8] = (m_CaLVA(u[1]) - u[8]) / t_m_CaLVA(u[1]) ;
	du[9] = (h_CaLVA(u[1]) - u[9]) / t_h_CaLVA(u[1]) ;
	du[10] = (z_Sk(u[11]) - u[10]) / t_z_Sk(u[11]) ;
	du[11] = (k_Ca_d / comp.Ca_vol) * p_CaHVA_s * u[7]^3 * z_CaHVA^2 * F^2 * u[1] * (Ca_in - Ca_out * exp(-z_CaHVA * F * u[1] / (R * Temp))) / 
			(R * Temp * (1.0 - exp(-z_CaHVA * F * u[1] / (R * Temp)))) * comp.area - (u[11] - Ca_base) / t_Ca ;

	return du

end

function compartment_dynamics(comp::dist_dend_t, u::Array{Float64}, V_connect::Array{Float64})

	du = Array{Float64}(n_eq_dd) ;
	V_rep_v = repmat([u[1]], length(V_connect))

	du[1] = (  - g_h_dd * u[2]^2 * (u[1] - E_h) * comp.area - 
				 p_CaHVA_dd * u[3]^3 * z_CaHVA^2 * F^2 * u[1] * (Ca_in - Ca_out * exp(-z_CaHVA * F * u[1] / (R * Temp))) / 
						(R * Temp * (1.0 - exp(-z_CaHVA * F * u[1] / (R * Temp)))) * comp.area -
				 g_CaLVA_dd * u[4]^2 * u[5] * (u[1] - E_Ca) * comp.area - 
				 g_Sk_dd * u[6] * (u[1] - E_K) * comp.area -
				 (u[1] - Em) / (RM_d / comp.area) -
				 1.0./(comp.Ra_connect_v') * (V_rep_v - V_connect)) / (CM * comp.area);

	du[2] = (m_h(u[1]) - u[2]) / t_m_h ;
	du[3] = (m_CaHVA(u[1]) - u[3]) / t_m_CaHVA(u[1]) ;
	du[4] = (m_CaLVA(u[1]) - u[4]) / t_m_CaLVA(u[1]) ;
	du[5] = (h_CaLVA(u[1]) - u[5]) / t_h_CaLVA(u[1]) ;
	du[6] = (z_Sk(u[7]) - u[6]) / t_z_Sk(u[7]) ;
	du[7] = (k_Ca_d / comp.Ca_vol) * p_CaHVA_s * u[3]^3 * z_CaHVA^2 * F^2 * u[1] * (Ca_in - Ca_out * exp(-z_CaHVA * F * u[1] / (R * Temp))) / 
			(R * Temp * (1.0 - exp(-z_CaHVA * F * u[1] / (R * Temp)))) * comp.area - (u[7] - Ca_base) / t_Ca ;

	return du

end