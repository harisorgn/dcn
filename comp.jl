
abstract type abstract_comp end

struct soma_t <: abstract_comp
	area ::Float64
	Ca_vol ::Float64
	Ra_connect_v ::Array{Float64}
	V_connect_idx ::Array{Int64}
	n_eq ::Int64
	reduce_fact ::Float64
end

struct axon_hill_t <: abstract_comp
	area ::Float64
	Ca_vol ::Float64
	Ra_connect_v ::Array{Float64}
	V_connect_idx ::Array{Int64}
	n_eq ::Int64
	reduce_fact ::Float64
end

struct axon_is_t <: abstract_comp
	area ::Float64
	Ca_vol ::Float64
	Ra_connect_v ::Array{Float64}
	V_connect_idx ::Array{Int64}
	n_eq ::Int64
	reduce_fact ::Float64
end

struct axon_in_t <: abstract_comp
	area ::Float64
	Ca_vol ::Float64
	Ra_connect_v ::Array{Float64}
	V_connect_idx ::Array{Int64}
	n_eq ::Int64
	reduce_fact ::Float64
end

struct prox_dend_t <: abstract_comp
	area ::Float64
	Ca_vol ::Float64
	Ra_connect_v ::Array{Float64}
	V_connect_idx ::Array{Int64}
	n_eq ::Int64
	reduce_fact ::Float64
end

struct dist_dend_t <: abstract_comp
	area ::Float64
	Ca_vol ::Float64
	Ra_connect_v ::Array{Float64}
	V_connect_idx ::Array{Int64}
	n_eq ::Int64
	reduce_fact ::Float64
end

function init_comp(reduced_comp_v ::Array{reduced_comp_t, 1})

	comp_v = Array{abstract_comp}(length(reduced_comp_v)) ;

	area_v = Array{Float64}(length(reduced_comp_v)) ;
	Ca_vol_v = Array{Float64}(length(reduced_comp_v)) ;
	V_idx = Array{Int64}(length(reduced_comp_v)) ;

	idx_stride = 1 ;
	for i = 1 : length(reduced_comp_v)

		if reduced_comp_v[i].ctype == :s

			area_v[i] = pi * reduced_comp_v[i].d^2.0 ;
			Ca_vol_v[i] = pi / 3.0 * (3.0 * reduced_comp_v[i].d^2.0 * Ca_thick - 
									6.0 * reduced_comp_v[i].d * Ca_thick^2.0 + 
									4.0 * Ca_thick^3.0) ;
			V_idx[i] = idx_stride ;

			idx_stride += n_eq_s ;

		elseif reduced_comp_v[i].ctype == :axhill

			area_v[i] = pi * reduced_comp_v[i].d * reduced_comp_v[i].l ;
			Ca_vol_v[i] = pi * reduced_comp_v[i].l * (reduced_comp_v[i].d * Ca_thick - Ca_thick^2.0) ;
			V_idx[i] = idx_stride ;

			idx_stride += n_eq_ax ;

		elseif reduced_comp_v[i].ctype == :axis

			area_v[i] = pi * reduced_comp_v[i].d * reduced_comp_v[i].l ;
			Ca_vol_v[i] = pi * reduced_comp_v[i].l * (reduced_comp_v[i].d * Ca_thick - Ca_thick^2.0) ;
			V_idx[i] = idx_stride ;

			idx_stride += n_eq_ax ;

		elseif reduced_comp_v[i].ctype == :axin

			area_v[i] = pi * reduced_comp_v[i].d * reduced_comp_v[i].l ;
			Ca_vol_v[i] = pi * reduced_comp_v[i].l * (reduced_comp_v[i].d * Ca_thick - Ca_thick^2.0) ;
			V_idx[i] = idx_stride ;

			idx_stride += n_eq_axin ;

		elseif reduced_comp_v[i].ctype == :pdend

			area_v[i] = pi * reduced_comp_v[i].d * reduced_comp_v[i].l ;
			Ca_vol_v[i] = pi * reduced_comp_v[i].l * (reduced_comp_v[i].d * Ca_thick - Ca_thick^2.0) ;
			V_idx[i] = idx_stride ;

			idx_stride += n_eq_pd ;

		elseif reduced_comp_v[i].ctype == :ddend

			area_v[i] = pi * reduced_comp_v[i].d * reduced_comp_v[i].l ;
			Ca_vol_v[i] = pi * reduced_comp_v[i].l * (reduced_comp_v[i].d * Ca_thick - Ca_thick^2.0) ;
			V_idx[i] = idx_stride ;

			idx_stride += n_eq_dd ;
			
		else
			println("Invalid compartment type \n")
			return
		end
	end

	for i = 1 : length(reduced_comp_v)
		Ra_connect_v = Array{Float64,1}() ;
		V_connect_idx = Array{Int64,1}() ;

		for j = 1 : length(reduced_comp_v)
			if reduced_comp_v[j].connect_to == reduced_comp_v[i].id
				push!(Ra_connect_v, reduced_comp_v[j].Ra / reduced_comp_v[j].mrg_n_comp) ;
				push!(V_connect_idx, V_idx[j]) ;
			elseif reduced_comp_v[j].id == reduced_comp_v[i].connect_to
				push!(Ra_connect_v, reduced_comp_v[i].Ra / reduced_comp_v[i].mrg_n_comp) ;
				push!(V_connect_idx, V_idx[j]) ;
			end
		end
		
		if reduced_comp_v[i].ctype == :s

			comp_v[i] = soma_t(area_v[i], Ca_vol_v[i], Ra_connect_v, V_connect_idx, n_eq_s, 
								reduced_comp_v[i].area / (pi * reduced_comp_v[i].d^2.0)) ;
			#comp_v[i] = soma_t(area_v[i], Ca_vol_v[i], Ra_connect_v, V_connect_idx, n_eq_s, 
			#					reduced_comp_v[i].full_area_sum / (pi * reduced_comp_v[i].d^2.0)) ;

		elseif reduced_comp_v[i].ctype == :axhill

			comp_v[i] = axon_hill_t(area_v[i], Ca_vol_v[i], Ra_connect_v, V_connect_idx, n_eq_ax,
									reduced_comp_v[i].area / (pi * reduced_comp_v[i].l * reduced_comp_v[i].d)) ;
			#comp_v[i] = axon_hill_t(area_v[i], Ca_vol_v[i], Ra_connect_v, V_connect_idx, n_eq_ax,
			#						reduced_comp_v[i].full_area_sum / (pi * reduced_comp_v[i].l * reduced_comp_v[i].d)) ;

		elseif reduced_comp_v[i].ctype == :axis

			comp_v[i] = axon_is_t(area_v[i], Ca_vol_v[i], Ra_connect_v, V_connect_idx, n_eq_ax,
								reduced_comp_v[i].area / (pi * reduced_comp_v[i].l * reduced_comp_v[i].d)) ;
			#comp_v[i] = axon_is_t(area_v[i], Ca_vol_v[i], Ra_connect_v, V_connect_idx, n_eq_ax,
			#					reduced_comp_v[i].full_area_sum / (pi * reduced_comp_v[i].l * reduced_comp_v[i].d)) ;

		elseif reduced_comp_v[i].ctype == :axin

			comp_v[i] = axon_in_t(area_v[i], Ca_vol_v[i], Ra_connect_v, V_connect_idx, n_eq_axin,
								reduced_comp_v[i].area / (pi * reduced_comp_v[i].l * reduced_comp_v[i].d)) ;
			#comp_v[i] = axon_in_t(area_v[i], Ca_vol_v[i], Ra_connect_v, V_connect_idx, n_eq_axin,
			#					reduced_comp_v[i].full_area_sum / (pi * reduced_comp_v[i].l * reduced_comp_v[i].d)) ;

		elseif reduced_comp_v[i].ctype == :pdend

			comp_v[i] = prox_dend_t(area_v[i], Ca_vol_v[i], Ra_connect_v, V_connect_idx, n_eq_pd,
									reduced_comp_v[i].area / (pi * reduced_comp_v[i].l * reduced_comp_v[i].d)) ;
			#comp_v[i] = prox_dend_t(area_v[i], Ca_vol_v[i], Ra_connect_v, V_connect_idx, n_eq_pd,
			#						reduced_comp_v[i].full_area_sum / (pi * reduced_comp_v[i].l * reduced_comp_v[i].d)) ;

		elseif reduced_comp_v[i].ctype == :ddend

			comp_v[i] = dist_dend_t(area_v[i], Ca_vol_v[i], Ra_connect_v, V_connect_idx, n_eq_dd,
									reduced_comp_v[i].area / (pi * reduced_comp_v[i].l * reduced_comp_v[i].d)) ;
			#comp_v[i] = dist_dend_t(area_v[i], Ca_vol_v[i], Ra_connect_v, V_connect_idx, n_eq_dd,
			#						reduced_comp_v[i].full_area_sum / (pi * reduced_comp_v[i].l * reduced_comp_v[i].d)) ;

		end		
	end
	return comp_v
end

function compartment_dynamics(comp::soma_t, u::Array{Float64}, V_connect::Array{Float64}, Ie::Float64)
	
	I_connect = 0.0 ;
	for i = 1 : length(comp.V_connect_idx)
		I_connect += 1.0/comp.Ra_connect_v[i] * (u[1] - V_connect[i]) ;
	end

	dv = (Ie - g_Naf_s * u[2]^3 * u[3] * (u[1] - E_Na) * comp.area * comp.reduce_fact - 
				 g_Nap_s * u[4]^3 * u[5] * (u[1] - E_Na) * comp.area * comp.reduce_fact -
				 g_TNC_s * (u[1] - E_TNC) * comp.area * comp.reduce_fact -
				 g_h_s * u[6]^2 * (u[1] - E_h) * comp.area * comp.reduce_fact - 
				 g_fKdr_s * u[7]^4 * (u[1] - E_K) * comp.area * comp.reduce_fact - 
				 g_sKdr_s * u[8]^4 * (u[1] - E_K) * comp.area * comp.reduce_fact -
				 p_CaHVA_s * u[9]^3 * z_CaHVA^2 * F^2 * u[1] * (u[13] - Ca_out * exp(-z_CaHVA * F * u[1] / (R * Temp))) / 
						(R * Temp * (1.0 - exp(-z_CaHVA * F * u[1] / (R * Temp)))) * comp.area * comp.reduce_fact-
				 g_CaLVA_s * u[10]^2 * u[11] * (u[1] - E_Ca) * comp.area * comp.reduce_fact - 
				 g_Sk_s * u[12] * (u[1] - E_K) * comp.area * comp.reduce_fact-
				 (u[1] - Em) / (RM_s / (comp.reduce_fact * comp.area)) -
				 I_connect) / (CM * comp.reduce_fact * comp.area);

	dm_Naf = (m_Naf(u[1]) - u[2]) / t_m_Naf(u[1]) ;
	dh_Naf = (h_Naf(u[1]) - u[3]) / t_h_Naf(u[1]) ;
	dm_Nap = (m_Nap(u[1]) - u[4]) / t_m_Nap ;
	dh_Nap = (h_Nap(u[1]) - u[5]) / t_h_Nap(u[1]) ;
	dm_h = (m_h(u[1]) - u[6]) / t_m_h ;
	dm_fKdr = (m_fKdr(u[1]) - u[7]) / t_m_fKdr(u[1]) ;
	dm_sKdr = (m_sKdr(u[1]) - u[8]) / t_m_sKdr(u[1]) ;
	dm_CaHVA = (m_CaHVA(u[1]) - u[9]) / t_m_CaHVA(u[1]) ;
	dm_CaLVA = (m_CaLVA(u[1]) - u[10]) / t_m_CaLVA(u[1]) ;
	dh_CaLVA = (h_CaLVA(u[1]) - u[11]) / t_h_CaLVA(u[1]) ;
	dz_Sk = (z_Sk(u[13]) - u[12]) / t_z_Sk(u[13]) ;
	dCa = -(k_Ca_s / comp.Ca_vol) * p_CaHVA_s * u[9]^3 * z_CaHVA^2 * F^2 * u[1] * (u[13] - Ca_out * exp(-z_CaHVA * F * u[1] / (R * Temp))) / 
			(R * Temp * (1.0 - exp(-z_CaHVA * F * u[1] / (R * Temp)))) * comp.area * comp.reduce_fact - 
			(u[13] - Ca_base) / t_Ca ;

	return @SVector [dv, dm_Naf, dh_Naf, dm_Nap, dh_Nap, dm_h, dm_fKdr, dm_sKdr, dm_CaHVA, dm_CaLVA, dh_CaLVA, dz_Sk, dCa]
end

function compartment_dynamics(comp::axon_hill_t, u::Array{Float64}, V_connect::Array{Float64}, Ie::Float64)
	
	I_connect = 0.0 ;
	for i = 1 : length(comp.V_connect_idx)
		I_connect += 1.0/comp.Ra_connect_v[i] * (u[1] - V_connect[i]) ;
	end

	dv = (  - g_Naf_ax * u[2]^3 * u[3] * (u[1] - E_Na) * comp.area * comp.reduce_fact - 
				 g_TNC_ax * (u[1] - E_TNC) * comp.area * comp.reduce_fact -
				 g_fKdr_ax * u[4]^4 * (u[1] - E_K) * comp.area * comp.reduce_fact - 
				 g_sKdr_ax * u[5]^4 * (u[1] - E_K) * comp.area * comp.reduce_fact -
				 (u[1] - Em) / (RM_ax / (comp.reduce_fact * comp.area)) -
				 I_connect) / (CM * comp.reduce_fact * comp.area);

	dm_Naf = (m_Naf(u[1]) - u[2]) / t_m_Naf(u[1]) ;
	dh_Naf = (h_Naf(u[1]) - u[3]) / t_h_Naf(u[1]) ;
	dm_fKdr = (m_fKdr(u[1]) - u[4]) / t_m_fKdr(u[1]) ;
	dm_sKdr = (m_sKdr(u[1]) - u[5]) / t_m_sKdr(u[1]) ;
	
	return @SVector [dv, dm_Naf, dh_Naf, dm_fKdr, dm_sKdr]
end

function compartment_dynamics(comp::axon_is_t, u::Array{Float64}, V_connect::Array{Float64}, Ie::Float64)
	
	I_connect = 0.0 ;
	for i = 1 : length(comp.V_connect_idx)
		I_connect += 1.0/comp.Ra_connect_v[i] * (u[1] - V_connect[i]) ;
	end

	dv = (  - g_Naf_ax * u[2]^3 * u[3] * (u[1] - E_Na) * comp.area * comp.reduce_fact - 
				 g_TNC_ax * (u[1] - E_TNC) * comp.area * comp.reduce_fact -
				 g_fKdr_ax * u[4]^4 * (u[1] - E_K) * comp.area * comp.reduce_fact - 
				 g_sKdr_ax * u[5]^4 * (u[1] - E_K) * comp.area * comp.reduce_fact -
				 (u[1] - Em) / (RM_ax / (comp.reduce_fact * comp.area)) -
				 I_connect) / (CM * comp.reduce_fact * comp.area);

	dm_Naf = (m_Naf(u[1]) - u[2]) / t_m_Naf(u[1]) ;
	dh_Naf = (h_Naf(u[1]) - u[3]) / t_h_Naf(u[1]) ;
	dm_fKdr = (m_fKdr(u[1]) - u[4]) / t_m_fKdr(u[1]) ;
	dm_sKdr = (m_sKdr(u[1]) - u[5]) / t_m_sKdr(u[1]) ;
	
	return @SVector [dv, dm_Naf, dh_Naf, dm_fKdr, dm_sKdr]
end

function compartment_dynamics(comp::axon_in_t, u::Array{Float64}, V_connect::Array{Float64}, Ie::Float64)
	
	I_connect = 0.0 ;
	for i = 1 : length(comp.V_connect_idx)
		I_connect += 1.0/comp.Ra_connect_v[i] * (u[1] - V_connect[i]) ;
	end

	dv = (  (u[1] - Em) / (RM_my / (comp.reduce_fact * comp.area)) -
			I_connect) / (CM_my * comp.reduce_fact * comp.area);
	
	return @SVector [dv]
end

function compartment_dynamics(comp::prox_dend_t, u::Array{Float64}, V_connect::Array{Float64}, Ie::Float64)
	
	I_connect = 0.0 ;
	for i = 1 : length(comp.V_connect_idx)
		I_connect += 1.0/comp.Ra_connect_v[i] * (u[1] - V_connect[i]) ;
	end

	dv = (  - g_Naf_pd * u[2]^3 * u[3] * (u[1] - E_Na) * comp.area * comp.reduce_fact - 
				 g_TNC_pd * (u[1] - E_TNC) * comp.area * comp.reduce_fact -
				 g_h_pd * u[4]^2 * (u[1] - E_h) * comp.area * comp.reduce_fact - 
				 g_fKdr_pd * u[5]^4 * (u[1] - E_K) * comp.area * comp.reduce_fact - 
				 g_sKdr_pd * u[6]^4 * (u[1] - E_K) * comp.area * comp.reduce_fact -
				 p_CaHVA_pd * u[7]^3 * z_CaHVA^2 * F^2 * u[1] * (u[11] - Ca_out * exp(-z_CaHVA * F * u[1] / (R * Temp))) / 
						(R * Temp * (1.0 - exp(-z_CaHVA * F * u[1] / (R * Temp)))) * comp.area * comp.reduce_fact -
				 g_CaLVA_pd * u[8]^2 * u[9] * (u[1] - E_Ca) * comp.area * comp.reduce_fact - 
				 g_Sk_pd * u[10] * (u[1] - E_K) * comp.area * comp.reduce_fact -
				 (u[1] - Em) / (RM_d / (comp.reduce_fact * comp.area)) -
				 I_connect) / (CM * comp.reduce_fact * comp.area);

	dm_Naf = (m_Naf(u[1]) - u[2]) / t_m_Naf(u[1]) ;
	dh_Naf = (h_Naf(u[1]) - u[3]) / t_h_Naf(u[1]) ;
	dm_h = (m_h(u[1]) - u[4]) / t_m_h ;
	dm_fKdr = (m_fKdr(u[1]) - u[5]) / t_m_fKdr(u[1]) ;
	dm_sKdr = (m_sKdr(u[1]) - u[6]) / t_m_sKdr(u[1]) ;
	dm_CaHVA = (m_CaHVA(u[1]) - u[7]) / t_m_CaHVA(u[1]) ;
	dm_CaLVA = (m_CaLVA(u[1]) - u[8]) / t_m_CaLVA(u[1]) ;
	dh_CaLVA = (h_CaLVA(u[1]) - u[9]) / t_h_CaLVA(u[1]) ;
	dz_Sk = (z_Sk(u[11]) - u[10]) / t_z_Sk(u[11]) ;
	dCa = -(k_Ca_d / comp.Ca_vol) * p_CaHVA_s * u[7]^3 * z_CaHVA^2 * F^2 * u[1] * (u[11] - Ca_out * exp(-z_CaHVA * F * u[1] / (R * Temp))) / 
			(R * Temp * (1.0 - exp(-z_CaHVA * F * u[1] / (R * Temp)))) * comp.area * comp.reduce_fact - 
			(u[11] - Ca_base) / t_Ca ;

	return @SVector [dv, dm_Naf, dh_Naf, dm_h, dm_fKdr, dm_sKdr, dm_CaHVA, dm_CaLVA, dh_CaLVA, dz_Sk, dCa]

end

function compartment_dynamics(comp::dist_dend_t, u::Array{Float64}, V_connect::Array{Float64}, Ie::Float64)
	
	I_connect = 0.0 ;
	for i = 1 : length(comp.V_connect_idx)
		I_connect += 1.0/comp.Ra_connect_v[i] * (u[1] - V_connect[i]) ;
	end

	dv = (  - g_h_dd * u[2]^2 * (u[1] - E_h) * comp.area * comp.reduce_fact - 
				 p_CaHVA_dd * u[3]^3 * z_CaHVA^2 * F^2 * u[1] * (u[7] - Ca_out * exp(-z_CaHVA * F * u[1] / (R * Temp))) / 
						(R * Temp * (1.0 - exp(-z_CaHVA * F * u[1] / (R * Temp)))) * comp.area * comp.reduce_fact -
				 g_CaLVA_dd * u[4]^2 * u[5] * (u[1] - E_Ca) * comp.area * comp.reduce_fact - 
				 g_Sk_dd * u[6] * (u[1] - E_K) * comp.area * comp.reduce_fact -
				 (u[1] - Em) / (RM_d / (comp.reduce_fact * comp.area)) -
				 I_connect) / (CM * comp.reduce_fact * comp.area);

	dm_h = (m_h(u[1]) - u[2]) / t_m_h ;
	dm_CaHVA = (m_CaHVA(u[1]) - u[3]) / t_m_CaHVA(u[1]) ;
	dm_CaLVA = (m_CaLVA(u[1]) - u[4]) / t_m_CaLVA(u[1]) ;
	dh_CaLVA = (h_CaLVA(u[1]) - u[5]) / t_h_CaLVA(u[1]) ;
	dz_Sk = (z_Sk(u[7]) - u[6]) / t_z_Sk(u[7]) ;
	dCa = -(k_Ca_d / comp.Ca_vol) * p_CaHVA_s * u[3]^3 * z_CaHVA^2 * F^2 * u[1] * (u[7] - Ca_out * exp(-z_CaHVA * F * u[1] / (R * Temp))) / 
			(R * Temp * (1.0 - exp(-z_CaHVA * F * u[1] / (R * Temp)))) * comp.area * comp.reduce_fact - 
			(u[7] - Ca_base) / t_Ca ;

	return @SVector [dv, dm_h, dm_CaHVA, dm_CaLVA, dh_CaLVA, dz_Sk, dCa]

end

function compartment_currents(comp::soma_t, u::Array{Float64}, V_connect::Array{Float64})

	V_rep_v = repmat([u[1]], length(V_connect))

	current_v = [ g_Naf_s * u[2]^3 * u[3] * (u[1] - E_Na) * comp.area,  
					g_Nap_s * u[4]^3 * u[5] * (u[1] - E_Na) * comp.area,
					g_TNC_s * (u[1] - E_TNC) * comp.area,
			 		g_h_s * u[6]^2 * (u[1] - E_h) * comp.area,
					g_fKdr_s * u[7]^4 * (u[1] - E_K) * comp.area, 
					g_sKdr_s * u[8]^4 * (u[1] - E_K) * comp.area,
			 		p_CaHVA_s * u[9]^3 * z_CaHVA^2 * F^2 * u[1] * (u[13] - Ca_out * exp(-z_CaHVA * F * u[1] / (R * Temp))) / 
					(R * Temp * (1.0 - exp(-z_CaHVA * F * u[1] / (R * Temp)))) * comp.area,
					g_CaLVA_s * u[10]^2 * u[11] * (u[1] - E_Ca) * comp.area,
					g_Sk_s * u[12] * (u[1] - E_K) * comp.area,
			 		(u[1] - Em) / (RM_s / (comp.reduce_fact * comp.area)),
					1.0./(comp.Ra_connect_v') * (V_rep_v - V_connect)] ;

	return current_v
end

