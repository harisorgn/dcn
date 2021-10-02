struct gate_t
	Vh::Float64		# V
	k::Float64		# V
end

struct tau_t
	A::Float64		# s
	B::Float64		# V
	C::Float64		# V
	D::Float64		# V
	E::Float64		# V
	F::Float64		# s
end

struct tau_h_Nap_t
	A::Float64		# s
	D::Float64		# V
	E::Float64		# V
	F::Float64		# s
end

function (x::gate_t)(V::Float64)

	return 1.0 / (1.0 + exp((V - x.Vh)/x.k))

end

function (t::tau_t)(V::Float64)

	return t.A / (exp((V - t.B) / t.C) + exp((V - t.D) / t.E)) + t.F

end

function (t::tau_h_Nap_t)(V::Float64)

	return t.A / (1.0 + exp((V - t.D) / t.E)) + t.F

end

function t_m_CaHVA(V::Float64)

	return 1.0e-3 / (31.746 * (exp((V - 5.0e-3) / -13.89e-3) + 1.0)^(-1.0) + 3.97e-1 * (V + 8.9e-3) * (exp((V + 8.9e-3) / 5.0e-3) - 1.0)^(-1.0))

end

function t_h_CaLVA(V::Float64)

	if V < -81.0e-3 	# V
		return 0.333e-3 * exp((V + 466.0e-3) / 66.0e-3)
	else
		return 0.333e-3 * exp((V + 21.0e-3) / -10.5e-3) + 9.32e-3
	end

end

function z_Sk(Ca::Float64)

	return (Ca^4.0) / (Ca^4.0 + (3.0e-4)^4.0) 

end

function t_z_Sk(Ca::Float64)

	if Ca < 0.005	# mM
		return 0.001 - 0.1867 * Ca  
	else
		return 6.667e-5
	end

end

function f(Vh::Float64, k::Float64)

	function g(V::Float64)

		return 1.0 / (1.0 + exp((V - Vh)/k))
	end
	return g
end

m_Naf(V::Float64) = 1.0 / (1.0 + exp((V + 45.0e-3)/-7.3e-3)) ;
h_Naf(V::Float64) = 1.0 / (1.0 + exp((V + 42.0e-3)/5.9e-3)) ;
t_m_Naf(V::Float64) = 5.83e-3 / (exp((V - 6.4e-3) / -9.0e-3) + exp((V + 97.0e-3) / 17.0e-3)) + 0.025e-3 ;
t_h_Naf(V::Float64) = 16.67e-3 / (exp((V - 8.3e-3) / -29.0e-3) + exp((V + 66.0e-3) / 9.0e-3)) + 0.2e-3 ;

m_Nap(V::Float64) = 1.0 / (1.0 + exp((V + 70.0e-3)/-4.1e-3)) ;
h_Nap(V::Float64) = 1.0 / (1.0 + exp((V + 80.0e-3)/4.0e-3)) ;
t_h_Nap(V::Float64) = 1750.0e-3 / (1.0 + exp((V + 65.0e-3) / -8.0e-3)) + 250.0e-3 ;

m_h(V::Float64) = 1.0 / (1.0 + exp((V + 80.0e-3)/5.0e-3)) ;

m_fKdr(V::Float64) = 1.0 / (1.0 + exp((V + 40.0e-3)/-7.8e-3)) ;
t_m_fKdr(V::Float64) = 13.9e-3 / (exp((V + 40.0e-3) / 12.0e-3) + exp((V + 40.0e-3) / -13.0e-3)) + 0.1e-3 ;

m_sKdr(V::Float64) = 1.0 / (1.0 + exp((V + 50.0e-3)/-9.1e-3)) ;
t_m_sKdr(V::Float64) = 14.95e-3 / (exp((V + 50.0e-3) / 21.74e-3) + exp((V + 50.0e-3) / -13.91e-3)) + 0.05e-3 ;

m_CaHVA(V::Float64) = 1.0 / (1.0 + exp((V + 34.5e-3)/-9.0e-3)) ;

m_CaLVA(V::Float64) = 1.0 / (1.0 + exp((V + 56.0e-3)/-6.2e-3)) ;
t_m_CaLVA(V::Float64) = 0.333e-3 / (exp((V + 131.0e-3) / -16.7e-3) + exp((V + 15.8e-3) / 18.2e-3)) + 0.204e-3 ;

h_CaLVA(V::Float64) = 1.0 / (1.0 + exp((V + 80.0e-3)/4.0e-3)) ;