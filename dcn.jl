
struct gate_t
	Vh::Float64		# mV
	k::Float64		# mV
end

struct tau_t
	A::Float64		# ms
	B::Float64		# mV
	C::Float64		# mV
	D::Float64		# mV
	E::Float64		# mV
	F::Float64		# ms
end

struct tau_h_Nap_t
	A::Float64		# ms
	D::Float64		# mV
	E::Float64		# mV
	F::Float64		# ms
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

	return 1.0e-3 / (31.746 * (exp((V - 5.0e-3) / -13.89e-3) + 1.0)^(-1) + 3.97e-1 * (V + 8.9e-3) * (exp((V + 8.9e-3) / 5.0e-3) - 1.0)^(-1))

end

function t_h_CaLVA(V::Float64)

	if V < -81.0e-3 	# V
		return 0.333e-3 * exp((V + 466.0e-3) / 66.0e-3)
	else
		return 0.333e-3 * exp((V + 21.0e-3) / -10.5e-3) + 9.32e-3
	end

end

function z_Sk(Ca::Float64)

	return (Ca^4) / (Ca^4 + (3.0e-7)^4) 

end

function t_z_Sk(Ca::Float64)

	if Ca < 0.005e-3	# M
		return 60.0e-3 - 11.2e3 * Ca 
	else
		return 4.0e-3
	end

end

function dcn_dyn(t, u, p, du)

	du[1] = (p - g_Naf * u[2]^3 * u[3] * (u[1] - E_Na) - 
					g_Nap * u[4]^3 * u[5] * (u[1] - E_Na) -
					g_TNC * (u[1] - E_TNC) -
					g_h * u[6]^2 * (u[1] - E_h) - 
					g_fKdr * u[7]^4 * (u[1] - E_K) - 
					g_sKdr * u[8]^4 * (u[1] - E_K) -
					p_CaHVA * u[9]^3 * z_CaHVA^2 * F^2 * u[1] * (u[13] - u0[13] * exp(-z_CaHVA * F * u[1] / (R * Temp))) / 
						(R * Temp * (1.0 - exp(-z_CaHVA * F * u[1] / (R * Temp)))) -
					g_CaLVA * u[10]^2 * u[11] * (u[1] - E_Ca) - 
					g_Sk * u[12] * (u[1] - E_K)) / Cm ;

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
	du[13] = B_Ca * p_CaHVA * u[9]^3 * z_CaHVA^2 * F^2 * u[1] * (u[13] - u0[13] * exp(-z_CaHVA * F * u[1] / (R * Temp))) / 
			(R * Temp * (1.0 - exp(-z_CaHVA * F * u[1] / (R * Temp)))) - (u[13] - Ca_base) / t_Ca ;

	println(u)
	println(du)
	#println(p_CaHVA * u[9]^3 * z_CaHVA^2 * F^2 * u[1] * (u[13] - Ca0 * exp(-z_CaHVA * F * u[1] / (R * Temp))))
	#println(exp(-z_CaHVA * F * u[1] / (R * Temp)))
	#println(-z_CaHVA * F * u[1])
	#println((R * Temp))
end

# Fast sodium current

g_Naf = 250.0 ;		# S/m^2
E_Na = 71.0e-3 ;	# V

m_Naf = gate_t( -45.0e-3, -7.3e-3 ) ;
t_m_Naf = tau_t( 5.83e-3, 6.4e-3, -9.0e-3, -97.0e-3, 17.0e-3, 0.025e-3 ) ;

h_Naf = gate_t( -42.0e-3, 5.9e-3 ) ;
t_h_Naf = tau_t( 16.67e-3, 8.3e-3, -29.0e-3, -66e-3, 9.0e-3, 0.2e-3 ) ;

# Persistent sodium current

g_Nap = 8.0 ;		# S/m^2

m_Nap = gate_t( -70.0e-3, -4.1e-3 ) ;
t_m_Nap = 50.0e-3 ;	# s

h_Nap = gate_t( -80.0e-3, 4.0e-3 ) ;
t_h_Nap = tau_h_Nap_t( 1750.0e-3, -65.0e-3, -8.0e-3, 250.0e-3 ) ;

# Tonic non-specific cation current

g_TNC = 0.3	;		# S/m^2
E_TNC = -35.0e-3 ;	# V

# Hyperpolarisation-activated cyclic nucleotide current

g_h = 2.0	;		# S/m^2
E_h = -45.0e-3 ;	# V

m_h = gate_t( -80.0e-3, 5.0e-3 ) ;
t_m_h = 400.0e-3 ;	# s

# Fast delayed rectifier fKdr

g_fKdr = 150.0 ; 	# S/m^2
E_K = -90.0e-3 ;	# V

m_fKdr = gate_t( -40.0e-3, -7.8e-3 ) ;
t_m_fKdr = tau_t( 13.9e-3, -40.0e-3, 12.0e-3, -40.0e-3, -13.0e-3, 0.1e-3 ) ;

# Slow delayed rectifier sKdr

g_sKdr = 125.0 ;	# S/m^2

m_sKdr = gate_t( -50.0e-3, -9.1e-3 ) ;
t_m_sKdr = tau_t( 14.95e-3, -50.0e-3, 21.74e-3, -50.0e-3, 13.91e-3, 0.05e-3 ) ;

# High voltage activated calcium current CaHVA

p_CaHVA = 7.5e-8 ;	# m/s

m_CaHVA = gate_t( -34.5e-3, -9.0e-3) ;

z_CaHVA = 2.0 ;			# valence of Ca 
R = 8.3145 ;			# J/(K*mol) , gas constant
F = 96480.0 ;			# C/mol , Faraday constant
Temp = 32.0 + 273.15 ;	# K

# Low voltage activated calcium current CaLVA / CaT

g_CaLVA = 1.5 ;		# S/m^2 
E_Ca = 139.0e-3 ;	# V

m_CaLVA = gate_t( -56.0e-3, -6.2e-3 ) ;
t_m_CaLVA = tau_t( 0.333e-3, -131.0e-3, -16.7e-3, -15.8e-3, 18.2e-3, 0.204e-3 ) ;

h_CaLVA = gate_t( -80.0e-3, 4.0e-3 ) ;

# Small conductance calcium dependent potassium current Sk

g_Sk = 2.2 ;		# S/m^2

# Calcium concentration

Ca_base = 50e-9 ;	# M
k_Ca = 3.45e-7 ;	# mol/C

radius_soma = 50.0e-6 ;	# m
thick_soma = 200.0e-9 ;	# m
B_Ca = k_Ca / (4.0 * pi * thick_soma * radius_soma^2) ;

t_Ca = 70.0e-3 ; 		# s 


#-------------------------
# Solving the ODE system 
#-------------------------

using DifferentialEquations 
using ODEInterface 	
using Plots

Cm = 0.0156 ;		# F/m^2
I =  100.0e-12 ;	# A 
pf = ParameterizedFunction(dcn_dyn, I) ;

V0 = -60.0e-3 ;		# V
Ca0 = 2.0e-3 ;		# M
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
	Ca0 ] ;

tspan = (0.0, 1.0) ;

prob = ODEProblem(pf, u0, tspan) ;

sol = solve(prob, alg=:radau, reltol=1e-8, abstol=1e-8);

plotlyjs(size = (700,500))

plt = plot(sol, vars = (0,1), legend = false)
title!("Soma voltage response with I = $I pA and V0 = $V0 mV")
xlabel!("t [s]")
ylabel!("V [V]")

#gui(plt)
#display(plt)

savefig(plt,"voltage.png")
 		
