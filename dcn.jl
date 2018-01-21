
struct gate_t
	Vh::Float64		# mV
	k::Float64		# mV
end

struct tau_coeff_t
	A::Float64		# ms
	B::Float64		# mV
	C::Float64		# mV
	D::Float64		# mV
	E::Float64		# mV
	F::Float64		# ms
end

struct tau_nap_coeff_t
	A::Float64		# ms
	D::Float64		# mV
	E::Float64		# mV
	F::Float64		# ms
end

function (x::gate_t)(V::Float64)

	return 1 / (1 + exp((V - x.Vh)/x.k))

end

function (t::tau_coeff_t)(V::Float64)

	return t.A / (exp((V - t.B) / t.C) + exp((V - t.D) / t.E)) + t.F

end

function (t::tau_nap_coeff_t)(V::Float64)

	return t.A / (1 + exp((V - t.D) / t.E)) + t.F

end

function dcn_dyn(t, u, p, du)

	du[1] = (p - g_Naf * u[2]^3 * u[3] * (u[1] - E_Na) - 
					g_Nap * u[4]^3 * u[5] * (u[1] - E_Na) -
					g_TNC * (u[1] - E_TNC) -
					g_h * u[6]^2 * (u[1] - E_h) - 
					g_fKdr * u[7]^4 * (u[1] - E_K) - 
					g_sKdr * u[8]^4 * (u[1] - E_K) ) / Cm ;
	du[2] = (m_Naf(u[1]) - u[2]) / t_m_Naf(u[1]) ;
	du[3] = (h_Naf(u[1]) - u[3]) / t_h_Naf(u[1]) ;
	du[4] = (m_Nap(u[1]) - u[4]) / t_m_Nap ;
	du[5] = (h_Nap(u[1]) - u[5]) / t_h_Nap(u[1]) ;
	du[6] = (m_h(u[1]) - u[6]) / t_m_h ;
	du[7] = (m_fKdr(u[1]) - u[7]) / t_m_fKdr(u[1]) ;
	du[8] = (m_sKdr(u[1]) - u[8]) / t_m_sKdr(u[1]) ;

end

# Fast sodium current

g_Naf = 250.0 ;		# S/m^2
E_Na = 71.0 ;		# mV

m_Naf = gate_t( -45.0, -7.3 ) ;
t_m_Naf = tau_coeff_t( 5.83, 6.4, -9.0, -97.0, 17.0, 0.025 ) ;

h_Naf = gate_t( -42.0, 5.9 ) ;
t_h_Naf = tau_coeff_t( 16.67, 8.3, -29.0, -66, 9.0, 0.2 ) ;

# Persistent sodium current

g_Nap = 8.0 ;		# S/m^2

m_Nap = gate_t( -70.0, -4.1 ) ;
t_m_Nap = 50.0 ;	# ms

h_Nap = gate_t( -80.0, 4.0 ) ;
t_h_Nap = tau_nap_coeff_t( 1750.0, -65.0, -8.0, 250.0 ) ;

# Tonic non-specific cation current

g_TNC = 0.3	;		# S/m^2
E_TNC = -35.0 ;		# mV

# Hyperpolarisation-activated cyclic nucleotide current

g_h = 2.0	;		# S/m^2
E_h = -45.0 ;		# mV

m_h = gate_t( -80.0, 5 ) ;
t_m_h = 400.0 ;		# ms

# Fast delayed rectifier fKdr

g_fKdr = 150.0 ; 
E_K = -90.0 ;			# mV

m_fKdr = gate_t( -40.0, -7.8 ) ;
t_m_fKdr = tau_coeff_t( 13.9, -40.0, 12.0, -40.0, -13.0, 0.1 ) ;

# Slow delayed rectifier sKdr

g_sKdr = 125 ;		# S/m^2

m_sKdr = gate_t( -50.0, -9.1 ) ;
t_m_sKdr = tau_coeff_t( 14.95, -50.0, 21.74, -50.0, 13.91, 0.05 ) ;


using DifferentialEquations 
using ODEInterface 	
using Plots

Cm = 0.0156 ;		# F/m^2
I = 0.0 ;			# pA
pf = ParameterizedFunction(dcn_dyn, I) ;

V0 = -50.0 ;		# mV
u0 = [V0 ; m_Naf(V0) ; h_Naf(V0) ; m_Nap(V0) ; h_Nap(V0) ; m_h(V0) ; m_fKdr(V0) ; m_sKdr(V0)] ;
tspan = (0.0, 1.0) ;

prob = ODEProblem(pf, u0, tspan) ;

sol = solve(prob, alg=:radau );

plotlyjs()
plt = plot(sol, vars = (0,1), show=true) 

#gui(plt)

#gui()
#display(plt)
savefig(plt,"voltage.png")
 		
