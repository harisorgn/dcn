
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

	return (Ca^4) / (Ca^4 + (3.0e-4)^4) 

end

function t_z_Sk(Ca::Float64)

	if Ca < 0.005	# mM
		#return 60.0e-3 - 11.2 * Ca
		return 0.001 - 0.1867 * Ca  
	else
		#return 4.0e-3
		return 6.667e-5
	end

end


# Passive properties

const CM = 0.0157 ;			# F/m^2
const RM_s = 3.556 ;		# Ω m^2
const RM_d = 3.556 ;		# Ω m^2
const RM_ax = 3.556 ;		# Ω m^2
const RA = 2.353 ;			# Ω m^2

const Em = -66e-3 ;			# V

# Fast sodium current

const g_Naf_s = 250.0 ;		# S/m^2
const g_Naf_pd = 100.0 ;	# S/m^2
const g_Naf_ax = 500.0 ;	# S/m^2

const E_Na = 71.0e-3 ;	# V

m_Naf = gate_t( -45.0e-3, -7.3e-3 ) ;
t_m_Naf = tau_t( 5.83e-3, 6.4e-3, -9.0e-3, -97.0e-3, 17.0e-3, 0.025e-3 ) ;

h_Naf = gate_t( -42.0e-3, 5.9e-3 ) ;
t_h_Naf = tau_t( 16.67e-3, 8.3e-3, -29.0e-3, -66e-3, 9.0e-3, 0.2e-3 ) ;

# Persistent sodium current

#const g_Nap_s = 8.0 ;		# S/m^2
const g_Nap_s = 0.0 ;

m_Nap = gate_t( -70.0e-3, -4.1e-3 ) ;
t_m_Nap = 50.0e-3 ;	# s

h_Nap = gate_t( -80.0e-3, 4.0e-3 ) ;
t_h_Nap = tau_h_Nap_t( 1750.0e-3, -65.0e-3, -8.0e-3, 250.0e-3 ) ;

# Tonic non-specific cation current

const g_TNC_s = 0.3	;		# S/m^2
const g_TNC_pd = 0.06	;	# S/m^2
const g_TNC_ax = 0.35	;	# S/m^2

const E_TNC = -35.0e-3 ;	# V

# Hyperpolarisation-activated cyclic nucleotide current

#const g_h_s = 2.0	;		# S/m^2
const g_h_s = 0.0 ;
const g_h_pd = 4.0	;		# S/m^2
const g_h_dd = 6.0	;		# S/m^2

const E_h = -45.0e-3 ;	# V

m_h = gate_t( -80.0e-3, 5.0e-3 ) ;
t_m_h = 400.0e-3 ;	# s

# Fast delayed rectifier fKdr

const g_fKdr_s = 150.0 ; 	# S/m^2
const g_fKdr_pd = 90.0 ; 	# S/m^2
const g_fKdr_ax = 300.0 ; 	# S/m^2

const E_K = -90.0e-3 ;	# V

m_fKdr = gate_t( -40.0e-3, -7.8e-3 ) ;
t_m_fKdr = tau_t( 13.9e-3, -40.0e-3, 12.0e-3, -40.0e-3, -13.0e-3, 0.1e-3 ) ;

# Slow delayed rectifier sKdr

const g_sKdr_s = 125.0 ;	# S/m^2
const g_sKdr_pd = 75.0 ;	# S/m^2
const g_sKdr_ax = 250.0 ;	# S/m^2

m_sKdr = gate_t( -50.0e-3, -9.1e-3 ) ;
t_m_sKdr = tau_t( 14.95e-3, -50.0e-3, 21.74e-3, -50.0e-3, -13.91e-3, 0.05e-3 ) ;

# High voltage activated calcium current CaHVA

const p_CaHVA_s = 7.5e-8 ;	# m/s
const p_CaHVA_pd = 5.0e-8 ;	# m/s
const p_CaHVA_dd = 5.0e-8 ;	# m/s

m_CaHVA = gate_t( -34.5e-3, -9.0e-3) ;

const z_CaHVA = 2.0 ;	# valence of Ca 
const R = 8.3145 ;			# J/(K*mol) , gas constant
const F = 96480.0 ;			# C/mol , Faraday constant
const Temp = 32.0 + 273.15 ;# K

# Low voltage activated calcium current CaLVA / CaT

#const g_CaLVA_s = 1.5 ;		# S/m^2
const g_CaLVA_s = 0.0 ;
const g_CaLVA_pd = 3.0 ;		# S/m^2
const g_CaLVA_dd = 3.0 ;		# S/m^2

const E_Ca = 139.0e-3 ;	# V

m_CaLVA = gate_t( -56.0e-3, -6.2e-3 ) ;
t_m_CaLVA = tau_t( 0.333e-3, -131.0e-3, -16.7e-3, -15.8e-3, 18.2e-3, 0.204e-3 ) ;

h_CaLVA = gate_t( -80.0e-3, 4.0e-3 ) ;

# Small conductance calcium dependent potassium current Sk

const g_Sk_s = 2.2 ;		# S/m^2
const g_Sk_pd = 0.66 ;		# S/m^2
const g_Sk_dd = 0.66 ;		# S/m^2

# Calcium concentration

const Ca_in = 50.0e-6 ;		# mM
const Ca_out = 2.0 ;		# mM
const Ca_base = 50.0e-6 ;	# mM

const k_Ca_s = 3.45e-7 ;	# mol/C
const k_Ca_d = 1.0364e-6 ;	# mol/C
const t_Ca = 70.0e-3 ; 		# s 

const Ca_thick = 200.0e-9 ;

const n_eq_s = 13 ;
const n_eq_ax = 5 ;
const n_eq_pd = 11 ;
const n_eq_dd = 7 ;

const V0 = -70e-3 ;
const Ca0 = 50.0e-6 ;

const Ie = 50.0e-12 ;
const t_pulse_on = 0.5 ;
const t_pulse_dur = 1.5 ;
