
# Passive properties

const CM = 0.0157 ;			# F/m^2
const CM_my = 0.01 * CM ;	# F/m^2
const RM_s = 3.556 ;		# Ω m^2
const RM_d = 3.556 ;		# Ω m^2
const RM_ax = 3.556 ;		# Ω m^2
const RM_my = 10.0 ;		# Ω m^2
const RA = 2.353 ;			# Ω/m

const Em = -66e-3 ;			# V

# Fast sodium current

const g_Naf_s = 250.0 ;		# S/m^2
const g_Naf_pd = 100.0 ;	# S/m^2
const g_Naf_ax = 500.0 ;	# S/m^2

const E_Na = 71.0e-3 ;	# V


# Persistent sodium current

const g_Nap_s = 0.0 ;

const t_m_Nap = 50.0e-3 ;	# s

# Tonic non-specific cation current

const g_TNC_s = 0.3	;		# S/m^2
const g_TNC_pd = 0.06	;	# S/m^2
const g_TNC_ax = 0.35	;	# S/m^2

const E_TNC = -35.0e-3 ;	# V

# Hyperpolarisation-activated cyclic nucleotide current

const g_h_s = 1.0	;		# S/m^2
#const g_h_s = 0.0 ;
const g_h_pd = 4.0	;		# S/m^2
const g_h_dd = 6.0	;		# S/m^2

const E_h = -45.0e-3 ;	# V

#m_h = gate_t( -80.0e-3, 5.0e-3 ) ;
const t_m_h = 400.0e-3 ;	# s

# Fast delayed rectifier fKdr

const g_fKdr_s = 150.0 ; 	# S/m^2
const g_fKdr_pd = 90.0 ; 	# S/m^2
const g_fKdr_ax = 300.0 ; 	# S/m^2

const E_K = -90.0e-3 ;	# V

# Slow delayed rectifier sKdr

const g_sKdr_s = 125.0 ;	# S/m^2
const g_sKdr_pd = 75.0 ;	# S/m^2
const g_sKdr_ax = 250.0 ;	# S/m^2

# High voltage activated calcium current CaHVA

const p_CaHVA_s = 7.5e-8 ;	# m/s
const p_CaHVA_pd = 5.0e-8 ;	# m/s
const p_CaHVA_dd = 5.0e-8 ;	# m/s

const z_CaHVA = 2.0 ;	# valence of Ca 
const R = 8.3145 ;			# J/(K*mol) , gas constant
const F = 96480.0 ;			# C/mol , Faraday constant
const Temp = 32.0 + 273.15 ;# K

# Low voltage activated calcium current CaLVA / CaT

const g_CaLVA_s = 0.0 ;
const g_CaLVA_pd = 3.0 ;		# S/m^2
const g_CaLVA_dd = 3.0 ;		# S/m^2

const E_Ca = 139.0e-3 ;	# V

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
const n_eq_axin = 1 ;
const n_eq_pd = 11 ;
const n_eq_dd = 7 ;

const V0 = -70e-3 ;
const Ca0 = 50.0e-6 ;

const Ie = 00.0e-12 ;
const t_pulse_on = 0.5 ;
const t_pulse_dur = 1.5 ;

