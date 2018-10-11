module param
	use nrtype
	real(dp),parameter :: C_light=2.99792458E+10 
	real(dp),parameter :: hplanck_erg=6.6260695729E-27
	real(dp),parameter :: hplanck_ev=4.13566751691E-15
	real(dp),parameter :: erg_ev=6.2415E+11
	real(dp),parameter :: E_HI=13.6d0	
	real(dp),parameter :: E_lyalpha=10.198834d0		
	real(dp),parameter :: E_HeI=24.59d0		
	real(dp),parameter :: E_HeII=54.42d0	
	real(dp),parameter :: E_max=1000.d0!10000.0d0
	real(dp),parameter :: nu_HI=e_hi/hplanck_ev
	real(dp),parameter :: nu_HeI=e_hei/hplanck_ev
	real(dp),parameter :: nu_HeiI=e_heii/hplanck_ev
	real(dp),parameter :: nu_max=e_max/hplanck_ev
	real(dp),parameter :: nu_alpha=E_lyalpha/hplanck_ev
	real(dp),parameter :: nu_beta=c_light/1026.d-8
	real(dp),parameter :: nu_limit=c_light/912.d-8
	real(dp),parameter :: charge_e=1.60217656e-19
	real(dp),parameter :: Mass_e=9.10938188E-28			!mass of electron in gram
	real(dp),parameter :: Kboltz_ev=8.62911325E-5 
	real(dp),parameter :: Kboltz_erg=1.380648813e-16	
	real(dp),parameter :: solar_mass=1.9891d+33
	real(dp),parameter :: year_sec=3.15569d+7
	real(dp),parameter :: Megayear=3.1536E+13
	real(dp),parameter :: Stefan_Boltzmann=5.6704d-5 !erg cm^-2 s^-1 k^-4		
	real(dp),parameter :: Megaparsec=3.08568025E+24                 !megaparsec in cm
	real(dp),parameter :: kiloparsec=3.08568025E+21			!kiloparsec in cm	
	real(dp),parameter :: G = 6.67259e-8 !cm^3 g^-1 s^-2
	real(dp),parameter :: m_p = 1.6726231e-24 ! proton mass (g) 
	real(dp),parameter :: A10_HYPERFINE=2.85d-15
	real(dp),parameter :: Alpha_A=4.2e-13	!case A recombination
	real(dp),parameter :: sig_s=6.6524E-25
	real(dp),parameter :: Y_He =0.24d0



!!Planck..

	real(dp),parameter::omega_l=0.6911d0
	real(dp),parameter::omega_b=0.049d0
	real(dp),parameter::omega_r=0.0d-0
	real(dp),parameter::omega_m=0.3089d0
	real(dp),parameter::omega_k=1.d0-(omega_l+omega_m+omega_r)
	real(dp),parameter::sigma_8=0.81d0
	real(dp),parameter:: hlittle =0.678d0 ! little hubble h
	real(dp),parameter::ns = 0.968d0
	real(dp),parameter::k0=0.01027d0
	real(dp),parameter::dn_dlnk=-0.003d0
	real(dp),parameter::cmb_temp=2.726d0




	real(dp),parameter:: delta_c=1.686d0
	real(dp),parameter:: Ho  = (hlittle*3.2407d-18) ! s^-1 at z=0 
	real(dp),parameter:: RHOcrit_cgs =(3.0*Ho*Ho / (8.0*PI*G)) ! g pcm^-3 at z=0  9.2057E-030
	real(dp),parameter:: H_No = (RHOcrit_cgs*omega_b*(1.d0-Y_He)/m_p)  ! current hydrogen number density estimate  cm^3)  ~1.92e-7
	real(dp),parameter:: He_No = (RHOcrit_cgs*omega_b*Y_He/(4.0*m_p)) !  current helium number density estimate
	real(dp),parameter:: N_b0 = (H_No+He_No) ! present-day baryon num density, H + He 
	real(dp),parameter:: f_H = (H_No/(H_No+He_No))  ! hydrogen number fraction 
	real(dp),parameter:: f_He = (He_No/(H_No+He_No))  ! helium number fraction 
	real(dp),parameter:: T21 =0.0628d0 ! in K !temperature corresponding to the 21cm photon 
	real(dp),parameter:: rho_c=RHOcrit_cgs*(megaparsec**3.d0)/solar_mass!2.7755d+11! in M0/mpc3
	real(dp),parameter:: rho_c_0=rho_c*omega_m !!in M0/(Mpc)^3
	real(dp),parameter:: rho_c_h=rho_c_0/hlittle**2.d0 !!in M0/h/(Mpc/h)^3h)^3
	real(dp),parameter:: mui_molw=1.22d0!f_h+4.d0*f_he
	real(dp),parameter:: p_alpha=0.79d0	!fraction of X-ray energy goes into Lya photons
	real(dp),parameter:: AHe = H_No/N_b0	!!correction factor in presence of helium to the total uv photons 	
	

!! source properties for popII and popIII

	real(dp), parameter :: n_gamma_popII=4000.d0
	real(dp) :: popII_frac
	real(dp), parameter :: n_lya_lyb_popII=6520.d0
	real(dp), parameter :: n_lya_lyl_popII=9690.d0
	real(dp), parameter :: alpha_lya_lyb_popII=0.14d0



	real(dp), parameter ::	f_esc_fiducial = 0.1d0
	real(dp), parameter ::	f_star_fiducial = 0.05d0
	real(dp), parameter ::	f_X_fiducial = 1.d0
	real(dp), parameter ::	alpha_X_fiducial = 0.5d0
	real(dp), parameter ::  t_vir_fiducial = 1d4
	real(dp):: n_rec_comm, Tk_comm, redshift_comm, xi_comm, xe_comm


!!additional sim parameters
	real(dp),parameter:: C_clumping=2.d0
	real(dp), parameter :: ion_cut=0.8d0	!!below which the analytical model deviates..
	real(dp),parameter:: z_star=35.d0 
	real(dp),parameter:: z_end = 6.d0

	 character(*), parameter :: output_path    = '../outputs/'



!!other common declaration


	real(dp) :: f_esc_param, N_ion_param, f_star_param,  Alpha_s_param, L_0_param
	real(dp) :: f_X_param,  t_vir_param!, vcmin
	REAL(DP)::Jalpha_star,Jalpha_x, delz_fcoll, Q_global
	integer, parameter :: ndim_fcoll = 100
	real(dp), dimension(2, ndim_fcoll)::arr_fcoll, arr_corel_dd, av_bias_arr
	character(12)::f_star_s, z_s, f_esc_s, alpha_star_pop_s, nlya_pop_s, N_ion_s
	character(12):: Alpha_s_s, f_X_s,  t_vir_s
	character(180) :: filename_comm





	integer, parameter :: sheth_tormen = 0	!!1 will set ST, 0 will set PS

	integer, parameter :: nkbin=10
	real(dp), dimension(nkbin)::kbin_comm, dndmarr_cdm, hmassarr, dimless_PS
	integer :: k01id

	real(dp), dimension(nkbin) :: Wxk_comm, Walphak_comm, gT_comm


!! number of parameters for parameter space study



	real(dp), parameter :: z_cosmic_dawn = 30.d0
	real(dp), parameter :: z_eor_end = 6.d0
	real(dp), parameter :: delz_eor = 0.2d0
	integer, parameter :: nzloop = int((z_cosmic_dawn-z_eor_end)/delz_eor)



!!make faster xhi

	real(dp) :: sigmamin_comm, k_zeta_comm, mminz_comm, fcollz_comm, B0z_comm, B1z_comm, ion_ef_comm, test_comm, mass_av_biasz
	integer, parameter :: nrdim_comm = 100
	integer, parameter :: massdim_comm = 100
	real(dp), parameter ::mmin_bin = 1d6
	real(dp), parameter ::mmax_bin = 5d16

	real(dp), parameter ::rmin_bin = 1d-2
	real(dp), parameter ::rmax_bin = 1d4
	real(dp), dimension(massdim_comm) :: M_arr_comm, dx_arr_comm, V_arr_comm, sigma_arr_comm
	real(dp), dimension(massdim_comm) :: dndlnm_arr_comm, halo_bias_arr_comm, R_m_arr_comm
	real(dp), dimension(nrdim_comm) :: r_arr_comm, cor_dd_arr_comm, cor_dx_arr_comm, cor_xx_arr_comm
	real(dp), dimension(massdim_comm) :: BSD_arr, BSD_R_arr



	integer, parameter :: nparam = 4	!!number of possible parameters , fx, al, n_gamma
	logical, parameter :: write_output_file = .false.
	logical, parameter :: wr_messages = .false.




end module param
