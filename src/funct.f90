module funct      
	use param
	use adaptint
	use cosmo
	use press_sch

	real(dp)::z_comm,r_comm,d2fact
	real(dp)::z_Wa_comm,nu_Wa_comm,k_Wa_comm
	real(dp)::z_Wx_comm,E_Wx_comm,k_Wx_comm
contains


!***************************TEMERATURE FUNCTIONS***************************

!Compton heating term/baryon
FUNCTION dcompton_heating_dt(z,tk,xe) 
	implicit none
	REAL(DP), INTENT(IN) :: z,tk,xe
	REAL(DP)::dcompton_heating_dt

	dcompton_heating_dt=(xe/(1.d0+f_he+xe))*(tcmb(z)-Tk)*(8.d0*sig_s/3.d0/mass_e/c_light)*u_gamma(z)

END FUNCTION dcompton_heating_dt


!energy deposited per unit vol as heat by X rays
SUBROUTINE Xray_heating(z,xe, Eps_xz, nbin, kbin, Wxk)			!one may apply log scale
	implicit none					!adaptive method takes long time
	REAL(DP), INTENT(IN) :: z,xe
	REAL(DP), INTENT(out) :: Eps_xz
	INTEGER, INTENT(IN) :: nbin
	REAL(DP), DIMENSION(nbin), INTENT(in) :: kbin
	REAL(DP), DIMENSION(nbin), INTENT(out) :: Wxk

	REAL(DP), DIMENSION(nbin) :: Wxk_tmp
	REAL(DP)::nu1,nu2,delnu,res,  WX, reswx, jxnuz, nu, q_hi, q_hei, q_heii, nhi, nhei, nheii, gamma_av
	real(DP), parameter :: increase_f = 1.1d0

	delnu=100.d0/hplanck_ev
	nu1=100.d0/hplanck_ev
	Eps_xz=0.d0
	gamma_av=0.d0
	Wxk=0.d0

	do !nu integration
		nu2=nu1+delnu
		nu=(nu1+nu2)*0.5d0
		if(nu .gt. nu_max) exit

		call sub_JX_nuz(nu, z, xe, Jxnuz, nbin, kbin,  Wxk_tmp)
		nhi=n_hi(z,xe)
		nhei=n_hei(z,xe)
		nheii=n_heii(z,xe)
		q_hi=(nu-nu_hi)*sigma_hi(nu)*delnu
		q_hei=(nu-nu_hei)*sigma_hei(nu)*delnu
		q_heii=(nu-nu_heii)*sigma_heii(nu)*delnu

		res=Jxnuz*(q_hi*nhi + q_hei*nhei + q_heii*nheii)

		Eps_xz=Eps_xz + res
		gamma_av = gamma_av + Jxnuz*(q_hi + q_hei + q_heii)

		Wxk=Wxk + Wxk_tmp*(q_hi + q_hei + q_heii)
		!if(res<1d-2*Eps_xz) exit	!!newly added
		nu1=nu2
		delnu=delnu*increase_f
	end do
	Eps_xz=4.d0*pi*hplanck_erg*f_heat(xe)*Eps_xz
	gamma_av=4.d0*pi*gamma_av
	Wxk=Wxk*4.d0*pi/gamma_av



END SUBROUTINE Xray_heating

subroutine sub_JX_nuz(nu, z, xe, Jxnuz, nbin, kbin,  Wxk)
	implicit none					!adaptive method takes long time
	REAL(DP), INTENT(IN) :: nu, z,xe
	REAL(DP), INTENT(out) ::Jxnuz
	INTEGER, INTENT(IN) :: nbin
	REAL(DP), DIMENSION(nbin), INTENT(in) :: kbin
	REAL(DP), DIMENSION(nbin), INTENT(out) :: Wxk

	REAL(DP) :: k, r
	REAL(DP)::J_X,nu_p,z1,z2,delz,res, res1, zp, q1, q2, gfc, norm
	real(dp), dimension(nkbin) :: J0arr, J2arr
	real(DP), parameter :: increase_f = 1.1d0
	INTEGER :: i

	delz=0.05d0
	z1=z
	Jxnuz=0.d0
	Wxk=0.d0
	gfc = D(z)

	DO !redshift integration
		z2=z1+delz
		zp=(z1+z2)*0.5d0
		if(zp .gt. z_star) exit
		nu_p=nu*(1.d0+zp)/(1.d0+z)
		res=(1.d0/hubble_constant(zp))*delz*epsilon_cap_X(zp,nu_p)*exp(-(tau(nu,z,zp,xe)))
		Jxnuz=Jxnuz+res

		r=comoving_distance(zp,z)	!!Mpc

		r=r*hlittle		!Mpc/h
		DO i=1, nbin
			k=kbin(i)	!h/Mpc
			J0arr(i) = J0(k*r)
			J2arr(i)=J2(k*r)
		END DO

		q1=D(zp)/gfc
		q2=1.d0+Bias(zp)
		Wxk=Wxk + q1*res*(q2*j0arr - 2.d0/3.d0*J2arr)

		if(res<1d-2*Jxnuz) exit	!!newly added
		z1=z2
		delz=delz*increase_f
	end do

	norm=((1.d0+z)**2.d0)*c_light/4.d0/pi
	Jxnuz=Jxnuz*norm
	Wxk=Wxk*norm


END subroutine sub_JX_nuz



!Optical depth, z'>z
function tau(nu,z,z_prime,xe) 
	implicit none
	real(dp),intent(in)::nu,z,z_prime,xe
	real(dp)::tau,z1,z2,delz,res,nu_pp

	delz=(z_prime-z)/100.d0!0.05d0
	z1=z
	res=0.d0

	do
		z2=z1+delz
		if(z2 .gt. z_prime) exit
		nu_pp=nu*(1.d0+z2)/(1.d0+z)
		res=res+(sigma_HI(nu_PP)*n_HI(z2,xe)+sigma_HeI(nu_PP)*n_HeI(z2,xe)+sigma_HeiI(nu_PP)*n_HeiI(z2,xe))*delz*(-dt_dz(z2))
		z1=z2
	end do
	tau=res*c_light

end function tau


function epsilon_cap_X(z,nu)
	implicit none
	REAL(DP), INTENT(IN) :: z,nu
	real(dp)::epsilon_cap_X

	epsilon_cap_X=epsilon_cap_X_nu(nu)*sfrd(z)
	
end function epsilon_cap_X


!! this estimate number of X-ray photons per unit frequency
!! 
function epsilon_cap_X_nu(nu)
	implicit none
	REAL(DP), INTENT(IN) :: nu
	real(dp)::epsilon_cap_X_nu,nu_0,nx,norm, al, nu30kev, nu100ev

	L_0_param=f_X_param*1.d+40*year_sec/solar_mass


	nu_0=1000.d0/hplanck_ev

	nu30kev=30000.d0/hplanck_ev
	nu100ev=100.d0/hplanck_ev 

	norm=(alpha_s_param)/nu_0	
!!This is an approximation, see also https://arxiv.org/pdf/1310.0465.pdf !!assuming lng span of x-ray band

!	norm=1.d0/(nu30kev-nu100ev) !! frequency difference
	epsilon_cap_X_nu=(L_0_param/hplanck_erg/nu_0)*(nu/nu_0)**(-alpha_s_param-1.d0)*norm		!per unit frequency

end function epsilon_cap_X_nu

!star formation rate density
function sfrd(z)
	implicit none
	REAL(DP), INTENT(IN) :: z
	real(dp)::sfrd

	sfrd=RHOcrit_cgs*omega_b*f_star_param*dfcoll_dt(z)

end function sfrd


!rate of production of ionising photon by UV rays per baryon
function gamma_i(z)
	implicit none
	real(dp),intent(in)::z
	real(dp)::gamma_i
	gamma_i=xi(z)*dfcoll_dt(z)
end function gamma_i



function xi(z)
	implicit none
	real(dp),intent(in)::z
	real(dp)::xi
	xi=Ahe*f_star_param*f_esc_param*N_ion_param
end function xi



function dfcoll_dt(z)
	implicit none
	real(dp),intent(in)::z
	real(dp)::dfcoll_dt
	real(dp)::z1,z2,fc1,fc2,delz

	delz=delz_fcoll+0.05d0
	z1=z+delz
	z2=z-delz
	fc1=fcoll(z1)
	fc2=fcoll(z2)
	dfcoll_dt=(fc1-fc2)/2.d0/delz/dt_dz(z)

end function dfcoll_dt


function dfcoll_dz(z)
	implicit none
	real(dp),intent(in)::z
	real(dp)::dfcoll_dz
	real(dp)::z1,z2,fc1,fc2,delz

	delz=delz_fcoll+0.05d0
	z1=z+delz
	z2=z-delz
	fc1=fcoll(z1)
	fc2=fcoll(z2)
	dfcoll_dz=(fc1-fc2)/2.d0/delz

end function dfcoll_dz

!mass fraction
function fcoll(z)
implicit none
real(dp),intent(in)::z
real(dp)::fcoll,mass1,r1,sig,nu, zp, zmid
integer::i, imin
real(dp), dimension(ndim_fcoll) :: x, y

x=arr_fcoll(1,:)
y=arr_fcoll(2,:)
fcoll=linear_interpol(z,x,y)


end function fcoll



!collapsed fraction from PS..
!!
function fcoll_analytical(z)
	implicit none
	real(dp),intent(in)::z
	real(dp)::fcoll_analytical,mass1,r1,sig,nu

	mass1=mass_min_cdm(t_vir_param,z) !in M0/h
	r1=MtoRadius(mass1)	!!Mpc/h
	sig=sigma_cdm(r1)
	nu=delta_c/(d(z)*sig)/sqrt(2.d0) 
	fcoll_analytical=erfc(nu)

end function fcoll_analytical





!********************************TEMPERATURE ENDS*****************************


!********************************SPIN TEMPERATURE START***************************

!!calculate spin temperature..
!checked
function get_Ts(Trad,Tk,Tceff,xa_tilde,xc)
	implicit none
	real(dp),intent(in)::Trad,Tk,Tceff,xa_tilde,xc
	real(dp)::GET_TS,TSinv

 	TSinv = ((1.0d0/Trad) + (xa_tilde/Tceff) + (xc/TK))/(1.0d0+xa_tilde+xc)
 	get_ts = 1.0d0/TSinv
 
end function get_ts

!collisional coupling coefficient
!https://arxiv.org/abs/astro-ph/0607234 equation 3
!checked
function xcoll(z,TK,delta,xe)
	implicit none
	real(dp),intent(in)::z,TK,delta,xe
	real(dp)::xcoll
	real(dp):: k_hh,k_eh,nH,Trad,ne

	Trad = tcmb(z)
	nH = H_No*((1.0+z)**3.)*(1.0+delta)
	k_hh = kappa_10(TK,0)
	ne = xe*N_b0*((1.0+z)**3.)*(1.0+delta)
	k_eh = kappa_10_elec(TK,0)
	xcoll= 4.d0*T21/Trad/A10_HYPERFINE/3.d0*(nh*k_hh + ne*k_eh)

end function xcoll

!collisional coupling .. only electron collision
!checked
function xcoll_e(z,TK,delta,xe)
	implicit none
	real(dp),intent(in)::z,TK,delta,xe
	real(dp)::xcoll_e
	real(dp):: k_hh,k_eh,nH,Trad,ne

	Trad = tcmb(z)
	ne = xe*N_b0*((1.0+z)**3)*(1.0+delta)
	k_eh = kappa_10_elec(TK,0)
	xcoll_e= 4.d0*T21/Trad/A10_HYPERFINE/3.d0*( ne*k_eh)
end function xcoll_e

!!collisional coupling .. only h atom collision
!checked
function xcoll_h(z,TK,delta,xe)
	implicit none
	real(dp),intent(in)::z,TK,delta,xe
	real(dp)::xcoll_h
	real(dp):: k_hh,k_eh,nH,Trad,ne

	Trad = tcmb(z)
	nH = H_No*((1.0+z)**3)*(1.0+delta)
	k_hh = kappa_10(TK,0)

	xcoll_h= 4.d0*T21/Trad/A10_HYPERFINE/3.d0*(nh*k_hh)
end function xcoll_h



!from mesinger et al 2010 equ 7
!see also https://arxiv.org/pdf/astro-ph/0508381.pdf
!checked
function xalpha_tilde(z,Jalpha,TK,TS,delta,xe)
	implicit none
	real(dp),intent(in)::z,Jalpha,TK,TS,delta,xe
	real(dp):: tgp,Stilde,xalpha_tilde,Trad
	Trad = tcmb(z)
  	tgp = taugp(z,delta,xe)
  	Stilde = Salpha_tilde(TK,TS,tgp)
  	xalpha_tilde = 1.72d11/(1.0+z)*Stilde*Jalpha 
end function xalpha_tilde


!Hirata 2006 equ 35
!http://adsabs.harvard.edu/abs/2006MNRAS.367..259H
!checked
function taugp(z, delta,  xe)
	implicit none
	real(dp),intent(in)::z,delta,xe
	real(dp)::taugp
	taugp= 1.3485342d-7 / hubble_constant(z)*H_No*(1.0+z)**3. * (1.0+delta)*(1.0-xe)
end function taugp

!Hirata 2006 equ 40 and 41
!http://adsabs.harvard.edu/abs/2006MNRAS.367..259H
!checked
function Salpha_tilde(TK,TS,tauGP)
	implicit none
	real(dp),intent(in)::tk,ts,taugp
	real(dp)::xi,ans,Salpha_tilde

  	xi = (1.0e-7*tauGP/TK/TK)**(0.333d0)
  	ans = 1.0d0 - 0.0631789d0/TK + 0.115995d0/TK/TK - 0.401403d0/TS/TK
  	ans =ans+ 0.336463d0/TS/TK/TK
  	ans =ans/( 1.0d0 + 2.98394d0*xi + 1.53583d0*xi*xi + 3.85289d0*xi*xi*xi)
  	salpha_tilde=ans
end function Salpha_tilde

!Hirata 2006 equ 42
!checkd
function Tc_eff(TK,TS)
	implicit none
	real(dp),intent(in)::TK,TS
	real(dp)::ans,Tc_eff

  	ans = 1.0d0/TK + 0.405535d0/TK*(1.0d0/TS - 1.0d0/TK)
  	Tc_eff= 1.0d0/ans
end function Tc_eff


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!LYMAN ALPHA COUPLING!!!!!!!!!!!!!!!!

function frecycle(n)
integer,intent(in)::n
real(dp):: frecycle

  if (n== 0 ) then 
    frecycle= 1.d0
end if
  if (n== 1 ) then 
    frecycle= 1.d0
end if
  if (n== 2 ) then 
    frecycle= 1.d0
end if
  if (n== 3 ) then 
    frecycle= 0.d0
end if
  if (n== 4 ) then 
    frecycle= 0.2609d0
end if
  if (n==  5) then 
    frecycle= 0.3078d0
end if
  if (n== 6 ) then 
    frecycle= 0.3259d0
end if
  if (n== 7 ) then 
    frecycle= 0.3353d0
end if
  if (n== 8 ) then 
    frecycle= 0.3410d0
end if
  if (n== 9 ) then 
    frecycle= 0.3448d0
end if
  if (n== 10 ) then 
    frecycle= 0.3476d0
end if
  if (n==  11) then 
    frecycle= 0.3496d0
end if
  if (n== 12 ) then 
    frecycle= 0.3512d0
end if
  if (n==  13) then 
    frecycle= 0.3524d0
end if
  if (n==  14) then 
    frecycle= 0.3535d0
end if
  if (n== 15 ) then 
    frecycle= 0.3543d0
end if 
 if (n== 16 ) then 
    frecycle= 0.3550d0
end if
  if (n== 17 ) then 
    frecycle= 0.3556d0
end if
  if (n== 18 ) then 
    frecycle= 0.3561d0
end if
  if (n== 19 ) then 
    frecycle= 0.3565d0
end if
  if (n== 20 ) then 
    frecycle= 0.3569d0
end if
  if (n== 21 ) then 
    frecycle= 0.3572d0
end if
  if (n== 22 ) then 
    frecycle= 0.3575d0
end if
  if (n==  23) then 
    frecycle= 0.3578d0
end if
  if (n== 24 ) then 
    frecycle= 0.3580d0
end if
  if (n== 25 ) then 
    frecycle= 0.3582d0
end if
  if (n== 26 ) then 
    frecycle= 0.3584d0
end if
  if (n== 27 ) then 
    frecycle= 0.3586d0
end if
  if (n== 28 ) then 
    frecycle= 0.3587d0
end if
  if (n== 29 ) then 
    frecycle= 0.3589d0
end if
  if (n== 30 ) then 
    frecycle= 0.3590d0
end if  
end function frecycle



!lyman alpha flux from stars
!checked
subroutine J_alpha_star(z, Jalphastar, nbin, kbin, Walphak)
	implicit none					!adaptive method takes long time
	REAL(DP), INTENT(IN) :: z
	REAL(DP), INTENT(out) ::Jalphastar
	INTEGER, INTENT(IN) :: nbin
	REAL(DP), DIMENSION(nbin), INTENT(in) :: kbin
	REAL(DP), DIMENSION(nbin), INTENT(out) :: Walphak

	real(dp), dimension(nkbin) :: J0arr, J2arr
	real(dp)::res,delz,z1,z2,nu1,z_max,nu_p, zp, fn_f, k, r, q1, q2
	integer::nmax,i,n
	real(DP), parameter :: increase_f = 1.1d0

	nmax=23
	delz=0.05d0
	Jalphastar=0.d0
	Walphak=0.d0

	DO n=2,nmax
		res=0.d0
		z1=z
		z_max=zmax(z,n)
		nu1=nu_n(n)*nu_alpha
		fn_f=frecycle(n)

		DO 
			z2=z1+delz
			zp=(z1+z2)*0.5d0
			if(zp >= z_max) exit

			nu_p=nu1*(1.d0+zp)/(1.d0+z)

			res=delz*epsilon_cap_star(zp,nu_p)/hubble_constant(zp)*fn_f
			Jalphastar=Jalphastar+res

			r=comoving_distance(zp,z)	!Mpc
			r=r*hlittle	!Mpc/h
			DO i=1, nbin
				k=kbin(i)	!h/Mpc
				J0arr(i) = J0(k*r)
				J2arr(i)=J2(k*r)
			END DO

			q1=D(zp)/D(z)
			q2=1.d0+Bias(zp)
			Walphak=Walphak + q1*res*(q2*j0arr - 2.d0/3.d0*J2arr)


			z1=z2
			delz=delz*increase_f

		END DO

	END DO
	Walphak=Walphak/4.d0/pi*c_light*(1.d0+z)**2.d0
	Jalphastar=Jalphastar/4.d0/pi*c_light*(1.d0+z)**2.d0
	Walphak=Walphak/Jalphastar


end SUBROUTINE J_alpha_star


!from Pritchard & Furlanetto 2006 !barkana & loeb 2005b
function zmax(z,n)
	implicit none
	real(dp),intent(in)::z
	integer,intent(in)::n
	real(dp):: num, denom,zmax

  	num = 1.d0 - (n+1)**(-2.d0)
  	denom = 1.d0 - (n)**(-2.d0)
 	zmax=(1.d0+z)*num/denom - 1.d0

end function zmax


! Returns frequency of Lyman-n, in units of Lyman-alpha !barkana & loeb 2005b
function nu_n(n)
	implicit none
	integer,intent(in)::n
	real(dp):: ans,nu_n

  	ans = 1.0d0 - (n)**(-2.0d0)
  	nu_n =ans/0.75d0
 
end function nu_n

!comoving photon emissivity(# of photon emitted per unit comoving volume, per proper time and frequency) 
!from Pritchard & Furlanetto 2006
function epsilon_cap_star(z,nu)
	implicit none
	real(dp),intent(in)::z,nu
	real(dp)::epsilon_cap_star
	epsilon_cap_star=n_b0*f_star_param*dfcoll_dt(z)*epsilon_nu_star(nu)
end function epsilon_cap_star



!this model assumes mix of popii and popiii stars
function epsilon_nu_star(nu)
	implicit none
	real(dp),intent(in)::nu
	real(dp)::epsilon_nu_star,alpha_star,norm,a1, nlya, alpha_lyl
	real(dp) :: norm_popiii, norm_popii, res_popiii, res_popii


	norm_popii=n_lya_lyb_popII*alpha_lya_lyb_popII/(nu_beta**alpha_lya_lyb_popII-nu_alpha**alpha_lya_lyb_popII)
	if((nu >= nu_alpha) .and. (nu <= nu_beta)) then
		res_popii=norm_popii*nu**(alpha_lya_lyb_popII-1.d0)
	elseif(nu> nu_beta) then

		nlya=n_lya_lyl_popII-n_lya_lyb_popII
		nlya=(nlya*alpha_lya_lyb_popII+norm_popii*nu_beta**alpha_lya_lyb_popII)/norm_popii
		alpha_lyl=log(nlya)/log(nu_limit)	!!spectral index in the lyman limit energy range
		res_popii=norm_popii*nu**(alpha_lyl-1.d0)
	else
		res_popii=0.d0
	END IF


	epsilon_nu_star=res_popii
	
end function epsilon_nu_star



!offset of 21 cm brightness from cmb brightness temperature
!ok
function get_tb(z,ts,tgam,delta,xe)
	implicit none
	real(dp),intent(in)::z,ts,tgam,delta,xe
	real(dp)::get_tb
	get_tb=27.d0*((ts-tgam)/ts)*sqrt(0.15d0*(1.d0+z)/omega_m/10.d0/hlittle/ &
	hlittle)*(omega_b*hlittle*hlittle/0.023d0)*(1.d0+delta)*(1.d0-xe)
end function get_tb

!************************************SPIN TEMP END****************************

!***********************************FLUCTUATION************************

function Q_I(z,xe,gammaX, H)
	implicit none
	real(dp),intent(in)::z,xe,gammaX, H
	real(dp)::Q_I
	Q_I=(1.d0-xe)*gammaX/xe/(1.d0+z)/H
end function Q_I

function Q_R(z,xe,T, H)
	implicit none
	real(dp),intent(in)::z,xe,T, H
	real(dp)::Q_R
	Q_R=alpha_b(T)*C_clumping*xe*n_hi(z,0.d0)/(1.d0+z)/H
end function Q_R

function Q_C(z,xe,tk, H)
	implicit none
	real(dp),intent(in)::z,xe,tk, H
	real(dp)::Q_C
	Q_C=xe*((1.d0+z)**3.d0)*tcmb(z)/(1.d0+f_he+xe)/H/tk/(1.e+13*year_sec/8.55d0)
end function Q_C

function Q_X(z,gamma_heat,tk, H)
	implicit none
	real(dp),intent(in)::z,gamma_heat,tk, h
	real(dp)::Q_X
	Q_X=2.d0*gamma_heat/3.d0/Kboltz_erg/tk/(1.d0+z)/H
end function Q_X


!!!!!!!!!!!!!!!!!!LYMAN ALPHA part!!!!!!!

!sperical bessel function
!!
function j0(kr)
	implicit none
	real(dp),intent(in)::kr
	real(dp)::j0
	j0=sin(kr)/kr
end function j0

function j2(kr)
	implicit none
	real(dp),intent(in)::kr
	real(dp)::j2
	j2=(3.d0/kr/kr-1.d0)*sin(kr)/kr-3.d0*cos(kr)/kr/kr
end function j2




! *******************************POWER SPECTRUM******************

!power spectrum of tb angular averaged P(k)..karr in h/Mpc
!equation 39 of prichard..https://arxiv.org/abs/astro-ph/0607234
!ok
SUBROUTINE cal_PS(ndim, karr,z,beta,beta_x,beta_alpha,beta_t,xx, PS_Tb3d, PS_T, PS_tb2)
	implicit none
	integer, intent(in) :: ndim
	real(dp), dimension (ndim), intent(in)::karr
	real(dp),intent(in)::z,beta,beta_x,beta_alpha,beta_t
	real(dp),dimension(:),intent(in)::xx
	real(dp),dimension(ndim),intent(out):: PS_Tb3d, Ps_T,  PS_tb2
	real(dp)::Pk, k,x_e, growth, con
	real(dp), dimension(ndim) :: g_e,g_T,  beta_p
	integer:: ii 
	real(dp),dimension(ndim) :: W_a, con_arr, Pk_arr, beta_p_arr, norm
	integer :: n1, n2
	
	n1=3+nkbin
	n2=n1+nkbin


	x_e=xx(1)
	g_T=xx(4:n1)	
	g_e=xx(n1+1:n2)


	PS_tb3d=0.0
	PS_T=0.0
	PS_tb2=0.0
	beta_p_arr=0.0
	W_a=0.0
	Pk_arr=0.0
	con_arr=0.0


	DO ii=1, ndim
		k=karr(ii)


		Pk_arr(ii) = pspec_cdm(k)
		norm = norm_cdm(sigma_8)

		con=k*k*k
		con_arr(ii)=con
		W_a(ii) = (Walphak_comm(ii)*Jalpha_star+Wxk_comm(ii)*Jalpha_x)/(Jalpha_star+Jalpha_x)
	END DO
	con_arr=con_arr/2.d0/pi/pi
	growth=d(z)
	Pk_arr=Pk_arr*norm*growth**2.d0

	beta_p=beta - beta_x*x_e*g_e(:)/(1.d0+x_e) &
	       + beta_T*g_t(:) 

	beta_p_arr=W_a*beta_alpha + beta_p 


	PS_tb3d=(beta_p_arr)**2.d0
	PS_tb3d=PS_tb3d+ 2.d0*beta_p_arr/3.d0
	PS_tb3d=PS_tb3d+1.d0/5.d0
	PS_Tb3d=PS_tb3d*Pk_arr*con_arr
	!PS_Tb3d=sqrt(PS_Tb3d)

	PS_T=Pk_arr*con_arr*g_t*g_t
	!PS_T=sqrt(PS_T)

	PS_tb2=2.d0*beta_p_arr*Pk_arr*con_arr
	!PS_tb2=sqrt(PS_tb2)

end subroutine cal_PS


!https://arxiv.org/abs/astro-ph/0607234
FUNCTION Cal_beta(x_c,x_a)
	IMPLICIT NONE
	real(dp),intent(in)::x_c,x_a
	real(dp)::Cal_beta,x_total
	x_total=x_c+x_a
	Cal_beta=1.d0+ x_c/x_total/(1.d0+x_total)
END FUNCTION Cal_beta

!https://arxiv.org/abs/astro-ph/0607234
FUNCTION Cal_beta_x(x_c,x_a,xc_h,xc_e)
	IMPLICIT NONE
	real(dp),intent(in)::x_c,x_a,xc_h,xc_e
	real(dp)::Cal_beta_x,x_total
	x_total=x_c+x_a
	Cal_beta_x=1.d0+(xc_h-xc_e)/x_total/(1.d0+x_total)
END FUNCTION Cal_beta_x

!https://arxiv.org/abs/astro-ph/0607234 eq 6
!
FUNCTION Cal_beta_a(x_c,x_a)
	IMPLICIT NONE
	real(dp),intent(in)::x_c,x_a
	real(dp)::Cal_beta_a,x_total
	x_total=x_c+x_a
	Cal_beta_a=x_a/x_total/(1.d0+x_total)
END FUNCTION Cal_beta_a

!https://arxiv.org/abs/astro-ph/0607234 equation  6
!
FUNCTION Cal_beta_T(x_c,x_a,z,xc_h,xc_e,T)
	IMPLICIT NONE
	real(dp),intent(in)::x_c,x_a,z,xc_h,xc_e,T
	real(dp)::Cal_beta_T,T_g,x_total,a1,a2, Tp
	x_total=x_c+x_a
	T_g=tcmb(z)

if(abs(T-T_g)<1.d0) then
if(T>T_g) Tp=T_g+1.d0
if(T<T_g) Tp=T_g-1.d0
else
Tp=T
end if

if(Tp>5d3) Tp=5d3
	a1=dLogKe_dLogTk(Tp)
	a2=dLogKh_dLogTk(Tp)



	Cal_beta_T=T_g/(Tp-T_g) + (1.d0/x_total/(1.d0+x_total))*(xc_e*a1 + xc_h*a2)


END FUNCTION Cal_beta_T

FUNCTION dLogKe_dLogTk(T)
	implicit none
	real(dp),intent(in)::T
	real(dp)::dLogKe_dLogTk,a1,a2,T1,T2,delt
	delt=T*0.1d0
	T1=T+ delt
	T2=T- delT
	a1=kappa_10_elec(T1,0)-kappa_10_elec(T2,0)
	a2=2.d0*delT
	dLogKe_dLogTk=(a1/a2)*(T/kappa_10_elec(T,0))
END FUNCTION dLogKe_dLogTk

FUNCTION dLogKh_dLogTk(T)
	implicit none
	real(dp),intent(in)::T
	real(dp)::dLogKh_dLogTk,a1,a2,T1,T2,delt
	delt=T*0.1d0
	T1=T+ delt
	T2=T- delT
	a1=kappa_10(T1,0)-kappa_10(T2,0)
	a2=2.d0*delT
	dLogKh_dLogTk=(a1/a2)*(T/kappa_10(T,0))
END FUNCTION dLogKh_dLogTk

!*******************GALAXY BIAS***********************
!!Interplated galaxy bias 
function bias(z)
	implicit none
	real(dp),intent(in)::z
	real(dp)::bias,mass1,r1,sig,nu, zp, zmid
	integer::i, imin
	real(dp), dimension(ndim_fcoll) :: x, y

	x=av_bias_arr(1,:)
	y=av_bias_arr(2,:)
	bias=linear_interpol(z,x,y)

end function bias

!!Halo Mass averaged  bias
! 
FUNCTION mass_av_bias(z)
	IMPLICIT NONE
	real(dp), intent(in) :: z
	real(dp) :: mass_av_bias, M1, M2, M, res, dndm, res1, bb, delm


	M1=mass_min_cdm(t_vir_param,z)	!!M0/h


	res=0.d0
	res1=0.d0
	DO 
		M2=M1*1.2d0
		M=(M1+M2)*0.5d0
		delm=M2-M1
		if(M>1d14) exit
		dndm=numdenm_cdm(M,z)
		bb=bias_m(M,z)
		res=res+dndm*bb*M*delm
		res1=res1+dndm*M*delm
		M1=M2
	END DO
	if(res1<=0.d0) then
		mass_av_bias=0.d0
	else
		mass_av_bias=res/res1
	end if
END FUNCTION mass_av_bias


!!m in M0/h
!http://iopscience.iop.org/article/10.1086/319797/pdf footnote 3
FUNCTION bias_m(m,z)
	implicit none
	real(dp),intent(in)::m,z
	real(dp)::bias_m,nu,sig,r1,delc
	r1=MtoRadius(m)
	sig=sigma_cdm(r1)
	nu=delta_c/d(z)/sig
	bias_m=1.d0+(nu*nu-1.d0)/delta_c

END FUNCTION bias_m


FUNCTION bias_m_cdm(m,z)
	implicit none
	real(dp),intent(in)::m,z
	real(dp)::bias_m_cdm,nu,sig,r1,delc
	r1=MtoRadius(m)

	sig=sigma_cdm(r1)
	delc=delta_c
	nu=delc/d(z)/sig
	bias_m_cdm=1.d0+(nu*nu-1.d0)/delc

END FUNCTION bias_m_cdm


end module funct
