module subr
	use nrtype
	use param
	use funct
	use cosmo

contains

SUBROUTINE declare
	IMPLICIT NONE

	write(*,*) 'Cosmological parameters:: omega_l, omega_m, omega_b, h, sigma_8', omega_l, omega_m, omega_b, hlittle, sigma_8
	write(*,*) 'X-ray properies :: f_X, Alpha_s', f_X_param, Alpha_s_param
	write(*,*) 'UV properties :: f_esc, N_ion, f_star', f_esc_param, N_ion_param, f_star_param
	if(sheth_tormen==1) then
	write(*,*) 'Using ST halo mass function'
	elseif(sheth_tormen==0) then
	write(*,*) 'Using PS halo mass function'
	else
	write(*,*) 'Choose right halo mass function'
	stop
	end if

END SUBROUTINE declare

SUBROUTINE  update_variables(var_arr)
	IMPLICIT NONE
	real(dp), dimension(:), intent(inout) :: var_arr

	xe_comm = var_arr(1)
	xi_comm = var_arr(2)
	if(var_arr(3)>1d4) var_arr(3)=1d4	!!by hand included as T is going very high for large fx and alpha. cooling should be incorporated
	TK_comm = var_arr(3)

END SUBROUTINE  update_variables

SUBROUTINE det_filename
use param
IMPLICIT NONE


	write(f_esc_s,'(f12.3)') f_esc_param
	f_esc_s=adjustl(f_esc_s)

	write(N_ion_s,'(f12.3)') N_ion_param
	N_ion_s=adjustl(N_ion_s)

	!write(f_star_s,'(f12.3)') f_star_param
	!f_star_s=adjustl(f_star_s)

	write(f_X_s,'(f12.3)') f_X_param
	f_X_s=adjustl(f_X_s)

	write(Alpha_s_s,'(f12.3)') ALpha_s_param
	Alpha_s_s=adjustl(Alpha_s_s)



	write(alpha_star_pop_s,'(f12.3)') alpha_lya_lyb_popII
	alpha_star_pop_s=adjustl(alpha_star_pop_s)

	write(nlya_pop_s,'(f12.3)') n_lya_lyb_popII
	nlya_pop_s=adjustl(nlya_pop_s)

	write(t_vir_s,'(f12.3)') t_vir_param
	t_vir_s=adjustl(t_vir_s)

	filename_comm="Nion"//trim(N_ion_s)//"_tvir"//trim(t_vir_s) &
	& //"_fx"//trim(f_X_s)//"_alx"//trim(Alpha_s_s)


END SUBROUTINE det_filename



!!Used to calculate the density threshold in equation 4 of FZH04
!
FUNCTION k_factor(ion_ef)
	IMPLICIT NONE
	real(dp), intent(in) :: ion_ef
	real(dp) :: k_factor, x, iof

	x=1.d0-(1.d0+n_rec_comm)/ion_ef
	if(x>1.d0 .or. x<-1.d0) then
		!write(*,*) 'x is outside limit -1, 1 for inverf(x), ion_ef', x, iof
		k_factor=-1.d0	!!by hand putted

	else
		k_factor=inverf(x)
	end if

END FUNCTION k_factor

!!Inverse error function..
!
FUNCTION inverf(x)
	IMPLICIT NONE
	real(dp), intent(in) ::x
	real(dp) ::yii, yi, dely, accu, inverf
	integer::cn
	integer, parameter :: cnmax = 200000

	accu=1d-6*abs(x)
	yi=x

	cn=0
	DO 
		cn=cn+1
		if(cn>cnmax) then
			write(*,*) 'Error!!..loop overshoot in inverf', x, yi, erf(yi)
			stop
		end if
		yii=erf(yi)
		if(abs(yii-x)<accu) exit
		dely=(yii-x)
		yi=yi-dely

	END DO

	inverf=yi


END FUNCTION inverf

!!calculate collapsed fraction in grid..that will be used for calculating the integrals late.
!
subroutine init_arrfcoll
	implicit none
	integer::i
	real(dp)::z, fcoll, delz, fc, m, sig, r, cor, avbi, x

	
	delz=(z_star+1.d0-z_end+1.d0)/(dble(ndim_fcoll)-2.d0)
	delz_fcoll=delz
	z=z_end-1.d0
	DO i=1, ndim_fcoll
		arr_fcoll(1,i)=z

		fc=fcoll_cdm(z)

		arr_fcoll(2,i) = fc

!!average bias
		avbi=mass_av_bias(z)
		av_bias_arr(1,i)=z
		av_bias_arr(2,i)=avbi

		z=z+delz
	END DO

!! mass bin, sigma(M) bin
	call logspace_dp(massdim_comm, mmin_bin, mmax_bin, M_arr_comm)

	DO i=1, massdim_comm
		M=M_arr_comm(i)
		sig=sigma_m_CDM(M)
		sigma_arr_comm(i)=sig
	    	x=3.d0*M/(4.d0*pi*rho_c_h)
	    	x=x**(1.d0/3.d0)
		BSD_R_arr(i) = x
	END DO

!! r bin

	call logspace_dp(nrdim_comm, rmin_bin, rmax_bin, r_arr_comm)

!!
	k_zeta_comm=k_factor(ion_ef_comm)

end subroutine init_arrfcoll





!!Initial conditions of the differential equation solver..
!!
SUBROUTINE initialize(z,x)
	IMPLICIT NONE
	real(dp),intent(in)::z
	real(dp),dimension(:),intent(inout)::x
	real(dp)::zdec
	integer :: n1, n2
	
	n1=3+nkbin
	n2=n1+nkbin

	zdec=z_decoupling()

	x(1)=1.e-4
	x(2)=1.e-4
	x(3)=2.725d0*(1.0d0+z)*(1.0d0+z)/(1.d0+zdec)
	x(4:n1)=2.d0/3.d0	!!g_T(k)
	x(n1+1:n2)=-1.d0	!!g_e(k)

END SUBROUTINE initialize




!!Right hand sides of the differential equations..
!!
SUBROUTINE derivs(xd,yd,dydxd)
	IMPLICIT NONE
	REAL(dP), INTENT(IN) :: xd
	REAL(dP), DIMENSION(:), INTENT(IN) :: yd
	REAL(dP), DIMENSION(:), INTENT(OUT) :: dydxd

	real(dp)::n_B,n_H, fex, fheat, xe, xi, Tk, dtdz, z, hub
	real(dp)::x_heat,jalpha
	real(dp)::QI,QR,QC,QX,gammaX,gamma_heat,gam_e
	real(dp)::WX
	integer :: n1, n2
	
	n1=3+nkbin
	n2=n1+nkbin

	z=xd
	xe=yd(1)	!!partial
	xi=yd(2)
	Tk=yd(3)


	n_B=n_b0*(1.d0+z)**3.d0 
	n_H=f_H*n_B 
	hub=hubble_CONSTANT(z)

	fex=f_ex(xe)
	fheat=f_heat(xe)
	dtdz=dt_dz(z)

	call  Xray_heating(z,xe, x_heat, nkbin, kbin_comm, Wxk_comm)
	Jalpha_x=c_light/4.d0/pi/hplanck_erg/nu_alpha/nu_alpha/hub*x_heat*fex*p_alpha/fheat


	call J_alpha_star(z, Jalpha_star, nkbin, kbin_comm, Walphak_comm)
	jalpha=Jalpha_x+Jalpha_star


	gammaX=x_heat*erg_ev*fex/e_hi/fheat/n_b0/(1.d0+z)**3.d0

	dydxd(1)=((1.d0-xe)*gammaX -alpha_B(Tk)*C_clumping*xe*xe*n_H)*dtdz
	if(xi >= 1.d0) then
		dydxd(2)=0.d0
	else
	dydxd(2)=((1.d0-xe)*gamma_i(z) -alpha_A*C_clumping*xi*xi*n_H)*dtdz
	end if
	dydxd(3)=(dcompton_heating_dt(z,Tk,xe)+x_heat*2.d0/3.d0/Kboltz_erg/(1.d0+xe)/n_b)*dtdz &
	+2.0/3.0*Tk*dydxd(1)/(1.d0+xe)- 2.*Tk*hub*dtdz !2.*Tk/(1.d0+z)	!!last term is hubble cooling


	gamma_heat=x_heat/n_b/(1.d0+xe)


	QI=Q_I(xd,xe,gammaX, hub)
	QR=Q_R(xd,xe,Tk, hub)
	QC=Q_C(xd,xe,Tk, hub)
	QX=Q_X(xd,gamma_heat,Tk, hub)


	dydxd(4:n1)=(yd(4:n1)-2.d0/3.d0)/(1.d0+xd)-QX*(Wxk_comm(:)-yd(4:n1))-QC*yd(4:n1)
	dydxd(n1+1:n2)=yd(n1+1:n2)/(1.d0+xd)-QI*(Wxk_comm(:)-yd(n1+1:n2))+QR*(1.d0+yd(n1+1:n2))

END SUBROUTINE derivs

SUBROUTINE show_time(ifile)
	IMPLICIT NONE
	character(*), intent(in) :: ifile
	real(4) :: time

	call cpu_time(time)
	write(*,*) ' Time:', time, ifile

END SUBROUTINE show_time


!!https://arxiv.org/abs/astro-ph/0408408
SUBROUTINE nrec_MH(z, Tvir, nrec)
	implicit none
	real(dp), intent(in) :: z, Tvir 
	real(dp), intent(out) :: nrec
	real(dp) :: Mj, Mmax, M1, M2, M, delM, res, crec, dndlnm
	real(dp), parameter::increase_f=1.2d0

	Mj=jeans_mass_cdm(z)
	Mmax = mass_min_cdm(Tvir,z)


	M1=MJ
	res=0.d0

	DO
		M2=M1*increase_f
		M=(M1+M2)*0.5d0
		delM=M2-M1
		crec=conrate_nrec(z, M)

		dndlnm=int_fcoll_cdm(M,z)	!!keep this to cdm.. changing fdm is not working prperly

		res=res+dndlnm*delM*crec
		if(M2>Mmax) exit
		M1=M2
	END DO

	res=res/rho_c_h
	nrec=res

END SUBROUTINE nrec_MH


!https://arxiv.org/abs/astro-ph/0408408
!M in M0/h
FUNCTION conrate_nrec(z, M)
implicit none
real(dp), intent(in) :: z, M
real(dp) :: conrate_nrec, A, B, C, D, E, F, G, M7, F0
	A=4.4
	B=0.334
	C=0.023
	D=0.199
	E=-0.042
	F=0.
	G=1.
	F0=1.
	M7=(M/hlittle/1d7)
	conrate_nrec = A*M7**(B+C*log10(M7))*F0**(D+E*log10(F0))*(F+G*(1.+z)/10.)

END FUNCTION conrate_nrec



end module subr
