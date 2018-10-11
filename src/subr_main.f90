module subr_main
contains
subroutine linspace_dp(n,xmin,xmax,x)
    	implicit none
	integer, intent(in) :: n
    	real(8),intent(in) :: xmin,xmax
    	real(8),dimension(n),intent(out) :: x
    	integer :: i

    	if (n == 1) then
	       	if(xmin /= xmax) then
		  	write(0,'("ERROR: Cannot call linspace with n=1 and xmin /= xmax")')
		  	stop
	       	else
        		x = xmin
       		end if
    	else
       		do i=1,n
          		x(i) = (xmax-xmin) * dble(i-1) / dble(n-1) + xmin
       		end do
    	end if
end subroutine linspace_dp

subroutine logspace_dp(n,xmin,xmax,x)
    	implicit none
	integer, intent(in) :: n
    	real(8),intent(in) :: xmin,xmax
    	real(8), dimension(n),intent(out) :: x
   	 if (size(x) == 1 .and. xmin /= xmax) then
       		write(0,'("ERROR: Cannot call logspace with n=1 and xmin /= xmax")')
       		stop
    		end if
    	call linspace_dp(n,log10(xmin),log10(xmax),x)
    	x = 10.d0**x
end subroutine logspace_dp



SUBROUTINE reionization_21cm(arr, zcut, z_arr, GLB_arr, k_arr, PSKZ_arr, BSD_R_out, BSDZ_arr)
	USE nrtype
	USE nr , ONLY : odeint,rkqs
	USE ode_path 
	use subr					
	use funct
	use param
	use subr_FZH
	use subr_plot
	IMPLICIT NONE
	real(dp), dimension(nparam), intent(in) ::arr
	real(dp), intent(in) :: zcut
	real(dp), dimension(nzloop), intent(out) :: z_arr
	real(dp), dimension(nzloop, 4), intent(out) :: GLB_arr
	real(dp), dimension(nkbin), intent(out) :: k_arr
	real(dp), dimension(nzloop, nkbin), intent(out) :: PSKZ_arr
	real(dp), dimension(massdim_comm), intent(out) :: BSD_R_out
	real(dp), dimension(nzloop, massdim_comm), intent(out) :: BSDZ_arr
	real(dp)::  tbmin, zmintb, psmax, zpsmax
	real(dp)::  z, zp, z2, tbp
	real(dp), dimension(nkbin) :: psk1, psk2, pskmiss
	REAL(dP):: eps, hmin, h
	real(dp),dimension(3+2*nkbin):: var_arr	!!3 global variable, rest are k dependent 
	real(dp):: Trad, xc, Tceff, xa, xc_h, xc_e, T_spin, T_B, Jalpha, psz9, psz8, psz10, Qp
	real(dp):: beta, beta_x, beta_a, beta_T, kvec, Q, avtb_inhomo, fc, pt, abtbmax
	integer :: i, iloop, ii1, ii2, set_con
	character(180) ::  ifile

	real(dp), dimension(nkbin):: ptot, PS_Tb3d, PS_T, PS_tb2

	integer, dimension(nzloop) :: track_missps

	double precision::z1,tau2, tau1, nh, nhe, con
	double precision::RELACC,ABSACC,acc
	integer:: MAXRUL, IPARM, ifail, N
	double precision,dimension(390)::alpha

	z_arr=0.d0
	GLB_arr=0.d0
	PSKZ_arr=0.d0
	BSDZ_arr = 0.d0
	k_arr=0.d0


	relacc=1.0D-9
	absacc=0.0d0
	maxrul=9
	iparm=0

!! Imitializing the common parameters

	f_X_param = arr(1)
	Alpha_s_param=arr(2)
	N_ion_param=arr(3)*4000.d0
	t_vir_param = arr(4)*10000.d0

	f_esc_param = f_esc_fiducial
	f_star_param = f_star_fiducial
	!t_vir_param = t_vir_fiducial


	ion_ef_comm=	f_esc_param *N_ion_param *f_star_param 
	call det_filename

!!!

	call init_kbin(nkbin, kbin_comm)
	k_arr=kbin_comm
	call declare
	call Initialize(z_star,var_arr)
	call init_arrfcoll


!! these files will contain global variables, contribution from different fluctuation and power spectrum at k=0.1 h/Mpc


!!Simulation starts here..


	Qp=0.d0
	h=-1.d0*delz_eor/3.d0
	zp=z_cosmic_dawn
	z=zp-delz_eor
	eps=1d-4
	hmin=0.0d0

	xi_comm=0.d0
	track_missps=0

	T_spin=tcmb(z)	!!initializing with cmb
	tbp=0.d0

	DO iloop=1, nzloop	!!redshift lop

	z_arr(iloop)= z
if((xi_comm+xe_comm)>=1.d0) then
	GLB_arr(iloop,:) =(/var_arr(1)+var_arr(2), var_arr(3),t_spin,0.d0/)
	PSKZ_arr(iloop,:)=0.d0
else

		call init_xhii_inhomogeneity(z)
		if(zp <= z_eor_end) exit

	
		call odeint(var_arr,zp,z,eps,h,hmin,derivs,rkqs)
		call update_variables(var_arr)


		Jalpha=Jalpha_star+Jalpha_x
		Trad=tcmb(z)
		xc = xcoll(z, TK_comm, 0.d0, xe_comm)
		xa = xalpha_tilde(z,Jalpha,TK_comm,T_spin,0.d0,xe_comm)
		xc_h = xcoll_h(z,TK_comm,0.d0,xe_comm)
		xc_e = xcoll_e(z,TK_comm,0.d0,xe_comm)
		Tceff = Tc_eff(TK_comm,T_spin)


		T_spin = get_Ts(Trad,TK_comm,Tceff,xa,xc)
		t_b = get_tb(z,t_spin,Trad,0.d0,xi_comm)


		write(27,'(17f20.4)') z,var_arr(1:3),t_spin,tcmb(z),t_b
		flush(27)


		GLB_arr(iloop,:) = (/var_arr(1)+var_arr(2), var_arr(3),t_spin,t_b/)


!!fluctuations..


	if(abs(TK_comm-tcmb(z)) > 1.d0) then
		gT_comm = var_arr(4:3+nkbin)
		xa = xalpha_tilde(z,Jalpha,TK_comm,T_Spin,0.d0,xe_comm)
		beta = Cal_beta(xc,xa)
		beta_x = Cal_beta_x(xc,xa,xc_h,xc_e)
		beta_a = Cal_beta_a(xc,xa)
		beta_T = Cal_beta_T(xc,xa,z,xc_h,xc_e,TK_comm)


		call cal_PS(nkbin, kbin_comm,z,beta,beta_x,beta_a,beta_t,var_arr, PS_Tb3d, PS_T, PS_tb2)
		dimless_PS = T_b*T_b*PS_Tb3d

!!xhii part

			if(xi_comm<1.d0) then
				 avtb_inhomo=abs(t_b)/(1.d0-xi_comm)
			else
				avtb_inhomo=0.d0
			end if

			call  PS_TB_ion(nkbin, kbin_comm, z, ion_ef_comm, Ptot, Q)

			if(Q<Qp) then
				Ptot=0.d0	!condition applied.. this is due to the fact that Q from furlanetto goes down above xhii>.8
			else
				call cal_VdndlnR_hii(z, ion_ef_comm, BSD_arr)
			end if
		dimless_PS = dimless_PS + avtb_inhomo**2.*Ptot



		if(Q>Qp) Qp=Q	!!this prevent wrong PS at end stages of EoR


	else
		if(tbp<t_b) track_missps(iloop) = 1
		tbp=t_b

	end if


		PSKZ_arr(iloop,:)=dimless_PS(:)
		BSDZ_arr(iloop,:)=BSD_arr

	

		zp=z
		z=zp-delz_eor

if(zp<=zcut) then
if(wr_messages) write(*,*) 'Exit limit reached, zcut=', zcut
exit
end if

		!inter_arr=(/tbmin, zmintb, psmax, zpsmax, psz8, psz9, psz10/)
end if
	end do

!!Now interpolate the missed redshifts
	set_con=0
	ii1=nzloop
	DO i=1, nzloop
	!if(track_missps(i)==0 .and. set_con==0) set_con=1
	if(track_missps(i)==1 ) then
	ii1=i-1
	exit
	end if
	END DO
	ii2=ii1+sum(track_missps)+1

if(ii1>=1 .and. ii2<nzloop) then



	z1=z_arr(ii1)
	psk1=PSKZ_arr(ii1,:)
	z2=z_arr(ii2)
	psk2=PSKZ_arr(ii2,:)

	DO i=1, nzloop
	if(track_missps(i)==1) then
	z=z_arr(i)
	pskmiss=psk2 + (psk1-psk2)*(z-z2)/(z1-z2)
	PSKZ_arr(i,:)=pskmiss

	end if
	END DO
END IF

BSD_R_out = BSD_R_arr

!!writting part

if(write_output_file) then
	ifile = output_path//trim(filename_comm)//".bin"
if(wr_messages)	write(*,*) ifile
	open(unit=27,file=ifile, status='replace', form='unformatted')
	write(27) nzloop, nkbin, 4, massdim_comm
	write(27) real(z_arr)
	write(27) real(kbin_comm)
	write(27) real(BSD_R_arr)
	write(27) real(GLB_arr)
	write(27) real(PSKZ_arr)
	write(27) real(BSDZ_arr)
	close(27)
end if
!!

end subroutine reionization_21cm
end module subr_main


