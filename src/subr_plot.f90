module subr_plot
use param
use subr_FZH
contains

subroutine Plot_sig_dndm_fcoll(filename)
IMPLICIT NONE
	character(*), intent(in) :: filename
	integer :: i
	character(180) :: ifile
	real(dp) :: M, z

!!write collapsed fraction.. 

	ifile = output_path//trim(filename)//"fcoll.dat"
	write(*,*) ifile
	open(unit=2,file=ifile)
	write(2,*) '#z, fcoll'

	DO i=1, ndim_fcoll
	z=arr_fcoll(1,i)
	write(2,'(17e20.4)') arr_fcoll(:,i), fcoll_analytical(z)
	END DO
	close(2)




!!dndlm for the scenario..
 
	ifile = output_path//trim(filename)//'z_6dndlnm.dat'
	open(unit=1, file=ifile)
	write(1,*) '#M(M0/h), dndlnm'
	write(*,*) ifile
	M=1d8
	z=6.0d0
	DO 

	write(1,*) M, int_fcoll_cdm(M, z)

	M=M*1.5d0
	if(M>1d16) exit
	END DO
	close(1)

!!sigma

	call plot_sigmam(filename)



end subroutine Plot_sig_dndm_fcoll



!!
SUBROUTINE plot_sigmam(filename)
	IMPLICIT NONE
	character(*), intent(in) :: filename
	real(dp) ::  M, M1, sig, res, ion_ef, M2, B, B0, B1, sigmin
	real(dp), parameter::increase_f=1.2d0
	character(12) ::zt_s
	character(180) :: ifile




	M=1d2 !in M0/h
	M1=M

	ifile = output_path//trim(filename)//"sigmam.dat"
	write(*,*) ifile
	open(unit=219,file=ifile)
	write(219,*) "#M/hlittle, sig**2.0, sig"

	DO
		M2=M1*increase_f
		sig=sigma_m_cdm(M2)

		write(219,'(17e20.4)') m2, sig**2.0, sig

		if(M2>1d14) exit
		M1=M2
	END DO
	close(219)

END SUBROUTINE  plot_sigmam

subroutine plot_window(z, filename)
	IMPLICIT NONE
	real(dp), intent(in) ::z
	character(*), intent(in) :: filename
	integer :: i
	character(12) ::zt_s
	character(180) :: ifile

	write(zt_s,'(f12.3)') z
	zt_s=adjustl(zt_s)

	ifile = output_path//trim(filename)//'_z'//trim(zt_s)//"windowf.dat"
	!write(*,*) ifile
	open(unit=215,file=ifile)
	write(215,*) "#k(h.Mpc), WX, W_a, gT"

	DO i=1, nkbin
	write(215, '(17e20.4)') kbin_comm(i), Wxk_comm(i), Walphak_comm(i), GT_comm(i)
	END DO
	close(215)

end subroutine plot_window

!!
SUBROUTINE cal_VdndlnR_hii(z, ion_ef, arr)
	IMPLICIT NONE
	real(dp), intent(in) ::z, ion_ef
	real(dp), dimension(massdim_comm), intent(out) :: arr
	real(dp) :: M, M1, sig, res, M2, h,  x, dndm, Q, mmin
	real(dp), parameter::increase_f=1.2d0
	character(12) ::zt_s
	character(180) :: ifile
	integer :: i


	Q=vol_frac_Q(z, ion_ef)


	Mmin=mminz_comm

    	h=hlittle
    	!rho_0=rho_c*h*h*omega_m  !rho_c in M0/mpc3 ! m is in M0/h, x is in mpc/h,rho_0 is in (M0/h)/(mpc/h)^3

	arr=0.d0

	DO i=1, massdim_comm
	M=M_arr_comm(i)
	if(M>Mmin) then
	    		x=3.d0*M/(4.d0*pi*rho_c_h)
	    		x=x**(1.d0/3.d0)
			dndm= dndlnm_arr_comm(i)
			res=dndm*3.d0*M**2.d0/rho_c_h/Q
			arr(i) = res
			!BSD_R_arr(i) = x
	ENDif
	END DO


!write(*,*) dndlnm_arr_comm, q
!STOP


END SUBROUTINE  cal_VdndlnR_hii

!!
SUBROUTINE plot_VdndlnR_hii(z, ion_ef, filename)
	IMPLICIT NONE
	real(dp), intent(in) ::z, ion_ef
	character(*), intent(in) :: filename
	real(dp) :: M, M1, sig, res, M2, h,  x, dndm, Q, mmin
	real(dp), parameter::increase_f=1.2d0
	character(12) ::zt_s
	character(180) :: ifile
	integer :: i

	write(zt_s,'(f12.3)') z
	zt_s=adjustl(zt_s)

	Q=vol_frac_Q(z, ion_ef)


	Mmin=mminz_comm

    	h=hlittle
    	!rho_0=rho_c*h*h*omega_m  !rho_c in M0/mpc3 ! m is in M0/h, x is in mpc/h,rho_0 is in (M0/h)/(mpc/h)^3

	ifile = output_path//trim(filename)//'_z'//trim(zt_s)//"VdndlnR.dat"
	!write(*,*) ifile
	open(unit=219,file=ifile)
	write(219,*) "#M(M0/h), x(Mpc/h), VdndlnR/Q , Q"

	DO i=1, massdim_comm-1
M=M_arr_comm(i)
if(M>Mmin) then
    		x=3.d0*M/(4.d0*pi*rho_c_h)
    		x=x**(1.d0/3.d0)
		dndm= dndlnm_arr_comm(i)
		res=dndm*3.d0*M**2.d0/rho_c_h/Q

		write(219,'(17e20.4)') M, x,  res, Q, dndm, B0z_comm, B1z_comm

ENDif
	END DO
	close(219)

!write(*,*) dndlnm_arr_comm, q
!STOP


END SUBROUTINE  plot_VdndlnR_hii

!!
SUBROUTINE plot_Qz(ion_ef, filename)
	IMPLICIT NONE
	real(dp), intent(in) :: ion_ef
	character(*), intent(in) :: filename
	real(dp) ::  z_max, z_min, del_z, x1, x2, Q, cxx, Qp, Q1, av_halo_bias
	real(dp), parameter::increase_f=1.2d0
	character(12) ::zt_s
	character(180) :: ifile

	z_max=25.d0
	z_min=6.d0
	del_z=-0.5d0
	x1=z_max
	x2=x1+del_z


	ifile = output_path//trim(filename)//"Q_z.dat"
	write(*,*) ifile
	open(unit=219,file=ifile)
	write(219,*) "#z, Q "

	do 
		if(x1 <= z_min) exit
		Q=vol_frac_Q(x1, ion_ef)

		!if(Qp>Q) exit
		write(219,'(17f20.4)') x1, Q 


		Qp=Q
		x1=x2
		x2=x2+del_z

	end do
	close(219)



END SUBROUTINE  plot_Qz

!!
SUBROUTINE plot_correl(ion_ef, z, filename)
	IMPLICIT NONE
	real(dp), intent(in) ::z, ion_ef
	character(*), intent(in) :: filename
	real(dp) ::  cxx, r12, cdd, cdx, Q, rmin, rmax
	real(dp), parameter::increase_f=2.d0
	character(12) ::zt_s
	character(180) :: ifile
	integer, parameter :: nrdim = 50
	real(dp), dimension(nrdim) ::rarr, xi_dd, xi_xx, xi_xd
	integer :: i


	Q=vol_frac_Q(z, ion_ef)


	rmin=1d-2
	rmax=1d2

	call logspace_dp(nrdim,rmin,rmax,rarr)	!!Mpc/h
	call Initialize_correlation_dd_xx_xd(z, ion_ef, Q, nrdim, rarr, xi_dd, xi_xx, xi_xd)

	ifile = output_path//trim(filename)//"cor_dd.dat"
	write(*,*) ifile
	open(unit=219,file=ifile)
	write(219,*) "#R (Mpc/h), xi_dd, xi_xx, xi_xd"

	DO i=1, nrdim
		write(219,'(17f20.4)') rarr(i), xi_dd(i), xi_xx(i), xi_xd(i)
		flush(219)
	END DO
	close(219)



END SUBROUTINE  plot_correl

subroutine plot_nfw(z, M)
implicit none
	real(dp), intent(in) ::z, M
	real(dp) ::  cxx, r12, cdd, cdx, Q, rmin, rmax, k
	real(dp), parameter::increase_f=1.2d0
	character(12) ::zt_s
	character(180) :: ifile
	integer, parameter :: nrdim = 50
	real(dp), dimension(nrdim) ::rarr, xi_dd, xi_xx, xi_xd
	integer :: i

	rmin=1d-5
	rmax=1d2


	ifile = output_path//"nfw.dat"
	write(*,*) ifile
	open(unit=219,file=ifile)
write(219,*) "#R (kpc), rho_nfw (M0/pc**3)"

	DO 
		write(219,'(17e20.4)') rmin*1000.0/hlittle, nfw_profile(M, rmin, z)*hlittle**2.0/1d18
rmin=rmin*increase_f
		flush(219)
if(rmin>rmax) exit
	END DO
	close(219)

!stop
!-----


k=0.1d0
	ifile = output_path//"u_nfw.dat"
	write(*,*) ifile
	open(unit=219,file=ifile)
write(219,*) "#k (h/Mpc), u_nfw"

	DO 
		write(219,'(17e20.4)') k, u_NFW(k, m, 0.d0)
k=k*increase_f
		flush(219)
if(k>1d3) exit
	END DO
	close(219)

stop
end subroutine plot_nfw





end module subr_plot
