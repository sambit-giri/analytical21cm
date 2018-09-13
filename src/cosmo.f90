module cosmo
use param
implicit none
contains

!h/Mpc
SUBROUTINE init_kbin(n,karr)
	IMPLICIT NONE
	INTEGER, intent(in) :: n
	real(dp), dimension(n), intent(out) :: karr
	integer::i
	real(dp) :: kmin, delk, kmax, k

	kmin=0.01d0
	kmax=10.d0


	call logspace_dp(n,kmin,kmax,karr)

	DO i=1, n
	k=karr(i)
	if(k>0.1d0) then
	k01id=i-1
	exit
	end if
	END DO
	write(*,*) 'kmin (h/Mpc), kmax', minval(karr), maxval(karr)
!stop

END SUBROUTINE init_kbin

SUBROUTINE mean_nh_nhe(z, nh, nhe)
	IMPLICIT NONE
	real(dp)::z, nh, nhe

	nh=H_No*(1.d0+z)**3.d0
	nhe=He_No*(1.d0+z)**3.d0

	write(*,*) 'nh, nhe, z', nh, nhe, z

END SUBROUTINE mean_nh_nhe


FUNCTION tau_func1(z,x)
	implicit none
	real(dp),intent(in)::z,x
	real(dp)::tau_func1

	tau_func1=x*(1.d0+z)**2.d0/SQRT(omega_k*(1.d0+z)**2.d0+omega_m*(1.d0+z)**3.d0+&
         &omega_r*(1.d0+z)**4.d0+omega_l)

END FUNCTION



FUNCTION tau_func2(z)
	implicit none
	real(dp),intent(in)::z
	real(dp)::tau_func2

	tau_func2=(1.d0+z)**2.d0/SQRT(omega_k*(1.d0+z)**2.d0+omega_m*(1.d0+z)**3.d0+&
         &omega_r*(1.d0+z)**4.d0+omega_l)

END FUNCTION

subroutine linspace_dp(n,xmin,xmax,x)
    	implicit none
	integer, intent(in) :: n
    	real(dp),intent(in) :: xmin,xmax
    	real(dp),dimension(n),intent(out) :: x
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
    	real(dp),intent(in) :: xmin,xmax
    	real(dp), dimension(n),intent(out) :: x
   	 if (size(x) == 1 .and. xmin /= xmax) then
       		write(0,'("ERROR: Cannot call logspace with n=1 and xmin /= xmax")')
       		stop
    		end if
    	call linspace_dp(n,log10(xmin),log10(xmax),x)
    	x = 10.d0**x
end subroutine logspace_dp

!growth function at redshift z
!!
DOUBLE PRECISION FUNCTION d(z)
	IMPLICIT NONE
	DOUBLE PRECISION, INTENT(in) :: z
	DOUBLE PRECISION :: dinv

    	IF ((omega_l .EQ. 0.d0).AND.(omega_m .EQ. 1.d0)) THEN
       		dinv=1.d0+z
    	ELSEIF ((omega_l .EQ. 0.d0).AND.(omega_m .NE. 1.d0)) THEN
       		dinv=1.d0+2.5d0*omega_m*z/(1.d0+1.5d0*omega_m)
    	ELSE
       		dinv=(1.d0+((1.d0+z)**3.d0-1.d0)/&
            	&(1.d0+0.4545d0*(omega_l/omega_m)))**(1.d0/3.d0)
    	ENDIF
    	d=1.d0/dinv

END FUNCTION d



!!R in Mpc/h, M in M0/h
FUNCTION RtoM(R)
	IMPLICIT NONE
	real(dp), Intent(in) :: R
	real(dp) :: RtoM

    	RtoM=(4.d0*pi*rho_c_h)/3.d0*R**3.d0
END FUNCTION RtoM

!m in M0 , R in mpc
function MtoR(m)
	implicit none
	real(dp),intent(in)::m
	real(dp)::MtoR
	MtoR=(3.d0*m/(4.d0*PI*omega_m*RHO_c))**(1.d0/3.d0)
end function MtoR


!!M in M0/h, V in (Mpc/h)^3
!
FUNCTION MtoV(m)
	IMPLICIT NONE
	real(dp), intent(in) ::m
	real(dp) :: MtoV

	MtoV=M/rho_c_h

END FUNCTION MtoV

!!M in M0/h, V in (Mpc/h)^3
!
FUNCTION MtoRadius(M)
	IMPLICIT NONE
	real(dp), intent(in) :: M
	real(dp) :: MtoRadius,  x

    	x=3.d0*m/(4.d0*pi*rho_c_h)
    	MtoRadius=x**(1./3.)

END FUNCTION MtoRadius



!growth function at redshift z
!!
DOUBLE PRECISION FUNCTION d1(z)
	USE adaptint
	IMPLICIT NONE
	DOUBLE PRECISION, INTENT(in) :: z
	DOUBLE PRECISION :: dinv

	double precision::RELACC,ABSACC,acc, res, amax, a
	integer::MAXRUL,IPARM,ifail,N
	double precision,dimension(390)::alpha
	relacc=1.0E-9
	absacc=0.0d0
	maxrul=9
	iparm=0

	amax=1.d0/(1.d0+z)
	a=amax
	call D01ARF(0.d0,amax,dfun,RELACC,ABSACC,MAXRUL,IPARM,ACC,res,N,&
		     &ALPHA,IFAIL)

	d1=(omega_l*a**3.d0 + omega_k*a + omega_m)**0.5d0/a**1.5d0*res


END FUNCTION d1

DOUBLE PRECISION FUNCTION dfun(a)
	IMPLICIT NONE
	real(dp), intent(in) :: a
	dfun= a**1.5d0/(omega_l*a**3.d0 + omega_k*a + omega_m)**1.5d0
END FUNCTION dfun


!Recombination coefficient from Abel et al. 1997
function alpha_b(T)
	real(dp),intent(in)::T
	real(dp):: logT, ans,alpha_b
  	logT = log(T/1.1604505e4)
  	ans = Exp(-28.6130338d0 - 0.72411256*logT - 2.02604473e-2*logT**2.d0 &
	     - 2.38086188e-3*logT**3.d0 - 3.21260521e-4*logT**4.d0 &
	     - 1.42150291e-5*logT**5.d0 + 4.98910892e-6*logT**6.d0 &
	    + 5.75561414e-7*logT**7.d0 - 1.85676704e-8*logT**8.d0 &
	    - 3.07113524e-9 *logT**9.d0)
	alpha_b=ans
end function alpha_b


!gas and cmb temp decoupling redshift
function z_decoupling()
	real(dp)::z_decoupling
	z_decoupling= 137.d0*(omega_b*hlittle*hlittle/0.022d0)**0.4d0-1.d0
end function z_decoupling

!time-redhsift derivative
function dt_dz(z) 
	implicit none
	real(dp),intent(in)::z
	real(dp)::dt_dz
	dt_dz=-1.d0/(hubble_constant(z)*(1.d0+z))  
end function dt_dz



!!! return hubble costant at redshift z in sec-1
function hubble_constant(z) 
	implicit none
	real(dp),intent(in)::z
	real(dp)::hubble_constant
	hubble_constant=Ho*SQRT(omega_k*(1.d0+z)**2+omega_m*(1.d0+z)**3+&
         & omega_r*(1.d0+z)**4+omega_l)
end function hubble_constant

!cmb temp at z
function tcmb(z) 
	implicit none
	real(dp),intent(in)::z
	real(dp)::tcmb
	tcmb=cmb_temp*(1.d0+z)
end function tcmb

!cmb energy density at redshift z
FUNCTION u_gamma(z) 
	implicit none
	REAL(DP), INTENT(IN) :: z
	REAL(DP)::u_gamma
	u_gamma=4.d0*Stefan_Boltzmann/c_light*(tcmb(z))**4 
END FUNCTION u_gamma

!comoving ! in cm arXiv:0809.1326v1 
!!Distance measures in cosmology--david w hogg
function drdz(z) 
	implicit none
	real(dp),intent(in)::z
	real(dp)::drdz
	drdz=- c_light/hubble_constant(z) 
end function drdz

!gives the minimum halo mass correspond to a T_vir in M0 h-1 unit 
!https://arxiv.org/abs/astro-ph/0010468
function mass_min_cdm(temp,zcoll) !in M0 h^-1
	implicit none
	real(dp),intent(in)::temp,zcoll
	real(dp)::mass_min_cdm,hsq,hsq1
	hsq=((1.98e+4)*mui_molw*(1.d0+zcoll)/0.6d0/10.d0)**3.d0
	hsq1=omega_m*delvir(zcoll)/omega_z(zcoll)/18.d0/pi/pi
	mass_min_cdm=1.e+8*sqrt(temp**3.d0/hsq/hsq1)
end function mass_min_cdm

!!jeans mass
function jeans_mass_cdm(zcoll) !in M0 h^-1
	implicit none
	real(dp),intent(in)::zcoll
	real(dp)::jeans_mass_cdm
jeans_mass_cdm=5.73d3*(omega_m*hlittle**2./0.15)**(-0.5)*(omega_b**hlittle**2./0.022)**(-3./5.)*&
((1.+zcoll)/10.)**1.5*hlittle
end function jeans_mass_cdm

!https://arxiv.org/abs/astro-ph/0010468
FUNCTION delvir(zcoll)
	IMPLICIT NONE
	DOUBLE PRECISION :: delvir,zcoll
	DOUBLE PRECISION :: x
	x=omega_z(zcoll)-1.d0
	delvir=18.d0*pi*pi+82.d0*x-39.d0*x*x
END FUNCTION delvir

function mass_min_Lya_cool(temp,zcoll) !in M0 h^-1
	implicit none
	real(dp),intent(in)::temp,zcoll
	real(dp)::mass_min_lya_cool
	mass_min_Lya_cool=6.3e+7*(temp/1d4)**(-1.5d0)*(omega_m*hlittle**2./0.14d0)**(-0.5)* &
	&(delvir(zcoll)/18.d0/pi/pi)**(-0.5d0)*(mui_molw/0.6d0)**(-1.5d0)*((1.d0+zcoll)/13.d0)**(-1.5d0)
end function mass_min_Lya_cool

DOUBLE PRECISION FUNCTION omega_z(z)
	IMPLICIT NONE
	DOUBLE PRECISION, INTENT(in)  :: z
	omega_z=hubbledist(z)*hubbledist(z)*(1.d0+z)**3*omega_m
END FUNCTION omega_z

DOUBLE PRECISION FUNCTION hubbledist(z)
	IMPLICIT NONE
	DOUBLE PRECISION, INTENT(in) :: z
	hubbledist=1.d0/SQRT(omega_k*(1.d0+z)**2+omega_m*(1.d0+z)**3+&
		 &omega_r*(1.d0+z)**4+omega_l)

END FUNCTION hubbledist


FUNCTION n_hi(z,xion)
	implicit none
	REAL(DP), INTENT(IN) :: z,xion
	REAL(DP)::n_hi
	n_hi=(1.d0-xion)*n_b0*f_h*(1.d0+z)**3 
END FUNCTION n_hi

FUNCTION n_hei(z,xion)
	implicit none
	REAL(DP), INTENT(IN) :: z,xion
	REAL(DP)::n_hei
	n_hei=(1.d0-xion)*n_b0*f_he*(1.d0+z)**3 
END FUNCTION n_hei

FUNCTION n_heii(z,xion)
	implicit none
	REAL(DP), INTENT(IN) :: z,xion
	REAL(DP)::n_heii
	n_heii=xion*n_b0*f_he*(1.d0+z)**3
END FUNCTION n_heii


function rho_B(z)
	implicit none
	REAL(DP), INTENT(IN) :: z
	real(dp)::rho_B
	rho_B=omega_b*RHOcrit_cgs*(1.d0+z)**3
end function rho_B


!in Mpc
function comoving_distance(z1,z2)
	USE adaptint
	implicit none
	real(dp),intent(in)::z1,z2
	real(dp)::comoving_distance,tH,dH,res,res2
	double precision::RELACC,ABSACC,acc
	integer::MAXRUL,IPARM,ifail,N
	double precision,dimension(390)::alpha
	relacc=1.0E-9
	absacc=0.0d0
	maxrul=9
	iparm=0

	call D01ARF(z1,z2,drdz,RELACC,ABSACC,MAXRUL,IPARM,ACC,res2,N,&
	       &ALPHA,IFAIL)
	comoving_distance=res2/megaparsec!res*dH  !both are same
end function comoving_distance

FUNCTION sigma_HI(nu)
    USE param
    implicit none
    DOUBLE PRECISION :: sigma_HI,nu
    DOUBLE PRECISION :: sigma_0,nu_0,P,x

    sigma_0=5.475d-14
    nu_0=4.298d-1/hPlanck_eV
    P=2.963d0
    x=nu/nu_0
    sigma_HI=sigma_0*((x-1.d0)**2)*(x**(0.5*P-5.5))/(1+SQRT(x/32.88d0))**P

  END FUNCTION sigma_HI

  FUNCTION sigma_HeI(nu)
    USE param
    IMPLICIT NONE
    DOUBLE PRECISION :: sigma_HeI,nu
    DOUBLE PRECISION :: sigma_0,nu_0,P,x,yy,yw,y0,y1

    sigma_0=9.492d-16
    nu_0=1.361d1/hPlanck_eV
    P=3.188
    yw=2.039d0
    y0=0.4434d0
    y1=2.136d0

    x=nu/nu_0-y0
    yy=SQRT(x*x+y1*y1)
    sigma_HeI=sigma_0*(yw**2+(x-1.d0)**2)*(yy**(0.5*P-5.5))/(1+SQRT(yy/1.469d0))**P

  END FUNCTION sigma_HeI

  FUNCTION sigma_HeII(nu)
    USE param
    IMPLICIT NONE
    DOUBLE PRECISION :: sigma_HeII,nu
    DOUBLE PRECISION :: sigma_0,nu_0,P,x

    sigma_0=1.369d-14
    nu_0=1.72d0/hPlanck_eV
    P=2.963
    x=nu/nu_0
    sigma_HeII=sigma_0*((x-1.d0)**2)*(x**(0.5*P-5.5))/(1+SQRT(x/32.88d0))**P

  END FUNCTION sigma_HeII


function kappa_10(tki,a)
	implicit none
	real(dp),intent(in)::tki
	integer,intent(in)::a
	real(dp)::kappa_10,ans,tk
	real(dp),dimension(27)::tkin,kap
	integer::i  
	tk=tki
	 
	tkin(1) = 1.0d0
	kap(1) = 1.38e-13
	tkin(2) = 2.0d0
	kap(2) = 1.43e-13
	tkin(3) = 4.0d0
	 kap(3) = 2.71e-13
	tkin(4) = 6.0d0
	 kap(4) = 6.60e-13
	tkin(5) = 8.0d0
	 kap(5) = 1.47e-12
	tkin(6) = 10.0d0
	 kap(6) = 2.88e-12
	tkin(7) = 15.0d0
	 kap(7) = 9.10e-12
	tkin(8) = 20.0d0
	 kap(8) = 1.78e-11
	tkin(9) = 25.0d0
	 kap(9) = 2.73e-11
	tkin(10) = 30.0d0
	 kap(10) = 3.67e-11
	tkin(11) = 40.0d0
	 kap(11) = 5.38e-11
	tkin(12) = 50.0d0
	 kap(12) = 6.86e-11
	tkin(13) = 60.0d0
	 kap(13) = 8.14e-11 
	tkin(14) = 70.0d0
	kap(14) = 9.25e-11
	tkin(15) = 80.0d0
	 kap(15) = 1.02e-10
	tkin(16) = 90.0d0
	 kap(16) = 1.11e-10
	tkin(17) = 100.0d0
	 kap(17) = 1.19e-10
	tkin(18) = 200.0d0
	 kap(18) = 1.75e-10
	tkin(19) = 300.0d0
	 kap(19) = 2.09e-10
	tkin(20) = 501.0d0
	kap(20) = 2.565e-10
	tkin(21) = 701.0d0
	 kap(21) = 2.91e-10
	tkin(22) = 1000.0d0
	 kap(22) = 3.31e-10
	tkin(23) = 2000.0d0
	 kap(23) = 4.27e-10
	tkin(24) = 3000.0d0
	 kap(24) = 4.97e-10
	tkin(25) = 5000.0d0
	 kap(25) = 6.03e-10
	tkin(26) = 7000.0d0
	 kap(26) = 6.87e-10
	tkin(27) = 50000.0d0
	 kap(27) = 7.87e-10

	do i=1,27
	tkin(i)=log(tkin(i))
	kap(i)=log(kap(i))
	end do

	   


	  if (log(TK) < tkin(1)) then
	    ans = kap(1)
	  else if (log(TK) > tkin(27)) then
	    ans = log(exp(kap(27))*(TK/exp(tkin(27))**0.381d0))
	  else 
	    TK = log(TK)
	    ans = linear_interpol(TK,tkin,kap)
	 end if
	kappa_10= exp(ans)   
end function kappa_10

FUNCTION linear_interpol(x1,x,y)
	implicit none
	real(dp),intent(in)::x1
	real(dp),dimension(:),intent(in)::x,y
	integer::i,j,k
	real(dp)::a,b,linear_interpol

	k=size(x)
	do i=1,k
	a=x(i)
	 if(a >= x1) then
		 j=i
		 exit
	 end if
	end do
	if(j>1) then
		linear_interpol=(y(j)-y(j-1))/(x(j)-x(j-1))*(x1-x(j-1)) + y(j-1)
	else
		linear_interpol=y(j)-(y(j+1)-y(j))/(x(j+1)-x(j))*(x(j)-x1)
	end if
END FUNCTION linear_interpol

function kappa_10_elec(tki,a)
	implicit none
	real(dp),intent(in)::tki
	integer,intent(in)::a
	real(dp)::kappa_10_elec,kappa_10,ans,tk
	real(dp),dimension(17)::tkin,kap
	integer::i
	tk=tki

	tkin(1) = 1.0d0
	kap(1) = 0.239e-9
	tkin(2) = 2.0d0
	kap(2) = 0.337e-9
	tkin(3) = 5.0d0
	 kap(3) = 0.53e-9
	tkin(4) = 10.0d0
	 kap(4) = 0.746e-9
	tkin(5) = 20.0d0
	 kap(5) = 1.05e-9
	tkin(6) = 50.0d0
	 kap(6) = 1.63e-9
	tkin(7) = 100.0d0
	 kap(7) = 2.26e-9
	tkin(8) = 200.0d0
	 kap(8) = 3.11e-9
	tkin(9) = 500.0d0
	 kap(9) = 4.59e-9
	tkin(10) = 1000.0d0
	 kap(10) = 5.92e-9
	tkin(11) = 2000.0d0
	 kap(11) = 7.15e-9
	tkin(12) = 3000.0d0
	 kap(12) = 7.71e-9
	tkin(13) = 5000.0d0
	 kap(13) =8.17e-9
	tkin(14) = 7000.0d0
	kap(14) = 8.32e-9
	tkin(15) = 10000.0d0
	 kap(15) = 8.37e-9
	tkin(16) = 15000.0d0
	 kap(16) = 8.29e-9
	tkin(17) = 50000.0d0
	 kap(17) = 8.11e-9

	do i=1,17
	tkin(i)=log(tkin(i))
	kap(i)=log(kap(i))
	end do

	   


	  if (log(TK) < tkin(1)) then
	    ans = kap(1)
	  else if (log(TK) > tkin(17)) then
	    ans = kap(17) ! need to change
	  else 
	    TK = log(TK)
	    ans = linear_interpol(TK,tkin,kap)
	 end if

	kappa_10_elec=exp(ans) 
	!kappa_10_elec=x_e*No*((1.d0+redshift)**3)*exp(-9.607d0+0.5d0*log(tk) &
	!*exp(-(log(tk)))**(4.5)/1800)
end function kappa_10_elec

!!Fraction of X-ray energy that goes into heating..
!!
function f_heat(xe)
	IMPLICIT NONE
	real(dp),intent(in)::xe
	real(dp)::f_heat
	if(xe > 1E-4) then
	f_heat=0.9971d0*(1.0d0-(1.0d0-xe**(0.2663d0))**(1.3163d0))
	else
	f_heat=0.15d0
	end if
end function f_heat

!!Fraction of X-ray energy that goes into ionization..
!!
function f_ex(xe)
	IMPLICIT NONE
	real(dp),intent(in)::xe
	real(dp)::f_ex

	f_ex=0.4766d0*(1.0d0-xe**(0.2735d0))**(1.5221d0)
end function f_ex

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!EXTRA PRT !!!!!!!!!!!!!!!!!!!!NOT INCLUDED IN THE PAPER!!!!!!!!!!
!----------------------------------------------------------
!!NFW profile
!!r in comoving Mpc/h, M in M0/h
function nfw_profile(M, r, z)
	implicit none
	real(dp), intent(in) ::M, r, z
	real(dp) :: nfw_profile, rs, rhos

	rs=rs_char(M, z)
	rhos=rhos_char(M, z)
	nfw_profile=rhos/(r/rs)/(1.d0+r/rs)**2.d0
!write(*,*) 1.33*pi*rs**3.*rhos, m, rs, rhos, virial_radius_m(M, z)
!stop
!	Mp=M/hlittle
!	CN=concentrationn_parameter(Mp,z)
!	rvir=virial_radius_m(M, z)
!dc=178.d0/3.d0*CN**3/(log(1.d0+CN)-CN/(1.d0+CN))
!x=r/rvir
!omega_mz=omega_z(z)
!nfw_profile=rho_c/hlittle**2*(1.d0+z)**3*omega_m/omega_mz*dc/CN/x/(1.d0+CN*x)/(1.d0+CN*x)	!!M0/h/(Mpc/h)^3
end function nfw_profile

!!M in M0/h, rs_char in comoving Mpc/h
function rs_char(M, z)
implicit none
real(dp), intent(in) :: M, z
real(dp):: rs_char, Mp, CN, rvir
	Mp=M/hlittle
	CN=concentrationn_parameter(Mp,z)

	rvir=virial_radius_m(M, z)
	rs_char = rvir/CN
end function rs_char

function rhos_char(M, z)
implicit none
real(dp), intent(in) :: M, z
real(dp):: rhos_char, dc, Mp, CN, omega_mz

	Mp=M/hlittle
	CN=concentrationn_parameter(Mp,z)
	dc=178.d0/3.d0*CN**3/(log(1.d0+CN)-CN/(1.d0+CN))
	omega_mz=omega_z(z)
	rhos_char=rho_c/hlittle**2*(1.d0+z)**3*omega_m/omega_mz*dc

!
end function rhos_char


!!M in M0/h, virial_radius_m in comoving Mpc/h
function virial_radius_m(M, z)
implicit none
real(dp), intent(in) :: m, z
real(dp) :: virial_radius_m, omega_mz

omega_mz=omega_z(z)


virial_radius_m = 0.784d0*(m*1d-8)**(1.d0/3.d0)*(omega_m*178.d0/omega_mz/18.d0/pi/pi)**(-0.333d0)*(10.d0/(1.d0+z))*1d-3*(1.d0+z)	!!Mpc/h
!write(*,*) omega_mz, virial_radius_m
!stop

end function virial_radius_m

!!M in M0

function concentrationn_parameter(M,z)
implicit none
real(dp), intent(in) :: M, z
real(dp) :: concentrationn_parameter

concentrationn_parameter= (M*1d-9)**(-0.1d0)*(25.d0/(1.d0+z))

end function concentrationn_parameter

!----------------------------
!!M in M0/h, k in Mpc/h
FUNCTION u_NFW(k, m, z)
IMPLICIT NONE
real(dp), intent(in) :: k, m, z
real(dp) :: rs, rhos, Mp, CN, krs, u_NFW, SI, CI

	rs=rs_char(M, z)
	rhos=rhos_char(M, z)
	Mp=M/hlittle
	CN=concentrationn_parameter(Mp,z)
krs=k*rs
u_NFW=sin(krs)*(Si_fun((1.d0+CN)*krs) -Si_fun(krs)) -sin(CN*krs)/(1.d0+CN)/krs &
& + cos(krs)*(Ci_fun((1.d0+CN)*krs) -Ci_fun(krs))

!write(*,*) u_NFW, 4.d0*pi*rhos*rs**3., m
!stop

u_NFW = u_NFW*4.d0*pi*rhos*rs**3./m
END FUNCTION u_NFW



function Si_fun(x)
USE adaptint
implicit none
real(dp), intent(in) :: x
real(dp) :: Si_fun

	double precision::RELACC,ABSACC,acc, res, amax, a
	integer::MAXRUL,IPARM,ifail,N
	double precision,dimension(390)::alpha
	relacc=1.0E-9
	absacc=0.0d0
	maxrul=9
	iparm=0


	a=amax
	call D01ARF(0.d0,x,Si_fun_integral,RELACC,ABSACC,MAXRUL,IPARM,ACC,res,N,&
		     &ALPHA,IFAIL)
	Si_fun=res

end function Si_fun

function Ci_fun(x)
USE adaptint
implicit none
real(dp), intent(in) :: x
real(dp) :: Ci_fun

	double precision::RELACC,ABSACC,acc, res, amax, a
	integer::MAXRUL,IPARM,ifail,N
	double precision,dimension(390)::alpha
	relacc=1.0E-9
	absacc=0.0d0
	maxrul=9
	iparm=0


	call D01ARF(0.d0,x,Ci_fun_integral,RELACC,ABSACC,MAXRUL,IPARM,ACC,res,N,&
		     &ALPHA,IFAIL)

	Ci_fun=0.577d0 + log(x) + res	!http://mathworld.wolfram.com/CosineIntegral.html

end function Ci_fun


function Si_fun_integral(x)
implicit none
real(dp), intent(in) :: x
real(dp) :: Si_fun_integral

Si_fun_integral=sin(x)/x

end function Si_fun_integral

function Ci_fun_integral(x)
implicit none
real(dp), intent(in) :: x
real(dp) :: Ci_fun_integral

	Ci_fun_integral=(cos(x)-1.d0)/x

end function Ci_fun_integral




end module cosmo
