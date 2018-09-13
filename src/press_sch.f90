module press_sch
use cosmo

DOUBLE PRECISION :: rcomm

contains


!!RMS of mass variance, M in M0/h.. at z=0.. This is due to the fact del_c(z) is used in fcoll calculation and sigma(m) at present epoch
!!!checked with online tool HMFcalc
DOUBLE PRECISION FUNCTION sigma_m_CDM(M)
IMPLICIT NONE
DOUBLE PRECISION, INTENT(in) :: M
DOUBLE PRECISION :: X
	x=MtoRadius(M)
    	sigma_m_CDM=sigma_CDM(x)
END FUNCTION sigma_m_CDM


!!RMS of mass variance, X in Mpc/h.. at z=0.. This is due to the fact del_c(z) is used in fcoll calculation and sigma(m) at present epoch
!!!
DOUBLE PRECISION FUNCTION sigma_CDM(x)
IMPLICIT NONE
DOUBLE PRECISION, INTENT(in) :: x
    sigma_CDM=SQRT(sigmasq_CDM(x))
END FUNCTION sigma_CDM

!normalised mass variance at z=0.. x in Mpc/h
!!!This works both for CDM and FDM
DOUBLE PRECISION FUNCTION sigmasq_CDM(x)
USE adaptint, ONLY : d01amf
IMPLICIT NONE
    DOUBLE PRECISION, INTENT(in) :: x
    INTEGER :: ifail,inf
    INTEGER, PARAMETER :: lw=2000,liw=lw/4
    DOUBLE PRECISION :: epsabs,epsrel,res,abserr
    DOUBLE PRECISION, DIMENSION(lw) :: w
    INTEGER, DIMENSION(liw) :: iw

    rcomm=x
    inf=1
    epsabs=0.
    epsrel=1.0d-3
    ifail=-1
    CALL d01amf(integrand_sigmasq_CDM,0.d0,inf,epsabs,epsrel,res,abserr,w,&
         &lw,iw,liw,ifail)
    sigmasq_CDM=res*norm_CDM(sigma_8)/(2.d0*pi*pi) 

END FUNCTION sigmasq_CDM

!!not normalized..r in in Mpc/h.. k is h/Mpc
!
DOUBLE PRECISION FUNCTION integrand_sigmasq_CDM(k)
	IMPLICIT NONE
	DOUBLE PRECISION, INTENT(in) :: k
	DOUBLE PRECISION::r
	r=rcomm
    	IF (r .LT. 1d-7) THEN
       		integrand_sigmasq_CDM=(1.d0-(k*k*r*r/5.d0))*k*k*pspec_cdm(k)
    	ELSE
       		integrand_sigmasq_CDM=k*k*pspec_CDM(k)*window(k*r)
    	ENDIF

END FUNCTION integrand_sigmasq_CDM

!top hat window function
!
FUNCTION window(x)
	IMPLICIT NONE
	REAL(DP), INTENT(in) :: x
	REAL(DP) :: window

    	window=3.d0*((SIN(x)/x**3.d0)-(COS(x)/x**2.d0))
    	window=window*window

END FUNCTION window

!matter power spectrum padmanabhan structure formation page320..not normalized.. at z=0.. k in h/Mpc
! checked with online plots
DOUBLE PRECISION FUNCTION pspec_CDM(k)
	IMPLICIT NONE
	DOUBLE PRECISION, INTENT(in) :: k
	DOUBLE PRECISION :: gam_h,pspecdm3d,n0,gama

  	gama=omega_m*hlittle*exp(-omega_b-sqrt(2.d0*hlittle)*omega_b/omega_m)
  	gam_h=gama*hlittle

  	IF (ABS(k) < 1.d-40) THEN
     		pspecdm3d=0.d0
  	ELSE
     		n0=ns+0.5d0*dn_dlnk*LOG(k/k0)
     		pspecdm3d=(k**n0)*transfun_CDM(k)*transfun_CDM(k)
  	END IF

	  	pspec_CDM=pspecdm3d

END FUNCTION pspec_CDM



!transfer function.. k in h/mpc..https://arxiv.org/pdf/astro-ph/9709112.pdf
!!
FUNCTION transfun_CDM(k)
	IMPLICIT NONE
	REAL(DP), INTENT(in) :: k
	REAL(DP) :: transfun_cdm, q, s, alpha_gamma, gamma_eff, theta, L0, C0


	s=44.5d0*log(9.83d0/omega_m/hlittle/hlittle)/sqrt(1.d0+10.d0*(omega_b*hlittle*hlittle)**(3.d0/4.d0))	!!in Mpc
	s=s*hlittle !Mpc/h
	alpha_gamma=1.d0-0.328d0*log(431.d0*omega_m*hlittle*hlittle)*omega_b/omega_m &
	& +0.38d0*log(22.3*omega_m*hlittle*hlittle)*(omega_b/omega_m)**2.d0
	gamma_eff=omega_m*hlittle*(alpha_gamma + (1.d0-alpha_gamma)/(1.d0+(0.43d0*k*s)**4.d0))
	theta=cmb_temp/2.7d0
	q=k*theta**2.d0/gamma_eff
	L0=log(2.d0*2.71828d0 + 1.8d0*q)
	C0= 14.2d0 + 731.d0/(1.d0+62.5d0*q)

  	transfun_cdm=L0/(L0+C0*q*q)



END FUNCTION transfun_CDM


!*******************************NORMALISATION******************************  
!!normalisation to sigma and matter power spectrum using..sigma 8
!
FUNCTION norm_CDM(sig)
	USE adaptint, ONLY : d01amf
	IMPLICIT NONE
    	REAL(DP), INTENT(in) :: sig
    	REAL(DP) :: norm_CDM
    	INTEGER, PARAMETER :: lw=2000,liw=lw/4
    	INTEGER :: inf,ifail
    	REAL(DP) :: epsabs,epsrel,res,abserr
    	INTEGER, DIMENSION(liw) :: iw
    	REAL(DP), DIMENSION(lw) :: w

    	inf=1
    	epsabs=0.d0
    	epsrel=1.0d-4
    	ifail=-1
    	CALL d01amf(pswin_CDM,0.d0,inf,epsabs,epsrel,res,abserr,w,lw,iw,liw,ifail)

    	norm_CDM=sig*sig*2.d0*pi*pi/res

END FUNCTION norm_CDM

FUNCTION pswin_CDM(k)
	IMPLICIT NONE
    	REAL(DP), INTENT(in) :: k
    	REAL(DP) :: pswin_CDM
    	REAL(DP) :: r1

    	r1=8.d0	!in Mpc/h	!!This is for intergraion k space
    	pswin_CDM=pspec_CDM(k)*k*k*window(k*r1)

END FUNCTION pswin_CDM  

!**********************************normalization ends**************************************
 
!!dndm for cdm..M is M0/h.. number per M0/h per (Mpc/h)^3
!! tested with HMFcalc
DOUBLE PRECISION FUNCTION numdenm_cdm(m,z)
    	IMPLICIT NONE
    	DOUBLE PRECISION, INTENT(in) :: m,z
    	DOUBLE PRECISION :: x


	x=MtoRadius(M)	!Mpc/h
    	numdenm_cdm=numdenx_cdm(x,z)/(4.d0*pi*x*x*rho_c_h)

END FUNCTION numdenm_cdm

DOUBLE PRECISION FUNCTION numdenx_cdm(x,z)
    	IMPLICIT NONE
    	DOUBLE PRECISION, INTENT(in) :: x,z
    	DOUBLE PRECISION :: sig,dlogsigmadlogx,step,err,nu

    	sig=sigma_cdm(x)
    	step=0.1d0
    	DO
       		dlogsigmadlogx=dfridr(logsigma_cdm,LOG(x),step,err)
       		IF (ABS(err).LT.1.d-5) exit
       		step=step/2.d0
    	END DO
    	nu=delta_c/(d(z)*sig)
    	numdenx_cdm=-probdist(nu)*(3.d0/(4.d0*pi*x**3))*dlogsigmadlogx*nu/x

END FUNCTION numdenx_cdm


DOUBLE PRECISION FUNCTION dfridr(func,x,h,err)
    IMPLICIT NONE

    	DOUBLE PRECISION, INTENT(in) :: x,h
    	DOUBLE PRECISION, INTENT(out) :: err
    	INTEGER, parameter :: NTAB=10
    	DOUBLE PRECISION, PARAMETER :: CON=1.4d0,CON2=CON*CON,BIG=1.d30,SAFE=2.d0
    	INTEGER :: ierrmin,i,j
    	INTEGER, DIMENSION(1) :: imin
    	DOUBLE PRECISION :: hh
    	DOUBLE PRECISION, DIMENSION(NTAB-1) :: errt,fac
    	DOUBLE PRECISION, DIMENSION(NTAB,NTAB) :: a
    	INTERFACE
       		DOUBLE PRECISION FUNCTION func(x)
         	DOUBLE PRECISION, intent(in) :: x
       		END FUNCTION func
    	END INTERFACE

    IF(h.EQ.0.) STOP 'h must be nonzero in dfridr'
    hh=h
    a(1,1)=(func(x+hh)-func(x-hh))/(2.d0*hh)
    err=BIG
    fac(1)=CON
    DO i=2,NTAB-1
       fac(i)=fac(i-1)*CON
    ENDDO
    DO i=2,NTAB
       hh=hh/CON
       a(1,i)=(func(x+hh)-func(x-hh))/(2.d0*hh)
       DO j=2,i
          a(j,i)=(a(j-1,i)*fac(j-1)-a(j-1,i-1))/(fac(j-1)-1.d0)
       ENDDO
       errt(1:i-1)=MAX(ABS(a(2:i,i)-a(1:i-1,i)),ABS(a(2:i,i)-a(1:i-1,i-1)))
       imin=MINLOC(errt(1:i-1))
       ierrmin=imin(1)
       IF (errt(ierrmin).LE.err) THEN
          err=errt(ierrmin)
          dfridr=a(1+ierrmin,i)
       ENDIF
       IF(ABS(a(i,i)-a(i-1,i-1)).GE.SAFE*err) RETURN
    ENDDO
END FUNCTION dfridr

DOUBLE PRECISION FUNCTION probdist(nu)
    	IMPLICIT NONE
    	DOUBLE PRECISION, INTENT(in) :: nu
    	DOUBLE PRECISION :: k

	DOUBLE PRECISION, parameter :: Ac=0.3222d0
	DOUBLE PRECISION, parameter :: a=0.707d0
	DOUBLE PRECISION, parameter :: p = 0.3d0

    	k=SQRT(2.d0/pi)
	if(sheth_tormen==1) then
    		probdist=Ac*sqrt(a)*k*(1.d0  + (1.d0/a/nu/nu)**p)*EXP(-nu*nu*a/2.d0)
	elseif(sheth_tormen==0) then
    		probdist=k*EXP(-nu*nu/2.d0)
	else
		stop
	end if

END FUNCTION probdist

DOUBLE PRECISION FUNCTION logsigma_cdm(logx)
    IMPLICIT NONE
    DOUBLE PRECISION, INTENT(in) :: logx
    DOUBLE PRECISION :: x

    x=EXP(logx)
    logsigma_cdm=LOG(sigma_cdm(x))

END FUNCTION logsigma_cdm


!!!!!!!!!!!!!


!!This is using the dndm relation..
!!M in M0/h
double precision function fcoll_cdm(z)
	implicit none
	real(dp), intent(in)::z
	real(dp):: M1, M2, delM, res, M, itg
	real(dp), parameter::increase_f=1.2d0


	M=mass_min_cdm(t_vir_param,z) !in M0/h

	M1=M
	res=0.d0

	DO
		M2=M1*increase_f
		M=(M1+M2)*0.5d0
		delM=M2-M1
		itg = int_fcoll_cdm(M,z)
		res=res+itg*delM
		if(M>1d13 .or. itg<1d-10) exit
		
		M1=M2
	END DO

!!res is amount os mass in M0/h in (Mpc/h)^3 volume

	fcoll_cdm=res/rho_c_h

end function fcoll_cdm



!!M in M0/h.. 
double precision function int_fcoll_cdm(M, z)
	implicit none
	double precision, intent(in) ::z, M

	int_fcoll_cdm=numdenm_cdm(M,z)*M

end function int_fcoll_cdm

end module press_sch
