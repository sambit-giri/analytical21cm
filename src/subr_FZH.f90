module subr_FZH
	use nrtype
	use subr					
	use funct
	use param
	use adaptint
	use cosmo

	DOUBLE PRECISION :: r_correl_dd

contains




subroutine init_xhii_inhomogeneity(z)
implicit none
	real(dp), intent(in) :: z
	real(dp) ::  dc, gf, M, sig
	integer :: i
    	DOUBLE PRECISION :: x, sfac, Tmz, S, Sm, V, haloB, r, cor, B, numdenm_hii
    	DOUBLE PRECISION :: T1, T2, T3, T4, T5

dndlnm_arr_comm=0.d0
cor_dd_arr_comm=0.d0
dx_arr_comm=0.d0
halo_bias_arr_comm=0.d0

	gf=d(z)
	dc=delta_c/gf


!!dd correlation
	DO i=1, nrdim_comm
		r=r_arr_comm(i)
		cor=correlation_func_dd(r, z)
		cor_dd_arr_comm(i)=cor	
	END DO



	mminz_comm=mass_min_cdm(t_vir_param,z) !in M0/h
    	sigmamin_comm=sigma_m_cdm(mminz_comm)


	B0z_comm=dc-sqrt(2.d0)*k_zeta_comm*sigmamin_comm
	if(B0z_comm<0.d0) B0z_comm=0.d0
	B1z_comm=k_zeta_comm/sqrt(2.d0)/sigmamin_comm


!open(unit=57, file='test.dat')

!!sigma array
DO i=1, massdim_comm
	M=M_arr_comm(i)
	sig=sigma_arr_comm(i)
	B=B0z_comm+ B1z_comm*sig**2.d0
!if(B<0.d0) B=0.d0

if(sigmamin_comm>sig) then
	dx_arr_comm(i) = dc - sqrt(2.d0)*k_zeta_comm*sqrt(sigmamin_comm**2.d0 - sig**2.d0)
else
	dx_arr_comm(i)=0.d0
end if
!!

	haloB=1.d0+ (B/sig/sig -1.d0/B0z_comm)/gf	!!ZFH05
	if(haloB<0.d0) haloB=0.d0
	halo_bias_arr_comm(i) = haloB

!!
	V=MtoV(M)
	V_arr_comm(i) = V

 

!!
	x=MtoRadius(M)
	R_m_arr_comm(i) = x
	numdenm_hii=numdenx_hii_fun(x,z, B, B0z_comm)/(4.d0*pi*x*x*rho_c_h)



if(sig>=sigmamin_comm) then
	dndlnm_arr_comm(i)=0.d0
!write(*,*) sigmamin_comm, sig, 'ra'
else
!if(sig>sigmamin_comm*1d-2) then
	S = sig**2.d0
	Sm = sigmamin_comm**2.d0
	sfac = Sm - S
	T1=-S/2.d0/sqrt(sfac)
	T2=S*S/2.d0/4.d0/sfac**(1.5d0)
	T3=-S**3./6.d0/8.d0*3.d0/sfac**(2.5d0)
	T4=S**4./24.d0/16.d0*15.d0/sfac**(3.5d0)
	T5=- S**5./120.d0/32.d0*105.d0/sfac**(4.5d0)
	Tmz=B + sqrt(2.d0)*k_zeta_comm*(T1+T2+T3+T4+T5) 

	if(Tmz<0.d0) Tmz=0.d0


!write(*,*) 'raghu',M,tmz, B ,sqrt(2.d0)*k_zeta_comm*(T1+T2+T3+T4+T5) , s, sm, B1z_comm


	dndlnm_arr_comm(i) = numdenm_hii*Tmz


end if


!if(z==8.d0) then
!write(*,*) M, B0z_comm, B1z_comm, numdenm_hii, Tmz, sigmamin_comm, sig
!end if

!write(57,'(17e20.4)') x, numdenm_hii, tmz, dndlnm_arr_comm(i), test_comm, numdenx_hii_fun(x,z, B, B0z_comm), 2.71828d0**(-0.5d0*test_comm*test_comm)

END DO
!	close(57)

!write(*,*) dndlnm_arr_comm, mminz_comm
!stop


!if(z==8.d0) stop

end subroutine init_xhii_inhomogeneity



!!this scheme will not work as at late stages, the bubble distribution will ve very narrow. potentially miss those distribution

!!Volume fraction ..This is using the analytical formula Mcquin05
!!
function vol_frac_Q(z, ion_ef)
	implicit none
	real(dp), intent(in)::z, ion_ef
	real(dp)::rho_0, M1, M2, delM, res, M, dndm, V, B, B0, B1, sigmasq, delxmz, sigmin, sig, mmin, resp, vol_frac_Q
	real(dp), parameter::increase_f=2.d0
	integer :: i


	Mmin=mminz_comm


	M1=Mmin
	res=0.0

	DO i=1, massdim_comm-1
	M1=M_arr_comm(i)
	if(M1>Mmin) then
		M2=M_arr_comm(i+1)
		delM=M2-M1
		M=M1
		V=V_arr_comm(i)
		dndm=dndlnm_arr_comm(i)

if(dndm>0.d0)		res=res+dndm*delM*V


	end  if

	END DO


	vol_frac_Q=res
!if(vol_frac_Q>=1.d0) vol_frac_Q=0.9999d0

!if(sum(dx_arr_comm) .ne. 0.d0 .and. sum(dndlnm_arr_comm)==0.d0) vol_frac_Q=0.9999d0

end function vol_frac_Q



!!
DOUBLE PRECISION FUNCTION numdenx_hii_fun(x,z, B, B0)
    	IMPLICIT NONE
    	DOUBLE PRECISION, INTENT(in) :: x,z
    	DOUBLE PRECISION :: sig,dlogsigmadlogx,step,err,nu, B, B0, nu1, Tmz


    	sig=sigma_cdm(x)
    	step=0.1d0
    	DO
       		dlogsigmadlogx=dfridr(logsigma_cdm,LOG(x),step,err)
       		IF (ABS(err).LT.1.d-5) exit
       		step=step/2.d0
    	END DO


    	nu=B/sig
	nu1=1.d0/sig	!!changes 
if(nu<15.d0 .and. nu>0.d0) then
    	numdenx_hii_fun=-probdist(nu)*(3.d0/(4.d0*pi*x**3.d0))*dlogsigmadlogx*nu1/x
else
    	numdenx_hii_fun=0.d0
end if


END FUNCTION numdenx_hii_fun





!!!!!!!! Correlation part..!!!!


!!overlapped volume in FZH04
!!M in M0/h and r12 in Mpc/h
!!
FUNCTION Vop_fn(R, r12)
	IMPLICIT NONE
	real(dp), intent(in) :: R, r12
	real(dp):: vop_fn, x, rho_0, h



	if(r12>2.d0*R) then
		Vop_fn=0.d0
	elseif(r12<2.d0*R) then
		Vop_fn = pi*r12*(R**2.d0 - r12**2.d0/12.d0) 
	else
		Vop_fn=4.d0*pi*R**3.d0/3.d0
	end if

END FUNCTION Vop_fn


!!volume encompass in mcquinn
!!M in M0/h and r12 in Mpc/h
!!
FUNCTION Vo_fn(R, r12)
	IMPLICIT NONE
	real(dp), intent(in) :: R, r12
	real(dp):: vo_fn,  x, rho_0, h




	if(r12>=2.d0*R) then
		Vo_fn=0.d0
	else
		Vo_fn =4.d0*pi*R**3.d0/3.d0-4.d0*pi*r12**3.d0/3.d0/8.d0

	end if

END FUNCTION Vo_fn


!!r1 r2 are the radius of the spheres and d is the distance between the centers..
!will return overlap volume
!!
FUNCTION overlap_sphere(r1, r2, d)
	IMPLICIT NONE
	real(dp), intent(in) :: r1, r2, d
	real(dp) :: overlap_sphere

	if(d>(r1+r2) .or. r1<=0.0 .or. r2<=0.0) then
		overlap_sphere=0.d0
	else
		overlap_sphere = pi*(r1+r2-d)**2.d0*(d**2.d0+2.d0*d*(r1+r2)-3.d0*(r1**2.0+r2**2.0)+6.d0*r1*r2)/12.d0/d
	end if

	if(overlap_sphere <0.0) then
		write(*,*) 'not a realistic scenario to calculate overlap'
		write(*,*) overlap_sphere, r1, r2, d
		stop
	end if

END FUNCTION overlap_sphere

!!This is using the dndm_hii relation..
!!M in M0/h
!!
SUBROUTINE Initialize_correlation_dd_xx_xd(z, ion_ef, Q, nrdim, rarr, xi_dd, xi_xx, xi_xd)
	implicit none
	real(dp), intent(in)::z, ion_ef, Q
	INTEGER, INTENT(in) :: nrdim
	real(dp), dimension(nrdim), intent(in) :: rarr
	real(dp), dimension(nrdim), intent(out) :: xi_dd, xi_xx, xi_xd
	real(dp), dimension(nrdim) :: resp, vop, P1, P2,  Pout, Pin
	real(dp)::M1, M2, delM, res, M, dndm, V, Mmin,   av_halo_bias, R1, r2, bias, r12, r
	real(dp) ::  m1p, m2p, delmp, mp, dndmp, vp, co_dd, sigmin, sig, delxmz, b, b0, b1, bxsq,  over, rmin, gfactor
	real(dp), parameter::increase_f=2.d0
	integer :: i, ii, jj

!	call init_arr_corel_dd(z, nrdim, rarr, xi_dd)
gfactor=d(z)
xi_dd = cor_dd_arr_comm

!	r12=rin

	av_halo_bias=av_halo_bias_B(z, ion_ef)
	bxsq=av_halo_bias**2.d0

	P1=0.d0
	P2=0.d0
	Pin=0.d0
	Pout=0.d0


	Mmin=mminz_comm


	M1=Mmin
	res=0.d0

	DO ii=1, massdim_comm-1

M1=M_arr_comm(ii)
if(M1>=Mmin) then
		M2=M_arr_comm(ii+1)
		delM=M2-M1
		M=M1
		dndm=dndlnm_arr_comm(ii)
		V=V_arr_comm(ii)
		R1=R_m_arr_comm(ii) 

		DO i=1, nrdim
			r12=rarr(i)
			Vop(i)=Vop_fn(R1, r12)	
			P1(i)=P1(i)+(dndm)*delM*Vop(i)

			if(r12>R1) then	!!rh is outside the bubble at Rb
				r=max(r12, R1)
				if(r==r12) then
					co_dd=xi_dd(i)
				else
					co_dd=linear_interpol(r,rarr,xi_dd)
				end if
				Pout(i)=Pout(i)+ delm*dndm*V*co_dd*av_halo_bias
			end if
		END DO





!!This is for the two point term.. xx



		resp=0.d0
		M1p=Mmin
		DO jj=1,  massdim_comm-1
M1p=M_arr_comm(jj)
if(M1p>=Mmin) then
			M2p=M_arr_comm(jj+1)
			delMp=M2p-M1p
			Mp=M1p
			dndmp=dndlnm_arr_comm(jj)
			Vp=V_arr_comm(jj)
			R2=R_m_arr_comm(jj) 

		DO i=1, nrdim
			r12=rarr(i)
			r=max(r12, (R1+R2))

			if(r==r12) then
				co_dd=xi_dd(i)
			else
				co_dd=linear_interpol(r,rarr,xi_dd) 
			end if
			if(r12>=(R1+r2)) then
				resp(i)=resp(i)+ delmp*dndmp*(1.d0+bxsq*co_dd)*V*Vp
			else
				over=overlap_sphere(r1, r2, max(r1,r2))
				resp(i)=resp(i)+ delmp*dndmp*(1.d0+bxsq*co_dd)*(V-over)*(Vp-over)
			end if
		END DO

end if
		END DO

		P2=P2+resp*delm*dndm

!!P2 ends..and Pin starts..



		delxmz=dx_arr_comm(ii)
		Pin=Pin+(dndm)*delM*Vop*(1.d0+delxmz*gfactor)		!!d(z) as different convension in Mcquinn than Furlanett
!if(z==8.d0) then
! write(*,*) pin, Vop
!write(*,*) dndm, delM,  delxmz, gfactor
!!stop
!end if

end if
	END DO

!!calculating correlation functions..

!!xi_xx..

	if(Q>ion_cut) then
		xi_xx = (1.d0-Q)*P1 + Q**2.d0
	else
		xi_xx = P1 + P2
	end if
	xi_xx = xi_xx -Q**2.d0
	!if(xi_xx<0.0) xi_xx=0.0

!xi_xd..

	if(Q>ion_cut) then
		xi_xd = p1-pin
	else
		xi_xd = p1-pin-pout
	end if



!if(z==8.d0) then
!write(*,*) pin, dx_arr_comm!Q, z, P1, P2, pin, pout,  xi_dd, xi_xx, xi_xd
!stop
!end if

END SUBROUTINE Initialize_correlation_dd_xx_xd

!!Mcquinn..checked values..
!!
FUNCTION av_halo_bias_B(z, ion_ef)
	IMPLICIT NONE
	real(dp), intent(in) :: z, ion_ef
	real(dp) ::  av_halo_bias_B,  B, B0, B1, sigsq
	integer :: i
	real(dp):: M1, M2, delM, res, M, dndm, V, res1, bias, Mmin, Q
	real(dp), parameter::increase_f=2.d0


	Mmin=mminz_comm !in M0/h
	res=0.d0
	res1=0.d0


	DO i=1, massdim_comm-1
		M1=M_arr_comm(i)

if(M1>=Mmin) then
		M2=M_arr_comm(i+1)
		delM=M2-M1
		M=M1
		dndm=dndlnm_arr_comm(i)
		V=V_arr_comm(i)
		res=res+(dndm)*delM*V
		bias=halo_bias_arr_comm(i)
		res1=res1+(dndm)*delM*V*bias
end if
	END DO


	av_halo_bias_B=res1/res

END FUNCTION av_halo_bias_B





!!!two point correlation of dark matter density..
! r12 in Mpc/h
!http://www.astro.caltech.edu/~george/ay21/eaa/eaa-powspec.pdf
!
FUNCTION correlation_func_dd(r12, z)
	IMPLICIT NONE
	real(dp), intent(in) :: r12, z
	real(dp) ::  correlation_func_dd, Q, Qo, x1x2, correc, xi, gf, norm

    	INTEGER :: ifail,inf
    	INTEGER, PARAMETER :: lw=2000,liw=lw/4
    	DOUBLE PRECISION :: epsabs,epsrel,res,abserr
   	 DOUBLE PRECISION, DIMENSION(lw) :: w
    	INTEGER, DIMENSION(liw) :: iw

    	r_correl_dd=r12
    	inf=1
   	epsabs=0.
   	epsrel=1.0d-3
    	ifail=-1
    	CALL d01amf(integrand_correl_dd,0.d0,inf,epsabs,epsrel,res,abserr,w,&
         &lw,iw,liw,ifail)

	gf=d(z)	!!as pspec(k) corresponds to z=0 
	norm=norm_cdm(sigma_8)


	correlation_func_dd = res*norm/2.d0/pi**2.d0*gf**2.d0

	if(correlation_func_dd<0.0) correlation_func_dd=0.0


END FUNCTION correlation_func_dd



!!
DOUBLE PRECISION FUNCTION integrand_correl_dd(k)
	IMPLICIT NONE
	DOUBLE PRECISION, INTENT(in) :: k
	DOUBLE PRECISION :: r

	r = r_correl_dd
	integrand_correl_dd=k*k*pspec_cdm(k)*window1(k*r)

END FUNCTION integrand_correl_dd

!!
FUNCTION window1(x)
	IMPLICIT NONE
	REAL(DP), INTENT(in) :: x
	REAL(DP) :: window1

    	IF (x .LT. 1d-7) THEN
    		window1=1.d0-(x*x/6.d0)
	ELSEIF(x>1d3) THEN
    		window1=0.d0
	ELSE
    		window1=SIN(x)/x
	END IF


END FUNCTION window1

!!This need to be improved.. results significant ringing..
!! k in h/Mpc
!!
SUBROUTINE PS_TB_ion(n,karr, z, ion_ef, Ptot, Q)
	IMPLICIT NONE
	integer, intent(in) ::n
	real(dp), dimension(n), intent(in) ::karr
	real(DP), intent(in) :: z, ion_ef
	real(dp), dimension(n), intent(out)::Ptot
	real(DP), intent(out) :: Q
	real(dp) :: r1, r2 , cxx, cdd, cdx, res,  ctot, delr,  xhi, r,  avtb, aa
	real(dp) ::  delr_max,  rr, dr, drtot, cdd1, cxx1, cdx1, k, win, rmin, rmax
	integer::cn, i, ii
	real(dp), dimension(n) :: winarr, coefarr, bbarr

	integer, parameter :: nrdim = nrdim_comm
	real(dp), dimension(nrdim) ::rarr, xi_dd, xi_xx, xi_xd

	Q=vol_frac_Q(z, ion_ef)


	!call logspace_dp(nrdim,rmin,rmax,rarr)	!!Mpc/h
	rarr = r_arr_comm
	call Initialize_correlation_dd_xx_xd(z, ion_ef, Q, nrdim, rarr, xi_dd, xi_xx, xi_xd)


	xhi=1.0-Q
	res=0.d0

	ptot=0.0


	DO ii=1,nrdim-1
		r1=rarr(ii)
		r2=rarr(ii+1)
		delr=r2-r1
		r=(r1+r2)*0.5d0

		winarr=0.
		DO i=1, n
			k=karr(i)
			bbarr(i)=k**3.d0/2.d0/pi/pi
			delr_max=0.1/k

		if(k*r>5000.)  then
			win=0.0
		else

			if(delr>delr_max) then
				dr=delr_max
				rr=r1
				!if(k*rr>100.) exit

				win=0.0
				drtot=0.0
				DO
					win=win+window1(k*rr)*dr
					drtot=drtot+dr
					rr=rr+dr
					if(rr>r2) exit
				end do
				win=win/drtot
			else	
				win=window1(k*r)
			end if
		end if
			winarr(i) =win
		end do


		if(r>1000.d0 .or. (k*r >10000.)) exit
		aa=r**2.d0*delr*4.d0*pi
		coefarr=winarr*aa
		if(coefarr(1) .ne. 0.0) then
			cdd=xi_dd(ii)
			cxx=xi_xx(ii)
			cdx=xi_xd(ii)
			ctot=xhi**2.d0*cdd + xhi*(1.d0-xhi)*(cdx+cxx) 	!!zaldarriaga 2004


			Ptot=Ptot+ctot*coefarr


		end if

		r1=r2

	END DO

DO i=1, nrdim
if(ptot(i)<0.d0) ptot(i) = 0.d0
END DO 

!write(*,*) Ptot(1), ctot, coefarr(1), Ptot(k01id)

!!dimensionless power spectrum..


	Ptot=Ptot*bbarr 


	!Ptot=sqrt(Ptot)

!write(*,*) Ptot(1), ctot, coefarr(1), Ptot(k01id), bbarr(k01id)
!stop
END SUBROUTINE PS_TB_ion




SUBROUTINE cal_PS_ion(ion_ef, z, ndim, karr, Ptot, Q)
	IMPLICIT NONE
	real(dp), intent(in) ::ion_ef, z
	integer, intent(in) :: ndim
	real(dp), dimension(ndim), intent(in)::karr
	real(dp), dimension(ndim), intent(out):: Ptot
	real(dp), intent(out) :: Q
	real(dp), dimension(ndim):: Pdd, Pxx, Pdx,  psdm

	real(dp) ::  k,  cxx, avtb

	integer ::i


	call  PS_TB_ion(ndim,karr, z, ion_ef, Ptot, Q)


END SUBROUTINE cal_PS_ion





end module subr_FZH


