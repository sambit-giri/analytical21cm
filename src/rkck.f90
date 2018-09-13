	SUBROUTINE rkck(y,dydx,x,h,yout,yerr,derivs)
	USE nrtype; USE nrutil, ONLY : assert_eq
	IMPLICIT NONE
	REAL(dp), DIMENSION(:), INTENT(IN) :: y,dydx
	REAL(dp), INTENT(IN) :: x,h
	REAL(dp), DIMENSION(:), INTENT(OUT) :: yout,yerr
	INTERFACE
		SUBROUTINE derivs(x,y,dydx)
		USE nrtype
		IMPLICIT NONE
		REAL(dp), INTENT(IN) :: x
		REAL(dp), DIMENSION(:), INTENT(IN) :: y
		REAL(dp), DIMENSION(:), INTENT(OUT) :: dydx
		END SUBROUTINE derivs
	END INTERFACE
	INTEGER(I4B) :: ndum
	REAL(dp), DIMENSION(size(y)) :: ak2,ak3,ak4,ak5,ak6,ytemp
	REAL(dp), PARAMETER :: A2=0.2d0,A3=0.3d0,A4=0.6d0,A5=1.0d0,&
		A6=0.875d0,B21=0.2d0,B31=3.0d0/40.0d0,B32=9.0d0/40.0d0,&
		B41=0.3d0,B42=-0.9d0,B43=1.2d0,B51=-11.0d0/54.0d0,&
		B52=2.5d0,B53=-70.0d0/27.0d0,B54=35.0d0/27.0d0,&
		B61=1631.0d0/55296.0d0,B62=175.0d0/512.0d0,&
		B63=575.0d0/13824.0d0,B64=44275.0d0/110592.0d0,&
		B65=253.0d0/4096.0d0,C1=37.0d0/378.0d0,&
		C3=250.0d0/621.0d0,C4=125.0d0/594.0d0,&
		C6=512.0d0/1771.0d0,DC1=C1-2825.0d0/27648.0d0,&
		DC3=C3-18575.0d0/48384.0d0,DC4=C4-13525.0d0/55296.0d0,&
		DC5=-277.0d0/14336.0d0,DC6=C6-0.25d0
	ndum=assert_eq(size(y),size(dydx),size(yout),size(yerr),'rkck')
	ytemp=y+B21*h*dydx
	call derivs(x+A2*h,ytemp,ak2)
	ytemp=y+h*(B31*dydx+B32*ak2)
	call derivs(x+A3*h,ytemp,ak3)
	ytemp=y+h*(B41*dydx+B42*ak2+B43*ak3)
	call derivs(x+A4*h,ytemp,ak4)
	ytemp=y+h*(B51*dydx+B52*ak2+B53*ak3+B54*ak4)
	call derivs(x+A5*h,ytemp,ak5)
	ytemp=y+h*(B61*dydx+B62*ak2+B63*ak3+B64*ak4+B65*ak5)
	call derivs(x+A6*h,ytemp,ak6)
	yout=y+h*(C1*dydx+C3*ak3+C4*ak4+C6*ak6)
	yerr=h*(DC1*dydx+DC3*ak3+DC4*ak4+DC5*ak5+DC6*ak6)
	END SUBROUTINE rkck
