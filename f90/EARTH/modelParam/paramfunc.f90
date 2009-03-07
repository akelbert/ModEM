! *****************************************************************************
module paramfunc
	! This module contains subroutines that can be passed by reference to the
	! implementation of the parametrization from the original model parameters
	! to the resistivity on the grid

  use polpak, only: legendre_associated,spherical_harmonic,d_factorial
  implicit none

Contains

  ! ***************************************************************************
  ! * the following are the subroutines for the *unnormalized* expressions for
  ! * the spherical harmonics with *no radial dependence*; they are used in the
  ! * (externally defined) parametrizations for the inversion. The reason we
  ! * need them here is that for the operators P and Pt we need the mapping from
  ! * the original model parameters to the resistivity on the grid and backwards.
  ! * For the log10 spherical harmonic parametrization, we compute the resistivity
  ! * from the coefficient through the following expression:
  ! * $\log_10(\rho) = \sum_\i a_i \tau_i(\phi,\th,r)$,
  ! * where $\tau_i$ are the unnormalized spherical harmonics
  ! * sin(m*phi)*P_l^m(cos(theta)) and cos(m*phi)*P_l^m(cos(theta))
  ! * The coefficients and the number of harmonics can be different layer to layer.

  function plgndr(l,m,x)
  
	integer, intent(in)		:: l,m
	real(8), intent(in)		:: x
	real(8)					:: plgndr
	real(8), dimension(0:l) :: P_lm

	call legendre_associated(l,m,x,P_lm)

	plgndr = P_lm(l)

  end function plgndr 

  ! ***************************************************************************
  ! * computes the value for the (l,m) spherical harmonic term in the expansion
  real(8) function SphHarm(l,m,phi,theta) result (Y)

	integer, intent(in)				:: l,m
	real(8), intent(in)				:: phi, theta
	real(8), dimension(0:l)				:: legendre
        real(8)                                         :: factor
        integer                                         :: m_abs

        m_abs = iabs(m)

	! I am not making efficient use of this subroutine, that outputs the first
	! l+1 values. However, this usage simplifies the parametrization by allowing
	! to immediately discard the coefficients which would otherwise be multiplied
	! by a zero-valued spherical harmonic term (with l<m)
	call legendre_associated(l,m_abs,cos(theta),legendre)

 !    We are using Schmidt Seminormalized Associated Legendre Functions, defined as:
 !    Q_n^m(x) = P_n(x) for m=0
 !    Q_n^m(x) = (-1)^m sqrt ( 2 * (n-m)! / (n+m)! ) P_n^m for m>0

        factor = ( (-1)**m_abs ) * &
             sqrt ( 2.0d0 * ( d_factorial(l - m_abs) ) / ( d_factorial(l + m_abs) ) )

	Y = legendre(l)

	if (m>0) then

	  Y = Y * factor * cos(m*phi)

	else if(m<0) then

	  Y = Y * factor * sin(iabs(m)*phi)

	end if

  end function SphHarm	! unnormalized spherical harmonics


  ! ***************************************************************************
  ! * to test whether xval belongs to the interval [xmin,xmax)
  integer function I1(xmin,xmax,xval)

	real(8),intent(in)	  :: xmin,xmax,xval

	if((xval>=xmin).and.(xval<xmax)) then

	  I1=1

	else

	  I1=0

	end if

  end function I1 ! identity function

  ! ***************************************************************************
  ! * to test whether (xval,yval) belongs to the rectangle area specified 
  integer function I2(xmin,xmax,xval,ymin,ymax,yval)

	real(8),intent(in)	  :: xmin,xmax,xval
	real(8),intent(in)	  :: ymin,ymax,yval

!	if((xval>=xmin).and.(xval<xmax).and.(yval>=ymin).and.(yval<ymax)) then
	I2 = I1(xmin,xmax,xval) * I1(ymin,ymax,yval)


  end function I2 ! identity function


end module paramfunc
