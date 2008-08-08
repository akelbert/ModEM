module NLCG

!use math_constants
!use utilities
use sensmatrix
   ! inherits meascomb,  dataspace, dataFunc, solnrhs, 
   !            modelspace, soln2d

implicit none

public  :: NLCGsolver

! iteration control for the NLCG solver is initialized once
! and saved in the module to be used by most subroutines

  type  :: NLCGiterControl_t
     ! maximum number of iterations in one call to iterative solver
     integer					:: maxIter
     ! convergence criteria: return from solver if rms < rmsTol
     real (kind=selectedPrec)	:: rmsTol
     ! the condition to identify when the inversion stalls
     real (kind=selectedPrec)   :: fdiffTol
     ! exit if lambda < lambdaTol approx. 1e-4
     real (kind=selectedPrec)   :: lambdaTol
     ! set lambda_i = k * lambda_{i-1} when the inversion stalls
     real (kind=selectedPrec)   :: k
     ! the factor that ensures sufficient decrease in the line search
     real (kind=selectedPrec)   :: c
     ! restart CG every nCGmax iterations to ensure conjugacy
     integer                    :: nCGmax
     ! restart CG if orthogonality is lost (not necessarily needed)
     ! real (kind=selectedPrec)   :: delta ! 0.5
     ! the starting step for the line search
     real (kind=selectedPrec)   :: alpha_1
     ! if alpha_{i+1} < alpha_i * k_{alpha}, set alpha_{i+1} = alpha_i/2
     ! real (kind=selectedPrec)   :: alpha_k ! 0.1
     ! if alpha_{i+1} - alpha_i < tol_{alpha}, set alpha_{i+1} = alpha_i/2
     ! real (kind=selectedPrec)   :: alpha_tol ! 1.0e-2
  end type NLCGiterControl_t

  type(NLCGiterControl_t), private, save :: iterControl

Contains

!**********************************************************************
   subroutine set_NLCGiterControl(iterControl)
   !  at present this is a simple routine to setup for
   !  iteration control; allowing parameters to be read from
   !   a file would be desireable.

   type(NLCGiterControl_t), intent(inout)	:: iterControl
 
     ! maximum number of iterations in one call to iterative solver
     iterControl%maxIter = 120
     ! convergence criteria: return from solver if rms < rmsTol
     iterControl%rmsTol  = 1.05
     ! the condition to identify when the inversion stalls approx. 1e-3
     iterControl%fdiffTol = 1.0e-3
     ! exit if lambda < lambdaTol approx. 1e-4
     iterControl%lambdaTol = 1.0e-4
     ! set lambda_i = k * lambda_{i-1} when the inversion stalls
     iterControl%k = 0.1
     ! the factor that ensures sufficient decrease in the line search >=1e-4
     iterControl%c = 1.0e-4
     ! restart CG every nCGmax iterations to ensure conjugacy
     iterControl%nCGmax = 8
     ! the starting step for the line search
     iterControl%alpha_1 = 0.001

   end subroutine set_NLCGiterControl


!**********************************************************************
   subroutine printf(comment,lambda,alpha,f,rms)
   
   ! Compute the full penalty functional F
   ! Also output the predicted data and the EM solution
   ! that can be used for evaluating the gradient
   
	 character(*), intent(in)               :: comment
   real(kind=selectedPrec), intent(in)  :: lambda, alpha, f, rms

		write(*,'(a10)',advance='no') trim(comment)//':'
		write(*,'(a3,es12.6)',advance='no') ' f=',f
		write(*,'(a5,f10.6)',advance='no') ' rms=',rms
		write(*,'(a8,f9.6)',advance='no') ' lambda=',lambda
		write(*,'(a7,f9.6)') ' alpha=',alpha
      
   end subroutine printf


!**********************************************************************
   subroutine func(lambda,d,m0,mHat,F,dHat,eAll,RMS)
   
   ! Compute the full penalty functional F
   ! Also output the predicted data and the EM solution
   ! that can be used for evaluating the gradient
   
   real(kind=selectedPrec), intent(in)  :: lambda
   type(dvecMTX), intent(in)              :: d
   type(modelParam_t), intent(in)           :: m0
   type(modelParam_t), intent(in)           :: mHat
   real(kind=selectedPrec), intent(out) :: F
   type(dvecMTX), optional, intent(out)   :: dHat
   type(EMsolnMTX), optional, intent(out) :: eAll
   real(kind=selectedPrec), optional, intent(out) :: RMS

   !  local variables
   type(dvecMTX)    :: res,Nres
   type(modelParam_t) :: m,JTd
   real(kind=selectedPrec) :: SS,mNorm
   integer :: Ndata
   
   ! compute the smoothed model parameter vector	
   call CmSqrtMult(mHat,m)
   call linComb_modelParam(ONE,m,ONE,m0,m)
   
   !  create data structure to save solutions for all transmitters
   ! call create_EMSolnMTX(d,eAll)
   
   ! initialize dHat
   dHat = d
   
   !  compute predicted data for current model parameter m 
   !   also sets up forward solutions for all transmitters in eAll
   call fwdPred(m,dHat,eAll)
   
!	call write_Z_ascii(fidWrite,cfile,nPer,periods,modes, &
!			nSites,sites,allData)
 
       
   ! compute residual: res = d-dHat
   call linComb_DvecMTX(ONE,d,MinusONE,dHat,res)

   ! normalize residuals, compute sum of squares
   call CdInvMult(res,Nres)
   SS = dotProd(res,Nres)
   
   ! compute the model norm
   mNorm = dotProd_modelParam(mHat,mHat)
   
   ! penalty functional = sum of squares + scaled model norm
   F = SS + (lambda * mNorm)
   
   ! if required, compute the Root Mean Squared misfit
   if (present(RMS)) then
   	Ndata = count_DvecMTX(res)
   	RMS = sqrt(SS/Ndata)
   end if
      
   end subroutine func

!**********************************************************************
   subroutine gradient(lambda,d,m0,mHat,grad,dHat,eAll)

   !  Computes the gradient of the penalty functional,
   !  using EM solution (eAll) and the predicted data (dHat)
   !  Here, mHat denotes the non-regularized model parameter that
   !  is normally referred to as \tilde{m} = C_m^{-1/2}(m - m_0),
   !  and the gradient is computed with respect to \tilde{m}.
   !  Before calling this routine, the forward solver must be run:
   !  call create_EMsolnMTX(d,eAll) 
   !  call CmSqrtMult(mHat,m)
   !  call linComb_modelParam(ONE,m,ONE,m0,m)
   !  call fwdPred(m,dHat,eAll)

   real(kind=selectedPrec), intent(in)  :: lambda
   type(dvecMTX), intent(in)              :: d
   type(modelParam_t), intent(in)           :: m0
   type(modelParam_t), intent(in)           :: mHat
   type(modelParam_t), intent(out)          :: grad
   type(dvecMTX), intent(in)              :: dHat
   type(EMsolnMTX), intent(in)            :: eAll

   !  local variables
   type(dvecMTX)    :: res
   type(modelParam_t) :: m,JTd,CmJTd
 
   ! integer :: j, Ny, NzEarth
   
   ! compute the smoothed model parameter vector	
   call CmSqrtMult(mHat,m)
   call linComb_modelParam(ONE,m,ONE,m0,m)
   
   ! compute residual: res = d-dHat
   call linComb_DvecMTX(ONE,d,MinusONE,dHat,res)
   
   ! loop over transmitters:
   !do j = 1,res%nTx
   !   call getSize_modelParam(eAll%solns(j)%sigma,Ny,NzEarth)
   !   print *, Ny, NzEarth
   !end do
   
   ! multiply by J^T
   call CdInvMult(res)
   call JmultT(m,res,JTd,eAll)
   call CmSqrtMult(JTd,CmJTd)
   ! multiply by 2 (to be consistent with the formula)
   ! and add the gradient of the model norm
   call linComb_modelParam(MinusTWO,CmJTd,TWO*lambda,mHat,grad)
   
   end subroutine gradient

!**********************************************************************
   subroutine update_damping_parameter(lambda,mHat,F,grad)

   real(kind=selectedPrec), intent(inout)  :: lambda
   type(modelParam_t), intent(in)              :: mHat
   real(kind=selectedPrec), intent(inout)  :: F
   type(modelParam_t), intent(inout)             :: grad

   real(kind=selectedPrec) :: SS, mNorm
	 type(modelParam_t)          :: dSS

   ! compute the model norm
   mNorm = dotProd_modelParam(mHat,mHat)
   
   ! sum of squares = penalty functional - scaled model norm
   SS = F - (lambda * mNorm)

	 ! subtract the model norm derivative from the gradient of the penalty functional
   call linComb_modelParam(ONE,grad,MinusTWO*lambda,mHat,dSS)

	 ! update the damping parameter lambda
	 lambda = iterControl%k * lambda
   
   ! penalty functional = sum of squares + scaled model norm
   F = SS + (lambda * mNorm)

	 ! add the model norm derivative to the gradient of the penalty functional
   call linComb_modelParam(ONE,dSS,TWO*lambda,mHat,grad)

   end subroutine update_damping_parameter
	    
!**********************************************************************
   subroutine CdInvMult(d_in,d_out)

   ! Divides by the data covariance C_d, which is a diagonal
   ! operator. Divides by the variances (squared error bars)
   ! and scales by the number of data (degrees of freedom).

   type(dvecMTX), intent(inout)           :: d_in
   type(dvecMTX), optional, intent(out)   :: d_out
   type(dvecMTX)                          :: d
   !integer                                :: Ndata
   
    d = d_in
   
    ! divide each data component by its variance
    call normalize2_DvecMTX(d)
   
    ! divide by the number of data
    !Ndata = count_DvecMTX(d)
    !call scDivide_DvecMTX(ONE*Ndata,d)
   
   	if (present(d_out)) then
   		d_out = d
   	else
   	    d_in = d
   	end if

   end subroutine CdInvMult

   
!**********************************************************************
   subroutine CmSqrtMult(m_in,m_out)

   ! Multiplies by the square root of the model covariance,
   ! which is viewed as a smoothing operator. Intended
   ! to be used to compute m = C_m^{1/2} \tilde{m} + m_0.
   ! This is a dummy subroutine at present.

   type(modelParam_t), intent(in)              :: m_in
   type(modelParam_t), intent(out)             :: m_out
   type(modelParam_t)                          :: m
   
    m = m_in

	! apply the operator Cm^(1/2) here
	  call modelCov(m)
   
	  m_out = m

   end subroutine CmSqrtMult

!**********************************************************************
   subroutine NLCGsolver(d,lambda,m0,m,alpha)

   ! computes inverse solution minimizing penalty functional
   !   for fixed value of regularization parameter, using
   !   a variant of non-linear conjugate gradient search.
   !   Various flavours of the algorithm and of the line search
   !   can be called from this routine
   !
   !  Note about the starting model: 
   !  The starting model has to be in the smoothed model space,
   !  i.e. of the form m = C_m^{1/2} \tilde{m} + m_0.
   !  In order to compute \tilde{m} from the starting model,
   !  C_m^{-1/2} has to be implemented. To avoid this issue
   !  altogether, we are always starting from the prior,
   !  with \tilde{m} = 0. However, in general we could also
   !  start with the result of a previous search.
   
   !  d is data; on output it contains the responses for the inverse model
   !  NOTE: trying to set d = dHat on exit results in a corrupted data structure
   !  that is not readable by Matlab. Have to figure out why!..
   type(dvecMTX), intent(in)		       :: d
   !  lambda is regularization parameter
   real(kind=selectedPrec), intent(inout)  :: lambda
   !  m0 is prior model parameter
   type(modelParam_t), intent(in)		       :: m0
   !  m is solution parameter ... on input m contains starting guess
   type(modelParam_t), intent(inout)	       :: m
   !  alpha is the initial step size
   real(kind=selectedPrec), intent(inout), optional  :: alpha
   !  flavor is a string that specifies the algorithm to use
   ! character(80), intent(in)               :: flavor
   character(80)                           :: flavor = 'Cubic'

   !  local variables
   type(dvecMTX)			:: dHat, res
   type(modelParam_t)			:: mHat, m_minus_m0, grad, g, h, gPrev
   !type(NLCGiterControl_t)			:: iterControl
   !real(kind=selectedPrec)		:: value, valuePrev, rms, rmsPrev, alpha, beta
   real(kind=selectedPrec)		:: value, valuePrev, rms, rmsPrev, beta
   real(kind=selectedPrec)      :: grad_dot_h, g_dot_g, g_dot_gPrev, gPrev_dot_gPrev, g_dot_h
   integer				:: iter, nCG, nLS, nfunc
   type(EMsolnMTX)      :: eAll

   call set_NLCGiterControl(iterControl)

   ! initialize the line search
   !if (.not.present(flavor)) flavor = 'Cubic'
   if (.not.present(alpha)) then
      alpha = iterControl%alpha_1
   end if

   ! starting from the prior hardcoded by setting mHat = 0 and m = m0
   m = m0
   mHat = m0
   call zero_modelParam(mHat)
     

   !  compute the penalty functional and predicted data
   call func(lambda,d,m0,mHat,value,dHat,eAll,rms)
   call printf('START',lambda,alpha,value,rms)
	 nfunc = 1
      
   ! compute gradient of the full penalty functional
   call gradient(lambda,d,m0,mHat,grad,dHat,eAll)
   ! call writeModelParam(grad)
   
   ! initialize CG: g = - grad; h = g
   nCG = 0
   iter = 0
   call linComb_modelParam(MinusONE,grad,R_ZERO,grad,g)
   h = g

   do
      !  test for convergence ...
      if((rms.lt.iterControl%rmsTol).or.(iter.ge.iterControl%maxIter)) then
         exit
      end if
   
	  iter = iter + 1

	  ! save the values of the functional and the directional derivative
		rmsPrev = rms
	  valuePrev = value
	  grad_dot_h = dotProd_modelParam(grad,h)

	  ! at the end of line search, set mHat to the new value
	  ! mHat = mHat + alpha*h  and evaluate gradient at new mHat
	  write(*,*) 'Starting line search...'
	  select case (flavor)
	  case ('Cubic')
	  	call lineSearchCubic(lambda,d,m0,h,alpha,mHat,value,grad,rms,nLS)
	  case ('Quadratic')
	  	call lineSearchQuadratic(lambda,d,m0,h,alpha,mHat,value,grad,rms,nLS)
	  case default
        call errStop('Unknown line search requested in NLCG')
	  end select
		nfunc = nfunc + nLS
	  gPrev = g
	  call linComb_modelParam(MinusONE,grad,R_ZERO,grad,g)
	  
	  ! compute the starting step for the next line search
	  alpha = 2*(value - valuePrev)/grad_dot_h
	  
	  ! adjust the starting step to ensure superlinear convergence properties
	  alpha = min(ONE,(ONE+0.01)*alpha)
		write(*,'(a25,i5)') 'Completed NLCG iteration ',iter
    call printf('with',lambda,alpha,value,rms)
	  
	  ! if alpha is too small, we are not making progress: update lambda
	  if (abs(rmsPrev - rms) < iterControl%fdiffTol) then
	  	! update lambda, penalty functional and gradient
      call update_damping_parameter(lambda,mHat,value,grad)
			call linComb_modelParam(MinusONE,grad,R_ZERO,grad,g)
			! check that lambda is still at a reasonable value
			if (lambda < iterControl%lambdaTol) then
				write(*,'(a55)') 'Unable to get out of a local minimum. Exiting...'
				! multiply by C^{1/2} and add m_0
                call CmSqrtMult(mHat,m_minus_m0)
                call linComb_modelParam(ONE,m_minus_m0,ONE,m0,m)
                !d = dHat	
				return
			end if
	  	! restart
			write(*,'(a55)') 'Restarting NLCG with the damping parameter updated'
			call printf('to',lambda,alpha,value,rms)
	  	h = g
	  	nCG = 0
	  	cycle
	  end if	  
	  
	  g_dot_g = dotProd_modelParam(g,g)
	  g_dot_gPrev = dotProd_modelParam(g,gPrev)
	  gPrev_dot_gPrev = dotProd_modelParam(gPrev,gPrev)
	  g_dot_h = dotProd_modelParam(g,h)	   
	  
	  ! Polak-Ribiere variant
	  beta = ( g_dot_g - g_dot_gPrev )/gPrev_dot_gPrev 
	  
	  ! restart CG if the orthogonality conditions fail. Using the fact that 
		! h_{i+1} = g_{i+1} + beta * h_i. In order for the next directional 
		! derivative = -g_{i+1}.dot.h_{i+1} to be negative, the condition 
		! g_{i+1}.dot.(g_{i+1}+beta*h_i) > 0 must hold. Alternatively, books
		! say we can take beta > 0 (didn't work as well)
	  if ((g_dot_g + beta*g_dot_h > 0).and.(nCG < iterControl%nCGmax)) then
      	call linComb_modelParam(ONE,g,beta,h,h)
      	nCG = nCG + 1
	  else
   	    ! restart
				write(*,'(a45)') 'Restarting NLCG to restore orthogonality'
        h = g
        nCG = 0
   	  end if

   end do

   ! multiply by C^{1/2} and add m_0
   call CmSqrtMult(mHat,m_minus_m0)
   call linComb_modelParam(ONE,m_minus_m0,ONE,m0,m)
   !d = dHat
   write(*,'(a25,i5,a25,i5)') 'NLCG iterations:',iter,' function evaluations:',nfunc
   
   ! cleaning up
   call deall_dvecMTX(dHat)
   call deall_dvecMTX(res)
   call deall_modelParam(mHat)
   call deall_modelParam(m_minus_m0)
   call deall_modelParam(grad)
   call deall_modelParam(g)
   call deall_modelParam(h)
   call deall_modelParam(gPrev)
   call deall_EMsolnMTX(eAll)  

   end subroutine NLCGsolver

!**********************************************************************
  subroutine lineSearchQuadratic(lambda,d,m0,h,alpha,mHat,f,grad,rms,niter)

   ! Line search that imitates the strategy of Newman & Alumbaugh (2000),
   ! except without the errors. In particular, we only test the sufficient
   ! decrease (Armijo) condition (ignoring the curvature condition) and
   ! we use quadratic (not cubic) interpolation for backtracking.
   ! This strategy only requires one gradient evaluation (but so does
   ! the cubic interpolation method described in the Numerical Recipes).
   ! This is likely to be less efficient than the cubic interpolation,
   ! but it is simple to implement, and in some cases will work just as
   ! well (assuming an adequate initial step size has been chosen).
   !
   ! The initial step size is set outside of this routine (in the NLCG)
   ! but these are the good choices (ref. Michael Ferris, Chapter 3, p 59):
   ! alpha_1 = alpha_{k-1} dotProd(grad_{k-1},h_{k-1})/dotProd(grad_k,h_k})
   ! or interpolate the quadratic to f(m_{k-1}), f(m_k) and 
   ! dotProd(grad_{k-1},h_{k-1}) and find the minimizer
   ! alpha_1 = 2(f_k - f_{k-1})/dotProd(grad_{k-1},h_{k-1}),
   ! the update alpha_1 <- min(1.00,1.01 * alpha_1).
   !
   ! Set f(alpha) = func(mHat + alpha*h). Fit the quadratic
   !     f_q(alpha) = a alpha^2 + b alpha + f(0)
   ! using the information f(0), f'(0) and f(alpha_1) to obtain
   ! a = (f(alpha_1) - f(0) - f'(0) alpha_1)/(alpha_1 * alpha_1),
   ! b = f'(0).
   ! Then, the minimum point of the quadratic is alpha_q = -b/(2a),
   ! assuming that a > 0. If this try is not successful, fit another
   ! quadratic using f(0), f'(0) and f(alpha_q). The new quadratic
   ! is not identical to the previous quadratic since f_q is only
   ! an approximation to f: in general, f(alpha_q) /= f_q(alpha_q),
   ! hence the new point does not lie on the same quadratic curve.
   !
   ! Our solution has to satisfy the sufficient decrease condition
   !     f(alpha) < f(0) + c alpha f'(0).
   
   real(kind=selectedPrec), intent(in)     :: lambda
   type(dvecMTX), intent(in)		       :: d
   type(modelParam_t), intent(in)		       :: m0
   type(modelParam_t), intent(in)            :: h  ! search direction
   real(kind=selectedPrec), intent(inout)  :: alpha ! step size
   type(modelParam_t), intent(inout)         :: mHat
   real(kind=selectedPrec), intent(inout)  :: f
   type(modelParam_t), intent(inout)         :: grad
   real(kind=selectedPrec), intent(out)    :: rms
   integer,intent(out)                     :: niter

   ! local variables
   real(kind=selectedPrec)                 :: alpha_1,alpha_i
   logical                                 :: starting_guess
   real(kind=selectedPrec)                 :: eps,k,c,a,b
   real(kind=selectedPrec)                 :: g_0,f_0,f_1,f_i,rms_1
   type(modelParam_t)                        :: mHat_0,mHat_1
   type(dvecMTX)                           :: dHat,dHat_1
   type(EMsolnMTX)                         :: eAll,eAll_1

   ! parameters
   c = iterControl%c
   !k = iterControl%alpha_k
   !eps = iterControl%alpha_tol

   ! initialize the line search     
   niter = 0
   mHat_0 = mHat
   f_0 = f
   starting_guess = .false.
   
   ! rescale the search direction
   !h_dot_h = dotProd_modelParam(h,h)
   !h = scMult_modelParam(ONE/sqrt(h_dot_h),h)
   
   ! g_0 is the directional derivative f'(0) = (df/dm).dot.h 
   g_0 = dotProd_modelParam(grad,h)
   
   ! alpha_1 is the initial step size, which will in future be set in NLCG
   !alpha_1 = ONE/maxNorm_modelParam(h)
   alpha_1 = alpha
      
   !  compute the trial parameter mHat_1
   call linComb_modelParam(ONE,mHat_0,alpha_1,h,mHat_1)
 
   !  compute the penalty functional and predicted data at mHat_1
   call func(lambda,d,m0,mHat_1,f_1,dHat_1,eAll_1,rms_1)
   call printf('STARTLS',lambda,alpha,f_1,rms_1)
   niter = niter + 1

	 if (f_1 - f_0 >= LARGE_REAL) then
		print *, 'Try a smaller starting value of alpha.'
		print *, 'Exiting...'
		stop
	 end if
   
   f_i = f_1
   alpha_i = alpha_1
   
   fit_quadratic: do
    a = (f_i - f_0 - g_0*alpha_i)/(alpha_i**2)
    b = g_0
    alpha = - b/(TWO*a) ! the minimizer of the quadratic
    ! if the quadratic has negative curvature & no minimum, exit
    if (a < 0) then
    	starting_guess = .true.
    	exit
    end if
	!	The step size alpha should not be adjusted manually at all!
	! Even when it is too small or too close to the previous try, 
	! adjusting it won't result in an improvement (it's better to exit 
	! the line search in that case, if anything)...
  !  if ((alpha_i - alpha < eps).or.(alpha < k*alpha_i)) then
  !  	alpha = alpha_i/TWO ! reset alpha to ensure progress
  !  end if 
    call linComb_modelParam(ONE,mHat_0,alpha,h,mHat)  
    call func(lambda,d,m0,mHat,f,dHat,eAll,rms)
    call printf('QUADLS',lambda,alpha,f,rms)
    niter = niter + 1
    ! check whether the solution satisfies the sufficient decrease condition
    if (f < f_0 + c * alpha * g_0) then
    	exit
    end if
    ! if not, iterate, using the most recent value of f & alpha
    alpha_i = alpha
    f_i = f 
   end do fit_quadratic
   
   ! if the initial guess was better than what we found, take it
   if (f_1 < f) then
   	starting_guess = .true.
   end if
   
   if (starting_guess) then
   	alpha = alpha_1
   	dHat = dHat_1
   	eAll = eAll_1
   	mHat = mHat_1
   	rms = rms_1
   	f = f_1
   end if
   
   ! compute gradient of the full penalty functional and exit
   call gradient(lambda,d,m0,mHat,grad,dHat,eAll)
	 print *, 'Gradient computed, line search finished'
     
  end subroutine lineSearchQuadratic
  
  
  !**********************************************************************
  subroutine lineSearchCubic(lambda,d,m0,h,alpha,mHat,f,grad,rms,niter)

   ! Line search that is based on the Numerical Recipes and on the
   ! text by Michael Ferris, Chapter 3, p 59. We only test the sufficient
   ! decrease (Armijo) condition (ignoring the curvature condition).
   ! We first interpolate using a quadratic approximation; if the
   ! solution does not satisfy the condition, we backtrack using
   ! cubic interpolation. This strategy only requires one gradient 
   ! evaluation and is very efficient when computing gradients is
   ! expensive.
   !
   ! The initial step size is set outside of this routine (in the NLCG)
   ! but these are the good choices (ref. Michael Ferris, Chapter 3, p 59):
   ! alpha_1 = alpha_{k-1} dotProd(grad_{k-1},h_{k-1})/dotProd(grad_k,h_k})
   ! or interpolate the quadratic to f(m_{k-1}), f(m_k) and 
   ! dotProd(grad_{k-1},h_{k-1}) and find the minimizer
   ! alpha_1 = 2(f_k - f_{k-1})/dotProd(grad_{k-1},h_{k-1}),
   ! the update alpha_1 <- min(1.00,1.01 * alpha_1).
   !
   ! Set f(alpha) = func(mHat + alpha*h). Fit the quadratic
   !     f_q(alpha) = a alpha^2 + b alpha + f(0)
   ! using the information f(0), f'(0) and f(alpha_1) to obtain
   ! a = (f(alpha_1) - f(0) - f'(0) alpha_1)/(alpha_1 * alpha_1),
   ! b = f'(0).
   ! Then, the minimum point of the quadratic is alpha_q = -b/(2a),
   ! assuming that a > 0. If this try is not successful, fit a cubic
   !     f_c(alpha) = a alpha^3 + b alpha^2 + f'(0) alpha + f(0)
   ! using f(0), f'(0), f(alpha_1) and f(alpha_q). Repeat as necessary.
   ! Here, a and b are as described in the code.
   ! A new cubic is not identical to a previous curve since f_c is only
   ! an approximation to f: in general, f(alpha_c) /= f_c(alpha_c),
   ! hence the new point does not lie on the approximating curve.
   !
   ! Our solution has to satisfy the sufficient decrease condition
   !     f(alpha) < f(0) + c alpha f'(0).
   
   real(kind=selectedPrec), intent(in)     :: lambda
   type(dvecMTX), intent(in)		       :: d
   type(modelParam_t), intent(in)		       :: m0
   type(modelParam_t), intent(in)            :: h  ! search direction
   real(kind=selectedPrec), intent(inout)  :: alpha ! step size
   type(modelParam_t), intent(inout)         :: mHat
   real(kind=selectedPrec), intent(inout)  :: f
   type(modelParam_t), intent(inout)         :: grad
   real(kind=selectedPrec), intent(out)    :: rms
   integer, intent(out)                    :: niter

    ! local variables
   real(kind=selectedPrec)                 :: alpha_1,alpha_i,alpha_j
   logical                                 :: starting_guess
   real(kind=selectedPrec)                 :: eps,k,c,a,b,q1,q2,q3
   real(kind=selectedPrec)                 :: g_0,f_0,f_1,f_i,f_j,rms_1
   type(modelParam_t)                        :: mHat_0,mHat_1
   type(dvecMTX)                           :: dHat,dHat_1
   type(EMsolnMTX)                         :: eAll,eAll_1

   ! parameters
   c = iterControl%c
   !k = iterControl%alpha_k
   !eps = iterControl%alpha_tol

   ! initialize the line search     
   niter = 0
   mHat_0 = mHat
   f_0 = f
   starting_guess = .false.
   
   ! rescale the search direction
   !h_dot_h = dotProd_modelParam(h,h)
   !h = scMult_modelParam(ONE/sqrt(h_dot_h),h)
   
   ! g_0 is the directional derivative f'(0) = (df/dm).dot.h 
   g_0 = dotProd_modelParam(grad,h)
   
   ! alpha_1 is the initial step size, which will in future be set in NLCG
   alpha_1 = alpha
     
   ! compute the trial mHat, f, dHat, eAll, rms
   call linComb_modelParam(ONE,mHat_0,alpha_1,h,mHat_1)
   call func(lambda,d,m0,mHat_1,f_1,dHat_1,eAll_1,rms_1)
   call printf('STARTLS',lambda,alpha,f_1,rms_1)
   niter = niter + 1
   
	 if (f_1 - f_0 >= LARGE_REAL) then
		print *, 'Try a smaller starting value of alpha.'
		print *, 'Exiting...'
		stop
	 end if

   ! try fitting a quadratic
   a = (f_1 - f_0 - g_0*alpha_1)/(alpha_1**2)
   b = g_0
   ! if the curvature is -ve, there is no minimum; take the initial guess
   if (a < 0) then
	starting_guess = .true.
  	alpha = alpha_1
   	dHat = dHat_1
   	eAll = eAll_1
   	mHat = mHat_1
   	rms = rms_1
   	f = f_1
    ! compute the gradient and exit
    call gradient(lambda,d,m0,mHat,grad,dHat,eAll)
		print *, 'Gradient computed, exiting line search'
   	return
   end if
   
   ! otherwise compute the functional at the minimizer of the quadratic
   alpha = - b/(TWO*a)
   call linComb_modelParam(ONE,mHat_0,alpha,h,mHat)  
   call func(lambda,d,m0,mHat,f,dHat,eAll,rms)
   call printf('QUADLS',lambda,alpha,f,rms)
   niter = niter + 1
   ! check whether the solution satisfies the sufficient decrease condition
   if (f < f_0 + c * alpha * g_0) then
    ! if the initial guess was better than what we found, take it
   	if (f_1 < f) then
   		starting_guess = .true.
   		alpha = alpha_1
   		dHat = dHat_1
   		eAll = eAll_1
   		mHat = mHat_1
   		rms = rms_1
   		f = f_1
    end if
    ! compute the gradient and exit
    call gradient(lambda,d,m0,mHat,grad,dHat,eAll)
	  print *, 'Gradient computed, exiting line search'
   	return
   end if    


   ! fit a cubic and backtrack (initialize)
   alpha_i = alpha_1
   f_i = f_1
   alpha_j = alpha
   f_j = f
      
   fit_cubic: do
    ! compute the minimizer
   	q1 = f_i - f_0 - g_0 * alpha_i
   	q2 = f_j - f_0 - g_0 * alpha_j
   	q3 = alpha_i**2 * alpha_j**2 * (alpha_j - alpha_i)
   	a = (alpha_i**2 * q2 - alpha_j**2 * q1)/q3
   	b = (alpha_j**3 * q1 - alpha_i**3 * q2)/q3
   	alpha = (- b + sqrt(b*b - 3*a*g_0))/(3*a)
 	! if alpha is too close or too much smaller than its predecessor
  !  if ((alpha_j - alpha < eps).or.(alpha < k*alpha_j)) then
  !  	alpha = alpha_j/TWO ! reset alpha to ensure progress
  !  end if 
	! compute the penalty functional
    call linComb_modelParam(ONE,mHat_0,alpha,h,mHat)  
    call func(lambda,d,m0,mHat,f,dHat,eAll,rms)
    call printf('CUBICLS',lambda,alpha,f,rms)
    niter = niter + 1
    ! check whether the solution satisfies the sufficient decrease condition
    if (f < f_0 + c * alpha * g_0) then
    	exit
    end if
    ! if not, iterate, using the two most recent values of f & alpha
    alpha_i = alpha_j
    f_i = f_j
    alpha_j = alpha
    f_j = f
   end do fit_cubic
   
   if (f_1 < f) then
   	starting_guess = .true.
   end if
   
   ! if the initial guess was better than what we found, take it
   if (starting_guess) then
   	alpha = alpha_1
   	dHat = dHat_1
   	eAll = eAll_1
   	mHat = mHat_1
   	rms = rms_1
   	f = f_1
   end if
   
   ! compute gradient of the full penalty functional and exit
   call gradient(lambda,d,m0,mHat,grad,dHat,eAll)
	 print *, 'Gradient computed, line search finished'

     
  end subroutine lineSearchCubic

end module NLCG
