module DCG

use math_constants
use utilities
use sensmatrix
   ! inherits datasens,  dataspace, dataFunc, solnrhs,
   !            modelspace, soln2d

! iteration control for PCG solver

  type  :: iterControl_t
     ! maximum number of iterations in one call to iterative solver
     integer					:: maxIt
     ! convergence criteria: return from solver if relative error < tol
     real (kind=prec)			:: tol
     ! actual number of iterations before return
     integer					:: niter
     ! relative error for each iteration
     real (kind=prec), pointer, dimension(:)	:: rerr
     ! logical variable indicating if algorithm "failed"
     logical					:: failed = .false.
  end type iterControl_t

 save
   type(EMsolnMTX_t), private		:: eAll
   real(kind=prec), private	:: lambda
   type(modelParam_t), private		:: sigma, JTd

Contains

!**********************************************************************
      subroutine multA(d,Ad)

   !  Multiplies by matrix A = ( J Cm J^T + lambda I )
   !  needed for data space CG algorithm
   !  eAll, lambda, sigma must be set before using
   !  JTd must be allocated before calling

   !use wscovar, only: CmMult => solveDiff
	 use modelparameter, only: CmMult => solveDiff
   type(dataVecMTX_t), intent(in)        :: d
   type(dataVecMTX_t), intent(out)       :: Ad

   !  local variables
   integer      :: i,j

   ! multiply by J Cm J^T
   call JmultT(sigma,d,JTd,eAll)
   call CmMult(JTd)
   call Jmult(JTd,sigma,Ad,eAll)

   ! add lambda to result  (need to know about details of dataVec)
   do j = 1,d%nTx
     do i = 1,d%d(j)%nDt
       Ad%d(j)%data(i)%value = Ad%d(j)%data(i)%value + lambda
     enddo
   enddo

   end subroutine multA

!**********************************************************************
   subroutine Minv(d,Md)
   !  precondtioner for PCG routine
   !   dummy routine ... no preconditioner at present

   type(dataVecMTX_t), intent(in)        :: d
   type(dataVecMTX_t), intent(out)       :: Md

   Md = d

   end subroutine Minv
!**********************************************************************
   subroutine setIterControl(PCGiter)
   !  at present this is a simple routine to setup for
   !  iteration control; allowing parameters to be read from
   !   a file would be desireable.

   type(iterControl_t), intent(inout)	:: PCGiter

   !  create data structure to save solutions for all transmitters
   call CreateSolnMTX(d,m%grid,eAll)

   PCGiter%maxit = 20
   PCGiter%tol = .001
   PCGiter%niter = 0
   allocate(PCGiter%rerr(PCGiter%maxit))

   end subroutine setIterControl


!**********************************************************************
   subroutine dataSpaceCG(d,lambda,m0,m)

   ! computes inverse solution minimizing penalty functional
   !   for fixed value of regularization parameter, by
   !   Gauss-Newton iteration.  Inner loop Gauss-Newton
   !   equations are solved in the data space with CG
   ! iteration control of solver/output of diagnostics is through
   !   module data block

   !  d is data
   type(dataVecMTX_t), intent(in)		:: d
   !  lambda is regularization parameter
   real(kind=prec)		:: lambda
   !   m0 is prior model parameter
   type(modelParam_t), intent(in)		:: m0
   !   m is solution parameter ... on input m contains starting guess
   type(modelParam_t), intent(inout)	:: m

   !  local variables
   type(dataVecMTX_t)			:: dHat, b, res
   type(modelParam_t)			:: m_minus_m0
   type(iterControl_t)			:: PCGiter
   real(kind=prec)		:: rms
   integer				:: iter, ndata

   ! these copies are just used to create data vectors with the
   !   proper structure
   dhat = d
   b = d

   !  create data structure to save solutions for all transmitters
   call create_EMSolnMTX(d,eAll)

   call setIterControl(PCGiter)

   iter = 0
   do
      !  compute predicted data for current model parameter m
      !   also sets up forward solutions for all transmitters in eAll
      Call fwdPred(m,res,eAll)

      ! compute residual: res = d-dHat
      call linComb_dataVecMTX(ONE,d,MinusONE,dHat,res)

      ! normalize data, compute rms
      call normalizeData(res,b)
      ndata = countData(d)
      rms = sqrt((b.dot.b)/ndata)

      !  test for convergence ...
      if((rms.lt.tol).or.(iter.ge.maxIter)) then
         exit
      endif

      iter = iter + 1

      ! compute dHat for next iteration
      call linComb_modelParam(ONE,m,MinusONE,m0,m_minus_m0)
      call Jmult(m_minus_m0,m,b,eAll)
      call linComb_dataVecMTX(ONE,res,ONE,b,dHat)

      ! normalize
      call normalizeData(dHat)

      ! solve data space equations with CG
      call PCG(dHat,b,PCGiter)

      ! multiply b by C_d^-1
      b%normalized = .false.
      call normalizeData(b)

      !  compute J^T b
      call JmultT(m,b,m_minus_m0,eAll)

      ! add m_0
      call linComb_modelParam(ONE,m_minus_m0,ONE,m0,m)

   enddo

   end subroutine dataSpaceCG

!**********************************************************************

subroutine PCG_dataVecMTC(b,x, PCGiter)
! Quasi-generic pre-conditioned conjugate gradient
! Actual code is generic, but the interface is not, and must be edited
!    to create a PCG which will work for different data types
!   parts that must be edited are marked by ===>
!  Also need to provide multA (to multiply by the coefficient matrix
!    for the system that is being solved with PCG) and Minv (to
!     precondition)
!  b is input rhs, and x is solution/input first guess
!  x must be allocated and initialized before calling PCG

! ===> define names of procedures for deletion, linear combinations, dot products
!        use copy (i.e., assume that we have overloaded = with this
!        procedure) to create input objects
!        dot product for input/output data type is assumed to overload .dot.
    use dataspace, only: delete => deall_dataVecMTX, &
        		   linComb => linComb_dataVecMTX, &
                           scMultAdd => scMultAdd_dataVecMTX

! ===> define names of procedures for multiplication by A
!         (solving Ax = b), and for preconditioning
!  I am actually setting the names of multiplication by A
!    and Minv inside this module   ... otherwise we seem to end up
!    with circular use arguments

  implicit none
  ! ====> Define type of input and outputs for soln, rhs (must be same type)
  type (dataVecMTX_t), intent(in)            :: b
  type (dataVecMTX_t), intent(inout)         :: x

  type (iterControl_t), intent(inout)     :: PCGiter

  ! local variables
  !======>  these local variables must be declared as same
  !           type as input (b) and output (x)
  type (dataVecMTX_t)        		:: r,z,p,q

  !  also need to change types on these for complex case
  !   We assume here that data objects are real, so result of
  !     dot product is also real
  real (kind=prec)      	:: beta,alpha,rho,rhoOld
  real (kind=prec)      	:: bnorm, rnorm
  integer               		:: i

  ! Allocation of r, z, p, q using =
  !  NOTE: "=" is overloaded by appropriate copy routine ...
  !    otherwise this won't work (it will compile, but may not run correctly!)
  r = x
  z = x
  p = x
  q = x

  call multA(x,r)
  call linComb(ONE,b,MinusONE,r,r)
  bnorm = b.dot.b
  rnorm = r.dot.r
  i = 1
  PCGiter%rerr(i) = real(rnorm/bnorm)

  loop: do while ((PCGiter%rerr(i).gt.PCGiter%tol).and.(i.lt.PCGiter%maxIt))
     call Minv(r,z)
     rho = r.dot.z
     if(i.eq.1) then
        p = z
     else
        beta = rho/rhoOld
        call linComb(ONE,z,beta,p,p)
     end if

     call multA(p,q)
     alpha = rho/(p.dot.q)
     call linComb(ONE,x,alpha,p,x)
     call linComb(ONE,r,-alpha,q,r)
     rhoOld = rho
     i = i + 1
     rnorm = r.dot.r
     PCGiter%rerr(i) = rnorm/bnorm


  end do loop

  PCGiter%niter = i

  ! deallocate all the work arrays
  call delete(r)
  call delete(z)
  call delete(p)
  call delete(q)

end subroutine PCG ! PCG

end module DCG
