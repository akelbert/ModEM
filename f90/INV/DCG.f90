module DCG

	use math_constants
 	use utilities
    use senscomp
     use main
#ifdef MPI
	Use MPI_main
	use MPI_sub
#endif
implicit none
! iteration control for CG solver
  type  :: iterControl_t
     ! maximum number of iterations in one call to iterative solver
     integer					:: maxIt
     ! convergence criteria: return from solver if relative error < tol
      real(kind=prec) 			:: tol
     ! actual number of iterations before return
     integer					:: niter
     ! relative error for each iteration
      real(kind=prec) , pointer, dimension(:)	:: rerr
     ! logical variable indicating if algorithm "failed"
     logical					:: failed = .false.
  end type iterControl_t
  


public  :: DCGsolver
   ! type(EMsolnMTX_t),save    :: eAll

Contains

!**********************************************************************
   subroutine setIterControl(CGiter)
   !  at present this is a simple routine to setup for
   !  iteration control; allowing parameters to be read from
   !   a file would be desireable.

   type(iterControl_t), intent(inout)	:: CGiter
   
   CGiter%maxit = 10
   CGiter%tol = 0.001
   CGiter%niter = 0
   allocate(CGiter%rerr(CGiter%maxit))

   end subroutine setIterControl
!**********************************************************************

  subroutine DCGsolver(d,m0,m,lambda)
  
  ! Subroutine to solve the inverse problem in data space using conjugate gradients (CG)  
   
   type(dataVecMTX_t), intent(inout)		       ::d
   !  m0 is prior model parameter
   type(modelParam_t), intent(in)		       ::m0
   !  m is solution parameter ... on input m contains starting guess
   type(modelParam_t), intent(inout)	       ::m
   !  lambda is regularization parameter
   real(kind=prec) , intent(in)							   ::lambda
   
!  local variables
   type(dataVecMTX_t)			:: dHat, b,dx,d_Pred,res,Nres,JmHat
   type(modelParam_t)			:: mHat,CmJTd,Cm_mHat,Cm_mHat1
   real(kind=prec)		  		:: rms,value,rms_old
   integer						:: iter, ndata,DS_iter,CG_iter,i,j,k,iDt
   character(3)        			:: iterChar
   character(100)       		:: modelFile,dataFile
   type(iterControl_t)			:: CGiter

   
   m    	=m0
   mHat 	=m
   Cm_mHat  =m
   Cm_mHat1 =m

call zero(mHat)
   
   JmHat	=d
   dx		=d
   b		=d
   d_Pred	=d
   res		=d
   Nres		=d
call zero_dataVecMTX(JmHat)
call zero_dataVecMTX(b)

! initialize the CG control parameters
		call setIterControl(CGiter)
   

! Compute the predicted data for the current model m
#ifdef MPI
        call Master_Job_fwdPred(m,d_Pred,eAll)
#else
        call fwdPred(m,d_Pred,eAll)
#endif
      ! compute residuals: res = d-d_Pred
        call linComb(ONE,d,MinusONE,d_Pred,res)

      ! normalize residuals with the error and compute rms
        Nres=res
        call normalize_dataVecMTX(Nres,2)
        Ndata = countData(d)
        rms = sqrt(dotProd(res,Nres)/Ndata)
        write(iterChar,'(i3.3)') DS_iter
        call printf('DCG_Start:',lambda,rms)
        
do DS_iter=1,30
	! Compute the right hand side vector (b) for the CG solver.
	! b= (d-dPred)+ J(m-m0)
	
	        if (DS_iter .gt. 1 )then	    
#ifdef MPI
            JmHat=d
            call zero_dataVecMTX(JmHat)
	        call Master_job_Jmult(mHat,m,JmHat,eAll) 
#else
	        call Jmult(mHat,m,JmHat,eAll)
#endif
	
	        end if
	        
	        call linComb(ONE,res,ONE,JmHat,b)
	        call normalize_dataVecMTX(b,1)   
	        call CG_DS(b,dx,m,d,lambda,CGiter)
	
	
	        !call normalize_dataVecMTX(dx,2)
	            do i=1,dx%nTx
	             do iDt=1,dx%d(i)%nDt
	              do j=1,dx%d(i)%data(iDt)%nSite
	               do k=1,dx%d(i)%data(iDt)%nComp
	                      dx%d(i)%data(iDt)%value(k,j)=  (dx%d(i)%data(iDt)%value(k,j)/d%d(i)%data(iDt)%error(k,j))
	                      dx%d(i)%data(iDt)%errorBar=.true.
	                end do      
	              end do
	            end do                                       
	        end do      
	 
#ifdef MPI           
	                call Master_job_JmultT(m,dx,mHat,eAll)              
#else
	                call JmultT(m,dx,mHat,eAll)
#endif
	
	
	    !  Cm_mHat=multBy_Cm(mHat)	  
	    !  Cm_mHat1= multBy_CmSqrt(mHat) 
	      Cm_mHat=  multBy_Cm(mHat) 
	      mHat=Cm_mHat
	       
	     
	     call linComb_modelParam(ONE,m0,ONE,mHat,m)
	   
	   	  write(iterChar,'(i3.3)') DS_iter
	   	  modelFile = 'DCG_Model_'//iterChar//'_MPI.cpr'
	      call write_modelParam(m,trim(modelFile))
rms_old=rms	      
	! Compute the predicted data for the current model m
#ifdef MPI
            call Master_Job_fwdPred(m,d_Pred,eAll)
#else
	        call fwdPred(m,d_Pred,eAll)
#endif
	      ! compute residuals: res = d-d_Pred
	        call linComb(ONE,d,MinusONE,d_Pred,res)
	
	      ! normalize residuals with the error and compute rms
	        Nres=res
	        call normalize_dataVecMTX(Nres,2)
	        Ndata = countData(d)
	        rms = sqrt(dotProd(res,Nres)/Ndata)
	        write(iterChar,'(i3.3)') DS_iter
	        call printf('DCG_Iter: '//iterChar,lambda,rms)
    	   dataFile ='DCG_FWD_'//iterChar//'.imp'
           call write_dataVecMTX(d_Pred,trim(dataFile))


           
 ! Clean temp vectors
end do

d=d_Pred
 
end subroutine DCGsolver
 
subroutine CG_DS(b,x,m,d,lambda,CGiter)


  type (dataVecMTX_t), intent(in)	 	::b
  type (dataVecMTX_t), intent(inout) 	::x
  type(modelParam_t),  intent(in)       ::m
  type (dataVecMTX_t), intent(in)	    ::d
  real(kind=prec),     intent(in)       ::lambda
  type(iterControl_t), intent(inout)	:: CGiter
  character(3)         					::iterChar
  
  !Local
    type (dataVecMTX_t)              	:: r,p,Ap
    real(kind=prec)					 	::alpha,beta,r_norm_pre,r_norm,b_norm,error,delta_new,delta_zero,delta_old
    integer                          	::cg_iter,i,j,k,ii,iDt
    
     
 
r=b
p=r
Ap=d
b_norm=dotProd(b,b)
call zero_dataVecMTX(x)
r_norm=dotProd(r,r)

ii = 1
CGiter%rerr(ii) = r_norm/b_norm

loop: do while ((CGiter%rerr(ii).gt.CGiter%tol).and.(ii.lt.CGiter%maxIt))

             
! Compute matrix-vector product A*p and save the result in Ap  
       call MultA(p,m,d,lambda,Ap)   

         do i=1,x%nTx 
          do iDt=1, x%d(i)%nDt 
	          r%d(i)%data(iDt)%errorBar= .false.
	          p%d(i)%data(iDt)%errorBar= .false.
	          x%d(i)%data(iDt)%errorBar= .false.
	          Ap%d(i)%data(iDt)%errorBar= .false.
          end do
          
         end do  
                       
! Compute alpha: alpha= (r^T r) / (p^T Ap)    
       alpha = r_norm/dotProd(p,Ap)
       
! Compute new x: x = x + alpha*p         
       Call scMultAdd_dataVecMTX(alpha,p,x)  
                                 
! Compute new r: r = r - alpha*Ap   
       Call scMultAdd_dataVecMTX(-alpha,Ap,r) 

        
                
        r_norm_pre=r_norm
        r_norm=dotProd(r,r)
! Compute beta: beta= r_norm /r_norm_previous           
        beta=r_norm/r_norm_pre
   
! Compute new p: p = r + beta*p    
          call linComb(ONE,r,beta,p,p)
          
       ii=ii+1
       CGiter%rerr(ii) = r_norm/b_norm 
       write(6,*) 'CG-error',ii, r_norm/b_norm
  end do loop

CGiter%niter = ii

! deallocate the help vectors
    call deall_dataVecMTX(r)
    call deall_dataVecMTX(p)
    call deall_dataVecMTX(Ap)
    
end subroutine CG_DS

!###################################################################################
subroutine MultA(p,m,d,lambda,Ap)
   type(dataVecMTX_t), intent(in)          ::p
   type(dataVecMTX_t), intent(out)         ::Ap
   type(modelParam_t), intent(in)          ::m
   type (dataVecMTX_t), intent(in)	       ::d
   real(kind=prec), intent(in)             ::lambda
!Local parameters
   type(modelParam_t)                      ::JTp,CmJTp,CmJTp1
   type(dataVecMTX_t)                      ::lambdaP,p_temp
   integer                                 ::i,j,k,iDt

JTp		=m
CmJTp	=m
CmJTp1	=m
p_temp	=p
lambdaP	=p


           !call normalize_dataVecMTX(p,1)
!  Mult C^(-1/2) p              
          do i=1,p_temp%nTx
            do iDt=1,p_temp%d(i)%nDt
             do j=1,p_temp%d(i)%data(iDt)%nSite
               do k=1,p%d(i)%data(iDt)%nComp
                      p_temp%d(i)%data(iDt)%value(k,j)=  (p_temp%d(i)%data(iDt)%value(k,j)/d%d(i)%data(iDt)%error(k,j))
              end do
            end do
           end do                                        
        end do 
        
         call linComb(R_ZERO,d,ONE,p_temp,p_temp) 
! Compute   J^T  C^(-1/2) p                   
#ifdef MPI
            call Master_job_JmultT(m,p_temp,JTp,eAll)
#else
            call JmultT(m,p_temp,JTp,eAll)
#endif
! Compute  Cm  J^T  C^(-1/2) p 
            !CmJTp=multBy_Cm(JTp)
            CmJTp= multBy_Cm(JTp) 
            !CmJTp= multBy_CmSqrt(CmJTp)
            
! Compute J Cm  J^T  C^(-1/2) p = Ap 
Ap=d    
#ifdef MPI
            call Master_job_Jmult(CmJTp,m,Ap,eAll)
#else
            call Jmult(CmJTp,m,Ap,eAll)
#endif

            call scMult_dataVecMTX(lambda,p,lambdaP)
            
!Normalize: C^(-1/2)*Ap
         do i=1,Ap%nTx
            do iDt=1,Ap%d(i)%nDt
             do j=1,Ap%d(i)%data(iDt)%nSite
               do k=1,Ap%d(i)%data(iDt)%nComp
                      Ap%d(i)%data(iDt)%value(k,j)=  (Ap%d(i)%data(iDt)%value(k,j)/d%d(i)%data(iDt)%error(k,j))
              end do
            end do
           end do                                       
        end do 
!Add Cd^(-1/2)*Ap*Cd^(-1/2) to lambda*p
         do i=1,lambdaP%nTx
           do iDt=1,lambdaP%d(i)%nDt
             lambdaP%d(i)%data(iDt)%errorBar= .false.
            end do
            
         end do
          
           call linComb_dataVecMTX(ONE,Ap,ONE,lambdaP,Ap)      
        
        
!Deallocate help vectors
    call deall_modelParam(CmJTp)
    call deall_modelParam(CmJTp1)
    call deall_dataVecMTX(p_temp)
    call deall_dataVecMTX(lambdaP)
    
    
                
end subroutine MultA
!###################################################################################
   subroutine printf(comment,lambda,rms)

   ! print some comments, rms, f and lambda
  character(*), intent(in)               :: comment
  real(kind=prec), intent(in)  :: lambda, rms

		write(*,'(a20)',advance='no') trim(comment)//':'
		write(*,'(a5,f11.6)',advance='no') ' rms=',rms
		write(*,'(a8,f11.6)',advance='no') ' lambda=',lambda

   end subroutine printf
end module DCG
