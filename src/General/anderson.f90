! **********************************************
! Anderson acceleration for fixed point iterations
! 
! Reference: https://www.nag.com/numeric/FL/manual/pdf/F08/f08zaf.pdf
!
! Needs DGGLSE to solve the Least Square method problem 
!
! min ||c-A.x|| = 0, constrained to b.x = d
! with
!
! A(i) = residue(i)
! b(i) = 1
! c(i) = 0
! d=1
!
! @author: Abhinav Jha (Continued from MooNMD)
! @date: 23.07.2018
! **********************************************

!** ************************************************************************************
! @param[in]:
! ndim = Dimension of the solution
! kdim = Number of past values needed (dimension of min_solution)
! sol_old = value of previous solution sol_old(nn)
! fpast = array previous values: fpast(:,i) = F(sol_old^(m-k))
! delta = tilde u^m - u(m-1-i)
! work = Minimum value of lwork required for optimum perfomance
! lwork = Dimension of work as declared in subprogram from which this is called
! info_error= Detects error (if any)

! @param[out]
! sol_old= Rewrite the computed solution

! ************************************************************************************** **

subroutine anderson_acceleration(&
 ndim,      &
 kdim,   &
 sol_old,      &
 fpast,   &
 delta, &
 work,    &
 lwork,   &
 info_error,  &
 alphas_x_i) BIND(C)
  USE, INTRINSIC :: ISO_C_BINDING
  IMPLICIT NONE
  ! ldb: first dimension of b_vector
  integer(4),parameter:: ldb=1
  integer(C_INT),VALUE,intent(in) :: ndim,kdim,lwork
  real(C_DOUBLE),intent(inout) :: sol_old(ndim)
  real(C_DOUBLE),intent(inout) :: alphas_x_i(kdim)
  real(C_DOUBLE),intent(inout) :: delta(ndim,kdim),fpast(ndim,kdim)
  real(C_DOUBLE),intent(inout) :: work(lwork)
  integer(C_INT),intent(out) :: info_error

! ==== end ====
  
  integer(4) :: i,j;
  
  ! minimization problem for anderson acceleration
  ! n_rows_b: Numbver of rows in vector b
  integer(4) :: n_rows_b,m_k
  REAL(8) :: sum
  REAL(8) :: b_vector(1,kdim),c_vector(1,ndim), d_vector(1),min_solution(kdim)

   
  n_rows_b=1

  ! for the first call, just compute the optimal work array size
  if(lwork.eq.-1) then 
  !    M, N, P, A, LENGHT_A, B, LENGHT_B, C, D, X, WORK, LWORK,INFO
     call DGGLSE(ndim,kdim,n_rows_b,delta,&
          ndim,b_vector,ldb,c_vector,d_vector,min_solution,work,lwork,info_error)
     return
  end if
  
  ! for successive calls
  ! set constraints
  d_vector=1.0d0  ! constrained rhs
  b_vector=1.0d0  ! constrained matrix b_vector.alphav=d_vector=1
  c_vector=0.0d0  ! linear system rhs

  m_k=kdim ! note: we assume that the arrays have length kdim

  ! solve the constrained problem

  call DGGLSE(ndim,m_k,n_rows_b,delta,&
  ndim,b_vector,ldb,c_vector,d_vector,min_solution,work,lwork,info_error)
  
  ! calculate a new guess x_{k+1}
  ! x = sum_i=1^k w_i fpast_i
  ! print *, " GGLSE solved, min_solution: "
  ! do i=1,kdim
  !   print *, i,min_solution(i)
  ! end do

  do i=1,ndim
    sum=0.0d0
    do j=1,kdim
      sum=sum+fpast(i,j)*min_solution(j)
      !print *, i,j,fpast(i,j)
    end do
    sol_old(i)=sum
  end do
  do i=1,ndim
    sum=0.0d0
    !x_i= g(x_i)-f_i
    do j=1,kdim
      sum=sum+(fpast(i,j)-delta(i,j))*min_solution(j)
      !print *, i,j,fpast(i,j)
    end do
    alphas_x_i(i)=sum
  end do
  
end subroutine






