program pcg
  implicit none
  integer,parameter :: max_iter = 10000
  double precision,parameter :: epsilon = 1e-7

  integer :: i, iold=0, j, k, iter, n, nn, nzero
  integer :: row_ptr_index=1, argc
  integer,allocatable,dimension(:) :: row_ptr, col_idx
  double precision :: val, r_norm, b_norm, roh, roh_old, alpha, beta
  double precision,allocatable,dimension(:) :: a, b, p, q, r, x, diag
  double precision,allocatable,dimension(:) :: a_mod, b_mod, diag_inv
  character(50) :: s
  logical :: debug
  real :: t1, t2

  include "initialize.f95"
  allocate(p(n))
  allocate(q(n))

  ! Ax = b => A'x = b' (A' = D-1 A, b' = D-1 b)
  allocate(a_mod(nzero))
  allocate(b_mod(n))
  allocate(diag_inv(n))

  x = 0d0

  do i=1,n
     diag_inv(i) = 1/diag(i)

     do k = row_ptr(i), row_ptr(i+1)-1
        ! j = col_idx(k)
        a_mod(k) = diag_inv(i) * a(k)
     end do

     b_mod(i) = diag_inv(i) * b(i)
     r(i) = b_mod(i) - a_mod(i) * x(i)

     ! if (debug) then
     !    print *, "(i, a, a~, diag, diag_inv) = ", i, a(i), a_mod(i), diag(i), diag_inv(i)
     ! end if
  end do

  ! こっから先は col_idx, row_ptr不要

  ! cg
  ! do i=1,n
  !    r(i) = b_mod(i) - a_mod(i) * x(i)
  ! end do

  do iter=1, max_iter
     roh = dot_product(r,r)

     if( iter==1 ) then
        p = r
     else
        beta = roh/roh_old
        p = r + beta*p
     end if

     do i=1,n
!        q(i) = a_mod(i) * p(i)
        do k = row_ptr(i), row_ptr(i+1)-1
           j = col_idx(k)
           q(i) = q(i) + a(k) * p(j)
        end do
     end do

     alpha = roh / dot_product(p,q)     
     x = x + alpha*p
     r = r - alpha*q

     roh_old = roh

     !! 残差評価
     ! 残差ベクトル r = b - Ax
     r_norm = sqrt(dot_product(r,r))

     if (debug) then
        print *, "iter=",iter, "r=", r_norm
     endif

     if( r_norm < epsilon*b_norm ) then
        exit ! break do
     end if

  end do ! iter

  include "finalize.f95"

end program pcg
