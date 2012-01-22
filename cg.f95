program cg
  implicit none
  integer,parameter :: max_iter = 10000
  double precision,parameter :: epsilon = 1e-7

  integer :: i, iold=0, j, k, iter, n, nn, nzero
  integer :: row_ptr_index=1, argc
  integer,allocatable,dimension(:) :: row_ptr, col_idx
  double precision :: val, r_norm, b_norm, roh, roh_old, alpha, beta
  double precision,allocatable,dimension(:) :: a, b, p, q, r, x, diag
  character(50) :: s
  logical :: debug
  real :: t1, t2

  include "initialize.f95"
  allocate(p(n))
  allocate(q(n))

  x = 0d0

  ! cg
  do i=1,n
     r(i) = b(i)
     do k = row_ptr(i), row_ptr(i+1)-1
        j = col_idx(k)
        r(i) = r(i) - a(k) * x(j)
     end do
  end do

  do iter=1, max_iter
     roh = dot_product(r,r)

     if( iter==1 ) then
        p = r
     else
        beta = roh/roh_old
        p = r + beta*p
     end if

     do i=1,n
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

end program cg

