program jacobi2
  implicit none
  integer,parameter :: max_iter = 100000
  double precision,parameter :: epsilon = 1e-7

  integer :: i,iold=0, j, k, iter, n, nn, nzero, argc
  integer :: row_ptr_index=1
  integer,allocatable,dimension(:) :: row_ptr, col_idx
  double precision :: res, val, r_norm, b_norm
  double precision,allocatable,dimension(:) :: a, b, x, xold, r, diag
  character(50) :: s
  logical :: debug
  real :: t1, t2

  include "initialize.f95"

  allocate(xold(n))
  xold = 0d0

  !! start main loop
  do iter=1, max_iter

     do i=1,n

        x(i) = b(i)

        ! k : 非ゼロ要素のval中の序列
        do k = row_ptr(i), row_ptr(i+1)-1
           j = col_idx(k)
           if ( i .ne. j ) then
              x(i) = x(i) - a(k)*xold(j)
           end if
        end do
        
        x(i) = x(i) / diag(i)

     end do

     !! 残差評価
     ! 残差ベクトル r = b - Ax
     do i=1,n
        r(i) = b(i)

        do k = row_ptr(i), row_ptr(i+1)-1
           j = col_idx(k)
           r(i) = r(i) - a(k) * x(j)
        enddo
     end do

     r_norm = sqrt(dot_product(r,r))

     if ( debug .and. mod(iter,100) == 1 ) then
        print *, "iter=",iter, "r", r_norm
     end if

     if( r_norm < epsilon*b_norm ) then
        exit ! break do
     end if

     xold = x

  end do ! iter
  
  include "finalize.f95"

end program jacobi2
