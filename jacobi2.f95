program jacobi2
  implicit none
  integer,parameter :: max_iter = 100000
  double precision,parameter :: epsilon = 1e-7

  integer :: i,iold=0, j, k, iter, n, nzero
  integer :: row_ptr_index=1
  integer,allocatable,dimension(:) :: row_ptr, col_idx
  double precision :: res, val, r_norm, b_norm
  double precision,allocatable,dimension(:) :: a, b, x, xold, r, diag

  ! initialize n, nzero
  open (10,file='poisson.matrix.900.data')
  read(10,*) n,nzero                !900, 3924

  allocate(a(nzero))
  allocate(col_idx(nzero))

  allocate(row_ptr(n+1))
  allocate(b(n))
  allocate(x(n))
  allocate(xold(n))
  allocate(r(n))
  allocate(diag(n))

  xold = 0d0

  !! CRS
  do k=1,nzero
     read(10,*) i,j,val
     a(k) = val
     col_idx(k) = j

     if ( i == j ) then
        diag(i) = val
     end if

     if (i>iold) then
        iold = i
        row_ptr(row_ptr_index) = k
        row_ptr_index = row_ptr_index + 1
     endif
  enddo

  ! trick
  row_ptr(row_ptr_index) = nzero+1

  open (20,file='poisson.rhs.900.data.mod')
  do i = 1,n
     read(20,*) val
     b(i) = val
     !   print *, val !! debug
  enddo

  b_norm = sqrt(dot_product(b,b))

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

     if ( mod(iter,100) == 1 ) then
        print *, "iter=",iter, "r", r_norm
     end if

     if( r_norm < epsilon*b_norm ) then
        exit ! break do
     end if

     xold = x

  end do ! iter

  !! print
  do i=1,n
     if (x(i) > 0.01) then
        print *, i, x(i)
     end if
  end do
  
  print *, "iter=",iter

end program jacobi2
