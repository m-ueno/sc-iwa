program cg
  implicit none
  integer,parameter :: max_iter = 10000
  double precision,parameter :: epsilon = 1e-7

  integer :: i, iold=0, j, k, iter, n, nn, nzero
  integer :: row_ptr_index=1, argc
  integer,allocatable,dimension(:) :: row_ptr, col_idx
  double precision :: val, r_norm, b_norm, roh, roh_old, alpha, beta
  double precision,allocatable,dimension(:) :: a, b, p, q, r, x
  character(50) :: s
  logical :: debug

  argc = iargc()
  print *, "argc", argc
  if ( argc > 1 ) then
     debug = .true.
  end if

  if ( argc > 0 ) then
     call getarg(1, s)
     read (s,*) n
  end if

  ! make file name str
  write(s, '("poisson.matrix.",I0,".data")') n
  if (debug) then
     print *, "s=", s
  end if

  ! initialize n, nzero
  open (10,file=s)
  read(10,*) nn,nzero                !900, 3924

  if ( n .ne. nn ) then
     return
  endif

  allocate(a(nzero))
  allocate(col_idx(nzero))

  allocate(row_ptr(n+1))
  allocate(b(n))
  allocate(p(n))
  allocate(q(n))
  allocate(r(n))
  allocate(x(n))

  x = 0d0

  !! CRS
  do k=1,nzero
     read(10,*) i,j,val
     a(k) = val
     col_idx(k) = j

     if (i>iold) then
        iold = i
        row_ptr(row_ptr_index) = k
        row_ptr_index = row_ptr_index + 1
     endif
  enddo

  ! trick
  row_ptr(row_ptr_index) = nzero+1

  write(s, '("poisson.rhs.",I0,".data.mod")') n

  open (20,file=s)
  do i = 1,n
     read(20,*) val
     b(i) = val
     !   print *, val !! debug
  enddo

  b_norm = sqrt(dot_product(b,b))

  ! cg
  do i=1,n
     r(i) = b(i)
     do k = row_ptr(i), row_ptr(i+1)-1
        j = col_idx(k)
        r(i) = r(i) - a(k) * x(j)
     end do
  end do

  !! ================================================================
  do iter=1, max_iter
     roh = dot_product(r,r)

     if( iter==1 ) then
        p = r
     else
        beta = roh/roh_old
        p = r + beta*p
     end if

     q = 0d0 ! vector
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

!     if ( mod(iter,100) == 1 ) then
        print *, "iter=",iter, "r=", r_norm
!     end if

     if( r_norm < epsilon*b_norm ) then
        exit ! break do
     end if

  end do ! iter

  !! print
  if (debug) then
     do i=1,n
        if (x(i) > 0.01) then
           print *, i, x(i)
        end if
     end do
  end if

  print *, "iter=",iter

end program cg

