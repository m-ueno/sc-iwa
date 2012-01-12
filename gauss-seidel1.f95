program gauss_seidel1
  implicit none
  integer,parameter :: n = 900, max_iter = 100000
  integer :: i,j,k,nn,nzero,debug=0
  double precision,parameter :: epsilon = 1e-7
  double precision :: val, res, r_norm, b_norm
  double precision :: a(n,n),b(n),x(n),xold(n),r(n)

  debug = iargc()

  a = 0d0
  open (10,file='poisson.matrix.900.data')
  read(10,*) nn,nzero             !900,3924
  do k=1,nzero
     read(10,*) i,j,val
     a(j,i)=val
  enddo

  open (20,file='poisson.rhs.900.data.mod')
  do i = 1,n
     read(20,*) val
     b(i) = val
!     print *, "b(i)",val !! debug
     xold(i) = 0d0
  enddo

  b_norm = sqrt(dot_product(b,b))

  !! start main loop
  do k=1, max_iter

     do i=1,n

        x(i) = b(i)

        do j=1,i-1
           x(i) = x(i) - a(j,i)*x(j)
        end do

        do j=i+1,n
           x(i) = x(i) - a(j,i)*xold(j)
        enddo

        x(i) = x(i) / a(i,i)

     end do

     !! 残差評価
     ! 残差ベクトル r(i) = -Ax + b ( = b - Ax )
     do i=1,n
        r(i) = b(i)
        do j=1,n
           r(i) = r(i) - a(j,i)*x(j)
        end do
     end do

     ! 相対残差ノルム
     r_norm = sqrt(dot_product(r,r))

     if ( mod(k,100) == 1 ) then
        print *, "k=",k, "r=", r_norm
     end if
    
     if( r_norm < epsilon*b_norm ) then
        exit ! break do
     end if

     xold = x

  end do ! k-loop

  !! print
  if (debug>0) then
     do i=1,n
        if (x(i) > 0.01) then
           print *, i, x(i)
        end if
     end do
  end if

  print *, "k=",k

end program gauss_seidel1
