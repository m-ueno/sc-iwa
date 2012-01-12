program jacobi1
  implicit none
  integer,parameter :: n = 900, max_iter = 100000
  integer :: i,j,k,nn,nzero
  double precision,parameter :: epsilon = 1e-7
  double precision :: a(n,n),b(n),x(n),xold(n),val,res

  a = 0d0
  open (10,file='poisson.matrix.900.data')
  read(10,*) nn,nzero             !900,3924
  do k=1,nzero
     read(10,*) i,j,val
     a(j,i)=val                 !caution
  enddo

  open (20,file='poisson.rhs.900.data.mod')
  do i = 1,n
     read(20,*) val
     b(i) = val
     !   print *, val !! debug
     xold(i) = 0d0
  enddo

  !! start main loop
  do k=1, max_iter

     res = 0d0                  ! residual

     do i=1,n

        x(i) = b(i)

        do j=1,n
           if ( i .ne. j ) then
              x(i) = x(i) - a(j,i)*xold(j)
           end if
        end do

        x(i) = x(i) / a(i,i)
        res = res + abs( x(i) - xold(i) )

     end do

     if ( mod(k,100) == 1 ) then
        print *, "res", res        !debug
     end if

     ! 
     
     ! if( res < epsilon ) then
     !    exit ! break do
     ! end if

     xold = x

  end do ! k-loop

  !! print
  do i=1,n
     if (x(i) > 0.01) then
        print *, i, x(i)
     end if
  end do
  
  print *, "k=",k

end program jacobi1
