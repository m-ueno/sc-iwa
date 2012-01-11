program jacobi2
  implicit none
  integer,parameter :: n = 900, max_iter = 100000, nzero=3924
  integer :: i,iold=0,j,k
  double precision,parameter :: epsilon = 1e-7
  integer :: col_ind(nzero), row_ptr(n), row_ptr_index=1, jj
  double precision :: a(nzero),b(n),x(n),xold(n),res,l,val

  !! CSR
  open (10,file='poisson.matrix.900.data')
  read(10,*) i             !900, 3924
  do k=1,nzero
     read(10,*) i,j,val
     a(k) = val
     col_ind = j
     if (i>iold) then
        iold = i
        row_ptr(row_ptr_index) = k
        row_ptr_index = row_ptr_index + 1
     endif
  enddo

  print *, "row_ptr", row_ptr
  ! ok ¤Ã¤Ý¤¤

  open (20,file='poisson.rhs.900.data.mod')
  do i = 1,n
     read(20,*) l
     b(i) = l
     !   print *, val !! debug
  enddo

  xold = 0d0

  !! start main loop
  do k=1, max_iter

     res = 0d0                  ! residual

     do i=1,n

        x(i) = b(i)

!        do j=1,n
        do j = row_ptr(i), row_ptr(i+1)-1
           jj = col_ind(j)
           if ( i .ne. j ) then
!              x(i) = x(i) - a(j,i)*xold(j)
              x(i) = x(i) - a()*xold(j)
           end if
        end do

!        x(i) = x(i) / a(i,i)
        x(i) = x(i) / a()    ! ??
        res = res + abs( x(i) - xold(i) )

     end do

     if ( mod(k,100) == 1 ) then
        print *, "res", res        !debug
     end if

     if( res < epsilon ) then
        exit ! break do
     end if

     xold = x

  end do

  !! print
  do i=1,n
     if (x(i) > 0.01) then
        print *, i, x(i)
     end if
  end do
  
  print *, "k=",k

end program jacobi2
