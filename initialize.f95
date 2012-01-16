  ! integer, external :: iargc

  ! argc = iargc()

  debug = .false.

  ! if ( argc > 1 ) then
  !    debug = .true.
  !    print *, "argc=", argc
  ! end if

  if ( argc > 0 ) then
     call getarg(1, s)
     read (s,*) n
  else
     n = 900
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
     stop
  endif

  allocate(a(nzero))
  allocate(col_idx(nzero))

  allocate(row_ptr(n+1))
  allocate(b(n))
  allocate(x(n))
  allocate(r(n))
  allocate(diag(n))

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

  write(s, '("poisson.rhs.",I0,".data.mod")') n
  open (20,file=s)

  do i = 1,n
     read(20,*) val
     b(i) = val
     !   print *, val !! debug
  enddo

  b_norm = sqrt(dot_product(b,b))

  call cpu_time( t1 )
