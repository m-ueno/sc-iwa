  call cpu_time( t2 )

  !! print
  if ( debug ) then
     do i=1,n
        if (x(i) > 0.01) then
           print *, i, x(i)
        end if
     end do
  end if

  print *, "n: ", n
  print *, "iter: ", iter
  print *, "cpu_time: ", t2-t1
