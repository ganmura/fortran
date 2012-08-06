program renritu
	implicit none
	integer i, j
	integer, parameter :: m = 2
	integer, parameter :: dp = kind(1.0d0)
	real(dp) :: a(m,m) = reshape( (/2., 4., 3., 1./), (/m,m/) )
	real(dp) p(m),aa(m,m),b(m),bb(m), x(m)
    b(1) = 4; b(2) = 3

	  p(2) = a(2,1) / a(1,1)
	  print *, a(1,1),  a(1,2)
	  print *, a(2,1),  a(2,2)
	  print *, b, p(2)
	  
	  aa(2,2) = a(2,2) - p(2)*a(1,2)
	  bb(2) = b(2) - p(2)*b(1)
	  
	  print *, aa(2,2), bb(2) ! -5, -5
	
	  x(2) = bb(2) / aa(2,2)
	  x(1) = (b(1) - a(1,2)*x(2)) / a(1,1)
	  
	  print *, x(1), x(2) ! 0.5, 1.0
 
end program renritu
