program dif

	implicit none
	integer, parameter :: dp = kind(1.0d0 )

	!実数の配列を宣言	
	real(dp) f(4), dfdx(4), error(4)
	
	real(dp) x, dx
	real(dp) pi
	integer i
	pi = acos( -1.0d0 ) !cosの逆関数
	
	x = pi * 0.25d0
	
	open(unit=11, file='output1.txt', status='unknown')
	open(unit=12, file='output2.txt', status='unknown')
	open(unit=14, file='output3.txt', status='unknown')
		
	do i = 1, 15
		dx = pi * 0.1d0 * (0.5d0 ** i) ! 0.5^i 半分づつdxを取っていく log座標に対応
	
		dfdx(1) = ( sin (x + dx) - sin( x ) ) / dx
		dfdx(2) = ( sin (x + dx) - sin( x - dx ) ) / dx * 0.5d0 ! 2Δx : 0.5d0
		dfdx(4) = 4.d0 * dfdx(2) - ( sin (x + dx * 2.d0) - sin( x - dx * 2.d0) ) / dx * 0.25d0
		dfdx(4) = dfdx(4) / 3.0d0             ! 
		error(1) = abs(cos(x) - dfdx(1))
		error(2) = abs(cos(x) - dfdx(2))
		error(4) = abs(cos(x) - dfdx(4))
	
!		print *, dfdx(1), error(1)
!		print *, dx, error(1), dfdx(1) 
		
		write(11, '(3e16.6)')dx, error(1), dfdx(1) !書式(3e16.6)
		write(12, '(3e16.6)')dx, error(2), dfdx(2)
		write(14, '(3e16.6)')dx, error(4), dfdx(4)
	end do

	close(11)
	close(12)
	close(14)
end program dif
