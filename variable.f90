program variable
	implicit none
	! 倍精度の実数にするためのおまじない
	integer, parameter :: dp = kind(1.0d0 )
	! integer
	integer i, j ! int 4byte
	! 文字(文字数) h
	character(len=12) h, h2
	! 実数
	!real*8 a ! 倍数度の実数
	!real a    ! 単精度の実数
	real(dp) a
	h = 'hello world'
	h2 = 'Not so bad'
	
	write(*, '(a)') 'input i ='
	read(*, *) i
!	i = 1
!	j = 6
	
!	print *, 'i / j = ', i / j 
	
	if(i>1)then
		print *, h
	else
		print *, h2
	endif

!	a = 3.14198572673294726529762962622659846609178275283
	
!	print *, a
end program variable
