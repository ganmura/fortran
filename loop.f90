program loop

	implicit none  !
	integer sum, i, n, j
	integer sum0
	integer lstat
	
	n = 1000

! 

! -------------- ファイルを書き込む処理 -----------------
	open(unit=10, file='output1.txt', status='unknown')
	open(unit=11, file='output2.txt', status='unknown')

	j = 3
!	do i = j, n, 2
	
	do n = 10, 200, 10
		sum = 0
		sum0 = n*(n+1)/2
		do i = 1, n
			sum = sum + i
		end do 
		
		print *, n, sum, sum0
		
		write(10, *) n, sum
		write(11, *)n, sum0
	end do
	
	close(10)
	close(11)
!	print *, 'sum = ', sum
!	print *, n

!	if(sum==sum0)then
!		print *, 'correct'
!	else
!		print *, 'wrong'
!	end if

! ---------- ファイルを読み込む ---------------------
	open(20, file='output1.txt', status='old')
	do
		read(20, *, iostat = lstat) n, sum
		if(lstat==0)then
			print *, n, sum
		else
			print *, 'end read'
			exit
		end if
	end do
	close(20)
		
end program loop
