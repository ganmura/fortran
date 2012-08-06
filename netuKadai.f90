! --------------------------------
! Created by M. Minowa
! Lastupdate 2012/7/25
! --------------------------------
program netuKadai
	implicit none
	integer, parameter :: dp = kind(1.0d0)
	integer, parameter :: m = 100
	
	! 変数の定義
	integer i, j, time
	real(dp) T(m), Tn(m)
	real(dp) dt, dx 
	real(dp) T0, Tm, T1
	real(dp) k, h
	real(dp) r, f
	real(dp) a(m), b(m), c(m), bb(m)
	real(dp) d(m), dd(m)
	!
	k = 1.0d0-3    ! 熱伝導率
	h = 10 	 !
	T0 = 0.0d0 	 !
	T1 = 1.0d0
	T = T0
	dx = h / (m - 1)
	dt = 0.5d0 * dx * dx / k	! 
	Tm = 1
	time = 5000
	! 熱伝導係数
	r = k * dt / dx**2
  do i = 1, m
    ! 係数の初期条件
    if (i==1) then
      a(i) = 0
    else
	    a(i) = -r
	  end if
	  b(i) = 1.0d0 + 2.0d0 * r
	  if (i==m-1) then
	    c(i) = 0
	    !print *, "c : ", c(i)
	  else
	    c(i) = -r
	    !print *, "c : ", c(i)
	  end if
    dd(i) = 0.0d0
	  d(i) = 0.0d0  ! 初期温度
	end do

	! 時間微分の始まり
	do j = 1, time
		! -------- 前進消去 -------------
		! i = 2 のとき
    bb(2) = b(2)
    dd(2) = d(2) - a(2) * d(1)
    !print *, "dd(2) : ", dd(2)
    !
		do i = 3, m-1
		  f = a(i) / bb(i-1)
		    !print *, "f : ", f
		  bb(i) = b(i) - c(i-1) * f
		    !print *, "bb : ", bb(i)
		  dd(i) = d(i) - dd(i-1) * f
			 ! print *, "dd : ", dd(i)
		end do
		! 
	  Tn(m-1) = dd(m-1) / bb(m-1)
		  !print *, "Tn(m-1)", Tn(m-1)
	  ! 初期条件と境界条件
	  Tn(1) = T1
  	Tn(m) = T0
		! -------- 後退代入 ------------
    do i= m-2, 2, -1
      Tn(i) = (dd(i) - c(i) * Tn(i+1)) / bb(i)
    end do
    ! 温度を更新する
    do i = 1, m
      d(i) = Tn(i)
    end do
		if (mod(j,500)==0) then ! 出力の間引き
			do i = 1, m
				print *, i, Tn(i) ! データの出力
			end do
		print *	    ! 空行 データ表示のため
		end if
	end do
end program netuKadai
