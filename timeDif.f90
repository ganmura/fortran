program timeDif

! ------------------------------------------
!
! 時間微分の数値計算のサンプルコード
! Created by M. Minowa
! Last Update, 2012/6/14
!
!
! ------------------------------------------
! 
! 初期条件
! f = -1 at t = 0
! 0 =< t =< 2pi
!
! -------------------------------------------

implicit none
integer, parameter :: dp = kind(1.0d0)
integer i, t, n
real(dp) pi, rad
real(dp), Dimension (0:361) :: f, f2, f3 ! 関数fを代入する配列
real(dp) dt

! 初期値の代入
f(0) = -1.0d0

  open(unit=11, file='output1.txt', status='unknown')
  open(unit=12, file='output2.txt', status='unknown')
  pi = acos(-1.0d0) ! πを求める
  ! print *, pi
  ! print *, 2*pi
	
	n = 20
	
  do t = 0, 360 ! 0 =< t =< 2pi で繰り返す
    
    dt = 2 * pi / n
    rad = pi * t / 180 ! pi/180 = x/y
    print *, dt
    print *, f(t)
    
    f(t+1) = f(t) + dt * sin(rad) ! sin には実数のラジアンが与えられる
!    f3(t) = (f(t+1) + f(t)) /  dt
    f(t) = f(t+1)
    f2(t) = -cos(rad)
  end do

  !ファイルに書きだすためのループ
  do t = 0, 360
    write(11, *)t, f(t)
    write(12, *)t, f2(t)
  end do
  
	close(11)
	close(12)	
	
end program timeDif
