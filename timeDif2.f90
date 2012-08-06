program timeDif2

! ------------------------------------------
!
! 時間微分の数値計算のサンプルコード
! Created by M. Minowa
! Last Update, 2012/6/18
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
integer i, n
real(dp) pi, f, f0
real(dp) t0, dt, t

  open(unit=11, file='output1.txt', status='unknown')
  open(unit=12, file='output2.txt', status='unknown')
  pi = acos(-1.0d0) ! πを求める
  ! print *, pi
  ! print *, 2*pi
	
	n = 20
	dt = 2*pi / n
	! 初期値の代入
	f0 = -1.0d0
  
  do i = 1, n	! 0 =< t =< 2pi で繰り返す
   
    t0 = (i-1) * dt
    f = f0 + dt * sin(t0)
    t = i * dt
    ! 後方差分(未知の t を使う)
    ! f = f0 + dt * sin(t)
    ! 中央差分(既知と未知を半分ずつ使う)
    !f = f0 + dt * (sin(t0) + sin(t)) * 0.5d0 ! 精度を向上する
    ! 前方差分(既知 t0 を使う)
     f = f0 + dt * sin(t0)
    
    f0 = f
		print *, t, f, -cos( t )
  end do

end program timeDif2
