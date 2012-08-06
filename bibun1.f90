! bibun1.f90
! du/dz=z を初期条件u(0)=3 で解く
! 2012/6/14 M.Minowa

parameter (NX=10) ! 
dimension u(0:NX) ! 配列

dz=0.1 
u(0)=3  ! 初期条件

! 微分の計算
do i=0,NX-1
         z=i*dz
         du=z*dz
         u(i+1)=u(i)+du
end do

! 書きだす
do i=0,NX
         write(*,*)i*dz,u(i)
end do

end program bibun1
