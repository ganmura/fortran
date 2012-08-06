************************************************************************
* DIFFU2.F
*  ONE-DIMENSIONAL DIFFUSION EQUATION.
*  BACKWARD SCHEME
************************************************************************
      PARAMETER(MZ=50,ONEDAY=24.0*3600.0)
      PARAMETER(IUT=20)
      REAL KH ! DIFFUSIVITY
      COMMON /PARA/DZ,DT,KH

      DIMENSION FB(MZ),F(MZ)

      DZ=1.0
      DT=30.0
      KH=1.0E-4

      TTIME=ONEDAY*5.0
      NSTEP=INT(TTIME/DT)
      PTIME=3600.0*4.0
      IPRT=INT(PTIME/DT)
      STIME=ONEDAY
      IPRS=INT(STIME/DT)

      WRITE(*,*)'TTIME=',TTIME
      WRITE(*,*)'NSTEP=',NSTEP
      WRITE(*,*)'PTIME=',PTIME
      WRITE(*,*)'IPRT=',IPRT

      OPEN(IUT,FILE='TSERIES.DAT',STATUS='UNKNOWN')
      TIME=0.0
      CALL INIT(FB,MZ)
      CALL PRITS(IUT,TIME,FB,MZ)
      CALL PRISH(TIME,NFL,FB,MZ)

      DO N=1,NSTEP
        CALL STEP(FB,F,MZ)
        DO K=1,MZ !UPDATE
          FB(K)=F(K)
        END DO
        IF(MOD(N,IPRT).EQ.0)THEN !TIME SERIES
          TIME=FLOAT(N)*DT/ONEDAY
          CALL PRITS(IUT,TIME,FB,MZ)
        END IF
        IF(MOD(N,IPRS).EQ.0)THEN !SNAPSHOT
          NFL=NFL+1
          TIME=FLOAT(N)*DT/ONEDAY
          CALL PRISH(TIME,NFL,FB,MZ)
        END IF
      END DO

      CLOSE(IUT)
      STOP
      END

************************************************************************
* STEP
* TASK
*  STEP FORWARD IN TIME
*  GOVERNING EQUATION:
*   dF/dt = d(KH*(dF/dz))/dz.
* METHOD
*  BACKWARD SCHEME
* REMARK
CJ 鉛直1次元乱流モデルを使うための準備
CJ 後の回で使う乱流モデルのものと同一の解法を採用している。
* REFERENCES
*  Mellor, G.L. (1996): User's guide for a three-dimensional, primitive
*    equation, numerical ocean model. Program in atmospheric and oceanic
*    sciences, Princeton University, 41pp.[Available online at http://
*    www.aos.princeton.edu/WWWPUBLIC/htdocs.pom/]
*  Richtmyer,R.R.,and K.W.Morton (1967) Difference Methods for Initial-
*     Value Problems, 2nd Ed. Interscience, New York, 405pp.
************************************************************************
      SUBROUTINE STEP(FB,F,MZ)

      DIMENSION FB(MZ),F(MZ)
      REAL KH ! DIFFUSIVITY
      COMMON /PARA/DZ,DT,KH
      
      DIMENSION A(MZ),C(MZ),D(MZ),EE(MZ),GG(MZ)

      DO 20 K=2,MZ-1
      A(K-1)=-DT*KH/DZ**2
      C(K)=-DT*KH/DZ**2
      D(K)=-FB(K)
   20 CONTINUE

C BOUNDARY CONDITIONS
C K=1
      EE(1)=0.
      FSURF=0.0
      GG(1)=FSURF
C K=MZ
      EE(MZ)=0.
      FBTTM=0.0
      GG(MZ)=FBTTM

      DO 101 K=2,MZ-2
      GG(K)=1./(A(K)+C(K)*(1.-EE(K-1))-1.)
      EE(K)=A(K)*GG(K) ! Eq.(9-5a) of Mellor(1996)
      GG(K)=(C(K)*GG(K-1)+D(K))*GG(K) ! Eq.(9-5b) of Mellor(1996)
  101 CONTINUE

      DO 105 K=2,MZ-1
      KI=MZ-K
      F(KI)=(EE(KI)*F(KI+1)+GG(KI)) ! Eq.(9-4) of Mellor(1996)
  105 CONTINUE

      RETURN
      END

************************************************************************
* INIT
*  SET INITIAL CONDITIONS
************************************************************************
      SUBROUTINE INIT(FB,MZ)
      DIMENSION FB(MZ)

      DO K=1,MZ
        FB(K)=0.0
      END DO
      FB(MZ/2-1)=1.0
      FB(MZ/2)=1.0
      FB(MZ/2+1)=1.0
      RETURN
      END

************************************************************************
* PRITS
*  TIME SERIES
************************************************************************
      SUBROUTINE PRITS(IUT,TIME,F,MZ)

      DIMENSION F(MZ)
C      WRITE(*,*)time
      WRITE(IUT,'(2F10.5)')TIME,F(MZ/2)
      RETURN
      END

************************************************************************
* PRISH
*  SNAPSHOT
************************************************************************
      SUBROUTINE PRISH(TIME,NFL,F,MZ)
      DIMENSION F(MZ)
      COMMON /PARA/DZ,DT,KH
      CHARACTER SNOUT*3

      j1 = mod(NFL,10)
      j2 = mod(int(NFL/10),10)
      j3 = int(NFL/100)

      SNOUT=CHAR(J3+ICHAR('0'))//CHAR(J2+ICHAR('0'))
     &//CHAR(J1+ICHAR('0'))
     
      OPEN(30,FILE='PROF'//SNOUT//'.DAT',STATUS='UNKNOWN')
      DO K=1,MZ
        Z=DZ*(FLOAT(K)-0.5)
        WRITE(30,'(2F10.5)')Z,F(K)
      END DO
      CLOSE(30)
      RETURN
      END
