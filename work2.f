      SUBROUTINE UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,
     1 RPL,DDSDDT,DRPLDE,DRPLDT,
     2 STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,
     3 NDI,NSHR,NTENS,NSTATV,PROPS,NPROPS,COORDS,DROT,PNEWDT,
     4 CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,KSTEP,KINC)
C
      INCLUDE 'ABA_PARAM.INC'
C
      CHARACTER*80 CMNAME
      DIMENSION STRESS(NTENS),STATEV(NSTATV),
     1 DDSDDE(NTENS,NTENS),DDSDDT(NTENS),DRPLDE(NTENS),
     2 STRAN(NTENS),DSTRAN(NTENS),TIME(2),PREDEF(1),DPRED(1),
     3 PROPS(NPROPS),COORDS(3),DROT(3,3),DFGRD0(3,3),DFGRD1(3,3)
C DEFINITIONS
C       -----------------------------------------
C       ROMIL KADIA(16105045)
C       ANKUR MAURYA(13124)
C       ROHIT KUMAVAT(13587)
C       -----------------------------------------
C GENERATING RIGHT CAUCHY-GREEN TENSOR:
      DIMENSION BB(6)
      PARAMETER(ZERO=0.0D0, ONE=1.0D0, TWO=2.0D0)
      MU=PROPS(1)
      LAMB=PROPS(2)
C       -----------------------------------------
C XJ IS DETERMINENT OF (F)
      XJ=DFGRD1(1, 1)*DFGRD1(2, 2)*DFGRD1(3, 3)
     1 -DFGRD1(1, 2)*DFGRD1(2, 1)*DFGRD1(3, 3)
     2 +DFGRD1(1, 2)*DFGRD1(2, 3)*DFGRD1(3, 1)
     3 +DFGRD1(1, 3)*DFGRD1(3, 2)*DFGRD1(2, 1)
     4 -DFGRD1(1, 3)*DFGRD1(3, 1)*DFGRD1(2, 2)
     5 -DFGRD1(2, 3)*DFGRD1(3, 2)*DFGRD1(1, 1)
C      ------------------------------------------
      BB(1)=DFGRD1(1, 1)**2+DFGRD1(1, 2)**2+DFGRD1(1, 3)**2
      BB(2)=DFGRD1(2, 1)**2+DFGRD1(2, 2)**2+DFGRD1(2, 3)**2
      BB(3)=DFGRD1(3, 1)**2+DFGRD1(3, 2)**2+DFGRD1(3, 3)**2
      BB(4)=DFGRD1(1, 1)*DFGRD1(2, 1)+DFGRD1(1, 2)*DFGRD1(2, 2)
     1 +DFGRD1(1, 3)*DFGRD1(2, 3)
      BB(5)=DFGRD1(1, 1)*DFGRD1(3, 1)+DFGRD1(1, 2)*DFGRD1(3, 2)
     1 +DFGRD1(1, 3)*DFGRD1(3, 3)
      BB(6)=DFGRD1(2, 1)*DFGRD1(3, 1)+DFGRD1(2, 2)*DFGRD1(3, 2)
     1 +DFGRD1(2, 3)*DFGRD1(3, 3)
C       -----------------------------------------
C STRESS UPDATION
      DO I=1, 3
            STRESS(I)=BB(I)*MU/XJ+((LAMB*LOG(XJ)-MU)/XJ)
      END DO
      DO I=4, 6
            STRESS(I)=BB(I)*MU/XJ
      END DO
C       -----------------------------------------
      DO I=1, 3
            DDSDDE(I, I)=(LAMB+TWO*MU*BB(I))/XJ
      END DO
      DDSDDE(1, 2)=LAMB/XJ
      DDSDDE(1, 3)=LAMB/XJ
      DDSDDE(2, 3)=LAMB/XJ
      DDSDDE(1, 6)=ZERO
      DDSDDE(2, 5)=ZERO
      DDSDDE(3, 4)=ZERO
      DDSDDE(1, 4)=BB(4)*MU/XJ
      DDSDDE(2, 4)=BB(4)*MU/XJ
      DDSDDE(1, 5)=BB(5)*MU/XJ
      DDSDDE(3, 5)=BB(5)*MU/XJ
      DDSDDE(2, 6)=BB(6)*MU/XJ
      DDSDDE(3, 6)=BB(6)*MU/XJ
      DDSDDE(4, 5)=BB(6)*MU/(TWO*XJ)
      DDSDDE(4, 6)=BB(5)*MU/(TWO*XJ)
      DDSDDE(5, 6)=BB(4)*MU/(TWO*XJ)
      DDSDDE(4, 4)=(BB(1)+BB(2))*MU/(TWO*XJ)
      DDSDDE(5, 5)=(BB(1)+BB(3))*MU/(TWO*XJ)
      DDSDDE(6, 6)=(BB(2)+BB(3))*MU/(TWO*XJ)
      DO I=1, 6
            DO J=1, I-1
                  DDSDDE(I, J)=DDSDDE(J, I)
            END DO
      END DO
      RETURN
      END
