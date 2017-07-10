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
C BRCT = RIGHT CAUCHY-GREEN TENSOR
C
      DIMENSION BRCT(6)
C
      PARAMETER(ZERO=0.0D0, ONE=1.0D0, TWO=2.0D0)
C
      MU=PROPS(1)
      LAMB=PROPS(2)
C
C CALCULATION OF J=DETF
C
      DETF=DFGRD1(1, 1)*DFGRD1(2, 2)*DFGRD1(3, 3)
     1 -DFGRD1(1, 2)*DFGRD1(2, 1)*DFGRD1(3, 3)
     2 +DFGRD1(1, 2)*DFGRD1(2, 3)*DFGRD1(3, 1)
     3 +DFGRD1(1, 3)*DFGRD1(3, 2)*DFGRD1(2, 1)
     4 -DFGRD1(1, 3)*DFGRD1(3, 1)*DFGRD1(2, 2)
     5 -DFGRD1(2, 3)*DFGRD1(3, 2)*DFGRD1(1, 1)
C
      BRCT(1)=DFGRD1(1, 1)**2+DFGRD1(1, 2)**2+DFGRD1(1, 3)**2
      BRCT(2)=DFGRD1(2, 1)**2+DFGRD1(2, 2)**2+DFGRD1(2, 3)**2
      BRCT(3)=DFGRD1(3, 1)**2+DFGRD1(3, 2)**2+DFGRD1(3, 3)**2
      BRCT(4)=DFGRD1(1, 1)*DFGRD1(2, 1)+DFGRD1(1, 2)*DFGRD1(2, 2)
     1 +DFGRD1(1, 3)*DFGRD1(2, 3)
      BRCT(5)=DFGRD1(1, 1)*DFGRD1(3, 1)+DFGRD1(1, 2)*DFGRD1(3, 2)
     1 +DFGRD1(1, 3)*DFGRD1(3, 3)
      BRCT(6)=DFGRD1(2, 1)*DFGRD1(3, 1)+DFGRD1(2, 2)*DFGRD1(3, 2)
     1 +DFGRD1(2, 3)*DFGRD1(3, 3)
C
      C1=ONE/DETF
      C2=LAMB/DETF
      C3=MU/DETF
      C4=MU/(TWO*DETF)
C STRESS UPDATION
      DO I=1, 3
            STRESS(I)=C1*(LAMB*LOG(DETF)-MU)+C3*BRCT(I)
      END DO
      DO I=4, 6
            STRESS(I)=C3*BRCT(I)
      END DO
C
      DO I=1, 3
            DDSDDE(I, I)=C1*(LAMB+TWO*MU*BRCT(I))
      END DO
      DDSDDE(1, 2)=C2
      DDSDDE(1, 3)=C2
      DDSDDE(2, 3)=C2
      DDSDDE(1, 6)=ZERO
      DDSDDE(2, 5)=ZERO
      DDSDDE(3, 4)=ZERO
      DDSDDE(1, 4)=C3*BRCT(4)
      DDSDDE(2, 4)=C3*BRCT(4)
      DDSDDE(1, 5)=C3*BRCT(5)
      DDSDDE(3, 5)=C3*BRCT(5)
      DDSDDE(2, 6)=C3*BRCT(6)
      DDSDDE(3, 6)=C3*BRCT(6)
      DDSDDE(4, 5)=C4*BRCT(6)
      DDSDDE(4, 6)=C4*BRCT(5)
      DDSDDE(5, 6)=C4*BRCT(4)
      DDSDDE(4, 4)=C4*(BRCT(1)+BRCT(2))
      DDSDDE(5, 5)=C4*(BRCT(1)+BRCT(3))
      DDSDDE(6, 6)=C4*(BRCT(2)+BRCT(3))
      DO I=1, 6
            DO J=1, I-1
                  DDSDDE(I, J)=DDSDDE(J, I)
            END DO
      END DO
C
C STRESS UPDATION
C      DO I=1, NTENS
C            DO J=1, NTENS
C                  STRESS(J)=STRESS(J)+DDSDDE(J, I)*DSTRAN(I)
C            END DO
C      END DO
C
      RETURN
      END