      PROGRAM MAIN
      INCLUDE 'ABA_PARAM.INC'
 
      CHARACTER*80 CMNAME
C      REAL*16
      PARAMETER (NTENS=6,NSTATV=38,NPROPS=12)
      DIMENSION STRESS(NTENS),STATEV(NSTATV),
     1 DDSDDE(NTENS,NTENS),DDSDDT(NTENS),DRPLDE(NTENS),
     2 STRAN(NTENS),DSTRAN(NTENS),TIME(2),PREDEF(1),DPRED(1),
     3 PROPS(NPROPS),COORDS(3),DROT(3,3),DFGRD0(3,3),DFGRD1(3,3)

      CMNAME='ABAQUS_TEST' 
      DO I=1,NTENS
         DSTRAN(I) = 0.
         STRAN(I)  = 0.
         STRESS(I) = 0.
         DO J=1,6
            DDSDDE(I,J) = 0.
         END DO 
      ENDDO
      DO I=1,NSTATV
         STATEV(I)=0.
      enddo
      DSTRAN(3) = -0.00002
      PROPS(1)  = 6.8E10
      PROPS(2)  = 0.21
      PROPS(3)  = 1.26E-13
      PROPS(4)  = 3.94E-11
      PROPS(5)  = -1.26E-12
      PROPS(6)  = 2.51E-12
      PROPS(7)  = 0.11E6
      PROPS(8)  = 2.2E6
      PROPS(9)  = 0.231
      PROPS(10) = 0.
      PROPS(11) = 0.
      PROPS(12) = 0.

      ISTEP = 1000

      WRITE(*,*) '========================================'
      WRITE(*,*) '=                                      ='
      WRITE(*,*) '=                                      ='
      WRITE(*,*) '=         ABAQUS UMAT DRIVER           ='
      WRITE(*,*) '=                                      ='
      WRITE(*,*) '=                                      ='
      WRITE(*,*) '========================================'

      open(unit=666, file='result.out',status='unknown')
      WRITE(666,6660) 
      DO INC=1,ISTEP
         WRITE(6,*) '*** Loading increment ', INC, ' ***'
c         WRITE(*,*) 'STRESS = ',(STRESS(I),I=1,NTENS)
c         WRITE(*,*) 'CMNAME = ', CMNAME
c         WRITE(*,*) 'PROPS  = ',(PROPS(I),I=1,NPROPS) 
      CALL UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,
     1 RPL,DDSDDT,DRPLDE,DRPLDT,
     2 STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,
     3 NDI,NSHR,NTENS,NSTATV,PROPS,NPROPS,COORDS,DROT,PNEWDT,
     4 CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,KSTEP,KINC)
      DO I=1,NTENS
         STRAN(I)=STRAN(I)+DSTRAN(I)
      enddo
         WRITE(666,'(I5,2X,9(E10.4,2X))') INC, (STRAN(I),I=1,3),
     1 (STRESS(I),I=1,3),(STATEV(I),I=1,3) 
      END DO
      close(666)
6660  format ('# STEP   STRAN(1)    STRAN(2)    STRAN(3)   STRESS(1)   '
     1 'STRESS(2)   STRESS(1)    OMEGA(1)    OMEGA(2)    OMEGA(3)')
      END

