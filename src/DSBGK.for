
      PROGRAM DSBGK
C    See a detailed introduction to the DSBGK method at http://arxiv.org/abs/1207.1040
C    It solves the BGK equation with noise reduction technique for single-component flows
C    Can modify it for multi-components as discussed at https://arxiv.org/pdf/1710.07795
C    Can have relaxation time changing with molecular speed for variable Prandtl numbers
C    Current setting has two tanks to reduce the end effect of channel flows
C    Easy for users who are familiar with DSMC codes due to similarity
C    Easy to develope hybrid DSBGK-DSMC solver with straightforward info exchange
C    Usually only modify subroutines PARAMETERS and TWALLGET and the variables file

      IMPLICIT NONE
      include 'variables'
      INTEGER I,DAY1(3),NOW1(3),DAY2(3),NOW2(3)
      DOUBLE PRECISION COEFFI,RF

      CALL IDATE(DAY1)
      CALL ITIME(NOW1)
      WRITE (*,10001) DAY1(2),DAY1(1),DAY1(3),NOW1(1),NOW1(2),NOW1(3)
10001 FORMAT ( 'Start at MM/DD/YY ', i2.2, '/', i2.2, '/', i4.4, 
     &         ';   TIME ', i2.2, ':', i2.2, ':', i2.2)

      OPEN(30,FILE='FLUX.DAT')
      WRITE(30,*) 'TITLE="DSBGK TRANSIENT RESULT"'
      WRITE(30,*) 'VARIABLES="NDT","U",,"V","W"'

      WRITE(*,*) '            '
      WRITE(*,*) ' CALL DSBGK SINGLEENSEMBLE SIMULATION'
      WRITE(*,*) '            '
      DO NPR=1,1  
        CALL SAMPI

        COEFFI=RF(-1)  !negative to reinitialize the seed
        DO I=1,NPR*10  !IF NOT REINITIALIZE, SOME ENSEMBLES STOP EARLIER
          COEFFI=RF(0) !EACH ENSEMBLE NEEDS DIFFERENT RF SEQUENCE
        END DO

        CALL INIT
        DO NDT=1,NDT_T
          IF(MOD(NDT,100).EQ.0) WRITE(*,99001) NPR,NDT,NM
99001     FORMAT(' DSBGK LOOP-> NPR:',I6,' NDT:',I6,' NM:',I12,' Mol')
    
          IF (NDT.GT.NDT_S.AND.MOD(NDT,NDT_D).EQ.0) THEN
            ISMP=1
          ELSE
            ISMP=0
          END IF

          IF (NDT.EQ.NDT_S+1) THEN
            CALL ITIME(NOW2)
            WRITE (*,10003) NOW2(1),NOW2(2),NOW2(3)
10003       FORMAT ( 'Sample at ', i2.2, ':', i2.2, ':', i2.2)          
          END IF

          TIME=TIME+DT
          CALL MOVE
          CALL UPDATE  

          IF (ISMP.EQ.1) CALL SAMPLE 
          IF (MOD(NDT,NDT_O).EQ.0) THEN
            CALL OUT
            CALL SUM
            WRITE(30,993) NDT,FLUX_MEAN(1),FLUX_MEAN(2),FLUX_MEAN(3)
          END IF
        END DO
      END DO
      CLOSE(30)
993   FORMAT(' ',I5,3E15.6)

      CALL IDATE(DAY2)
      CALL ITIME(NOW2)
      WRITE (*,10001) DAY1(2),DAY1(1),DAY1(3),NOW1(1),NOW1(2),NOW1(3)
      WRITE (*,10002) DAY2(2),DAY2(1),DAY2(3),NOW2(1),NOW2(2),NOW2(3)
10002 FORMAT ( 'Close at MM/DD/YY ', i2.2, '/', i2.2, '/', i4.4, 
     &         ';   TIME ', i2.2, ':', i2.2, ':', i2.2)

      STOP
      END


      SUBROUTINE ENSEMBLE

      IMPLICIT NONE
      include 'variables'
      INTEGER I
      DOUBLE PRECISION COEFFI,RF

      WRITE(*,*) '            '
      WRITE(*,*) ' CALL DSBGK ENSEMBLE SIMULATION'
      WRITE(*,*)

*----NPT is the total ensemble number
      WRITE(*,*) ' INPUT ENSEMBLE NUMBER'
      READ(*,*) NPT
      WRITE(*,*)

      CALL SAMPI

      OPEN(24,FILE='FORCE.DAT')
      WRITE(24,*) 'TITLE="NORMALIZED FORCE ON PLATES INSIDE DOMAIN"'
      WRITE(24,*) 'VARIABLES="NSMP_F","F1","F2","F3"'

      DO NPR=1,NPT
        COEFFI=RF(-1)  !negative to reinitialize the seed
        DO I=1,NPR*10  !IF NOT REINITIALIZE, SOME ENSEMBLES STOP EARLIER
          COEFFI=RF(0) !EACH ENSEMBLE NEEDS DIFFERENT RF SEQUENCE
        END DO

        CALL INIT
        DO NDT=1,NDT_S
          WRITE(*,99001) NPR,NDT,NM
99001     FORMAT(' DSBGK LOOP-> NPR:',I6,' NDT:',I6,' NM:',I12,' Mol')

          IF (NDT.EQ.NDT_S) THEN
            ISMP=1
          ELSE
            ISMP=0
          END IF          

          TIME=TIME+DT
          CALL MOVE
          CALL UPDATE
          IF (ISMP.EQ.1) CALL SAMPLE 

          IF (ISMP.EQ.1) THEN
            COEFFI=SP(5)/(NSMP_F*DT)/(DENS0*BOLTZ*TEMP0*25.D-9)
            WRITE(*,*) 'TOTAL FORCE ON SEVERAL PLATES INSIDE DOMAIN:'
            WRITE(*,*) (FORCE(I)*COEFFI,I=1,3)
            WRITE(24,998) NSMP_F, (FORCE(I)*COEFFI,I=1,3)
          END IF
998       FORMAT(' ',I7,3E15.6)

        END DO
        CALL OUT
      END DO

      CLOSE(24)

      RETURN
      END


      SUBROUTINE INIT

      IMPLICIT NONE
      include 'variables'
      INTEGER MX,MY,MZ,NNC,I,MM,M,K
      DOUBLE PRECISION U(3),REM,VMP,A,RF,AA,UX,UY,UZ,DC,FEQ

      CALL PARAMETERS

      DX_MIN=(CB(2)-CB(1))/NCX
      DY_MIN=(CB(4)-CB(3))/NCY !THE SMALLEST CELL SIZE IF NOT UNIFORM CELLS
      DZ_MIN=(CB(6)-CB(5))/NCZ                                      

      CG(1,1)=CB(1)
      CG(3,1)=CB(3)
      CG(5,1)=CB(5)
      DO MZ=1,NCZ
        DO MY=1,NCY
          DO MX=1,NCX
            NNC=(MZ-1)*NCX*NCY+(MY-1)*NCX+MX
                     
            IF (MX.EQ.1) THEN
              CG(1,NNC)=CG(1,1)
            ELSE
              CG(1,NNC)=CG(2,NNC-1)
            END IF
            CG(2,NNC)=CG(1,NNC)+DX_MIN

            IF (MY.EQ.1) THEN 
              CG(3,NNC)=CG(3,1)
            ELSE IF (MX.EQ.1) THEN
              CG(3,NNC)=CG(4,NNC-1)
            ELSE
              CG(3,NNC)=CG(3,NNC-1)
            END IF
            CG(4,NNC)=CG(3,NNC)+DY_MIN

            IF (MZ.EQ.1) THEN
              CG(5,NNC)=CG(5,1)
            ELSE IF ((MX.EQ.1).AND.(MY.EQ.1)) THEN
              CG(5,NNC)=CG(6,NNC-1)
            ELSE
              CG(5,NNC)=CG(5,NNC-1)
            END IF
            CG(6,NNC)=CG(5,NNC)+DZ_MIN
          END DO
        END DO
      END DO

      TIME=0.D0
      NDT=0

      SUM_INIT=0.D0 
      NM=0
      REM=0.D0
      VMP=DSQRT(2.D0*BOLTZ*TEMP0/SP(5))
      DO MZ=1,NCZ
        DO MY=1,NCY
          DO MX=1,NCX
            NNC=(MZ-1)*NCX*NCY+(MY-1)*NCX+MX

            IF (OBST(NNC).NE.0) THEN !solid cells
              TEMP(NNC)=TEMP0        !can be special to differentiate
              DENS(NNC)=DENS0
              VELO(1,NNC)=0.D0
              VELO(2,NNC)=0.D0
              VELO(3,NNC)=0.D0
              GOTO 110
            END IF

            DENS(NNC)=DENS0*(1.+(DRATIO-1.)*(MX-1)/(NCX-1)) !by inlet/outlet ratio
            FNUMC(NNC)=FNUM0*DENS(NNC)/DENS0
C    DSBGK can let FNUM vary with density to make particle number per cell uniform,
C    very useful advantage for high contrast of density among cells or multi-components;
C    also let FNUM vary with cell size of high contrast, ..., in multi-scale simulations.

            SUM_INIT=SUM_INIT+DENS(NNC)*(CG(2,NNC)-CG(1,NNC))*
     &               (CG(4,NNC)-CG(3,NNC))*(CG(6,NNC)-CG(5,NNC))

            TEMP(NNC)=TEMP0
            VELO(1,NNC)=VELO0(1)
            VELO(2,NNC)=VELO0(2)
            VELO(3,NNC)=VELO0(3)
            TEMPA(NNC)=0.D0
            DENSA(NNC)=0.D0
            VELOA(1,NNC)=0.D0
            VELOA(2,NNC)=0.D0
            VELOA(3,NNC)=0.D0
            NIN(1,NNC)=0.D0
            NIN(2,NNC)=0.D0
            NIN(3,NNC)=0.D0
            NIN(4,NNC)=0.D0
            NIN(5,NNC)=0.D0
            NIN(6,NNC)=0.D0

            A=REM+DENS(NNC)*(CG(2,NNC)-CG(1,NNC))*
     &    (CG(4,NNC)-CG(3,NNC))*(CG(6,NNC)-CG(5,NNC))/FNUMC(NNC)
            IF (NNC.LT.MNC) THEN
              MM=A
              REM=(A-MM)
            ELSE
              MM=NINT(A)
            END IF

            IF (MM.GT.0) THEN
              DO M=1,MM
                IF (NM.LT.MNM) THEN
                  NM=NM+1
              PP(1,NM)=CG(1,NNC)+RF(0)*(CG(2,NNC)-CG(1,NNC))
              PP(2,NM)=CG(3,NNC)+RF(0)*(CG(4,NNC)-CG(3,NNC))
              PP(3,NM)=CG(5,NNC)+RF(0)*(CG(6,NNC)-CG(5,NNC))
                  ICELL(NM)=NNC
                  DO K=1,3
                    CALL RVELC(U(K),AA,VMP)
                    PV(K,NM)=U(K)+VELO0(K)
                  END DO
                  FNUM(NM)=FNUMC(NNC)
                  UX=U(1)
                  UY=U(2)
                  UZ=U(3)
                  DC=DENS(NNC)
                  CALL FEQGET(BOLTZ,PI,SP(5),DC,TEMP0,UX,UY,UZ,FEQ)
                  FD(NM)=FEQ
                ELSE 
                  WRITE(*,*)' ERROR IN INIT: NM.GT.MNM'
                  PAUSE
                END IF
              END DO
            END IF

110         CONTINUE
          END DO
        END DO
      END DO

      CALL WALLD

      WRITE(*,10001) DX_MIN,DY_MIN,DZ_MIN
      WRITE(*,10002) DT,DT*VMP,BOLTZ
      WRITE(*,10003) DENS0*DX_MIN*DY_MIN*DZ_MIN/FNUM0,NM/(MNC+0.),FNUM0
      WRITE(*,*) 
10001 FORMAT(' DX_MIN :',E14.6,' DY_MIN :',E14.6,' DZ_MIN  :',E14.6)
10002 FORMAT(' DT     :',E14.6,' DT*VMP :',E14.6,' BOLTZ   :',E14.6)
10003 FORMAT(' NM_MIN :',E14.6,' NM/MNC :',E14.6,' FNUM0   :',E14.6)
      WRITE(*,*) ' ---Initialization is done, now start simulation:---'
      WRITE(*,*) ' '

      RETURN
      END



      SUBROUTINE MOVE

      IMPLICIT NONE
      include 'variables'
      INTEGER IFT,N,NNC1,MX,MY,MZ,NNC2,JJ
      DOUBLE PRECISION U(3),X(3),DTR,RF,FD1,FNUM1
      DOUBLE PRECISION DTUSE,FD2,FNUM2,WEIGHT

      IF (ISMP.EQ.1) NSMP_F=NSMP_F+1

      IFT=-1
      N=0
100   N=N+1
      IF (N.LE.NM) THEN
        IF (IFT.LT.0) THEN
          DTR=DT
        ELSE
          DTR=DT*RF(0)
        END IF

        U(1)=PV(1,N)
        U(2)=PV(2,N)
        U(3)=PV(3,N)
        X(1)=PP(1,N)
        X(2)=PP(2,N)
        X(3)=PP(3,N)
        FD1=FD(N)
        FNUM1=FNUM(N)

        NNC1=ICELL(N)
        MZ=(NNC1-0.1)/NCX/NCY+1
        MY=(NNC1-(MZ-1)*NCX*NCY-0.1)/NCX+1
        MX=NNC1-(MZ-1)*NCX*NCY-(MY-1)*NCX

350     CONTINUE !Need: DTR,U,X,MX,MY,MZ,NNC1 & FD1,FNUM1

        CALL PPGET(DTR,U,X,MX,MY,MZ,NNC1,NNC2,DTUSE)

        CALL F2GET(FD1,U,DTUSE,NNC1,FD2) !FD2 AFTER INTER COLLISION
        FNUM2=FNUM1*FD2/FD1

        WEIGHT=FNUM1-FNUM2
        DENSA(NNC1)=DENSA(NNC1)+WEIGHT
        TEMPA(NNC1)=TEMPA(NNC1)+WEIGHT*(U(1)*U(1)+U(2)*U(2)+U(3)*U(3))
        VELOA(1,NNC1)=VELOA(1,NNC1)+WEIGHT*U(1)
        VELOA(2,NNC1)=VELOA(2,NNC1)+WEIGHT*U(2)
        VELOA(3,NNC1)=VELOA(3,NNC1)+WEIGHT*U(3)

        IF (NNC2.EQ.0) THEN
C    Moving inside NNC1 if NNC2=0, not yet reaching at the interfaces of NNC1
C    Update all molecular variables, but PV is unchanged or already updated with U elsewhere (efficient way). 
          PP(1,N)=X(1)
          PP(2,N)=X(2)
          PP(3,N)=X(3)
          ICELL(N)=NNC1
          FD(N)=FD2
          FNUM(N)=FNUM2
          GOTO 100
        ELSE IF (NNC2.GT.0) THEN
C    Moving into NNC2. In PPGET, MX, MY, MZ are already updated to represent NNC2 and X(3) is also updated
C    U remains unchanged when passing through interface between NNC1 and NNC2
          DTR=DTR-DTUSE
          NNC1=NNC2
          FD1=FD2   
          FNUM1=FNUM2
          GOTO 350 !Need: DTR,U,X,MX,MY,MZ,NNC1 & FD1,FNUM1
        ELSE
C    Running into the BC (IB=1,2,3,-1) of cell NNC1 if NNC2<0, then remove it or reflect back into NNC1
C    New X(3) from PPGET and MX, MY, MZ & NNC1 will remain unchanged, except for IB=-1
*----IB=1 stream (DB.GE.0, vacuum), 2 solid wall, 3 symmetry, 0 pseudo in 2D case, -1 periodic, 4 cell interface  
          IF (IB(-NNC2,NNC1).EQ.1) THEN
            CALL REMOVE(N)
            GOTO 100
          ELSE IF (IB(-NNC2,NNC1).EQ.2) THEN
            NIN(-NNC2,NNC1)=NIN(-NNC2,NNC1)+FNUM2
            IF (ISMP.EQ.1) THEN
              IF ((MX-1)*(MX-NCX)*(MY-1)*(MY-NCY).NE.0) THEN !for possible internal plates
                FORCE(1)=FORCE(1)+U(1)*FNUM2
                FORCE(2)=FORCE(2)+U(2)*FNUM2
                FORCE(3)=FORCE(3)+U(3)*FNUM2
              END IF
            END IF
            CALL BCWALL(U,FD1,NNC2,NNC1)
            PV(1,N)=U(1)
            PV(2,N)=U(2)
            PV(3,N)=U(3)
            IF (ISMP.EQ.1) THEN
              IF ((MX-1)*(MX-NCX)*(MY-1)*(MY-NCY).NE.0) THEN !for possible internal plates
                FORCE(1)=FORCE(1)-U(1)*FNUM2
                FORCE(2)=FORCE(2)-U(2)*FNUM2
                FORCE(3)=FORCE(3)-U(3)*FNUM2
              END IF
            END IF
          ELSE IF (IB(-NNC2,NNC1).EQ.-1) THEN
            IF (NNC2.EQ.-1) THEN
              X(1)=CB(2)
              MX=NCX
            ELSE IF (NNC2.EQ.-2) THEN
              X(1)=CB(1)
              MX=1
            ELSE IF (NNC2.EQ.-3) THEN
              X(2)=CB(4)
              MY=NCY
            ELSE IF (NNC2.EQ.-4) THEN
              X(2)=CB(3)
              MY=1
            ELSE IF (NNC2.EQ.-5) THEN
              X(3)=CB(6)
              MZ=NCZ
            ELSE IF (NNC2.EQ.-6) THEN
              X(3)=CB(5)
              MZ=1
            END IF
            NNC1=(MZ-1)*NCX*NCY+(MY-1)*NCX+MX

            FD1=FD2 !U(3) remains unchanged
          ELSE IF (IB(-NNC2,NNC1).EQ.3) THEN
            IF (NNC2.EQ.-1.OR.NNC2.EQ.-2) THEN
              JJ=1
            ELSE IF (NNC2.EQ.-3.OR.NNC2.EQ.-4) THEN
              JJ=2
            ELSE IF (NNC2.EQ.-5.OR.NNC2.EQ.-6) THEN
              JJ=3
            END IF

            U(JJ)=-U(JJ)
            PV(JJ,N)=U(JJ)
            FD1=FD2
          ELSE
            WRITE(*,*) ' ERROR IN MOVE: IB=?',IB(-NNC2,NNC1)
            PAUSE
          END IF 

          DTR=DTR-DTUSE
          FNUM1=FNUM2
          GOTO 350 !Need: DTR,U,X,MX,MY,MZ,NNC1 & FD1,FNUM1
        END IF
      ELSE IF (IFT.LT.0) THEN
        IFT=1
*----new molecules enter
        CALL ENTER
        N=N-1
        GOTO 100
      END IF

      RETURN
      END



      SUBROUTINE ENTER
*----IB=1 stream (DB.GE.0, vacuum), 2 solid wall, 3 symmetry, 0 pseudo in 2D case, -1 periodic, 4 cell interface  
C    For IB=1, we preset DB, FNUMC, TB, and VB(0); if VB(0)=1 (active), also preset VB(1:3)
C    Here, only deal with IB=1 of stream BC

      IMPLICIT NONE
      include 'variables'
      INTEGER MX,MY,MZ,NNC,I,MM,M
      DOUBLE PRECISION UFLOW(3),DX,DY,DZ,VMP,SC,A,B,ERF,RF,FS1,FS2
      DOUBLE PRECISION QA,U,UN,PA,V1,V2,V3,DBI,TBI,FEQ

      DO MZ=1,NCZ
        DO MY=1,NCY
          DO MX=1,NCX
            IF ((MX-1)*(MX-NCX)*(MY-1)*(MY-NCY)*
     &          (MZ-1)*(MZ-NCZ).NE.0) GOTO 120 !not if CONCAVE
  
            NNC=(MZ-1)*NCY*NCX+(MY-1)*NCX+MX
            DX=CG(2,NNC)-CG(1,NNC)
            DY=CG(4,NNC)-CG(3,NNC)
            DZ=CG(6,NNC)-CG(5,NNC)
  
            DO I=1,6
              MM=0
              IF (I.EQ.1) THEN
                IF (IB(I,NNC).EQ.1) THEN
                  VMP=DSQRT(2.D0*BOLTZ*TEMPB(I,NNC)/SP(5))    
                  IF (VELOB(0,I,NNC).GT.0.D0) THEN
                    UFLOW(1)=VELOB(1,I,NNC)
                    UFLOW(2)=VELOB(2,I,NNC)
                    UFLOW(3)=VELOB(3,I,NNC)
                  ELSE
                    UFLOW(1)=VELO(1,NNC)
                    UFLOW(2)=VELO(2,NNC)
                    UFLOW(3)=VELO(3,NNC)
                  END IF
                  SC=UFLOW(1)/VMP
                  A=(EXP(-SC*SC)+SPI*SC*(1.D0+ERF(SC)))/(2.D0*SPI)
                  B=(DENSB(I,NNC)/FNUMC(NNC))*A*VMP*DT*DY*DZ              
                  MM=B
                  IF ((B-MM).GT.RF(0)) MM=MM+1
                  GOTO 100
                END IF
                GOTO 110
              ELSE IF (I.EQ.2) THEN
                IF (IB(I,NNC).EQ.1) THEN
                  VMP=DSQRT(2.D0*BOLTZ*TEMPB(I,NNC)/SP(5))    
                  IF (VELOB(0,I,NNC).GT.0.D0) THEN
                    UFLOW(1)=VELOB(1,I,NNC)
                    UFLOW(2)=VELOB(2,I,NNC)
                    UFLOW(3)=VELOB(3,I,NNC)
                  ELSE
                    UFLOW(1)=VELO(1,NNC)
                    UFLOW(2)=VELO(2,NNC)
                    UFLOW(3)=VELO(3,NNC)
                  END IF
                  SC=0.D0-UFLOW(1)/VMP
                  A=(EXP(-SC*SC)+SPI*SC*(1.D0+ERF(SC)))/(2.D0*SPI)
                  B=(DENSB(I,NNC)/FNUMC(NNC))*A*VMP*DT*DY*DZ
                  MM=B
                  IF ((B-MM).GT.RF(0)) MM=MM+1
                  GOTO 100
                END IF
                GOTO 110    
              ELSE IF (I.EQ.3) THEN
                IF (IB(I,NNC).EQ.1) THEN
                  VMP=DSQRT(2.D0*BOLTZ*TEMPB(I,NNC)/SP(5))    
                  IF (VELOB(0,I,NNC).GT.0.D0) THEN
                    UFLOW(1)=VELOB(1,I,NNC)
                    UFLOW(2)=VELOB(2,I,NNC)
                    UFLOW(3)=VELOB(3,I,NNC)
                  ELSE
                    UFLOW(1)=VELO(1,NNC)
                    UFLOW(2)=VELO(2,NNC)
                    UFLOW(3)=VELO(3,NNC)
                  END IF
                  SC=UFLOW(2)/VMP
                  A=(EXP(-SC*SC)+SPI*SC*(1.D0+ERF(SC)))/(2.D0*SPI)
                  B=(DENSB(I,NNC)/FNUMC(NNC))*A*VMP*DT*DX*DZ              
                  MM=B
                  IF ((B-MM).GT.RF(0)) MM=MM+1
                  GOTO 100
                END IF
                GOTO 110
              ELSE IF (I.EQ.4) THEN
                IF (IB(I,NNC).EQ.1) THEN
                  VMP=DSQRT(2.D0*BOLTZ*TEMPB(I,NNC)/SP(5))    
                  IF (VELOB(0,I,NNC).GT.0.D0) THEN
                    UFLOW(1)=VELOB(1,I,NNC)
                    UFLOW(2)=VELOB(2,I,NNC)
                    UFLOW(3)=VELOB(3,I,NNC)
                  ELSE
                    UFLOW(1)=VELO(1,NNC)
                    UFLOW(2)=VELO(2,NNC)
                    UFLOW(3)=VELO(3,NNC)
                  END IF
                  SC=0.D0-UFLOW(2)/VMP
                  A=(EXP(-SC*SC)+SPI*SC*(1.D0+ERF(SC)))/(2.D0*SPI)
                  B=(DENSB(I,NNC)/FNUMC(NNC))*A*VMP*DT*DX*DZ              
                  MM=B
                  IF ((B-MM).GT.RF(0)) MM=MM+1                
                  GOTO 100
                END IF
                GOTO 110    
              ELSE IF (I.EQ.5) THEN
                IF (IB(I,NNC).EQ.1) THEN
                  VMP=DSQRT(2.D0*BOLTZ*TEMPB(I,NNC)/SP(5))    
                  IF (VELOB(0,I,NNC).GT.0.D0) THEN
                    UFLOW(1)=VELOB(1,I,NNC)
                    UFLOW(2)=VELOB(2,I,NNC)
                    UFLOW(3)=VELOB(3,I,NNC)
                  ELSE
                    UFLOW(1)=VELO(1,NNC)
                    UFLOW(2)=VELO(2,NNC)
                    UFLOW(3)=VELO(3,NNC)
                  END IF
                  SC=UFLOW(3)/VMP
                  A=(EXP(-SC*SC)+SPI*SC*(1.D0+ERF(SC)))/(2.D0*SPI)
                  B=(DENSB(I,NNC)/FNUMC(NNC))*A*VMP*DT*DX*DY              
                  MM=B
                  IF ((B-MM).GT.RF(0)) MM=MM+1
                  GOTO 100
                END IF
                GOTO 110
              ELSE IF (I.EQ.6) THEN
                IF (IB(I,NNC).EQ.1) THEN
                  VMP=DSQRT(2.D0*BOLTZ*TEMPB(I,NNC)/SP(5))    
                  IF (VELOB(0,I,NNC).GT.0.D0) THEN
                    UFLOW(1)=VELOB(1,I,NNC)
                    UFLOW(2)=VELOB(2,I,NNC)
                    UFLOW(3)=VELOB(3,I,NNC)
                  ELSE
                    UFLOW(1)=VELO(1,NNC)
                    UFLOW(2)=VELO(2,NNC)
                    UFLOW(3)=VELO(3,NNC)
                  END IF
                  SC=0.D0-UFLOW(3)/VMP
                  A=(EXP(-SC*SC)+SPI*SC*(1.D0+ERF(SC)))/(2.D0*SPI)
                  B=(DENSB(I,NNC)/FNUMC(NNC))*A*VMP*DT*DX*DY              
                  MM=B
                  IF ((B-MM).GT.RF(0)) MM=MM+1                
                  GOTO 100
                END IF
                GOTO 110    
              END IF
   
100           CONTINUE
              IF (MM.GT.0) THEN
                FS1=SC+DSQRT(SC*SC+2.D0)
                FS2=0.5D0*(1.D0+SC*(2.D0*SC-FS1))
*----the above constants are required for the entering distn. of eqn (12.5)
                DO M=1,MM
                  NM=NM+1
                  IF (NM.GT.MNM) THEN
                    WRITE(*,*)' ERROR IN ENTER: NM.GT.MNM'
                    PAUSE
                  END IF
*----NM is now the index of the new molecule
                      
                  IF ((I.EQ.1).OR.(I.EQ.2)) THEN
                    PP(1,NM)=CG(I,NNC)
                    PP(2,NM)=CG(3,NNC)+RF(0)*DY
                    PP(3,NM)=CG(5,NNC)+RF(0)*DZ 
                  ELSE IF ((I.EQ.3).OR.(I.EQ.4)) THEN                      
                    PP(1,NM)=CG(1,NNC)+RF(0)*DX
                    PP(2,NM)=CG(I,NNC)
                    PP(3,NM)=CG(5,NNC)+RF(0)*DZ
                  ELSE IF ((I.EQ.5).OR.(I.EQ.6)) THEN                      
                    PP(1,NM)=CG(1,NNC)+RF(0)*DX
                    PP(2,NM)=CG(3,NNC)+RF(0)*DY
                    PP(3,NM)=CG(I,NNC)
                  END IF

                  QA=3.D0
                  IF (SC.LT.(0.D0-3.D0)) QA=ABS(SC)+1.D0
2                 U=2.D0*QA*RF(0)-QA
*----U is a potential normalised thermal velocity component
                  UN=U+SC
*----UN is a potential inward velocity component
                  IF (UN.LT.0.D0) GOTO 2
                  PA=(2.D0*UN/FS1)*EXP(FS2-U*U)
                  IF (PA.LT.RF(0)) GOTO 2
*----the inward normalised vel. component has been selected (eqn (12.5))

                  IF ((I.EQ.1).OR.(I.EQ.2)) THEN
                    CALL RVELC(V2,V3,VMP)
                    IF (I.EQ.1) THEN
                      V1=UN*VMP-UFLOW(1)
                    ELSE
                      V1=(0.D0-UN*VMP)-UFLOW(1)
                    END IF
                  ELSE IF ((I.EQ.3).OR.(I.EQ.4)) THEN                      
                    CALL RVELC(V1,V3,VMP)
                    IF (I.EQ.3) THEN
                      V2=UN*VMP-UFLOW(2)
                    ELSE
                      V2=(0.D0-UN*VMP)-UFLOW(2)
                    END IF
                  ELSE IF ((I.EQ.5).OR.(I.EQ.6)) THEN                      
                    CALL RVELC(V1,V2,VMP)
                    IF (I.EQ.5) THEN
                      V3=UN*VMP-UFLOW(3)
                    ELSE
                      V3=(0.D0-UN*VMP)-UFLOW(3)
                    END IF
                  END IF

                  DBI=DENSB(I,NNC)
                  TBI=TEMPB(I,NNC)
                  CALL FEQGET(BOLTZ,PI,SP(5),DBI,TBI,V1,V2,V3,FEQ)
C    Here, V1, V2 and V3 are thermal velocity in global frame
                  FD(NM)=FEQ
                  ICELL(NM)=NNC
                  FNUM(NM)=FNUMC(NNC)

                  PV(1,NM)=V1+UFLOW(1)
                  PV(2,NM)=V2+UFLOW(2)
                  PV(3,NM)=V3+UFLOW(3)

                END DO
              END IF
110           CONTINUE
            END DO
120         CONTINUE
          END DO
        END DO
      END DO

      RETURN
      END


      SUBROUTINE REFLECT(BOLTZ,PI,TWALL,SP5,ALPHAN,ALPHAT,
     &                   U1,U2,U3,V1,V2,V3)
C    In surface frame, 2,3 for tangential and 1 for normal, incoming U1<0 and reflected V1>0

      IMPLICIT NONE
      DOUBLE PRECISION BOLTZ,PI,TWALL,SP5,ALPHAN,ALPHAT,U1,U2,U3,V1,V2
      DOUBLE PRECISION V3,VMP,U1_N,U2_N,U3_N,R,TH,UM,ANG,U23_N,V02,V03
      DOUBLE PRECISION RF

      VMP=DSQRT(2.D0*BOLTZ*TWALL/SP5)

      U1_N=ABS(U1)/VMP
      R=DSQRT(-ALPHAN*LOG(RF(0)))
      TH=2.D0*PI*RF(0)
      UM=DSQRT(1.D0-ALPHAN)*U1_N
      V1=VMP*DSQRT(R*R+UM*UM+2.D0*R*UM*COS(TH))

      IF (1.EQ.0) THEN               !orginal scheme as in arXiv paper
        U2_N=U2/VMP
        U3_N=U3/VMP
        ANG=ATAN2(U3_N,U2_N)         !azimuthal angle \theta
        U23_N=DSQRT(U2_N*U2_N+U3_N*U3_N)
        R=DSQRT(-ALPHAT*LOG(RF(0)))
        TH=2.D0*PI*RF(0)             !phi
        UM=DSQRT(1.D0-ALPHAT)*U23_N
        V02=VMP*(UM+R*COS(TH))
        V03=VMP*R*SIN(TH)
        V2=V02*COS(ANG)-V03*SIN(ANG)
        V3=V02*SIN(ANG)+V03*COS(ANG)
      ELSE                           !simplified scheme as in arXiv paper
        R=DSQRT(-ALPHAT*LOG(RF(0)))
        TH=2.D0*PI*RF(0)
        V2=VMP*R*COS(TH)+U2*DSQRT(1.D0-ALPHAT)
        V3=VMP*R*SIN(TH)+U3*DSQRT(1.D0-ALPHAT)
      END IF

      RETURN
      END



      SUBROUTINE REMOVE(N)
*----remove molecule N and replace it by molecule NM

      IMPLICIT NONE
      include 'variables'
      INTEGER N

      PP(1,N)=PP(1,NM)
      PP(2,N)=PP(2,NM)
      PP(3,N)=PP(3,NM)

      PV(1,N)=PV(1,NM)
      PV(2,N)=PV(2,NM)
      PV(3,N)=PV(3,NM)

      ICELL(N)=ICELL(NM)

      FD(N)=FD(NM)

      FNUM(N)=FNUM(NM)

      NM=NM-1
      N=N-1

      RETURN
      END



      SUBROUTINE SAMPI

      IMPLICIT NONE
      include 'variables'
      INTEGER NNC
      
      NSMP_C=0
      DO NNC=1,MNC
        TEMPS(NNC)=0.D0
        DENSS(NNC)=0.D0
        VELOS(1,NNC)=0.D0
        VELOS(2,NNC)=0.D0
        VELOS(3,NNC)=0.D0
      END DO
      NSMP_F=0
      FORCE(1)=0.D0
      FORCE(2)=0.D0
      FORCE(3)=0.D0

      RETURN
      END



      SUBROUTINE SAMPLE

      IMPLICIT NONE
      include 'variables'
      INTEGER NNC

      NSMP_C=NSMP_C+1
      DO NNC=1,MNC
        TEMPS(NNC)=TEMPS(NNC)+TEMP(NNC)
        DENSS(NNC)=DENSS(NNC)+DENS(NNC)
        VELOS(1,NNC)=VELOS(1,NNC)+VELO(1,NNC)
        VELOS(2,NNC)=VELOS(2,NNC)+VELO(2,NNC)
        VELOS(3,NNC)=VELOS(3,NNC)+VELO(3,NNC)
      END DO

      RETURN
      END



      SUBROUTINE OUT

      IMPLICIT NONE
      include 'variables'
      INTEGER MX,MY,MZ,NNC
      DOUBLE PRECISION VX1,VY1,VZ1,DCELL1,TCELL1,PCELL1
      DOUBLE PRECISION X,Y,Z,VX,VY,VZ,DCELL,TCELL,PCELL,VMP
      COMMON/FILE_NAME/FILENAME(1)
      CHARACTER *40 FILENAME

      VMP=DSQRT(2.D0*BOLTZ*TEMP0/SP(5))

      WRITE(FILENAME(1),20000) NDT
20000 FORMAT('DSBGK_NDT',I6,'.DAT')
      OPEN(17,FILE=FILENAME(1))
      WRITE(17,*) 'TITLE="DSBGK TRANSIENT RESULT"'
      WRITE(17,*) 'VARIABLES="X","Y","u","v","w","n","T","OBST"'
      WRITE(17,*) 'ZONE T="DSBGK",I=',NCX,',J=',NCY,',F=POINT'
      MZ=1
      DO MY=1,NCY
        DO MX=1,NCX
          NNC=(MZ-1)*NCX*NCY+(MY-1)*NCX+MX
          X=(CG(1,NNC)+0.5D0*(CG(2,NNC)-CG(1,NNC)))/CB(6)
          Y=(CG(3,NNC)+0.5D0*(CG(4,NNC)-CG(3,NNC)))/CB(6)
          VX=VELO(1,NNC)/VMP
          VY=VELO(2,NNC)/VMP
          VZ=VELO(3,NNC)/VMP
          DCELL=DENS(NNC)/DENS0
          TCELL=TEMP(NNC)/TEMP0
          PCELL=DCELL*TCELL
          WRITE(17,998) X,Y,VX,VY,VZ,DCELL-1,TCELL-1,OBST(NNC)
        END DO
      END DO
      CLOSE(17)
998   FORMAT(' ',7E15.6,I2)
999   FORMAT(' ',8E15.6,I2)

      IF (NDT.EQ.NDT_T) THEN
        WRITE(FILENAME(1),10000) NSMP_C
10000   FORMAT('DSBGK_3D_NSMP_C',I6,'.DAT')
        OPEN(17,FILE=FILENAME(1))
        WRITE(17,*) 'TITLE="DSBGK TIME AVERAGE RESULT"'
        WRITE(17,*) 'VARIABLES="X","Y","Z","u","v","w","n","T","OBST"'
        WRITE(17,*) 'ZONE T="3D",I=',NCX,',J=',NCY,',K=',NCZ,',F=POINT'
        DO MZ=1,NCZ
        DO MY=1,NCY
          DO MX=1,NCX
            NNC=(MZ-1)*NCX*NCY+(MY-1)*NCX+MX
            X=(CG(1,NNC)+0.5D0*(CG(2,NNC)-CG(1,NNC)))/CB(6)
            Y=(CG(3,NNC)+0.5D0*(CG(4,NNC)-CG(3,NNC)))/CB(6)
            Z=(CG(5,NNC)+0.5D0*(CG(6,NNC)-CG(5,NNC)))/CB(6)
            VX=VELOS(1,NNC)/NSMP_C/VMP
            VY=VELOS(2,NNC)/NSMP_C/VMP
            VZ=VELOS(3,NNC)/NSMP_C/VMP
            DCELL=DENSS(NNC)/NSMP_C/DENS0
            TCELL=TEMPS(NNC)/NSMP_C/TEMP0
            PCELL=DCELL*TCELL

            WRITE(17,999) X,Y,Z,VX,VY,VZ,DCELL-1,TCELL-1,OBST(NNC)

          END DO
        END DO
        END DO
        CLOSE(17)

        WRITE(FILENAME(1),30000) NSMP_C
30000   FORMAT('DSBGK_NSMP_C',I6,'.DAT')
        OPEN(17,FILE=FILENAME(1))
        WRITE(17,*) 'TITLE="DSBGK TIME AVERAGE RESULT"'
        WRITE(17,*) 'VARIABLES="X","Y","u","v","w","n","T","OBST"'
        WRITE(17,*) 'ZONE T="DSBGK",I=',NCX,',J=',NCY,',F=POINT'
        DO MY=1,NCY
          DO MX=1,NCX

            MZ=1 !layer 1
            NNC=(MZ-1)*NCX*NCY+(MY-1)*NCX+MX
            VX1=VELOS(1,NNC)/NSMP_C/VMP
            VY1=VELOS(2,NNC)/NSMP_C/VMP
            VZ1=VELOS(3,NNC)/NSMP_C/VMP
            DCELL1=DENSS(NNC)/NSMP_C/DENS0
            TCELL1=TEMPS(NNC)/NSMP_C/TEMP0
            PCELL1=DCELL1*TCELL1

            MZ=1 !could be a different layer for average
            NNC=(MZ-1)*NCX*NCY+(MY-1)*NCX+MX
            X=(CG(1,NNC)+0.5D0*(CG(2,NNC)-CG(1,NNC)))/CB(6)
            Y=(CG(3,NNC)+0.5D0*(CG(4,NNC)-CG(3,NNC)))/CB(6)
            VX=VELOS(1,NNC)/NSMP_C/VMP
            VY=VELOS(2,NNC)/NSMP_C/VMP
            VZ=VELOS(3,NNC)/NSMP_C/VMP
            DCELL=DENSS(NNC)/NSMP_C/DENS0
            TCELL=TEMPS(NNC)/NSMP_C/TEMP0
            PCELL=DCELL*TCELL

            VX=(VX+VX1)/2
            VY=(VY+VY1)/2
            VZ=(VZ+VZ1)/2
            DCELL=(DCELL+DCELL1)/2
            TCELL=(TCELL+TCELL1)/2
            PCELL=(PCELL+PCELL1)/2

            WRITE(17,998) X,Y,VX,VY,VZ,DCELL-1,TCELL-1,OBST(NNC)
          END DO
        END DO
        CLOSE(17)
      END IF

      RETURN
      END



      SUBROUTINE SUM
C    also can be applied to momentum summation for checking

      IMPLICIT NONE
      include 'variables'
      INTEGER N,NNC
      DOUBLE PRECISION VMP

      SUM_NM=0.D0
      SUM_CELL=0.D0
      FLUX_MEAN(1)=0.D0
      FLUX_MEAN(2)=0.D0
      FLUX_MEAN(3)=0.D0

      DO N=1,NM
        SUM_NM=SUM_NM+FNUM(N)
      END DO
      DO NNC=1,MNC
        IF (OBST(NNC).EQ.0) THEN
          SUM_CELL=SUM_CELL+DENS(NNC)*(CG(2,NNC)-CG(1,NNC))*
     &             (CG(4,NNC)-CG(3,NNC))*(CG(6,NNC)-CG(5,NNC))

          FLUX_MEAN(1)=FLUX_MEAN(1)+VELO(1,NNC)*DENS(NNC)
          FLUX_MEAN(2)=FLUX_MEAN(2)+VELO(2,NNC)*DENS(NNC)
          FLUX_MEAN(3)=FLUX_MEAN(3)+VELO(3,NNC)*DENS(NNC)
        END IF
      END DO

      VMP=DSQRT(2.D0*BOLTZ*TEMP0/SP(5))
      FLUX_MEAN(1)=FLUX_MEAN(1)/(MNC*DENS0*VMP)
      FLUX_MEAN(2)=FLUX_MEAN(2)/(MNC*DENS0*VMP)
      FLUX_MEAN(3)=FLUX_MEAN(3)/(MNC*DENS0*VMP)

      WRITE(*,*) 'SUM OVER MOLECULES AND CELLS TO CHECK CONVERGE:'
      WRITE(*,*) 'MOLE & CELL', SUM_NM/SUM_INIT, SUM_CELL/SUM_INIT

      RETURN
      END



      SUBROUTINE UPDATE

      IMPLICIT NONE
      include 'variables'
      INTEGER NNC
      DOUBLE PRECISION SUMN0,SUMU0,SUMV0,SUMW0,VV0,SUME0,VOLUM
      DOUBLE PRECISION SUMN1,SUMU1,SUMV1,SUMW1,VV1,SUME1
      DOUBLE PRECISION DENS_OLD,TEMP_OLD,VELO_OLD(3),AAA

      IF (NDT.LE.NDT_S) THEN !be careful if change the convergence speed by AAA
        AAA=1.D0-0.8D0*(NDT-1.D0)/NDT_S
      ELSE 
        AAA=1.D0-0.8D0
      END IF

      DO NNC=1,MNC
C       DENS_OLD=DENS(NNC)
C       TEMP_OLD=TEMP(NNC)
C       VELO_OLD(1)=VELO(1,NNC)
C       VELO_OLD(2)=VELO(2,NNC)
C       VELO_OLD(3)=VELO(3,NNC)

        VOLUM=(CG(2,NNC)-CG(1,NNC))*
     &        (CG(4,NNC)-CG(3,NNC))*(CG(6,NNC)-CG(5,NNC))

        SUMN0=DENS(NNC)*VOLUM
        SUMU0=SUMN0*VELO(1,NNC)
        SUMV0=SUMN0*VELO(2,NNC)
        SUMW0=SUMN0*VELO(3,NNC)
        VV0=VELO(1,NNC)**2.D0+VELO(2,NNC)**2.D0+VELO(3,NNC)**2.D0
        SUME0=SUMN0*(1.5D0*BOLTZ*TEMP(NNC)+0.5D0*SP(5)*VV0)

        SUMN1=SUMN0+DENSA(NNC)
        SUMU1=SUMU0+VELOA(1,NNC)
        SUMV1=SUMV0+VELOA(2,NNC)
        SUMW1=SUMW0+VELOA(3,NNC)
        SUME1=SUME0+TEMPA(NNC)*0.5D0*SP(5)
  
        DENS(NNC)=SUMN1/VOLUM
        VELO(1,NNC)=SUMU1/SUMN1
        VELO(2,NNC)=SUMV1/SUMN1
        VELO(3,NNC)=SUMW1/SUMN1
        VV1=VELO(1,NNC)**2.D0+VELO(2,NNC)**2.D0+VELO(3,NNC)**2.D0
        TEMP(NNC)=(SUME1/SUMN1-0.5D0*SP(5)*VV1)/1.5D0/BOLTZ
  
C       DENS(NNC)=(1.D0-AAA)*DENS_OLD+AAA*DENS(NNC)
C       TEMP(NNC)=(1.D0-AAA)*TEMP_OLD+AAA*TEMP(NNC)
C       VELO(1,NNC)=(1.D0-AAA)*VELO_OLD(1)+AAA*VELO(1,NNC)
C       VELO(2,NNC)=(1.D0-AAA)*VELO_OLD(2)+AAA*VELO(2,NNC)
C       VELO(3,NNC)=(1.D0-AAA)*VELO_OLD(3)+AAA*VELO(3,NNC)

        TEMPA(NNC)=0.D0
        DENSA(NNC)=0.D0
        VELOA(1,NNC)=0.D0
        VELOA(2,NNC)=0.D0
        VELOA(3,NNC)=0.D0
      END DO

      CALL WALLD

      RETURN
      END



      SUBROUTINE WALLD 
*----IB=1 stream (DB.GE.0, vacuum), 2 solid wall, 3 symmetry, 0 pseudo in 2D case, -1 periodic, 4 cell interface  
C    For IB=2, we preset VB(0:3) and TB; usually VB(0)=1 for active
C    Here, only deal with IB=2 of solid wall to determine the effective number density

      IMPLICIT NONE
      include 'variables'
      INTEGER MX,MY,MZ,NNC,I
      DOUBLE PRECISION DCELL,TCELL,VMP,SCOS,FLUXIN,ERF,DX,DY,DZ

      IF (0.EQ.0) THEN !smooth/stable but has small errors at high Kn > 1
C    The mass flux is computed using the density, velocity and temperature of adjacent cell 
C    in an equilibrium distribution, can be imporved using additional non-equilibrium terms;
C    also can use n, \vec u, T extrapolated at surface location for possible improvement. 

        DO MZ=1,NCZ
          DO MY=1,NCY
            DO MX=1,NCX  
              NNC=(MZ-1)*NCY*NCX+(MY-1)*NCX+MX
              DCELL=DENS(NNC)
              TCELL=TEMP(NNC)
              VMP=DSQRT(2.D0*BOLTZ*TCELL/SP(5))
  
              I=1
              IF (IB(I,NNC).EQ.2) THEN
                SCOS=0.D0-(VELO(1,NNC)-VELOB(1,I,NNC))/VMP
                FLUXIN=EXP(0.D0-SCOS**2.D0)+SPI*SCOS*(1.D0+ERF(SCOS))
                DENSB(I,NNC)=DCELL*DSQRT(TCELL/TEMPB(I,NNC))*FLUXIN
              END IF
  
              I=3
              IF (IB(I,NNC).EQ.2) THEN
                SCOS=0.D0-(VELO(2,NNC)-VELOB(2,I,NNC))/VMP
                FLUXIN=EXP(0.D0-SCOS**2.D0)+SPI*SCOS*(1.D0+ERF(SCOS))
                DENSB(I,NNC)=DCELL*DSQRT(TCELL/TEMPB(I,NNC))*FLUXIN
              END IF

              I=5
              IF (IB(I,NNC).EQ.2) THEN
                SCOS=0.D0-(VELO(3,NNC)-VELOB(3,I,NNC))/VMP
                FLUXIN=EXP(0.D0-SCOS**2.D0)+SPI*SCOS*(1.D0+ERF(SCOS))
                DENSB(I,NNC)=DCELL*DSQRT(TCELL/TEMPB(I,NNC))*FLUXIN
              END IF

              I=2
              IF (IB(I,NNC).EQ.2) THEN
                SCOS=(VELO(1,NNC)-VELOB(1,I,NNC))/VMP
                FLUXIN=EXP(0.D0-SCOS**2.D0)+SPI*SCOS*(1.D0+ERF(SCOS))
                DENSB(I,NNC)=DCELL*DSQRT(TCELL/TEMPB(I,NNC))*FLUXIN
              END IF

              I=4
              IF (IB(I,NNC).EQ.2) THEN
                SCOS=(VELO(2,NNC)-VELOB(2,I,NNC))/VMP
                FLUXIN=EXP(0.D0-SCOS**2.D0)+SPI*SCOS*(1.D0+ERF(SCOS))
                DENSB(I,NNC)=DCELL*DSQRT(TCELL/TEMPB(I,NNC))*FLUXIN
              END IF

              I=6
              IF (IB(I,NNC).EQ.2) THEN
                SCOS=(VELO(3,NNC)-VELOB(3,I,NNC))/VMP
                FLUXIN=EXP(0.D0-SCOS**2.D0)+SPI*SCOS*(1.D0+ERF(SCOS))
                DENSB(I,NNC)=DCELL*DSQRT(TCELL/TEMPB(I,NNC))*FLUXIN
              END IF
  
            END DO
          END DO
        END DO
      ELSE !statistically accurate but noisy 
        DO MZ=1,NCZ
          DO MY=1,NCY
            DO MX=1,NCX  
              NNC=(MZ-1)*NCY*NCX+(MY-1)*NCX+MX
              DX=CG(2,NNC)-CG(1,NNC)
              DY=CG(4,NNC)-CG(3,NNC)
              DZ=CG(6,NNC)-CG(5,NNC)
  
              I=1
              IF (IB(I,NNC).EQ.2) THEN
                FLUXIN=DSQRT(2.D0*PI*SP(5)/BOLTZ/TEMPB(I,NNC))
                DENSB(I,NNC)=FLUXIN*NIN(I,NNC)/DT/DY/DZ
                NIN(I,NNC)=0.D0
              END IF
  
              I=2
              IF (IB(I,NNC).EQ.2) THEN
                FLUXIN=DSQRT(2.D0*PI*SP(5)/BOLTZ/TEMPB(I,NNC))
                DENSB(I,NNC)=FLUXIN*NIN(I,NNC)/DT/DY/DZ
                NIN(I,NNC)=0.D0
              END IF

              I=3
              IF (IB(I,NNC).EQ.2) THEN
                FLUXIN=DSQRT(2.D0*PI*SP(5)/BOLTZ/TEMPB(I,NNC))
                DENSB(I,NNC)=FLUXIN*NIN(I,NNC)/DT/DX/DZ
                NIN(I,NNC)=0.D0
              END IF

              I=4
              IF (IB(I,NNC).EQ.2) THEN
                FLUXIN=DSQRT(2.D0*PI*SP(5)/BOLTZ/TEMPB(I,NNC))
                DENSB(I,NNC)=FLUXIN*NIN(I,NNC)/DT/DX/DZ
                NIN(I,NNC)=0.D0
              END IF

              I=5
              IF (IB(I,NNC).EQ.2) THEN
                FLUXIN=DSQRT(2.D0*PI*SP(5)/BOLTZ/TEMPB(I,NNC))
                DENSB(I,NNC)=FLUXIN*NIN(I,NNC)/DT/DX/DY
                NIN(I,NNC)=0.D0
              END IF

              I=6
              IF (IB(I,NNC).EQ.2) THEN
                FLUXIN=DSQRT(2.D0*PI*SP(5)/BOLTZ/TEMPB(I,NNC))
                DENSB(I,NNC)=FLUXIN*NIN(I,NNC)/DT/DX/DY
                NIN(I,NNC)=0.D0
              END IF
  
            END DO
          END DO
        END DO
      END IF

      RETURN
      END



      SUBROUTINE BCWALL(U,FD1,NNC2,NNC1)
C    Update FD1 and U(3) (incoming velocity and then reflected value in global frame)
C    U(3) is in global frame; incoming U123 and reflected V123 are in surface frame
C    In surface frame, 2,3 for tangential and 1 for normal, incoming U1<0 and reflected V1>0
C    Selection of axis directions of the surface frame doesn't satisfy the right-hand rule 
C    but doesn't matter as the two tangential components are equivalent.  

      IMPLICIT NONE
      include 'variables'
      INTEGER NNC2,NNC1
      DOUBLE PRECISION U(3),FD1,TWALL,VMP,RF,V1,V2,V3
      DOUBLE PRECISION D,FEQ,U1,U2,U3,FCLL

      TWALL=TEMPB(-NNC2,NNC1)

      IF (IREFLE.EQ.1) THEN !Maxwell diffuse reflection
        VMP=DSQRT(2.D0*BOLTZ*TWALL/SP(5))
        V1=DSQRT(-LOG(RF(0)))*VMP
        CALL RVELC(V2,V3,VMP)

        IF (NNC2.EQ.-1) THEN
          U(1)=V1
          U(2)=V2
          U(3)=V3
        ELSE IF (NNC2.EQ.-2) THEN
          U(1)=-V1
          U(2)=V2
          U(3)=V3
        ELSE IF (NNC2.EQ.-3) THEN
          U(1)=V2
          U(2)=V1
          U(3)=V3
        ELSE IF (NNC2.EQ.-4) THEN
          U(1)=V2
          U(2)=-V1
          U(3)=V3
        ELSE IF (NNC2.EQ.-5) THEN
          U(1)=V3
          U(2)=V2
          U(3)=V1
        ELSE IF (NNC2.EQ.-6) THEN
          U(1)=V3
          U(2)=V2
          U(3)=-V1
        END IF
      
        D=DENSB(-NNC2,NNC1)
        CALL FEQGET(BOLTZ,PI,SP(5),D,TWALL,U(1),U(2),U(3),FEQ)
C    Generally, we use the velocity \vec c_r in surface frame as FD1=f_B(\vec c_r)
C    For Maxwell reflection model, we can use the components of \vec c_r in the global or surface frames. 
C    Here, we use the components of \vec c_r in the global frame. 
        FD1=FEQ

        U(1)=U(1)+VELOB(1,-NNC2,NNC1)
        U(2)=U(2)+VELOB(2,-NNC2,NNC1)
        U(3)=U(3)+VELOB(3,-NNC2,NNC1)
      ELSE                            !CLL reflection model
        U(1)=U(1)-VELOB(1,-NNC2,NNC1) !need relative velocity
        U(2)=U(2)-VELOB(2,-NNC2,NNC1)
        U(3)=U(3)-VELOB(3,-NNC2,NNC1)
  
        IF (NNC2.EQ.-1) THEN
          U1=U(1)
          U2=U(2)
          U3=U(3)
        ELSE IF (NNC2.EQ.-2) THEN
          U1=-U(1)
          U2=U(2)
          U3=U(3)
        ELSE IF (NNC2.EQ.-3) THEN
          U1=U(2)
          U2=U(1)
          U3=U(3)
        ELSE IF (NNC2.EQ.-4) THEN
          U1=-U(2)
          U2=U(1)
          U3=U(3)
        ELSE IF (NNC2.EQ.-5) THEN
          U1=U(3)
          U2=U(2)
          U3=U(1)
        ELSE IF (NNC2.EQ.-6) THEN
          U1=-U(3)
          U2=U(2)
          U3=U(1)
        END IF

        CALL REFLECT(BOLTZ,PI,TWALL,SP(5),ALPHAN,ALPHAT,
     &               U1,U2,U3,V1,V2,V3)
        D=DENSB(-NNC2,NNC1)
        CALL FCLLGET(BOLTZ,PI,TWALL,SP(5),ALPHAN,ALPHAT,
     &               U1,U2,U3,V1,V2,V3,D,FCLL) !?tentative algorithm
C    Here, incoming and reflected relative velocity are in surface frame rather than global frame
        FD1=FCLL

        IF (NNC2.EQ.-1) THEN
          U(1)=V1
          U(2)=V2
          U(3)=V3
        ELSE IF (NNC2.EQ.-2) THEN
          U(1)=-V1
          U(2)=V2
          U(3)=V3
        ELSE IF (NNC2.EQ.-3) THEN
          U(1)=V2
          U(2)=V1
          U(3)=V3
        ELSE IF (NNC2.EQ.-4) THEN
          U(1)=V2
          U(2)=-V1
          U(3)=V3
        ELSE IF (NNC2.EQ.-5) THEN
          U(1)=V3
          U(2)=V2
          U(3)=V1
        ELSE IF (NNC2.EQ.-6) THEN
          U(1)=V3
          U(2)=V2
          U(3)=-V1
        END IF

        U(1)=U(1)+VELOB(1,-NNC2,NNC1)
        U(2)=U(2)+VELOB(2,-NNC2,NNC1)
        U(3)=U(3)+VELOB(3,-NNC2,NNC1)
      END IF

      RETURN
      END



      SUBROUTINE F2GET(FD1,U,DTUSE,NNC1,FD2) 
C    (1) time_relax=DENS*BOLTZ*TEMP/VISCO_MU
C    (2) time_relax=(5*BOLTZ/2/MASS)*DENS*BOLTZ*TEMP/HEAT_K if we model monatomic molecules

      IMPLICIT NONE
      include 'variables'
      INTEGER NNC1
      DOUBLE PRECISION FD1,U(3),DTUSE,FD2,TCELL,DCELL,TIME_RELAX,FEQ

      TCELL=TEMP(NNC1)
      DCELL=DENS(NNC1)

      IF (IRELAX.EQ.1) THEN
        TIME_RELAX=DCELL*BOLTZ*TCELL/VISCO
      ELSE
        TIME_RELAX=DCELL*BOLTZ*TCELL*2.D0/VISCO/3.D0 !monatomic, useful for thermal transpiration 
      END IF

      FEQ=EXP(-(SP(5)/BOLTZ/2.D0/TCELL)*((U(1)-VELO(1,NNC1))**2.D0+
     &    (U(2)-VELO(2,NNC1))**2.D0+(U(3)-VELO(3,NNC1))**2.D0))*
     &    DCELL*((SP(5)/BOLTZ/2.D0/TCELL/PI)**1.5D0)

      FD2=FEQ+(FD1-FEQ)*EXP(0.D0-DTUSE*TIME_RELAX)
      
      RETURN
      END



      SUBROUTINE FCLLGET(BOLTZ,PI,TWALL,SP5,ALPHAN,ALPHAT,
     &                   U1,U2,U3,V1,V2,V3,D,FCLL)
C    Here, the tentative formula of f_CLL is not rigorous and valid only when alpha is close to 1 (see the arXiv article)
C    In surface frame, 2,3 for tangential and 1 for normal, incoming U1<0 and reflected V1>0

      IMPLICIT NONE
      INTEGER LOOP,L
      DOUBLE PRECISION BOLTZ,PI,TWALL,SP5,ALPHAN,ALPHAT,U1,U2,U3,V1,V2
      DOUBLE PRECISION V3,D,FCLL,VMP,U1_N,U2_N,U3_N,V1_N,V2_N,V3_N
      DOUBLE PRECISION FACT,QUAD,Q

      VMP=DSQRT(2.D0*BOLTZ*TWALL/SP5)
      U1_N=U1/VMP
      U2_N=U2/VMP
      U3_N=U3/VMP
      V1_N=V1/VMP
      V2_N=V2/VMP
      V3_N=V3/VMP

CCCCC FCLL=D/(VMP**3.D0)/(PI**2.5D0)/2.D0/ALPHAT/ALPHAN* !old version
      FCLL=(D/(VMP**3.D0)/DSQRT(PI)/2.D0                 !=a as in arXiv paper
     &     /PI/PI/ALPHAT/ALPHAN)*
     &     EXP(-((V2_N-DSQRT(1.D0-ALPHAT)*U2_N)**2.D0)/ALPHAT)*
     &     EXP(-((V3_N-DSQRT(1.D0-ALPHAT)*U3_N)**2.D0)/ALPHAT)*
     &     EXP(-(V1_N**2.D0+(1.D0-ALPHAN)*U1_N**2.D0)/ALPHAN)

      FACT=2.D0*DSQRT(1.D0-ALPHAN)*V1_N*ABS(U1_N)/ALPHAN
      QUAD=0.D0
      LOOP=100                                           !CHECK if large enough
      DO L=1,LOOP
        Q=(L-0.5D0)*2.D0*PI/LOOP
        QUAD=QUAD+EXP(FACT*COS(Q))
      END DO
      QUAD=QUAD*2.D0*PI/LOOP

      FCLL=FCLL*QUAD

      RETURN
      END



      SUBROUTINE FEQGET(BOLTZ,PI,SP5,D,T,UX,UY,UZ,FEQ) 

      IMPLICIT NONE
      DOUBLE PRECISION BOLTZ,PI,SP5,D,T,UX,UY,UZ,FEQ

      FEQ=D*EXP(-SP5*(UX*UX+UY*UY+UZ*UZ)/2.D0/BOLTZ/T)/     
     &    ((2.D0*PI*BOLTZ*T/SP5)**1.5D0)

      RETURN
      END



      SUBROUTINE PPGET(DTR,U,X,MX,MY,MZ,NNC1,NNC2,DTUSE)

      IMPLICIT NONE
      include 'variables'
      INTEGER MX,MY,MZ,NNC1,NNC2,JJ,I
      DOUBLE PRECISION DTR,U(3),X(3),DTUSE,DTNEED(6)

      IF (IB(1,NNC1).NE.0) THEN
        IF (U(1).LT.0.D0) THEN
          DTNEED(1)=(CG(1,NNC1)-X(1))/U(1)
          IF (DTNEED(1).LT.0.D0) THEN
            WRITE(*,*) ' ERROR IN PPGET: 1'
            PAUSE
          END IF
          DTNEED(2)=-1.D0
        ELSE
          DTNEED(1)=-1.D0
          DTNEED(2)=(CG(2,NNC1)-X(1))/U(1)
          IF (DTNEED(2).LT.0.D0) THEN
            WRITE(*,*) ' ERROR IN PPGET: 2'
            PAUSE
          END IF
        END IF
      ELSE
        DTNEED(1)=-1.D0
        DTNEED(2)=-1.D0
      END IF

      IF (IB(3,NNC1).NE.0) THEN
        IF (U(2).LT.0.D0) THEN
          DTNEED(3)=(CG(3,NNC1)-X(2))/U(2)
          IF (DTNEED(3).LT.0.D0) THEN
            WRITE(*,*) ' ERROR IN PPGET: 3'
            PAUSE
          END IF
          DTNEED(4)=-1.D0
        ELSE
          DTNEED(3)=-1.D0
          DTNEED(4)=(CG(4,NNC1)-X(2))/U(2)
          IF (DTNEED(4).LT.0.D0) THEN
            WRITE(*,*) ' ERROR IN PPGET: 4'
            PAUSE
          END IF 
        END IF
      ELSE
        DTNEED(3)=-1.D0
        DTNEED(4)=-1.D0
      END IF

      IF (IB(5,NNC1).NE.0) THEN
        IF (U(3).LT.0.D0) THEN
          DTNEED(5)=(CG(5,NNC1)-X(3))/U(3)
          IF (DTNEED(5).LT.0.D0) THEN
            WRITE(*,*) ' ERROR IN PPGET: 5'
            PAUSE
          END IF
          DTNEED(6)=-1.D0
        ELSE
          DTNEED(5)=-1.D0
          DTNEED(6)=(CG(6,NNC1)-X(3))/U(3)
          IF (DTNEED(6).LT.0.D0) THEN
            WRITE(*,*) ' ERROR IN PPGET: 6'
            PAUSE
          END IF 
        END IF
      ELSE
        DTNEED(5)=-1.D0
        DTNEED(6)=-1.D0
      END IF

      JJ=0
      DTUSE=DTR
      DO I=1,6
        IF ((DTNEED(I).GT.0).AND.(DTNEED(I).LE.DTUSE)) THEN
          DTUSE=DTNEED(I)
          JJ=I
        END IF
      END DO

*----IB=1 stream (DB.GE.0, vacuum), 2 solid wall, 3 symmetry, 0 pseudo in 2D case, -1 periodic, 4 cell interface  

      IF (JJ.EQ.0) THEN
        NNC2=0
      ELSE IF (IB(JJ,NNC1).EQ.4) THEN !cell-cell interface always inside domain
        IF (JJ.EQ.1) THEN
          MX=MX-1
        ELSE IF (JJ.EQ.2) THEN
          MX=MX+1
        ELSE IF (JJ.EQ.3) THEN
          MY=MY-1
        ELSE IF (JJ.EQ.4) THEN
          MY=MY+1
        ELSE IF (JJ.EQ.5) THEN
          MZ=MZ-1
        ELSE IF (JJ.EQ.6) THEN
          MZ=MZ+1
        END IF
        NNC2=(MZ-1)*NCX*NCY+(MY-1)*NCX+MX
      ELSE !IB=1,2,3,-1, never be 0
        NNC2=-JJ
      END IF

      IF (IB(1,NNC1).NE.0) THEN
        X(1)=X(1)+DTUSE*U(1)
      END IF

      IF (IB(3,NNC1).NE.0) THEN
        X(2)=X(2)+DTUSE*U(2)
      END IF

      IF (IB(5,NNC1).NE.0) THEN
        X(3)=X(3)+DTUSE*U(3)
      END IF

      RETURN
      END



      SUBROUTINE RVELC(U,V,VMP)
*----generates two random velocity components U an V in an equilibrium
*----gas with most probable speed VMP  (based on eqns (C10) and (C12))

      IMPLICIT NONE
      DOUBLE PRECISION U,V,VMP,A,B,RF

      A=DSQRT(-LOG(RF(0)))
      B=6.283185308*RF(0)
      U=A*SIN(B)*VMP
      V=A*COS(B)*VMP

      RETURN
      END


      FUNCTION GAM(X)

      IMPLICIT NONE
      DOUBLE PRECISION X,A,Y,GAM

      A=1.
      Y=X
      IF (Y.LT.1.) THEN
        A=A/Y
      ELSE
40      Y=Y-1
        IF (Y.GE.1.) THEN
          A=A*Y
          GOTO 40
        END IF
      END IF
      GAM=A*(1.-0.5748646*Y+0.9512363*Y**2-0.6998588*Y**3+
     &    0.4245549*Y**4-0.1010678*Y**5)
      RETURN
      END


      FUNCTION RF(IDUM) !use IMPLICIT NONE and DOUBLE PRECISION RF
*----generates a uniformly distributed random fraction between 0 and 1
*----IDUM will generally be 0, but can be negative to re-initialize the seed
      DOUBLE PRECISION RF

      SAVE MA,INEXT,INEXTP
      PARAMETER (MBIG=1000000000,MSEED=161803398,MZ=0,FAC=1.E-9)
      DIMENSION MA(55)
      DATA IFF/0/
      IF (IDUM.LT.0.OR.IFF.EQ.0) THEN
        IFF=1
        MJ=MSEED-IABS(IDUM)
        MJ=MOD(MJ,MBIG)
        MA(55)=MJ
        MK=1
        DO 50 I=1,54
          II=MOD(21*I,55)
          MA(II)=MK
          MK=MJ-MK
          IF (MK.LT.MZ) MK=MK+MBIG
          MJ=MA(II)
50      CONTINUE
        DO 100 K=1,4
          DO 60 I=1,55
            MA(I)=MA(I)-MA(1+MOD(I+30,55))
            IF (MA(I).LT.MZ) MA(I)=MA(I)+MBIG
60        CONTINUE
100     CONTINUE
        INEXT=0
        INEXTP=31
      END IF
200   INEXT=INEXT+1
      IF (INEXT.EQ.56) INEXT=1
      INEXTP=INEXTP+1
      IF (INEXTP.EQ.56) INEXTP=1
      MJ=MA(INEXT)-MA(INEXTP)
      IF (MJ.LT.MZ) MJ=MJ+MBIG
      MA(INEXT)=MJ
      RF=MJ*FAC
      IF (RF.GT.1.E-8.AND.RF.LT.0.99999999) RETURN
      GOTO 200
      END


      FUNCTION ERF(S)

      IMPLICIT NONE
      DOUBLE PRECISION S,B,D,C,T,ERF

      B=ABS(S)
      IF (B.GT.4.) THEN
        D=1.
      ELSE
        C=EXP(-B*B)
        T=1./(1.+0.3275911*B)
        D=1.-(0.254829592*T-0.284496736*T*T+1.421413741*T*T*T-
     &    1.453152027*T*T*T*T+1.061405429*T*T*T*T*T)*C
      END IF
      IF (S.LT.0.) D=-D
      ERF=D

      RETURN
      END



      SUBROUTINE PARAMETERS

      IMPLICIT NONE
      include 'variables'
      INTEGER MX,MY,MZ,NNC,I,POROSITY
      DOUBLE PRECISION DENSM,PRES0,VMP,LAMBDA,KN,SIZE,TWALLGET

      PI=3.141592654
      SPI=DSQRT(PI)
      BOLTZ=1.380622E-23

C     AR MOLECULE
C     SP(1)=4.17D-10
C     SP(2)=273.D0
C     SP(3)=0.81D0
C     SP(4)=1.D0
      SP(5)=66.3D-27 !only need MOLECULAR MASS in solving BGK equation
      VISCO=2.117D-5 !VISCO is dynamic viscosity not kinematic viscosity
C    relaxation time = DENS * BOLTZ * TEMP / VISCO
C    relaxation time = 5 * BOLTZ * BOLTZ * DENS * TEMP / 2 / SP(5) / Heat conductivity

      IRELAX=1     !=1 matching dynamic viscosity and =2 for heat conductivity
      IREFLE=1     !=1 for Maxwell boundary (=2 for CLL but not ready)
      ALPHAN=1.0D0 !ENERGY ACCOMMODATION COEFFICIENT OF NORMAL COMPONENT 
      ALPHAT=1.0D0 !ENERGY ACCOMMODATION COEFFICIENT OF TANGENTIAL ...

      PRES0=101325.D0
      TEMP0=273D0       !INITIAL TEMPERATURE
      DENS0=PRES0/(BOLTZ*TEMP0)
      DRATIO=0.5D0      !outlet/inlet density ratio
      VELO0(1)=0.D0     !INITIAL VELOCITY
      VELO0(2)=0.D0
      VELO0(3)=0.D0

      VMP=DSQRT(2.D0*BOLTZ*TEMP0/SP(5))
      LAMBDA=VISCO*SPI*VMP/(2*PRES0) !Bird

      NDT_T=5000
*----NDT_T is the total time-step number in each sigle run
      NDT_D=1
*----NDT_D is the time-step number between two consecutive samples
      NDT_O=500
*----NDT_O is the time-step number between two consecutive OUT
      NDT_S=4000
*----NDT_S is the time-step before starting the sampling process

      IF (IREFLE.EQ.2) WRITE(*,10000) ALPHAN,ALPHAT
      WRITE(*,10001) IRELAX,IREFLE,NPT
      WRITE(*,10003) NDT_T,NDT_D,NDT_S
      WRITE(*,10004) NCX,NCY,NCZ
10000 FORMAT(' ALPHAN :',E14.6,' ALPHAT :',E14.6)
10001 FORMAT(' IRELAX :',I14  ,' IREFLE :',I14  ,' NPT     :',I14)
10003 FORMAT(' NDT_T  :',I14  ,' NDT_D  :',I14  ,' NDT_S   :',I14)
10004 FORMAT(' NCX    :',I14  ,' NCY    :',I14  ,' NCZ     :',I14)
      WRITE(*,*)
      WRITE(*,10005) VELO0(1),VELO0(2),VELO0(3)
      WRITE(*,10006) DENS0,TEMP0,PRES0
      WRITE(*,10007) SP(5),VISCO,LAMBDA
10005 FORMAT(' VELO(1):',E14.6,' VELO(2):',E14.6,' VELO(3) :',E14.6)
10006 FORMAT(' DENS0  :',E14.6,' TEMP0  :',E14.6,' PRES0   :',E14.6)
10007 FORMAT(' MASS   :',E14.6,' VISCO  :',E14.6,' LAMBDA  :',E14.6)
10008 FORMAT(' UWALL  :',E14.6,' KN     :',E14.6,' SIZE    :',E14.6)

      CB(1)=0.0D0
      CB(2)=1.1D-6
      CB(3)=0.0D0
      CB(4)=1.0D-7
      CB(5)=0.0D0
      CB(6)=1.0D-7

      DT=2.0*CB(2)/VMP/NCX  !can be large as it is divided into smaller segments in MOVE

      SUM_INIT=DENS0*(CB(2)-CB(1))*(CB(4)-CB(3))*(CB(6)-CB(5))
      FNUM0=SUM_INIT/(NCX*NCY*NCZ*10.D0) !depends on the molecule number per cell

*----IB=1 stream (DB.GE.0, vacuum), 2 solid wall, 3 symmetry, 0 pseudo in 2D case, -1 periodic, 4 cell interface  

*----set arbitrary voxelized geometry by OBST distribution
      POROSITY=0
      DO MZ=1,NCZ
        DO MY=1,NCY
          DO MX=1,NCX
            NNC=(MZ-1)*NCX*NCY+(MY-1)*NCX+MX
            OBST(NNC)=0       !void cells

            IF (MX.GT.50.AND.MX.LE.NCX-50) THEN
              IF (MY.LE.40.OR.MY.GT.NCY-40) THEN
                OBST(NNC)=1   !solid cells
              END IF
            END IF

            IF (OBST(NNC).EQ.0) POROSITY=POROSITY+1
          END DO
        END DO
      END DO
      WRITE(*,*) 'Void cells number', POROSITY

      UWALL=0.D0        !speed of possible moving walls
      SIZE=(CB(4)-CB(3))*(NCY-40-40)/NCY
      KN=LAMBDA/SIZE    !use actual channel size < tank size = CB(4)-CB(3)
      WRITE(*,10008) UWALL,KN,SIZE
      WRITE(*,*) 

      DO MZ=1,NCZ
        DO MY=1,NCY
          DO MX=1,NCX
            NNC=(MZ-1)*NCX*NCY+(MY-1)*NCX+MX
      
            DO I=1,6
              IB(I,NNC)=4               !FIRST SET ALL SIDES AS cell-cell INTERFACE 
            END DO

            IF (MX.EQ.1)   IB(1,NNC)=1  !STREAM at left  inlet
            IF (MX.EQ.NCX) IB(2,NNC)=1  !STREAM at right outlet

            IF (MY.EQ.1)   IB(3,NNC)=1  !STREAM at lateral side of big tank
            IF (MY.EQ.NCY) IB(4,NNC)=1  !STREAM at lateral side of big tank

            IF (MZ.EQ.1)   IB(5,NNC)=-1 !periodic
            IF (MZ.EQ.NCZ) IB(6,NNC)=-1 !periodic

C           IB(5,NNC)=0                 !pseudo=0 BC sets all cells to 0 
C           IB(6,NNC)=0                 !periodic BC sets only boundary cells to -1

            IF (OBST(NNC).NE.0) THEN    !modify default IB by solid cells
              IB(:,NNC)=2               !solid walls for solid cells
            ELSE
              IF (MX.GT.1.AND.
     &            OBST(NNC-1).NE.0) IB(1,NNC)=2 !solid wall due to solid neighbor
              IF (MX.LT.NCX.AND.
     &            OBST(NNC+1).NE.0) IB(2,NNC)=2

              IF (MY.GT.1.AND.
     &            OBST(NNC-NCX).NE.0) IB(3,NNC)=2
              IF (MY.LT.NCY.AND.
     &            OBST(NNC+NCX).NE.0) IB(4,NNC)=2

              IF (MZ.GT.1.AND.
     &            OBST(NNC-NCX*NCY).NE.0) IB(5,NNC)=2
              IF (MZ.LT.NCZ.AND.
     &            OBST(NNC+NCX*NCY).NE.0) IB(6,NNC)=2
            END IF

            DO I=1, 6
              IF (IB(I,NNC).EQ.1) THEN !STREAM
                VELOB(0,I,NNC)=-1.D0   !+/- FOR ACTIVE/INACTIVE
C               VELOB(1,I,NNC)=0.D0    !specify density ratio below, cannot also specify VB 
C               VELOB(2,I,NNC)=0.D0
C               VELOB(3,I,NNC)=0.D0
                TEMPB(I,NNC)=TWALLGET(I,MX,MY,MZ)

                IF (MX.LT.NCX/2) THEN  !rough estimation for inlet area, not always valid 
                  DENSB(I,NNC)=DENS0
                ELSE
                  DENSB(I,NNC)=DENS0*DRATIO
                END IF
              ELSE IF (IB(I,NNC).EQ.2) THEN !SOLID WALL
                VELOB(0,I,NNC)=1.D0
                VELOB(1,I,NNC)=0.D0    !can add nonzero for each moving wall
                VELOB(2,I,NNC)=0.D0
                VELOB(3,I,NNC)=0.D0
                TEMPB(I,NNC)=TWALLGET(I,MX,MY,MZ)
              END IF
            END DO

          END DO
        END DO
      END DO
      
      RETURN
      END

C*****

      FUNCTION TWALLGET(I,MX,MY,MZ) 

      include 'variables'
      INTEGER I,MX,MY,MZ
      DOUBLE PRECISION TWALLGET

      IF (MX*2.LE.NCX) THEN
        TWALLGET=1.0D0*TEMP0  !can change with cell location for thermal transpiration
      ELSE
        TWALLGET=1.0D0*TEMP0  !can change with cell side I for more precise setting
      END IF

      RETURN
      END
