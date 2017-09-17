         PROGRAM ARC_CH4 

C-AD Modified on Jun, 2017 for mass independent fractionation of oxygen by
C-AD Alexander Dimoff. Deleted parts of original photochemical code that
C-AD are not needed to perform the box model of the solar nebula.
C
C       THIS PROGRAM IS A BOX MODEL OF THE SOLAR NEBULA. CHEMICAL SPECIES
C     ARE INJECTED INTO THE BOX, AND REACTIONS TAKE PLACE, WITH THE INTENT 
C     TO SEE WHETHER OR NOT THERE IS MASS INDEPENDT FRACTIONATION OF OXYGEN
C     IN THE REACTIONS.
C  
C     THE MIXING RATIO/NHS OF THE LONG-LIVED SPECIES  
C     ARE CALCULATED FROM THE EQUATION
C   
C         DF/DT  =  (1/N)*D/DZ(KN*DF/DZ + WNF) + P/N - LF 
C   
C     WHERE
C     F = MIXING RATIO (USOL)
C     N = TOTAL NUMBER DENSITY (DEN) = SUM(USOL(NQ))
C     L = CHEMICAL LOSS FREQUENCY (XL)
C     P = CHEMICAL PRODUCTION RATE (XP)
C 
C          THE SYSTEM OF PDES IS SOLVED USING THE REVERSE EULER      
C     METHOD.  LETTING THE TIME STEP GO TO INFINITY GIVES YOU NEWTONS
C     METHOD, E.G. IT REVERTS TO AN EFFICIENT STEADY-STATE SOLVER.  
C
C          THE LIST OF SUBROUTINES IS AS FOLLOWS:
C     (1) RATES  - DEFINES CHEMICAL REACTION RATES AND RAINOUT RATE
C
C     (2) OUTPUT - PRINTS OUT RESULTS
C     (3) DOCHEM - DOES CHEMISTRY FOR ALL SPECIES AT ALL GRID
C                  POINTS BY CALLING CHEMPL
C     (4) CHEMPL - COMPUTES CHEMICAL PRODUCTION AND LOSS RATES
C                  FOR ONE SPECIES AT ALL GRID POINTS
C     (5) LUDCMP - MATRIX SOLVING VIA DECMPOSITION, USED IN ACCORDANCE
C                  WITH THE FOLLOWING SUBROUTINE
C     (6) LUBKSB - MATRIX SOLVING VIA BACK SUBSTITUTION
C
C          OTHER DEFINED FUNCTIONS INCLUDE:
C     (1) TBDY   - COMPUTES 3-BODY REACTION RATES
C
C ***** REACTION LIST *****
C ***** FIND REACTIONS IN CHEM.DAT.CH4_OMIF *****
C     1)  SiO + O2 = SiO2 + O
C     2)  SiO + H2 = SiOH + H
C     3)  SiO + OH = SiO2 + H
C     4)  SiO + HO2 = SiO2 + OH
C     5)  SiO + O = SiO2
C     6)  SiO + H = SiOH
C     7)  SiO + O3 = SiO2 + O2
C     8)  O + O2 + M = O3 + M
C     9)  H2 + O = OH + H
C    10)  OH + O = H2O + O
C    11)  OH + OH = H2O + O
C    12)  H + O = OH
C    13)  H + O2 = OH + O
C    14)  O + HO2 = OH + O2
C    15)  O + H2O2 = OH + HO2
C    16)  OH + H2 = H2O + H
C    17)  OH + HO2 = H2O + O2
C    18)  OH + H2O2 = H2O + HO2
C    19)  HO2 + HO2 = H2O2 + O2
C    20)  OH + O3 = HO2 + O2
C    21)  HO2 + O3 = OH + 2O2
C    22)  MgO + SiO2 = MgSiO3
C    23)  Mg + O = MgO
C    24)  MgO + MgSiO3 = Mg2SiO4
C
C ***** THREE BODY RXN ***** ELIMINATE COMMON TERMS 
C    25)  H + O2 + H2O = HO2 + H2O
C    26)  H + O2 + H2 = HO2 + H2
C    27)  H + H + H2 = H2 + H2
C    28)  H + OH + H2O = H2O + H2O
C
C        THIS PROGRAM DOES THE CHEMISTRY AUTOMATICALLY.  THE CHEMICAL
C     REACTIONS ARE ENTERED ON DATA CARDS (.DAT) IN FIVE 10-DIGIT 
C     COLUMNS STARTING IN COLUMN 11, I.E.                         
C
C         REAC1     REAC2     PROD1     PROD2     PROD3
C
C     THE IMPORTANT PARAMETERS DESCRIBING THE CHEMISTRY ARE
C        NR   = NUMBER OF REACTIONS
C        NSP  = NUMBER OF CHEMICAL SPECIES
C        NQ   = NUMBER OF SPECIES FOR WHICH A DIFFUSION EQUATION
C               IS SOLVED AND WHICH ARE IN THE BIG, REVERSE EULER MATR
C        NMAX = MAXIMUM NUMBER OF REACTIONS IN WHICH AN INDIVIDUAL
C               SPECIES PARTICIPATES
C
C     THREE-BODY REACTIONS ARE WRITTEN
C     IN TWO-BODY FORM, SO THE DENSITY FACTOR MUST BE INCLUDED IN
C     THE RATE CONSTANT.
C
C        THE CHEMICAL REACTION SCHEME IS STORED IN THE FOLLOWING MATRICE
C
C     ISPEC(NSP) = VECTOR CONTAINING THE HOLLERITH NAMES OF THE
C                  CHEMICAL SPECIES.
C     JCHEM(5,NR) = MATRIX OF CHEMICAL REACTIONS.  THE FIRST TWO ARE
C                   REACTANTS, THE LAST THREE ARE PRODUCTS.
C     ILOSS(2,NSP,NMAX) = MATRIX OF LOSS PROCESSES.  ILOSS(1,I,L)
C                         HOLDS REACTION NUMBER J, ILOSS(2,I,L) HOLDS
C                         REACTANT NUMBER.
C     IPROD(NSP,NMAX) = MATRIX OF PRODUCTION PROCESSES.  IPROD(I,L)
C                       HOLDS REACTION NUMBER J.
C     NUML(NSP) = NUMBER OF NON-ZERO ELEMENTS FOR EACH ROW OF ILOSS
C     NUMP(NSP) = NUMBER OF NON-ZERO ELEMENTS FOR EACH ROW OF IPROD
C
      PARAMETER(NQ=14, NR=28, NSP=16, NMAX=10)
      DIMENSION FVAL(NQ),FV(NQ),DJAC(NQ,NQ),RHS(NQ),REL(NQ),IPVT(NQ)
     2  ,USAVE(NQ),USOL(NQ),ISPEC(NSP),TP(NQ),TL(NQ),D(NQ)
      DIMENSION RAT(NR)
      CHARACTER*30 CHEM(5,NR),PRODRX(NSP,NR),LOSSRX(NSP,NR)
      CHARACTER*4,DIRDATA
     
      COMMON/NBLOK/LO,LO2,LO3,LH,LH2,LOH,LH2O,LHO2,LH2O2,
     2  LSiO,LSiO2,LSiOH,LMg,LMgO,LMgSiO3,LMg2SiO4
      COMMON/RBLOK/A(NR),ILOSS(2,NSP,NMAX),IPROD(NSP,NMAX),
     2  JCHEM(5,NR),NUML(NSP),NUMP(NSP),T,DEN
C-AP *************************************************************

      DATA LO,LO2,LO3,LH,LH2,LOH,LH2O,LHO2,LH2O2,
     2  LSiO,LSiO2,LSiOH,LMg,LMgO,LMgSiO3,LMg2SiO4/
     3  1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16/
C
C ***** SPECIES DEFINITIONS *****
C   LONG-LIVED - ALL OUR SPECIES ARE "LONG-LIVED"
      ISPEC(1) = 1HO
      ISPEC(2) = 2HO2
      ISPEC(3) = 2HO3
      ISPEC(4) = 1HH
      ISPEC(5) = 2HH2
      ISPEC(6) = 2HOH
      ISPEC(7) = 3HH2O
      ISPEC(8) = 3HHO2
      ISPEC(9) = 4HH2O2
      ISPEC(10) = 3HSiO
      ISPEC(11) = 4HSiO2
      ISPEC(12) = 4HSiOH
      ISPEC(13) = 2HMg
      ISPEC(14) = 3HMgO
C   ROCK SPECIES (INERT SPECIES)
      ISPEC(15) = 6HMgSiO3
      ISPEC(16) = 7HMg2SiO4

C ****************************************************************
      DIRDATA='DATA'
C Input files
      OPEN(UNIT=7,FILE='species+T_in.dat')
      OPEN(UNIT=9,FILE=DIRDATA//'/CHEM.DAT.CH4_OMIF')
C
C Output files
      OPEN(UNIT=98,FILE='MIFprintout.dat')
C
      OPEN(UNIT=8,FILE='species+T_out.dat')
C
      OPEN(UNIT=15,FILE='OUT/int.rates.out.dat')
C
C ***** READ THE CHEMISTRY DATA CARDS *****
      READ(9,200)JCHEM
 200  FORMAT(10X,A8,2X,A8,2X,A8,2X,A8,2X,A8)
C     print 201,(J,(JCHEM(M,J),M=1,5),J=1,NR)
C201  FORMAT(1X,I3,1H),5X,A8,4H +  ,A8,7H  =    ,A8,4H +  ,A8,4X,A8)
C
C     READS JCHEM AND PRINTS WELL.
C
C ***** REPLACE HOLLERITH LABELS WITH SPECIES NUMBERS IN JCHEM *****
      DO 5 J=1,NR
      DO 5 M=1,5
      IF(JCHEM(M,J).EQ.1H ) GO TO 5
      DO 6 I=1,NSP
      IF(JCHEM(M,J).NE.ISPEC(I)) GO TO 6
      JCHEM(M,J) = I
      GO TO 5
   6  CONTINUE
      IERR = J
      GO TO 25
   5  CONTINUE
      PRINT *
C
C     print 401,(J,(JCHEM(M,J),M=1,5),J=1,NR)
C401  FORMAT(1X,I3,1H),5X,I2,4H +  ,I2,7H  =    ,I2,4H +  ,I2,4X,I2)
C     LABLES REPLACED, PRINT FOR PROOFIES
C
C ***** Read character array for P&L tables, "int.rates.out.dat"
C     REWIND 9
C     READ(9,200)CHEM
C
C ***** FILL UP CHEMICAL PRODUCTION AND LOSS MATRICES *****
      DO 7 M=1,2
      N = 3-M
      DO 7 J=1,NR
      I = JCHEM(M,J)
      IF(I.LT.1.OR.I.GT.NSP) GO TO 7
      NUML(I) = NUML(I) + 1
      IF(NUML(I).GT.NMAX) GO TO 20
      K = NUML(I)
      ILOSS(1,I,K) = J
      ILOSS(2,I,K) = JCHEM(N,J)
   7  CONTINUE
C
      DO 8 M=3,5
      DO 8 J=1,NR
      I = JCHEM(M,J)
      IF(I.LT.1.OR.I.GT.NSP) GO TO 8
      NUMP(I) = NUMP(I) + 1
      IF(NUMP(I).GT.NMAX) GO TO 20
      K = NUMP(I)
      IPROD(I,K) = J
   8  CONTINUE
C
C-AD ***** PRINT PROD/LOSS MATRICIES *****
C     print *
C     DO I=1,NSP
C     print 689,I,(IPROD(I,K),K=1,NMAX)
C689  FORMAT(1X,I3,3X,10(I2,1X))
C     ENDDO
C
C     print *,' '
C     DO I=1,NSP
C     print 689,I,(ILOSS(1,I,K),K=1,NMAX)
C     ENDDO
C
C     print *,' '
C     DO I=1,NSP
C     print 689,I,(ILOSS(2,I,K),K=1,NMAX)
C     ENDDO
C
C     PROD/LOSS MATRICIES WORK AND PRINT
C
C ***** READ THE INPUT DATAFILE species+T_in.dat *****
      READ(7,500) USOL,T
 500  FORMAT(1P1E10.3)
C
      print 1899,USOL,T
 1899 FORMAT(1X,1P1E10.3)
C
      DEN = 0.
      DO I=1,NQ
      DEN = DEN + USOL(I)
      ENDDO
C
C ***** SET MODEL PARAMETERS *****
C     DT = INITIAL TIME STEP
C     TSTOP = TIME AT WHICH CALCULATION IS TO STOP
C     NSTEPS = NUMBER OF TIME STEPS TO RUN (IF TSTOP IS NOT REACHED)
C 
      EPSJ = 1.E-7
C
C-AD ****** OUTPUT FILE HEADER ******
c-as   print the head of the output file
        write(98, 700) 
  700   format('********************************************',
     & /2x,'OUTPUT FOR BOX MODEL OF SOLAR NEBULA',f6.2,/,
     & '********************************************')
      write(98, 202) NQ !,KJAC
 202  FORMAT(//1X,'NQ = ',I2) !,5X,'KJAC=',I7  
C
C-AD *****************************

      CALL RATES(T,DEN) !DEFINE T AND D(USOL) species+T_in.dat

C     DO I=1,NR
C     print 534,I,A(I)
C534  FORMAT(1X,I3,3X,I2.3)
C     ENDDO
C     print *
C LOTSA PROBLEMS HERE IN RATES AND PRINTING RATE MATRIX
C
C ***** PRINT OUT INITIAL DATA *****
      CALL OUTPUT(0,NSTEPS,0.)
C
C   SET JACOBIAN PARAMETERS
      KD = 2*NQ + 1
      KU = KD - NQ
      KL = KD + NQ
C
C   PRINT OUT RESULTS EVERY NPR TIME STEPS
      NPR = NSTEPS
      PRN = NPR
C
C ***** START THE TIME-STEPPING LOOP *****
      TIME = 0.
      NORDER = 1
      STEPCOUTNT = 0
      TSTOP = 1.E17

      DO 1 N=1,NSTEPS

      STEPCOUNT = STEPCOUNT + 1

      IF(NORDER.EQ.1) THEN
     	DTINV = 1./DT

      ELSE IF(NORDER.EQ.2) THEN
      	DTINV = (3./2.)/DT
      END IF

      TIME = TIME + DT
C-AD
C-AD
C ***** SET UP THE JACOBIAN MATRIX AND RIGHT-HAND SIDE *****
      DO 17 J=1,NQ
      DO 17 K=1,NQ
  17  DJAC(J,K) = 0.
      DO 19 K=1,NQ
  19  RHS(K) = 0.
C
C   (DJAC IS EQUAL TO (1/DT)*I - J, WHERE J IS THE JACOBIAN MATRIX)
C
C   COMPUTE CHEMISTRY TERMS AT ALL GRID POINTS
      IDO = 0
      IF (NN.EQ.NSTEPS) IDO = 1
      CALL DOCHEM(FVAL,IDO)
C
      DO 9 I=1,NQ
      RHS(I) = FVAL(I)
   9  USAVE(I) = USOL(I)
C
      DO 3 I=1,NQ
      R = EPSJ * ABS(USOL(I))
  11  USOL(I) = USAVE(I) + R
      CALL DOCHEM(FV,0)
C
      DO 12 K=1,NQ  !CHANGED M INDEX TO K INDEX
      M = M - I + KD
  12  DJAC(K,I) = (FVAL(I) - FV(I))/R
C
      USOL(I) = USAVE(I)
      DJAC(I,I) = DJAC(I,I) + DTINV
   3  CONTINUE
C
C-AD
C ***** FACTOR THE JACOBIAN AND SOLVE THE LINEAR SYSTEM *****
C 
      RCOND = 0.
      CALL LUDCMP(DJAC,N,NP,INDX,D)
C     write(99,152) N,RCOND,DT,USOL
C152  FORMAT(I4,2X,1P16E22.15,2X,1P16E22.15,2X,1P16E22.15)

      IF(INDEX.NE.0) write(98, 103) N,INDEX
 103  FORMAT(/1X,'N =',I3,5X,'INDEX =',I9)
      CALL LUBKSB(DJAC,N,NP,INDX,RHS)

C   COMPUTE NEW CONCENTRATIONS
      EMAX = 0.
      DO 26 I=1,NQ
C
      REL(I) = RHS(I)/USOL(I)
      EREL = ABS(REL(I))
      EMAX = AMAX1(EMAX,EREL)!RAT
      IF(EREL.LT.EMAX) GO TO 26
      IS = I
      UMAX = USOL(I)
      RMAX = RHS(I)
  26  USOL(I) = USOL(I) + RHS(I) ! X = X + DX
C
C
      KCOUNT = 1
C
ccccccccccccccccccccccccc
c INCREASE STEP SIZE 2x c
ccccccccccccccccccccccccc
C
C   AUTOMATIC TIME STEP CONTROL
      DTSAVE = DT
      IF(EMAX.GT.0.15)  DT = 0.9*DTSAVE
      IF(EMAX.GT.0.20)  DT = 0.7*DTSAVE
      IF(EMAX.LT.0.10)  DT = 1.1*DTSAVE
      IF(EMAX.LT.0.05)  DT = 1.3*DTSAVE
      IF(EMAX.LT.0.03)  DT = 1.5*DTSAVE
      IF(EMAX.LT.0.01)  DT = 2.0*DTSAVE
      IF(EMAX.LT.0.003) DT = 5.0*DTSAVE
      IF(EMAX.LT.0.001) DT = 10.*DTSAVE
      DTINV = 1./DT
C		
      ISP = ISPEC(I)
      IF(SM-MS.GT.0.01) GO TO 317
      write(98, 100) N,EMAX,ISP,UMAX,RMAX,DT,TIME
 100  FORMAT(1X,'N =',I4,2X,'EMAX =',1PE9.2,
     2  1X,'U =',E9.2,1X,'RHS =',E9.2,
     3  2X,'DT =',E9.2,2X,'TIME =',E9.2)
 317  CONTINUE
C
      DEN = 0.
      DO I=1,NQ
      DEN = DEN + USOL(I)
      ENDDO
C      
      IF (EMAX.LT.0.5) GO TO 28
      print *, 'Decreasing Step'
      DT = 0.5*DTSAVE     
      TIME = TIME - DTSAVE

      DO 27 I=1,NQ
  27  USOL(I) = USAVE(I)
  28  CONTINUE
C
      NS = N/NPR
      SN = N/PRN
C
      IF(NN.EQ.NSTEPS) CALL OUTPUT(NN,NSTEPS,TIME)
  37  CONTINUE
      IF(INDEX.NE.0) STOP
      IF(NN.EQ.NSTEPS) GO TO 22
      IF(TIME.GT.TSTOP) NN = NSTEPS - 1
  22  CONTINUE
  1   CONTINUE
C ***** END THE TIME-STEPPING LOOP *****
C
      WRITE(8,501) USOL,T
 501  FORMAT(1P1E10.3)

C print out P&L tables with integrated rxn rates, "int.rates.out.dat"
      DO 702 I=1,NSP
         ISP = ISPEC(I)
         WRITE(15,703) ISP,TP(I)
 703     FORMAT(/A8,12X,'PRODUCTION RXS',14X,'INT RX RATE',4X,
     2      'TP = ',1PE9.2)
       DO 704 N=1,NR 
          IF(JCHEM(3,N).EQ.I .OR. JCHEM(4,N).EQ.I .OR. 
     2       JCHEM(5,N).EQ.I)THEN
           IF(RAT(N).NE.0.) WRITE(15,705) N,(CHEM(J,N),J=1,5),RAT(N)
 705       FORMAT(1X,I3,1H),1X,A7,3H + ,A7,3H = ,A7,3H + ,A6,2X,A4,
     2      1PE10.3)
          ENDIF
 704   CONTINUE
C
         WRITE(15,706) ISP,TL(I)
 706     FORMAT(/A8,15X,'LOSS RXS',16X,'INT RX RATE',4X,'TL = ',1PE9.2)
       DO 707 N=1,NR 
          IF(JCHEM(1,N).EQ.I .OR. JCHEM(2,N).EQ.I)THEN
             IF(RAT(N).NE.0.) WRITE(15,705) N,(CHEM(J,N),J=1,5),RAT(N)
          ENDIF
 707   CONTINUE
 702  CONTINUE
C
      GO TO 21
  20  write(98, 300) I
 300  FORMAT(//1X,'NMAX EXCEEDED FOR SPECIES ',I3)
      GO TO 21
  25  write(98, 301) IERR
 301  FORMAT(//1X,'ERROR IN REACTION ',I3)

  21  CONTINUE

      STOP
      END PROGRAM ARC_CH4

C-PK ********************************
      SUBROUTINE RATES(T,DEN)
      PARAMETER(NQ=14, NR=28, NSP=16, NMAX=10)
      COMMON/RBLOK/A(NR),ILOSS(2,NSP,NMAX),IPROD(NSP,NMAX),
     2  JCHEM(5,NR),NUML(NSP),NUMP(NSP)
C
C
C ***** FILL UP RATE MATRIX *****
C     DO 4 I=1,NR !NR FOR NUMBER OF REACTIONS
C   A MATRIX DIMENSIONED NR, K COEFFS
C   GO DIRECTLY INTO RATE MATRIX A
C   RATE CONSTANTS TAKEN FROM THIEMENS_13_SUPP_INFO TABLE S5
      A(1) = 1.0E-11*EXP(110/T) !ASSUMED
      A(2) = 3.0E-15            !ASSUMED
      A(3) = 1.0E-12            !ASSUMED
      A(4) = 5.0E-12            !ASSUMED
      A(5) = 1.0E-12            !ASSUMED
      A(6) = 1.0E-12            !ASSUMED
      A(7) = 1.0E-12            !ASSUMED
      A(8) = 6.0E-34
      A(9) = 8.5E-20*(T**2.67)*EXP(-3160./T)
      A(10) = 2.4E-11*EXP(-352./T)
      A(11) = 2.5E-15*(T**1.14)*EXP(-50./T)
      A(12) = 1.0E-12
      A(13) = 3.3E-10*EXP(-8460./T)
      A(14) = 5.3E-11
      A(15) = 1.1E-12*EXP(-2000./T)
      A(16) = 1.7E-16*(T**1.6)*EXP(-1660./T)
      A(17) = 4.8E-11*EXP(250./T)
      A(18) = 1.3E-11*EXP(-670./T)
      A(19) = 3.1E-12*EXP(-775./T)
      A(20) = 1.9E-12*EXP(-1000./T)
      A(21) = 1.4E-14*EXP(-600./T)
      A(22) = 1.0E-12           !MADE-UP KASTING
      A(23) = 1.0E-12           !MADE-UP KASTING
      A(24) = 1.0E-12           !MADE-UP KASTING
C
C
C ***** THREE-BODY COEFFICIENTS *****
C     A(I) = TBDY(K0,KI,N,M,T,DEN)
C
C  H + O2 + H2O = HO2 + H2O
      A(25) = TBDY(4.3E-30,1.E-10,-0.8,0.,T,DEN)
C 
C  H + O2 + H2 = HO2 + H2
      A(26) = TBDY(5.8E-30,1.E-10,-0.8,0.,T,DEN)
C 
C  H + H + H2 = H2 + H2
      A(27) = TBDY(2.7E-31,1.E-10,-0.6,0.,T,DEN)
C 
C  H + OH + H2O = H2O + H2O
      A(28) = TBDY(3.9E-25,1.E-10,-2.,0.,T,DEN)
C
      print *, 'A = '
      print 100, A
 100  FORMAT (1P6E10.3)

      RETURN
      END

C-AD ***************************************************
      FUNCTION TBDY(A0,AI,CN,CM,T,D)                      
      B0 = A0*(T**CN)                                
      BI = AI*(T**CM)                                
      Y = ALOG10(B0*D/BI)                                 
      X = 1./(1. + Y**2)                                  
      TBDY = B0*D/(1. + B0*D/BI) * 0.6**X
      RETURN
      END

C-AD ***************************************************
      SUBROUTINE OUTPUT(N,NSTEPS,TIME)
      PARAMETER (NQ=14, NR=28, NSP=16, NMAX=10)
      DIMENSION TP(NSP),TL(NSP),SL(NSP)
      DIMENSION USOL(NQ)
      COMMON/NBLOK/LO,LO2,LO3,LH,LH2,LOH,LH2O,LHO2,LH2O2,
     2  LSiO,LSiO2,LSiOH,LMg,LMgO,LMgSiO3,LMg2SiO4
C
C   THIS SUBROUTINE PRINTS OUT ALL THE DATA. THE VARIABLE ISKIP 
C   SAYS HOW MANY POINTS YOU WANT TO LOOK AT.
C
C   SET ISKIP=1 FOR MIXING RATIO USOL (NUMBER DENSITY) AT EVERY ITERATION
      ISKIP = 4
      JSKIP = ISKIP
      IF(ISKIP.GT.1 .AND. N.EQ.NSTEPS) ISKIP = 2
      TIMEY = TIME/3600./24./365.25
      write(98, 100) TIME,TIMEY
 100  FORMAT(/1X,'TIME =', 1PE9.2,5X,'TIMEY =',E9.2,1X,'YEARS')
      write(98, 101) USOL
 101  FORMAT(/1X,'USOL =',I3)
C
C ***** PRINT ON LAST ITERATION ONLY *****
      write(98, 125)
 125  FORMAT(/1X,'NUMBER DENSITIES OF SPECIES'/)

      DO 1 K=1,NQ
      SL(K) = USOL(K)*DEN
   1  CONTINUE

      RETURN
      END

C-AD **********************************
      SUBROUTINE DOCHEM(FVAL,N)                           
      PARAMETER NQ=14                  
      PARAMETER(NR=28, NSP=16, NMAX=10)                                           
      DIMENSION FVAL(NQ),D(NQ),YP(NQ),YL(NQ),USOL(NQ)
      DIMENSION RXTOT(NR),RAT(NR)
      COMMON/NBLOK/LO,LO2,LO3,LH,LH2,LOH,LH2O,LHO2,LH2O2,
     2  LSiO,LSiO2,LSiOH,LMg,LMgO,LMgSiO3,LMg2SiO4
      COMMON/RBLOK/A(NR),ILOSS(2,NSP,NMAX),IPROD(NSP,NMAX),
     2  JCHEM(5,NR),NUML(NSP),NUMP(NSP)
C                                                         
C   THIS SUBROUTINE DOES THE CHEMISTRY BY CALLING CHEMPL.  PHOTO-         
C   CHEMICAL EQUILIBRIUM SPECIES ARE DONE FIRST.  THESE MUST CON-         
C   TAIN NO NONLINEARITIES AND MUST BE DONE IN THE PROPER ORDER 
C   (I.E. IF SPECIES A REACTS TO FORM B, THEN A MUST BE FOUND FIRST).
C   LONG-LIVED SPECIES CAN BE DONE IN ANY ORDER.
C
C ***** LONG-LIVED SPECIES CHEMISTRY ***** IN JAC MATRIX
      DO 4 I=1,NQ
      D(I) = USOL(I) !NUMBER DENSITY
      CALL CHEMPL(D,XP,XL,I)
      FVAL = XP - XL*USOL(I)
      YP(I) = XP
      YL(I) = XL
      IF (N.LT.1) RETURN

C ***** CALCULATE PRODUCTION AND LOSS *****
      DO 10 L=1,NR
  10  RAT(L) = 0.
C
      DO 12 L=1,NR
      M = JCHEM(1,L)
      K = JCHEM(2,L)
      RXTOT(L) = A(L)*USOL(M)*USOL(K)
  12  RAT(L) = RAT(L) + A(L)*USOL(M)*USOL(K)
   4  CONTINUE
      RETURN

      END
C-PK *******************************
      SUBROUTINE CHEMPL(D,XP,XL,K)
      PARAMETER NQ=14
      PARAMETER(NR=28, NSP=16, NMAX=10)
      DIMENSION D(NQ)
      COMMON/RBLOK/A(NR),ILOSS(2,NSP,NMAX),IPROD(NSP,NMAX),
     2  JCHEM(5,NR),NUML(NSP),NUMP(NSP)
C
C   THIS SUBROUTINE CALCULATES CHEMICAL PRODUCTION AND LOSS RATES
C   USING THE INFORMATION IN THE MATRICES JCHEM, ILOSS, AND IPROD.
C   CALLED BY SUBROUTINE DOCHEM.
C
      XP = 0.
      XL = 0.
C
C   LOSS FREQUENCY XL
      NL = NUML(K)
      DO 2 L=1,NL
      J = ILOSS(1,K,L)
      M = ILOSS(2,K,L)
   2  XL = XL + A(J)*D(M)
C
C   PRODUCTION RATE XP
      NP = NUMP(K)
      DO 3 L=1,NP
      J = IPROD(K,L)
      M = JCHEM(1,J)
      N = JCHEM(2,J)
   3  XP = XP + A(J)*D(M)*D(N)
C
      RETURN
      END

C-AD ****************************************
C     SUBROUTINE ISOTOPE(NQ,NR,NSP)
C     PARAMETER (NQ=25,NR=95,NSP=29,NMAX=10)
C
C     COMMON/NBLOK/LO,LXO,LO2,LOXO,LO3,LO2XO,LH,LH2,LOH,LXOH,LH2O,LH2XO,
C    2  LHO2,LHOXO,LH2O2,LH2OXO,LSiO,LSiXO,LSiO2,LSiOXO,LSiOH,LSiXOH,
C    3  LMg,LMgO,LMgXO,LMgSiO3,LMgSiO2XO,LMg2SiO4,LMg2SiO3XO
C
C     DATA LO,LXO,LO2,LOXO,LO3,LO2XO,LH,LH2,LOH,LXOH,LH2O,LH2XO,
C    2  LHO2,LHOXO,LH2O2,LH2OXO,LSiO,LSiXO,LSiO2,LSiOXO,LSiOH,LSiXOH,
C    3  LMg,LMgO,LMgXO,LMgSiO3,LMgSiO2XO,LMg2SiO4,LMg2SiO3XO/
C    4  1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23
C    5  24,25,26,27,28,29/
C
C ***** SPECIES DEFINITIONS *****
C   LONG-LIVED - ALL OUR SPECIES ARE "LONG-LIVED"
C     ISPEC(1) = 1HO
C     ISPEC(2) = 2HXO
C     ISPEC(3) = 2HO2
C     ISPEC(4) = 3HOXO
C     ISPEC(5) = 2HO3
C     ISPEC(6) = 4HO2XO
C     ISPEC(7) = 1HH
C     ISPEC(8) = 2HH2
C     ISPEC(9) = 2HOH
C     ISPEC(10) = 3HXOH
C     ISPEC(11) = 3HH2O
C     ISPEC(12) = 4HH2XO
C     ISPEC(13) = 3HHO2
C     ISPEC(14) = 4HHOXO
C     ISPEC(15) = 4HH2O2
C     ISPEC(16) = 5HH2OXO
C     ISPEC(17) = 3HSiO
C     ISPEC(18) = 4HSiXO
C     ISPEC(19) = 4HSiO2
C     ISPEC(20) = 5HSiOXO
C     ISPEC(21) = 4HSiOH
C     ISPEC(22) = 5HSiXOH
C     ISPEC(23) = 2HMg
C     ISPEC(24) = 3HMgO
C     ISPEC(25) = 4HMgXO
C   ROCK SPECIES (INERT SPECIES)
C     ISPEC(26) = 6HMgSiO3
C     ISPEC(27) = 8HMgSiO2XO
C     ISPEC(28) = 7HMg2SiO4
C     ISPEC(29) = 9HMg2SiO3XO
C
C ***** NEW REACTION LIST IN CHEM.DAT.CH4_ISOTOPE*****
C

C-PK ****************************************
      SUBROUTINE LUDCMP(A,N,NP,INDX,D)
      PARAMETER (NMAX=10,TINY=1.0E-20)
      DIMENSION A(NP,NP), INDX(N), VV(NMAX)
      D=1.
      DO 42 I=1,N
        AAMAX=0.
        DO 41 J=1,N
          IF (ABS(A(I,J)).GT.AAMAX) AAMAX=ABS(A(I,J))
41      CONTINUE
        IF (AAMAX.EQ.0.) PAUSE 'Singular matrix.'
        VV(I)=1./AAMAX
42    CONTINUE
      DO 49 J=1,N
        IF (J.GT.1) THEN
          DO 44 I=1,J-1
            SUM=A(I,J)
            IF (I.GT.1)THEN
              DO 43 K=1,I-1
                SUM=SUM-A(I,K)*A(K,J)
43            CONTINUE
              A(I,J)=SUM
            ENDIF
44        CONTINUE
        ENDIF
        AAMAX=0.
        DO 46 I=J,N
          SUM=A(I,J)
          IF (J.GT.1)THEN
            DO 45 K=1,J-1
              SUM=SUM-A(I,K)*A(K,J)
45          CONTINUE
            A(I,J)=SUM
          ENDIF
          DUM=VV(I)*ABS(SUM)
          IF (DUM.GE.AAMAX) THEN
            IMAX=I
            AAMAX=DUM
          ENDIF
46      CONTINUE
        IF (J.NE.IMAX)THEN
          DO 47 K=1,N
            DUM=A(IMAX,K)
            A(IMAX,K)=A(J,K)
            A(J,K)=DUM
47        CONTINUE
          D=-D
          VV(IMAX)=VV(J)
        ENDIF
        INDX(J)=IMAX
        IF(J.NE.N)THEN
          IF(A(J,J).EQ.0.)A(J,J)=TINY
          DUM=1./A(J,J)
          DO 48 I=J+1,N
            A(I,J)=A(I,J)*DUM
48        CONTINUE
        ENDIF
49    CONTINUE
      IF(A(N,N).EQ.0.)A(N,N)=TINY
      RETURN
      END SUBROUTINE LUDCMP

C-PK ****************************************
      SUBROUTINE LUBKSB(A,N,NP,INDX,B)
      DIMENSION A(NP,NP),INDX(N),B(N)

      II=0
      DO 112 I=1,N
        LL=INDX(I)
        SUM=B(LL)
        B(LL)=B(I)
        IF (II.NE.0)THEN
          DO 111 J=II,I-1
            SUM=SUM-A(I,J)*B(J)
111        CONTINUE
        ELSE IF (SUM.NE.0.) THEN
          II=I
        ENDIF
        B(I)=SUM
112    CONTINUE
      DO 114 I=N,1,-1
        SUM=B(I)
        IF(I.LT.N)THEN
          DO 113 J=I+1,N
            SUM=SUM-A(I,J)*B(J)
113        CONTINUE
        ENDIF
        B(I)=SUM/A(I,I)
114    CONTINUE
      RETURN
      END SUBROUTINE LUBKSB

