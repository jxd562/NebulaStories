         PROGRAM ARC_CH4 

C-JKD June 26, 2017 removed Mg species--current error: singular matrix

C-JKD Modified on June 21, 2017 for mass independent fractionation of oxygen.
C-AD Deleted parts of original photochemical code that are not needed
C-to perform the box model of the solar nebula.
C
C       THIS PROGRAM IS A BOX MODEL OF THE SOLAR NEBULA. CHEMICAL SPECIES
C     ARE INJECTED INTO THE BOX, AND REACTIONS TAKE PLACE, WITH THE INTENT 
C     TO SEE WHETHER OR NOT THERE IS MASS INDEPENDT FRACTIONATION OF OXYGEN
C     IN THE REACTIONS.
C
C     THE IMPORTANT PARAMETERS DESCRIBING THE CHEMISTRY ARE
C        NR   = NUMBER OF REACTIONS
C        NSP  = NUMBER OF CHEMICAL SPECIES
C        NQ   = NUMBER OF SPECIES FOR WHICH A DIFFUSION EQUATION
C               IS SOLVED AND WHICH ARE IN THE BIG, REVERSE EULER MATR
C        NMAX = MAXIMUM NUMBER OF REACTIONS IN WHICH AN INDIVIDUAL
C               SPECIES PARTICIPATES
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
C     1)  SiO + O2 = SiO2 + O       !Thiemens#1
C     2)  SiO + H2 = SiOH + H       !Thiemens#2
C     3)  SiO + OH = SiO2 + H       !Thiemens#3
C     4)  SiO + HO2 = SiO2 + H      !Thiemens#4
C     5)  SiO + O = SiO2            !Thiemens#5
C     6)  SiO + H = SiOH            !Thiemens#6
C     7)  SiO + O3 = SiO2 + O2      !Thiemens#7
C     8)  H2 + O = OH + O           !Thiemens#9
C     9)  OH + O = O2 + H          !Thiemens#10
C    10)  OH + OH = H2O + O         !Thiemens#11
C    11)  H + O = OH                !Thiemens#12
C    12)  H + O2 = OH + O           !Thiemens#13
C    13)  O + HO2 = OH + O2         !Thiemens#14
C    14)  O + H2O2 = OH + HO2       !Thiemens#15
C    15)  OH + H2 = H2O + H         !Thiemens#16
C    16)  OH + HO2 = H2O + O2       !Thiemens#17
C    17)  OH + H2O2 = H2O + HO2     !Thiemens#18
C    18)  HO2 + HO2 = H2O2 + O2     !Thiemens#24
C    19)  OH + O3 = HO2 + O2        !Thiemens#25
C    20)  HO2 + O3 = OH + O2 + O2 !(2O2)    !Thiemens#1
C
C ***** THREE BODY RXN ***** ELIMINATE COMMON TERMS 
C    21)  O + O2 + M = O3 + M           !Thiemens#8
C    22)  H + O2 + H2O(M) = HO2 + H2O(M)!Thiemens#20
C    23)  H + O2 + H2(M) = HO2 + H2(M)  !Thiemens#21
C    24)  H + H + H2(M) = H2 + H2(M)    !Thiemens#22
C    25)  H + OH + H2O(M) = H2O + H2O(M)!Thiemens#23
C
C    26)  MgO + SiO2 = MgSiO3       !check this
C REMOVED   27)  Mg + O = MgO              !check this
C    28)  MgO + MgSiO3 = Mg2SiO4    !check this
C
C        THIS PROGRAM DOES THE CHEMISTRY AUTOMATICALLY.  THE CHEMICAL
C     REACTIONS ARE ENTERED ON DATA CARDS (.DAT) IN FIVE 10-DIGIT 
C     COLUMNS STARTING IN COLUMN 11, I.E.                         
C
C         REAC1     REAC2     PROD1     PROD2     PROD3
C
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
      PARAMETER(NQ=13, NR=27, NSP=15, NMAX=10)
      DIMENSION FVAL(NQ),FV(NQ),DJAC(NQ,NQ),RHS(NQ),REL(NQ),IPVT(NQ)
     2  ,USAVE(NQ),USOL(NQ),TP(NQ),TL(NQ),D(NQ),INDX(NQ)
C     DIMENSION RAT(NR)
      CHARACTER*30 CHEM(5,NR),PRODRX(NSP,NR),LOSSRX(NSP,NR)
      CHARACTER*4,DIRDATA
     
      COMMON/NBLOK/LO,LO2,LO3,LH,LH2,LOH,LH2O,LHO2,LH2O2,
     2  LSiO,LSiO2,LSiOH,LMg,LMgO,LMgSiO3,LMg2SiO4
      COMMON/RBLOK/A(NR),ILOSS(2,NSP,NMAX),IPROD(NSP,NMAX),
     2  JCHEM(5,NR),NUML(NSP),NUMP(NSP),ISPEC(NSP),T,DEN
      COMMON/ZBLOK/RAT(NR)
C-AP *************************************************************

      DATA LO,LO2,LO3,LH,LH2,LOH,LH2O,LHO2,LH2O2,
     2  LSiO,LSiO2,LSiOH,LMgO,LMgSiO3,LMg2SiO4/
     3  1,2,3,4,5,6,7,8,9,10,11,12,13,14,15/
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
C      ISPEC(13) = 2HMg
      ISPEC(13) = 3HMgO
C   ROCK SPECIES (INERT SPECIES)
      ISPEC(14) = 6HMgSiO3
      ISPEC(15) = 7HMg2SiO4

C ****************************************************************
      DIRDATA='DATA'
C Input files
      OPEN(UNIT=7,FILE='species+T_in.dat') !these are mixing ratios +T
      OPEN(UNIT=9,FILE=DIRDATA//'/CHEM.DAT.CH4_OMIF')!reactions 
C
C Output files
      OPEN(UNIT=98,FILE='MIFprintout.dat')!output file
C
      OPEN(UNIT=8,FILE='species+T_out.dat')!resulting mixing ratios +T
C
      OPEN(UNIT=15,FILE='OUT/int.rates.out.dat')
C
C ***** READ THE CHEMISTRY DATA CARDS *****
      READ(9,200)JCHEM !reads the names, matched to hollerith
 200  FORMAT(10X,A8,2X,A8,2X,A8,2X,A8,2X,A8)
      print 201,(J,(JCHEM(M,J),M=1,5),J=1,NR)
 201  FORMAT(1X,I3,1H),5X,A8,4H +  ,A8,7H  =    ,A8,4H +  ,A8,4X,A8)
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
      PRINT *,'Reactant #'
C
      print 401,(J,(JCHEM(M,J),M=1,5),J=1,NR)
 401  FORMAT(1X,I3,1H),5X,I2,4H +  ,I2,7H  =    ,I2,4H +  ,I2,4X,I2)
C     JCHEM PRINTS WELL, MATCHES REACTION LIST
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
      print *, 'PROD'
      DO I=1,NSP
      print 689,I,(IPROD(I,K),K=1,NMAX)
 689  FORMAT(1X,I3,3X,10(I2,1X))
      ENDDO
C
      print *,'LOSS 1'
      DO I=1,NSP
      print 689,I,(ILOSS(1,I,K),K=1,NMAX)
      ENDDO
C
      print *,'LOSS 2'
      DO I=1,NSP
      print 689,I,(ILOSS(2,I,K),K=1,NMAX)
      ENDDO
C
C     PROD/LOSS MATRICIES WORK
C
C ***** READ THE INPUT DATAFILE species+T_in.dat *****
      READ(7,500) USOL,T
 500  FORMAT(1P1E10.3)
C
C     print 1899,USOL,T
C1899 FORMAT(1X,1P1E10.3) !Reads species+T_in.dat 
 
      DEN = 1E15 !5.17598E14 FROM P=NKT

      DO I=1,NQ
      D(I) = USOL(I)*(DEN)
C individual number densities D(I) for each species (#NQ)
      ENDDO
C
C     print *, 'D ='
C     print 154, D
C154  FORMAT (1P1E10.3) 
C
C   **INDIVIDUAL NUMBER DENSITIES ARE GOOD**
C
C ***** SET MODEL PARAMETERS *****
C     DT = INITIAL TIME STEP
C     TSTOP = TIME AT WHICH CALCULATIONmak IS TO STOP
C     NSTEPS = NUMBER OF TIME STEPS TO RUN (IF TSTOP IS NOT REACHED)
C
C     EPSJ = AMOUNT BY WHICH TO PERTURB THINGS FOR JACOBIAN CALCULATION

      EPSJ = 1.E-7 
C
C-AD ****** OUTPUT FILE HEADER ******
c-as   print the head of the output file
        write(98, 700) 
  700   format('********************************************',
     & /2x,'OUTPUT FOR BOX MODEL OF SOLAR NEBULA',f6.2,/,
     & '********************************************')
      write(98, 202) NQ 
 202  FORMAT(//1X,'NQ = ',I2)   
C
C-AD *****************************

      CALL RATES(D)
C
C ***** PRINT OUT INITIAL DATA *****
      CALL OUTPUT(USOL,0,NSTEPS,0.)
C
C   PRINT OUT RESULTS EVERY NPR TIME STEPS
C      NPR = NSTEPS
C      PRN = NPR
C
C ***** START THE TIME-STEPPING LOOP *****
      TIME = 0.
      DT = 1.E-8    !SIZE OF TIME STEP
      STEPCOUNT = 0
      TSTOP = 1.E17
      NSTEPS = 1    !~1000 when we increase Time

      DO 1 N=1,NSTEPS

      STEPCOUNT = STEPCOUNT + 1

      DTINV = 1./DT

      TIME = TIME + DT
C-AD
C ***** SET UP THE JACOBIAN MATRIX AND RIGHT-HAND SIDE *****
      DO 17 J=1,NQ
      DO 17 K=1,NQ
  17  DJAC(J,K) = 0.
      DO 19 K=1,NQ
  19  RHS(K) = 0.
C
C     print *
C     print *, 'DJAC BEFORE DOCHEM'
C     print 670, ((DJAC(I,J),J=1,NQ),I=1,NQ)
C670  FORMAT (1P14E8.1)
C     ALL ZEROES
C
C   (DJAC IS EQUAL TO (1/DT)*I - J, WHERE J IS THE JACOBIAN MATRIX)
C
C   COMPUTE CHEMISTRY TERMS
C
      IDO = 0 !this is a diagnostic for last time step
      IF (NN.EQ.NSTEPS) IDO = 1
      CALL DOCHEM(USOL,FVAL,IDO) !okay to have USOL b/c not in cmn blk
C
C     print *, 'D = '
C     print 155, D
C155  FORMAT (1P1E10.3)
C
      print *
      print *, 'FVAL'
      print 668, FVAL
 668  FORMAT (1P7E10.1)
C
      DO 9 I=1,NQ
      RHS(I) = FVAL(I)
   9  USAVE(I) = USOL(I)
C
      DO 3 I=1,NQ
      R = EPSJ * ABS(USOL(I)) !reminder, EPSJ-perturbation to Jac.
      USOL(I) = USAVE(I) + R
      CALL DOCHEM(USOL,FV,0)
C
      DO 12 K=1,NQ 
  12  DJAC(K,I) = (FVAL(K) - FV(K))/R
C
      USOL(I) = USAVE(I)
      DJAC(I,I) = DJAC(I,I) + DTINV
   3  CONTINUE
C
C-AD
      print *
      print *, 'DJAC BEFORE MATRIX SOLVER'
      print 666, ((DJAC(I,J),J=1,NQ),I=1,NQ)
 666  FORMAT (1P14E8.1)

      print *, 'D ='
      print 157, D
 157  FORMAT (1P1E10.3) 
C   **D VALUES LOOK GOOD BEFORE MATRIX SOLVER**
C
C ***** FACTOR THE JACOBIAN AND SOLVE THE LINEAR SYSTEM *****
C 
      CALL LUDCMP(DJAC,NQ,NQ,INDX,DP) !N, NP GOES TO NQ, NQ,, D -> DP 
      CALL LUBKSB(DJAC,NQ,NQ,INDX,RHS)

C   COMPUTE NEW CONCENTRATIONS
      EMAX = 0.
      DO 26 I=1,NQ
C
      REL(I) = RHS(I)/USOL(I)
      EREL = ABS(REL(I))
      EMAX = AMAX1(EMAX,EREL)

      IF(EREL.LT.EMAX) GO TO 26

      UMAX = USOL(I)
      RMAX = RHS(I)
  26  USOL(I) = USOL(I) + RHS(I) ! X = X + DX
C 
      print *
      print *, '  USAVE,      RHS,        USOL'
      DO I=1,NQ
      print 667, USAVE(I), RHS(I), USOL(I)
 667  FORMAT (1P3E12.3)
      ENDDO
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

      write(98, 100) N,EMAX,UMAX,RMAX,DT,TIME
 100  FORMAT(1X,'N =',I4,2X,'EMAX =',1PE9.2,
     2  1X,'U =',1PE9.2,1X,'RHS =',1PE9.2,
     3  2X,'DT =',1PE9.2,2X,'TIME =',1PE9.2)
      CONTINUE
C      
      IF (EMAX.LT.0.5) GO TO 28
C     IF EMAX > 0.5, continue and reject current time step
C     Stability of Reverse Euler--maintain between 0.1 and 0.2
      print *, 'Decreasing Step'
      DT = 0.5*DTSAVE     
      TIME = TIME - DTSAVE

      DO 27 I=1,NQ
  27  USOL(I) = USAVE(I) !to return to previous time step
  28  CONTINUE
C
C      NS = N/NPR !relic code(?)
C      SN = N/PRN
C
      CALL OUTPUT(USOL,NN,NSTEPS,TIME)
  37  CONTINUE
      IF(INDEX.NE.0) STOP
      IF(NN.EQ.NSTEPS) GO TO 22
      IF(TIME.GT.TSTOP) NN = NSTEPS - 1
  22  CONTINUE
  1   CONTINUE
C ***** END THE TIME-STEPPING LOOP *****
C
C PRINT REACTION RATES (wrt ambient density)
      print *
      print *, 'RAT'
      print 502,RAT
 502  FORMAT(1P10E10.3)
C
C PRINT INDIVIDUAL DENSITIES
C      print *
C      print *, 'D'
C      print 503,D
C 503  FORMAT(1P10E10.3)
C   **FIRST VALUE OF D IS 1, INSTEAD OF 5E+12**
C
C PRINT ORIGINAL RXN RATES (INPUT MATRIX)
C
C     print *
C     print *, 'A'
C     print 504,A
C504  FORMAT(1P10E10.3)
C   **RATE VECTOR LOOKS GOOD**
C
C
      WRITE(8,501) USOL,T
 501  FORMAT(1P1E10.3)

C print out P&L tables with integrated rxn rates, "int.rates.out.dat"
C NEED TO *REDO* AT SOME POINT, AD/JKD
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
C
      STOP
      END PROGRAM ARC_CH4

C-PK ********************************
      SUBROUTINE RATES(D)
      PARAMETER(NQ=13, NR=27, NSP=15, NMAX=10)
      COMMON/RBLOK/A(NR),ILOSS(2,NSP,NMAX),IPROD(NSP,NMAX),
     2  JCHEM(5,NR),NUML(NSP),NUMP(NSP),ISPEC(NSP),T,DEN
C
C ***** FILL UP RATE MATRIX *****
C     DO 4 I=1,NR !NR FOR NUMBER OF REACTIONS
C   A MATRIX DIMENSIONED NR (NUMBER OF REACTIONS), K COEFFS
C   GO DIRECTLY INTO MATRIX A
C   RATE CONSTANTS TAKEN FROM THIEMENS_13_SUPP_INFO TABLE S5
      A(1) = 1.0E-11*EXP(110/T) !Thiemens_assumed #1
      A(2) = 3.0E-15            !ibid #2
      A(3) = 1.0E-12            !ibid #3
      A(4) = 5.0E-12            !ibid #4
      A(5) = 1.0E-12            !ibid #5
      A(6) = 1.0E-12            !ibid #6
      A(7) = 1.0E-12            !ibid #7
      A(8) = 8.5E-20*(T**2.67)*EXP(-3160./T) !HCP #9
      A(9) = 2.4E-11*EXP(-352./T) !HCP #10
      A(10) = 2.5E-15*(T**1.14)*EXP(-50./T) !HCP #11
      A(11) = 1.0E-12 !HCP #12
      A(12) = 3.3E-10*EXP(-8460./T) !HCP #13
      A(13) = 5.3E-11 !HCP #14
      A(14) = 1.1E-12*EXP(-2000./T) !HCP #15
      A(15) = 1.7E-16*(T**1.6)*EXP(-1660./T) !HCP #16
      A(16) = 4.8E-11*EXP(250./T) !HCP #17
      A(17) = 1.3E-11*EXP(-670./T) !HCP #18
      A(18) = 3.1E-12*EXP(-775./T) !HCP #19
      A(19) = 1.9E-12*EXP(-1000./T) !HCP #24
      A(20) = 1.4E-14*EXP(-600./T) !HCP #25
C  THREE BODY REACTION K COEFFICIENTS (ABOVE)
C      A(21) = 6.0E-34 !JPL #8
C      A(22) = 4.3E-30*(T**-0.8) !HCP #20
C      A(23) = 5.8E-30*(T**-0.8) !HCP #21
C      A(24) = 2.7E-31*(T**-0.6) !HCP #22
C      A(25) = 3.9E-25*(T**-2.)  !HCP #23
C  ROCK REACTIONS
      A(26) = 1.0E-12           !MADE-UP
C      A(27) = 1.0E-12           !MADE-UP
      A(27) = 6.0E-34           !MADE-UP
C
C
C ***** THREE-BODY COEFFICIENTS *****
C     A(I) = TBDY(K0,KI,N,M,T,DEN)

C Three-body reactions should have input DEN (total num dens)
C This is because it is based on [M], which is f(DEN)
C
C    21)  O + O2 + M = O3 + M           !Thiemens#8
      A(21) = TBDY(6.0E-34,1.E-10,0,0.,T,DEN)
C 
C    22)  H + O2 + H2O(M) = HO2 + H2O(M)!Thiemens#20
      A(22) = TBDY(4.3E-30,1.E-10,-0.8,0.,T,DEN)
C 
C    23)  H + O2 + H2(M) = HO2 + H2(M)  !Thiemens#21
      A(23) = TBDY(5.8E-30,1.E-10,-0.8,0.,T,DEN)
C 
C    24)  H + H + H2(M) = H2 + H2(M)    !Thiemens#22
      A(24) = TBDY(2.7E-31,1.E-10,-.6,0.,T,DEN)

C    25)  H + OH + H2O(M) = H2O + H2O(M)!Thiemens#23
      A(25) = TBDY(3.9E-25,1.E-10,-2.,0.,T,DEN)

C High-pressure estimate for TBDY--gas kinetic rate
C
C
C     print *, 'A = '
C     print 100, A
C100  FORMAT (1P5E10.3)

      RETURN
      END

C-AD ***************************************************
      FUNCTION TBDY(A0,AI,CN,CM,T,D)                
C AO--low pressure limit, AI--high pressure limit
C CN--exponent low pressure, CM--exponent high pressure      
      B0 = A0*(T**CN) !low-pressure                               
      BI = AI*(T**CM) !high pressure                               
      Y = ALOG10(B0*D/BI)                                 
      X = 1./(1. + Y**2)                                  
      TBDY = B0*D/(1. + B0*D/BI) * 0.6**X
C     Equation above is based on three-body w/in gas phase
C     Complicated kinetics, maybe Marcus (?)
      RETURN
      END

C-AD ***************************************************
      SUBROUTINE OUTPUT(USOL,N,NSTEPS,TIME)
      PARAMETER (NQ=13, NR=27, NSP=15, NMAX=10)
      DIMENSION USOL(NQ),D(NSP)
      COMMON/ZBLOK/RAT(NR),TP(NQ),TL(NQ)
      COMMON/NBLOK/LO,LO2,LO3,LH,LH2,LOH,LH2O,LHO2,LH2O2,
     2  LSiO,LSiO2,LSiOH,LMg,LMgO,LMgSiO3,LMg2SiO4
      COMMON/RBLOK/A(NR),ILOSS(2,NSP,NMAX),IPROD(NSP,NMAX),
     2  JCHEM(5,NR),NUML(NSP),NUMP(NSP),ISPEC(NSP),T,DEN
C
C   THIS SUBROUTINE PRINTS OUT ALL THE DATA. THE VARIABLE ISKIP
C   SAYS HOW MANY POINTS YOU WANT TO LOOK AT.
C
C   SET ISKIP=1 FOR MIXING RATIO USOL AT EVERY ITERATION
      ISKIP = 4
      IF(ISKIP.GT.1 .AND. N.EQ.NSTEPS) ISKIP = 2
C
      write(98, 149)
 149  FORMAT('----------------------------------------')
      TIMEY = TIME/3600./24./365.25
C
      write(98, 100) TIME,TIMEY
 100  FORMAT(/1X,'TIME =', 1PE9.2,5X,'TIMEY =',E9.2,1X,'YEARS')
C
      write(98, 105)
 105  FORMAT(/1X,'MIXING RATIOS OF LONG-LIVED SPECIES')
C
      IROW = 7
      LR = NQ/IROW + 1         !changing nq to nsp makes tl and tp
      RL = FLOAT(NQ)/IROW + 1  !write 3 times instead of 2
      DIF = RL - LR
      IF (DIF.LT.0.001) LR = LR - 1
C
      DO 8 L=1,LR
      K1 = 1 + (L-1)*IROW
      K2 = K1 + IROW - 1
      IF (L.EQ.LR) K2 = NQ
C
C     print *, K1,K2
C
      write(98, 110) (ISPEC(K),K=K1,K2)
 110  FORMAT(/5X,7(A8,1X))     !THIS FORMAT STATEMENT IS THE CULPRIT
C
      write(98, 120) (USOL(K),K=K1,K2)!SPLITS USOL INTO SEPARATED LINES
C     write(98, 120) (USOL(K),K=K1,K2)
C
 120  FORMAT(1X,1P7E9.2)
C     IF (N.EQ.0) GO TO 8             !N = 0, SKIPS TP AND TL PRINT
C
      write(98, 140)
 140  FORMAT(1X,'TP, TL')
      write(98, 145) (TP(K),K=K1,K2)
      write(98, 145) (TL(K),K=K1,K2)
 145  FORMAT(1X,1P7E9.2)
   8  CONTINUE
C
C ***** PRINT REACTION RATES IN MIFprintout.dat *****
      write(98, 179)
  179 FORMAT(/1X,'INTEGRATED REACTION RATES'/)
C
      write(98, 181)
      IROW = 10
      LR = NR/IROW + 1
      RL = FLOAT(NR)/IROW + 1
      DIF = RL - LR
      IF (DIF.LT.0.001) LR = LR - 1
C
      DO 17 L=1,LR
      K1 = 1 + (L-1)*IROW
      K2 = K1 + IROW - 1
C
      IF (L.EQ.LR) THEN
        K2 = NR
        write(98, 186) K1,(RAT(K),K=K1,K2),K2
  186   FORMAT(I3,2X,1P8E10.3,22X,I3) !EDIT 22X FOR NUMER OF RXNS
        GO TO 17
      ENDIF
      write(98, 180) K1,(RAT(K),K=K1,K2),K2
  180 FORMAT(I3,2X,1P10E10.3,2X,I3)
   17 CONTINUE
      write(98, 181)
  181 FORMAT(9X,'1',9X,'2',9X,'3',9X,'4',9X,'5',9X,'6',9X,'7',9X,
     2    '8',9X,'9',8X,'10')     
C
C ***** PRINT ON LAST ITERATION ONLY *****
      write(98, 125)
 125  FORMAT(/1X,'NUMBER DENSITIES OF LONG-LIVED SPECIES'/)

      DO 1 K=1,NQ
      D(K) = USOL(K)*DEN
   1  CONTINUE

      write(98, 119), D
 119  FORMAT(1X,1P7E9.2)

      RETURN
      END

C      SUBROUTINE OUTPUT(USOL,N,NSTEPS,TIME)
C      PARAMETER (NQ=14, NR=28, NSP=16, NMAX=10)
C      DIMENSION D(NQ)
C      DIMENSION USOL(NQ)
C      COMMON/NBLOK/LO,LO2,LO3,LH,LH2,LOH,LH2O,LHO2,LH2O2,
C     2  LSiO,LSiO2,LSiOH,LMg,LMgO,LMgSiO3,LMg2SiO4
C
C   THIS SUBROUTINE PRINTS OUT ALL THE DATA. THE VARIABLE ISKIP 
C   SAYS HOW MANY POINTS YOU WANT TO LOOK AT.
C
C   SET ISKIP=1 FOR MIXING RATIO USOL AT EVERY ITERATION
C      ISKIP = 4
C      IF(ISKIP.GT.1 .AND. N.EQ.NSTEPS) ISKIP = 2
C      TIMEY = TIME/3600./24./365.25
C      write(98, 100) TIME,TIMEY
C TIMEY--time in 
C 100  FORMAT(/1X,'TIME =', 1PE9.2,5X,'TIMEY =',E9.2,1X,'YEARS')
C      write(98, 101) USOL
C 101  FORMAT(/1X,'USOL =',1PE10.3)
C
C     write(98, 105)
C105  FORMAT(/1X,'MIXING RATIOS OF LONG-LIVED SPECIES'/)
C     DO 8 L=1,LR
C     K1 = 1 + (L-1)*IROW
C     K2 = K1 + IROW - 1
C     IF (L.EQ.LR) K2 = NQ
C     write(98, 110) (ISPEC(K),K=K1,K2)
C110  FORMAT(/5X,'Z',7X,13(A8,1X))
C     DO 20 I=1,3
C 20  write(98, 120) Z(I),(USOL(K,I),K=K1,K2)
C     DO 21 I=4,NZ,ISKIP
C 21  write(98, 120) Z(I),(USOL(K,I),K=K1,K2)
C120  FORMAT(1X,1P13E9.2)
C     IF (N.EQ.0) GO TO 8
C     write(98, 140)
C140  FORMAT(/1X,'TP, TL')
C     write(98, 145) (TP(K),K=K1,K2)
C     write(98, 145) (TL(K),K=K1,K2)
C145  FORMAT(10X,1P12E9.2)
C  8  CONTINUE
C
C ***** PRINT ON LAST ITERATION ONLY *****
C      write(98, 125)
C 125  FORMAT(/1X,'NUMBER DENSITIES OF SPECIES'/)

C      DO 1 K=1,NQ
C      D(K) = USOL(K)*DEN
C   1  CONTINUE

C      RETURN
C      END

C-AD **********************************
      SUBROUTINE DOCHEM(USOL,FVAL,N)                                             
      PARAMETER(NQ=13, NR=27, NSP=15, NMAX=10)                                           
      DIMENSION FVAL(NQ),D(NQ),USOL(NQ)
      COMMON/ZBLOK/RAT(NR),TP(NQ),TL(NQ)
C TP--total production, TL--totl loss (indexed by species)
      COMMON/NBLOK/LO,LO2,LO3,LH,LH2,LOH,LH2O,LHO2,LH2O2,
     2  LSiO,LSiO2,LSiOH,LMg,LMgO,LMgSiO3,LMg2SiO4
      COMMON/RBLOK/A(NR),ILOSS(2,NSP,NMAX),IPROD(NSP,NMAX),
     2  JCHEM(5,NR),NUML(NSP),NUMP(NSP),ISPEC(NSP),T,DEN
C                                                         
C   THIS SUBROUTINE DOES THE CHEMISTRY BY CALLING CHEMPL.
C   THESE MUST CONTAIN NO NONLINEARITIES AND MUST BE DONE 
C   IN THE PROPER ORDER (I.E. IF SPECIES A REACTS TO FORM B,
C   THEN A MUST BE FOUND FIRST). LONG-LIVED SPECIES CAN BE 
C   DONE IN ANY ORDER.

C *****SHORT-LIVED SPECIES CHEMISTRY*****
      I = 15 !species 15 is MgSiO3 (no loop!)
      CALL CHEMPL(D,XP,XL,I)
      D(I)=XP/XL

C ***** LONG-LIVED SPECIES CHEMISTRY ***** 
      DO 4 I=1,NSP
   4  D(I) = USOL(I)*DEN
C
      DO 5 I=1,NQ
      CALL CHEMPL(D,XP,XL,I)
      FVAL(I) = XP/DEN - XL*USOL(I)
C line above should work, XP, XL, DEN are 1-D b/c no Z-term
      TP(I) = XP !prod rate (mol/cm^3/sec)
      TL(I) = XL*D(I) !loss rate (mol/cm^3/sec)
   5  CONTINUE
C     IF (N.LT.1) RETURN

C *****CALCULATE MgO INPUT ******

C MADD - MgO ADDITION variable (c.f Lightning / V-outgassing)
      
      MADD =  1.000E10 ! mol / cm^3 sec^-1 (density)
      FVAL(LMgO) = FVAL(LMgO) + MADD/DEN
      TP(LMgO) = TP(LMgO) + MADD

C *****CALCULATE H2O SINK  ******
C HLOSS - H2O LOSS variable (c.f.Rainout)

      HLOSS = 1.500E3 !units of sec^-1 
      TL(LH2O) = TL(LH2O) + HLOSS
      FVAL(LH2O) = TP(LH2O)/DEN - HLOSS*USOL(LH2O)

C ***** CALCULATE PRODUCTION AND LOSS *****
      DO 10 L=1,NR
  10  RAT(L) = 0.
C
      DO 12 L=1,NR
      M = JCHEM(1,L)
      K = JCHEM(2,L)
C 12  RXTOT(L) = A(L)*D(M)*D(K) !*USOL(M)*USOL(K) 
  12  RAT(L) = RAT(L) + A(L)*D(M)*D(K)
C
C      DO 8 I=1,NQ
C      TP(I) = TP(I) + YP(I)
C      TL(I) = TL(I) + YL(I)*D(I)
C   8  CONTINUE
C PROB DON'T NEED STUFF ABOVE, KEPT FOR FUN

      RETURN

      END
C-PK *******************************
      SUBROUTINE CHEMPL(D,XP,XL,K)
      PARAMETER(NQ=13, NR=27, NSP=15, NMAX=10)
      DIMENSION D(NQ)
      COMMON/RBLOK/A(NR),ILOSS(2,NSP,NMAX),IPROD(NSP,NMAX),
     2  JCHEM(5,NR),NUML(NSP),NUMP(NSP),ISPEC(NSP),T,DEN
C
C   THIS SUBROUTINE CALCULATES CHEMICAL PRODUCTION AND LOSS RATES
C   USING THE INFORMATION IN THE MATRICES JCHEM, ILOSS, AND IPROD.
C   CALLED BY SUBROUTINE DOCHEM.
C   
C   XL and XP have no dimenions, because no Z(height) term
C
      XL = 0. !loss frequency 
      XP = 0. !prod rate 
C
C   LOSS FREQUENCY XL (RATE)
      NL = NUML(K)
      DO 2 L=1,NL
      J = ILOSS(1,K,L) !3d matrix loss
      M = ILOSS(2,K,L)
   2  XL = XL + A(J)*D(M) !loss frequency (sec-1)
C
C   PRODUCTION FREQUENCY XP (RATE)
      NP = NUMP(K)
      DO 3 L=1,NP
      J = IPROD(K,L)
      M = JCHEM(1,J)
      N = JCHEM(2,J)
   3  XP = XP + A(J)*D(M)*D(N) 
C XP--production rate b/c based on number densities
C XP: mol/cm^3/sec
C
c     print 748, XP, XL
c748  FORMAT (1P2E10.3)
      
      RETURN
      END



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

