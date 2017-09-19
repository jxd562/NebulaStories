         PROGRAM ARC_CH4 

C-JK Modified on Jun 15, 2017 for mass independent fractionation of oxygen by
C- Alexander Dimoff and Jasmeet K Dhaliwal. 
C Deleted parts of original photochemical code that are not needed to 
C perform the box model of the solar nebula.
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
C     N = TOTAL NUMBER DENSITY (DEN)
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
C ***** FIND REACTIONS IN CHEM.DAT.CH4_OMIF2 *****
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
C    27)  Mg + O = MgO              !check this
C    28)  MgO + MgSiO3 = Mg2SiO4    !check this
C
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
      PARAMETER(NQ=13, NR=28, NSP=15, NMAX=10)
      DIMENSION FVAL(NQ),FV(NQ),DJAC(NQ,NQ),RHS(NQ),REL(NQ),IPVT(NQ)
     2  ,USAVE(NQ),USOL(NQ),TP(NQ),TL(NQ),SPEC(NSP)
      DIMENSION RAT(NR)
      CHARACTER*30 CHEM(5,NR),PRODRX(NSP,NR),LOSSRX(NSP,NR)
      CHARACTER*12 FILSPEC
      CHARACTER*4,DIRDATA
     
      COMMON/NBLOK/LO,LO2,LO3,LH,LH2,LOH,LH2O,LHO2,LH2O2,
     2  LSiO,LSiO2,LSiOH,LMg,LMgO,LMgSiO3,LMg2SiO4
      COMMON/RBLOK/AR(NR),ILOSS(2,NSP,NMAX),IPROD(NSP,NMAX),
     2  JCHEM(5,NR),NUML(NSP),NUMP(NSP)
      COMMON/DBLOK/D(NQ)
C-AP *************************************************************

      DATA LO,LO2,LH,LOH,LH2O,LHO2,LH2O2,
     2  LSiO,LSiO2,LSiOH,LMgO,LMgSiO3,LMg2SiO4,LO3,LH2/
     3  1,2,3,4,5,6,7,8,9,10,11,12,13,14,15/
C      DATA LO,LO2,LO3,LH,LH2,LOH,LH2O,LHO2,LH2O2,
C     2  LSiO,LSiO2,LSiOH,LMg,LMgO,LMgSiO3,LMg2SiO4/
C     3  1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16/
C
C ***** SPECIES DEFINITIONS *****
C   LONG-LIVED - ALL OUR SPECIES ARE "LONG-LIVED"
C      ISPEC(1) = 1HO
C      ISPEC(2) = 2HO2
C SHORT-LIVED
C      ISPEC(3) = 2HO3
C
C      ISPEC(4) = 1HH
C INERT
C      ISPEC(5) = 2HH2
C
C      ISPEC(6) = 2HOH
C      ISPEC(7) = 3HH2O
C      ISPEC(8) = 3HHO2
C      ISPEC(9) = 4HH2O2
C      ISPEC(10) = 3HSiO
C      ISPEC(11) = 4HSiO2
C      ISPEC(12) = 4HSiOH
C      ISPEC(13) = 2HMg
C      ISPEC(14) = 3HMgO
C   ROCK SPECIES (LONG-LIVED)
C      ISPEC(15) = 6HMgSiO3
C      ISPEC(16) = 7HMg2SiO4

C ****************************************************************
      DIRDATA='DATA'
C Input files
      OPEN(UNIT=1,FILE='species_matrix.dat')

C
C ORIGINAL FILES
C
      OPEN(UNIT=7,FILE='species+T_in.dat')

      OPEN(UNIT=9,FILE=DIRDATA//'/CHEM.DAT.CH4_OMIF2')
C
C Output files
      OPEN(UNIT=98,FILE='Photoprintout.dat')
C
      OPEN(UNIT=8,FILE='species+T_out.dat')
C
      OPEN(UNIT=15,FILE='OUT/int.rates.out.dat')
C
C ***** READ THE CHEMISTRY DATA CARDS *****
      READ(9,200)JCHEM
 200  FORMAT(10X,A8,2X,A8,2X,A8,2X,A8,2X,A8)
C      print 201,(J,(JCHEM(M,J),M=1,5),J=1,NR)
C 201  FORMAT(1X,I3,1H),5X,A8,4H +  ,A8,7H  =    ,A8,4H +  ,A8,4X,A8)

C ***** READ SPECIES FILE ****

C      READ(1,*)SPEC,USOL

      READ(1,201)ISPEC
 201  FORMAT(A9)

      PRINT*,'ISPEC'
      PRINT 201,ISPEC

C      READ(1,*)FILSPEC,USOL

C      READ(1,203)USOL
C 203  FORMAT(1P1E10.3)

C
C ***** REPLACE HOLLERITH LABELS WITH SPECIES NUMBERS IN JCHEM *****
C EQ = equal, NE = not equal

      DO 5 J=1,NR
      DO 5 M=1,5 
C NR - number of rxns, M to 5 b/c max prod / reactant species

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

C      print 401,(J,(JCHEM(M,J),M=1,5),J=1,NR)
C 401  FORMAT(1X,I3,1H),5X,I2,4H +  ,I2,7H  =    ,I2,4H +  ,I2,4X,I2)

C Read character array for P&L tables, "int.rates.out.dat"
C      REWIND 9
C      READ(9,200)CHEM
C
C ***** FILL UP CHEMICAL PRODUCTION AND LOSS MATRICES *****
C LT = less than, OR = or, GT = greater than

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
      print *
      DO I=1,NSP
      print 689,I,(IPROD(I,K),K=1,NMAX)
 689  FORMAT(1X,I3,3X,10(I2,1X))
      ENDDO
C
      print *,' '
      DO I=1,NSP
      print 689,I,(ILOSS(1,I,K),K=1,NMAX)
      ENDDO
C
      print *,' '
      DO I=1,NSP
      print 689,I,(ILOSS(2,I,K),K=1,NMAX)
      ENDDO

C
C ***** READ THE INPUT DATAFILE species+T_in.dat *****
      print *,' '
      READ(7,500) USOL,T
 500  FORMAT(1P1E10.3)
C Reading exponential 10 characters, 3 decimal points
C Format 1P1 allows for leading zero
C
      print 1899,USOL,T
 1899 FORMAT(1X,1P1E10.3)
C species+T_in.dat: mixing ratios of species based on
C Lodders (2003) and total number density 10^15
C T at end of file

C 
      DEN = 1.00e15
C      DO I=1,NQ
C      DEN = DEN + USOL(I) !only if USOL in mol/cm^3
C      ENDDO

C     print 1899,USOL,T   !The read works, print for proofies
C1899 FORMAT(1X,1P2E10.1)
C 
C ***** SET MODEL PARAMETERS *****
C     DT = INITIAL TIME STEP
C     TSTOP = TIME AT WHICH CALCULATION IS TO STOP
C     NSTEPS = NUMBER OF TIME STEPS TO RUN (IF TSTOP IS NOT REACHED)
C 
      EPSJ = 1.E-7

C-AD ****** OUTPUT FILE HEADER ******
c-as   print the head of the output file
        write(98, 700) 
  700   format('********************************************',
     & /2x,'OUTPUT FOR BOX MODEL OF SOLAR NEBULA',f6.2,/,
     & '********************************************')
      write(98, 202) NQ!,KJAC
 202  FORMAT(//1X,'NQ = ',I2)!,5X,'KJAC=',I7  
C
C-AD *****************************

      CALL RATES(T,DEN) !DEFINE TEMPERATURE AND NUMBER DENSITY
      NSTEPS = 1
C ***** PRINT OUT INITIAL DATA *****
      CALL OUTPUT(0,NSTEPS,0.)
C
C   SET JACOBIAN PARAMETERS
      KD = 2*NQ + 1
      KU = KD - NQ
      KL = KD + NQ
C
      NPR = NSTEPS
      PRN = NPR
C
C ***** START THE TIME-STEPPING LOOP *****
      TIME = 0.
      NORDER = 1
      STEPCOUNT = 0
      TSTOP = 1.E17

      print *,'nsteps =',NSTEPS

      DO 1 N=1,NSTEPS

      STEPCOUNT = STEPCOUNT + 1

      IF(NORDER.EQ.1) THEN
     	DTINV = 1./DT

      ELSE IF(NORDER.EQ.2) THEN
      	DTINV = (3./2.)/DT
      END IF

      TIME = TIME + DT
C      TOTTIME(STEPCOUNT) = TIME	
      print *,'TIME =', TIME
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

      print *

      print *,'after 1st DOCHEM call, D='
      print 1100,D
 1100 FORMAT(1P10E10.3)
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
      write(99,152) N,RCOND,DT,USOL
 152  FORMAT(I4,2X,1P16E22.15,2X,1P16E22.15,2x,1P16E22.15)

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
 501  FORMAT(1P2E10.3)

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
      PARAMETER NQ=14
      PARAMETER(NR=28, NSP=16, NMAX=7)
      COMMON/RBLOK/A(NR),ILOSS(2,NSP,NMAX),IPROD(NSP,NMAX),
     2  JCHEM(5,NR),NUML(NSP),NUMP(NSP)
C
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
      A(27) = 1.0E-12           !MADE-UP
      A(28) = 6.0E-34           !MADE-UP
C   4  CONTINUE

C      RETURN
C      END
      print *, 'DEN =', DEN
C ***** THREE-BODY COEFFICIENTS *****
C     A(I) = TBDY(K0,KI,N,M,T,DEN)
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
C
C   3  CONTINUE

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
      SUBROUTINE OUTPUT(USOL,NSTEPS,TIME)
      PARAMETER NQ = 16
      PARAMETER (NR = 28, NSP = 16, NMAX = 7)
      DIMENSION TP(NSP),TL(NSP),SL(NSP)
      DIMENSION USOL(NQ)
      COMMON/NBLOK/LO,LO2,LO3,LH,LH2,LOH,LH2O,LHO2,LH2O2,
     2  LSiO,LSiO2,LSiOH,LMg,LMgO,LMgSiO3,LMg2SiO4
      COMMON/DBLOK/D(NQ)
C
C   THIS SUBROUTINE PRINTS OUT ALL THE DATA. THE VARIABLE ISKIP 
C   SAYS HOW MANY POINTS YOU WANT TO LOOK AT.
C
C   SET ISKIP=1 FOR MIXING RATIO (NUMBER DENSITY) AT EVERY ITERATION
      ISKIP = 5
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
C     SL(K) = USOL(K)/(DEN - DEN*USOL(K)) !trying this 
      D(K) = USOL(K)*DEN !original, should be D(I)?
   1  CONTINUE

      RETURN
      END

C-AD **********************************
      SUBROUTINE DOCHEM(FVAL,N)                           
      PARAMETER NQ=14                  
      PARAMETER(NR=28, NSP=16, NMAX=7)                                           
      DIMENSION FVAL(NQ),YP(NQ),YL(NQ),USOL(NQ)
      DIMENSION RXTOT(NR),RAT(NR)
      COMMON/NBLOK/LO,LO2,LO3,LH,LH2,LOH,LH2O,LHO2,LH2O2,
     2  LSiO,LSiO2,LSiOH,LMg,LMgO,LMgSiO3,LMg2SiO4
      COMMON/RBLOK/A(NR),ILOSS(2,NSP,NMAX),IPROD(NSP,NMAX),
     2  JCHEM(5,NR),NUML(NSP),NUMP(NSP)
      COMMON/DBLOK/D(NQ)
C                                                         
C   THIS SUBROUTINE DOES THE CHEMISTRY BY CALLING CHEMPL.  PHOTO-         
C   CHEMICAL EQUILIBRIUM SPECIES ARE DONE FIRST.  THESE MUST CON-         
C   TAIN NO NONLINEARITIES AND MUST BE DONE IN THE PROPER ORDER 
C   (I.E. IF SPECIES A REACTS TO FORM B, THEN A MUST BE FOUND FIRST).
C   LONG-LIVED SPECIES CAN BE DONE IN ANY ORDER.
C
C ***** LONG-LIVED SPECIES CHEMISTRY ***** IN JAC MATRIX
      DO 4 I=1,NQ
      D(I) = USOL(I)*DEN !NUMBER DENSITY
C converting to number density above using mixing ratio
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
      PARAMETER(NR=28, NSP=16, NMAX=7)
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
c      SUBROUTINE ISOTOPE(NQ,NR,NSP)
c      PARAMETER (NQ=25,NR=,NSP=29,NMAX=10)

c      COMMON/NBLOK/LO,LXO,LO2,LOXO,LO3,LO2XO,LH,LH2,LOH,LXOH,LH2O,LH2XO,
c     2  LHO2,LHOXO,LH2O2,LH2OXO,LSiO,LSiXO,LSiO2,LSiOXO,LSiOH,LSiXOH,
c     3  LMg,LMgO,LMgXO,LMgSiO3,LMgSiO2XO,LMg2SiO4,LMg2SiO3XO

c      DATA LO,LXO,LO2,LOXO,LO3,LO2XO,LH,LH2,LOH,LXOH,LH2O,LH2XO,
c     2  LHO2,LHOXO,LH2O2,LH2OXO,LSiO,LSiXO,LSiO2,LSiOXO,LSiOH,LSiXOH,
c     3  LMg,LMgO,LMgXO,LMgSiO3,LMgSiO2XO,LMg2SiO4,LMg2SiO3XO/
c     4  1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23
c     5  24,25,26,27,28,29/
C
C ***** SPECIES DEFINITIONS *****
C   LONG-LIVED - ALL OUR SPECIES ARE "LONG-LIVED"
c      ISPEC(1) = 1HO
c      ISPEC(2) = 2HXO
c      ISPEC(3) = 2HO2
c      ISPEC(4) = 3HOXO
c      ISPEC(5) = 2HO3
c      ISPEC(6) = 4HO2XO
c      ISPEC(7) = 1HH
c     ISPEC(8) = 2HH2
c      ISPEC(9) = 2HOH
c      ISPEC(10) = 3HXOH
c      ISPEC(11) = 3HH2O
c      ISPEC(12) = 4HH2XO
c      ISPEC(13) = 3HHO2
c      ISPEC(14) = 4HHOXO
c      ISPEC(15) = 4HH2O2
c      ISPEC(16) = 5HH2OXO
c      ISPEC(17) = 3HSiO
c      ISPEC(18) = 4HSiXO
c      ISPEC(19) = 4HSiO2
c      ISPEC(20) = 5HSiOXO
c      ISPEC(21) = 4HSiOH
c      ISPEC(22) = 5HSiXOH
c      ISPEC(23) = 2HMg
c      ISPEC(24) = 3HMgO
c      ISPEC(25) = 4HMgXO
C   ROCK SPECIES (INERT SPECIES)
c     ISPEC(26) = 6HMgSiO3
c     ISPEC(27) = 8HMgSiO2XO
c      ISPEC(28) = 7HMg2SiO4
c      ISPEC(29) = 9HMg2SiO3XO
C
C ***** NEW REACTION LIST *****
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

