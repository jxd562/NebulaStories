         PROGRAM ARC_CH4

C-AD Modified on Jun, 2017 for mass independent fractionation of oxygen by
C-AD Alexander Dimoff. Used structure of original photochemical code that
C-AD are needed to perform the box model of the solar nebula, and added
C    remaining subroutines and chemistry for our reactions.
C
C    The overall goal is to get this code running and converging with good
C    numbers, and then adapt to an isotopic environment to test fractionation.
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
C
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
C ***** FIND REACTION NETWORK IN CHEM.DAT.CH4_OMIF *****
C     1)  SiO + O2 = SiO2 + O
C     2)  SiO + H2 = SiOH + H
C     3)  SiO + OH = SiO2 + H
C     4)  SiO + HO2 = SiO2 + OH
C     5)  SiO + O = SiO2
C     6)  SiO + H = SiOH
C     7)  SiO + O3 = SiO2 + O2
C     8)  H2 + O = OH + H
C     9)  OH + O = O2 + H
C    10)  OH + OH = H2O + O
C    11)  H + O = OH
C    12)  H + O2 = OH + O
C    13)  O + HO2 = OH + O2
C    14)  O + H2O2 = OH + HO2
C    15)  OH + H2 = H2O + H
C    16)  OH + HO2 = H2O + O2
C    17)  OH + H2O2 = H2O + HO2
C    18)  HO2 + HO2 = H2O2 + O2
C    19)  OH + O3 = HO2 + O2
C    20)  HO2 + O3 = OH + 2O2
C    21)  SiOH + H = SiO + H2
C    22)  SiOH + OH = SiO + H2O
C
C ***** THREE BODY RXN *****
C    23)  O + O2 + M = O3 + M 
C    24)  H + O2 + H2O = HO2 + H2O
C    25)  H + O2 + H2 = HO2 + H2
C    26)  H + H + H2 = H2 + H2
C    27)  H + OH + H2O = H2O + H2O
C
C ***** ADDED REACTIONS JUN 23, 2017 *****
C    28)  MgO + SiO2 = MgSiO3
C    29)  MgO + MgSiO3 = Mg2SiO4
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
C               IS SOLVED AND WHICH ARE IN THE BIG, REVERSE EULER MATRIX
C        NMAX = MAXIMUM NUMBER OF REACTIONS IN WHICH AN INDIVIDUAL
C               SPECIES PARTICIPATES
C
C     THREE-BODY REACTIONS ARE WRITTEN IN TWO-BODY FORM, SO THE DENSITY 
C     FACTOR MUST BE INCLUDED IN THE RATE CONSTANT.
C
C     THREE-BODY REACTIONS THAT DO NOT FALL INTO THAT CATEGORY (EX. En, Fo)
C     ARE CALCULATED SEPARATELY, ACCORDING TO THEIR REACTION KINEMATICS.
C
C     THE CHEMICAL REACTION SCHEME IS STORED IN THE FOLLOWING MATRICES
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
      INCLUDE 'header.inc'
      DIMENSION FVAL(NQ),FV(NQ),DJAC(NQ,NQ),RHS(NQ),REL(NQ),IPVT(NQ)
     2  ,USAVE(NQ),USOL(NQ),D(NSP),INDX(NQ)
      CHARACTER*30 CHEM(5,NR),PRODRX(NSP,NR),LOSSRX(NSP,NR)
      CHARACTER*4,DIRDATA
     
      COMMON/NBLOK/LO,LO2,LH,LOH,LH2O,LHO2,LH2O2,
     2  LSiO,LSiO2,LSiOH,LMgO,LO3,LH2,LMgSiO3,LMg2SiO4
      COMMON/RBLOK/A(NR),ILOSS(2,NSP,NMAX),IPROD(NSP,NMAX),
     2  JCHEM(5,NR),NUML(NSP),NUMP(NSP),ISPEC(NSP),T,DEN
      COMMON/ZBLOK/RAT(NR),TP(NQ),TL(NQ),SL,H2,DH2,USOLTOT,
     2  SADD,AMADD,OADD,CON(NQ),CONO,CONSi,CONMg,SiO2SMR,
     3  EnSMR,FoSMR,PRESS,SiO2COND,EnCOND,FoCOND,H2OLOSS
C
C-AP *************************************************************

      DATA LO,LO2,LH,LOH,LH2O,LHO2,LH2O2,
     2  LSiO,LSiO2,LSiOH,LMgO,LMgSiO3,LMg2SiO4,LO3,LH2/
     3  1,2,3,4,5,6,7,8,9,10,11,12,13,14,15/
C
C ***** SPECIES DEFINITIONS *****
C   LONG-LIVED SPECIES
      ISPEC(1) = 1HO
      ISPEC(2) = 2HO2
      ISPEC(3) = 1HH
      ISPEC(4) = 2HOH
      ISPEC(5) = 3HH2O
      ISPEC(6) = 3HHO2
      ISPEC(7) = 4HH2O2
      ISPEC(8) = 3HSiO
      ISPEC(9) = 4HSiO2
      ISPEC(10) = 4HSiOH
      ISPEC(11) = 3HMgO
C   ROCK SPECIES (LONG-LIVED)
      ISPEC(12) = 6HMgSiO3
      ISPEC(13) = 7HMg2SiO4
C   SHORT-LIVED SPECIES
      ISPEC(14) = 2HO3
C    INERT SPECIES
      ISPEC(15) = 2HH2
C
C ****************************************************************
      DIRDATA='DATA'
C Input files
      OPEN(UNIT=7,FILE='species+T_in.dat')
C
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
C
C     print 201,(J,(JCHEM(M,J),M=1,5),J=1,NR)
C201  FORMAT(1X,I3,1H),5X,A8,4H +  ,A8,7H  =    ,A8,4H +  ,A8,4X,A8)
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
C
C ***** Read character array for P&L tables for "int.rates.out.dat"
      REWIND 9
      READ(9,200)CHEM
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
C     print *,'PROD'
C     DO I=1,NSP
C     print 689,I,(IPROD(I,K),K=1,NMAX)
C689  FORMAT(1X,I3,3X,12(I2,1X))
C     ENDDO
C
C     print *,'LOSS 1'
C     DO I=1,NSP
C     print 689,I,(ILOSS(1,I,K),K=1,NMAX)
C     ENDDO
C
C     print *,'LOSS 2'
C     DO I=1,NSP
C     print 689,I,(ILOSS(2,I,K),K=1,NMAX)
C     ENDDO
C
C ***** READ THE INPUT DATAFILE species+T_in.dat *****
      READ(7,500) USOL,T
 500  FORMAT(1P1E10.3)
C
C ***** TOTAL DENSITY *****
      DEN = 1E15 !5.17598E14 FROM P=nKT
C
C ***** CALCULATE INDIVIDUAL NUMBER DENSITY OF LONG-LIVED SPECIES*****
      DO I=1,NQ
      D(I) = USOL(I)*(DEN)
      ENDDO
C
C ***** SET MODEL PARAMETERS *****
C     DT = INITIAL TIME STEP
C     TSTOP = TIME AT WHICH CALCULATION IS TO STOP
C     NSTEPS = NUMBER OF TIME STEPS TO RUN (IF TSTOP IS NOT REACHED)
C     EPSJ = PERTURBATION AMOUNT FOR JACOBIAN CALCULATION
C
      EPSJ = 1.E-7 
C
C-AD ****** OUTPUT FILE HEADER ******
C-AD  PRINT THE HEAD OF THE OUTPUT FILE, MIFprintout
       write(98, 700) 
 700   FORMAT('********************************************',
     & /2x,'OUTPUT FOR BOX MODEL OF SOLAR NEBULA'/,
     & '********************************************')
      write(98, 202) NQ 
 202  FORMAT(/1X,'NQ = ',I2/)
C
C-AD ***** CALL RATES FOR RATE MATRIX *****
      CALL RATES(D)
C
C ***** PRINT OUT INITIAL DATA *****
C     CALL OUTPUT(USOL,0,NSTEPS,0.)
C  ALL PRINTED VALUES SHOULD BE ZERO,
C  NOTHING HAS BEEN DEFINED OR CALCULATED YET.
C
C ***** START THE TIME-STEPPING LOOP *****
      TIME = 0.
      DT = 1.E-8    !SIZE OF TIME STEP
      STEPCOUNT = 0
      TSTOP = 1.E17
      NSTEPS = 1000
C
      DO 1 N=1,NSTEPS
C
      STEPCOUNT = STEPCOUNT + 1
C
      DTINV = 1./DT
C
      TIME = TIME + DT
C-AD
C ***** SET UP THE JACOBIAN MATRIX AND RIGHT-HAND SIDE *****
      DO 17 J=1,NQ
      DO 17 K=1,NQ
  17  DJAC(J,K) = 0.
      DO 19 K=1,NQ
  19  RHS(K) = 0.
C
C     DJAC = (1/DT)*I - J, WHERE J IS THE JACOBIAN MATRIX
C
C ***** COMPUTE CHEMISTRY TERMS *****
C
      IDO = 0
      IF (NN.EQ.NSTEPS) IDO = 1
      CALL DOCHEM(USOL,FVAL,IDO)
C
C     print *, 'D = '
C     print 155, D
C155  FORMAT (1P1E10.3)
C
C     print *
C     print *, 'FVAL'
C     print 668, FVAL
C668  FORMAT (1P7E10.1)
C
      DO 9 I=1,NQ
      RHS(I) = FVAL(I)
   9  USAVE(I) = USOL(I)
C
      DO 3 I=1,NQ
      R = EPSJ * ABS(USOL(I))
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
C     print *
C     print *, 'DJAC BEFORE MATRIX SOLVER'
C     print 666, ((DJAC(I,J),J=1,NQ,I=1,NQ)
C666  FORMAT (1P14E8.1)
C   **DJAC PRINTS WELL**
C
C     print *, 'D ='
C     print 157, D
C157  FORMAT (1P1E10.3) 
C   **D VALUES LOOK GOOD BEFORE MATRIX SOLVER**
C
C ***** FACTOR THE JACOBIAN AND SOLVE THE LINEAR SYSTEM *****
C 
      CALL LUDCMP(DJAC,NQ,NQ,INDX,DP)  !N, NP GOES TO NQ, NQ,, D -> DP
      CALL LUBKSB(DJAC,NQ,NQ,INDX,RHS) !N, NP GOES TO NQ, NQ

C ***** COMPUTE NEW CONCENTRATIONS *****
      EMAX = 0.
      DO 26 I=1,NQ
C
      REL(I) = RHS(I)/USOL(I)
      EREL = ABS(REL(I))
      EMAX = AMAX1(EMAX,EREL)
C
      IF(EREL.LT.EMAX) GO TO 26
C
      ISP = ISPEC(I)
      UMAX = USOL(I)
      RMAX = RHS(I)
  26  USOL(I) = USOL(I) + RHS(I) ! X = X + DX
C
C
C ***** AUTOMATIC TIME STEP CONTROL *****
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
      write(98, 100) N,EMAX,ISP,UMAX,RMAX,DT,TIME
 100  FORMAT(1X,'N =',I4,2X,'EMAX =',1PE9.2,' FOR ',
     2  A8,'U =',1PE9.2,1X,'RHS =',1PE9.2,
     3  2X,'DT =',1PE9.2,2X,'TIME =',1PE9.2)
      CONTINUE
C
      IF (EMAX.LT.0.5) GO TO 28
      print *, 'Decreasing Step'
      DT = 0.5*DTSAVE     
      TIME = TIME - DTSAVE
C
      DO 27 I=1,NQ
  27  USOL(I) = USAVE(I)
  28  CONTINUE
C
      IF(N.EQ.1) CALL OUTPUT(USOL,NN,NSTEPS,TIME)
C
      IF(INDEX.NE.0) STOP
C
      IF(NN.EQ.NSTEPS) GO TO 22
C
      IF(TIME.GT.TSTOP) NN = NSTEPS - 1
C
      IF(TIME.GT.1E18) EXIT !when model runs to older than age of universe, exit loop
C
  22  CONTINUE
  1   CONTINUE
C ***** END THE TIME-STEPPING LOOP *****
C
      CALL OUTPUT(USOL,NN,NSTEPS,TIME)
C
C SHOW THE RATES FOR EACH REACTION 
C     print *
C     print *, 'RAT'
C     print 502,RAT
C502  FORMAT(1P10E10.3)
C
C     print *
C     print *, 'A'
C     print 504,A
C504  FORMAT(1P10E10.3)
C
C ** WRITE ENDING MIXING RATIOS AND TEMPERATURE TO species+T_out.dat
C    can use this to run final results in new iterations
      write(8,501) USOL,T
 501  FORMAT(1P1E10.3)
C
C print out P&L tables with integrated rxn rates, "int.rates.out.dat"
      DO 702 I=1,NSP
         ISP = ISPEC(I)
         WRITE(15,703) ISP,TP(I)
 703     FORMAT(/A8,12X,'PRODUCTION RXS',13X,'INT RX RATE',4X,
     2      'TP = ',1PE9.2)
 
       DO 704 N=1,NR
          IF(JCHEM(3,N).EQ.I .OR. JCHEM(4,N).EQ.I .OR.
     2       JCHEM(5,N).EQ.I)THEN
           IF(RAT(N).NE.0.) WRITE(15,705) N,(CHEM(J,N),J=1,5),RAT(N)
 705       FORMAT(1X,I3,')',1X,A7,' + ',A7,' = ',A7,' + ',A7,2X,A4,
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
C
      STOP
      END PROGRAM ARC_CH4

C-PK ********************************
      SUBROUTINE RATES(D)
      INCLUDE 'header.inc'
      COMMON/RBLOK/A(NR),ILOSS(2,NSP,NMAX),IPROD(NSP,NMAX),
     2  JCHEM(5,NR),NUML(NSP),NUMP(NSP),ISPEC(NSP),T,DEN
C
C ***** FILL UP RATE MATRIX *****
C   A MATRIX DIMENSIONED NR, K COEFFS
C   GO DIRECTLY INTO RATE MATRIX A
C   RATE CONSTANTS TAKEN FROM CHAKORBORDY_13_SUPP_INFO TABLE S5
C
      A(1) = 1.0E-11*EXP(110/T)                 !Assumed Chakorbordy 
      A(2) = 3.0E-12      !edited by AD 7/10    !Assumed Chakorbordy 
      A(3) = 1.0E-09      !edited by AD 7/10    !Assumed Chakorbordy 
      A(4) = 5.0E-12                            !Assumed Chakorbordy 
      A(5) = 1.0E-12                            !Assumed Chakorbordy 
      A(6) = 1.0E-18      !edited by AD 7/10    !Assumed Chakorbordy
      A(7) = 1.0E-12                            !Assumed Chakorbordy 
      A(8) = 8.5E-20*(T**2.67)*EXP(-3160./T)    !HCP from Chakorbordy
      A(9) = 2.4E-11*EXP(-352./T)               !HCP "
      A(10) = 2.5E-15*(T**1.14)*EXP(-50./T)     !HCP "
      A(11) = 1.0E-12                           !HCP "
      A(12) = 3.3E-10*EXP(-8460./T)             !HCP "
      A(13) = 5.3E-11                           !HCP "
      A(14) = 1.1E-12*EXP(-2000./T)             !HCP "
      A(15) = 1.7E-16*(T**1.6)*EXP(-1660./T)    !HCP "
      A(16) = 4.8E-11*EXP(250./T)               !HCP "
      A(17) = 1.3E-11*EXP(-670./T)              !HCP "
      A(18) = 3.1E-12*EXP(-775./T)              !HCP "
      A(19) = 1.9E-12*EXP(-1000./T)             !HCP "
      A(20) = 1.4E-14*EXP(-600./T)              !HCP "
      A(21) = 1.0E-10                           !MADE-UP, KASTING
      A(22) = 1.0E-10                           !MADE-UP, KASTING
C
C HCP = CRC HANDBOOK OF CHEMISTRY AND PHYSICS
C JPL = CHEMICAL KINETICS AND PHOTOCHEMICAL DATA FOR USE IN ATMOSPHERIC
C       STUDIES EVALUATION NUMBER 15, NASA, JPL (NASA, JPL, 2006)
C
C ***** THREE-BODY REACTION COEFFICIENTS *****
C     A(I) = TBDY(K0,KI,N,M,T,D)
C
C       O + O2 + M = O3 + M                          
      A(23) = TBDY(6.0E-34,1.E-10,2.3,0.,T,DEN) !JPL from Chakorbordy paper
C       H + O2 + H2O = HO2 + H2O
      A(24) = TBDY(4.3E-30,1.E-10,-0.8,0.,T,DEN)!HCP "
C       H + O2 + H2 = HO2 + H2
      A(25) = TBDY(5.8E-30,1.E-10,-0.8,0.,T,DEN)!HCP "
C       H + H + H2 = H2 + H2
      A(26) = TBDY(2.7E-31,1.E-10,-0.6,0.,T,DEN)!HCP "
C       H + OH + H2O = H2O + H2O
      A(27) = TBDY(3.9E-25,1.E-10,-2.,0.,T,DEN) !HCP "
C
C ***** ROCK REACTIONS *****
      A(28) = 1.0E-12                           !MADE-UP, KASTING
      A(29) = 1.0E-10                           !MADE-UP, KASTING
C
      RETURN
      END
C
C-AD ***************************************************
      FUNCTION TBDY(A0,AI,CN,CM,T,D)                      
C     A0 -- LOW PRESSURE LIMIT
C     AI -- HIGH PRESSURE LIMIT
C     CN -- EXPONENT LOW PRESSURE
C     CM -- EXPONENT HIGH PRESSURE
      B0 = A0*(T**CN) !LOW PRESSURE                                
      BI = AI*(T**CM) !HIGH PRESSURE                 
      Y = ALOG10(B0*D/BI)                                 
      X = 1./(1. + Y**2)                                  
      TBDY = B0*D/(1. + B0*D/BI) * 0.6**X
C
      RETURN
      END
C
C-AD ***************************************************
      SUBROUTINE OUTPUT(USOL,N,NSTEPS,TIME)
      INCLUDE 'header.inc'
      DIMENSION USOL(NQ),D(NSP)
      COMMON/NBLOK/LO,LO2,LH,LOH,LH2O,LHO2,LH2O2,
     2  LSiO,LSiO2,LSiOH,LMgO,LO3,LH2,LMgSiO3,LMg2SiO4
      COMMON/RBLOK/A(NR),ILOSS(2,NSP,NMAX),IPROD(NSP,NMAX),
     2  JCHEM(5,NR),NUML(NSP),NUMP(NSP),ISPEC(NSP),T,DEN
      COMMON/ZBLOK/RAT(NR),TP(NQ),TL(NQ),SL,H2,DH2,USOLTOT,
     2  SADD,AMADD,OADD,CON(NQ),CONO,CONSi,CONMg,SiO2SMR,
     3  EnSMR,FoSMR,PRESS,SiO2COND,EnCOND,FoCOND,H2OLOSS
C
C   THIS SUBROUTINE PRINTS OUT ALL THE DATA. THE VARIABLE ISKIP 
C   SAYS HOW MANY POINTS YOU WANT TO LOOK AT.
C   PRINTS TO MIFprintout.dat FILE.
C
C ** SET ISKIP=1 FOR MIXING RATIO USOL AT EVERY ITERATION **
c     ISKIP = 2
c     IF(ISKIP.GT.1 .AND. N.EQ.NSTEPS) ISKIP = 2
C
      TIMEY = TIME/3600./24./365.25
C     
C
      write(98, 100) TIME,TIMEY
 100  FORMAT(/1X,'TIME =', 1PE9.2,5X,'TIMEY =',E9.2,1X,'YEARS')
C
C
C ** WRITE MIXING RATIOS **
      write(98, 105)
 105  FORMAT(/1X,'MIXING RATIOS OF LONG-LIVED SPECIES,')
      write(98, 106)
 106  FORMAT(1X,'TOTAL PRODUCTION AND LOSS, CONSERVATION')
C
      IROW = 7                      !FORMATTING
      LR = NQ/IROW + 1              !LOGIC
      RL = FLOAT(NQ)/IROW + 1       
      DIF = RL - LR                 
      IF (DIF.LT.0.001) LR = LR - 1 
C
      DO 8 L=1,LR
      K1 = 1 + (L-1)*IROW
      K2 = K1 + IROW - 1
      IF (L.EQ.LR) K2 = NQ
C
      write(98, 110) (ISPEC(K),K=K1,K2)
 110  FORMAT(/5X,7(A8,1X))     
C
      write(98, 120) (USOL(K),K=K1,K2)
 120  FORMAT(1X,1P7E9.2)
C
      write(98, 140)
 140  FORMAT(1X,'TP, TL, CONSERV')
      write(98, 145) (TP(K),K=K1,K2)
      write(98, 145) (TL(K),K=K1,K2)
 145  FORMAT(1X,1P8E9.2)
      write(98,147) (CON(I),I=K1,K2)
 147  FORMAT(1X,1P8E9.2)
   8  CONTINUE
C
C
C ** WRITE MIXING RATIO OF H2, SAT. MIXING RATIO FOR SiO2, En **
      USOLH2 = DH2/DEN !SAME AS (1-USOLTOT)
C
      write(98,148)
 148  FORMAT(/1X,'MIXING RATIOS OF H2 AND CONDENSATES')
      write(98,149) USOLH2,SiO2SMR,EnSMR,FoSMR
 149  FORMAT(/3X,'H2MR',6X,'SiO2SMR',3X,'EnSMR',5X,'FoSMR'/1P4E10.2)
C
      write(98,160)
 160  FORMAT(/1X,'CONDENSATION AND INJECTION RATES')
      write(98,161) H2OLOSS,SiO2COND,EnCOND,FoCOND
 161  FORMAT(/3X,'H2OLOSS'3x'SiO2COND',2X,'EnCOND',4X,
     2  'FoCOND'/1X,1P4E10.3)
      write (98,162)
 162  FORMAT(/4X,'OADD',5X,'SADD',5X,'AMADD')
      WRITE(98,163) OADD, SADD, AMADD
 163  FORMAT(1X,1P3E10.3)
C
C
C ** ELEMENTAL MASS CONSERVATION **
      write(98,151)
 151  FORMAT(/1X,'FRACTIONAL ELEMENTAL MASS CONSERVATION')
C
      write(98,152)
 152  FORMAT(/5X,'O',7X,'Si',7X,'Mg')
C
      write(98,153)CONO,CONSi,CONMg
 153  FORMAT(1X,1P3E9.2)
C
C
C ** WRITE REACTION RATES **
      write(98, 179)
 179  FORMAT(/1X,'INTEGRATED REACTION RATES'/)
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
 186    FORMAT(I3,2X,1P9E10.3,12X,I3) !EDIT 12X FOR NUMER OF RXNS
        GO TO 17
      ENDIF
      write(98, 180) K1,(RAT(K),K=K1,K2),K2
 180  FORMAT(I3,2X,1P10E10.3,2X,I3)
   17 CONTINUE
      write(98, 181)
 181  FORMAT(9X,'1',9X,'2',9X,'3',9X,'4',9X,'5',9X,'6',9X,'7',9X,
     2    '8',9X,'9',8X,'10')
C
C
C ** WRITE NUMBER DENSITIES **
      write(98, 125)
 125  FORMAT(/1X,'NUMBER DENSITIES OF LONG-LIVED SPECIES'/)
C
      DO 1 K=1,NQ
      D(K) = USOL(K)*DEN
    1 CONTINUE
C
      IROW = 6
      LR = NQ/IROW + 1         
      RL = FLOAT(NQ)/IROW + 1  
      DIF = RL - LR
      IF (DIF.LT.0.001) LR = LR - 1
C
      DO 7 L=1,LR
      K1 = 1 + (L-1)*IROW
      K2 = K1 + IROW - 1
      IF (L.EQ.LR) K2 = NQ
      write(98, 119) (ISPEC(K),K=K1,K2)
 119  FORMAT(5X,7(A8,1X))     
C
      write(98, 118) (D(K),K=K1,K2)
 118  FORMAT(1X,1P7E9.2)
   7  CONTINUE
C  
      write(98,182)
 182  FORMAT(/1X,'NUMBER DENSITIES OF SHORT-LIVED SPECIES'/)
C
      write(98, 121) ISPEC(LO3)
 121  FORMAT(5X,1(A8,1X))     
C     
      write(98, 127) SL
 127  FORMAT(1X,1PE9.2)
   6  CONTINUE
C
      write(98,184)
 184  FORMAT(/1X,'NUMBER DENSITIES OF INERT SPECIES'/)
C
      write(98, 123) ISPEC(LH2)
 123  FORMAT(5X,1(A8,1X))     
C
      write(98, 129) DH2
 129  FORMAT(1X,1PE9.2/)
C
      RETURN
      END

C-AD **********************************
      SUBROUTINE DOCHEM(USOL,FVAL,N)                                             
      INCLUDE 'header.inc'
      DIMENSION FVAL(NQ),D(NSP),USOL(NQ)
      COMMON/NBLOK/LO,LO2,LH,LOH,LH2O,LHO2,LH2O2,
     2  LSiO,LSiO2,LSiOH,LMgO,LO3,LH2,LMgSiO3,LMg2SiO4
      COMMON/RBLOK/A(NR),ILOSS(2,NSP,NMAX),IPROD(NSP,NMAX),
     2  JCHEM(5,NR),NUML(NSP),NUMP(NSP),ISPEC(NSP),T,DEN
      COMMON/ZBLOK/RAT(NR),TP(NQ),TL(NQ),SL,H2,DH2,USOLTOT,
     2  SADD,AMADD,OADD,CON(NQ),CONO,CONSi,CONMg,SiO2SMR,
     3  EnSMR,FoSMR,PRESS,SiO2COND,EnCOND,FoCOND,H2OLOSS
C                                                         
C   THIS SUBROUTINE DOES THE CHEMISTRY BY CALLING CHEMPL.
C   THESE MUST CONTAIN NO NONLINEARITIES AND MUST BE DONE 
C   IN THE PROPER ORDER (I.E. IF SPECIES A REACTS TO FORM B,
C   THEN A MUST BE FOUND FIRST). LONG-LIVED SPECIES CAN BE 
C   DONE IN ANY ORDER.
C
      PRESS = DEN*1.38E-16*T          !calculate total pressure of system
C                                     !193.200 [DYNES CM^-2] using T=1400 and DEN=1E15
      CONFAC = 1.                     !condensation factor
C
      DO 1 I=1,NQ                                !loop over species
   1  D(I) = USOL(I)*DEN                         !calculate densities
C
      USOLTOT = 0.                               !create USOLTOT
      DO 2 I=1,NQ
   2  USOLTOT = USOLTOT + USOL(I)                !calculate USOLTOT
C
C
C ***** SPECIES CHEMISTRY *****
C THIS SECTION CALCUATES DENSITIES AND PRODUCTION AND LOSS RATES
C
C ***** INERT SPECIES CHEMISTRY ***** (H2 IS ONLY INERT SPECIES)
C FOR MORE INERT SPECIES, MAKE IL A VECTOR
      D(LH2) = (1 - USOLTOT)*DEN
      DH2 = D(LH2)
C
C NO TP/TL FOR H2, FOR IT IS ETERNAL AND CONSTANT   
C
C ***** SHORT-LIVED SPECIES CHEMISTRY ***** (O3 IS ONLY SHORT-LIVED SPECIES)
C FOR MORE SHORT-LIVED SPECIES, MAKE SL A VECTOR
      I = LO3
      CALL CHEMPL(D,XP,XL,I)
      D(I) = XP/XL
      SL = D(I)
      TP(I) = XP
      TL(I) = XL*SL
C
C -- DIAGNOSTIC PRINT STATEMENTS FOR O3 --
C     PRINT 141, D
C 141 FORMAT(1P10E10.3)
C     STOP
C
C ***** LONG-LIVED SPECIES CHEMISTRY ***** 
      DO 3 I=1,NQ
      CALL CHEMPL(D,XP,XL,I)
      FVAL(I) = XP/DEN - XL*USOL(I)
      TP(I) = XP                     !TOTAL PRODUCTION RATE [CM^-3 S^-1]
      TL(I) = XL*D(I)                !TOTAL LOSS RATE       [CM^-3 S^-1]
   3  CONTINUE
C
C
C ***** INJECTION (ENTRY) SOURCES *****
C -- MADD - MgO (Mg) ADDITION RATE VARIABLE
      AMADD = 1.020E11 !EARLY SOLAR NEBULA ABUNDANCE RATIO (FIXED)
      FVAL(LMgO) = FVAL(LMgO) + AMADD/DEN
      TP(LMgO) = TP(LMgO) + AMADD
C
C -- SADD - SiO (Si) ADDITION VARIABLE (FIXED)
      SADD = 1.000E11 !EARLY SOLAR NEBULA ABUNDANCE RATIO (FIXED)
      FVAL(LSiO) = FVAL(LSiO) + SADD/DEN ![S ^-1]
      TP(LSiO) = TP(LSiO) + SADD         ![CM^-3 SEC^-1]
C
C -- OADD - O ADDITION VARIABLE (TUNEABLE)
C     OADD = 1.413E12 !EARLY SOLAR NEBULA ABUNDANCE RATIO
      OADD = 1.413E12
      FVAL(LO) = FVAL(LO) + OADD/DEN
      TP(LO) = TP(LO) + OADD
C
C -- HADD - H2O (O/H) ADDITION VARIABLE
C     HADD = 1.25E14
C     FVAL(LH2O) = FVAL(LH2O) + HADD/DEN
C     TP(LH2O) = TP(LH2O) + HADD
C  WHEN WE USE THIS, DONT FORGET TO USE HADD PARAMETER IN CONSERVATION
C  OUR REACTION NETWORK DOES NOT DISSOCIATE H2O
C
C
C ***** CONDENSATION (EXIT) SOURCES *****
C
C  H2O ADVECTION SINK
C -- H2OLOSS - H2O LOSS VARIABLE
      H2OLFREQ = 1.  !H2O LOSS FREQUENCY [SEC^-1] (TUNABLE PARAM)
      H2OLOSS = H2OLFREQ*D(LH2O) !LOSS RATE [CM^-3 SEC^-1]
C    INCORPORATE H2OLOSS INTO CHEMICAL SCHEME
      TL(LH2O) = TL(LH2O) + H2OLOSS
      FVAL(LH2O) = TP(LH2O)/DEN - H2OLFREQ*USOL(LH2O)
C
C  SiO2 CONDENSATION 
      SiO2VAP = 4.39311833E-03 ![DYNES CM^-2], SAT. VAP. PRES. OF SILICATE
C        !FROM JASMEET'S SAT. VAP. PRES. DATA AT T = 1400 K
      SiO2SMR = SiO2VAP / PRESS !SATURATION MIXING RATIO
      CONDSiO2 = CONFAC * (USOL(LSiO2) - SiO2SMR)
      CONDSiO2 = AMAX1(CONDSiO2,0.)
C
C   INCORPORATE CONSiO2 INTO CHEMISTRY
      FVAL(LSiO2) = FVAL(LSiO2) - CONDSiO2
      TL(LSiO2) = TL(LSiO2) + CONDSiO2*DEN
      TP(LSiO2) = TP(LSiO2)
C
C  MgSiO3 (En) CONDENSTAION 
      EnVAP = 1.03108108E-04 ![DYNES CM^-2]
C        !FROM JASMEET'S SAT. VAP. PRES. DATA
      EnSMR = EnVAP / PRESS !SATURATION MIXING RATIO
      CONDEn = CONFAC * (USOL(LMgSiO3) - EnSMR)
      CONDEn = AMAX1(CONDEn,0.)
C
C   INCORPORATE CONDEn INTO CHEMISTRY
      FVAL(LMgSiO3) = FVAL(LMgSiO3) - CONDEn
      TL(LMgSiO3) = TL(LMgSiO3) + CONDEn*DEN
      TP(LMgSiO3) = TP(LMgSiO3)
C
C  Mg2SiO4 (Fo) CONDENSATION
      FoVAP = 3.339622114E-04 ![DYNES CM^-2]
C        !FROM JASMEET'S SAT. VAP. PRES. DATA
      FoSMR = FoVAP / PRESS !SATURATION MIXING RATIO
      CONDFo = CONFAC * (USOL(LMg2SiO4) - FoSMR)
      CONDFo = AMAX1(CONDFo,0.)
C
C   INCORPORATE CONDFo INTO CHEMISTRY
      FVAL(LMg2SiO4) = FVAL(LMg2SiO4) - CONDFo
      TL(LMg2SiO4) = TL(LMg2SiO4) + CONDFo*DEN
      TP(LMg2SiO4) = TP(LMg2SiO4)
C
      SiO2COND = CONDSiO2 * DEN 
      EnCOND = CONDEn * DEN
      FoCOND = CONDFo * DEN
C
C
C ***** SPECIES MASS CONSERVATION *****
      DO 6 I = 1,NQ
      CON(I) = TP(I) - TL(I)
   6  CONTINUE
C
C
C ***** FRACTIONAL ELEMENTAL MASS CONSERVATION *****
C   OXYGEN CONSERVATION
      CONO = (AMADD + SADD + OADD - SiO2COND*2. - EnCOND*3.
     2  - FoCOND*4. - H2OLOSS) / (SADD + AMADD + OADD)
C
C
C   SILICON CONSERVATION
      CONSi = (SADD - SiO2COND - EnCOND - FoCOND) / SADD
C
C   MAGNESIUM CONSERVATION
      CONMg = (AMADD - EnCOND - FoCOND*2.) / AMADD
C  
C
C ***** CALCULATE PRODUCTION (AND LOSS) RATES *****
      DO 10 L=1,NR
  10  RAT(L) = 0.
C
      DO 12 L=1,NR
      M = JCHEM(1,L)
      K = JCHEM(2,L)
  12  RAT(L) = RAT(L) + A(L)*D(M)*D(K)
C    
C
      RETURN
      END
C-PK *******************************
      SUBROUTINE CHEMPL(D,XP,XL,K)
      INCLUDE 'header.inc'
      DIMENSION D(NSP)
      COMMON/RBLOK/A(NR),ILOSS(2,NSP,NMAX),IPROD(NSP,NMAX),
     2  JCHEM(5,NR),NUML(NSP),NUMP(NSP),ISPEC(NSP),T,DEN
C
C   THIS SUBROUTINE CALCULATES CHEMICAL PRODUCTION AND LOSS RATES
C   USING THE INFORMATION IN THE MATRICES JCHEM, ILOSS, AND IPROD.
C   CALLED BY SUBROUTINE DOCHEM.
C
C   INITIALIZE VALUES AS ZEROS
      XL = 0. !LOSS FREQUENCY  [S^-1]
      XP = 0. !PRODUCTION RATE [CM^-3 S^-1]
C
C ** LOSS FREQUENCY XL (1/S)
      NL = NUML(K)
      DO 2 L=1,NL
      J = ILOSS(1,K,L)
      M = ILOSS(2,K,L)
   2  XL = XL + A(J)*D(M)
C     PRINT *, A(J),D(M)
C
C ** PRODUCTION RATE XP (MOL/CM3/S)
      NP = NUMP(K)
      DO 3 L=1,NP
      J = IPROD(K,L)
      M = JCHEM(1,J)
      N = JCHEM(2,J)
   3  XP = XP + A(J)*D(M)*D(N)
C
C ** DIAGNOSTIC PRINT STATEMENTS **
C     PRODO = A(J)*D(M)*D(N)
C     PRINT *,'L =',L,' J =',J,' PRODO =',PRODO
C     PRINT *,'M =',M,' N =',N
C     PRINT *,'A(J) =',A(J),' D(M) =',D(M),'D(N) =',D(N)C
C     STOP
C
      RETURN
      END
C
C-PK ****************************************
      SUBROUTINE LUDCMP(A,N,NP,INDX,D)
      PARAMETER (NMAX=12,TINY=1.0E-20)
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
