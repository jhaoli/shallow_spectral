      SUBROUTINE SPCANL(PHISC,ZETASC,DIVSC) 
C                                                                              
C THIS ROUTINE COMPUTES THE SPECTRUM FOR VARIOUS QUANTITIES
C AND PLOTS THEM. CURRENTLY THE CODE DOES SPECTRAL ANALYSIS FOR
C THE HEIGHT AND  KINETIC ENERGY.
C                                                                              
C CALLED BY: STSWM
C CALLS: AGSETI, AGSETR, ANOTAT, EZXY
C
C---- Model Parameters -------------------------------------------------
C
      INCLUDE 'params.i'
C
C---- Common Blocks ----------------------------------------------------
C
C CONSTANTS
      INCLUDE 'consts.i'
C TRANSFORM ARRAYS
      INCLUDE 'trnsfm.i'
C
C------------ Arguments ------------------------------------------------
C
C     Input
C
C COMPUTED HEIGHT SPECTRAL COEFFICIENTS
      COMPLEX PHISC(NALP)
C COMPUTED VORTICITY SPECTRAL COEFFICIENTS
      COMPLEX ZETASC(NALP)
C COMPUTED DIVERGENCE SPECTRAL COEFFICIENTS
      COMPLEX DIVSC(NALP)
C
C------ Local Parameters -----------------------------------------------
C
C     SPECIAL VALUE FOR SPECTRUM PLOTS
C
      REAL SPVAL
      PARAMETER (SPVAL=1.0E36)
C
C------ Local Variables ------------------------------------------------
C                                                                               
C     STORAGE ARRAYS FOR SPECTRUM
C                                                                               
      REAL PHIG(KK),KEG(KK),WAVE(KK)
C
C     TEMPORARY VARIABLES FOR SPECTRAL ANALYSIS
C
      INTEGER M,N,K,IS,JN,JM
      REAL FAC
      COMPLEX ETASC
C
C     CHARACTER ARRAYS FOR PLOT TITLES
C
      CHARACTER*64 GPH
      CHARACTER*40 XAXIS
      CHARACTER*1 DASH(1)
C
C----- Initialized Variables -------------------------------------------
C
C     FULL LINE GRAPHS
C
      DATA DASH(1) / '$' /
C                                                                               
C----- Executable Statements -------------------------------------------
C
C     INITIALIZE ARRAYS
C
      DO 100 K=1,KK
         WAVE(K) = REAL(K)
         PHIG(K) = 0.0
         KEG(K) = 0.0
  100 CONTINUE
C
C     ADD MERIDIONAL WAVES ALONG DIAGONALS
C
      DO 120 JN=0,NN                                                           
        DO 110 JM=1,LDIAG(JN,1)
C
C       USE SYMMETRY OF M
C
        IF (JM .EQ. 1) THEN
           FAC = 1.0 
        ELSE
           FAC = 2.0 
        ENDIF
C
C          COMPUTE SPECTRAL INDEX
C
           M = JM - 1
           N = JM - 1 + JN
           IS = LDIAG(JN,2)
C
C          COMPUTE REALTIVE VORTICITY BY SUBTRACTING
C          SPECTRAL CORIOLIS TERM
C
           IF ((N .EQ. 1) .AND. (M .EQ. 0)) THEN
              ETASC = ZETASC(IS+JM) - CORSC1
           ELSEIF ((N .EQ. 1) .AND. (M .EQ. 1)) THEN
              ETASC = ZETASC(IS+JM) - CORSC2
           ELSE
              ETASC = ZETASC(IS+JM)
           ENDIF
C
C          IGNORE GLOBAL MEAN !
C
           IF (N .NE. 0) THEN
              KEG(N) = KEG(N) + FAC *
     $                 (DIVSC(IS+JM)*CONJG(DIVSC(IS+JM)) + 
     $                  ETASC*CONJG(ETASC))/(4.0*A2NNP1(N))
              PHIG(N) = PHIG(N) + FAC *
     $                  PHISC(IS+JM)*CONJG(PHISC(IS+JM))*0.5
           ENDIF
  110   CONTINUE
C
  120 CONTINUE                                                                  
C
C     PRINT RESULTS
C
      WRITE (6,987) NSTEP, TAU
  987 FORMAT (/, ' SPCANL: SPECTRAL ANALYSIS FOR NSTEP = ', I4,
     $           ', TAU = ', 0PF6.2, ' HRS')
C
      DO 200 K=1,KK
         WRITE (6,124) K,PHIG(K),KEG(K)
  124    FORMAT (' WAVE K = ',I4,' PHI/KE AMPLITUDE = ',1PE16.9,
     $           '/',1PE16.9)
C
C        CHECK FOR ZERO VALUES (CAUSES PROBLEMS WITH LOGARITHMIC PLOTS)
C
         IF (PHIG(K) .EQ. 0.0) THEN
            PHIG(K) = SPVAL       
         ENDIF
         IF (KEG(K) .EQ. 0.0) THEN
            KEG(K) = SPVAL
         ENDIF
  200 CONTINUE
C
C     LOGARITHMIC SPECTRUM PLOTS AS FUNCTION OF WAVENUMBER
C
      IF (LGPHS) THEN
C
    5    FORMAT(A21,I2,',EXP.',A4,',',A6,',TAU = ',0PF6.2,' HRS$')
C
C        GRAPH AXIS TEXT
C
         CALL AGSETI('LINE/MAXIMUM.', 64)
         CALL AGSETR('X/LOGARITHMIC.', -1.)
         CALL AGSETR('Y/LOGARITHMIC.', -1.)
         XAXIS = 'WAVENUMBER N (LOGARITHMIC)$'
C
C        PLOT HEIGHT ERROR OVER TIME
C
         WRITE(GPH,5) 'PHI SPECTRUM,   TEST ',ICOND,CHEXP,
     $                STRUNC, TAU
C
         CALL ANOTAT(XAXIS,'M^2/S^2 (LOGARITHMIC)$',1,1,0,DASH)
         CALL EZXY(WAVE, PHIG, KK, GPH)
C
         WRITE(GPH,5) 'KE  SPECTRUM,   TEST ',ICOND,CHEXP,
     $                STRUNC, TAU
C
         CALL ANOTAT(XAXIS,'J/KG (LOGARITHMIC)$',1,1,0,DASH)
         CALL EZXY(WAVE, KEG, KK, GPH)
C
         CALL AGSETR('X/LOGARITHMIC.', 0.)
         CALL AGSETR('Y/LOGARITHMIC.', 0.)
C
C     GRAPHICS DONE
C                                                                               
      ENDIF
C
      RETURN                                                                    
      END                                                                       
