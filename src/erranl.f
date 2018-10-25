      SUBROUTINE ERRANL(PHI,U,V,DIV,ZETA,LN,GPHS,L2CTR,NFLAG)
C
C THIS ROUTINE COMPUTES VARIOUS FORMS OF SCALAR ERROR ESTIMATES          
C BY COMPARISON WITH THE ANALYTIC SOLUTIONS
C THE ROUTINE THEN PLOTS THE ANALYTIC SOLUTION AND THE DEPARTURE            
C (THE ERROR) BY CALLING THE STANDARD PLOTTING ROUTINE WITH 
C VARIABLES PLACED INTO THE WORKSPACE (IF GPHS = .TRUE.).        
C                                                                              
C CALLED BY: STSWM
C CALLS: AGSETI, ANLYTC, ANOTAT, EZXY, GLAT, PLOTS, WEIGHT
C
C REVISIONS:
C 7-10-92 CHANGE TO CCM CODING CONVENTIONS (R. JAKOB)
C
C---- Model Parameters -------------------------------------------------
C
      INCLUDE 'params.i'
C
C---- Common Blocks ----------------------------------------------------
C
C Constants
      INCLUDE 'consts.i'
C Initial Conditions (test case 4 special treatment)
      INCLUDE 'finit.i'
C
C------------ Arguments ------------------------------------------------
C
C     Input
C
C HEIGHT
      REAL PHI(NLON+2,NLAT)
C U WIND VELOCITY
      REAL U(NLON+2,NLAT)
C V WIND VELOCITY
      REAL V(NLON+2,NLAT)
C DIVERGENCE
      REAL DIV(NLON+2,NLAT,LVLS)
C VORTICITY
      REAL ZETA(NLON+2,NLAT,LVLS)
C CURRENT TIMELEVEL OF DIVERGENCE,VORTICITY
      INTEGER LN
C IF .TRUE. PLOT CURRENT ANALYTIC SOLUTION AND ERROR
      LOGICAL GPHS
C ARRAY LOCATION TO STORE ERRORS FOR LATER PLOTTING
      INTEGER L2CTR
C IF .TRUE. CAUSES PLOTS OF TIME HISTORY OF ERRORS
      LOGICAL NFLAG
C
C------ Local Variables ------------------------------------------------
C
C     STORAGE ARRAYS FOR ERROR MEASURES
C     SUFFICIENT FOR UP TO NGRPHS ENTRIES                
C     SAVE ERROR MEASURES IN ARRAY FOR PLOTTING
C                                                                               
C     L1, L2 AND L(INFINITY) ERRORS FOR HEIGHT AND
C     VECTOR VELOCITY, TIME
C
      REAL PL1G(NGRPHS),VL1G(NGRPHS),VL2G(NGRPHS), PL2G(NGRPHS),
     $     PINFG(NGRPHS),VINFG(NGRPHS),TME(NGRPHS)
      SAVE PL1G,VL1G,PL2G,VL2G,PINFG,VINFG,TME
C
C     MINIMUM, MAXIMUM, MEAN AND VARIANCE
C
      REAL PMING(NGRPHS),PMAXG(NGRPHS),PAVGG(NGRPHS),PVARG(NGRPHS)
      SAVE PMING,PMAXG,PAVGG,PVARG
C
C     GRID WEIGHTS
C
      REAL WTS
C
C     NORMALIZATION CONSTANTS (FOR COMPARABLE GRAPHS)
C     SAVES INITIAL VALUES OF ANALYTIC SOLUTION
C     AND THEN PLOTS ERRORS AS RELATIVE ERRORS
C
      REAL MINN,MAXN,MEANN,VARN
      SAVE MINN,MAXN,MEANN,VARN
C
C     FIELDS FOR ANALYTIC SOLUTION 
C
      REAL UANL(NLON+2,NLAT),VANL(NLON+2,NLAT),PANL(NLON+2,NLAT),
     $     DANL(NLON+2,NLAT),ZANL(NLON+2,NLAT)
C
C     TEMPORARY VARIABLES FOR COMPUTED SOLUTION
C
      REAL MINP,MAXP,MAXV,MEAN,VAR,
     $     VMAX,PINF,VL1,PL1,VL2,PL2,VLEN
C
C     TEMPORARY VARIABLES FOR ANALYTIC SOLUTION
C
      REAL MEANT,VART,MINT,MAXT,
     $     VMAXT,PINFT,VL1T,PL1T,VL2T,PL2T,VLENT,PABS
C
C     TEMPORARY VARIABLES FOR ERROR MEASURES
C
      INTEGER I,J
      REAL VDIFF
C
C     CHARACTER ARRAYS FOR PLOT TITLES
C
      CHARACTER*64 GPH
      CHARACTER*40 XAXIS
      CHARACTER*1 DASH(1)
      SAVE DASH
C
C----- External Functions ----------------------------------------------
C
C     GAUSSIAN GRID AND WEIGHTS
C
      EXTERNAL GLAT, WEIGHT
      REAL GLAT, WEIGHT
C
C----- Initialized Variables -------------------------------------------
C
C     FULL LINE GRAPHS
C
      DATA DASH(1) / '$' /
C                                                                               
C----- Executable Statements -------------------------------------------        
C
      IF (NFLAG) THEN                                                    
         IF (.NOT. LGPHS) THEN
C
C           NO PLOTING OUTPUT -> RETURN
C
         ELSE
C                                                                               
C        PLOT TIME DEPENDENCE OF ERRORS                             
C                                                                               
C        GRAPH AXIS TEXT 
C
         CALL AGSETI('LINE/MAXIMUM.', 64)
         XAXIS = 'TIME IN HOURS$'
C
C        PLOT HEIGHT ERROR OVER TIME
C
         WRITE(GPH,5) 'HEIGHT ERROR,TEST ',ICOND,CHEXP,STRUNC,INT(DT)
    5    FORMAT(A18,I2,',EXP.',A4,',',A6,',DT = ',I4,' SEC$')
C
         CALL ANOTAT(XAXIS,'L1-NORM$',1,1,0,DASH)
         CALL EZXY(TME, PL1G, L2CTR, GPH)
C
         CALL ANOTAT(XAXIS,'L2-NORM$',1,1,0,DASH)
         CALL EZXY(TME, PL2G, L2CTR, GPH)
C
         CALL ANOTAT(XAXIS,'L(INF)-NORM$',1,1,0,DASH)
         CALL EZXY(TME, PINFG, L2CTR, GPH)
C
C        PLOT WIND ERROR OVER TIME
C
         WRITE(GPH,5) 'WIND ERROR,  TEST ',ICOND,CHEXP,STRUNC,INT(DT)

         CALL ANOTAT(XAXIS,'L1-NORM$',1,1,0,DASH)
         CALL EZXY(TME,VL1G, L2CTR, GPH)
C
         CALL ANOTAT(XAXIS,'L2-NORM$',1,1,0,DASH)
         CALL EZXY(TME,VL2G, L2CTR, GPH)
C
         CALL ANOTAT(XAXIS,'L(INF)-NORM$',1,1,0,DASH)
         CALL EZXY(TME, VINFG, L2CTR, GPH)
C
C        PLOT MINIMUM/MAXIMUM
C
         WRITE(GPH,5) 'HEIGHT MIN., TEST ',ICOND,CHEXP,STRUNC,INT(DT)
C
         CALL ANOTAT(XAXIS,'RELATIVE ERROR$',1,1,0,DASH)
         CALL EZXY(TME, PMING, L2CTR, GPH)
C
         WRITE(GPH,5) 'HEIGHT MAX., TEST ',ICOND,CHEXP,STRUNC,INT(DT)
C
         CALL ANOTAT(XAXIS,'RELATIVE ERROR$',1,1,0,DASH)
         CALL EZXY(TME, PMAXG, L2CTR, GPH)
C
C        PLOT NORMALIZED MEAN/VARIANCE
C
         WRITE(GPH,5) 'HEIGHT AVG., TEST ',ICOND,CHEXP,STRUNC,INT(DT)
C
         CALL ANOTAT(XAXIS,'RELATIVE ERROR$',1,1,0,DASH)
         CALL EZXY(TME, PAVGG, L2CTR, GPH)
C
         WRITE(GPH,5) 'HEIGHT VAR., TEST ',ICOND,CHEXP,STRUNC,INT(DT)
C
         CALL ANOTAT(XAXIS,'RELATIVE ERROR$',1,1,0,DASH)
         CALL EZXY(TME, PVARG, L2CTR, GPH)
C
C        DONE TIME SERIES PLOTS -> RETURN
C
         ENDIF
C                                                                               
      ELSE
C
C     DO ERROR ANALYSIS ON MODEL GRID
C                                                                               
C     INITIALIZE ERROR VARIABLES
C                                                                               
      VMAX  = 0.0                                                               
      VMAXT = 0.0
      PINF  = 0.0                                                               
      PINFT = 0.0
      VL1   = 0.0                                                               
      VL1T  = 0.0
      PL1   = 0.0                                                               
      PL1T  = 0.0
      VL2   = 0.0                                                               
      VL2T  = 0.0
      PL2   = 0.0                                                               
      PL2T  = 0.0
      MINP  = PHI(1,1)
      MAXP  = PHI(1,1)
      MAXV  = 0.0
      MEAN  = 0.0
      MEANT = 0.0
      VAR   = 0.0
      VART  = 0.0
C
C     GET ANALYTIC SOLUTION
C     CALL ROUTINE FOR TIME DEPENDENT ANALYTIC SOLUTIONS WHICH WILL BE
C     STUCK IN UANL, VANL, PANL, DANL AND ZANL
C
      CALL ANLYTC(1,1,NSTEP*DT,ICOND,UANL,VANL,PANL,DANL,ZANL)
C                                                                               
      DO 120 J=1,NLAT                                                           
C                                                                               
C        WEIGHT INDEPENDENT OF LONGITUDE
C
         WTS = WEIGHT(1,J)
C
C        CALCULATE COMPONENTS OF L2 NORM, AND OTHER MEASURES OF ERROR
C        FOR SMALL U, SMALL V, AND HEIGHT
C
         DO 10 I=1,NLON                                                         
C
C           WEIGH DATA DEPENDENT ON LATITUDE & LONGITUDE
C
C           WTS   = WEIGHT(I,J)
C
            VLEN  = SQRT(U(I,J)**2+V(I,J)**2)
            VDIFF = SQRT((U(I,J)-UANL(I,J))**2 +
     $              (V(I,J)-VANL(I,J))**2)
C
C           L1 ERRORS
C
            VL1   = VL1 + VDIFF*WTS
            PL1   = PL1 + ABS(PHI(I,J)-PANL(I,J))*WTS
C
C           L2 ERRORS
C
            VL2   = VL2 + VDIFF**2*WTS
            PL2   = PL2 + (PHI(I,J)-PANL(I,J))**2*WTS
C
C           L(INFINITY) ERRORS
C
            IF (VDIFF .GT. VMAX) THEN
               VMAX  = VDIFF
            ENDIF
            IF (ABS(PANL(I,J)-PHI(I,J)) .GT. PINF) THEN
               PINF  = ABS(PANL(I,J)-PHI(I,J))
            ENDIF
C
C           MAXIMUM WIND COMPUTED
C
            IF (PHI(I,J) .LT. MINP) THEN
               MINP = PHI(I,J)
            ENDIF
C
            IF (PHI(I,J) .GT. MAXP) THEN
               MAXP = PHI(I,J)
            ENDIF
            IF (VLEN .GT. MAXV) THEN
               MAXV = VLEN
            ENDIF
C
C           MINIMUM/MAXIMUM ANALYTIC SOLUTION
C
            IF ((I .EQ. 1) .AND. (J .EQ. 1)) THEN
               MINT = PANL(1,1)
               MAXT = PANL(1,1)
            ELSE
               IF (PANL(I,J) .LT. MINT) THEN
                  MINT = PANL(I,J)
               ENDIF
               IF (PANL(I,J) .GT. MAXT) THEN
                  MAXT = PANL(I,J)
               ENDIF
            ENDIF
C
C           MEAN VALUES
C
            MEAN = MEAN + PHI(I,J)*WTS
            MEANT = MEANT + PANL(I,J)*WTS
C
C           ERROR NORMALIZATION
C
            IF (ICOND .NE. 4) THEN
               PABS  = ABS(PANL(I,J))
               VLENT = SQRT(UANL(I,J)**2+VANL(I,J)**2)
            ELSE
C
C     SPECIAL HANDLING FOR CASE 4: SUBTRACT MEAN FLOW TO COMPUTE
C     ERROR RELATIVE TO FORCING
C     SEE PAPER BY BROWNING ET. AL.: A COMPARISON OF THREE
C     NUMERICAL METHODS FOR ... P. 1072 !
C       
               PABS  = ABS(PANL(I,J) - PHICON(J))
               VLENT = SQRT((UANL(I,J) - UCON(J))**2 
     $                 + (VANL(I,J) - VCON(J))**2)
            ENDIF
            PL1T  = PL1T + PABS*WTS
            VL1T  = VL1T + VLENT*WTS
            PL2T  = PL2T + PABS**2*WTS
            VL2T  = VL2T + VLENT**2*WTS
C
            IF (VLENT .GT. VMAXT) THEN 
               VMAXT  = VLENT
            ENDIF                                                               
            IF (PABS .GT. PINFT) THEN 
               PINFT  = PABS 
            ENDIF                                                               
   10    CONTINUE                                                               
C                                                                               
C        NEXT LATITUDE
C                                                                               
  120 CONTINUE                                                                  
C
C     COMPUTE VARIANCE
C
      DO 70 J=1,NLAT
         WTS = WEIGHT(1,J)
         DO 60 I=1,NLON
C           WTS = WEIGHT(I,J)
            VAR  = VAR + (PHI(I,J)-MEAN)**2*WTS
            VART = VART + (PANL(I,J)-MEANT)**2*WTS
   60    CONTINUE
   70 CONTINUE
C
C     CALCULATE L1 ERRORS
C
      IF (VL1T .NE. 0.0) THEN
         VL1G(L2CTR) = VL1/VL1T
      ELSE
         VL1G(L2CTR) = VL1
      ENDIF
C
      IF (PL1T .NE. 0.0) THEN
         PL1G(L2CTR) = PL1/PL1T
      ELSE
         PL1G(L2CTR) = PL1
      ENDIF
C                                                                               
C     CALCULATE L2 ERRORS 
C                                                                               
      IF (VL2T .NE. 0.0) THEN                                                  
         VL2G(L2CTR) = SQRT(VL2/VL2T) 
      ELSE                                                                      
         VL2G(L2CTR) = SQRT(VL2)
      ENDIF                                                                     
C                                                                               
      IF (PL2T .NE. 0.0) THEN
         PL2G(L2CTR) = SQRT(PL2/PL2T)
      ELSE                                                                      
         PL2G(L2CTR) = SQRT(PL2)
      ENDIF                                                                     
C
C     CALCULATE L(INFINITY) ERRORS
C
      IF (VMAXT .NE. 0.0) THEN
         VINFG(L2CTR) = VMAX/VMAXT
      ELSE
         VINFG(L2CTR) = VMAX
      ENDIF
C
      IF (PINFT .NE. 0.0) THEN
         PINFG(L2CTR) = PINF/PINFT
      ELSE
         PINFG(L2CTR) = PINF
      ENDIF
C
C     HEIGHT MIN,MAX,MEAN AND VARIANCE
C     SAVE INITIAL VALUES FOR TIME-INDEPENDENT NORMALIZATION
C
      IF (L2CTR .EQ. 1) THEN
         MINN  = MINT
         MAXN  = MAXT
         MEANN = MEANT
         VARN  = VART
         WRITE(6,986) MINN,MAXN,MEANT,VART
  986    FORMAT (/, ' ERRANL: INITIAL VALUES FOR NORMALIZ',
     $              'ATION OF RELATIVE ERRORS',/,
     $              ' (SLIGHTLY GRID DEPENDENT!)',/,
     $              ' HEIGHT MIN./MAX. = ',1PE16.9,'/',1PE16.9,/,
     $              ' HEIGHT AVG./VAR. = ',1PE16.9,'/',1PE16.9)
C
      ELSEIF (L2CTR .GT. 1) THEN
C                                                                               
C     PRINT RESULTS
C
         WRITE (6,987) NSTEP, TAU
  987    FORMAT (/, ' ERRANL: ERROR ESTIMATES FOR NSTEP = ', I4, 
     $           ', TAU = ', 0PF6.2, ' HRS') 
C
         WRITE (6,122) PL1G(L2CTR),PL2G(L2CTR),PINFG(L2CTR)
  122    FORMAT (' HEIGHT ERROR',/,
     $           ' L1 = ', 1PE16.9, ', L2 = ', 1PE16.9,
     $           ' L(INF) = ', 1PE16.9)
C                                                                               
         WRITE (6,123) VL1G(L2CTR),VL2G(L2CTR),VINFG(L2CTR)
  123    FORMAT (' VECTOR WIND ERROR',/,
     $           ' L1 = ', 1PE16.9, ', L2 = ', 1PE16.9,
     $           ' L(INF) = ', 1PE16.9)
C
         WRITE (6,124) MINP,MAXP
  124    FORMAT (' HEIGHT MIN./MAX.  = ',1PE16.9,'/',1PE16.9)
C
         WRITE (6,126) MEAN,VAR
  126    FORMAT (' HEIGHT AVG./VAR.  = ',1PE16.9,'/',1PE16.9)
C
C        CHECK FOR ARRAY OVERFLOW
C                                                                               
         IF (L2CTR .GE. NGRPHS) THEN 
            WRITE (0,905) L2CTR 
  905       FORMAT(/,' STSWM: FATAL ERROR IN SUBROUTINE ERRANL ',  
     $         /,' ARGUMENT L2CTR EXCEEDS ALLOCATED ',          
     $           'STORAGE FOR VARIABLES L2CTR = ',I4, /)           
            STOP
         ENDIF
      ENDIF
C
C     SAVE RELATIVE ERRORS FOR LATER PLOTTING
C
      PMING(L2CTR) = (MINP-MINT)/(MAXN-MINN)
      PMAXG(L2CTR) = (MAXP-MAXT)/(MAXN-MINN)
      PAVGG(L2CTR) = (MEAN-MEANT)/MEANN
      PVARG(L2CTR) = (VAR-VART)/VARN
C
C     SAVE TIME POINT
C
      TME(L2CTR) = TAU
C                                                                               
C     2D PLOTS OF ANALYTIC SOLUTION AND ERROR
C                                                                               
      IF (GPHS) THEN
C                                                                               
C        PLOTS OF ANALYTIC SOLUTION
C        CALLING ARGUMENT = 2
C                                                                               
         CALL PLOTS (2,PANL,PHI,UANL,VANL,ZANL,DANL,1,1)
C                                                                               
C        NEXT CALL TO GRAPHICS PROCEDURE PLOTS ERRORS IN SOLUTION 
C        STICK ERROR FIELDS INTO UANL, VANL, ZANL, DANL, AND PANL 
C                                                                               
         DO 260 J=1,NLAT
C                                                                               
C           ANALYTIC SOLUTION IS STILL IN xANL FROM PREVIOUS CALL TO
C           ANLYTC; SUBTRACT FROM COMPUTED SOLUTION FOR ERROR
C                                                                               
            DO 250 I=1,NLON 
               UANL(I,J) = U(I,J)-UANL(I,J)
               VANL(I,J) = V(I,J)-VANL(I,J)
               PANL(I,J) = PHI(I,J)-PANL(I,J)
               ZANL(I,J) = ZETA(I,J,LN)-ZANL(I,J)  
               DANL(I,J) = DIV(I,J,LN)-DANL(I,J)
  250       CONTINUE 
  260    CONTINUE   
C                                                                               
C        PLOTS OF ERROR (ANALYTIC - COMPUTED)
C        CALLING ARGUMENT = 1 
C                                                                               
         CALL PLOTS (1,PANL,PANL,UANL,VANL,ZANL,DANL,1,1)
C                                                                               
C        END 2D PLOTS
C
      ENDIF
C                                                                               
C     END ERROR ANALYSIS ON GRID
C                                                                               
      ENDIF
C
      RETURN                                                                    
      END                                                                       
