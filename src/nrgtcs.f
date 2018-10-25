      SUBROUTINE NRGTCS(H, U, V, ZETA, DIV, MOUNT, LN, ICTR, NFLAG)
C                                                                              
C THIS ROUTINE COMPUTES A NUMBER OF ENERGETICS PARAMETERS WHEN              
C INVOKED, STORES THEM IN LOCAL ARRAYS SO THEY MAY BE PLOTTED AS            
C A FUNCTION OF TIME, AND PLOTS THEM WHEN INVOKED WITH NFLAG=.TRUE.
C                                                                              
C CALLED BY: STSWM
C CALLS: AGSETI, ANOTAT, EZXY, GLAT, GLON, WEIGHT
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
C
C------------ Arguments ------------------------------------------------
C
C     Input
C
C HEIGHT FIELD
      REAL H(NLON+2,NLAT)
C EASTWARD WIND
      REAL U(NLON+2,NLAT)
C NORTHWARD WIND
      REAL V(NLON+2,NLAT)
C VORTICITY FIELD
      REAL ZETA(NLON+2,NLAT,LVLS)
C DIVERGENCE FIELD
      REAL DIV(NLON+2,NLAT,LVLS)
C MOUNTAIN HEIGHT
      REAL MOUNT(NLON+2,NLAT)
C CURRENT TIMELEVEL OF VORTICITY/DIVERGENCE
      INTEGER LN
C ARRAY INDEX FOR ENERGETICS DATA STORE
      INTEGER ICTR
C IF .TRUE., PLOT TIME EVOLUTION OF DATA
      LOGICAL NFLAG
C
C------ Local Parameters -----------------------------------------------
C
      REAL PI
      PARAMETER (PI=3.141592653589793)
C
C------ Local Variables ------------------------------------------------
C
C     STORAGE FOR ENERGETICS PARAMETERS FOR UP TO NGRPHS ENTRIES                
C     SAVE ANALYSIS DATA FOR SUBSEQUENT PLOTTING
C                                                                               
      REAL TE(NGRPHS), PE(NGRPHS), APE(NGRPHS),              
     $          AAM(NGRPHS), TMSS(NGRPHS), TME(NGRPHS),
     $          HTS(NGRPHS)
      SAVE TE, PE, APE, AAM, TMSS, TME, HTS
C
C     SAVE INITIAL VALUES FOR NORMALIZED  GRAPHS
C
      REAL TOTEN,PEN,TMSSN
      SAVE TOTEN,PEN,TMSSN
C                                                                               
C     TEMPORARY ARRAYS
C
      REAL TMP1(NLON), TMP2(NLON), TMP3(NLON), TMP4(NLON),
     $          TMP5(NLON)
C
C     TEMPORARY SUMMATION VARIABLES
C
      REAL TPE, TAM, TAPE, TMASS, TOTE
C
C     TEMPORARIES
C
      INTEGER I,J
      REAL HFLUID, WTS 
C
C     CHARACTER ARRAYS FOR PLOT TITLES
C
      CHARACTER*64 GPH
      CHARACTER*40 XAXIS
      CHARACTER*1 DASH(1)
      SAVE DASH
C
C     HEIGHT TIMESERIES FOR POINT (LONGITUDE,LATITUDE)
C     (IN DEGREES EAST OF GREENWICH, NORTH OF EQUATOR)
C
      INTEGER ILON,JLAT
      SAVE ILON,JLAT
      REAL PLON,PLAT,PLONG,PLATG
      SAVE PLON,PLAT
C
C----- External Functions ----------------------------------------------
C
      EXTERNAL GLAT, GLON
      REAL GLAT,GLON
C
C     RELATIVE WEIGHTS
C
      EXTERNAL WEIGHT
      REAL WEIGHT
C
C----- Initialized Variables -------------------------------------------
C
C     HEIGHT TIMESERIES FOR POINT (LONGITUDE,LATITUDE)
C     (NEAR BOULDER, COLORADO)
C
      DATA PLON,PLAT / 255.0, 40.0 /
C
C     GRAPH LINE TYPE (=SOLID)
C
      DATA DASH(1) / '$' /
C
C----- Executable Statements -------------------------------------------        
C                                                                               
      IF (NFLAG) THEN  
         IF (.NOT. LGPHS) THEN
C
C           NO TIME SERIES PLOTS -> RETURN
C
         ELSE
C                                                                               
C        PLOT TIME DEPENDENCE OF ENERGETICS PARAMETERS 
C                                                                               
    5    FORMAT(A21,I2,',EXP.',A4,',',A6,',DT = ',I4,' SEC$') 
C                                                                               
C        TITLE LENGTH
C
         CALL AGSETI('LINE/MAXIMUM.', 64)
         XAXIS = 'TIME IN HOURS$'
C
C        PLOTS OVER TIME
C
         WRITE(GPH,5) 'MASS,           TEST ',ICOND,CHEXP,
     $                STRUNC, INT(DT)
         CALL ANOTAT(XAXIS,'RELATIVE ERROR$',1,1,0,DASH)
         CALL EZXY(TME,TMSS, ICTR, GPH)
C
         WRITE(GPH,5) 'TOTAL ENERGY,   TEST ',ICOND,CHEXP,
     $                STRUNC, INT(DT)
         CALL ANOTAT(XAXIS,'RELATIVE ERROR$',1,1,0,DASH)
         CALL EZXY(TME,TE, ICTR, GPH)                                       
C
         WRITE(GPH,5) 'VORTICITY,      TEST ',ICOND,CHEXP,
     $                STRUNC, INT(DT)
         CALL ANOTAT(XAXIS,'1/SECOND$',1,1,0,DASH)
         CALL EZXY(TME, APE, ICTR, GPH)                                       
C
         WRITE(GPH,5) 'DIVERGENCE,     TEST ',ICOND,CHEXP,
     $                STRUNC, INT(DT)
         CALL ANOTAT(XAXIS,'1/SECOND$',1,1,0,DASH)
         CALL EZXY(TME, AAM, ICTR, GPH)                                       
C
         WRITE(GPH,5) 'POT. ENSTROPHY, TEST ',ICOND,CHEXP,
     $                STRUNC, INT(DT)
         CALL ANOTAT(XAXIS,'RELATIVE ERROR$',1,1,0,DASH)
         CALL EZXY(TME,  PE, ICTR, GPH)
C
         WRITE(GPH,5) 'POINT HEIGHT,   TEST ',ICOND,CHEXP,
     $                STRUNC, INT(DT)
         CALL ANOTAT(XAXIS,'METER$',1,1,0,DASH)
         CALL EZXY(TME,HTS, ICTR, GPH)
C                                                                               
C        END TIME SERIES PLOTS -> RETURN
C                                                                               
         ENDIF
C
      ELSE
C
C-----ENERGETICS ANALYSIS ON MODEL GRID-------------------------------
C                                                                               
C     INITIALIZE LOCAL SUMMATION VARIABLES FOR ENERGY TOTALS
C                                                                               
      TOTE  = 0.0
      TPE   = 0.0                                                               
      TAM   = 0.0                                                               
      TAPE  = 0.0                                                               
      TMASS = 0.0                                                               
C
      DO 100 I=1,NLON                                                           
         TMP1(I) = 0.0                                                          
         TMP2(I) = 0.0                                                          
         TMP3(I) = 0.0                                                          
         TMP4(I) = 0.0
         TMP5(I) = 0.0
  100 CONTINUE                                                                  
C                                                                               
C     AVERAGE KINETIC ENERGY, POTENTIAL ENERGY, AND PSEUDO-MASS                 
C                                                                               
      DO 120 J=1,NLAT                                                           
C
C        WEIGHT INDEPENDENT OF LONGITUDE
C
         WTS = WEIGHT(1,J)
         DO 110 I=1,NLON                                                        
C
C           WEIGHT DEPENDENT ON LATITUDE & LONGITUDE
C
C           WTS = WEIGHT(I,J)
            IF (FTOPO) THEN
               HFLUID = H(I,J) - MOUNT(I,J)
               TMP1(I) = TMP1(I) + 0.5*ZETA(I,J,LN)**2/HFLUID*WTS
               TMP2(I) = TMP2(I) + HFLUID*WTS
               TMP5(I) = TMP5(I) + ((U(I,J)**2+V(I,J)**2)*HFLUID 
     $              + (H(I,J)**2 - MOUNT(I,J)**2)*GRAV)/2.0*WTS
            ELSE
               HFLUID = H(I,J) 
               TMP1(I) = TMP1(I) + 0.5*ZETA(I,J,LN)**2/HFLUID*WTS
               TMP2(I) = TMP2(I) + HFLUID*WTS
               TMP5(I) = TMP5(I) + (U(I,J)**2+V(I,J)**2 +
     $              GRAV*HFLUID)*HFLUID/2.0*WTS
            ENDIF
            TMP3(I) = TMP3(I) + ZETA(I,J,LN)*WTS
            TMP4(I) = TMP4(I) + DIV(I,J,LN)*WTS
  110    CONTINUE                                                               
  120 CONTINUE                                                                  
C                                                                               
      DO 130 I=1,NLON                                                           
         TPE   = TPE   + TMP1(I)
         TMASS = TMASS + TMP2(I)                                                
         TAPE  = TAPE  + TMP3(I)
         TAM   = TAM   + TMP4(I)
         TOTE  = TOTE  + TMP5(I)
  130 CONTINUE                                                                  
C                                                                               
C     WRITE OUT ENERGETICS INFORMATION                                 
C                                                                               
      WRITE(6,199) NSTEP, TAU
  199 FORMAT (/, ' CONSERVATION ANALYSIS FOLLOWS FOR NSTEP = ', I4,
     $           ', TAU = ', 0PF6.2, ' HRS')
C
C     CHECK FOR DATA ARRAY OVERFLOW
C                                                                               
      IF (ICTR .GE. NGRPHS) THEN  
         WRITE (0,905) NGRPHS  
  905    FORMAT(/,' STSWM: FATAL ERROR IN SUBROUTINE NRGTCS:',/,
     $       ' STORAGE FOR CONSERVATION PLOTS INSUFFICIENT',/,
     $       ' INCREASE PARAMETER NGRPHS = ', I4)               
         STOP 
      ENDIF 
C
C     SAVE INITIAL VALUES FOR SUBSEQUENT NORMALIZATION
C
      IF (ICTR .EQ. 1) THEN
         TOTEN = TOTE
         PEN   = TPE
         TMSSN = TMASS
C
C        FIND CLOSEST GRIDPOINT FOR HEIGHT TIMESERIES
C        AND CHECK FOR CORRECTNESS
C
         ILON = 1 + NINT(PLON/360.0*NLON)
         IF ((ILON .LT. 1) .OR. (ILON .GT. NLON)) THEN
            WRITE(0,*) ' NRGTCS: OUT OF BOUND LONGITUDE'
            STOP 1
         ENDIF
         PLONG = GLON(ILON)/PI*180.0
         IF (ABS(PLONG-PLON) .GT. 360.0/NLON) THEN
            WRITE(0,*) ' NRGTCS: ERROR IN GRIDPOINT CHOICE'
            STOP 1
         ENDIF
         JLAT = NINT((90.0-PLAT)/180.0*NLAT)
         IF ((JLAT .LT. 1) .OR. (JLAT .GT. NLAT)) THEN
            WRITE(0,*) ' NRGTCS: OUT OF BOUND LATITUDE'
            STOP 1
         ENDIF
         PLATG = GLAT(JLAT)/PI*180.0
         IF (ABS(PLATG-PLAT) .GT. 180.0/NLAT) THEN
            WRITE(0,*) ' NRGTCS: ERROR IN GRIDPOINT CHOICE'
            STOP 1
         ENDIF
         WRITE(6,906) PLATG, PLONG
  906 FORMAT (/, ' INITIAL DATA FOR NORMALIZATION',/,
     $           ' GRIDPOINT LATITUDE    = ',F8.2,/,
     $           ' GRIDPOINT LONGITUDE   = ',F8.2)
C
      ENDIF
C
      WRITE(6,986) TMASS,TOTE,TAPE,TAM,TPE,H(ILON,JLAT)
  986 FORMAT (/, 
     $           ' MASS                  = ',1PE16.9,/,
     $           ' TOTAL ENERGY          = ',1PE16.9,/,
     $           ' VORTICITY             = ',1PE16.9,/,
     $           ' DIVERGENCE            = ',1PE16.9,/,
     $           ' POTENTIAL ENSTROPHY   = ',1PE16.9,/,
     $           ' GRIDPOINT HEIGHT      = ',1PE16.9)
C                                                                               
C     STORE RELATIVE ERRORS FOR LATER PLOTTING
C
      IF (TOTEN .NE. 0.0) THEN
         TE  (ICTR) = (TOTE-TOTEN)/TOTEN 
      ELSE
         TE  (ICTR) = TOTE
      ENDIF
      IF (PEN .NE. 0.0) THEN
         PE  (ICTR) = (TPE-PEN)/PEN
      ELSE
         PE  (ICTR) = TPE
      ENDIF
C
C     GLOBAL VORTICITY APPROXIMATELY ZERO
C
      APE (ICTR) = TAPE
C
C     GLOBAL DIVERGENCE APPROXIMATELY ZERO
C
      AAM (ICTR) = TAM
C
      IF (TMSSN .NE. 0.0) THEN
         TMSS(ICTR) = (TMASS-TMSSN)/TMSSN
      ELSE
         TMSS(ICTR) = TMASS
      ENDIF
C
C     GRIDPOINT FOR HEIGHT TIMESERIES
C
      HTS(ICTR) = H(ILON,JLAT)
C
C     SAVE MODEL TIME
C
      TME (ICTR) = TAU                                                          
C                                                                               
C     END ENERGETICS ANALYSIS ON GRID
C
      ENDIF
C
      RETURN                                                                    
      END                                                                       
