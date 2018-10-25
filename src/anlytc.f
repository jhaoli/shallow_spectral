      SUBROUTINE ANLYTC(L3D,LN,TIME,CASE,
     $                  UICLL,VICLL,PICLL,DICLL,EICLL)
C
C THIS PROCEDURE RETURNS THE ANALYTIC SOLUTION FOR THE TEST CASES.
C IT IS CALLED FOR THE INITIALIZATION OF
C THE PROGNOSTIC AND ANALYTIC FIELDS AND DURING ERROR ANALYSIS. 
C THE RESULTS ARE RETURNED IN THE ARRAYS UICLL,VICLL,PICLL,DICLL, 
C AND EICLL. THE ROUTINE USES MANY VALUES FROM THE COMMON
C BLOCKS /CONST/ AND /CONST2/ WHICH WERE COMPUTED IN ROUTINES 
C INPUT AND INIT.
C                                                                              
C CALLED BY: ERRANL, STSWM
C CALLS: DBUBF, EVAL, GLAT, GLON, INPTP, ROTATE
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
C Plotting Data
      INCLUDE 'complt.i'
C Constants  
      INCLUDE 'consts.i'
C Initial Conditions
      INCLUDE 'finit.i'
      INCLUDE 'case4.i'
C
C------------ Arguments ------------------------------------------------
C
C     Input
C
C NUMBER OF TIMELEVELS OF DIVERGENCE,VORTICITY
C          (=1 WHEN CALLED FROM SUBROUTINE ERRANL)
C          (=LVLS WHEN CALL FROM MAIN PROGRAM)
      INTEGER L3D
C CURRENT TIMELEVEL OF DIVERGENCE,VORTICITY
      INTEGER LN
C MODEL TIME IN SECONDS
      REAL TIME
C TEST CASE NUMBER
      INTEGER CASE
C
C     Output
C
C EASTWARD WIND VELOCITY ON GRID
      REAL UICLL(NLON+2,NLAT)
C NORTHWARD WIND VELOCITY ON GRID
      REAL VICLL(NLON+2,NLAT)
C HEIGHT ON GRID
      REAL PICLL(NLON+2,NLAT)
C DIVERGENCE ON GRID
      REAL DICLL(NLON+2,NLAT,L3D)
C VORTICITY ON GRID
      REAL EICLL(NLON+2,NLAT,L3D)
C
C------ Local Parameters -----------------------------------------------
C
      REAL PI
      PARAMETER (PI=3.141592653589793)
C
C------ Local Variables ------------------------------------------------
C
      REAL DIV,ZETA,PHI,U,V
C     TEST CASE 1 TEMPORARIES
      REAL RADIUS, DIST
C     LONGITUDES
      REAL RLONA, RRLONA, RLON
      INTEGER I
C     LATITUDES
      INTEGER J
      REAL RLATA, RRLATA, RLAT
C     TEST CASE 4 TEMPORARIES
      REAL TMSHFT, AI, A2I, SNJ, CSJ, SRCSJ, TMPRY, TMPRY2, DEN,
     $     AACSJI, CORR, BIGUBR, DBUB, C, PSIB, DCDM, DCDL, D2CDM,
     $     D2CDL, TMP1, TMP2, DKDM, DKDL, D2KDM, D2KDL, DLON, SINT,
     $     COST
C
C     SPECTRAL COEFFICIENTS FROM NETCDF FILE
      COMPLEX ZETASC(NALPH),DIVSC(NALPH),PHISC(NALPH)
C
C     SPECTRAL TRUNCATION PARAMETERS OF FILE
C
      INTEGER MMH,NNH,KKH
      INTEGER LDIAG(0:MAXH,2)
C
C----- External Functions ----------------------------------------------
C
C     ZONAL FLOW FUNCTIONS
C
      EXTERNAL DBUBF
      REAL DBUBF
C
C     LAT./ LONG. GRID
C
      EXTERNAL GLAT, GLON
      REAL GLAT, GLON
C
C----- Executable Statements -------------------------------------------
C                                                                               
      IF (CASE .EQ. 1) THEN
C
C-------- INITIAL CONDITION #1 ---------------------------------------
C
C     COPY DATA FOR U,V,ZETA,DIV SINCE STEADY FLOW FIELD
C     COMPUTE LOCATION OF ADVECTED HEIGHT
C
      RLATA = RLAT0
      RLONA = RLON0 + SU0*COS(RLAT0)*TIME/A
      CALL ROTATE(RLONA,RLATA,-ALPHA,RRLONA,RRLATA)
C
C     SET VIEWPOINT FOR GRAPHICS
C
      POLAT = RRLATA*180.0/PI
      POLNG = RRLONA*180.0/PI
      POROT = ALPHA*180.0/PI-90.0
C
C     SIZE OF FEATURE
C
      RADIUS = A/3.0
C
      DO 41 J=1,NLAT
         RLAT = GLAT(J)
         DO 40 I=1,NLON
            RLON = GLON(I)
C
C           LONGITUDE = RLON = GLON(I) 
C           LATITUDE = RLAT = GLAT(J)
C
            DICLL(I,J,LN) = DIC12(I,J)
            EICLL(I,J,LN) = EIC12(I,J)
            UICLL(I,J)    = UIC12(I,J)
            VICLL(I,J)    = VIC12(I,J)
C
C        CONSTRUCT ADVECTED LOW
C
            DIST = A*ACOS(SIN(RRLATA)*SIN(RLAT) + COS(RRLATA)
     $          *COS(RLAT)*COS(RLON-RRLONA))
            IF (DIST .LE. RADIUS) THEN
               PICLL(I,J) = PHI0/2.0*(1.0 + COS(PI*DIST/RADIUS))
            ELSE
               PICLL(I,J) = 0.0
            ENDIF
C
   40    CONTINUE
   41 CONTINUE
C
      ELSEIF ((CASE .EQ. 2) .OR. (CASE .EQ. 3)) THEN
C
C-------- INITIAL CONDITION #2, #3 ------------------------------------
C                                                                               
C     SET VIEWPOINT FOR GRAPHICS
C
      IF (CASE .EQ. 2) THEN
         POLAT = 90.0 - ALPHA/PI*180.0
         POLNG = 180.0
      ENDIF
C
C     COPY INITIAL DATA, SINCE STEADY STATE SOLUTION
C
      DO 51 J=1,NLAT
         DO 50 I=1,NLON  
C
C        LONGITUDE = RLON = GLON(I)
C        LATITUDE = RLAT = GLAT(J)
C
            DICLL(I,J,LN) = DIC12(I,J)
            EICLL(I,J,LN) = EIC12(I,J)
            UICLL(I,J)    = UIC12(I,J)
            VICLL(I,J)    = VIC12(I,J)
            PICLL(I,J)    = PIC12(I,J)
   50    CONTINUE 
   51 CONTINUE
C
      ELSEIF (CASE .EQ. 4) THEN
C
C-------- INITIAL CONDITION #4 ---------------------------------------
C                                                                               
C     CALCULATE ANALYTIC SOLUTION TO FORCED NONLINEAR PROBLEM                   
C                                                                               
C
C     LONGITUDINAL CHANGE OF LOW IN BASIC FLOW
C
      TMSHFT =  SU0*TIME/A                                                      
C
C     SET VIEWPOINT FOR GRAPHICS
C
      POLNG = (RLON0+TMSHFT)*180.0/PI
      POLAT = RLAT0*180.0/PI
      DO 101 J=1,NLAT
C
C        GAUSSIAN LATITUDE
C
         RLAT = GLAT(J)
C
C        TEMPORARY VARIABLES INDEPENDENT OF LONGITUDE
C
         AI     = 1.0/A  
         A2I    = 1.0/(A*A) 
         SNJ    = SIN(RLAT)
         CSJ    = COS(RLAT)*COS(RLAT)
         SRCSJ  = COS(RLAT)
         TMPRY  = TAN(RLAT)
         TMPRY2 = TMPRY*TMPRY
         DEN    = 1.0/COS(RLAT)
         AACSJI = 1.0/(A*A*CSJ)
         CORR   = 2.0*OMEGA*SNJ
C
C        NONLINEAR STEADY ZONAL FLOW
C
         BIGUBR = UCON(J)*SRCSJ
         DBUB   = DBUBF(RLAT)
C                                                                               
         DO 100 I=1,NLON  
            RLON = GLON(I)
C
C           LONGITUDE = RLON = GLON(I)
C           LATITUDE = RLAT = GLAT(J)
C
C           COMPUTE LOCATION OF TRANSLATING LOW
C
            C      = SIN(RLAT0)*SNJ + COS(RLAT0)*SRCSJ* 
     $                                   COS(RLON-TMSHFT-RLON0)             
            PSIB   = ALFA*EXP(-SIGMA*((1.0-C)/(1.0+C)))   
C                                                                               
C        COMPUTE PARTIAL DERIVATIVES OF C
C
            DCDM   = SIN(RLAT0) - COS(RLON-TMSHFT-RLON0)* 
     $               COS(RLAT0)*TMPRY
            DCDL   = -COS(RLAT0)*SRCSJ*SIN(RLON-TMSHFT-RLON0)                  
            D2CDM  = -COS(RLAT0)*COS(RLON-TMSHFT-RLON0)*                       
     $                      (1.0 + TMPRY2)/SRCSJ                                
            D2CDL  = -COS(RLAT0)*SRCSJ*COS(RLON-TMSHFT-RLON0)                  
C                                                                               
C        COMPUTE PARTIAL DERIVATIVES OF PSI BAR
C
            TMP1   = 2.0*SIGMA*PSIB/((1.0 + C)**2) 
            TMP2   = (SIGMA - (1.0 + C))/((1.0 + C)**2) 
            DKDM   = TMP1*DCDM   
            DKDL   = TMP1*DCDL  
            D2KDM  = TMP1*(D2CDM + 2.0*(DCDM**2)*TMP2) 
            D2KDL  = TMP1*(D2CDL + 2.0*(DCDL**2)*TMP2) 
C                                                                               
C           ANALYTIC SOLUTIONS 
C
            DICLL(I,J,LN) = 0.0
            EICLL(I,J,LN) = D2KDL*AACSJI + CORR - DBUB*AI
     $                  + (CSJ*D2KDM - 2.0*SNJ*DKDM)*A2I 
            UICLL(I,J)    = BIGUBR*DEN - SRCSJ*AI*DKDM 
            VICLL(I,J)    = (DKDL*AI)*DEN 
            PICLL(I,J)    = PHICON(J)+CORR*PSIB/GRAV
  100    CONTINUE 
  101 CONTINUE
C                                                                               
      ELSEIF ((CASE .EQ. 5) .AND. (TIME .EQ. 0.0)) THEN
C
C-------- INITIAL CONDITION #5 ---------------------------------------
C
C     ZONAL FLOW OVER ISOLATED MOUNTAIN
C
C     COPY INITIAL DATA 
C
      DO 121 J=1,NLAT
         DO 120 I=1,NLON
C
C           LONGITUDE = RLON = GLON(I)
C           LATITUDE = RLAT = GLAT(J)
C
            DICLL(I,J,LN) = DIC12(I,J)
            EICLL(I,J,LN) = EIC12(I,J)
            UICLL(I,J)    = UIC12(I,J)
            VICLL(I,J)    = VIC12(I,J)
            PICLL(I,J)    = PIC12(I,J) 
  120    CONTINUE
  121 CONTINUE
C
      ELSEIF ((CASE .EQ. 6) .AND. (TIME .EQ. 0.0)) THEN
C
C-------- INITIAL CONDITION #6 ---------------------------------------
C
C     ROSSBY-HAURWITZ WAVE
C
C     LONGITUDINAL CHANGE OF FEATURE
C
      DLON = (R*(3+R)*OMG - 2.0*OMEGA)/((1+R)*(2+R))*TIME
C
      DO 151 J=1,NLAT
         RLAT = GLAT(J)
         DO 150 I=1,NLON
            RLON = GLON(I)
C
C           LONGITUDE = RLON = GLON(I)
C           LATITUDE = RLAT = GLAT(J)
C
            SINT = SIN(RLAT)
            COST = COS(RLAT)
            DICLL(I,J,LN) = 0.0
            EICLL(I,J,LN) = 2.0*(OMG+OMEGA)*SINT - (1+R)*(2+R)*SINT*
     $                   K*COST**R*COS(R*(RLON-DLON))
            UICLL(I,J)    = A*OMG*COST + A*K*COST**(R-1)*
     $                   (R*SINT*SINT-COST*COST)*COS(R*(RLON-DLON))
            VICLL(I,J)    = -A*K*R*COST**(R-1)*SINT
     $                    * SIN(R*(RLON-DLON))
            PICLL(I,J)    = PHI0 + (A*A*(PHIA(J)+PHIB(J)
     $                    * COS(R*(RLON-DLON))+PHIC(J)
     $                    * COS(2*R*(RLON-DLON))))/GRAV
  150    CONTINUE
  151 CONTINUE
C
      ELSEIF ((CASE .EQ. 7) .OR. (TIME .GT. 0.0)) THEN
C
C-------- INITIAL CONDITION #7 ---------------------------------------
C
C     USE HIGH RESOLUTION SOLUTION FROM FILE
C
C     READ SPECTRAL COEFFICIENTS BACK IN
C
      CALL INPTP(FNIN,TIME,MAXH,LDIAG,ZETASC,DIVSC,PHISC,
     $              MMH,NNH,KKH)
C
C     PROJECT ONTO IRREGULAR GRID
C
      DO 451 J=1,NLAT
         RLAT = GLAT(J)
         DO 450 I=1,NLON
            RLON = GLON(I)
C
C           LONGITUDE = RLON = GLON(I)
C           LATITUDE = RLAT = GLAT(J)
C
C           COMPUTE VALUES AT LOCATION
C
            CALL EVAL(DIVSC,ZETASC,PHISC,MMH,NNH,KKH,A,ALPHA,
     $             OMEGA,RLON,RLAT,DIV,ZETA,PHI,U,V)
C
            DICLL(I,J,LN) = DIV
            EICLL(I,J,LN) = ZETA
            UICLL(I,J)    = U
            VICLL(I,J)    = V
            PICLL(I,J)    = PHI/GRAV
            IF (FTOPO) THEN
               PICLL(I,J) = PICLL(I,J) + MOUNT (I,J)
            ENDIF
C
  450    CONTINUE
  451 CONTINUE
C
C--------------------LAST CASE ----------------------------------
C
      ELSE
C
         WRITE (0,300) CASE
  300    FORMAT (/,' STSWM: FATAL ERROR IN ANLYTC: ',/,
     $      ' NO ANALYTIC SOLUTION PROVIDED FOR TEST CASE ',I2,/)
         STOP
      ENDIF
C
C----------------------------------------------------------------
C
      RETURN
C
      END                                                                       
