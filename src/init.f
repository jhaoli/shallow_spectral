      SUBROUTINE INIT
C                                                                             
C THIS ROUTINE INITIALIZES VARIABLES FOR THE VARIOUS TEST CASES.
C FOR STEADY-STATE FLOWS, THE INFORMATION IS STORED IN COMMON
C BLOCK /CONST2/, ARRAYS UIC12,VIC12,PIC12,DIC12 AND EIC12.
C THE CASE IS DETERMINED BY THE VARIABLE ICOND:
C    CASE 1: ADVECTION EQUATION FOR SOLID BODY FLOW
C    CASE 2: SOLID BODY ROTATION STEADY STATE FLOW
C    CASE 3: JETSTREAM STEADY STATE FLOW
C    CASE 4: FORCED LOW IN JETSTREAM
C    CASE 5: ZONAL FLOW OVER ISOLATED MOUNTAIN
C    CASE 6: ROSSBY-HAURWITZ WAVE
C    CASE 7: REAL DATA (500MB) TEST CASE
C                                                                             
C CALLED BY: STSWM
C CALLS: BUBFNC, D01AHE, FU, FUNC2, GLAT, GLON, ROTATE, 
C        SHTRNS, US, ZD
C
C REVISIONS:
C 7-13-92 CHANGE TO CCM CODING CONVENTIONS (R. JAKOB)
C
C---- Model Parameters -------------------------------------------------
C
      INCLUDE 'params.i'
C
C---- Common Blocks ----------------------------------------------------
C
C CONSTANTS & TIMESTEPS
      INCLUDE 'consts.i'
C INITIAL CONDITIONS
      INCLUDE 'finit.i'
      INCLUDE 'case4.i'
C
C------ Local Parameters -----------------------------------------------
C
      REAL PI
      PARAMETER (PI=3.141592653589793)
C
C------ Local Variables ------------------------------------------------
C
C     TEMPORARY ARRAYS FOR COMPUTING INITIAL DIVERGENCE,VORTICITY
C
      COMPLEX DIVSC(NALP), ZETASC(NALP)
C
C     TEMPORARIES
C
      INTEGER I,J
      REAL SINA, COSA, RLON, COSL, RLAT, ETAAMP, UBAR, PHIAMP,
     $     SINL, SINT, COST, ROTLON, ROTLAT, PHITMP, RELERR,
     $     RADIUS, DIST, MOUNTA
C
C     NAG INTEGRATION ROUTINE ARGUMENTS
C
      INTEGER NLIMIT, IFAIL, NPTS
C
C----- External Functions ----------------------------------------------
C
C     Actual function arguments for balance equations
C
      EXTERNAL BUBFNC, FU, FUNC2, US
      REAL BUBFNC, FU, FUNC2, US
C
C     GRID LATITUDE AND LONGITUDE
C
      EXTERNAL GLAT, GLON
      REAL GLAT, GLON
C
C     QUADRATURE ROUTINE
C
      EXTERNAL D01AHE
      REAL D01AHE
C                                                                               
C----- Executable Statements -------------------------------------------
C
C     DETERMINE WHICH INITIAL CONDITION TO USE                                  
C                                                                               
      IF (ICOND .EQ. 1) THEN
C
C-------------------INITIAL CONDITION 1 ---------------------------
C
C     ANALYTIC SPECIFICATION FOR ADVECTION OF COSINE BELL
C     SOLUTION AS SPECIFIED BY WILLIAMSON AND
C     RASCH, 1989, MON. WEATHER REVIEW
C
      WRITE (6,70) ALPHA
   70 FORMAT(/,' TEST CASE #1: ADVECTION OF COSINE BELL',
     $       /,' ROTATED BY AN ANGLE ALPHA = ',F5.3)
C
C     INITIAL LOCATION AND AMPLITUDE
C
      RLON0 = -90.0*(PI/180.0)
      RLAT0 = + 0.0*(PI/180.0)
      PHI0  = 1000.0
C
C     ONLY EXPLICIT TIMESTEPPING WITH ADVECTION EQUATION !
C
      IF (SITS) THEN
         WRITE(0,75)
   75 FORMAT(/,' STSWM: CANNOT RUN TEST 1 WITH SEMI-IMPLICIT',
     $       /,' TIMESTEPPING. CHANGE PARAMETER SITS = .FALSE.',
     $       /,' IN NAMELIST INPUT FILE')
         STOP
      ENDIF
C
C     CHOOSE VELOCITY FOR ONE ROTATION IN 12 DAYS
C
      SU0  = 2.0*PI*A/(3600.0*24*12)
      SINA = SIN(ALPHA)
      COSA = COS(ALPHA)
      ETAAMP = 2.0*(SU0/A + OMEGA)
      DO 85 I=1,NLON
         RLON = GLON(I)
         DO 80 J=1,NLAT
            RLAT = GLAT(J)
C
C           LONGITUDE = RLON = GLON(I)
C           LATITUDE = RLAT = GLAT(J)
C
C           SET UP STEADY FLOW FIELD
C
            IF (ALPHA .NE. 0.0) THEN
               UIC12(I,J) = SU0*(COS(RLAT)*COSA+COS(RLON)*
     $                      SIN(RLAT)*SINA)
               VIC12(I,J) = -SU0*SIN(RLON)*SINA
            ELSE
              UIC12(I,J) = SU0*COS(RLAT)
              VIC12(I,J) = 0.0
            ENDIF
C
            DIC12(I,J) = 0.0
            EIC12(I,J) = ETAAMP*(-COS(RLON)*COS(RLAT)*SINA
     $                   +SIN(RLAT)*COSA)
   80    CONTINUE
   85 CONTINUE
C
      RETURN
C
      ELSEIF (ICOND .EQ. 2) THEN
C
C-------------------INITIAL CONDITION 2 ---------------------------        
C                                                                               
C     ANALYTIC SPECIFICATION OF U, V, PHI, DIV, ZETA FIELD ON GAUSSIAN          
C     GRID. STEADY STATE SOLUTION (WILLIAMSON AND BROWNING, 1973 JAM)           
C                                                                               
      WRITE (6,90) ALPHA
   90 FORMAT(/,' TEST CASE #2: STEADY STATE NONLINEAR GEOSTROPHIC FLOW',
     $       /,' ROTATED BY AN ANGLE ALPHA = ',F5.3)
C
C     INITIAL LOCATION AND AMPLITUDE
C
      RLON0 = -90.0*(PI/180.0)
      RLAT0 = + 0.0*(PI/180.0)
C
      UBAR = (2.0*PI*A)/(12.0*24.0*3600.0)
      PHI0 = 2.94E4
C
      SINA = SIN(ALPHA)
      COSA = COS(ALPHA)
      ETAAMP = 2.0*(UBAR/A + OMEGA)
      PHIAMP = A*OMEGA*UBAR + (UBAR**2)/2.0
      DO 100 I=1,NLON                                                           
         RLON = GLON(I)
         SINL = SIN(RLON)
         COSL = COS(RLON)
         DO 95 J=1,NLAT                                                         
            RLAT = GLAT(J)
C
C           LONGITUDE = RLON = GLON(I)
C           LATITUDE = RLAT = GLAT(J)
C
            SINT = SIN(RLAT)
            COST = COS(RLAT)
            UIC12(I,J) = UBAR*(COST*COSA + COSL*SINT*SINA)
            VIC12(I,J) = - UBAR*SINL*SINA
            PIC12(I,J) = (PHI0-PHIAMP*(- COSL*COST*SINA + 
     $                   SINT*COSA)**2)/GRAV
            DIC12(I,J) = 0.0
            EIC12(I,J) = ETAAMP*(- COSL*COST*SINA + SINT*COSA)
   95    CONTINUE                                                               
  100 CONTINUE                                                                  
C                                                                               
      RETURN
C
      ELSE IF (ICOND .EQ. 3) THEN
C                                                                               
C--------------------INITIAL CONDITION 3 --------------------------        
C                                                                               
C     INITIAL U SPECIFIED AS BUMP (INFINITELY DIFFERENTABLE) FUNCTION           
C     V=0; SOLVE FOR PHI BY INTEGRATING 1 DIMENSIONAL BALANCE EQUATION          
C     SEE PAPER BY BROWNING ET. AL., (MONTHLY WEATHER REVIEW, 1989)
C                                                                               
      WRITE (6,120) ALPHA
  120 FORMAT(/,' TEST CASE #3: STEADY STATE NONLINEAR GEOSTROPHIC FLOW',
     $       /,' WITH COMPACT SUPPORT',
     $       /,' ROTATED BY AN ANGLE ALPHA = ',F5.3)
C
      PHI0   = 2.94E4
      NLIMIT = -1                                                               
      IFAIL  =  0                                                               
      DO 150 J=1,NLAT                                                           
         RLAT = GLAT(J)
         DO 145 I = 1, NLON
            RLON = GLON(I)
C
C           LONGITUDE = RLON = GLON(I)
C           LATITUDE = RLAT = GLAT(J)
C
C           COMPUTE COORDINATES IN ROTATED SYSTEM
C
            CALL ROTATE(RLON,RLAT,ALPHA,ROTLON,ROTLAT)
C                                                                               
C        CALL NUMERICAL INTEGRATION PROCEDURE FROM NAG FOR PHI                  
C        WHAT FOLLOWS IS A QUICK DESCRIPTION OF THE ARGUMENT LIST SO THE        
C        USER CAN REPLACE THIS ROUTINE WITH ANOTHER ONE IF NECESSARY            
C                                                                               
C        FUNCTION D01AHE (A, B, EPSR, NPTS, RELERR, F, NLIMIT, IFAIL)           
C           A      - SPECIFIES THE LOWER LIMIT OF INTEGRATION                   
C           B      - SPECIFIES THE UPPER LIMIT OF INTEGRATION                   
C           EPSR   - SPECIFIES THE RELATIVE ACCURACY REQUIRED                   
C           NPTS   - NUMBER OF POINTS AT WHICH TO EVALUATE THE INTEGRAL         
C           RELERR - CONTAINS ROUGH ESTIMATE OF RELATIVE ERROR ON EXIT          
C           F      - REAL FUNCTION, SUPPLIED BY THE USER                        
C           NLIMIT - SPECIFIES A LIMIT TO NUMBER OF FUNCTION EVALUATIONS        
C                       NLIMIT.LE.0 => LIMIT OF 10,000                          
C           IFAIL  - MUST BE PREASSIGNED WHEN ROUTINE IS CALLED                 
C                       CONTAINS 0 ON OUTPUT IF NO ERROR OCCURRED               
C                                                                               
            PHITMP = D01AHE(-0.5*PI, ROTLAT, 100.0*EPS, NPTS, 
     $                RELERR, FU, NLIMIT, IFAIL)                                
            IF (IFAIL .NE. 0) THEN  
               WRITE (0,140) IFAIL 
  140          FORMAT(/,' STSWM: FATAL ERROR IN ROUTINE INIT:',/,
     $             ' FAILURE IN NAG INTEGRATION ROUTINE D01AHE',/,
     $             ' IFAIL = ',I4)
               STOP 
            ENDIF 
C                                                                               
C           ROTATE FIELD VARIABLES
C
            IF (ALPHA .NE. 0.0) THEN
               UIC12(I,J) = US(ROTLAT)*(COS(ALPHA)*SIN(ROTLON)
     $                      *SIN(RLON)+COS(RLON)*COS(ROTLON))
               VIC12(I,J) = US(ROTLAT)*(COS(ALPHA)*COS(RLON)
     $                      * SIN(ROTLON)*SIN(RLAT) 
     $                      - COS(ROTLON)*SIN(RLON)*SIN(RLAT) 
     $                      - SIN(ALPHA)*SIN(ROTLON)*COS(RLAT))
               PIC12(I,J) = (PHI0-PHITMP)/GRAV
            ELSE
C
C              NO ROTATION -> ORIGINAL FIELD
C
               UIC12(I,J) = US(RLAT)
               VIC12(I,J) = 0.0
               PIC12(I,J) = (PHI0-PHITMP)/GRAV
            ENDIF
  145    CONTINUE                                                               
  150 CONTINUE                                                                  
C                                                                               
C     GET SPECTRAL COEFFICIENTS FOR ZETA AND DIV                  
C     (STEADY STATE)
C                                                                               
      CALL ZD(UIC12,VIC12,DIVSC,ZETASC)
C                                                                               
C     INVERSE TRANSFORM ZETA,DIV FOR INITIAL CONDITION INFORMATION 
C     DESTROYS ZETASC, DIVSC !
C                                                                               
      CALL SHTRNS (1,1,+1,EIC12,ZETASC)        
      CALL SHTRNS (1,1,+1,DIC12,DIVSC)
C                                                                               
      RETURN
C
      ELSE IF (ICOND .EQ. 4) THEN
C                                                                               
C-----------------------INITIAL CONDITION 4 -----------------------        
C
C     STEADY STATE NONLINEAR ZONAL GEOSTROPHIC FLOW WITH COMPACT
C     SUPPORT, SEE BROWNING ET AL, MONTHLY WEATHER REVIEW, 1989
C
      WRITE(6,200) 
  200 FORMAT(/,' TEST CASE #4: FORCED NONLINEAR SYSTEM WITH',
     $       ' ADVECTING LOW')
C
      IF (.NOT. FORCED) THEN
         WRITE(0,205)
  205    FORMAT(/,' STSWM: FORCING TERMS MUST BE INCLUDED FOR',
     $          /,' TEST 4. SET PARAMETER FORCED = .TRUE. IN',
     $          /,' NAMELIST INPUT FILE')
         STOP
      ENDIF
C                                                                               
C     FORCED NONLINEAR SOLUTION                                                 
C                                                                               
C     DETERMINE INITIAL CONDITION ON TIME DEPENDENT VARIABLES FOR FIRST         
C     TWO TIME STEPS (NO TIME TRUNCATION ERROR TO START LEAPFROG PROC.)         
C                                                                               
C     CONSTANTS FOR ANALYTIC STREAM FUNCTION (FORCED CASE)                      
C
C     BASIC ZONAL FLOW AMPLITUDE
C                                                                               
      SU0    = 20.0                                                             
      PHI0   = 1.0E5
C
C     USE A TRANSLATING LOW INSTEAD OF A TRANSLATING HIGH 
C     (PRODUCE PLOTS AS IN PAPER)
C
C
C     AMPLITUDE OF LOW
C     TO REMOVE LOW, SET TO ZERO
C
      ALFA   = -0.03*(PHI0/(2.0*OMEGA*SIN(PI/4.0))) 
C*****ALFA   = 0.0
C
C     INITIAL POSITION OF LOW (GREENWICH)
C
      RLON0  = 0.0*(PI/180.0) 
      RLAT0  = +45.0*(PI/180.0)                                                 
C
C     AREAL EXTENT OF LOW
C
      SIGMA  = (2.0*A/1.0E6)**2                                                 
      NPWR   =  14                                                              
C                                                                               
C     BALANCED PHI (PHICON), UCON, VCON IS USED IN ANLYTC FOR 
C     GEOPOTENTIAL THAT BALANCES STEADY ZONAL FLOW
C                                                                               
      NLIMIT = -1                                                               
      IFAIL  = 0   
      DO 210 J=1,NLAT                                                           
         RLAT = GLAT(J)
C
C        LONGITUDE = RLON = GLON(I)
C        LATITUDE = RLAT = GLAT(J)
C
C        SOLVE NON-LINEAR BALANCE EQUATION FOR PHI;
C                                                                               
         PHICON(J) = D01AHE(-0.5*PI, RLAT, 100*EPS,           
     $               NPTS, RELERR, FUNC2, NLIMIT, IFAIL)                        
C
         IF (IFAIL .NE. 0) THEN 
            WRITE (0,140) IFAIL
            STOP                                                                
         ENDIF                                                                  
         PHICON(J) = (PHI0 - PHICON(J))/GRAV
         UCON(J) = BUBFNC(RLAT)
         VCON(J) = 0.0
  210 CONTINUE                                                                  
C                                                                               
      RETURN
C
      ELSE IF (ICOND .EQ. 5) THEN
C
C-------------------INITIAL CONDITION 5 ---------------------------
C
C     ZONAL FLOW OVER AN ISOLATED MOUNTAIN AS USED BY TAKACS.
C     ANALYTIC SPECIFICATION OF U, V, PHI, DIV, ZETA FIELD ON GAUSSIAN
C     GRID. STEADY STATE SOLUTION (WILLIAMSON AND BROWNING, 1973 JAM)
C
      WRITE (6,280) 
  280 FORMAT(/,' TEST CASE #5: ZONAL FLOW OVER AN ISOLATED MOUNTAIN')
C
C     SET MOUNTAIN SURFACE
C
      FTOPO = .TRUE.
      MOUNTA = 2000.0
      RADIUS = PI/9.0
      DO 50 J = 1, NLAT
         RLAT = GLAT(J)
         DO 40 I = 1, NLON
            RLON = GLON(I)
            DIST = SQRT((RLON -  1.5*PI)**2 + (RLAT - PI/6.0)**2)
            IF (DIST .LT. RADIUS) THEN
               MOUNT(I,J) = MOUNTA*(1.0 - DIST/RADIUS)
            ELSE
               MOUNT(I,J) = 0.0
            ENDIF
C
C           TEMPORARY COPY FOR SPECTRAL TRANSFORM
C
            PIC12(I,J)=MOUNT(I,J)  
   40    CONTINUE
   50 CONTINUE
C
C     COMPUTE SPECTRAL COEFFICIENTS
C
      CALL SHTRNS(1,1,-1,PIC12,TOPOSC)
C
C     INITIAL CONDITIONS
C
      PHI0 = 5960.0
      UBAR = 20.0
C
      SINA = 0.0
      COSA = 1.0
      ETAAMP = 2.0*(UBAR/A + OMEGA)
      PHIAMP = A*OMEGA*UBAR + (UBAR**2)/2.0
      DO 295 I=1,NLON
         RLON = GLON(I)
         SINL = SIN(RLON)
         COSL = COS(RLON)
         DO 290 J=1,NLAT
            RLAT = GLAT(J)
C
C           LONGITUDE = RLON = GLON(I)
C           LATITUDE = RLAT = GLAT(J)
C
            SINT = SIN(RLAT)
            COST = COS(RLAT)
            UIC12(I,J) = UBAR*COST
            VIC12(I,J) = 0.0
C
C           FREE SURFACE HEIGHT (INCLUDE MOUNTAINS)
C
            PIC12(I,J) = PHI0-PHIAMP*SINT**2/GRAV
            DIC12(I,J) = 0.0
            EIC12(I,J) = ETAAMP*SINT
  290    CONTINUE
  295 CONTINUE
C
      RETURN
C
      ELSE IF (ICOND .EQ. 6) THEN
C
C-----------------------INITIAL CONDITION 6 -----------------------
C
C     ROSSBY-HAURWITZ WAVE AS USED BY PHILIPS IN
C     MONTHLY WEATHER REVIEW, 1959
C
      R = 4
      K = 7.848E-6
      OMG = 7.848E-6
      PHI0 = 8000.0
      WRITE(6,300) R
  300 FORMAT(/,' TEST CASE #6: ROSSBY-HAURWITZ WAVE, WAVENUMBER ',
     $       I2)
C
C     COMPUTE LATITUDE-DEPENDENT FACTORS FOR GEOPOTENTIAL
C     PHIA(NLAT),PHIB(NLAT) AND PHIC(NLAT)
C
      DO 310 J=1,NLAT
         RLAT = GLAT(J)
C
C        LONGITUDE = RLON = GLON(I)
C        LATITUDE = RLAT = GLAT(J)
C
         COST = COS(RLAT)
         PHIA(J) = 0.5*OMG*(2.0*OMEGA+OMG)*COST*COST +
     $             0.25*K*K*COST**(2*R) *
     $             ((R+1)*COST*COST+(2*R*R-R-2) - 
     $             2.0*R*R/(COST*COST))
         PHIB(J) = (2.0*(OMEGA+OMG)*K)/((R+1)*(R+2))*
     $             COST**R*((R*R+2*R+2)-(R+1)**2*COST*COST)
         PHIC(J) = 0.25*K*K*COST**(2*R)*
     $             ((R+1)*COST*COST-(R+2))
  310 CONTINUE
C
      RETURN
C
      ELSE IF (ICOND .EQ. 7) THEN
C
C-----------------------INITIAL CONDITION 7 -----------------------
C
C     REAL DATA TEST CASE
C
      WRITE(6,400) 
  400 FORMAT(/,' TEST CASE #7: 500 MB GEOPOTENTIAL FROM ECMWF ',
     $       'ANALYSIS')
C
      RETURN
C                                                                               
C----------------------------------------------------------------------
C
      ELSE
C
         WRITE (0,900) ICOND
  900    FORMAT(/,' STSWM: FATAL ERROR IN SUBROUTINE INIT:',/,               
     $       ' MISSING SPECIFICATION OF ',                
     $       'INITIAL CONDITION CASE = ',I3)
         STOP                                                                   
      ENDIF                                                                     
C
      END                                                                       
