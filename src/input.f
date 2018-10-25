      SUBROUTINE INPUT                                                          
C                                                                              
C THIS ROUTINE INPUTS OR DETERMINES NECESSARY CONSTANTS, CALLS              
C THE ROUTINES WHICH CALCULATE GAUSSIAN LATITUDES, WEIGHTS, BASIS           
C FUNCTIONS (ASSOCIATED LEGENDRE POLYNOMIALS) AND DERIVATIVES, AND          
C CALLS ROUTINES TO DO THE INITIAL SETUP FOR THE FAST FOURIER               
C TRANSFORM PROCEDURE (IN THIS CASE SET99).                                
C
C CALLED BY: STSWM
C CALLS: CALP, CEPS, EPSLON, GLAT, GLON, GLATS, SET99
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
C TRANSFORM ARRAYS
      INCLUDE 'trnsfm.i'
C CONSTANTS & TIMESTEPS
      INCLUDE 'consts.i'
C PLOT DATA
      INCLUDE 'complt.i'
C
C------ Local Parameters -----------------------------------------------
C
      REAL PI
      PARAMETER (PI=3.141592653589793)
C
C------ Local Variables ------------------------------------------------
C
C*    The following arrays contain information required for calculation*        
C*    of the Legendre Polynomials and derivatives using the method     *        
C*    proposed by Belousov (1962).  The length of the three recurrence *        
C*    coefficient matrices used in this procedure is given by the      *        
C*    relation LRM + 1 = ((MM+1)*(KK+1)-(MM**2+MM)/2).  The length of  *        
C*    each row is stored in the array LROW(0:KK,1).  If the Belousov   *        
C*    algorithm is not selected, LRM can be set to 0 to save static    *        
C*    storage.  This is also where we store EPSIL which is used in a   *        
C*    simpler polynomial recurrence relationship as well as in the     *        
C*    calculation of the polynomial derivatives.                                
C*    The arrays are passed from subroutine CEPS to CALP.
C*                                                                     *        
      REAL CMN(LRM+1), DMN(LRM+1), EMN(LRM+1), EPSIL(NALP)        
      INTEGER LROW(0:KK,2)
C
C     TEMPORARIES
C
      INTEGER I,J,N,NL,IGROUP,LATMIN
      REAL RLAT
C
C----- External Functions ----------------------------------------------
C
C     MACHINE EPSILON FUNCTION
C
      EXTERNAL EPSLON
      REAL EPSLON
C
C     LAT./LONG. GRID
C
      EXTERNAL GLAT, GLON
      REAL GLAT, GLON
C                                                                     *        
C----- Namelist Input Blocks -------------------------------------------
C                                                                               
      NAMELIST /EXPDEF/ DT, SITS, AFC, TAUE, TAUO, GPHFRQ, SPCFRQ,
     $                  FORCED, MOMENT, EGYFRQ, ERRFRQ, ICOND
C
      NAMELIST /PHYVAR/ A, OMEGA, GRAV, HDC, ALPHA
C
      NAMELIST /PLTDEF/ LCONT, LPSP, LOP, LCP, LG, LU, LV, LZ, LD, LVV,         
     $                  LVVG, LGPHS, POLAT, POLNG, POROT
      NAMELIST /FNAMES/ FNIN, FNOUT
C                                                                               
C----- Executable Statements -------------------------------------------
C                                                                               
C     THE FOLLOWING PHYSICAL CONSTANTS ARE DEFINED AS:                          
C        A      => RADIUS OF THE EARTH                                          
C        OMEGA  => ANGULAR VELOCITY OF THE EARTH                                
C        GRAV   => GRAVITATIONAL CONSTANT                                       
C        HDC    => LINEAR HORIZONTAL DIFFUSION COEFFICIENT                      
C        ALPHA  => ROTATION ANGLE OF COORDINATE SYSTEM
C                                                                               
      A      = 6.37122E6                                                        
      OMEGA  = 7.292E-5                                                         
      GRAV   = 9.80616 
      HDC    = 0.0E00                                                           
      ALPHA  = ATAN(1.0)
C                                                                               
C     OTHER CONSTANTS/VARIABLES REQUIRED FOR MODEL INTEGRATION:                 
C        DT     => TIME STEP IN SECONDS                                         
C        AFC    => ASSELIN FILTER PARAMETER                                     
C        TAU    => CURRENT MODEL TIME (IN HOURS)                                
C        TAUO   => FREQUENCY OF PRINTED OUTPUT (IN HOURS)                       
C        TAUE   => END OF MODEL RUN (IN HOURS)                                  
C        SITS   => LOGICAL FLAG FOR SEMI-IMPLICIT TIME DIFFERENCING             
C        GPHFRQ => FREQUENCY OF DETAILED GRAPHICAL OUTPUT (IN HOURS)            
C        EPS    => MACHINE ACCURACY OF FLOATING POINT REPRESENTATION            
C        NSTEP  => CURRENT MODEL TIME STEP (NSTEP=0 INITIALLY)                  
C        FORCED => LOGICAL FLAG .TRUE. IF CALL TO FORCING IN ROUTINE COMP1
C        MOMENT => LOGICAL FLAG .TRUE. FOR MOMENTUM FORCING
C        EGYFRQ => FREQUENCY OF MODEL ENERGETICS EVALUATION (IN HOURS)          
C        ERRFRQ => FREQUENCY OF L2-ERROR ANALYSIS (IN HOURS)
C        SPCFRQ => FREQUENCY OF SPECTRAL ANALYSIS (IN HOURS)
C        ICOND  => TYPE OF INITIAL CONDITIONS (DEFAULT = 2)
C                  (ARE SET UP IN SUBROUTINE INIT)
C        FTOPO  => LOGICAL FLAG .TRUE. IF SURFACE TOPOGRAPHY (CASE 5)
C
      EPS    = EPSLON(1.0)                                                      
      DT     = 2400.0
      AFC    = 0.000                                                            
      TAUO   = 999.0                                                           
      TAUE   = 120.0                                                           
      SITS   = .TRUE.                                                          
      GPHFRQ = 24.0                                                          
      FORCED = .FALSE.                                                          
      MOMENT = .FALSE.
      EGYFRQ = 3.0                                                          
      ERRFRQ = 24.0
      SPCFRQ = 999.0
      ICOND  = 2
      FTOPO  = .FALSE.
C
C     GET EXPERIMENT NUMBER FROM EXECUTION ENVIRONMENT
C       (NONSTANDARD)
C
      CALL GETENV('EXPERIMENT',CHEXP)
      IF (CHEXP .EQ. '    ') THEN
         WRITE(0,590)
  590    FORMAT(/,' STSWM: ''EXPERIMENT'' HAS NOT BEEN DEFINED IN',/,
     $            ' ENVIRONMENT; ASSUMING EXPERIMENT=0000')
         CHEXP = '0000'
      ENDIF
      WRITE(6,670) CHEXP
  670 FORMAT(/,' EXPERIMENT ',A4)
C
C     CHECK TRUNCATION PARAMETERS FOR CONSISTENCY
C
      IF ((KK .LT. NN) .OR. (KK .LT. MM)) THEN
         WRITE(0,600) KK,NN,MM
  600    FORMAT (/,' STSWM: FATAL ERROR IN SUBROUTINE INPUT ',/,
     $          ' TRUNCATION PARAMETER KK MUST BE GREATER OR EQUAL'
     $          ' THAN NN AND MM IN FILE ''PARAMS.i'':',/,
     $          ' KK = ',I4,' NN = ',I4,' MM = ',I4)
         STOP
      ENDIF
      IF (KK .GT. MM+NN) THEN
         WRITE(0,605) KK,MM+NN
  605    FORMAT (/,' STSWM: FATAL ERROR IN SUBROUTINE INPUT ',/,
     $          ' TRUNCATION PARAMETER KK MUST BE LESS THAN OR'
     $          ' EQUAL MM+NN IN FILE ''PARAMS.i'':',/,
     $          ' KK = ',I4,' MM+NN = ',I4)
         STOP
      ENDIF
C
C     DETERMINE TRUNCATION TYPE
C
      IF (NN .EQ. KK) THEN
         IF (NN .EQ. MM) THEN
C
C           TRIANGULAR TRUNCATION MM=NN=KK
C
            WRITE(6,20) MM
            WRITE(STRUNC,'(''T-'',I4)') MM
         ELSE 
C
C           TRAPEZOIDAL TRUNCATION (NN .GT. MM)
C
            WRITE(6,24) MM,NN 
            STRUNC = 'TRAPEZ'
         ENDIF
      ELSE
C
C        NN .LT. KK 
C
         IF (KK .EQ. NN+MM) THEN
            IF (NN .EQ. MM) THEN
C
C              RHOMBOIDAL TRUNCATION
C
               WRITE(6,22) MM
               WRITE(STRUNC,'(''R-'',I4)') MM
            ELSE
C
C           PARALLELOGRAMIC TRUNCATION
C
               WRITE(6,23) MM,NN
               STRUNC = 'PARALL'
            ENDIF
         ELSE
C
C           PENTAGONAL TRUNCATION
C
            WRITE(6,26) MM,NN,KK
            STRUNC = 'PENTAG'
         ENDIF
      ENDIF 
   20 FORMAT (/,' SPECTRAL TRUNCATION TYPE: TRIANGULAR',
     $        /,' M = N = K = ',I4)
   22 FORMAT (/,' SPECTRAL TRUNCATION TYPE: RHOMBOIDAL',
     $        /,' K = N + M, M = N = ',I4)
   23 FORMAT (/,' SPECTRAL TRUNCATION TYPE: PARALLELOGRAMMIC',
     $        /,' K = N + M, M = ',I4,' N = ',I4)
   24 FORMAT (/,' SPECTRAL TRUNCATION TYPE: TRAPEZOIDAL',
     $        /,' N = K > M, M = ',I4,' N = ',I4)
   26 FORMAT (/,' SPECTRAL TRUNCATION TYPE: PENTAGONAL',
     $        /,' M = ',I4,' N = ',I4,' K = ',I4)
C
C     CHECK GRIDPOINT RESOLUTION
C 
      IF (NLON .LT. 3*MM+1) THEN
         WRITE(0,610) NLON,MM
  610    FORMAT (/,' STSWM: FATAL ERROR IN SUBROUTINE INPUT ',/,
     $          ' UNALIASED EVALUATION OF QUADRATIC TERMS REQUIRES',
     $          ' (NLON .GE. 3*MM+1) IN FILE ''PARAMS.i'':',/,
     $          ' NLON = ',I4,' MM = ',I4)
         STOP
      ENDIF
      IF (NN .EQ. KK) THEN
         LATMIN = (3*NN+1)/2
      ELSE
         LATMIN = (3*NN+2*MM+1)/2
      ENDIF
      IF (NLAT .LT. (3*NN+1)/2) THEN
         WRITE(0,620) NLAT,LATMIN
  620    FORMAT (/,' STSWM: FATAL ERROR IN SUBROUTINE INPUT ',/,
     $          ' UNALIASED EVALUATION OF QUADRATIC TERMS REQUIRES',
     $          ' (NLAT .GE. LATMIN) IN FILE ''PARAMS.i'':',/,
     $          ' NLAT = ',I4,' LATMIN = ',I4)
         STOP
      ENDIF
      IF (MOD(NLAT,2) .EQ. 1) THEN
         WRITE(0,630) NLAT
  630    FORMAT (/,' STSWM: FATAL ERROR IN SUBROUTINE INPUT ',/,
     $          ' NLAT MUST BE EVEN TO TAKE ADVANTAGE OF HEMISPHERIC',
     $          ' SYMMETRY PROPERTIES IN FILE ''PARAMS.i'':',/,
     $          ' NLAT = ',I4)
         STOP
      ENDIF
      WRITE (6,640) NLAT, NLON
  640 FORMAT (/,' NUMBER OF GRIDPOINTS IN MODEL'
     $        /,' NORTH-SOUTH GAUSSIAN GRID: NLAT = ',I4,
     $        /,' EAST-WEST EQUIDISTANT GRID: NLON = ',I4)
C
C     PRINT MACHINE EPSILON
C
      WRITE (6,645) EPS
  645 FORMAT (/,' MACHINE EPSILON (1.0 + EPS > 1.0) = ',1PE16.9)
C
C     DETERMINE GUASSIAN LATITUDES AND ASSOCIATED WEIGHTS                       
C                                                                               
C     GET LATITUDES AND WGHTS        
C                                                                               
      CALL GLATS(NLAT, EPS, THTA, WTS)                                      
C                                                                               
C     CALCULATE THE STRUCTURE OF ASSOCIATED LEGENDRE POLYNOMIAL ARRAY,          
C     EPSILON ARRAY, AND OTHER REQUIRED RECURRENCE COEFFICIENT MATRICES         
C                                                                               
      CALL CEPS(CMN,DMN,EMN,EPSIL,LROW,LDIAG)
C                                                                               
C     CALCULATE THE ASSOCIATED LEGENDRE POLYNOMIALS AND DERIVATIVES             
C                                                                               
      CALL CALP(CMN,DMN,EMN,EPSIL,LROW,LDIAG,ALP,DALP)
C                                                                               
C     CALCULATE ANNP1(KK), A2NNP1(0:KK), AND WTACSJ(1:NLAT)                   
C     FOR LATER USE IN VARIOUS TRANFORM PROCEDURES.                             
C                                                                               
      A2NNP1(0) = 0.0                                                           
      DO 70 N=1,KK                                                              
         ANNP1(N)  = A/REAL(N*(N+1))                                           
         A2NNP1(N) = REAL(N*(N+1))/A**2                                        
   70 CONTINUE                                                                  
C                                                                               
      DO 80 NL=1,NLAT                                                           
         RLAT = GLAT(NL)
         WTACSJ(NL) = 1.0/(A*COS(RLAT)**2)  
   80 CONTINUE                                                                  
C                                                                               
C     CALL ROUTINE TO FACTOR TRANSFORM LENGTH NLON AND CREATE TRIG TABLE        
C                                                                               
C     CALL FFTFAX (NLON,IFAX,TRIGS)                                             
C     REPLACE BY NEW INITIALIZATION:
      CALL SET99(TRIGS,IFAX,NLON)
      IF (IFAX(1) .EQ. 99) THEN                                                 
         WRITE (0,100)                                                          
  100    FORMAT (/,' STSWM: FATAL ERROR IN ROUTINE INPUT:',/,
     $    ' FAILURE OF ECMFFT LIBRARY ROUTINE SET99 (FFT SETUP)',/,
     $    ' PARAMETER NLON IN FILE ''PARAMS.i'' MUST HAVE ONLY',
     $    ' PRIME FACTORS 2,3 OR 5 !') 
         STOP                                                                
      ENDIF                                                                     
C                                                                               
C     DETERMINE LATITUDE AND LONGITUDE OF GAUSSIAN GRID,
C     RLOND AND RLATD IN DEGREES, FOR PLOTTING ROUTINES.           
C     EXTRA POINT RLOND(NLON+1) IS REQUIRED FOR PLOTS                           
C     ***** N.B. MUCH OF THE FOLLOWING CODE IS NOT RELEVENT *****               
C     *****      WITHOUT THE NCAR SYSTEM PLOT PACKAGE       *****               
C                                                                               
      DO 150 I=1,NLON                                                           
         RLOND(I) = GLON(I)*180.0/PI
  150 CONTINUE                                                                  
      RLOND(NLON+1) = RLOND(NLON) + (RLOND(NLON)-RLOND(NLON-1))                 
     $                              - 100.*EPS                                  
C                                                                               
      DO 160 J=1,NLAT                                                           
         RLATD(J) = GLAT(J)*180./PI                                 
  160 CONTINUE                                                                  
C
C     SPECIAL VALUE FEATURE OF CONTOUR PLOT PACKAGE
C
      SPVAL  = 1.0E36
      SPV(1) = 1.0E36
      SPV(2) = 1.0E36
C
C     SCALING PARAMETER FOR VECTOR PLOTS IN CPMVXY
C
      VVSF   = 1.0E-3
C                                                                               
C     INITIALIZE USER SELECTABLE CONTROL VARIABLES FOR PLOTTING ROUTINES        
C     SEE COMMON BLOCK /COMPLT/ FOR DEFINITIONS                                 
C                                                                               
      LCONT  = .TRUE.                                                          
      LPSP   = .FALSE.                                                          
      LOP    = .TRUE.                                                          
      LCP    = .FALSE.                                                          
      LG     = .TRUE.                                                          
      LU     = .FALSE.                                                          
      LV     = .FALSE.                                                          
      LZ     = .FALSE.                                                          
      LD     = .FALSE.                                                          
      LVV    = .TRUE.                                                          
      LVVG   = .FALSE.                                                          
      LGPHS  = .TRUE.                                                          
C                                                                               
C     CENTER LATITUDE, LONGITUDE AND ROTATION ANGLE FOR PLOTS
C     (DEFAULT VALUES IN DEGREES)
C                                                                               
      POLAT = 0.0
      POLNG = 0.0 
      POROT = 0.0
C                
C     DEFAULT FILENAMES FOR INPUT/OUTPUT OF SPECTRAL COEFFICIENTS
C
      FNIN  = 'VDGDATA.cdf'
      FNOUT = 'REFDATA.cdf'
C                                                                               
C     READ NAMELIST INFORMATION                                                 
C     (COMMON BUT NOT UNIVERSAL FORTRAN EXTENSION)                            
C                                                                               
      WRITE(6,200)
  200 FORMAT(/,' READING NAMELIST PARAMETERS FROM STANDARD INPUT:',/)
C
      IGROUP = 1                                                                
      READ (5,PHYVAR,ERR=1000)                                                  
C                                                                               
      IGROUP = 2                                                                
      READ (5,EXPDEF,ERR=1000)                                                  
C                                                                               
      IGROUP = 3                                                                
      READ (5,PLTDEF,ERR=1000)                                                  
C                                                                               
      IGROUP = 4
      READ (5,FNAMES,ERR=1000)
C
      GO TO 580                                                                 
C                                                                               
 1000 WRITE (0,1010) IGROUP                                                     
 1010 FORMAT(/, ' STSWM: FATAL ERROR IN SUBROUTINE INPUT ',/,          
     $      ' ERROR READING NAMELIST GROUP NUMBER ', I1)                 
      STOP                                                                     
C                                                                               
  580 CONTINUE
C
C     CHECK MAP PROJECTION FLAGS (ONLY ONE TRUE)
C
      IF ((LOP .AND. LPSP) .OR. (LOP .AND. LCP)
     $    .OR. (LCP .AND. LPSP)) THEN
         WRITE(0,650)
  650    FORMAT(/,' STSWM: ILLEGAL SPECIFICATION OF MAP PROJECTION',
     $          /,' IN NAMELIST INPUT FILE: AT MOST ONE OF LPSP,',
     $          /,' LOP, LCP MAY BE .TRUE.')
         STOP
      ENDIF
C
C     COMPUTE CORIOLIS PARAMETER FOR ROTATED COORDINATES
C
C     f = 2.0*OMEGA*(-COS(LAMBDA)*COS(THETA)*SIN(ALPHA)+
C                    SIN(THETA)*COS(ALPHA))
C     TRANSFORM INTO SPECTRAL SPACE:
C
C     WAVE M=0, N=1
      CORSC1 = SQRT(8.0/3.0)*OMEGA*COS(ALPHA)
C     WAVE M=1, N=1
      CORSC2 = - SQRT(4.0/3.0)*OMEGA*SIN(ALPHA)
C
      RETURN                                                                    
C                                                                               
      END                                                                       
