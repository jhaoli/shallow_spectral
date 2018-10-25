C***********************************************************************
C* INCLUDE FILE 'consts.i'                                             *
C***********************************************************************        
C*    Common Block /CONSTS/ contains constant information required by  *        
C*    the routines which integrate the spectral shallow water equations*        
C*    For initializtion of these variables see routine INPUT !         *
C*                                                                     *        
C*    A=RADIUS OF EARTH, OMEGA=ANGULAR VELOCITY, GRAV=GRAVITATIONAL    *
C*    ACCELERATION, DT=TIMESTEP,                                       *
C*    AFC=ASSELIN TIMEFILTER COEFFICIENT, EPS=MACHINE ACCURACY, SITS=  *
C*    .TRUE. FOR SEMI-IMPLICIT TIMESTEPPING, ELSE EXPLICIT TIMESTEPPING*
C*    NSTEP=TIMESTEP COUNTER, TAU=MODEL TIME IN HOURS, TAUO=FREQUENCY  *
C*    OF MODEL STATE OUTPUT TO FILE, TAUE=LENGTH OF MODEL RUN IN HOURS,*
C*    GPHFRQ=FREQUENCY OF PLOTS, HDC=HORIZONTAL DIFFUSION COEFFICEINT, *
C*    FORCED=.TRUE. FOR EXTERNAL FORCING (CASE 4), MOMENT=.TRUE. FOR   *
C*    EXTERNAL MOMENTUM FORCING, ELSE VORTICITY-DIVERGENCE FORCING,    *
C*    ALPHA= ROTATION                                                  *
C*    ANKLE IN RADIANS, EGYFRQ=FREQUENCY OF ANALYSIS FOR CONSERVATION  *
C*    PROPERTIES, ERRFRQ=FREQUENCY OF ANALYSIS FOR L-ERRORS,MIN,MAX,   *
C*    ICOND=NUMBER OF TEST CASE, LGPHS=.TRUE. FOR  PLOTTING OUTPUT,    *
C*    FTOPO=.TRUE. IF SURFACE TOPOGRAPHY (MOUNT .NE. 0),               *
C*    FNIN/FNOUT=FILENAMES FOR INPUT/OUTPUT OF SPECTRAL COEFFICIENTS,  *
C*    CHEXP=EXPERIMENT #,STRUNC=SPECTRAL TRUNCATION TYPE               *
C***********************************************************************
      REAL
     $      A, OMEGA, GRAV, DT, AFC, EPS, 
     $      TAU, TAUO, TAUE, GPHFRQ,
     $      HDC, ALPHA, EGYFRQ, ERRFRQ, SPCFRQ
      INTEGER
     $      NSTEP, ICOND
      LOGICAL
     $      SITS, FORCED, MOMENT, LGPHS, FTOPO
      CHARACTER*80 FNIN, FNOUT
      CHARACTER*4 CHEXP
      CHARACTER*6 STRUNC
      COMMON  / CONSTS / 
     $      A, OMEGA, GRAV, DT, AFC, EPS,  
     $      TAU, TAUO, TAUE, GPHFRQ, 
     $      HDC, ALPHA, EGYFRQ, ERRFRQ, SPCFRQ
      COMMON  / CONSTS / 
     $      NSTEP, ICOND
      COMMON  / CONSTS / 
     $      SITS, FORCED, MOMENT, LGPHS, FTOPO
      COMMON / FILES / FNIN, FNOUT, CHEXP, STRUNC
C***********************************************************************
C* END INCLUDE FILE                                                    *        
C***********************************************************************        
