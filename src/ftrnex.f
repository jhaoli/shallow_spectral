      SUBROUTINE FTRNEX(DTA)
C                                                                              
C THIS IS THE MAIN COMPUTATIONAL PROCEDURE FOR ONE EXPLICIT TIME-
C STEP. THIS TIMESTEPPING IS DONE FOR THE PROGNOSTIC VARIABLES
C GEOPOTENTIAL, DIVERGENCE AND VORTICITY (SPECTRAL TRANSFORM
C ALGORITHM), USING A LEAPFROG TIMESTEPPING ALGORITHM.
C THE OLD TIMELEVEL IS LNM1, THE DERIVATIVE IS EVALUATED AT
C TIMELEVEL LN, AND THE COMPUTED VALUES GO INTO TIMELEVEL LNP1.
C THE FORCING TERMS ARE ADDED FOR TEST CASE 4, USING
C THE ROUTINE FORCE TO OBTAIN THE GRID FORCING.
C THE NONLINEAR PRODUCTS ARE EVALUATED AT THE GRIDPOINTS, AND
C A FFT TRANSFORM (FFT991) IS USED FOR EACH LATITUDE. 
C THE NEW TIMELEVEL SPECTRAL COEFFICIENTS ARE COMPUTED BY THE
C ROUTINES
C VORTICITY: FTRNVE 
C DIVERGENCE: FTRNDE 
C GEOPOTENTIAL: FTRNPE 
C                                                                              
C CALLED BY: COMP1
C CALLS: FFT991, FORCE, FTRNDE, FTRNPE, FTRNVE, GLAT
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
C TRANSFORM ARRAYS
      INCLUDE 'trnsfm.i'
C TIME DEPENDENT FIELDS
      INCLUDE 'tdvars.i'
C CONSTANTS & TIMESTEPS
      INCLUDE 'consts.i'
C WORKSPACE
      INCLUDE 'wrkspc.i'
C MOUNTAIN DATA FOR CASE 5
      INCLUDE 'finit.i'
C
C------------ Arguments ------------------------------------------------
C
C     Input
C
C TIMESTEP IN SECONDS
      REAL DTA
C
C------ Local Variables ------------------------------------------------
C
C     EQUIVALENCE TO STORAGE POOL IN COMMON BLOCK WRKSPC WHERE POSSIBLE         
C                                                                               
      REAL TRMA(NLON+2,NLAT), TRMB(NLON+2,NLAT),
     $     TRMC(NLON+2,NLAT), TRMD(NLON+2,NLAT),
     $     TRME(NLON+2,NLAT), TRMF(NLON+2,NLAT),
     $     TRMG(NLON+2,NLAT), TRMH(NLON+2,NLAT)
      COMPLEX TRMAFC(NFC,NLAT), TRMBFC(NFC,NLAT),               
     $        TRMCFC(NFC,NLAT), TRMDFC(NFC,NLAT),
     $        TRMEFC(NFC,NLAT), TRMFFC(NFC,NLAT),
     $        TRMGFC(NFC,NLAT), TRMHFC(NFC,NLAT)
      EQUIVALENCE (TRMA(1,1), TRMAFC(1,1), WS1(1,1)),
     $            (TRMB(1,1), TRMBFC(1,1), WS2(1,1)),
     $            (TRMC(1,1), TRMCFC(1,1), WS3(1,1)),
     $            (TRMD(1,1), TRMDFC(1,1), WS4(1,1)),
     $            (TRME(1,1), TRMEFC(1,1), WS5(1,1)),
     $            (TRMF(1,1), TRMFFC(1,1), WS6(1,1)),
     $            (TRMG(1,1), TRMGFC(1,1), WS7(1,1)),
     $            (TRMH(1,1), TRMHFC(1,1), WS8(1,1))
C                                                                               
C     FORCING TERMS (PHYSICAL & SPECTRAL)
C                                                                               
      REAL ETAFCG(NLON+2,NLAT), DIVFCG(NLON+2,NLAT), 
     $     PHIFCG(NLON+2,NLAT)
C
C     TEMPORARIES
C
      REAL FAC,RLAT
      INTEGER I,J,L
C
C----- External Functions ----------------------------------------------
C
C     GAUSSIAN LATITUDES
C
      REAL GLAT
      EXTERNAL GLAT
C
C----- Executable Statements -------------------------------------------
C                                                                               
C     MAIN COMPUTATIONAL ROUTINE                                                
C                                                                               
C     ZERO TERMS FOR LATER ACCUMULATIONS                         
C                                                                            
      DO 815 L=1,NALP                                                        
         ZETASC(L) = (0.0,0.0)                                                 
         DIVSC(L) = (0.0,0.0)                                                 
         PHISC(L) = (0.0,0.0)                                                 
  815 CONTINUE                                                               
C
C     COMPUTE THE MOMENTUM FORCING TERMS (TEST CASE 4)
C
      IF (FORCED) THEN
         CALL FORCE(ETAFCG,DIVFCG,PHIFCG)
      ENDIF
C 
C     EVALUATE NON-LINEAR ADVECTION TERMS 
C     AND OLD TIMELEVEL OF VORTICITY, DIVERGENCE, GEOPOTENTIAL
C                                                                               
      DO 20 J=1,NLAT                                                            
         RLAT = GLAT(J)
C
C        LATITUDE = RLAT
C
         FAC = 1.0/(2.0*COS(RLAT)**2)
         DO 10 I=1,NLON                                                         
            TRMA(I,J) = UCOS(I,J)*ZETA(I,J,LN)
            TRMB(I,J) = VCOS(I,J)*ZETA(I,J,LN)

            IF (FORCED .AND. MOMENT) THEN
               TRMA(I,J) = TRMA(I,J)-DIVFCG(I,J)
               TRMB(I,J) = TRMB(I,J)+ETAFCG(I,J)
            ENDIF
            TRMC(I,J) = UCOS(I,J)*PHI(I,J,LN)
            TRMD(I,J) = VCOS(I,J)*PHI(I,J,LN)
            TRME(I,J) = ZETA(I,J,LNM1)
            IF (FORCED .AND. (.NOT. MOMENT)) THEN
               TRME(I,J) = TRME(I,J)+DTA*ETAFCG(I,J)
            ENDIF
            TRMF(I,J) = DIV(I,J,LNM1)
            IF (FORCED .AND. (.NOT. MOMENT)) THEN
               TRMF(I,J) = TRMF(I,J)+DTA*DIVFCG(I,J)
            ENDIF
            TRMG(I,J) = (UCOS(I,J)**2 + VCOS(I,J)**2)*FAC +
     $                  PHI(I,J,LN)
            IF (FTOPO) THEN
               TRMG(I,J) = TRMG(I,J) + GRAV*MOUNT(I,J)
            ENDIF
            TRMH(I,J) = PHI(I,J,LNM1) - DTA*PHIBAR*DIV(I,J,LN)
            IF (FORCED) THEN
               TRMH(I,J) = TRMH(I,J) + DTA*PHIFCG(I,J)
            ENDIF
   10    CONTINUE                                                               
   20 CONTINUE                                                                  
C                                                                               
C     FOURIER TRANSFORM NON-LINEAR TERMS IN PLACE                               
C                                                                               
      CALL FFT991(TRMA,FFTWS2,TRIGS,IFAX,1,NLON+2,NLON,NLAT,-1)
      CALL FFT991(TRMB,FFTWS2,TRIGS,IFAX,1,NLON+2,NLON,NLAT,-1)
      CALL FFT991(TRMC,FFTWS2,TRIGS,IFAX,1,NLON+2,NLON,NLAT,-1)
      CALL FFT991(TRMD,FFTWS2,TRIGS,IFAX,1,NLON+2,NLON,NLAT,-1)
      CALL FFT991(TRME,FFTWS2,TRIGS,IFAX,1,NLON+2,NLON,NLAT,-1)
      CALL FFT991(TRMF,FFTWS2,TRIGS,IFAX,1,NLON+2,NLON,NLAT,-1)
      CALL FFT991(TRMG,FFTWS2,TRIGS,IFAX,1,NLON+2,NLON,NLAT,-1)
      CALL FFT991(TRMH,FFTWS2,TRIGS,IFAX,1,NLON+2,NLON,NLAT,-1)
C
C     FWD LEGENDRE TRANSFORM & UPDATE OF VORTICITY SPECTRAL COEFFICIENTS
C
      CALL FTRNVE (DTA, TRMAFC, TRMBFC, TRMEFC, ZETASC) 
C
C     FWD LEGENDRE TRANSFORM NONLINEAR TERMS AND TIMESTEPPING
C     IN DIVERGENCE EQUATION
C
      CALL FTRNDE (DTA, TRMAFC, TRMBFC, TRMFFC, TRMGFC, DIVSC)        
C                                                                               
C     FWD LEGENDRE TRANSFORM NONLINEAR TERMS IN CONTINUITY EQUATION             
C     (BECOMES ADVECTION EQUATION FOR NONDIVERGENT FIELD)
C                                                                               
      CALL FTRNPE (DTA, TRMCFC, TRMDFC, TRMHFC, PHISC)
C
      RETURN                                                                    
C
      END                                                                       
