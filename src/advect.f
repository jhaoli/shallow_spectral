      SUBROUTINE ADVECT(DTA,FORCED)
C                                                                              
C THIS IS THE MAIN COMPUTATIONAL PROCEDURE FOR THE ADVECTION 
C EQUATION THAT ADVANCES ONE TIME-
C STEP. THIS TIMESTEPPING IS DONE FOR THE PROGNOSTIC VARIABLE
C GEOPOTENTIAL USING A LEAPFROG TIMESTEPPING ALGORITHM.
C DIVERGENCE AND VORTICITY ARE JUST COPIED FOR SIMPLICITY.
C THE OLD TIMELEVEL IS LNM1, THE DERIVATIVE IS EVALUATED AT
C TIMELEVEL LN, AND THE COMPUTED VALUES GO INTO TIMELEVEL LNP1.
C THE NONLINEAR PRODUCTS ARE EVALUATED AT THE GRIDPOINTS, AND
C A FFT TRANSFORM (FFT991) IS USED FOR EACH LATITUDE. 
C SHTRNS INVERSE TRANSFORMS TO GRIDPOINT SPACE.
C                                                                              
C CALLED BY: COMP1
C CALLS: FFT991, SHTRNS
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
C WORKSPACE
      INCLUDE 'wrkspc.i'
C
C------------ Arguments ------------------------------------------------
C
C     Input
C
C TIMESTEP IN SECONDS
      REAL DTA
C FLAG FOR FORCING TERMS
      LOGICAL FORCED
C
C------ Local Variables ------------------------------------------------
C
C     EQUIVALENCE TO STORAGE POOL IN COMMON BLOCK WRKSPC WHERE POSSIBLE         
C                                                                               
      REAL TRMC(NLON+2,NLAT), TRMD(NLON+2,NLAT),
     $     TRMG(NLON+2,NLAT)
      COMPLEX TRMCFC(NFC,NLAT), TRMDFC(NFC,NLAT),
     $        TRMGFC(NFC,NLAT)
      EQUIVALENCE (ETAFCG(1,1), WS1(1,1)),
     $            (DIVFCG(1,1), WS2(1,1)),
     $            (TRMC(1,1), TRMCFC(1,1), WS3(1,1)),
     $            (TRMD(1,1), TRMDFC(1,1), WS4(1,1)),
     $            (PHIFCG(1,1), WS5(1,1)),
     $            (TRMG(1,1), TRMGFC(1,1), WS6(1,1))
C                                                                               
C     FORCING TERMS (PHYSICAL & SPECTRAL)
C                                                                               
      REAL ETAFCG(NLON+2,NLAT), DIVFCG(NLON+2,NLAT), 
     $     PHIFCG(NLON+2,NLAT)
C
C     TEMPORARIES
C
      INTEGER I,J,L
C
C----- Executable Statements -------------------------------------------
C                                                                               
C     ADVECTION ROUTINE                                                
C                                                                               
C     ZERO TERMS FOR LATER ACCUMULATIONS                         
C                                                                            
      DO 815 L=1,NALP                                                        
         PHISC(L) = (0.0,0.0)                                                 
  815 CONTINUE                                                               
C
C     COMPUTE THE GEOPOTENTIAL FORCING TERM
C
      IF (FORCED) THEN
         CALL FORCE(ETAFCG,DIVFCG,PHIFCG)
      ENDIF
C 
C     EVALUATE NON-LINEAR ADVECTION TERMS 
C     AND OLD TIMELEVEL OF DIVERGENCE & GEOPOTENTIAL
C                                                                               
      DO 20 J=1,NLAT                                                            
         DO 10 I=1,NLON                                                         
            TRMC(I,J) = UCOS(I,J)*PHI(I,J,LN)
            TRMD(I,J) = VCOS(I,J)*PHI(I,J,LN)
            TRMG(I,J) = PHI(I,J,LNM1) - DTA*PHIBAR*DIV(I,J,LN)
            IF (FORCED) THEN
               TRMG(I,J) = TRMG(I,J) + DTA*PHIFCG(I,J)
            ENDIF
   10    CONTINUE                                                               
   20 CONTINUE                                                                  
C                                                                               
C     FOURIER TRANSFORM NON-LINEAR TERMS IN PLACE                               
C                                                                               
      CALL FFT991(TRMC,FFTWS2,TRIGS,IFAX,1,NLON+2,NLON,NLAT,-1)
      CALL FFT991(TRMD,FFTWS2,TRIGS,IFAX,1,NLON+2,NLON,NLAT,-1)
      CALL FFT991(TRMG,FFTWS2,TRIGS,IFAX,1,NLON+2,NLON,NLAT,-1)
C
C     FWD LEGENDRE TRANSFORM NONLINEAR TERMS IN CONTINUITY EQUATION             
C     (BECOMES ADVECTION EQUATION FOR NONDIVERGENT FIELD)
C                                                                               
      CALL FTRNPE (DTA, TRMCFC, TRMDFC, TRMGFC, PHISC)
C
C     TRANSFORM ADVECTED GEOPOTENTIAL BACK TO GRID SPACE 
C
      CALL SHTRNS(LVLS,LNP1,+1,PHI,PHISC)
C
C     COPY VORTICITY & DIVERGENCE FIELD
C     (SO THERE AREN'T MORE CHANGES TO PROGRAM STRUCTURE)
C
      DO 710 J=1,NLAT
         DO 700 I=1,NLON
            DIV(I,J,LNP1) = DIV(I,J,LNM1)
            ZETA(I,J,LNP1) = ZETA(I,J,LNM1)
  700    CONTINUE
  710 CONTINUE
C
C     FINISHED FOR ADVECTION EQUATION
C
      RETURN
C
      END                                                                       
