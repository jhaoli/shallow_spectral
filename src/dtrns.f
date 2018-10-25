      SUBROUTINE DTRNS (NNH,SCOEF,LDIAG,TRIGS,ALP,RLAT,PHI)                
C                                                                              
C THIS SUBROUTINE PERFORMS INVERSE SPHERICAL HARMONIC TRANSFORMS   
C INTO FOURIER SPACE USING THE ASSOCIATED LEGENDRE POLYNOMIALS.
C IT IS A REDUCED VERSION OF THE SHTRNS ROUTINE 
C                                                                              
C CALLED BY: EVAL
C CALLS: DFT991
C
C REVISIONS:
C 7-10-92 CHANGE TO CCM CODING CONVENTIONS (R. JAKOB)
C
C---- Model Parameters -------------------------------------------------
C
      INCLUDE 'params.i'
C
C
C------------ Arguments ------------------------------------------------
C
C     Input
C
C SPECTRAL TRUNCATION PARAMETER
      INTEGER NNH
C SPHERICAL HARMONIC COEFFICIENT ARRAY
      COMPLEX SCOEF(NALPH)
C INDEX INTO SPECTRAL ARRAY
      INTEGER LDIAG(0:MAXH,2)
C COEFFICIENTS FOR DISCRETE FOURIER TRANSFORM
      REAL TRIGS(MAXH+1,2)
C ASSOC. LEGENDRE POLYNOMIALS
      REAL ALP(NALPH)
C LATITUDE
      REAL RLAT
C
C     Output
C
C COMPUTED FUNCTION VALUE AT GRIDPOINT
      REAL PHI
C
C------ Local Variables ------------------------------------------------
C
      COMPLEX   CEVEN(MAXH+1), CODD(MAXH+1)            
C
C     TEMPORARIES
C
      INTEGER JM,JN,IS
C
C     FOURIER COEFFICIENTS
C
      COMPLEX X(MAXH+1)
C                                                                               
C----- Executable Statements -------------------------------------------
C                                                                               
C     DETERMINE FOURIER COEFFICIENTS BY INVERSE LEGENDRE TRANSFORM.          
C     VARY M AND N SO PROCEDURE MOVES ALONG DIAGONALS DENOTED BY             
C     INDEX JN.  M IS GIVEN BY (JM-1) WHILE N IS GIVEN BY (JN+M).            
C     FIRST ACCUMULATE EVEN WAVENUMBER CONTRIBUTION                          
C                                                                               
      DO 180 JN=0,NNH,2 
         IS = LDIAG(JN,2)                                                    
C                                                                               
C        THIS DETOUR HELPS AVOID THE NEED TO SEPARATELY ZERO CEVEN           
C                                                                               
         IF (JN .EQ. 0) THEN  
            DO 150 JM=1,LDIAG(0,1)                                           
               CEVEN(JM) = SCOEF(IS+JM)*ALP(IS+JM)               
  150       CONTINUE                                                         
C                                                                               
         ELSE                                                               
C                                                                               
            DO 170 JM=1,LDIAG(JN,1)
               CEVEN(JM) = CEVEN(JM) + SCOEF(IS+JM)*ALP(IS+JM) 
  170       CONTINUE 
C
         ENDIF
  180 CONTINUE                                                               
C                                                                               
C     NEXT ACCUMULATE ODD WAVENUMBER CONTRIBUTION                            
C                                                                               
      DO 215 JN=1,NNH,2 
         IS = LDIAG(JN,2)                                                    
C                                                                               
C        THIS DETOUR HELPS AVOID THE NEED TO SEPARATELY ZERO CODD            
C                                                                               
         IF (JN .EQ. 1) THEN
            DO 190 JM=1,LDIAG(1,1)                                           
                  CODD (JM) = SCOEF(IS+JM)*ALP(IS+JM)               
  190       CONTINUE                                                         
C                                                                               
C        ACCOUNT FOR THE FACT THAT THE FIRST ODD DIAGONAL MAY BE             
C        SHORTER THAN THE FIRST EVEN DIAGONAL (PART OF THE GAME              
C        TO AVOID EXPLICITLY ZEROING THE ENTIRE CODD ARRAY)                  
C                                                                               
            IF (LDIAG(1,1) .LT. LDIAG(0,1)) THEN                             
               DO 200 JM=LDIAG(1,1)+1, LDIAG(0,1)                            
                     CODD (JM) =  (0.0,0.0)                               
  200          CONTINUE                                                      
C                                                                               
            ENDIF                                                            
         ELSE
C                                                                               
            DO 210 JM=1,LDIAG(JN,1) 
               CODD (JM) = CODD (JM) + SCOEF(IS+JM)*ALP(IS+JM) 
  210       CONTINUE  
C
         ENDIF
  215 CONTINUE                                                               
C                                                                               
C     COMBINE CONTRIBUTIONS OF EVEN AND ODD WAVENUMBERS TO OBTAIN            
C     COMPLEX FOURIER COEFFICIENTS, FOLLOWED BY INVERSE FFT                  
C                                                                               
      DO 225 JM=1,LDIAG(0,1)                                              
         IF (RLAT .GE. 0.0) THEN
            X(JM) = CEVEN(JM) + CODD(JM) 
         ELSE
            X(JM) = CEVEN(JM) - CODD(JM)
         ENDIF
  225 CONTINUE                                                            
C                                                                               
C     INVERSE FOURIER TRANSFORMATION 
C                                                                               
      CALL DFT991(X,TRIGS,MAXH,LDIAG(0,1),PHI)
C
      RETURN                                                                 
C                                                                               
      END                                                                       
