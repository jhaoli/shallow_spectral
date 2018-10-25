      SUBROUTINE CALP(CMN,DMN,EMN,EPSIL,LROW,LDIAG,ALP,DALP)
C
C CALCULATES THE ASSOCIATED LEGENDRE POLYNOMIALS USING A HIGHLY        
C STABLE FOUR TERM RECURRENCE RELATION GIVEN BY BELOUSOV (1962).            
C THE SPECTRAL TRUNCATION PARAMETERS ARE DEFINED BY 3 PARAMETERS:           
C MM, THE LARGEST FOURIER WAVENUMBER; KK, THE HIGHEST DEGREE OF THE         
C ASSOCIATED LEGENDRE POLYNOMIALS, AND NN THE HIGHEST DEGREE OF THE         
C ASSOCIATED LEGENDRE POLYNOMIALS FOR M=0.  THE LENGTH OF THE               
C ASSOCIATED LEGENDRE POLYNOMIAL ARRAY ALP IS GIVEN BY THE RELATION         
C LEN = ((MM+1)*(NN+1) - (LMN**2 + LMN)/2)  WHERE LMN = MM + NN - KK        
C VARIABLES ARE STORED ALONG DIAGONALS STARTING WITH DIAGONAL M=N.          
C THE LENGTH OF EACH ROW IS STORED IN THE ARRAY LROW(0:KK,2).
C THE LENGTH OF EACH DIAGONAL IS STORED IN THE ARRAY LDIAG(0:NN,1)          
C AND IS EVALUATED IN SUBROUTINE CEPS AS (MM+1)-AMAX(MM+N-KK,0)             
C WHERE 0>N>NN.  CUMULATIVE DIAGONAL LENGTHS (CUMULATIVE                    
C DISPLACEMENTS) ARE ALSO STORED IN LDIAG(0:NN,2) SO THAT THE               
C ASSOCIATED LEGENDRE POLYNOMIAL OF ORDER M, DEGREE N, AND ARGUMENT         
C SNJ(NL) IS ADDRESSED AS ALP(1 + LDIAG(N-M,2)+M),NL), OR USING THE         
C STATEMENT FUNCTION IDSP AS ALP(IDSP(M,N),NL).  THE SAME FORM              
C APPLIES FOR ADDRESSING THE DERIVATIVES (DALP) AND RECURRENCE              
C COEFFICIENTS DEFINED IN THE EPSIL MATRIX.                                 
C   
C CALLED BY: INPUT
C CALLS: GLAT
C
C REVISIONS:
C 7-10-92 CHANGE TO CCM CODING CONVENTIONS (R. JAKOB)
C
C---- Model Parameters -------------------------------------------------
C
      INCLUDE 'params.i'
C
C------------ Arguments ------------------------------------------------
C
C     Input
C
C RECURRENCE COEFFICIENTS FOR ALP
      REAL CMN(LRM+1)
C RECURRENCE COEFFICIENTS FOR ALP
      REAL DMN(LRM+1)
C RECURRENCE COEFFICIENTS FOR ALP
      REAL EMN(LRM+1)
C RECURRENCE COEFFICIENTS FOR DALP
      REAL EPSIL(NALP)
C ROW LENGTH AND INDEX
      INTEGER LROW(0:KK,2)
C DIAGONAL LENGTH AND INDEX
      INTEGER LDIAG(0:NN,2)
C
C     Output
C
C ASSOC. LEGENDRE POLYNOMIALS
      REAL ALP(NALP,NLAT/2)
C DERIV. ASSOC. LEGENDRE POLYNOMIALS
      REAL DALP(NALP,NLAT/2)
C
C------ Local Parameters -----------------------------------------------
C
      REAL PI
      PARAMETER (PI=3.141592653589793)
C
C------ Local Variables ------------------------------------------------
C
      REAL CN2N1N(KK+1,3)
      REAL COL1(KK), COL2(KK), SNPSUM(KK), CSPSUM(KK), SQNP(KK),           
     $     AN(KK), CSFAC(KK), SNFAC(KK), COSTBL(KK), SINTBL(KK)            
C
C     TEMPORARIES
C
      REAL RNORM,COS2P,RLAT,SNJ,THTA
      INTEGER N,M,JK,JM,JN,MSTART,NSTART
      INTEGER MLIM,IS,ISM1,NL,LAT
C
C----- External Functions ----------------------------------------------
C
C     GAUSSIAN LATITUDES
C
      EXTERNAL GLAT
      REAL GLAT
C
C----- Statement Functions ---------------------------------------------
C
C     ADDRESS COMPUTATION OF ALP, DALP, EPSIL            
C     E.G., EPSIL(M,N) = EPSIL(IDSP(M,N))                                       
C                                                                               
      INTEGER MDUM, NDUM, IDSP, IDSPR
      IDSP(MDUM,NDUM) = 1 + LDIAG(NDUM-MDUM,2)+MDUM                             
C                                                                               
C     ADDRESS COMPUTATION OF BELOUSOV RECURRENCE         
C     COEFFICIENT MATRICES.  E.G., CMN(M,N) = CMN(IDSPR(M,N))                   
C                                                                               
      IDSPR(MDUM,NDUM) = 1 + LROW(NDUM,2)+MDUM                                  
C
C----- Executable Statements -------------------------------------------
C
C     COMPUTE ASSOCIATED LEGENDRE POLYNOMIALS AND THEIR DERIVATIVES FOR         
C     THE TRUNCATED WAVENUMBER SPACE DEFINED IN LDIAG(0:NN,1:2) WITH            
C     ARGUMENTS GIVEN BY GRID FUNCTION GLAT USING BELOUSOV'S ALGORITHM.         
C     THE PROCEDURE INVOLVES EXTRA WORK FOR ANY TRUNCATION OTHER THAN           
C     TRIANGULAR BECAUSE BELOUSOV'S RECURRENCE REQUIRES POLYNOMIAL              
C     INFORMATION IN PART OF THE REGION BETWEEN KK AND NN.  MORE                
C     ELABORATE BOOKKEEPING COULD REDUCE THIS EXTRA COMPUTATION, BUT BY         
C     LESS THAN A FACTOR OF TWO (PROBABLY LESS THAN 10% OF TOTAL WORK).         
C
C     PRECOMPUTE SOME FACTORS
C
      AN(1) = SQRT(0.75)
      DO 190 N=2,KK
         AN(N) = AN(N-1)*SQRT(1.0-(1.0/REAL(4*N*N)))
  190 CONTINUE
C                                                                               
C     BEGIN PROCEDURE ... OUTER LOOP OVER LATITUDE (ARGUMENT)                   
C                                                                               
      DO 290 NL=1,NLAT/2                                                        
C                                                                               
C     DETERMINE PROPER INDEX FOR LATITUDE DEPENDENT QUANTITIES                  
C     (PROVISION FOR DOING THIS COMPUTATION ON THE FLY, NOT IN USE)             
C                                                                               
         LAT = NL                                                               
         RLAT = GLAT(LAT)
C
C        LATITUDE = RLAT
C                                                                               
C        BEGIN BY COMPUTING 1ST TWO ELEMENTS IN EACH ROW (M=0:1;N=0:KK)         
C        EVALUATE THE SERIES EXPANSIONS (19) AND (21) IN BELOUSOV (1962)        
C        FINAL RESULTS ARE STORED IN WORK ARRAYS COL1 AND COL2                  
C                                                                               
C        COS2P = SQRT(CSJ(LAT))                                                 
         COS2P = COS(RLAT)
         SNJ   = SIN(RLAT)
         THTA  = PI/2.0 - RLAT
         CN2N1N(1,1)  = SQRT(0.50)                                              
         ALP(1,NL)  = CN2N1N(1,1)                                               
         DALP(1,NL) = 0.0                                                       
C                                                                               
C        INITIALIZE WORKING SPACE                                               
C                                                                               
         DO 205 N=1,KK                                                          
            SNPSUM(N) = 0.0                                                     
            CSPSUM(N) = 0.0                                                     
            SQNP(N)   = 1.0/SQRT(REAL(N*N + N))                                
            CSFAC(N)  = 1.0                                                     
            SNFAC(N)  = REAL(N)*SQNP(N)                                        
            COSTBL(N) = COS(REAL(N)*THTA)                                 
            SINTBL(N) = SIN(REAL(N)*THTA)                                 
  205    CONTINUE                                                               
C                                                                               
C        EACH INCREMENT IN JK EVALUATES AN ADDITIONAL TERM IN EXPANSIONS        
C                                                                               
         JK=1                                                                   
         DO 215 N=1,KK                                                          
            CSPSUM(N) = CSPSUM(N)+COSTBL(N-JK+1)*CSFAC(N)                       
            SNPSUM(N) = SNPSUM(N)+SINTBL(N-JK+1)*SNFAC(N)                       
  215    CONTINUE                                                               
C                                                                               
         DO 225 JK=3,KK+1,2                                                     
C                                                                               
            NSTART = MAX0(JK-1,1)                                               
            N = NSTART                                                          
            CSFAC(N)  = REAL(JK-2)/REAL(JK-1)*REAL(2*N-JK+3)                 
     $                  /REAL(2*N-JK+2)*CSFAC(N)                               
            CSPSUM(N) = CSPSUM(N) + CSFAC(N)*0.50                               
C                                                                               
            DO 220 N=NSTART+1,KK                                                
               CSFAC(N)  = REAL(JK-2)/REAL(JK-1)*REAL(2*N-JK+3)              
     $                     /REAL(2*N-JK+2)*CSFAC(N)                            
               SNFAC(N)  = CSFAC(N)*REAL(N-JK+1)*SQNP(N)                       
               CSPSUM(N) = CSPSUM(N)+COSTBL(N-JK+1)*CSFAC(N)                    
               SNPSUM(N) = SNPSUM(N)+SINTBL(N-JK+1)*SNFAC(N)                    
  220       CONTINUE                                                            
  225    CONTINUE                                                               
C                                                                               
         RNORM = 1.0/ALP(1,NL)                                                  
         DO 230 N=1,KK                                                          
            COL1(N) = AN(N)*CSPSUM(N)*RNORM                                     
            COL2(N) = AN(N)*SNPSUM(N)*RNORM                                     
  230    CONTINUE                                                               
C                                                                               
         DO 260 N=1,KK                                                          
C                                                                               
C        EVALUATE REMAINING POLYNOMIALS BY SWEEPING THROUGH ROWS N=1:KK         
C        FIRST TWO ELEMENTS OBTAINED FROM THE ABOVE SERIES EXPANSIONS           
C                                                                               
            CN2N1N(1,3) = COL1(N)                                               
            CN2N1N(2,3) = COL2(N)                                               
C                                                                               
C           NECESSARY DETOUR TO "PRIME THE PIPELINE" (FIRST PASS)               
C                                                                               
            IF (N .EQ. 1) THEN                                                     
               CN2N1N(1,2)   = CN2N1N(1,3)                                      
               CN2N1N(2,2)   = CN2N1N(2,3)                                      
               ALP(MM+2,NL)  = CN2N1N(1,2)                                      
               DALP(MM+2,NL) = SQRT(3.0)*CN2N1N(1,1)                            
     $                           - SNJ*CN2N1N(1,2)                         
               ALP(2,NL)     = CN2N1N(2,2)                                      
               DALP(2,NL)    = -SNJ*CN2N1N(2,2)                            
            ELSE                                                               
C                                                                               
C           EVALUATE THE REMAINDER OF THIS ROW (M = 2, 3, 4, ...)               
C           USING BELOUSOV'S RECURRENCE RELATION                                
C                                                                               
            MLIM = MIN0(MM,N-1)                                                 
            DO 235 M=2,MLIM                                                     
               JM = M+1                                                         
               CN2N1N(JM,3) = CMN(IDSPR(M,N))*CN2N1N(JM-2,1)                    
     $                        -SNJ*(DMN(IDSPR(M,N))                        
     $                        *CN2N1N(JM-2,2) - EMN(IDSPR(M,N))                 
     $                        *CN2N1N(JM,2))                                    
  235       CONTINUE                                                            
C                                                                               
C           PUT VALUES OF THE POLYNOMIALS CONTAINED IN CN2N1N(0:MLIM,3)         
C           INTO THE ASSOCIATED LEGENDRE POLYNOMIAL ARRAY ALP                   
C                                                                               
            MSTART = MAX0(N-NN,0)                                               
            DO 240 M=MSTART,MLIM                                                
               JM = M+1                                                         
               ALP(IDSP(M,N),NL) = CN2N1N(JM,3)                                 
  240       CONTINUE                                                            
C                                                                               
C           SPECIAL EVALUATION REQUIRED FOR DIAGONAL ELEMENT M=N                
C                                                                               
            IF (N .LE. MM) THEN
               CN2N1N(N+1,3)      = SQRT(1.0 + (1.0/REAL(2*N)))*COS2P          
     $                                *CN2N1N(N,2)                              
               ALP(IDSP(N,N),NL)  = CN2N1N(N+1,3)                               
               DALP(IDSP(N,N),NL) = -REAL(N)*SNJ*CN2N1N(N+1,3)            
            ENDIF
C                                                                               
C           MAKE ROOM FOR NEW POLYNOMIAL EVALUATION IN CN2N1N(0:MLIM,3)         
C                                                                               
            DO 250 JM=1,N+1                                                     
               CN2N1N(JM,1) = CN2N1N(JM,2)                                      
               CN2N1N(JM,2) = CN2N1N(JM,3)                                      
  250       CONTINUE                                                            
C
            ENDIF
C                                                                               
  260    CONTINUE                                                               
C                                                                               
C        EFFICIENTLY EVALUATE DERIVATIVES SEPARATELY (ALONG DIAGONALS)          
C                                                                               
         DO 280 JN = 1,NN                                                       
            IS   = LDIAG(JN,2)                                                  
            ISM1 = LDIAG(JN-1,2)                                                
            DO 270 M=0,LDIAG(JN,1)-1                                            
               JM = M+1                                                         
               N  = JN+M                                                        
               DALP(IS+JM,NL) = REAL(2*N+1)*EPSIL(IS+JM)*
     $                     ALP(ISM1+JM,NL)-REAL(N)*SNJ*ALP(IS+JM,NL)
  270       CONTINUE                                                            
  280    CONTINUE                                                               
C                                                                               
  290 CONTINUE                                                                  
C                                                                               
      RETURN                                                                    
      END                                                                       
