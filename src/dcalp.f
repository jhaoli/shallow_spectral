      SUBROUTINE DCALP(MMH,NNH,KKH,CMN,DMN,EMN,EPSIL,LROW,LDIAG,RLAT,
     $                 ALP,DALP)
C                                                                              
C CALCULATES THE ASSOCIATED LEGENDRE POLYNOMIALS USING A HIGHLY        
C STABLE FOUR TERM RECURRENCE RELATION GIVEN BY BELOUSOV (1962).            
C THE SPECTRAL TRUNCATION PARAMETERS ARE DEFINED BY 3 PARAMETERS:           
C MMH, THE LARGEST FOURIER WAVENUMBER; KKH, THE HIGHEST DEGREE OF THE 
C ASSOCIATED LEGENDRE POLYNOMIALS, AND NNH THE HIGHEST DEGREE OF THE 
C ASSOCIATED LEGENDRE POLYNOMIALS FOR M=0.  THE LENGTH OF THE               
C ASSOCIATED LEGENDRE POLYNOMIAL ARRAY ALP IS GIVEN BY THE RELATION         
C LEN = ((MMH+1)*(NNH+1) - (LMN**2 + LMN)/2)  WHERE LMN = MMH + NNH - KKH 
C VARIABLES ARE STORED ALONG DIAGONALS STARTING WITH DIAGONAL M=N.          
C THE LENGTH OF EACH ROW IS STORED IN THE ARRAY LROW(0:MAXH,2).
C THE LENGTH OF EACH DIAGONAL IS STORED IN THE ARRAY LDIAG(0:MAXH,1)
C AND IS EVALUATED IN SUBROUTINE CEPS AS (MMH+1)-AMAX(MMH+N-KKH,0)  
C WHERE 0>N>NNH.  CUMULATIVE DIAGONAL LENGTHS (CUMULATIVE  
C DISPLACEMENTS) ARE ALSO STORED IN LDIAG(0:MAXH,2) SO THAT THE  
C ASSOCIATED LEGENDRE POLYNOMIAL OF ORDER M, DEGREE N, AND ARGUMENT         
C SNJ IS ADDRESSED AS ALP(1 + LDIAG(N-M,2)+M)), OR USING THE         
C STATEMENT FUNCTION IDSP AS ALP(IDSP(M,N)).  THE SAME FORM              
C APPLIES FOR ADDRESSING THE DERIVATIVES (DALP) AND RECURRENCE              
C COEFFICIENTS DEFINED IN THE EPSIL MATRIX.                                 
C
C CALLED BY: EVAL
C CALLS: 
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
C TRUNCATION PARAMETERS
      INTEGER MMH,NNH,KKH
C RECURRENCE COEFFICIENTS FOR ALP
      REAL CMN(NALPH), DMN(NALPH), EMN(NALPH)
C RECURRENCE COEFFICIENTS FOR DALP
      REAL EPSIL(NALPH)
C ROW LENGTH AND INDEX
      INTEGER LROW(0:MAXH,2)
C DIAGONAL LENGTH AND INDEX
      INTEGER LDIAG(0:MAXH,2)
C LATITUDE FOR EVALUATION OF BASISFUNCTIONS
      REAL RLAT
C
C     Output
C
C ASSOC. LEGENDRE POLYNOMIALS
      REAL ALP(NALPH)
C DERIV. ASSOC. LEGENDRE POLYNOMIALS
      REAL DALP(NALPH)
C
C------ Local Parameters -----------------------------------------------
C
      REAL PI
      PARAMETER (PI=3.141592653589793)
C
C------ Local Variables ------------------------------------------------
C
      REAL CN2N1N(MAXH+1,3)
      REAL COL1(MAXH), COL2(MAXH), SNPSUM(MAXH), CSPSUM(MAXH), 
     $     SQNP(MAXH), AN(MAXH), CSFAC(MAXH), SNFAC(MAXH), 
     $     COSTBL(MAXH), SINTBL(MAXH)            
C
C     TEMPORARIES
C
      REAL RNORM,COS2P,SNJ,THTA
      INTEGER N,M,JK,JM,JN,MSTART,NSTART
      INTEGER MLIM,IS,ISM1
C
C----- Statement Function ----------------------------------------------
C
C     ADDRESS COMPUTATION OF ALP, DALP, EPSIL            
C     E.G., EPSIL(M,N) = EPSIL(IDSP(M,N))                                       
C                                                                               
      INTEGER MDUM,NDUM,IDSP,IDSPR
      IDSP(MDUM,NDUM) = 1 + LDIAG(NDUM-MDUM,2)+MDUM                             
C                                                                               
C     STATEMENT FUNCTION FOR ADDRESS COMPUTATION OF BELOUSOV RECURRENCE         
C     COEFFICIENT MATRICES.  E.G., CMN(M,N) = CMN(IDSPR(M,N))                   
C                                                                               
      IDSPR(MDUM,NDUM) = 1 + LROW(NDUM,2)+MDUM                                  
C                                                                               
C----- Executable Statements -------------------------------------------
C                                                                               
C     COMPUTE ASSOCIATED LEGENDRE POLYNOMIALS AND THEIR DERIVATIVES FOR         
C     THE TRUNCATED WAVENUMBER SPACE DEFINED IN LDIAG(0:MAXH,2) WITH            
C     ARGUMENTS GIVEN BY PARAMETER RLAT USING BELOUSOV'S ALGORITHM.         
C     THE PROCEDURE INVOLVES EXTRA WORK FOR ANY TRUNCATION OTHER THAN           
C     TRIANGULAR BECAUSE BELOUSOV'S RECURRENCE REQUIRES POLYNOMIAL              
C     INFORMATION IN PART OF THE REGION BETWEEN KKH AND NNH.  MORE  
C     ELABORATE BOOKKEEPING COULD REDUCE THIS EXTRA COMPUTATION, BUT BY         
C     LESS THAN A FACTOR OF TWO (PROBABLY LESS THAN 10% OF TOTAL WORK).         
C
C     LATITUDE = RLAT
C                                                                               
C     BEGIN BY COMPUTING 1ST TWO ELEMENTS IN EACH ROW (M=0:1;N=0:KKH)
C     EVALUATE THE SERIES EXPANSIONS (19) AND (21) IN BELOUSOV (1962)        
C     FINAL RESULTS ARE STORED IN WORK ARRAYS COL1 AND COL2                  
C                                                                               
      COS2P = COS(RLAT)
      SNJ   = SIN(RLAT)
      THTA  = PI/2.0 - RLAT
      CN2N1N(1,1)  = SQRT(0.50)                                              
      ALP(1)  = CN2N1N(1,1)                                               
      DALP(1) = 0.0                                                       
C                                                                               
C        INITIALIZE WORKING SPACE                                               
C                                                                               
      DO 205 N=1,KKH 
         SNPSUM(N) = 0.0                                                     
         CSPSUM(N) = 0.0                                                     
         SQNP(N)   = 1.0/SQRT(REAL(N*N + N))                                
         CSFAC(N)  = 1.0                                                     
         SNFAC(N)  = REAL(N)*SQNP(N)                                        
         COSTBL(N) = COS(REAL(N)*THTA)                                 
         SINTBL(N) = SIN(REAL(N)*THTA)                                 
  205 CONTINUE                                                               
C                                                                               
      AN(1) = SQRT(0.75)                                                     
      DO 210 N=2,KKH 
         AN(N) = AN(N-1)*SQRT(1.0-(1.0/REAL(4*N*N)))                        
  210 CONTINUE                                                               
C                                                                               
C     EACH INCREMENT IN JK EVALUATES AN ADDITIONAL TERM IN EXPANSIONS        
C                                                                               
      JK=1                                                                   
      DO 215 N=1,KKH 
         CSPSUM(N) = CSPSUM(N)+COSTBL(N-JK+1)*CSFAC(N)                       
         SNPSUM(N) = SNPSUM(N)+SINTBL(N-JK+1)*SNFAC(N)                       
  215 CONTINUE                                                               
C                                                                               
      DO 225 JK=3,KKH+1,2  
C                                                                               
         NSTART = MAX0(JK-1,1)                                               
         N = NSTART                                                          
         CSFAC(N)  = REAL(JK-2)/REAL(JK-1)*REAL(2*N-JK+3)                 
     $               /REAL(2*N-JK+2)*CSFAC(N)                               
         CSPSUM(N) = CSPSUM(N) + CSFAC(N)*0.50                               
C                                                                               
         DO 220 N=NSTART+1,KKH 
            CSFAC(N)  = REAL(JK-2)/REAL(JK-1)*REAL(2*N-JK+3)              
     $                  /REAL(2*N-JK+2)*CSFAC(N)                            
            SNFAC(N)  = CSFAC(N)*REAL(N-JK+1)*SQNP(N)                       
            CSPSUM(N) = CSPSUM(N)+COSTBL(N-JK+1)*CSFAC(N)                    
            SNPSUM(N) = SNPSUM(N)+SINTBL(N-JK+1)*SNFAC(N)                    
  220    CONTINUE                                                            
  225 CONTINUE                                                               
C                                                                               
      RNORM = 1.0/ALP(1)                                                  
      DO 230 N=1,KKH 
         COL1(N) = AN(N)*CSPSUM(N)*RNORM                                     
         COL2(N) = AN(N)*SNPSUM(N)*RNORM                                     
  230 CONTINUE                                                               
C                                                                               
      DO 260 N=1,KKH 
C                                                                            
C        EVALUATE REMAINING POLYNOMIALS BY SWEEPING THROUGH ROWS N=1:KKH 
C        FIRST TWO ELEMENTS OBTAINED FROM THE ABOVE SERIES EXPANSIONS           
C                                                                               
         CN2N1N(1,3) = COL1(N)                                               
         CN2N1N(2,3) = COL2(N)                                               
C                                                                               
C        NECESSARY DETOUR TO "PRIME THE PIPELINE" (FIRST PASS)               
C                                                                               
         IF (N .EQ. 1) THEN                                                     
            CN2N1N(1,2)   = CN2N1N(1,3)                                      
            CN2N1N(2,2)   = CN2N1N(2,3)                                      
            ALP(MMH+2)  = CN2N1N(1,2)                                      
            DALP(MMH+2) = SQRT(3.0)*CN2N1N(1,1)-SNJ*CN2N1N(1,2) 
            ALP(2)     = CN2N1N(2,2)                                      
            DALP(2)    = -SNJ*CN2N1N(2,2)                            
         ELSE                                                               
C                                                                               
C        EVALUATE THE REMAINDER OF THIS ROW (M = 2, 3, 4, ...)               
C        USING BELOUSOV'S RECURRENCE RELATION                                
C                                                                               
         MLIM = MIN0(MMH,N-1) 
         DO 235 M=2,MLIM                                                     
            JM = M+1                                                         
            CN2N1N(JM,3) = CMN(IDSPR(M,N))*CN2N1N(JM-2,1)                    
     $                     -SNJ*(DMN(IDSPR(M,N))                        
     $                     *CN2N1N(JM-2,2) - EMN(IDSPR(M,N))                 
     $                     *CN2N1N(JM,2))                                    
  235    CONTINUE                                                            
C                                                                               
C        PUT VALUES OF THE POLYNOMIALS CONTAINED IN CN2N1N(0:MLIM,3)         
C        INTO THE ASSOCIATED LEGENDRE POLYNOMIAL ARRAY ALP                   
C                                                                               
         MSTART = MAX0(N-NNH,0) 
         DO 240 M=MSTART,MLIM                                                
            JM = M+1                                                         
            ALP(IDSP(M,N)) = CN2N1N(JM,3)                                 
  240    CONTINUE                                                            
C                                                                               
C        SPECIAL EVALUATION REQUIRED FOR DIAGONAL ELEMENT M=N                
C                                                                               
         IF (N .LE. MMH) THEN
            CN2N1N(N+1,3)      = SQRT(1.0 + (1.0/REAL(2*N)))*COS2P          
     $                             *CN2N1N(N,2)                              
            ALP(IDSP(N,N))  = CN2N1N(N+1,3)                               
            DALP(IDSP(N,N)) = -REAL(N)*SNJ*CN2N1N(N+1,3)            
         ENDIF
C                                                                               
C        MAKE ROOM FOR NEW POLYNOMIAL EVALUATION IN CN2N1N(0:MLIM,3)         
C                                                                               
         DO 250 JM=1,N+1                                                     
            CN2N1N(JM,1) = CN2N1N(JM,2)                                      
            CN2N1N(JM,2) = CN2N1N(JM,3)                                      
  250    CONTINUE                                                            
         ENDIF
C                                                                               
  260 CONTINUE                                                               
C                                                                               
C     EFFICIENTLY EVALUATE DERIVATIVES SEPARATELY (ALONG DIAGONALS)          
C                                                                               
      DO 280 JN = 1,NNH 
         IS   = LDIAG(JN,2)                                                  
         ISM1 = LDIAG(JN-1,2)                                                
         DO 270 M=0,LDIAG(JN,1)-1                                            
            JM = M+1                                                         
            N  = JN+M                                                        
            DALP(IS+JM) = REAL(2*N+1)*EPSIL(IS+JM)*
     $                  ALP(ISM1+JM)-REAL(N)*SNJ*ALP(IS+JM)
  270    CONTINUE                                                            
  280 CONTINUE                                                               
C                                                                               
      RETURN                                                                    
      END                                                                       
