      SUBROUTINE CEPS(CMN,DMN,EMN,EPSIL,LROW,LDIAG)
C                                                                              
C CALCULATES THE STRUCTURE OF THE RECURRENCE COEFFICIENT MATRICES           
C USED IN THE BELOUSOV PROCEDURE, STORES THIS INFORMATION IN                
C LROW(0:KK,2), AND EVALUATES THE COEFFICIENT MATRICES                    
C THE SPECTRAL TRUNCATION PARAMETERS ARE DEFINED BY 3 PARAMETERS:           
C MM, THE LARGEST FOURIER WAVENUMBER; KK, THE HIGHEST DEGREE OF THE         
C ASSOCIATED LEGENDRE POLYNOMIALS, AND NN THE HIGHEST DEGREE OF THE         
C ASSOCIATED LEGENDRE POLYNOMIALS FOR M=0.  THE LENGTH OF THE EPSIL         
C RECURRENCE ARRAY (AS WITH THE ASSOC. LEGENDRE POLYNOMIAL ARRAY)           
C IS GIVEN BY THE RELATION LEN = ((MM+1)*(NN+1) - (LMN**2 + LMN)/2)         
C WHERE LMN = MM + NN - KK.  VARIABLES ARE STORED ALONG DIAGONALS           
C STARTING WITH THE DIAGONAL DEFINED BY M=N.  THE LENGTH OF EACH            
C DIAGONAL IS STORED IN THE ARRAY LDIAG(0:NN,1) AND IS EVALUATED            
C USING THE EXPRESSION (MM+1) - AMAX(MM + N - KK, 0) WHERE 0>N>NN.          
C CUMULATIVE DIAGONAL LENGTHS (CUMULATIVE DISPLACEMENTS) ARE STORED         
C IN LDIAG(0:NN,2) SO THAT EPSIL OF ORDER M AND DEGREE N IS                 
C ADDRESSED AS EPSIL(1 + LDIAG(N-M,2)+M)).                                  
C
C CALLED BY: INPUT
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
C     Output
C
C RECURRENCE COEFFICIENTS FOR ALP
      REAL CMN(LRM+1)
C RECURRENCE COEFFICIENTS FOR ALP
      REAL DMN(LRM+1)
C RECURRENCE COEFFICIENTS FOR ALP
      REAL EMN(LRM+1)
C RECURRENCE COEFFICIENTS FOR DALP
      REAL EPSIL(NALP)
C INDEX FOR ROWS
      INTEGER LROW(0:KK,2)
C INDEX FOR DIAGONALS
      INTEGER LDIAG(0:NN,2)
C
C------ Local Variables ------------------------------------------------
C
      INTEGER K,M,N,JM,JN,IS
      INTEGER NSQD,MSQD
C
C----- Executable Statements -------------------------------------------
C
C     FIRST DETERMINE LDIAG(0:NN,1:2), DIAGONAL LENGTHS AND CUMULATIVE          
C     DISPLACEMENTS FOR ALP, DALP, AND EPSIL ARRAYS                             
C                                                                               
      LDIAG(0,1) = MM+1
      LDIAG(0,2) = 0                                                            
C                                                                               
      DO 10 N=1,NN                                                              
         LDIAG(N,1) = (MM+1) - MAX0(MM+N-KK, 0)                                
         LDIAG(N,2) = LDIAG(N-1,1) + LDIAG(N-1,2)                               
   10 CONTINUE                                                                  
C                                                                               
C     OTHER STRUCTURE INFORMATION (ROWS)                            
C     BEGIN BY CALCULATING ROW LENGTHS AND CUMULATIVE DISPLACEMENTS             
C     AND STORE THE INFORMATION IN LROW(0:KK,1:2)                               
C     ROW INFO IS FROM M=0 TO M=MM FOR EACH 0 <= N <= KK (I.E., IT DOES        
C     NOT EXCLUDE THE EMPTY AREA BETWEEN (NN+M) AND KK WHEN KK > NN)            
C                                                                               
      LROW(0,1) = 1                                                             
      LROW(0,2) = 0                                                             
C                                                                               
      DO 15 K=1,KK                                                              
         LROW(K,1) = MIN0(K+1,MM+1)                                             
         LROW(K,2) = LROW(K-1,1) + LROW(K-1,2)                                  
   15 CONTINUE                                                                  
C                                                                               
C     COMPUTE EPSIL MATRIX                                                      
C     BEGIN WITH FIRST TWO DIAGONALS                                            
C                                                                               
      DO 25 M=0,MM                                                              
         JM = M+1                                                               
         IS = LDIAG(1,2)                                                        
         EPSIL(JM)    = 0.0                                                     
         EPSIL(IS+JM) = 1.0/SQRT(REAL(2*M+3))                                  
   25 CONTINUE                                                                  
C                                                                               
C     REMAINING DIAGONALS (JN = 2 THROUGH NN)                                   
C                                                                               
      DO 40 JN = 2,NN                                                           
         IS   = LDIAG(JN,2)                                                     
         DO 30 M=0,LDIAG(JN,1)-1                                                
            JM   = M+1                                                          
            N    = JN+M                                                         
            NSQD = N*N                                                          
            MSQD = M*M                                                          
            EPSIL(IS+JM) = SQRT(REAL(NSQD-MSQD)/REAL(4*NSQD-1))               
   30    CONTINUE                                                               
   40 CONTINUE                                                                  
C                                                                               
C     COMPUTE RECURRENCE COEFFICIENT MATRICES CMN, DMN, EMN                     
C     REASON FOR STORING BY ROW IS TO AVOID BANK CONFLICTS WHEN                 
C     EVALUATING THE POLYNOMIALS IN SUBROUTINE CALP                             
C                                                                               
C     CHECK FOR SUFFICIENT STORAGE (PARAMETER LRM CORRECTLY SPECIFIED)          
C                                                                               
      IF (LRM .LT. ((MM+1)*(KK+1)-(MM*MM+MM)/2-1)) THEN                         
         WRITE (0,2) LRM                                                        
    2    FORMAT(/,' STSWM: FATAL ERROR IN SUBROUTINE CEPS ',/,        
     $     ' LRM INCORRECTLY SPECIFIED FOR BELOUSOV RECURRENCE',        
     $     ', LRM = ', I10)                                                
         STOP                                                                  
      ENDIF                                                                     
C                                                                               
      DO 130 N = 2,KK                                                           
         IS   = LROW(N,2)                                                       
         DO 120 M=2,LROW(N,1)-1                                                 
            JM = M+1                                                            
            CMN(IS+JM) = SQRT(REAL((2*N+1)*(M+N-1)*(M+N-3))                    
     $                      /(REAL((2*N-3)*(M+N)*(M+N-2))))                    
            DMN(IS+JM) = SQRT(REAL((2*N+1)*(M+N-1)*(N-M+1))                    
     $                      /(REAL((2*N-1)*(M+N)*(M+N-2))))                    
            EMN(IS+JM) = SQRT(REAL((2*N+1)*(N-M))                              
     $                      /(REAL((2*N-1)*(M+N))))                            
  120    CONTINUE                                                               
  130 CONTINUE                                                                  
C                                                                               
      RETURN                                                                    
      END                                                                       
