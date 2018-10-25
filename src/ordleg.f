      SUBROUTINE ORDLEG(COA, IR, SX)
C
C THIS ROUTINE IS USED TO CALCULATE ORDINARY LEGENDRE POLYNOMIALS
C
C CALLED BY: GLATS
C CALLS:
C
C REVISIONS:
C 7-10-92 CHANGE TO CCM CODING CONVENTIONS (R. JAKOB)
C
C------------ Arguments ------------------------------------------------
C
C     Input
C
C COSINE OF COLATITUDE
      REAL COA
C WAVE NUMBER
      INTEGER IR
C
C     Output
C
C LEGENDRE POLYNOMIAL EVALUATED AT COA
      REAL SX
C
C------ Local Variables ------------------------------------------------
C
      REAL SQR2,DELTA,THETA,C1,FN,FN2,FN2SQ,ANG,S1,C4,A,B,FK
      INTEGER IRPP,IRPPM,N,N1,KK,K
C
C----- Executable Statements -------------------------------------------
C
      SQR2  = SQRT(2.)                                                          
      IRPP  = IR + 1                                                            
      IRPPM = IRPP - 1                                                          
      DELTA = ACOS(COA)                                                         
C                                                                               
      THETA = DELTA                                                             
      C1    = SQR2                                                              
C                                                                               
      DO 20 N = 1, IRPPM                                                        
         FN    = REAL(N)                                                       
         FN2   = 2.*FN                                                          
         FN2SQ = FN2*FN2                                                        
         C1    = C1* SQRT(1.0-1.0/FN2SQ)                                        
   20 CONTINUE                                                                  
C                                                                               
      N   = IRPPM                                                               
      ANG = FN*THETA                                                            
      S1  = 0.0                                                                 
      C4  = 1.0                                                                 
      A   = -1.0                                                                
      B   = 0.0                                                                 
      N1  = N+1                                                                 
C                                                                               
      DO 30 KK = 1, N1, 2                                                       
         K = KK-1                                                               
         IF (K .EQ. N) C4 = 0.5*C4   
         S1  = S1+C4* COS(ANG)                                                  
         A   = A+2.0                                                            
         B   = B+1.0                                                            
         FK  = REAL(K)                                                         
         ANG = THETA*(FN-FK-2.0)                                                
         C4  = (A*(FN-B+1.0)/(B*(FN2-A)))*C4                                    
   30 CONTINUE                                                                  
C                                                                               
      SX = S1*C1                                                                
C                                                                               
      RETURN                                                                    
      END                                                                       
