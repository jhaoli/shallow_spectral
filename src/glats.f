      SUBROUTINE GLATS(NLAT, EPS, F, WT)                                    
C                                                                              
C THIS ROUTINE CALCULATES THE ROOTS OF THE ORDINARY LEGENDRE                
C POLYNOMIALS OF ORDER NLAT WHICH CORRESPOND TO THE LATITUDE              
C POINTS USED IN A GAUSSIAN QUADRATURE PROCEDURE ON THE SPHERE.             
C THE ROUTINE USES A NEWTON-RAPHSON ITERATION PROCEDURE TO                  
C DETERMINE THE ROOTS TO A PRECISION GIVEN BY EPS.  A ROUTINE TO            
C CALCULATE THE ORDINARY LEGENDRE POLYNOMIALS, IN THIS CASE ORDLEG,         
C IS REQUIRED.                                                              
C                                                                              
C CALLED BY: INPUT
C CALLS: ORDLEG
C
C REVISIONS:
C 7-10-92 CHANGE TO CCM CODING CONVENTIONS (R. JAKOB)
C
C
C------------ Arguments ------------------------------------------------
C
C     Input
C
C NUMBER OF ZEROS BETWEEN POLES
      INTEGER NLAT
C CONVERGENCE LIMIT FOR THE ITERATION PROCEDURE
      REAL EPS
C
C     Output
C
C RESULTING ROOTS EXPRESSED IN LATITUDE
      REAL F(NLAT)
C CORRESPONDING GAUSSIAN WEIGHTS
      REAL WT(NLAT)
C
C------ Local Parameters -----------------------------------------------
C
      REAL PI
      PARAMETER (PI=3.141592653589793)
C
C------ Local Variables ------------------------------------------------
C
      REAL FI, FI1, FN, DOT, DN, DN1, A, B, G, GM, GP, GT,
     &     FTEMP, GTEMP
      INTEGER IR, NL, IRP, IRM, NITER
C
C----- Executable Statements -------------------------------------------
C
      IR  = NLAT                                                            
      FI  = REAL(IR)                                                           
      FI1 = FI+1.                                                               
      FN  = REAL(NLAT/2)                                                       
C                                                                               
C     FIRST GUESS AT ROOTS; F REPRESENTS COSINE OF COLATITUDE                   
C                                                                               
      DO 20 NL = 1, NLAT/2                                                      
         DOT   = REAL(NL-1)                                                    
         F(NL) = -PI*.5*(DOT+.5)/FN + PI*.5                                     
         F(NL) =  SIN(F(NL))                                                    
   20 CONTINUE                                                                  
C                                                                               
      DN  = FI/SQRT(4.*FI*FI-1.)                                                
      DN1 = FI1/SQRT(4.*FI1*FI1-1.)                                             
      A   = DN1*FI                                                              
      B   = DN*FI1                                                              
      IRP = IR + 1                                                              
      IRM = IR - 1                                                              
C                                                                               
      DO 50 NL = 1, NLAT/2                                                      
C                                                                               
C        NEWTON ITERATION LOOP                                                  
C                                                                               
         NITER = 0                                                              
C                                                                               
   30    NITER = NITER + 1                                                      
         IF(NITER .GE. 100) THEN                                                
            WRITE (0, 90) NL                                                    
   90       FORMAT(/,' STSWM: FATAL ERROR IN GLATS:',/,
     &        ' NO CONVERGENCE IN NEWTON ITERATION FOR ROOT',/, 
     &             ' LATITUDE  NL = ', I3)                                  
            STOP                                                                
         ENDIF                                                                  
C                                                                               
         CALL ORDLEG(F(NL), IR, G)                                            
         CALL ORDLEG(F(NL), IRM, GM)                                            
         CALL ORDLEG(F(NL), IRP, GP)                                            
         GT    = (A*GP-B*GM)/(F(NL)*F(NL)-1.)                                   
         FTEMP = F(NL) - G/GT                                                   
         GTEMP = F(NL) - FTEMP                                                  
         F(NL) = FTEMP                                                          
C
C        CONVERGENCE CRITERION
C
         IF (ABS(GTEMP) .GT. 2.0*EPS) GO TO 30  
C                                                                               
   50 CONTINUE                                                                  
C                                                                               
C     CALCULATE WEIGHTS AND EXPRESS ROOTS IN TERMS OF COLATITUDE                
C                                                                               
      DO 60 NL = 1, NLAT/2                                                      
         A      = 2.0*(1.0 - F(NL)**2)                                          
         CALL ORDLEG(F(NL), IRM, B)                                             
         B      = B*B*FI*FI                                                     
         WT(NL) = A*(FI - 0.5)/B                                                
         F (NL) = PI/2.0-ACOS(F(NL)) 
   60 CONTINUE                                                                  
C
C     FILL IN THE INFORMATION ON THE OTHER HALF OF THE SPHERE
C
      DO 70 NL = 1,NLAT/2
         F (NLAT-NL+1) = - F(NL)
         WT(NLAT-NL+1) = WT(NL)
   70 CONTINUE
C                                                                               
      RETURN                                                                    
      END                                                                       
