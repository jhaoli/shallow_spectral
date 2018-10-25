      SUBROUTINE FTRNPE (DT, ALPHA, BETA, GAMMA, PCOEF)
C                                                                              
C THIS SUBROUTINE PERFORMS A FORWARD TRANSFORM PROCEDURE USED           
C TO EVALUATE THE RIGHT HAND SIDE OF THE SHALLOW WATER EQUATIONS            
C FOR THE GEOPOTENTIAL PROGNOSTIC EQUATION.
C THE COMPLEX COEFFICIENT VECTOR RETURNED BY THIS ROUTINE          
C IS NOT ZEROED AT THE START OF THE FORWARD TRANSFORM PROCEDURE             
C (I.E., THE USER MUST SPECIFY INITIAL STATE).  THE ROUTINE MAKES           
C USE OF INFORMATION CONTAINED IN COMMON BLOCK TRNSFM (GUASSIAN             
C WEIGHTS ASSOCIATED LEGENDRE POLYNOMIALS AND DERIVATIVES, ETC.)            
C                                                                              
C CALLED BY: ADVECT, FTRNEX
C CALLS:
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
C
C------------ Arguments ------------------------------------------------
C
C     Input
C
C TIMESTEP
      REAL DT
C ARRAY OF FOURIER COEFFICIENTS (UCOS*PHI)^M
      COMPLEX ALPHA(NFC,NLAT)
C ARRAY OF FOURIER COEFFICIENTS (VCOS*PHI)^M
      COMPLEX BETA(NFC,NLAT)
C ARRAY OF FOURIER COEFFICIENTS (PHI^TAU-1)^M 
      COMPLEX GAMMA(NFC,NLAT)
C
C     Output
C
C COMPUTED GEOPOTENTIAL NEW TIMELEVEL
      COMPLEX PCOEF(NALP)
C
C------ Local Variables ------------------------------------------------
C                                                                               
      COMPLEX CTMP1(MM+1), CTMP2(MM+1), CTMP3(MM+1), CTMP4(MM+1) 
C
C     TEMPORARIES
C
      REAL FAC, FAC2
      INTEGER M,JM,JN,NL,IS,SL
C
C----- Executable Statements -------------------------------------------
C
C     VARY M AND N SO THAT PROCEDURE MOVES ALONG DIAGONALS DENOTED              
C     BY INDEX JN.  M IS GIVEN BY (JM-1); N IS GIVEN BY (JN+M).                 
C     TAKE ADVANTAGE OF SYMMETRIC CHARACTER OF LEGENDRE POLYNOMIALS             
C     PROCEDURE ASSUMES EVEN NUMBER OF GAUSSIAN LATITUDES ...                   
C                                                                               
      DO 150 NL=1,NLAT/2                                                        
         SL = NLAT-NL+1
         FAC = DT*WTACSJ(NL)
         FAC2 = FAC*WTS(NL)
         DO 105 M=0,MM
            JM = M+1                                                            
            CTMP1(JM) = (ALPHA(JM,NL) + ALPHA(JM,SL))
     $                * CMPLX(0.0,-REAL(M)*FAC)
     $                + (GAMMA(JM,NL) + GAMMA(JM,SL))
            CTMP1(JM) = CTMP1(JM) * WTS(NL)
            CTMP2(JM) = (BETA(JM,NL) - BETA(JM,SL)) * FAC2 
C
            CTMP3(JM) = (ALPHA(JM,NL) - ALPHA(JM,SL))
     $                * CMPLX(0.0,-REAL(M)*FAC)
     $                + (GAMMA(JM,NL) - GAMMA(JM,SL))
            CTMP3(JM) = CTMP3(JM) * WTS(NL)
            CTMP4(JM) = (BETA(JM,NL) + BETA(JM,SL)) * FAC2  
  105    CONTINUE                                                               
C                                                                               
         DO 120 JN=0,NN,2                                                       
            IS = LDIAG(JN,2)                                                    
            DO 110 JM=1,LDIAG(JN,1)                                             
               PCOEF(IS+JM) = PCOEF(IS+JM) 
     $                      + ALP(IS+JM,NL) * CTMP1(JM)
     $                      + DALP(IS+JM,NL) * CTMP2(JM) 
  110       CONTINUE                                                            
  120    CONTINUE                                                               
C                                                                               
         DO 140 JN=1,NN,2                                                       
            IS = LDIAG(JN,2)                                                    
            DO 130 JM=1,LDIAG(JN,1)                                             
               PCOEF(IS+JM) = PCOEF(IS+JM) 
     $                      + ALP(IS+JM,NL) * CTMP3(JM) 
     $                      + DALP(IS+JM,NL) * CTMP4(JM) 
  130       CONTINUE                                                            
  140    CONTINUE                                                               
C                                                                               
  150 CONTINUE                                                                  
C                                                                               
      RETURN                                                                    
C                                                                               
      END                                                                       
