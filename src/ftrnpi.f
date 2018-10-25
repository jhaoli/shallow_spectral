      SUBROUTINE FTRNPI (DT, F3, F4, F7, PCOEF)
C                                                                              
C THIS SUBROUTINE PERFORMS A FORWARD TRANSFORM PROCEDURE USED           
C TO EVALUATE THE EXPLICIT PART OF THE RIGHT HAND SIDE 
C FOR THE GEOPOTENTIAL PROGNOSTIC EQUATION, USING SEMI-IMPLICIT
C TIMESTEPPING. (THE TERM Q OF EQ. (9) IN H. RITCHIE'S PAPER)
C THE COMPLEX COEFFICIENT VECTOR RETURNED BY THIS ROUTINE          
C IS NOT ZEROED AT THE START OF THE FORWARD TRANSFORM PROCEDURE             
C (I.E., THE USER MUST SPECIFY INITIAL STATE).  THE ROUTINE MAKES           
C USE OF INFORMATION CONTAINED IN COMMON BLOCK TRNSFM (GUASSIAN             
C WEIGHTS ASSOCIATED LEGENDRE POLYNOMIALS AND DERIVATIVES, ETC.)            
C                                                                              
C CALLED BY: FTRNIM
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
C TIME STEP
      REAL DT
C ARRAY OF FOURIER COEFFICIENTS (UCOS*PHI)^M
      COMPLEX F3(NFC,NLAT)
C ARRAY OF FOURIER COEFFICIENTS (VCOS*PHI)^M
      COMPLEX F4(NFC,NLAT)
C ARRAY OF FOURIER COEFFICIENTS (PHI^(TAU-1))^M
      COMPLEX F7(NFC,NLAT)
C
C     Output
C
C COMPUTED GEOPOTENTIAL NEW TIMESTEP
      COMPLEX PCOEF(NALP)                        
C                                                                               
C------ Local Variables ------------------------------------------------
C                                                                               
      COMPLEX CTMP1(MM+1), CTMP3(MM+1), CTMP4(MM+1), CTMP6(MM+1)
C
C     TEMPORARIES
C
      REAL FAC,FAD
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
         FAC = WTS(NL)*WTACSJ(NL)*DT
         FAD = WTS(NL)
C
         DO 105 M=0,MM
            JM = M+1                                                            
C     EVEN TERMS
            CTMP1(JM) = (F3(JM,NL) + F3(JM,SL))
     $                * CMPLX(0.0,-REAL(M)*FAC)
            CTMP1(JM) = CTMP1(JM) + (F7(JM,NL)+F7(JM,SL))*FAD
            CTMP3(JM) = (F4(JM,NL) - F4(JM,SL)) * FAC
C     ODD TERMS
            CTMP4(JM) = (F3(JM,NL) - F3(JM,SL))
     $                * CMPLX(0.0,-REAL(M)*FAC)      
            CTMP4(JM) = CTMP4(JM) + (F7(JM,NL)-F7(JM,SL))*FAD
            CTMP6(JM) = (F4(JM,NL) + F4(JM,SL)) * FAC
  105    CONTINUE                                                               
C                                                                               
         DO 120 JN=0,NN,2                                                       
            IS = LDIAG(JN,2)                                                    
C
            DO 110 JM=1,LDIAG(JN,1)                                             
               PCOEF(IS+JM) = PCOEF(IS+JM) 
     $           + ALP(IS+JM,NL) * CTMP1(JM)
     $           + DALP(IS+JM,NL) * CTMP3(JM)
  110       CONTINUE                                                            
  120    CONTINUE                                                               
C                                                                               
         DO 140 JN=1,NN,2                                                       
            IS = LDIAG(JN,2)                                                    
C
            DO 130 JM=1,LDIAG(JN,1)                                             
               PCOEF(IS+JM) = PCOEF(IS+JM) 
     $           + ALP(IS+JM,NL) * CTMP4(JM)
     $           + DALP(IS+JM,NL) * CTMP6(JM)
  130       CONTINUE                                                            
  140    CONTINUE                                                               
C                                                                               
  150 CONTINUE                                                                  
C                                                                               
      RETURN                                                                    
C                                                                               
      END                                                                       
