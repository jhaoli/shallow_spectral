      SUBROUTINE FTRNDI (DT, F1, F2, F5, F6, DCOEF)
C                                                                              
C THIS ROUTINE PERFORMS A FORWARD TRANSFORM PROCEDURE USED           
C TO EVALUATE THE EXPLICIT PART OF THE RIGHT HAND SIDE
C FOR THE DIVERGENCE PROGNOSTIC EQUATION, USING SEMI-IMPLICIT
C TIMESTEPPING. (THE TERM M OF EQ. (8) IN RITCHIE'S PAPER)
C THE COMPLEX COEFFICIENT VECTOR RETURNED BY THIS ROUTINE          
C IS NOT ZEROED AT THE START OF THE FORWARD TRANSFORM PROCEDURE             
C (I.E., THE USER MUST SPECIFY INITIAL STATE).  THE ROUTINE MAKES           
C USE OF INFORMATION CONTAINED IN COMMON BLOCK TRNSFM (GAUSSIAN             
C WEIGHTS, ASSOCIATED LEGENDRE POLYNOMIALS AND DERIVATIVES, ETC.) 
C
C CALLED BY: FTRNIM
C CALLS:
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
C
C------------ Arguments ------------------------------------------------
C
C     Input
C
C TIME STEP
      REAL DT
C ARRAY OF FOURIER COEFFICIENTS (UCOS*ZETA)^M
      COMPLEX F1(NFC,NLAT)
C ARRAY OF FOURIER COEFFICIENTS (VCOS*ZETA)^M
      COMPLEX F2(NFC,NLAT)
C ARRAY OF FOURIER COEFFICIENTS (UCOS**2+VCOS**2)/FAC
      COMPLEX F5(NFC,NLAT)
C ARRAY OF FOURIER COEFFICIENTS (DIV^(TAU-1))^M
      COMPLEX F6(NFC,NLAT)
C
C     Ouput
C
C COMPUTED DIVERGENCE AT NEW TIMELEVEL
      COMPLEX DCOEF(NALP)
C
C------ Local Variables ------------------------------------------------
C                                                                               
      COMPLEX CTMP1(MM+1), CTMP2(MM+1), CTMP3(MM+1), 
     $        CTMP5(MM+1), CTMP6(MM+1), CTMP7(MM+1)
C
C     TEMPORARIES
C
      REAL FAC,FAD
      INTEGER M,N,JM,JN,NL,IS,SL
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
         FAD = WTS(NL)*DT
C
         DO 105 M=0,MM
            JM = M+1                                                            
C EVEN TERMS
            CTMP1(JM) = (F2(JM,NL) + F2(JM,SL))
     $                * CMPLX(0.0,REAL(M)*FAC)
            CTMP1(JM) = CTMP1(JM) + (F6(JM,NL)+F6(JM,SL))*WTS(NL)
            CTMP2(JM) = (F5(JM,NL) + F5(JM,SL))*FAD
            CTMP3(JM) = (F1(JM,NL) - F1(JM,SL))*FAC
C ODD TERMS
            CTMP5(JM) = (F2(JM,NL) - F2(JM,SL))
     $                * CMPLX(0.0,REAL(M)*FAC)
            CTMP5(JM) = CTMP5(JM) + (F6(JM,NL)-F6(JM,SL))*WTS(NL)
            CTMP6(JM) = (F5(JM,NL) - F5(JM,SL))*FAD
            CTMP7(JM) = (F1(JM,NL) + F1(JM,SL))*FAC
  105    CONTINUE                                                               
C                                                                               
         DO 120 JN=0,NN,2                                                       
            IS = LDIAG(JN,2)                                                    
C
            DO 110 JM=1,LDIAG(JN,1)                                             
               N = JM + JN - 1
               DCOEF(IS+JM) = DCOEF(IS+JM)  
     $          + ALP(IS+JM,NL)*(CTMP1(JM) + CTMP2(JM)*A2NNP1(N))
     $          + DALP(IS+JM,NL)*CTMP3(JM) 
  110       CONTINUE                                                            
  120    CONTINUE                                                               
C                                                                               
         DO 140 JN=1,NN,2                                                       
            IS = LDIAG(JN,2)                                                    
C
            DO 130 JM=1,LDIAG(JN,1)                                             
               N = JM + JN - 1
               DCOEF(IS+JM) = DCOEF(IS+JM) 
     $         + ALP(IS+JM,NL)*(CTMP5(JM) + CTMP6(JM)*A2NNP1(N))
     $         + DALP(IS+JM,NL)*CTMP7(JM) 
  130       CONTINUE                                                            
  140    CONTINUE                                                               
C                                                                               
  150 CONTINUE                                                                  
C                                                                               
      RETURN    
C                                                                               
      END                                                                       
