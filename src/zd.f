      SUBROUTINE ZD(UIC,VIC,DIVSC,ZETASC)
C                                                                              
C THIS SUBROUTINE OBTAINS VORTICITY AND DIVERGENCE SPECTRAL 
C COEFFICIENTS FROM THE GRIDPOINT VELOCITIES (UNSCALED BY COS(THETA)) 
C THE SUBROUTINE MAKES USE OF INFORMATION           
C CONTAINED IN COMMON BLOCKS (POLYNOMIALS, COEFFICIENTS, ETC.)              
C
C CALLED BY: INIT, STSWM
C CALLS: FFT991, GLAT
C
C---- Model Parameters -------------------------------------------------
C
      INCLUDE 'params.i'
C
C---- Common Blocks ----------------------------------------------------
C
C TRANSFORM ARRAYS
      INCLUDE 'trnsfm.i'
C WORKSPACE
      INCLUDE 'wrkspc.i'
C
C------------ Arguments ------------------------------------------------
C
C     Input
C
C EASTWARD WIND FIELD
      REAL UIC(NLON+2,NLAT)
C NORTHWARD WIND FIELD
      REAL VIC(NLON+2,NLAT)
C
C     Output
C
C COMPUTED DIVERGENCE COEFFICIENTS
      COMPLEX DIVSC(NALP)
C COMPUTED VORTICITY COEFFICIENTS
      COMPLEX ZETASC(NALP)
C
C------ Local Variables ------------------------------------------------
C                                                                               
C     EQUIVALENCE TO WORKSPACE STORAGE POOL /WRKSPC/                            
C                                                                               
      REAL UTMP(NLON+2,NLAT),VTMP(NLON+2,NLAT)
      COMPLEX UTMPFC(NFC,NLAT),VTMPFC(NFC,NLAT)
      EQUIVALENCE (UTMP(1,1),UTMPFC(1,1),WS1(1,1)),
     $            (VTMP(1,1),VTMPFC(1,1),WS2(1,1))
C
      REAL STMP1(NLON+2), STMP2(NLON+2)                                    
      EQUIVALENCE (STMP1(1),RVEC1(1)), (STMP2(1),RVEC2(1))                      
C
      COMPLEX CTMP1(MM+1), CTMP2(MM+1), CTMP3(MM+1), CTMP4(MM+1)              
C
C     TEMPORARIES
C
      REAL RLAT, FAC
      INTEGER I,J,JM,JN,M,NL,IS,INDEX,SL
C
C----- Statement Function ----------------------------------------------
C
C     GRID LATITUDE
C
      EXTERNAL GLAT
      REAL GLAT
C                                                                               
C----- Executable Statements -------------------------------------------
C
C     COMPUTE VORTICITY AND DIVERGENCE SPECTRAL COEFFICIENTS                    
C     FROM U AND V MOMENTUM COEFFICIENTS                                        
C                                                                               
C     FOURIER TRANSFORM U*COS(PHI) AND V*COS(PHI) FIELDS                
C                                                                               
      DO 20 J=1,NLAT
         RLAT = GLAT(J)
C
C        LATITUDE = RLAT = GLAT(J)
C
         DO 10 I=1,NLON
            UTMP(I,J)=UIC(I,J)*COS(RLAT)
            VTMP(I,J)=VIC(I,J)*COS(RLAT)
   10    CONTINUE
   20 CONTINUE
C
      CALL FFT991(UTMP,FFTWS2,TRIGS,IFAX,1,NLON+2,NLON,NLAT,-1)                 
      CALL FFT991(VTMP,FFTWS2,TRIGS,IFAX,1,NLON+2,NLON,NLAT,-1)                 
C                                                                               
C     COMPUTE VORTICITY AND DIVERGENCE COEFFICIENTS FROM U AND V FOURIER        
C     COEFFICIENTS. VARY M AND N SO THAT PROCEDURE MOVES ALONG DIAGONALS        
C     DENOTED BY INDEX JN.  M IS GIVEN BY (JM-1); N IS GIVEN BY (JN+M).         
C     TAKE ADVANTAGE OF SYMMETRIC CHARACTER OF LEGENDRE POLYNOMIALS             
C                                                                               
C     PROCEDURE FOR EVEN NUMBER OF GAUSSIAN LATITUDES ...                       
C                                                                               
      DO 210 NL=1,NLAT/2                                                        
         SL = NLAT-NL+1
         FAC = WTS(NL)*WTACSJ(NL)
C                                                                               
C        CALCULATE INTERMEDIATE QUANTITIES FOR COMPUTATIONL EFFICIENCY          
C                                                                               
C
         DO 100 M=0,MM
            JM=M+1                                                              
            STMP1(JM) = REAL(M)*FAC                                            
C                                                                               
            CTMP1(JM) = (UTMPFC(JM,NL) + UTMPFC(JM,SL)) * 
     $                  CMPLX(0.0,STMP1(JM))        
            CTMP2(JM) = (VTMPFC(JM,NL) - VTMPFC(JM,SL))*FAC              
C                                                                               
            CTMP3(JM) = (VTMPFC(JM,NL) + VTMPFC(JM,SL)) * 
     $                  CMPLX(0.0,STMP1(JM))        
            CTMP4(JM) = (UTMPFC(JM,NL) - UTMPFC(JM,SL))*FAC              
  100    CONTINUE                                                               
C                                                                               
C        SPECIAL HANDLING OF THE FIRST LATITUDE IN THE ACCUMULATIONS            
C        DETOUR ELIMINATES THE NEED TO SEPARATELY ZERO DIVSC AND ZETASC         
C                                                                               
         IF (NL .EQ. 1) THEN  
C                                                                               
C        EVEN DIAGONALS FIRST LATITUDE                                          
C                                                                               
            DO 120 JN=0,NN,2                                                    
               IS = LDIAG(JN,2)                                                 
C
               DO 110 M=0,LDIAG(JN,1)-1                                         
                  JM = M + 1                                                    
                  DIVSC(IS+JM) = CTMP1(JM)*ALP(IS+JM,NL)                     
     $                              - CTMP2(JM)*DALP(IS+JM,NL)                  
C                                                                               
                  ZETASC(IS+JM) = CTMP3(JM)*ALP(IS+JM,NL)                    
     $                               + CTMP4(JM)*DALP(IS+JM,NL)                 
  110          CONTINUE                                                         
  120       CONTINUE                                                            
C
         ELSE                                                                  
C                                                                               
C        EVEN DIAGONALS OTHER LATITUDES                                         
C                                                                               
            DO 140 JN=0,NN,2  
               IS = LDIAG(JN,2)   
C
               DO 130 M=0,LDIAG(JN,1)-1  
                  JM = M + 1 
                  DIVSC(IS+JM) = DIVSC(IS+JM)                                
     $                            + CTMP1(JM)*ALP(IS+JM,NL)                     
     $                            - CTMP2(JM)*DALP(IS+JM,NL)                    
C                                                                               
                  ZETASC(IS+JM) = ZETASC(IS+JM)                              
     $                             + CTMP3(JM)*ALP(IS+JM,NL)                    
     $                             + CTMP4(JM)*DALP(IS+JM,NL)                   
  130          CONTINUE  
  140       CONTINUE  
C
         ENDIF
C                                                                               
C        CALCULATE INTERMEDIATE QUANTITIES FOR COMPUTATIONAL EFFICIENCY
C                                                                               
         DO 160 M=0,MM
            JM=M+1                                                              
            STMP1(JM) = REAL(M)*FAC                                            
C                                                                               
            CTMP1(JM) = (UTMPFC(JM,NL) - UTMPFC(JM,SL)) * 
     $                  CMPLX(0.0,STMP1(JM))        
            CTMP2(JM) = (VTMPFC(JM,NL) + VTMPFC(JM,SL))*FAC              
C                                                                               
            CTMP3(JM) = (VTMPFC(JM,NL) - VTMPFC(JM,SL)) * 
     $                  CMPLX(0.0,STMP1(JM))        
            CTMP4(JM) = (UTMPFC(JM,NL) + UTMPFC(JM,SL))*FAC              
  160    CONTINUE                                                               
C                                                                               
C        SPECIAL HANDLING OF THE FIRST LATITUDE IN THE ACCUMULATIONS            
C        DETOUR ELIMINATES THE NEED TO SEPARATELY ZERO DIVSC AND ZETASC         
C                                                                               
         IF (NL .EQ. 1) THEN 
C                                                                               
C        ODD DIAGONALS FIRST LATITUDE                                           
C                                                                               
            DO 180 JN=1,NN,2                                                    
               IS = LDIAG(JN,2)                                                 
C
               DO 170 M=0,LDIAG(JN,1)-1                                         
                  JM = M + 1                                                    
                  DIVSC(IS+JM) = CTMP1(JM)*ALP(IS+JM,NL)                     
     $                              - CTMP2(JM)*DALP(IS+JM,NL)                  
C                                                                               
                  ZETASC(IS+JM) = CTMP3(JM)*ALP(IS+JM,NL)                    
     $                               + CTMP4(JM)*DALP(IS+JM,NL)                 
  170          CONTINUE                                                         
  180       CONTINUE                                                            
C
         ELSE                                                                  
C                                                                               
C        ODD DIAGONALS REMAING LATITUDES                                        
C                                                                               
            DO 200 JN=1,NN,2 
               IS = LDIAG(JN,2)  
C
               DO 190 M=0,LDIAG(JN,1)-1   
                  JM = M + 1 
                  DIVSC(IS+JM) = DIVSC(IS+JM)                                
     $                            + CTMP1(JM)*ALP(IS+JM,NL)                     
     $                            - CTMP2(JM)*DALP(IS+JM,NL)                    
C                                                                               
                  ZETASC(IS+JM) = ZETASC(IS+JM)                              
     $                             + CTMP3(JM)*ALP(IS+JM,NL)                    
     $                             + CTMP4(JM)*DALP(IS+JM,NL)                   
  190          CONTINUE   
  200       CONTINUE 
C
         ENDIF
C
  210 CONTINUE                                                                  
C                                                                               
C     ADD CORIOLIS TERM FOR ROTATED COORDINATE SYSTEM
C
C     WAVE M=0, N=1
      INDEX = LDIAG(1,2)+1
      ZETASC(INDEX) = ZETASC(INDEX) + CORSC1
C     WAVE M=1, N=1
      INDEX = LDIAG(0,2)+2
      ZETASC(INDEX) = ZETASC(INDEX) + CORSC2
C                                                                               
      RETURN                                                                    
      END
