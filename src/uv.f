      SUBROUTINE UV(DIVSC,ZETASC,UCOS,VCOS)  
C                                                                              
C THIS SUBROUTINE OBTAINS CAPITAL U AND V MOMENTUM COMPONENTS FROM THE  
C VORTICITY AND DIVERGENCE SPECTRAL COEFFICIENTS
C THE SUBROUTINE MAKES EXTENSIVE USE OF INFORMATION           
C CONTAINED IN COMMON BLOCKS (POLYNOMIALS, COEFFICIENTS, ETC.)              
C                                                                              
C CALLED BY: COMP1
C CALLS: FFT991
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
C WORKSPACE
      INCLUDE 'wrkspc.i'
C
C------------ Arguments ------------------------------------------------
C
C     Input
C
C DIVERGENCE SPECTRAL COEFFICIENTS
      COMPLEX DIVSC(NALP)
C VORTICITY SPECTRAL COEFFICIENTS
      COMPLEX ZETASC(NALP)
C
C     Output
C
C COMPUTED EASTWARD WIND FIELD
      REAL UCOS(NLON+2,NLAT)
C COMPUTED NORTHWARD WIND FIELD
      REAL VCOS(NLON+2,NLAT)
C
C------ Local Variables ------------------------------------------------
C
C     EQUIVALENCE TO WORKSPACE STORAGE POOL /WRKSPC/                            
C                                                                               
      COMPLEX CMAT1(NALP), CMAT2(NALP), CMAT3(NALP), CMAT4(NALP)              
      EQUIVALENCE (CMAT1(1),RMAT1(1)), (CMAT2(1),RMAT2(1)),                     
     $            (CMAT3(1),RMAT3(1)), (CMAT4(1),RMAT4(1))                      
      COMPLEX UTMP1(NLAT/2,NFC), UTMP2(NLAT/2,NFC),             
     $        VTMP1(NLAT/2,NFC), VTMP2(NLAT/2,NFC)              
      EQUIVALENCE (UTMP1(1,1),WS1(1,1)), (UTMP2(1,1),WS1(1,NLAT/2+1)),          
     $            (VTMP1(1,1),WS2(1,1)), (VTMP2(1,1),WS2(1,NLAT/2+1))
C
C     TEMPORARIES
C
      INTEGER M,N,JM,JN,NL,IS,INDEX,SL
C
C----- Executable Statements -------------------------------------------
C
C     COMPUTE U AND V MOMENTUM COMPONENTS FROM VORTICITY AND DIVERGENCE         
C                                                                               
C                                                                               
C     COMPUTE INTERMEDIATE QUANTITIES FOR COMPUTATIONAL EFFICIENCY           
C     INCLUDE CORIOLIS TERM HERE
C                                                                               
      DO 230 JN=0,NN                                                         
         IS = LDIAG(JN,2)                                                    
C
         DO 220 JM=1,LDIAG(JN,1)                                             
            M = JM-1                                                         
            N = M + JN                                                       
C
C           COEFFICIENT N=0 UNDEFINED
C
            IF (N .GT. 0) THEN
               CMAT1(IS+JM)  = CMPLX(0.0,ANNP1(N)*REAL(M)) * 
     $               DIVSC(IS+JM)
               CMAT2(IS+JM)  = ANNP1(N)*ZETASC(IS+JM)
               CMAT3(IS+JM)  = CMPLX(0.0,ANNP1(N)*REAL(M)) *
     $               ZETASC(IS+JM)
               CMAT4(IS+JM)  = ANNP1(N)*DIVSC(IS+JM)                         
             ENDIF
  220    CONTINUE                                                            
  230 CONTINUE                                                               
C
C     ADD CORIOLIS TERM FOR ROTATED COORDINATES
C
C     WAVE M=0, N=1
      INDEX = LDIAG(1,2)+1
      CMAT2(INDEX) = CMAT2(INDEX)-CMPLX(ANNP1(1)*CORSC1,0.0) 
C     WAVE M=1, N=1
      INDEX = LDIAG(0,2)+2
      CMAT2(INDEX) = CMAT2(INDEX)-CMPLX(ANNP1(1)*CORSC2,0.0)
      CMAT3(INDEX) = CMAT3(INDEX)-CMPLX(0.0,ANNP1(1)*CORSC2)
C                                                                               
C     EVEN DIAGONALS FIRST                                                      
C                                                                               
      DO 280 JN=0,NN,2                                                          
         IS = LDIAG(JN,2)                                                       
C                                                                               
C     THIS DETOUR ELIMINATES NEED TO EXPLICITLY ZERO UTMP1, UTMP2, ETC.         
C                                                                               
         IF (JN .EQ. 0) THEN                                                   
C
C           SPECIAL HANDLING M=N=0
C
            DO 245 NL=1,NLAT/2
               UTMP1(NL,1) = CMPLX(0.0,0.0)
               UTMP2(NL,1) = CMPLX(0.0,0.0)
               VTMP1(NL,1) = CMPLX(0.0,0.0)
               VTMP2(NL,1) = CMPLX(0.0,0.0)
  245       CONTINUE
            DO 250 JM=2,MM+1
               DO 240 NL=1,NLAT/2                                            
                  UTMP1(NL,JM) = - CMAT1(JM)*ALP(JM,NL)                
                  UTMP2(NL,JM) =   CMAT2(JM)*DALP(JM,NL)               
                  VTMP1(NL,JM) = - CMAT3(JM)*ALP(JM,NL)                
                  VTMP2(NL,JM) = - CMAT4(JM)*DALP(JM,NL)               
  240          CONTINUE                                                      
  250       CONTINUE                                                         
         ELSE
C
C           ALL OTHER EVEN DIAGONALS
C                                                                               
            DO 270 JM=1,LDIAG(JN,1) 
               DO 260 NL=1,NLAT/2 
                  UTMP1(NL,JM) = UTMP1(NL,JM) 
     $                           - CMAT1(IS+JM)*ALP(IS+JM,NL)         
                  UTMP2(NL,JM) = UTMP2(NL,JM) 
     $                           + CMAT2(IS+JM)*DALP(IS+JM,NL)        
                  VTMP1(NL,JM) = VTMP1(NL,JM) 
     $                           - CMAT3(IS+JM)*ALP(IS+JM,NL)         
                  VTMP2(NL,JM) = VTMP2(NL,JM) 
     $                           - CMAT4(IS+JM)*DALP(IS+JM,NL)        
  260          CONTINUE 
  270       CONTINUE   
         ENDIF
  280 CONTINUE                                                                  
C                                                                               
C     ODD DIAGONALS NEXT                                                        
C                                                                               
      DO 310 JN=1,NN,2                                                          
         IS = LDIAG(JN,2)                                                       
         DO 300 JM=1,LDIAG(JN,1)                                                
            DO 290 NL=1,NLAT/2                                                  
               UTMP2(NL,JM) = UTMP2(NL,JM) - CMAT1(IS+JM)*ALP(IS+JM,NL)         
               UTMP1(NL,JM) = UTMP1(NL,JM) + CMAT2(IS+JM)*DALP(IS+JM,NL)        
               VTMP2(NL,JM) = VTMP2(NL,JM) - CMAT3(IS+JM)*ALP(IS+JM,NL)         
               VTMP1(NL,JM) = VTMP1(NL,JM) - CMAT4(IS+JM)*DALP(IS+JM,NL)        
  290       CONTINUE                                                            
  300    CONTINUE                                                               
  310 CONTINUE                                                                  
C                                                                               
C        COMBINE CONTRIBUTIONS OF EVEN AND ODD WAVENUMBERS TO OBTAIN            
C        FOURIER COEFFICIENTS                                                   
C                                                                               
         DO 330 JM=1,MM+1
            DO 320 NL=1,NLAT/2                                                  
               SL = NLAT-NL+1
C              UCOSFC(JM  ,NL)   = UTMP1(NL,JM) + UTMP2(NL,JM) 
               UCOS(2*JM-1,NL)   = REAL(UTMP1(NL,JM) + UTMP2(NL,JM))
               UCOS(2*JM  ,NL)   =AIMAG(UTMP1(NL,JM) + UTMP2(NL,JM))
C              UCOSFC(JM  ,NLAT-NL+1) = UTMP1(NL,JM) - UTMP2(NL,JM) 
               UCOS(2*JM-1,SL)   = REAL(UTMP1(NL,JM) - UTMP2(NL,JM))
               UCOS(2*JM  ,SL)   =AIMAG(UTMP1(NL,JM) - UTMP2(NL,JM))
C              VCOSFC(JM,NL)     = VTMP1(NL,JM) + VTMP2(NL,JM)  
               VCOS(2*JM-1,NL)   = REAL(VTMP1(NL,JM) + VTMP2(NL,JM))
               VCOS(2*JM  ,NL)   =AIMAG(VTMP1(NL,JM) + VTMP2(NL,JM))
C              VCOSFC(JM,NLAT-NL+1)   = VTMP1(NL,JM) - VTMP2(NL,JM)  
               VCOS(2*JM-1,SL)   = REAL(VTMP1(NL,JM) - VTMP2(NL,JM))
               VCOS(2*JM  ,SL)   =AIMAG(VTMP1(NL,JM) - VTMP2(NL,JM))
  320       CONTINUE                                                            
  330    CONTINUE                                                               
C                                                                               
         DO 350 JM=MM+2,NFC
            DO 340 NL=1,NLAT                                                    
C              UCOSFC(JM,NL)   = (0.0,0.0)  
               UCOS(2*JM-1,NL) = 0.0
               UCOS(2*JM  ,NL) = 0.0
C              VCOSFC(JM,NL)   = (0.0,0.0)
               VCOS(2*JM-1,NL) = 0.0
               VCOS(2*JM  ,NL) = 0.0
  340       CONTINUE                                                            
  350    CONTINUE                                                               
C                                                                               
C                                                                               
C     CALL FFT ROUTINE FOR INVERSE TRANSFORM OF UCOS AND VCOS                   
C                                                                               
      CALL FFT991(UCOS,FFTWS2,TRIGS,IFAX,1,NLON+2,NLON,NLAT,+1)               
      CALL FFT991(VCOS,FFTWS2,TRIGS,IFAX,1,NLON+2,NLON,NLAT,+1)               
C                                                                               
      RETURN                                                                    
C                                                                               
      END                                                                       
