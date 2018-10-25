      SUBROUTINE DUV(NNH,DIVSC,ZETASC,LDIAG,TRIGS,ALP,DALP,ANNP1,
     $              CORSC1,CORSC2,RLAT,U,V)  
C                                                                               
C THIS SUBROUTINE OBTAINS CAPITAL U AND V MOMENTUM COMPONENTS FROM THE  
C VORTICITY AND DIVERGENCE SPECTRAL COEFFICIENTS
C IT IS A REDUCED VERSION OF ROUTINE UV FOR A SINGLE POINT
C                                                                              
C CALLED BY: EVAL
C CALLS: DFT991
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
C SPECTRAL TRUNCATION COEFFICIENT
      INTEGER NNH
C DIVERGENCE SPECTRAL COEFFICIENTS
      COMPLEX DIVSC(NALPH)
C VORTICITY SPECTRAL COEFFICIENTS
      COMPLEX ZETASC(NALPH)
C SPECTRAL ARRAY INDEX
      INTEGER LDIAG(0:MAXH,2)
C FOURIER COEFFICIENTS
      REAL TRIGS(MAXH+1,2)
C ASSOC. LEGENDRE POLYNOMIALS
      REAL ALP(NALPH)
C DERIV. ASSOC. LEGENDRE POLYNOMIALS
      REAL DALP(NALPH)
C LAPLACIAN SCALING FACTORS FOR U,V
      REAL ANNP1(MAXH)
C CORIOLIS COEFFICIENTS
      REAL CORSC1,CORSC2
C LATITUDE OF POINT
      REAL RLAT
C
C     Output
C
C COMPUTED EASTWARD WIND
      REAL U
C COMPUTED NORTHWARD WIND
      REAL V
C
C------ Local Variables ------------------------------------------------
C                                                                               
C                                                                               
      COMPLEX CMAT1(NALPH), CMAT2(NALPH), 
     $        CMAT3(NALPH), CMAT4(NALPH) 
      COMPLEX UTMP1(MAXH+1), UTMP2(MAXH+1),             
     $        VTMP1(MAXH+1), VTMP2(MAXH+1)              
      COMPLEX UCOSFC(MAXH+1),VCOSFC(MAXH+1)
C
C     TEMPORARIES
C
      INTEGER M,N,JM,JN,IS,INDEX
      REAL UCOS,VCOS,DIV
C
C----- Executable Statements -------------------------------------------
C                                                                               
C     COMPUTE U AND V MOMENTUM COMPONENTS FROM VORTICITY AND DIVERGENCE         
C                                                                               
C     COMPUTE INTERMEDIATE QUANTITIES FOR COMPUTATIONAL EFFICIENCY           
C                                                                               
      DO 230 JN=0,NNH 
         IS = LDIAG(JN,2)                                                    
C
         DO 220 JM=1,LDIAG(JN,1)                                             
            M = JM-1                                                         
            N = M + JN                                                       
C
C           TERMS N=0 UNDEFINED
C
            IF (N .GT. 0) THEN
               CMAT1(IS+JM)  = ANNP1(N)*REAL(M)*CMPLX(
     $                 -AIMAG(DIVSC(IS+JM)),REAL(DIVSC(IS+JM)))        
               CMAT2(IS+JM)  = ANNP1(N)*ZETASC(IS+JM)
               CMAT3(IS+JM)  = ANNP1(N)*REAL(M)*CMPLX( 
     $               -AIMAG(ZETASC(IS+JM)),REAL(ZETASC(IS+JM)))        
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
      DO 280 JN=0,NNH,2 
         IS = LDIAG(JN,2)                                                       
C                                                                               
C     THIS DETOUR ELIMINATES NEED TO EXPLICITLY ZERO UTMP1, UTMP2, ETC.         
C                                                                               
         IF (JN .EQ. 0) THEN                                                   
C
C           SPECIAL HANDLING FOR N=M=0
C
            UTMP1(1) = CMPLX(0.0,0.0)
            UTMP2(1) = CMPLX(0.0,0.0)
            VTMP1(1) = CMPLX(0.0,0.0)
            VTMP2(1) = CMPLX(0.0,0.0)
            DO 250 JM=2,LDIAG(0,1)                                           
               UTMP1(JM) = - CMAT1(JM)*ALP(JM)                
               UTMP2(JM) =   CMAT2(JM)*DALP(JM)               
               VTMP1(JM) = - CMAT3(JM)*ALP(JM)                
               VTMP2(JM) = - CMAT4(JM)*DALP(JM)               
  250       CONTINUE                                                         
         ELSE
C                                                                               
C           ALL OTHER EVEN DIAGONALS
C
            DO 270 JM=1,LDIAG(JN,1)  
               UTMP1(JM) = UTMP1(JM) - CMAT1(IS+JM)*ALP(IS+JM)         
               UTMP2(JM) = UTMP2(JM) + CMAT2(IS+JM)*DALP(IS+JM)        
               VTMP1(JM) = VTMP1(JM) - CMAT3(IS+JM)*ALP(IS+JM)         
               VTMP2(JM) = VTMP2(JM) - CMAT4(IS+JM)*DALP(IS+JM)        
  270       CONTINUE 
         ENDIF
  280 CONTINUE                                                                  
C                                                                               
C     ODD DIAGONALS NEXT                                                        
C                                                                               
      DO 310 JN=1,NNH,2 
         IS = LDIAG(JN,2)                                                       
         DO 300 JM=1,LDIAG(JN,1)                                                
            UTMP2(JM) = UTMP2(JM) - CMAT1(IS+JM)*ALP(IS+JM)         
            UTMP1(JM) = UTMP1(JM) + CMAT2(IS+JM)*DALP(IS+JM)        
            VTMP2(JM) = VTMP2(JM) - CMAT3(IS+JM)*ALP(IS+JM)         
            VTMP1(JM) = VTMP1(JM) - CMAT4(IS+JM)*DALP(IS+JM)        
  300    CONTINUE                                                               
  310 CONTINUE                                                                  
C                                                                               
C     COMBINE CONTRIBUTIONS OF EVEN AND ODD WAVENUMBERS TO OBTAIN            
C     FOURIER COEFFICIENTS                                                   
C                                                                               
      DO 330 JM=1,LDIAG(0,1)                                                 
         IF (RLAT .GE. 0.0) THEN
            UCOSFC(JM) = UTMP1(JM) + UTMP2(JM)               
            VCOSFC(JM) = VTMP1(JM) + VTMP2(JM)
         ELSE
            UCOSFC(JM) = UTMP1(JM) - UTMP2(JM)               
            VCOSFC(JM) = VTMP1(JM) - VTMP2(JM)
         ENDIF
  330 CONTINUE                                                               
C                                                                               
C     CALL DFT ROUTINE FOR INVERSE TRANSFORM OF UCOS AND VCOS                   
C                                                                               
      CALL DFT991(UCOSFC,TRIGS,MAXH,LDIAG(0,1),UCOS)
      CALL DFT991(VCOSFC,TRIGS,MAXH,LDIAG(0,1),VCOS)
C
C     RENORMALIZE WIND FIELD
C
      DIV = COS(RLAT)
      IF (DIV .NE. 0.0) THEN
         U = UCOS/DIV
         V = VCOS/DIV
      ELSE
         WRITE (0,350) 
  350    FORMAT(/,' STSWM: FATAL ERROR IN ROUTINE DUV',
     $          /,' VELOCITY COMPONENTS MULTIVALUED AT POLES')
         STOP
      ENDIF
C                                                                               
      RETURN                                                                    
C                                                                               
      END                                                                       
