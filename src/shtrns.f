      SUBROUTINE SHTRNS(L3D,LV,ITYPE,X,SCOEF)                
C                                                                              
C THIS SUBROUTINE PERFORMS SPHERICAL HARMONIC TRANSFORMS AND INVERSE  
C TRANFORMS OF ARBITRARY SCALAR DATA USING THE GAUSSIAN GRID (SEE           
C THE GUASSIAN LATITUDES AND WEIGHTS) AND THE ASSOCIATED LEGENDRE           
C POLYNOMIALS AND DERIVATIVES SPECIFIED IN COMMON BLOCK TRNSFM.             
C IT USES THE LIBRARY ROUTINE FFT991 FOR IN-PLACE FAST FOURIER
C TRANSFORMS.
C                                                                              
C CALLED BY: COMP1, INIT, STSWM
C CALLS: FFT991
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
C WORKSPACE
      INCLUDE 'wrkspc.i'
C
C------------ Arguments ------------------------------------------------
C
C     Input
C
C NUMBER OF TIME LEVELS OF X
      INTEGER L3D
C CURRENT TIME LEVEL FOR FIELD X
      INTEGER LV
C SPECIFIES TYPE OF TRANSFORM DESIRED 
C          -1 => FORWARD TRANSFORM X -> SCOEF
C          +1 => INVERSE TRANSFORM SCOEF -> X
      INTEGER ITYPE
C
C     Input/Output
C
C SCALAR DATA (ON GAUSSIAN GRID)
      REAL X(NLON+2,NLAT,L3D)
C SPHERICAL HARMONIC COEFFICIENT ARRAY
      COMPLEX SCOEF(NALP)
C                                                                               
C------ Local Variables ------------------------------------------------
C
C     EQUIVALENCE TO WORKSPACE IN STORAGE POOL /WRKSPC/, ETC.                   
C                                                                               
      COMPLEX   CEVEN(NLAT/2,NFC), CODD (NLAT/2,NFC)            
      EQUIVALENCE (CEVEN(1,1),WS4(1,1)), (CODD(1,1),WS4(1,NLAT/2+1))            
C
      COMPLEX   CTMP1(MM+1), CTMP2(MM+1)                                      
C
C     TEMPORARIES
C
      INTEGER NL, JM, JN, IS, SL
C                                                                               
C----- Executable Statements -------------------------------------------
C                                                                               
C     FIRST CHECK FOR VALID ARGUMENTS; INVALID ARGUMENTS=>FATAL ERROR           
C                                                                               
      IF ((ITYPE .NE. +1) .AND. (ITYPE .NE. -1)) THEN                           
         WRITE (0,900) ITYPE                                                    
  900    FORMAT(/,' STSWM: FATAL ERROR IN SUBROUTINE SHTRNS:',/,
     $            ' UNKNOWN TYPE OF TRANSFORM SPECIFIED',/, 
     $          ' ITYPE = ',I3)                                              
         STOP                                                                   
      ENDIF                                                                     
C                                                                               
      IF (ITYPE .EQ. -1) THEN                                                   
C                                                                               
C     FORWARD TRANSFORM FROM PHYSICAL TO WAVENUMBER SPACE (ITYPE=-1)         
C                                                                               
C     IN PLACE FFT FOR ALL LATITUDES OF TIMELEVEL LV
C
         CALL FFT991(X(1,1,LV),FFTWS2,TRIGS,IFAX,1,NLON+2,
     $                  NLON,NLAT,-1)              
C                                                                               
C        PROCEDURE FOR FORWARD GAUSS-LEGENDRE TRANSFORM.                        
C        VARY M AND N SO THAT PROCEDURE MOVES ALONG DIAGONALS DENOTED           
C        BY INDEX JN.  M IS GIVEN BY (JM-1) WHILE N IS GIVEN BY (JN+M).         
C        SINCE NUMBER OF GAUSSIAN LATITUDES IS EVEN, TAKE ADVANTAGE OF          
C        THE SYMMETRIC CHARACTER OF THE ASSOCIATED LEGENDRE POLYNOMIALS         
C                                                                               
         DO 100 NL=1,NLAT/2                                                     
            SL = NLAT-NL+1
C
            DO 30 JM=1,MM+1
C              CTMP1(JM) = (X(JM,NL,LV) + X(JM,SL,LV))*WTS(NL) 
               CTMP1(JM) = (CMPLX(X(2*JM-1,NL,LV),X(2*JM,NL,LV)) +
     $            CMPLX(X(2*JM-1,SL,LV),X(2*JM,SL,LV)))
     $            * WTS(NL)
C              CTMP2(JM) = (X(JM,NL,LV) - X(JM,SL,LV))*WTS(NL)
               CTMP2(JM) = (CMPLX(X(2*JM-1,NL,LV),X(2*JM,NL,LV)) -
     $            CMPLX(X(2*JM-1,SL,LV),X(2*JM,SL,LV)))
     $            * WTS(NL)

   30       CONTINUE                                                          
C                                                                               
C           SPECIAL HANDLING OF THE FIRST LATITUDE IN THE ACCUMULATION          
C           THIS DETOUR ELIMINATES THE NEED TO SEPARATELY ZERO SCOEF            
C                                                                               
            IF (NL .EQ. 1) THEN                                                 
C                                                                               
C           EVEN WAVENUMBERS FIRST                                              
C                                                                               
               DO 50 JN=0,NN,2                                                  
                  IS = LDIAG(JN,2)                                              
C
                  DO 45 JM=1,LDIAG(JN,1)                                        
                     SCOEF(IS+JM) = CTMP1(JM)*ALP(IS+JM,NL)                  
   45             CONTINUE                                                      
   50          CONTINUE                                                         
C                                                                               
C           ODD WAVENUMBERS NEXT                                                
C                                                                               
               DO 60 JN=1,NN,2                                                  
                  IS = LDIAG(JN,2)                                              
C
                  DO 55 JM=1,LDIAG(JN,1)                                        
                     SCOEF(IS+JM) = CTMP2(JM)*ALP(IS+JM,NL)                  
   55             CONTINUE                                                      
   60          CONTINUE                                                         
            ELSE                                                               
C                                                                               
C           REMAINING LATITUDES CONTRIBUTING TO SUM                             
C           EVEN WAVENUMBERS                                                    
C                                                                               
            DO 75 JN=0,NN,2                                                     
               IS = LDIAG(JN,2)                                                 
C
               DO 70 JM=1,LDIAG(JN,1)                                           
                  SCOEF(IS+JM) = SCOEF(IS+JM) +                           
     $                              CTMP1(JM)*ALP(IS+JM,NL)                     
   70          CONTINUE                                                         
   75       CONTINUE                                                            
C                                                                               
C           ODD WAVENUMBERS                                                     
C                                                                               
            DO 90 JN=1,NN,2                                                     
               IS = LDIAG(JN,2)                                                 
C
               DO 80 JM=1,LDIAG(JN,1)                                           
                  SCOEF(IS+JM) = SCOEF(IS+JM) +                           
     $                              CTMP2(JM)*ALP(IS+JM,NL)                     
   80          CONTINUE                                                         
   90       CONTINUE                                                            
C 
            ENDIF
  100    CONTINUE                                                               
C                                                                               
C        TRANSFORMATION TO WAVENUMBER SPACE (FORWARD TRANSFORM) COMPLETE        
C                                                                               
         RETURN                                                                 
C                                                                               
      ELSE                                                                      
C                                                                               
C        INVERSE TRANSFORM FROM WAVENUMBER TO PHYSICAL SPACE (ITYPE=+1)         
C        DETERMINE FOURIER COEFFICIENTS BY INVERSE LEGENDRE TRANSFORM.          
C        VARY M AND N SO PROCEDURE MOVES ALONG DIAGONALS DENOTED BY             
C        INDEX JN.  M IS GIVEN BY (JM-1) WHILE N IS GIVEN BY (JN+M).            
C        FIRST ACCUMULATE EVEN WAVENUMBER CONTRIBUTION                          
C                                                                               
         DO 180 JN=0,NN,2                                                       
            IS = LDIAG(JN,2)                                                    
C                                                                               
C           THIS DETOUR HELPS AVOID THE NEED TO SEPARATELY ZERO CEVEN           
C                                                                               
            IF (JN.EQ.0) THEN                                                   
               DO 150 JM=1,MM+1
                  DO 140 NL=1,NLAT/2                                            
                     CEVEN(NL,JM) = SCOEF(IS+JM)*ALP(IS+JM,NL)               
  140             CONTINUE                                                      
  150          CONTINUE                                                         
C                                                                               
            ELSE                                                               
C                                                                               
               DO 170 JM=1,LDIAG(JN,1) 
                  DO 160 NL=1,NLAT/2 
                     CEVEN(NL,JM) = CEVEN(NL,JM) +  
     $                             SCOEF(IS+JM)*ALP(IS+JM,NL)                
  160             CONTINUE 
  170          CONTINUE 
            ENDIF
  180    CONTINUE                                                               
C                                                                               
C        NEXT ACCUMULATE ODD WAVENUMBER CONTRIBUTION                            
C                                                                               
         DO 215 JN=1,NN,2                                                       
            IS = LDIAG(JN,2)                                                    
C                                                                               
C           THIS DETOUR HELPS AVOID THE NEED TO SEPARATELY ZERO CODD            
C                                                                               
            IF (JN.EQ.1) THEN                                                   
               DO 190 JM=1,LDIAG(1,1)                                           
                  DO 185 NL=1,NLAT/2                                            
                     CODD (NL,JM) = SCOEF(IS+JM)*ALP(IS+JM,NL)               
  185             CONTINUE                                                      
  190          CONTINUE                                                         
C                                                                               
C           ACCOUNT FOR THE FACT THAT THE FIRST ODD DIAGONAL MAY BE             
C           SHORTER THAN THE FIRST EVEN DIAGONAL (PART OF THE GAME              
C           TO AVOID EXPLICITLY ZEROING THE ENTIRE CODD ARRAY)                  
C                                                                               
               IF (LDIAG(1,1) .LT. MM+1) THEN                             
                  DO 200 JM=LDIAG(1,1)+1, MM+1
                     DO 195 NL=1,NLAT/2                                         
                        CODD (NL,JM) =  (0.0,0.0)                               
  195                CONTINUE                                                   
  200             CONTINUE                                                      
C                                                                               
               ENDIF                                                            
            ELSE                                                               
C                                                                               
               DO 210 JM=1,LDIAG(JN,1)  
                  DO 205 NL=1,NLAT/2 
                     CODD (NL,JM) = CODD (NL,JM) +  
     $                             SCOEF(IS+JM)*ALP(IS+JM,NL)                
  205             CONTINUE  
  210          CONTINUE 
C
            ENDIF
  215    CONTINUE                                                               
C                                                                               
C        COMBINE CONTRIBUTIONS OF EVEN AND ODD WAVENUMBERS TO OBTAIN            
C        COMPLEX FOURIER COEFFICIENTS, FOLLOWED BY INVERSE FFT                  
C                                                                               
         DO 225 JM=1,MM+1
            DO 220 NL=1,NLAT/2                                               
               SL = NLAT-NL+1
C              X(JM,NL,LV)        = CEVEN(NL,JM) + CODD(NL,JM) 
               X(2*JM-1,NL,LV) = REAL(CEVEN(NL,JM)+CODD(NL,JM))
               X(2*JM,NL,LV)   = AIMAG(CEVEN(NL,JM)+CODD(NL,JM))
C              X(JM,NLAT-NL+1,LV) = CEVEN(NL,JM) - CODD(NL,JM)
               X(2*JM-1,SL,LV) = REAL(CEVEN(NL,JM)-CODD(NL,JM))
               X(2*JM,SL,LV)   = AIMAG(CEVEN(NL,JM)-CODD(NL,JM))
  220       CONTINUE                                                         
  225    CONTINUE                                                            
C                                                                               
C        ZERO THE TAIL OF THE COMPLEX COEFFICIENT SPECTRUM                   
C                                                                               
         DO 235 JM=MM+2,NFC
            DO 230 NL=1,NLAT                                                 
C              X(JM,NL,LV) = (0.0,0.0)
               X(2*JM-1,NL,LV) = 0.0
               X(2*JM,NL,LV)   = 0.0
  230       CONTINUE                                                         
  235    CONTINUE                                                            
C                                                                               
C        INVERSE FAST FOURIER TRAANSFORM FOR ALL LATITUDES
C        AND TIMELEVEL LV
C
         CALL FFT991(X(1,1,LV),FFTWS2,TRIGS,IFAX,1,NLON+2,NLON,
     $                  NLAT,+1)              
C                                                                               
C        TRANSFORMATION TO PHYSICAL SPACE (INVERSE TRANSFORM) COMPLETE          
C                                                                               
         RETURN                                                                 
C                                                                               
      ENDIF                                                                     
C                                                                               
      END                                                                       
