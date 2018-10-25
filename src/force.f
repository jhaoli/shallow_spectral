      SUBROUTINE FORCE(ETAFCG,DIVFCG,PHIFCG)
C
C COMPUTES VORTICITY, DIVERGENCE AND GEOPOTENTIAL FORCING FOR
C TEST CASE 4 ON THE GRID
C ALTERNATIVELY, COMPUTES U/V/H-MOMENTUM FORCING
C WHICH FORCING IS DETERMINED BY LOGICAL VARIABLE MOMENT
C = .TRUE.  U/V MOMENTUM FORCING
C = .FALSE. VORTICITY/DIVERGENCE FORCING
C
C CALLED BY: COMP1
C CALLS: BUBFNC, DBUBF, D2BUBF, GLAT, GLON
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
C CONSTANTS & TIMESTEPS
      INCLUDE 'consts.i'
C INITIAL CONDITIONS CASE 4
      INCLUDE 'finit.i'
      INCLUDE 'case4.i'
C
C------------ Arguments ------------------------------------------------
C
C     Output
C                                                                               
C VORTICITY FORCING / U-MOMENTUM*COS(PHI) FORCING
      REAL ETAFCG(NLON+2,NLAT)
C DIVERGENCE FORCING / V-MOMENTUM*COS(PHI) FORCING
      REAL DIVFCG(NLON+2,NLAT)
C GEOPOTENTIAL FORCING / H-MOMENTUM*G FORCING
      REAL PHIFCG(NLON+2,NLAT)         
C
C------ Local Variables ------------------------------------------------
C
C      LOCAL VARIABLES REQUIRED TO EFFICIENTLY COMPUTE FORCING:  A LOT          
C      OF SPACE HERE, BUT NECESSARY FOR TRACTABLE AND DEBUGGABLE PROBLEM        
C                                                                               
      REAL BUB, DBUB, D2BUB	
      REAL C     (NLON), PSIB  (NLON)                        
      REAL DCDM  (NLON), DCDL  (NLON),                                     
     $     D2CDM (NLON), D2CDL (NLON),                                     
     $     D3CDM (NLON), D3CDL (NLON),                                     
     $     DMDCDL(NLON), DMD2CL(NLON),                                     
     $     DKDM  (NLON), DKDL  (NLON),                                     
     $     D2KDM (NLON), D2KDL (NLON),                                     
     $     D3KDM (NLON), D3KDL (NLON),                                     
     $     DLDKDM(NLON), DLD2KM(NLON),                                     
     $     DMD2KL(NLON), DLD2CM(NLON)                                      
      REAL TMP1  (NLON), TMP2  (NLON), TMP3  (NLON),                       
     $     UT    (NLON), VT    (NLON),                                     
     $     DUTDL (NLON), DUTDM (NLON),                                     
     $     DVTDL (NLON), DVTDM (NLON)                                      
C                                                                               
C     TEMPORARIES
C
      INTEGER I,J
      REAL TMSHFT, DFDM, AI, A2I, RLAT, SNJ, CSJ, COR, SRCSJ,
     $     TMPRY, CSJI, CSJ2I, ACSJI, AACSJI, ACSJ2I, RLON
C
C----- External Functions ----------------------------------------------
C
C     ZONAL FLOW FUNCTIONS AND DERIVATIVES
C
      EXTERNAL BUBFNC, DBUBF, D2BUBF
      REAL BUBFNC, DBUBF, D2BUBF
C
C     SPATIAL LAT./LONG. GRID
C
      EXTERNAL GLAT, GLON
      REAL GLAT, GLON
C
C----- Executable Statements -------------------------------------------
C
      TMSHFT = SU0*REAL(NSTEP-1)*DT/A                                          
      DFDM   = 2.0*OMEGA                                                        
      AI    = 1.0/A                                                             
      A2I   = 1.0/(A*A)                                                         
C                                                                               
      DO 510 J=1,NLAT                                                           
         RLAT = GLAT(J)
C
C        LATITUDE = RLAT = GLAT(J)
C
         SNJ  = SIN(RLAT)
         CSJ  = COS(RLAT)*COS(RLAT)
C
C        STEADY ZONAL FLOW BUB AND ITS TWO DERIVATIVES
C
         BUB    = BUBFNC(RLAT)*COS(RLAT)
         DBUB   = DBUBF(RLAT)
         D2BUB  = D2BUBF(RLAT)
C
C        CORIOLIS PARAMETER
C
         COR    = 2.0*OMEGA*SNJ
C
C        LONGITUDE INDEPENDENT FACTORS
C
         SRCSJ  = COS(RLAT)
         TMPRY  = TAN(RLAT)
C        TMPRY2 = TMPRY*TMPRY                                                   
         CSJI   = 1.0/CSJ                                                    
         CSJ2I  = 1.0/(CSJ**2)                                               
         ACSJI  = 1.0/(A*CSJ)                                                
         AACSJI = 1.0/(A*A*CSJ)                                              
         ACSJ2I = 1.0/((A*CSJ)**2)                                           
C                                                                               
         DO 100 I=1,NLON                                                        
            RLON = GLON(I)
C
C           LONGITUDE = RLON = GLON(I)
C           LATITUDE = RLAT = GLAT(J)
C
C           CALCULATE C AND ALL OF ITS DERIVATIVES FOR THIS TIME STEP
C
C (A.28)
            C(I)      = SIN(RLAT0)*SNJ + COS(RLAT0)*SRCSJ                    
     $                                     *COS(RLON-TMSHFT-RLON0)          
C (A.38)
            DCDM(I)   = SIN(RLAT0) - COS(RLON-TMSHFT-RLON0)*                
     $                      COS(RLAT0)*TMPRY                                    
C (A.39)
            DCDL(I)   = -COS(RLAT0)*SRCSJ*SIN(RLON-TMSHFT-RLON0)            
C (A.40)
C SIMPLIFIED TERM
            D2CDM(I)  = -COS(RLAT0)*COS(RLON-TMSHFT-RLON0)*                 
     $                      CSJI/SRCSJ                             
C (A.41)
            D2CDL(I)  = -COS(RLAT0)*SRCSJ*COS(RLON-TMSHFT-RLON0)            
C (A.42)
C SIMPLIFIED TERM
            D3CDM(I)  = D2CDM(I)*3.0*SNJ*CSJI
C
C (A.43)
            D3CDL(I)  = -DCDL(I)                                                
C (A.44)
            DMDCDL(I) = +COS(RLAT0)*SIN(RLON-TMSHFT-RLON0)*TMPRY            
C (A.45)
            DLD2CM(I) = +COS(RLAT0)*SIN(RLON-TMSHFT-RLON0)*                 
     $                       CSJI/SRCSJ
C (A.46)
            DMD2CL(I) = +COS(RLAT0)*COS(RLON-TMSHFT-RLON0)*TMPRY                                                
  100    CONTINUE                                                               
C                                                                               
C     CALCULATE PSI BAR AND ALL OF ITS DERIVATIVES (DKDM, ETC.)               
C                                                                               
         DO 120 I=1,NLON                                                        
C (A.27)
            PSIB(I)   = ALFA*EXP(-SIGMA*((1.0-C(I))/(1.0+C(I))))                
            TMP1(I)   = 2.0*SIGMA*PSIB(I)/((1.0 + C(I))**2)                     
            TMP2(I)   = (SIGMA - (1.0 + C(I)))/((1.0 + C(I))**2)                
            TMP3(I)   = (((1.0+C(I))**2)-2.0*SIGMA*(1.0+C(I)))                  
     $                      /((1.0 + C(I))**4)                                  
C (A.29)
            DKDM(I)   = TMP1(I)*DCDM(I)                                         
C (A.30)
            DKDL(I)   = TMP1(I)*DCDL(I)                                         
C (A.31)
            D2KDM(I)  = TMP1(I)*(D2CDM(I) + 2.0*(DCDM(I)**2)                    
     $                    *TMP2(I))                                             
C (A.32)
            D2KDL(I)  = TMP1(I)*(D2CDL(I) + 2.0*(DCDL(I)**2)                    
     $                    *TMP2(I))                                             
C (A.33)
            D3KDM(I)  = TMP1(I)*(D3CDM(I) + 2.0*(DCDM(I)**3)                    
     $                    *TMP3(I) + 2.0*DCDM(I)*TMP2(I)*(3.0                   
     $                    *D2CDM(I) + 2.0*(DCDM(I)**2)*TMP2(I)))                
C (A.34)
            D3KDL(I)  = TMP1(I)*(D3CDL(I) + 2.0*(DCDL(I)**3)                    
     $                    *TMP3(I) + 2.0*DCDL(I)*TMP2(I)*(3.0                   
     $                    *D2CDL(I) + 2.0*(DCDL(I)**2)*TMP2(I)))                
C (A.35)
            DLDKDM(I) = TMP1(I)*(DMDCDL(I) + 2.0*DCDL(I)                        
     $                    *DCDM(I)*TMP2(I))                                     
C (A.37)
            DMD2KL(I) = TMP1(I)*(DMD2CL(I) + 2.0*(DCDL(I)**2)                   
     $                    *DCDM(I)*TMP3(I) + 2.0*DCDM(I)*TMP2(I)                
     $                    *(D2CDL(I) + 2.0*(DCDL(I)**2)*TMP2(I))                
     $                    + 4.0*DCDL(I)*DMDCDL(I)*TMP2(I))                      
C (A.36)
            DLD2KM(I) = TMP1(I)*(DLD2CM(I) + 2.0*(DCDM(I)**2)                   
     $                    *DCDL(I)*TMP3(I) + 2.0*DCDL(I)*TMP2(I)                
     $                    *(D2CDM(I) + 2.0*(DCDM(I)**2)*TMP2(I))                
     $                    + 4.0*DCDM(I)*DMDCDL(I)*TMP2(I))                      
  120    CONTINUE                                                               
C                                                                               
C     COMPUTE COMMONLY UTILIZED TERMS IN FORCING (U AND V TILDE)                
C                                                                               
         DO 160 I=1,NLON                                                        
            UT(I)    = BUB - CSJ*DKDM(I)*AI                               
            VT(I)    = DKDL(I)*AI                                               
            DUTDL(I) = -CSJ*DLDKDM(I)*AI                                     
            DVTDL(I) = D2KDL(I)*AI                                              
            DUTDM(I) = DBUB - (CSJ*D2KDM(I) - 2.0*SNJ*DKDM(I))*AI               
            DVTDM(I) = DLDKDM(I)*AI                                             
  160    CONTINUE                                                               
C                                                                               
C     COMPUTE FORCING TERMS                                                     
C                                                                               
         DO 500 I=1,NLON                                                        
C                                                                               
            IF (MOMENT) THEN
C
C           U-MOMENTUM FORCING
C
            TMP1(I) = (CSJ*SU0)/(A*A)*DLDKDM(I)
            TMP1(I) = TMP1(I) + UT(I)*ACSJI*DUTDL(I)
            TMP1(I) = TMP1(I) + VT(I)*AI*DUTDM(I)
C
            ETAFCG(I,J) = TMP1(I) + COR*(AI*DKDL(I)-VT(I))
C
            ELSE
C
C           FORCING ON VORTICITY EQUATION                                       
C           COMPUTE LAMBDA DERIVATIVE 1ST (CAN BE USED FOR TIME DRVTVE)         
C                                                                               
C (A.16,A.17)
            TMP1(I)     = D3KDL(I)*AACSJI + (CSJ*DLD2KM(I)                   
     $                    - 2.0*SNJ*DLDKDM(I))*A2I                           
            ETAFCG(I,J) = (UT(I)*ACSJI - (SU0*AI))*TMP1(I)                      
C (A.18)
            TMP2(I)     = (DMD2KL(I)*CSJI + 2.0*SNJ                          
     $                    *D2KDL(I)*CSJ2I + CSJ*D3KDM(I)                     
     $                    - 4.0*SNJ*D2KDM(I) - 2.0*DKDM(I))*A2I              
     $                    + DFDM - D2BUB*AI                                  
C (A.13)
            ETAFCG(I,J) = ETAFCG(I,J) + TMP2(I)*VT(I)*AI                        
            ENDIF
C                                                                               
C           COMPUTE FORCING ON GEOPOTENTIAL EQUATION                            
C                                                                               
C (A.25)
            TMP3(I)     = COR*DKDL(I)
            PHIFCG(I,J) = -SU0*AI*TMP3(I) + UT(I)*TMP3(I)*ACSJI 
            PHIFCG(I,J) = PHIFCG(I,J) + AI*VT(I)*(PSIB(I)*DFDM+
     $                    COR*DKDM(I))
C
C     ADD TERM TO BALANCE ZONAL FLOW
C
            PHIFCG(I,J) = PHIFCG(I,J) - VT(I)*BUB*CSJI*
     $                    (COR + BUB*ACSJI*SNJ)
C                                                                               
            IF (MOMENT) THEN
C
C           V-MOMENTUM FORCING
C
            TMP3(I) = - SU0*A2I*D2KDL(I)
            TMP3(I) = TMP3(I) + UT(I)*ACSJI*DVTDL(I)
            TMP3(I) = TMP3(I) + VT(I)*AI*DVTDM(I)
            TMP3(I) = TMP3(I) + (UT(I)*UT(I)+VT(I)*VT(I))*
     $                SNJ*ACSJI
C
            DIVFCG(I,J) = TMP3(I) + CSJ*AI*(COR*DKDM(I)+
     $                  PSIB(I)*DFDM)
            DIVFCG(I,J) = DIVFCG(I,J) + COR*(UT(I)-BUB) -
     $                  ACSJI*SNJ*BUB*BUB
C
            ELSE
C           
C           COMPUTE FORCING ON DIVERGENCE EQUATION (MOST MESSY!)                
C           LAPLACIAN OF PHI FIRST                                              
C                                                                               
C (A.19)
            DIVFCG(I,J) = COR*D2KDL(I)*AACSJI + ( CSJ*
     $                    (COR*D2KDM(I) + 2.0*DFDM*DKDM(I)) - 2.0
     $                    * SNJ*(PSIB(I)*DFDM + COR*DKDM(I)))*A2I 
C                                                                               
C           CALCULATE ETA-F NEXT 
C                                                                               
            TMP3(I)     = ACSJI*DVTDL(I)-AI*DUTDM(I)
C                                                                               
            DIVFCG(I,J) = DIVFCG(I,J) - TMP3(I)*COR + UT(I)*AI*DFDM
C                                                                               
C           FINALLY, ADD LAPLACIAN OF U**2/(2*CSJ) AND V**2/(2*CSJ)             
C                                                                               
            DIVFCG(I,J) = DIVFCG(I,J) + 2.0*SNJ*AACSJI*(
     $                    UT(I)*CSJI*DVTDL(I) - VT(I)*CSJI*DUTDL(I) 
     $                    + UT(I)*DUTDM(I) + VT(I)*DVTDM(I) )
C                                                                               
            DIVFCG(I,J) = DIVFCG(I,J) + 2.0*AACSJI*(
     $                    DUTDM(I)*DVTDL(I) - DUTDL(I)*DVTDM(I) )
C
            DIVFCG(I,J) = DIVFCG(I,J) + (UT(I)*UT(I)+VT(I)*VT(I))
     $                    *(1.0+SNJ*SNJ)*ACSJ2I
C                                                                               
C     ADD TERM TO BALANCE ZONAL FLOW
C
            DIVFCG(I,J) = DIVFCG(I,J) - COR*AI*DBUB - BUB*AI*DFDM -
     $                    ACSJ2I*(BUB*BUB*(1+SNJ*SNJ)+
     $                    2.0*SNJ*CSJ*BUB*DBUB)
C
            ENDIF
C
  500    CONTINUE                                                               
  510 CONTINUE                                                                  
C                                                                               
      RETURN                                                                    
      END                                                                       
