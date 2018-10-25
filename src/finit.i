C***********************************************************************
C* INCLUDE FILE 'finit.i'                                              *
C***********************************************************************        
C*    COMMON BLOCK /CONST2/ CONTAINS ALL INFORMATION REQUIRED TO MAKE  *        
C*    USE OF INITIAL CONDITIONS PROVIDED AS TEST CASES IN SUBROUTINE   *
C*    INIT. THESE INITIAL CONDITIONS ARE THEN USED BY ANLYTC TO GET    *
C*    ANALYTIC SOLUTIONS.                                              *
C*    UIC12,VIC12: initial u,v wind                                    *
C*    PIC12: initial height                                            *
C*    DIC12,EIC12: initial divergence,vorticity                        *
C*    MOUNT: surface height (mountains) for CASE 5                     *
C*    TOPOSC: Spectral coefficients of mountains                       *
C***********************************************************************        
      INTEGER R
      REAL ALFA, SIGMA, RLAT0, RLON0, PHI0, K, OMG,
     $    PHICON(NLAT),PHIA(NLAT),PHIB(NLAT),PHIC(NLAT),
     $    UIC12(NLON+2,NLAT), VIC12(NLON+2,NLAT), PIC12(NLON+2,NLAT),
     $    DIC12(NLON+2,NLAT), EIC12(NLON+2,NLAT),
     $    MOUNT(NLON+2,NLAT), UCON(NLAT), VCON(NLAT)
      COMPLEX TOPOSC(NALP)
C
      COMMON  / CONST2 / R
      COMMON  / CONST2 / ALFA, SIGMA, RLAT0, RLON0, PHI0, K, OMG, 
     $    PHICON,PHIA,PHIB,PHIC,
     $    UIC12, VIC12, PIC12, DIC12, EIC12,
     $    MOUNT, UCON, VCON
      COMMON  / CONST2 / TOPOSC
C***********************************************************************
C* END INCLUDE FILE                                                    *        
C***********************************************************************
