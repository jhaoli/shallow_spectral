C***********************************************************************
C* INCLUDE FILE 'tdvars.i'                                             *
C***********************************************************************        
C*    Common Block /TDVARS/ contains the arrays for the prognostic     *        
C*    variables in both grid point space (LVLS levels in time) and in  *        
C*    wavenumber space (one level in time). The fields contain the     *
C*    vorticity (ZETA), divergence (DIV) and geopotential (PHI).       *
C*    PHIBAR is the global mean steady geopotential.                   *
C*    The diagnostic variables UCOS and VCOS for horizontal velocities *
C*    exist only in grid point space.                                  *
C*    The three variables LNM1, LN, LNP1 are index to the third        *
C*    dimension of the prognostic variables and implement a circular   *
C*    buffer with LNM1 = LN-1 and LNP1 = LN+1.                         *
C***********************************************************************
      INTEGER LNM1, LN, LNP1
      REAL  PHIBAR
      REAL  ZETA(NLON+2,NLAT,LVLS), DIV(NLON+2,NLAT,LVLS),
     $      PHI(NLON+2,NLAT,LVLS),
     $      UCOS(NLON+2,NLAT), VCOS(NLON+2,NLAT)
      COMPLEX ZETASC(NALP), DIVSC(NALP), PHISC(NALP)
      COMMON / TDVARS / LNM1, LN, LNP1, PHIBAR
      COMMON / TDVARS / ZETA, DIV, PHI, UCOS, VCOS
      COMMON / TDVARS / ZETASC, DIVSC, PHISC
C***********************************************************************
C* END INCLUDE FILE                                                    *        
C***********************************************************************        
