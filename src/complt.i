C***********************************************************************
C* INCLUDE FILE 'complt.i'                                             *
C***********************************************************************        
C*    COMMON BLOCK /COMPLT/ CONTAINS INFO REQUIRED BY PLOTTING ROUTINES*        
C*                                                                     *        
C*    LCONT  = .TRUE. => DRAW CONTINENTAL OUTLINES                     *        
C*    LPSP   = .TRUE. => DRAW POLAR STEREOGRAPHIC PROJECTION           *        
C*    LOP    = .TRUE. => DRAW ORTHOGRAPHIC PROJECTION                  *        
C*    LCP    = .TRUE. => DRAW CYLINDRICAL EQUIDISTANT PROJECTION       *        
C*    LNORTH = .TRUE. => NORTH/ SOUTH POLAR STEREOGRAPHIC PROJECTION   *
C*    LG     = .TRUE. => DRAW HEIGHT FIELD                             *        
C*    LU     = .TRUE. => DRAW SMALL U FIELD                            *        
C*    LV     = .TRUE. => DRAW SMALL V FIELD                            *        
C*    LZ     = .TRUE. => DRAW VORTICITY FIELD                          *        
C*    LD     = .TRUE. => DRAW DIVERGENCE FIELD                         *        
C*    LVV    = .TRUE. => DRAW VECTOR VELOCITY FIELD                    *        
C*    LVVG   = .TRUE. => DRAW VECTOR VELOCITY OVERLAID BY GEOPOTENTIAL *        
C*                                                                     *        
C*    POLAT, POLNG, AND POROT SPECIFY CENTER LATITUDE, LONGITUDE AND   *        
C*    ROTATION ANGLE FOR THE PLOTTING ROUTINE MAPROJ IN PLOTS.         *
C*    SPV & SPVAL HOLD THE CONREC SPECIAL VALUE FEATURE FOR MASKING    *        
C*    VVSF IS A SCALING PARAMETER FOR VECTOR PLOTS                     *
C*    RLOND AND RLATD CONTAIN LONGITUDE AND LATITUDE OF GAUSSIAN GRID  *
C*    ON THE SPHERE. THEY ARE USED IN THE MAPPING FUNCTIONS            *
C*    CPMPXY AND CPMVXY OF THE GRAPHICS PACKAGE                        *
C*    SEE NCAR SYSTEM PLOT PACKAGE DOCUMENTATION FOR EXACT DEFINITIONS *        
C*    THESE VARIABLES ARE INITIALIZED IN SUBROUTINE INPUT              *
C***********************************************************************        
      REAL
     $      SPV(2), SPVAL, VVSF,
     $      POLAT, POLNG, POROT,
     $      RLOND(NLON+1), RLATD(NLAT)
      LOGICAL
     $      LCONT, LPSP, LOP, LCP, LNORTH, LG, LU, LV, LZ, LD,
     $      LVV, LVVG
      COMMON /COMPLT/ 
     $      SPV, SPVAL, VVSF, POLAT, POLNG, POROT,
     $      RLOND, RLATD,
     $      LCONT, LPSP, LOP, LCP, LNORTH, LG, LU, LV, LZ, LD, 
     $      LVV, LVVG
C***********************************************************************
C* END INCLUDE FILE                                                    *        
C***********************************************************************        
