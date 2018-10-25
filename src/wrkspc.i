C***********************************************************************
C* INCLUDE FILE 'wrkspc.i'                                             *
C***********************************************************************        
C*    Common Block /WRKSPC/ contains a large number of scratch arrays  *        
C*    used by various routines in this spectral shallow water model.   *        
C*    It provides a pool of memory for temporary storage by model codes*        
C***********************************************************************        
      REAL RMAT1(2*NALP), RMAT2(2*NALP), RMAT3(2*NALP), RMAT4(2*NALP),
     $     WS1(NLON+2,NLAT), WS2(NLON+2,NLAT), WS3(NLON+2,NLAT),
     $     WS4(NLON+2,NLAT), WS5(NLON+2,NLAT), WS6(NLON+2,NLAT),
     $     WS7(NLON+2,NLAT), WS8(NLON+2,NLAT),
     $     RVEC1(NLON+2), RVEC2(NLON+2), RVEC3(NLON+2)
      COMMON  / WRKSPC / 
     $       RMAT1, RMAT2, RMAT3, RMAT4,        
     $       WS1, WS2, WS3, WS4, WS5, WS6, WS7, WS8,
     $       RVEC1, RVEC2, RVEC3                        
C***********************************************************************
C* END INCLUDE FILE                                                    *        
C***********************************************************************        
