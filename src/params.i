C***********************************************************************
C* INCLUDE FILE 'PARAMS.i'                                             *
C***********************************************************************        
C*                                                                     *        
C*    DEFINITION OF ARRAY SIZES                                        *        
C*    SPECTRAL SPACE TRUNCATION:                                       *
C*                                                                     *        
C*         |                                                           *        
C*         |                                                           *        
C*      KK |- - - - - - - - + + + + + + + + +                          *        
C*         |              + + + + + + + + + +   ^ NORTH-SOUTH          *        
C*         |            + + + + + + + + + + +   | LEGENDRE FUNCTION    *        
C*         |          + + + + + + + + + + + +                          *        
C*         |        + + + + + + + + + + + + .                          *        
C*         |      + + + + + + + + + + + +   .                          *        
C*         |    + + + + + + + + + + + +     .                          *        
C*         |  + + + + + + + + + + + +       .                          *        
C*   ^  NN |+ + + + + + + + + + + +         .                          *        
C*   |     |+ + + + + + + + + + +           .                          *        
C*   |     |+ + + + + + + + + +             .                          *        
C*   N     |+ + + + + + + + +               .                          *        
C*         |+ + + + + + + +                 .                          *        
C*         |+ + + + + + +                   .                          *        
C*       . |+ + + + + +                     .                          *        
C*       . |+ + + + +                       .                          *        
C*       . |+ + + +                         .                          *        
C*       2 |+ + +                           .                          *        
C*       1 |+ +                             .                          *        
C*       0 |+_______________________________._________________         *        
C*          0 1 2 ...                      MM                          *        
C*                                                                     *        
C*                                M --->                               *        
C*              (EAST-WEST FOURIER WAVENUMBER)                         *
C*                                                                     *        
C*    MM, NN, and KK are the spectral truncation parameters illustrated*        
C*    above.  NLAT specifies the number of gaussian latitudes, and NLON*        
C*    specifies the number of longitudes.  For reasons of efficiency   *        
C*    NLAT is required to be an even number. For the FFT transforms,   *
C*    NLON can have only prime factors of 2,3 and 5.                   *
      INTEGER MM,NN,KK,NLAT,NLON
C
C     COMMON RHOMBOIDAL TRUNCATION CASES (KK=MM+NN,MM=NN):
C
C     PARAMETER (MM=15, NN=15, KK=30, NLAT=38, NLON=48)
C
C     COMMON TRIANGULAR TRUNCATION CASES (MM=NN=KK):
C
      PARAMETER (MM=42, NN=42, KK=42, NLAT=64, NLON=128)
C     PARAMETER (MM=63, NN=63, KK=63, NLAT=96, NLON=192)
C     PARAMETER (MM=106, NN=106, KK=106, NLAT=160, NLON=320)
C     PARAMETER (MM=170, NN=170, KK=170, NLAT=256, NLON=512)
C     PARAMETER (MM=213, NN=213, KK=213, NLAT=320, NLON=640)
C
C*    MAXH is the triangular truncation limit for input of spectral    *
C*    coefficients for a previous high-resolution reference run        *
C*    solution stored on a file.                                       *
C
      INTEGER MAXH
C     PARAMETER (MAXH=106)
C
      PARAMETER (MAXH=42)
C     PARAMETER (MAXH=213)
C
C*    RTRUNC is the triangular truncation limit for output of spectral *
C*    coefficients as a reference solution.                            *
C
      INTEGER RTRUNC
C     PARAMETER (RTRUNC=106)
C
      PARAMETER (RTRUNC=42)
C     PARAMETER (RTRUNC=213)
C
C***********************************************************************
C
C     DERIVED PARAMETER VARIABLES - DO NOT CHANGE -
C
C***********************************************************************
C
C                                             NALP is the number of     
C     Associated Legendre Polynomials, and is calculated using the      
C     intermediate parameter LMN. NFC is the number of complex Fourier  
C     coefficients. 
C                                                                       
C     LRM is the length of the Belousov recurrence coefficient matrices 
C     NGRPHS is the amount of storage allocated for energetics graphs   
C     LVLS is the number of time levels for prognostic variables.        
C                                                                       
      INTEGER NFC,LMN,LRM,NALP
      PARAMETER (NFC=(NLON+2)/2)
      PARAMETER (LMN=MM+NN-KK, NALP=(MM+1)*(NN+1)-(LMN**2+LMN)/2)           
      PARAMETER (LRM=(MM+1)*(KK+1)-(MM*MM+MM)/2-1)                              
      INTEGER NGRPHS,LVLS
      PARAMETER (NGRPHS=1100)
      PARAMETER (LVLS=3)
      INTEGER NALPH
      PARAMETER (NALPH=(MAXH+1)*(MAXH+1)-(MAXH**2+MAXH)/2)
C***********************************************************************
C* END INCLUDE FILE                                                    *        
C***********************************************************************        
