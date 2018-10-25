C***********************************************************************
C* INCLUDE FILE 'trnsfm.i'                                             *
C***********************************************************************        
C*    Common Block /TRNSFM/ contains the fields THTA, the actual       *
C*    Gaussian co-latitude angle (radians): WTS, the Gaussian weights; *
C*    ALP, the Associated Legendre Polynomials; DALP, the derivatives  *
C*    of the Associated Legendre Polynomials; and LDIAG in which       *
C*    information on the structure of ALP, DALP, etc., is stored.      *
C*    IFAX contains the factorization of NLON for use by the FFT991    *
C*    routine. TRIGS contains trigonometric function values used by    *
C*    FFT991, and FFTWS2 is the FFT991 work array.  ANNP1 contains the *
C*    expression A/(N*(N+1)) and A2NNP1 contains the expression        *
C*    (N*(N+1))/A**2, where A is the radius of the earth and N is the  *
C*    2-dimensional wavenumber.  WTACSJ is 1.0/(A*COS(LAT)**2).        * 
C*    The spectral coriolis coefficients are CORSC1 and CORSC2.        *
C*    For initialization of these variables see routine INPUT !        *
C***********************************************************************
      INTEGER IFAX(13), LDIAG(0:NN,2)
      REAL
     $          THTA(NLAT), WTS(NLAT), TRIGS(3*NLON/2+1),
     $          ALP(NALP,NLAT/2), DALP(NALP,NLAT/2),
     $          FFTWS2(NLAT*(NLON+1)), ANNP1(KK), A2NNP1(0:KK),
     $          WTACSJ(NLAT),CORSC1,CORSC2
C
      COMMON  / TRNSFM / IFAX, LDIAG
      COMMON  / TRNSFM / 
     $          THTA, WTS, TRIGS, ALP, DALP, 
     $          FFTWS2, ANNP1, A2NNP1, WTACSJ,CORSC1,CORSC2
C***********************************************************************
C* END INCLUDE FILE                                                    * 
C***********************************************************************
