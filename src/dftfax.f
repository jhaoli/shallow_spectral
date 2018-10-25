      SUBROUTINE DFTFAX(MAXH,NFC,LAMBDA,TRIGS)
C
C COMPUTE COEFFICIENTS FOR INVERSE DISCRETE FOURIER TRANSFORM
C
C CALLED BY: EVAL
C CALLS:
C
C REVISIONS:
C 7-10-92 CHANGE TO CCM CODING CONVENTIONS (R. JAKOB)
C
C------------ Arguments ------------------------------------------------
C
C     Input
C
C ARRAY DIMENSION OF TRIGS
      INTEGER MAXH
C NUMBER OF FOURIER COEFFICIENTS
      INTEGER NFC
C LONGITUDE
      REAL LAMBDA
C
C     Output
C
C SINE AND COSINE COEFFICIENTS
      REAL TRIGS(MAXH+1,2)
C
C------ Local Variables ------------------------------------------------
C
      INTEGER I
C
C----- Executable Statements -------------------------------------------
C
      DO 100 I = 1,NFC
         TRIGS(I,1) = SIN((I-1)*LAMBDA)
         TRIGS(I,2) = COS((I-1)*LAMBDA)
  100 CONTINUE
C
      RETURN
      END
