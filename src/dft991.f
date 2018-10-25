      SUBROUTINE DFT991(X,TRIGS,MAXH,NFC,SUM)
C
C INVERSE DISCRETE FOURIER TRANSFORM 
C
C CALLED BY: DTRNS, DUV
C CALLS:
C
C REVISIONS:
C 7-10-92 CHANGE TO CCM CODING CONVENTIONS (R. JAKOB)
C
C------------ Arguments ------------------------------------------------
C
C     Input
C
C COMPLEX FOURIER COEFFICIENTS
      COMPLEX X(MAXH+1)
C SINE AND COSINE COEFFICIENTS
      REAL TRIGS(MAXH+1,2)
C DIMENSION OF TRIGS
      INTEGER MAXH
C NUMBER OF VALID COEFFICIENTS
      INTEGER NFC
C
C     Output
C
      REAL SUM
C
C------ Local Variables ------------------------------------------------
C
      INTEGER I
C
C----- Executable Statements -------------------------------------------
C
      SUM = REAL(X(1))
      DO 100 I = 2, NFC
         SUM = SUM + 2.0*(REAL(X(I))*TRIGS(I,2) -
     $                    AIMAG(X(I))*TRIGS(I,1))
  100 CONTINUE
C
      RETURN
C
      END
