      REAL FUNCTION GLON(I)
C
C THIS FUNCTION RETURNS THE LONGITUDE OF THE REGULAR
C GRID IN RADIANS
C
C CALLED BY: ANLYTC, ERRANL, FORCE, INIT, INPUT
C CALLS:
C
C REVISIONS:
C 7-13-92 CHANGE TO CCM CODING CONVENTIONS (R. JAKOB)
C
C
C---- Model Parameters -------------------------------------------------
C
      INCLUDE 'params.i'
C
C------------ Arguments ------------------------------------------------
C
C     Input
C
C GRID INDEX (LONGITUDE)
      INTEGER I
C
C------ Local Parameters -----------------------------------------------
C
      REAL PI
      PARAMETER (PI=3.141592653589793)
C
C----- Executable Statements -------------------------------------------
C
      GLON = (2.0*PI)/REAL(NLON)*REAL(I-1)
      RETURN
C
      END
