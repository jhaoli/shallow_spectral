      REAL FUNCTION GLAT(J)
C
C THIS FUNCTION RETURNS THE LATITUDE OF THE GAUSSIAN
C GRID IN RADIANS
C
C CALLED BY: ANLYTC, CALP, COMP1, ERRANL, FORCE, INIT, INPUT,
C            NRGTCS, OUTPTP, PLOTS, STSWM, ZD
C CALLS:
C
C REVISIONS:
C 7-13-92 CHANGE TO CCM CODING CONVENTIONS (R. JAKOB)
C
C---- Model Parameters -------------------------------------------------
C
      INCLUDE 'params.i'
C
C---- Common Blocks ----------------------------------------------------
C
C TRANSFORM ARRAYS
      INCLUDE 'trnsfm.i'
C
C------------ Arguments ------------------------------------------------
C
C     Input
C
C GRID INDEX (LATITUDE)
      INTEGER J
C
C----- Executable Statements -------------------------------------------
C
      GLAT = THTA(J)
      RETURN
C
      END
