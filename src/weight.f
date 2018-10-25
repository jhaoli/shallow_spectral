      REAL FUNCTION WEIGHT(I,J)
C
C THIS FUNCTION RETURNS A WEIGHT FOR NUMERIC INTEGRATION ON THE
C DEFINED GRID. IT IS NORMALIZED SUCH THAT THE SUM OVER ALL
C LATITUDE AND LONGITUDE POINTS IS 1.0
C
C#    DUMMY ARGUMENTS
C#
C#    INTEGER I,J: LATITUDE, LONGITUDE INDEX
C#
C#    CALLED BY: NRGTCS
C#    CALLS:
C#
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
C LONGITUDE INDEX
      INTEGER I
C LATITUDE INDEX
      INTEGER J
C
C----- Executable Statements -------------------------------------------
C
      WEIGHT = WTS(J)/(2.0*REAL(NLON))
      RETURN
C
      END
