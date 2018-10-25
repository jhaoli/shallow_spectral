      SUBROUTINE ROTATE(RLON,RLAT,ALPHA,ROTLON,ROTLAT)
C                                                                              
C THIS SUBROUTINE COMPUTES THE ROTATED COORDINATES
C ROTLON,ROTLAT FOR A ROTATION BY ANGLE ALPHA
C                                                                              
C CALLED BY: ANLYTC, INIT
C CALLS:
C
C REVISIONS:
C 7-13-92 CHANGE TO CCM CODING CONVENTIONS (R. JAKOB)
C
C
C------------ Arguments ------------------------------------------------
C
C     Input
C
C ORIGINAL LONGITUDE
      REAL RLON
C ORIGINAL LATITUDE
      REAL RLAT
C ROTATION ANGLE
      REAL ALPHA
C
C     Output
C
C ROTATED LONGITUDE
      REAL ROTLON
C ROTATED LATITUDE
      REAL ROTLAT
C
C------ Local Parameters -----------------------------------------------
C
      REAL PI
      PARAMETER (PI=3.141592653589793)
C 
C------ Local Variables ------------------------------------------------
C
      REAL TEST
C   
C----- Executable Statements -------------------------------------------
C                                                                               
      IF (ALPHA .EQ. 0.0) THEN
C
C        NO ROTATION
C
         ROTLON = RLON
         ROTLAT = RLAT
      ELSE
C
C        ROTATION BY ANGLE ALPHA
C
C        ROTATED LATITUDE
C
         TEST = SIN(RLAT)*COS(ALPHA)-
     $            COS(RLAT)*COS(RLON)*SIN(ALPHA)
         IF (TEST .GT. 1.0) THEN
            ROTLAT = PI/2.0
         ELSEIF (TEST .LT. -1.0) THEN
            ROTLAT = -PI/2.0
         ELSE
            ROTLAT = ASIN(TEST)
         ENDIF
C
C        ROTATED LONGITUDE
C
         TEST = COS(ROTLAT)
         IF (TEST .EQ. 0.0) THEN
            ROTLON = 0.0
         ELSE
            TEST = SIN(RLON)*COS(RLAT)/TEST
            IF (TEST .GT. 1.0) THEN
               ROTLON = PI/2.0
            ELSEIF (TEST .LT. -1.0) THEN
               ROTLON = -PI/2.0
            ELSE
               ROTLON = ASIN(TEST)
            ENDIF
         ENDIF
C
C        ADJUST FOR CORRECT BRANCH OF INVERSE SINE
C
         TEST = COS(ALPHA)*COS(RLON)*COS(RLAT)
     $          + SIN(ALPHA)*SIN(RLAT)
         IF (TEST .LT. 0.0) THEN
            ROTLON = PI - ROTLON
         ENDIF
      ENDIF
C
      RETURN
C
      END
