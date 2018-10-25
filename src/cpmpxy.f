C The contents of this file are copyright protected:
C --------------------------------------------------
C NCAR Graphics - UNIX Version 3.1.2
C Copyright (C) 1987, 1988, 1989, 1991
C University Corporation for Atmospheric Research
C All Rights Reserved
C
C Modified for high non-equidistant grids as part
C of the shallow water model STSWM by R. Jakob
C
C#################################################################
C#
C#    SUBROUTINE CPMPXY(IMAP,XINP,YINP,XOTP,YOTP)
C#
C#    THIS ROUTINE PROVIDES A MAPPING FOR NON-EQUALLY SPACED
C#    DATA IN THE CONPACK PACKAGE OF NCAR GRAPHICS. THE
C#    ROUTINE IS DESCRIBED IN VERSION 3.00, PAGES 3-55, 3-56.
C#    SEE ALSO MANUAL PAGES 3-15, 3-16 FOR CORRECT USE
C#
C#    DUMMY ARGUMENTS
C#
C#    INTEGER   IMAP = 0 NO MAPPING, CPMPXY NOT CALLED
C#              IMAP = 1 STANDARD SCALED MAPPING
C#              IMAP = 2 POLAR COORDINATE SYSTEM
C#              IMAP = 3 USER-DEFINED GRID MAPPING (GAUSSIAN)
C#              IMAP IS DEFINED BY THE VALUE OF 'MAP'
C#    REAL XINP,YINP: INPUT DATA
C#    REAL XOTP,YOTP: OUTPUT DATA
C#
C#    CALLED BY: CPCNRC, PLOTS, VELVCT
C#    CALLS: MAPTRN
C#
      SUBROUTINE CPMPXY(IMAP,XINP,YINP,XOTP,YOTP)
C
C*****************************************************************
C
C     DEFINE PARAMETERS
      INCLUDE 'params.i'
C     MAPPING DATA
      INCLUDE 'complt.i'
C
C     PARAMETER TYPES
C
      REAL XINP,YINP,XOTP,YOTP,XLON,YLAT
      INTEGER IMAP
C
C------ Local Parameters -----------------------------------------------
C
      REAL PI
      PARAMETER (PI=3.141592653589793)
C
C     HANDLE INTERMAEDIATE POINTS
C
      INTEGER M,N
      REAL XFRACT, YFRACT, TEST, PLAT, PLON
C
C*********************************************************************
C
      IF (IMAP .EQ. 1) THEN
         CALL MAPTRN(YINP,XINP,XOTP,YOTP)
      ELSE IF (IMAP .EQ. 2) THEN
         XOTP = XINP * COS(PI/180.0 * YINP)
         YOTP = XINP * SIN(PI/180.0 * YINP)
      ELSE IF (IMAP .EQ. 3) THEN
         M = INT(XINP)
         N = INT(YINP)
C
C        COMPUTE LONGITUDE
C
         XFRACT = XINP-REAL(M)
         IF (XFRACT .GT. 0.0) THEN
            XLON = RLOND(M) + (RLOND(M+1)-RLOND(M))*XFRACT
         ELSE
            XLON = RLOND(M)
         ENDIF
C
C        COMPUTE LATITUDE
C
         YFRACT = YINP-REAL(N)
         IF (YFRACT .GT. 0.0) THEN
            YLAT = RLATD(N) + (RLATD(N+1)-RLATD(N))*YFRACT
         ELSE
            YLAT = RLATD(N)
         ENDIF
         CALL MAPTRN(YLAT,XLON,XOTP,YOTP)
C
C        OUT-OF-RANGE VALUES FOR POLAR STEREOGRAPHIC PROJECTION
C
         IF (LPSP) THEN
            IF (LNORTH) THEN
C
C              NORTHERN HEMISPHERE
C
               IF (POROT .EQ. 0.0) THEN
                  TEST = YLAT
               ELSE
                  PLAT = YLAT/180.0*PI
                  PLON = XLON/180.0*PI
                  TEST = SIN(PI/2.0-POROT)*SIN(PLAT) +
     &               COS(PI/2.0-POROT)*COS(PLAT)*COS(-PI-PLON)
               ENDIF
               IF (TEST .LT. 0.0) THEN
                  XOTP = 1.0E12
               ENDIF
            ELSE
C
C              SOUTHERN HEMISPHERE
C
               IF (POROT .EQ. 0.0) THEN
                  TEST = -YLAT
               ELSE
                  PLAT = YLAT/180.0*PI
                  PLON = XLON/180.0*PI
                  TEST = SIN(-PI/2.0+POROT)*SIN(PLAT) +
     &               COS(-PI/2.0+POROT)*COS(PLAT)*COS(-PLON)
               ENDIF
               IF (TEST .LT. 0.0) THEN
                  XOTP = 1.0E12
               ENDIF
            ENDIF
         ENDIF
      ELSE
         XOTP = XINP
         YOTP = YINP
      END IF
      RETURN
      END

