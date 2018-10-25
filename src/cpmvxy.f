C##################################################################
C#
C#    SUBROUTINE CPMVXY(IMAP,XX,YY,UU,VV,SFX,SFY,MX,MY,MXF,MYF)
C#
C#    REPLACES THE TWO FUNCTIONS 
C#       MXF(XX,YY,UU,VV,SFX,SFY,MX,MY)  AND
C#       MYF(XX,YY,UU,VV,SFX,SFY,MX,MY)
C#       FUNCTION VALUES ARE REPLACED BY LAST TWO ARGUMENTS
C#
C#    DUMMY ARGUMENTS
C#
C#    INTEGER IMAP = 3 USER MAPPING FOR GAUSSIAN LATITUDES 
C#    REAL XX,YY: POSITION OF VECTOR
C#    REAL UU,VV: LENGTH OF VECTOR
C#    INTEGER MX,MY: VECTOR TAIL POSITION (PLOTTER COORDINATES)
C#    REAL SFX,SFY: SCALING PARAMETERS
C#    INTEGER MXF,MYF: VECTOR HEAD POSITION (PLOTTER COORDINATES)
C#
C#    CALLED BY: VELVCT
C#    CALLS: MAPTRN
C#
      SUBROUTINE CPMVXY(IMAP,XX,YY,UU,VV,SFX,SFY,MX,MY,MXF,MYF)
C
C------------------------------------------------------------------
C
C     DEFINE PARAMETERS
      INCLUDE 'params.i'
C     MAPPING DATA
      INCLUDE 'complt.i'
C
C     PARAMETER TYPES
C
      REAL XX,YY,UU,VV,SFX,SFY
      INTEGER IMAP,MX,MY,MXF,MYF
C
C     TEMPORARIES
C
      REAL X1,X2,Y1,Y2,VLEN,VELO,U,V,
     &     XLON,YLAT,XFRACT,YFRACT,CFCT
      INTEGER M,N
C-------------------------------------------------------------------
      IF (IMAP .EQ. 3) THEN
         M = INT(XX)
         N = INT(YY)
         XFRACT = XX-REAL(M)
         IF (XFRACT .GT. 0.0) THEN
            XLON = RLOND(M) + (RLOND(M+1)-RLOND(M))*XFRACT
         ELSE
            XLON = RLOND(M)
         ENDIF
         YFRACT = YY-REAL(N)
         IF (YFRACT .GT. 0) THEN
            YLAT = RLATD(N) + (RLATD(N+1)-RLATD(N))*YFRACT
         ELSE
            YLAT = RLATD(N)
         ENDIF
         IF (.NOT. LCP) THEN
            CFCT = COS(.017453292519943 * YLAT)
         ELSE
            CFCT = 1.0
         ENDIF
         CALL MAPTRN(YLAT,XLON,X1,Y1)
         IF (CFCT .NE. 0.0) THEN
            CALL MAPTRN(YLAT+VVSF*VV,XLON+VVSF*UU/CFCT,X2,Y2)
         ELSE
            CALL MAPTRN(YLAT+VVSF*VV,XLON,X2,Y2)
         ENDIF
         VLEN = SQRT((X2-X1)**2+(Y2-Y1)**2)
         VELO = SQRT(UU**2+VV**2)
         IF (VLEN .EQ. 0.0) THEN
            U = 0.0
            V = 0.0
         ELSE
            U = ((X2-X1)/VLEN)*VELO
            V = ((Y2-Y1)/VLEN)*VELO
         ENDIF
         MXF = MX + IFIX(SFX*U)
         MYF = MY + IFIX(SFY*V)
      ELSE
         WRITE(0,100) IMAP
  100    FORMAT(/,' STSWM: FATAL ERROR IN CPMVXY:',/,
     &      ' ILLEGAL MAPPING PARAMETER IMAP = ',I2)
         STOP
      END IF
      RETURN
      END
