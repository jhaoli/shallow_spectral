      SUBROUTINE SHADER (XCS,YCS,NCS,IAI,IAG,NAI)
C
C This version of SHADER shades the polygon whose edge is defined by
C the points ((XCS(I),YCS(I)),I=1,NCS) if and only, relative to edge
C group 3, its area identifier is a 1.  The package SOFTFILL is used
C to do the shading.
C
        INTEGER IAI, IAG, NCS, NAI
        REAL XCS, YCS
        DIMENSION XCS(*),YCS(*),IAI(*),IAG(*)
C
C Define workspaces for the shading routine.
C
        INTEGER NWORK
        PARAMETER (NWORK=10000)
        INTEGER IND
        REAL DST
        DIMENSION DST(NWORK),IND(NWORK)
C
C Temporaries
C
        INTEGER ISH, I 
C
C Turn off shading.
C
        ISH=0
C
C If the area identifier for group 3 is a 1, turn on shading.
C
        DO 101 I=1,NAI
          IF (IAG(I).EQ.3.AND.IAI(I).EQ.1) ISH=1
  101   CONTINUE

C
C If shading is turned on, shade the area.  The last point of the
C edge is redundant and may be omitted.
C
        IF (ISH.NE.0) THEN
          CALL SFSETI ('ANGLE',45)
          CALL SFSETR ('SPACING',0.01)
          CALL SFWRLD (XCS,YCS,NCS-1,DST,NWORK,IND,NWORK)
        END IF
C
C Done.
C
        RETURN
C
      END

