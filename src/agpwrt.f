      SUBROUTINE AGPWRT (XPOS, YPOS, CHRS, NCHS, ISIZ, IORI, ICEN)
C
C     HIGH QUALITY CHARACTER PLOT LABELS
C
      CHARACTER*(*) CHRS
      REAL XPOS, YPOS, CSFL
      INTEGER NCHS, ISIZ, IORI, ICEN
      CALL PCGETR ('CS - CONSTANT SPACING FLAG', CSFL)
      IF (ICEN.NE.0) THEN
        CALL PCSETR ('CS - CONSTANT SPACING FLAG', 1.25)
      ELSE
        CALL PCSETR ('CS - CONSTANT SPACING FLAG', 0.0 )
      ENDIF
      CALL PLCHHQ (XPOS, YPOS, CHRS(1:NCHS),
     +             .8*REAL(ISIZ), REAL(IORI), REAL(ICEN))
      CALL PCSETR ('CS - CONSTANT SPACING FLAG', CSFL)
C                                                                       
      RETURN                                                            
      END   
