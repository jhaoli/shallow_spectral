C The contents of this file are copyright protected:
C --------------------------------------------------
C NCAR Graphics - UNIX Version 3.1.2
C Copyright (C) 1987, 1988, 1989, 1991
C University Corporation for Atmospheric Research
C All Rights Reserved
C
C Modified for high resolution plots as part
C of the shallow water model STSWM by R. Jakob
C
C-----------------------------------------------------------------------
C
      SUBROUTINE CPCNRC (ZDAT,KZDT,MZDT,NZDT,FINC,NDSH)
C
      INTEGER KZDT,MZDT,NZDT,NDSH
      REAL ZDAT, FINC
      DIMENSION ZDAT(KZDT,*)
C
C This routine simulates the old routine CONREC.
C
C Define some needed dimensions.
C
      INTEGER LRWK, LIWK, LAWK
      PARAMETER (LRWK=4000,LIWK=40000,LAWK=300000)
      INTEGER LIRA, LCRA
      PARAMETER (LIRA=1000,LCRA=30000)
C
C     Flag to turn on labels
C
      LOGICAL PLABEL
      PARAMETER (PLABEL = .FALSE.)
C
C Declare arrays needed by ARSCAM
C
      INTEGER IARA,IGRA
      REAL XCRA,YCRA
      DIMENSION IARA(LIRA),IGRA(LIRA)
      DIMENSION XCRA(LCRA),YCRA(LCRA)
C
C Define required workspace arrays.
C
      REAL RWRK
      INTEGER IWRK,IAMA
      DIMENSION RWRK(LRWK),IWRK(LIWK),IAMA(LAWK)
C
C Contour level value
C
      REAL CLEV
C
C Temporaries
C
      INTEGER IDSH,ICFF,IWU,IRWU,NCLO,ICLO,ICLU
C
C     ROUTINE FOR DRAWING CONTOUR LINES
C
      EXTERNAL DRAWCL
C
C     ROUTINE FOR SHADING (NEGATIVE VALUES)
C
      EXTERNAL SHADER
C
C     BIT OPERATIONS
C
      INTEGER IAND, IOR, ISHIFT
      EXTERNAL IAND, IOR, ISHIFT
C
C GKS aspect source flags
C
      INTEGER IASF
      DIMENSION IASF(13)
      DATA IASF / 13*1 /
C
C-----------------------------------------------------------------------
C
C Turn off clipping, so that the informational label won't disappear.
C
      CALL GSCLIP (0)
C
C     Set GKS aspect source flags to 'individual'
C
      CALL GSASF (IASF)
C
C Arrange for the selection of contour levels as desired by the user.
C
      IF (FINC.LT.0.) THEN
        CALL CPSETI ('CLS - CONTOUR LEVEL SELECTOR',MAX(1,INT(-FINC)))
        CALL CPSETR ('CIS - CONTOUR INTERVAL SPECIFIER',0.)
      ELSE IF (FINC.EQ.0.) THEN
        CALL CPSETI ('CLS - CONTOUR LEVEL SELECTOR',16)
        CALL CPSETR ('CIS - CONTOUR INTERVAL SPECIFIER',0.)
      ELSE
        CALL CPSETI ('CLS - CONTOUR LEVEL SELECTOR',1)
        CALL CPSETR ('CIS - CONTOUR INTERVAL SPECIFIER',FINC)
      END IF
C
C     Turn on positioning of labels by penalty scheme
C
      CALL CPSETI ('LLP - LINE LABEL POSITIONING',2)
C
C Decide what dash pattern to use.
C
      IDSH=ABS(NDSH)
      IF ((IDSH .EQ. 0) .OR. (IDSH .EQ. 1) .OR. 
     &    (IDSH .EQ. 1023)) THEN
        IDSH=IOR(ISHIFT(32767,1),1)
      ELSE
        IDSH=IOR(ISHIFT(IDSH,6),IAND(ISHIFT(IDSH,-4),63))
      END IF
C
      IF (PLABEL) THEN
C
C Label highs and lows.
C
         CALL CPSETC ('HLT - HIGH/LOW LABEL TEXT',
     &               'H:B:$ZDV$:E:''L:B:$ZDV$:E:')
         CALL CPSETC ('HIT - HIGH/LOW LABEL TEXT',
     &               'H:B:$ZDV$:E:')
         CALL CPSETC ('LOT - HIGH/LOW LABEL TEXT',
     &               'L:B:$ZDV$:E:')
C
C And put box around it
C
         CALL CPSETI('HLB - HIGH/LOW LABEL BOX FLAG',1)
         CALL CPSETI('HLO - HIGH/LOW LABEL OVERLAP FLAG',11)
      ELSE
C
C No high/low labels
C
         CALL CPSETC ('HLT - HIGH/LOW LABEL TEXT','''')
      ENDIF
C
C Initialize CONPACK and give it all array dimensions.
C
      CALL CPRECT (ZDAT,KZDT,MZDT,NZDT,RWRK,LRWK,IWRK,LIWK)
C
C If the field was constant, skip some of the following code.
C
      CALL CPGETI ('CFF - CONSTANT FIELD FLAG',ICFF)
C
      IF (ICFF.EQ.0) THEN
C
C Pick contour levels.
C
        CALL CPPKCL (ZDAT,RWRK,IWRK)
C
C     Check work space used
C
      CALL CPGETI ('IWU',IWU)
      CALL CPGETI ('RWU',IRWU)
      IF ((IWU .GT. LIWK) .OR. (IRWU .GT. LRWK)) THEN
        WRITE(0,*)'CPPKCL WORKSPACE USED/AVAILABLE'
        WRITE(0,*)' INTEGER =',IWU,LIWK
        WRITE(0,*)' REAL =', IRWU,LRWK
      ENDIF
C
C Use fat lines each fifth contour interval, mark zero contour 
C for later shading, and set dash pattern for negative contours.
C
        CALL CPGETI ('NCL - NUMBER OF CONTOUR LEVELS',NCLO)
C
        DO 10001 ICLO=1,NCLO
          CALL CPSETI ('PAI - PARAMETER ARRAY INDEX',ICLO)
          CALL CPGETI ('CLU - CONTOUR LEVEL USE FLAG',ICLU)
          IF (ICLU .EQ. 3) THEN
             CALL CPSETI ('CLL - CONTOUR LINE WIDTH',2)
          ENDIF
          CALL CPGETR ('CLV - CONTOUR LEVEL VALUE',CLEV)
          IF (CLEV .NE. 0.0) THEN
             CALL CPSETI ('AIA - AREA IDENTIFIER ABOVE LINE',0)
             CALL CPSETI ('AIB - AREA IDENTIFIER BELOW LINE',0)
          ELSE
             CALL CPSETI ('AIA - AREA IDENTIFIER ABOVE LINE',2)
             CALL CPSETI ('AIB - AREA IDENTIFIER BELOW LINE',1)
          ENDIF
          IF (NDSH .GT. 0 .OR. 
     &        ((NDSH .LT. 0.0) .AND. (CLEV .LT. 0.0))) THEN
             CALL CPSETI ('CLD - CONTOUR LINE DASH PATTERN',IDSH)
          ENDIF
10001   CONTINUE
C
      END IF
C
C Draw the contour lines, masking them if necessary.
C
C     CALL CPCLDR (ZDAT,RWRK,IWRK)
C
C
C     INITIALIZE AREA MAP
C
      CALL ARINAM(IAMA,LAWK)
C
C     PUT LABEL BOXES INTO AREA MAP
C
      CALL CPLBAM(ZDAT,RWRK,IWRK,IAMA)
C
C     Check work space used
C
      CALL CPGETI ('IWU',IWU)
      CALL CPGETI ('RWU',IRWU)
      IF ((IWU .GT. LIWK) .OR. (IRWU .GT. LRWK)) THEN
        WRITE(0,*)'CPLBAM WORKSPACE USED/AVAILABLE'
        WRITE(0,*)' INTEGER =',IWU,LIWK
        WRITE(0,*)' REAL =', IRWU,LRWK
      ENDIF
C
C     DRAW CONTOUR LINES
C
      CALL CPCLDM(ZDAT,RWRK,IWRK,IAMA,DRAWCL)
C
C     Check work space used
C
      CALL CPGETI ('IWU',IWU)
      CALL CPGETI ('RWU',IRWU)
      IF ((IWU .GT. LIWK) .OR. (IRWU .GT. LRWK)) THEN
        WRITE(0,*)'CPCLDM WORKSPACE USED/AVAILABLE'
        WRITE(0,*)' INTEGER =',IWU,LIWK
        WRITE(0,*)' REAL =', IRWU,LRWK
      ENDIF
C
C Plot labels.
C
      CALL CPLBDR (ZDAT,RWRK,IWRK)
C
C     Check work space used
C
      CALL CPGETI ('IWU',IWU)
      CALL CPGETI ('RWU',IRWU)
      IF ((IWU .GT. LIWK) .OR. (IRWU .GT. LRWK)) THEN
        WRITE(0,*)'CPLBDR WORKSPACE USED/AVAILABLE'
        WRITE(0,*)' INTEGER =',IWU,LIWK
        WRITE(0,*)' REAL =', IRWU,LRWK
      ENDIF
C
C     Add contour lines to area map
C
      CALL CPCLAM(ZDAT,RWRK,IWRK,IAMA)
C
C     Check work space used
C
      CALL CPGETI ('IWU',IWU)
      CALL CPGETI ('RWU',IRWU)
      IF ((IWU .GT. LIWK) .OR. (IRWU .GT. LRWK)) THEN
        WRITE(0,*)'CPLBDR WORKSPACE USED/AVAILABLE'
        WRITE(0,*)' INTEGER =',IWU,LIWK
        WRITE(0,*)' REAL =', IRWU,LRWK
      ENDIF
C
C     SCAN AREA MAP FOR SHADING OF NEGATIVE VALUES
C 
      CALL ARSCAM(IAMA,XCRA,YCRA,LCRA,IARA,IGRA,LIRA,SHADER)
C
      RETURN
C
      END
