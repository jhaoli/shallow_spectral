C#######################################################################
C#
C#    SUBROUTINE ORTHO
C#                                                                              
C#    CHECKS ORTHOGONALITY OF BASISFUNCTIONS
C#                                                                              
C#    CALLED BY: STSWM
C#
      SUBROUTINE ORTHO
C                                                                               
C***********************************************************************        
C     DEFINE PARAMETERS
      INCLUDE 'params.i'
C***********************************************************************        
C                                                                               
C     GLOBAL QUANTITIES (COMMON BLOCKS)                                         
C                                                                               
C***********************************************************************        
C     TRANSFORM ARRAYS
      INCLUDE 'trnsfm.i'
C                                                                               
C***********************************************************************        
C                                                                               
C     DUMMY STORAGE (SUBROUTINE SHTRNS)                     
C
C     ORTHONORMALITY/ORTHOGONALITY CHECKSUMS
C
      REAL ALPON(NALP),ALPOG(NALP),TMP(0:NN)
      REAL PLOTON(0:MM,0:NN),PLOTOG(0:MM,0:NN)
C
C     TEMPORARIES
C
      INTEGER NL, JM, JN, IS, ISO
      REAL COEFF1,COEFF2, VMAX, VMIN
C                                                                               
C-----------------------------------------------------------------------        
C                                                                               
C     INITIALIZE ORTHONORMAL ARRAY
C
      DO 100 IS = 1, NALP
         ALPON(IS) = 0.0
         ALPOG(IS) = 0.0
  100 CONTINUE
C
C     COMPUTE ORTHONORMAL SUM
C
      DO 300 NL = 1, NLAT/2
         DO 200 IS = 1, NALP
            ALPON(IS) = ALPON(IS) + ALP(IS,NL)*ALP(IS,NL)*WTS(NL)
  200    CONTINUE 
  300 CONTINUE
C
      DO 390 JN = 0, NN
         IS = LDIAG(JN,2)
         DO 389 JM = 1,LDIAG(JN,1)
         IF (MOD(JN,2) .EQ. 0) THEN
            COEFF1 = 1.0
         ELSE
            COEFF1 = -1.0
         ENDIF
C
C        INITIALIZE
C
         DO 330 ISO = 0,JN-1
            TMP(ISO) = 0.0
  330    CONTINUE
C
C        INTEGRATE
C
         DO 380 NL = 1, NLAT/2
            DO 350 ISO = 0,JN-1
               IF (MOD(ISO,2) .EQ. 0) THEN
                  COEFF2 = 1.0
               ELSE
                  COEFF2 = -1.0
               ENDIF
               TMP(ISO) = TMP(ISO) + 
     &         (1.0+COEFF1*COEFF2)*ALP(LDIAG(ISO,2)+JM,NL)
     &                    *ALP(IS+JM,NL)* WTS(NL)
  350       CONTINUE
  380    CONTINUE
C
C        GET MAXIMUM
C
         VMAX = 0.0
         VMIN = 0.0
         DO 385 ISO = 0, JN-1
            VMIN = MIN(TMP(ISO),VMIN)
            VMAX = MAX(TMP(ISO),VMAX)
  385    CONTINUE
         IF (ABS(VMIN) .GT. ABS(VMAX)) THEN
            ALPOG(IS+JM) = VMIN
         ELSE
            ALPOG(IS+JM) = VMAX
         ENDIF
  389    CONTINUE
  390 CONTINUE
C
C     NORMALIZE
C
      DO 400 IS = 1, NALP
         ALPON(IS) = 2.0*ALPON(IS)
  400 CONTINUE
C
C     PLOT FIELDS & OUTPUT VALUES
C
      DO 600 JN=0,NN
         DO 550 JM=0,MM
            PLOTON(JM,JN) = 0.0
            PLOTOG(JM,JN) = 0.0
  550    CONTINUE
  600 CONTINUE

      DO 500 JN = 0,NN
         IS = LDIAG(JN,2)
         WRITE (6,999) (JN+JM-1,JM-1,ALPON(IS+JM)-1.0,ALPOG(IS+JM),
     &                  JM=1,LDIAG(JN,1))
         DO 525 JM=1,LDIAG(JN,1)
            PLOTON(JM-1,JN+JM-1) = ALPON(IS+JM)-1.0
            PLOTOG(JM-1,JN+JM-1) = ALPOG(IS+JM)
  525    CONTINUE
  500 CONTINUE
C
      CALL CPCNRC(PLOTON,MM+1,MM+1,NN+1,0.0,0)
      CALL LABTOP('Orthonormality Error',0.015)
      CALL FRAME
      CALL CPCNRC(PLOTOG,MM+1,MM+1,NN+1,0.0,0)
      CALL LABTOP('Orthogonality Error',0.015)
      CALL FRAME
C
      RETURN
C
  999 FORMAT('N = ',I4,' M = ',I4,' ALPON-1.0 = ',F18.16,
     &       ' ALPOG = ',F18.16)
      END
