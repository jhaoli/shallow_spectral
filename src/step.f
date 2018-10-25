      SUBROUTINE STEP                                                           
C
C THIS ROUTINE CONTROLS THE TIMESTEPPING PROCESS, AND CALLS ROUTINE
C COMP1 TO COMPUTE A SINGLE TIMESTEP. IT CHANGES THE INDEXES LNM1,
C LN AND LNP1 AS POINTERS INTO A CIRCULAR BUFFER OF TIMELEVELS OF
C THE PROGNOSTIC VARIABLES DIVERGENCE, VORTICITY AND GEOPOTENTIAL.
C INITIALIZATION IS BY TWO SEMI-IMPLICIT TIMESTEPS, NORMAL TIME-
C STEPPING BY A CENTERED LEAPFROG-SCHEME. AN ASSELIN FILTER CAN BE
C USED TO PREVENT MODAL SPLITTING.
C                                                                              
C CALLED BY: STSWM	
C CALLS: COMP1
C
C REVISIONS:
C 7-10-92 CHANGE TO CCM CODING CONVENTIONS (R. JAKOB)
C
C---- Model Parameters -------------------------------------------------
C
      INCLUDE 'params.i'
C
C---- Common Blocks ----------------------------------------------------
C
C Constants
      INCLUDE 'consts.i'
C Prognostic variables
      INCLUDE 'tdvars.i'
C
C------ Local Variables ------------------------------------------------
C
C    LATITUDE/LONGITUDE INDEX
C
      INTEGER I,J
C
C----- Executable Statements -------------------------------------------
C
C     MAIN LOOP OF TIME INTEGRATION                                             
C     CENTERED (LEAPFROG) PROCEDURE                                          
C                                                                               
      IF (NSTEP .EQ. 1) THEN                                    
C                                                                               
C        PERFORM INITIAL TIME STEPS BY COPYING PROGNOSTIC
C        VARIABLES INTO PREVIOUS TIMELEVEL FOR SMOOTH START
C
C        FIRST FORWARD TIMESTEP (DT/2.0)
C
         DO 70 J=1,NLAT                                                      
            DO 60 I=1,NLON
               ZETA(I,J,LVLS) = ZETA(I,J,1)
               DIV (I,J,LVLS) = DIV (I,J,1)
               PHI (I,J,LVLS) = PHI (I,J,1) 
   60       CONTINUE
   70    CONTINUE                                                            
C
         LNM1   = LVLS
         LN     = 1
         LNP1   = 2
         CALL COMP1(DT/2.0)
C
C        SECOND CENTERED TIMESTEP (DT)
C
C        DO 90 J=1,NLAT
C           DO 80 I=1,NLON
C              ZETA(I,J,LVLS) = ZETA(I,J,1)
C              DIV (I,J,LVLS) = DIV (I,J,1)
C              PHI (I,J,LVLS) = PHI (I,J,1)
C  80       CONTINUE
C  90    CONTINUE

         LNM1 = LVLS
         LN = 2
         LNP1 = 1
         CALL COMP1(DT)
C
C        SET BUFFER POINTERS FOR FIRST LEAPFROG TIMESTEP (2.0*DT)
C
         LNM1 = LVLS
         LN   = 1
         LNP1 = 2
C        
         RETURN                                                              
C
      ELSE
C                                                                               
C     ALL OTHER INSTANCES OF CENTERED DIFFERENCING                           
C                                                                               
         CALL COMP1(2.0*DT)
C                                                                               
C        ASSELIN FILTER (FILTERS COMPUTATIONAL MODE BY
C        DAMPING HIGH-FREQUENCY MODES)
C                                                                               
         IF (AFC .NE. 0.0) THEN
            DO 100 J=1,NLAT 
               DO 95 I=1,NLON
                  ZETA(I,J,LN) = ZETA(I,J,LN) + AFC*(ZETA(I,J,LNM1)
     $                          - 2.0*ZETA(I,J,LN) + ZETA(I,J,LNP1))            
                  DIV(I,J,LN)  = DIV(I,J,LN) + AFC*(DIV(I,J,LNM1)
     $                          - 2.0*DIV(I,J,LN) + DIV(I,J,LNP1))              
                  PHI(I,J,LN)  = PHI(I,J,LN) + AFC*(PHI(I,J,LNM1) 
     $                          - 2.0*PHI(I,J,LN) + PHI(I,J,LNP1))              
   95          CONTINUE
  100       CONTINUE   
         ENDIF
C                                                                               
C        ADVANCE CIRCULAR BUFFER POINTERS
C
         LNM1 = LN                                                              
         LN   = LNP1                                                            
         LNP1 = MOD(LN+1,LVLS)
         IF (LNP1 .EQ. 0) LNP1 = LVLS
C                                                                               
         RETURN                                                                 
C                                                                               
      ENDIF
C
      END                                                                       
