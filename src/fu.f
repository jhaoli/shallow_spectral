      REAL FUNCTION FU(RLATD)
C
C SPECIFIES RIGHT HAND SIDE OF BALANCE EQUATION FOR PURPOSES OF             
C NUMERICAL INTEGRATION; REQUIRES REAL FUNCTION US; 
C USED FOR TEST CASE 3 GEOPOTENTIAL FIELD IN SUBROUTINE INIT
C                                                                             
C CALLED BY: D01AHE
C CALLS: US
C
C REVISIONS:
C 7-10-92 CHANGE TO CCM CODING CONVENTIONS (R. JAKOB)
C
C---- Common Blocks ----------------------------------------------------
C
C Constants & Timesteps
      INCLUDE 'consts.i'
C
C------------ Arguments ------------------------------------------------
C
C     Input
C
C LATITUDE
      REAL RLATD
C
C----- External Functions ----------------------------------------------
C
C Zonal flow
C
      EXTERNAL US
      REAL US
C
C----- Executable Statements -------------------------------------------
C
      FU = (2.0*OMEGA*A*SIN(RLATD) + TAN(RLATD)*US(RLATD))*US(RLATD)            
      RETURN                                                                    
      END                                                                       
