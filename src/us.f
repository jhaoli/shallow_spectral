      REAL FUNCTION US(RLATD)  
C                                                                              
C SPECIFIES ZONALLY SYMMETRIC U FIELD AS A FUNCTION OF LATITUDE             
C (BUMP FUNCTION).  REQUIRES REAL FUNCTION BF2                              
C USED FOR TEST CASE 3 (SEE PAPER BY BROWNING ET AL.)
C
C CALLED BY: FU, INIT
C CALLS: BF2
C
C REVISIONS:
C 7-10-92 CHANGE TO CCM CODING CONVENTIONS (R. JAKOB)
C
C---- Common Blocks ----------------------------------------------------
C
C Constants
      INCLUDE 'consts.i'
C
C------------ Arguments ------------------------------------------------
C
C     Input
C
C LATITUDE
      REAL RLATD
C
C------ Local Parameters -----------------------------------------------
C
      REAL PI
      PARAMETER (PI=3.141592653589793)
C
C------ Local Variables ------------------------------------------------
C
C     NORTH-SOUTH EXTENT OF FLOW FIELD
      REAL RLATE,RLATB
C     FLOW PROFILE TEMPORARIES
      REAL XE,X
C     MAX AMPLITUDE OF FLOW (12 DAY ROTATION SPEED)
      REAL UBAR
C
C----- External Functions ----------------------------------------------
C
C     AUXILIARY FUNCTION (SEE BELOW)
C
      EXTERNAL BF2
      REAL BF2
C
C----- Executable Statements -------------------------------------------
C
      UBAR  = (2.0*PI*A)/(12.0*24.0*3600.0)
      RLATB = -PI/6.                                                            
      RLATE =  PI/2.                                                            
      XE    = 3.0E-1                                                            
      X     = XE*(RLATD-RLATB)/(RLATE-RLATB)                                    
      US    = UBAR*BF2(X)*BF2(XE-X)*EXP(4.0/XE)
      RETURN                                                                    
      END                                                                       
