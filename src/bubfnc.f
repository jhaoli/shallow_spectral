      REAL FUNCTION BUBFNC(RLAT)
C
C SPECIFIES FUNCTIONAL FORM OF U BAR
C (USED FOR FORCED SOLUTION ICOND = 4)
C SEE PAPER BY BROWNING ET. AL.:
C A COMPARISON OF THREE NUMERICAL METHODS ...
C (EQUATION 5.10a)*COS(RLAT)
C
C CALLED BY: ANLYTC, ERRANL, FORCE, FUNC2
C CALLS:
C
C REVISIONS:
C 7-10-92 CHANGE TO CCM CODING CONVENTIONS (R. JAKOB)
C
C---- Common Blocks ----------------------------------------------------
C
C     Case 4 variables
      INCLUDE 'case4.i'
C
C------------ Arguments ------------------------------------------------
C
C     Input
C
C LATITUDE
      REAL RLAT
C
C----- Executable Statements -------------------------------------------
C
      BUBFNC = SU0*(2.0*SIN(RLAT)*COS(RLAT))**NPWR
C                                                                               
      RETURN                                                                    
      END                                                                       
