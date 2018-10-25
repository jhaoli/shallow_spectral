      REAL FUNCTION DBUBF(RLAT)
C
C SPECIFIES FUNCTIONAL FORM OF 1ST DERIVATIVE OF BIG U BAR (BUBFNC)
C
C CALLED BY: ANLYTC, FORCE
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
C------ Local Variables ------------------------------------------------
C
      REAL RMU,COSLAT
C
C----- Executable Statements -------------------------------------------
C
      RMU    = SIN(RLAT)
      COSLAT = COS(RLAT)
      DBUBF  = 2.0*SU0*(2.0*RMU*COSLAT)**(NPWR-1)*
     $         (REAL(NPWR)-REAL(2*NPWR+1)*RMU*RMU)
C                                                                               
      RETURN                                                                    
      END                                                                       
