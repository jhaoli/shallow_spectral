      REAL FUNCTION D2BUBF(RLAT)
C
C SPECIFIES FUNCTIONAL FORM OF 2ND DERIVATIVE OF BIG U BAR (BUBFNC)
C (USED FOR FORCED CASE ICOND = 4)
C
C CALLED BY: FORCE
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
      D2BUBF = 8.0*SU0*(2.0*RMU*COSLAT)**(NPWR-3)*RMU*
     $         (REAL((NPWR-1)*NPWR)+RMU*RMU*(REAL(NPWR-1)-
     $         REAL(2*NPWR*(2*NPWR+1))*COSLAT*COSLAT))
C                                                                               
      RETURN                                                                    
      END                                                                       
