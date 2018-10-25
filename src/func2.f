      REAL FUNCTION FUNC2(THETA)  
C
C SPECIFIES RIGHT HAND SIDE OF NON-LINEAR BALANCE EQUATION FOR              
C PURPOSES OF NUMERICAL INTEGRATION; REQUIRES REAL FUNCTION BUBFNC          
C MAKE SURE THAT CONSTANTS ARE CONSISTANT WITH THOSE IN INPUT!!             
C USED FOR TEST CASE 4 IN ROUTINE INIT
C                                                                              
C CALLED BY: D01AHE
C CALLS: BUBFNC
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
      REAL THETA
C
C------ Local Variables ------------------------------------------------
C
      REAL SU
C
C----- External Functions ----------------------------------------------
C
C     ZONAL FLOW FUNCTION
C
      EXTERNAL BUBFNC
      REAL BUBFNC
C
C----- Executable Statements -------------------------------------------
C
      SU     = BUBFNC(THETA)
      FUNC2  = (2.0*OMEGA*A*SIN(THETA) + TAN(THETA)*SU)*SU                      
C                                                                               
      RETURN                                                                    
      END                                                                       
