      REAL FUNCTION BF2(X)
C
C USED IN CONJUNCTION WITH FINCTION US
C
C CALLED BY: US
C CALLS:
C
C REVISIONS:
C 7-10-92 CHANGE TO CCM CODING CONVENTIONS (R. JAKOB)
C
C------------ Arguments ------------------------------------------------
C
C     Input
C
C Scaled LATITUDE
      REAL X
C
C----- Executable Statements -------------------------------------------
C
      IF (X .LE. 0.0) THEN  
         BF2 = 0.0                                                              
      ELSE                                                                      
         BF2 = EXP(-1.0/X)                                                      
      ENDIF                                                                     
C                                                                               
      RETURN                                                                    
      END                                                                       
