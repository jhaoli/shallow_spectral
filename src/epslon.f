      REAL FUNCTION EPSLON (X)
C
C ESTIMATE MACHINE ARITHMETIC ROUNDOFF IN QUANTITIES OF SIZE X
C
C     THIS PROGRAM SHOULD FUNCTION PROPERLY ON ALL SYSTEMS
C     SATISFYING THE FOLLOWING TWO ASSUMPTIONS,
C        1.  THE BASE USED IN REPRESENTING FLOATING POINT
C            NUMBERS IS NOT A POWER OF THREE.
C        2.  THE QUANTITY  A  IN STATEMENT 10 IS REPRESENTED TO
C            THE ACCURACY USED IN FLOATING POINT VARIABLES
C            THAT ARE STORED IN MEMORY.
C     THE STATEMENT NUMBER 10 AND THE GO TO 10 ARE INTENDED TO
C     FORCE OPTIMIZING COMPILERS TO GENERATE CODE SATISFYING
C     ASSUMPTION 2.
C     UNDER THESE ASSUMPTIONS, IT SHOULD BE TRUE THAT,
C            A  IS NOT EXACTLY EQUAL TO FOUR-THIRDS,
C            B  HAS A ZERO FOR ITS LAST BIT OR DIGIT,
C            C  IS NOT EXACTLY EQUAL TO ONE,
C            EPS  MEASURES THE SEPARATION OF 1.0 FROM
C                 THE NEXT LARGER FLOATING POINT NUMBER.
C
C CALLED BY: INPUT
C CALLS:
C
C REVISIONS:
C 7-10-92 CHANGE TO CCM CODING CONVENTIONS (R. JAKOB)
C
C------------ Arguments ------------------------------------------------
C
C     Input
C
C APPROXIMATED VALUE
      REAL X                      
C                           
C------ Local Variables ------------------------------------------------
C
      REAL A,B,C,EPS
C
C----- Executable Statements -------------------------------------------
C
      A = 4.0E0/3.0E0                 
   10 B = A - 1.0E0                 
      C = B + B + B              
      EPS = ABS(C-1.0E0)             
      IF (EPS .EQ. 0.0E0) GO TO 10     
      EPSLON = EPS*ABS(X)                   
      RETURN                    
      END                     
