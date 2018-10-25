      SUBROUTINE INFLD(CDFID,VARID,STEP,MAXH,LDIAG,MM,NN,KK,VAR) 
C                                                                              
C READS A SPECTRAL COEFFICIENT FIELD
C IN NETCDF FORMAT
C
C CALLED BY: INPTP
C CALLS:
C
C REVISIONS:
C 7-13-92 CHANGE TO CCM CODING CONVENTIONS (R. JAKOB)
C
C---- Model Parameters -------------------------------------------------
C
C NETCDF DECLARATIONS
      INCLUDE 'netcdf.inc'
C
C------------ Arguments ------------------------------------------------
C
C     Input
C
C NETCDF FILE ID
      INTEGER CDFID
C NETCDF VARIABLE ID
      INTEGER VARID
C TIMESTEP
      INTEGER STEP
C DIMENSION OF LDIAG
      INTEGER MAXH
C STRUCTURE INFORMATION FOR ARRAY VAR
      INTEGER LDIAG(0:MAXH,2)                                      
C SPECTRAL TRUNCATION PARAMETERS
      INTEGER MM,NN,KK
C   
C     Output
C
C SPECTRAL ARRAY TO BE READ
      COMPLEX VAR(*)
C
C------ Local Variables ------------------------------------------------
C
      INTEGER M,N,MDUM,NDUM
      REAL DREAL,DIMAG
C
C     NETCDF VARIABLES
C
      INTEGER IRET
      INTEGER DIMS(4)
C
C----- Statement Function ----------------------------------------------
C
C     ADDRESS COMPUTATION OF ALP, DALP, EPSIL            
C                                                                               
      INTEGER IDSP
      IDSP(MDUM,NDUM) = 1 + LDIAG(NDUM-MDUM,2)+MDUM                             
C                                                                               
C----- Executable Statements -------------------------------------------
C
      DIMS(4) = STEP
      DO 1060 N=0,MAXH                                                   
         DIMS(3) = N+1
         DO 1040 M=0,N                                                      
C                                                                               
            IF ((M .LE. MM) .AND. (N .LE. MIN0(KK,M+NN)))
     $         THEN
C                                                                               
               DIMS(2) = M+1
C
C              READ REAL AND IMAGINARY PARTS
C
               DIMS(1) = 1
               CALL NCVGT1(CDFID,VARID,DIMS,DREAL,IRET)
               IF (IRET .NE. 0) THEN
                  WRITE(0,500) M,N
  500             FORMAT(/,' STSWM: FATAL ERROR IN ROUTINE INFLD:',/,
     $                 ' CANNOT READ FIELD VARIABLE (M=',I4,
     $                 ',N=',I4,')')
                  STOP
               ENDIF
               DIMS(1) = 2
               CALL NCVGT1(CDFID,VARID,DIMS,DIMAG,IRET)
               IF (IRET .NE. 0) THEN
                  WRITE(0,510) M,N
  510             FORMAT(/,' STSWM: FATAL ERROR IN ROUTINE INFLD:',/,
     $                 ' CANNOT READ FIELD VARIABLE (M=',I4,
     $                 ',N=',I4,')')
                  STOP
               ENDIF
               VAR(IDSP(M,N)) = CMPLX(DREAL,DIMAG)
            ELSE
C
C              UNDEFINED COEFFICIENT (ASSUMING 0.0)
C
               VAR(IDSP(M,N)) = CMPLX(0.0,0.0)
            ENDIF
C                                                                               
 1040    CONTINUE                                                               
 1060 CONTINUE                                                                  
C                                                                               
      RETURN                                                                    
      END                                                                       
