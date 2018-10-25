      SUBROUTINE PRNT(CDFID,VARID,STEP,VAR,LDIAG)                       
C                                                                              
C WRITES A SPECTRAL COEFFICIENT FIELD
C IN NETCDF FORMAT
C
C CALLED BY: OUTPTP
C CALLS:
C
C REVISIONS:
C 7-13-92 CHANGE TO CCM CODING CONVENTIONS (R. JAKOB)
C
C---- Model Parameters -------------------------------------------------
C
      INCLUDE 'params.i'
C NETCDF DECLARATIONS
      INCLUDE 'netcdf.inc'
C
C------------ Arguments ------------------------------------------------
C
C     Input
C
C NETCFD FILE ID
      INTEGER CDFID
C NETCDF VARIABLE ID
      INTEGER VARID
C TIMESTEP
      INTEGER STEP
C SPECTRAL ARRAY TO BE WRITTEN TO FILE
      COMPLEX VAR(NALP)
C STRUCTURE INFORMATION FOR ARRAY VAR
      INTEGER LDIAG(0:NN,2)                                      
C
C------ Local Variables ------------------------------------------------
C
      INTEGER M,N,IRET
      INTEGER DIMS(4)
      COMPLEX DATA
      REAL PAIR(2)
      EQUIVALENCE (DATA,PAIR(1))
C
C----- Statement Function ----------------------------------------------
C
C     ADDRESS COMPUTATION OF ALP, DALP, EPSIL            
C                                                                               
      INTEGER MDUM,NDUM,IDSP
      IDSP(MDUM,NDUM) = 1 + LDIAG(NDUM-MDUM,2)+MDUM                             
C                                                                               
C----- Executable Statements -------------------------------------------
C
C     SET TIMESTEP
C
      DIMS(4) = STEP
C
C     PUT OUT SPECTRAL COEFFICIENTS FOR TRIANGULAR TRUNCATION RTRUNC
C
      DO 1060 N=0,RTRUNC
         DIMS(3) = N+1
         DO 1040 M=0,N
            DIMS(2) = M+1
C                                                                               
            IF ((M .GT. MM) .OR. (N .GT. MIN0(KK,NN+M))) THEN
C
C              OUT OF SPECTRAL TRUNCATION RANGE
C
               DATA = CMPLX(0.0,0.0)
            ELSE
               DATA = VAR(IDSP(M,N))
            ENDIF
C                                                                               
C           WRITE REAL AND IMAGINARY PARTS
C
            DIMS(1) = 1
            CALL NCVPT1(CDFID,VARID,DIMS,PAIR(1),IRET)
            IF (IRET .NE. 0) THEN
               WRITE(0,600)
  600          FORMAT(/,' STSWM: FATAL ERROR IN ROUTINE PRNT:',/,
     $               ' CANNOT WRITE FIELD VARIABLE')
               STOP
            ENDIF
            DIMS(1) = 2 
            CALL NCVPT1(CDFID,VARID,DIMS,PAIR(2),IRET)
            IF (IRET .NE. 0) THEN
               WRITE(0,610)
  610          FORMAT(/,' STSWM: FATAL ERROR IN ROUTINE PRNT:',/,
     $               ' CANNOT WRITE FIELD VARIABLE')
               STOP
            ENDIF
C                                                                               
 1040    CONTINUE                                                               
 1060 CONTINUE                                                                  
C                                                                               
      RETURN                                                                    
      END                                                                       
