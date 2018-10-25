      SUBROUTINE OUTPTP(ZETASC,DIVSC,PHISC,PHIBAR,LDIAG,NEW,STEP) 
C                                                                              
C WRITES THE PROGNOSTIC FIELDS TO A NETCDF FILE
C                                                                              
C CALLED BY: STSWM
C CALLS: PRNT
C
C REVISIONS:
C 7-10-92 CHANGE TO CCM CODING CONVENTIONS (R. JAKOB)
C
C---- Model Parameters -------------------------------------------------
C
      INCLUDE 'params.i'
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
C VORTICITY COEFFICIENTS
      COMPLEX ZETASC(NALP)
C DIVERGENCE COEFFICIENTS
      COMPLEX DIVSC(NALP)
C GEOPOTENTIAL COEFFICIENTS
      COMPLEX PHISC(NALP)
C GLOBAL MEAN STEADY GEOPOTENTIAL
      REAL PHIBAR
C ARRAY INDICES
      INTEGER LDIAG(0:NN,2)
C .TRUE., IF NEW FILE
      LOGICAL NEW
C TIMESTEP INDEX
      INTEGER STEP
C
C------ Local Variables ------------------------------------------------
C
C     NETCDF DECLARATIONS
      INCLUDE 'netcdf.inc'
C
C     NETCDF VARIABLES
C
C     ERROR RETURN CODE
      INTEGER IRET
C     NETCDF ID
      INTEGER CDFID
C     DIMENSION IDS
      INTEGER CXDIM, MMDIM, NNDIM, TDIM
C     VARIABLE IDS
      INTEGER CASEID, ROTID, TIMEID, MMID,NNID, KKID,
     $        ZETAID, DIVID, PHIID
C     VARIABLE SHAPES
      INTEGER DIMS(4)
C
C----- Executable Statements -------------------------------------------
C
C     SET ERROR HANDLING FOR NETCDF
C     (PRINT MESSAGES, BUT DO NOT TERMINATE)
C
      CALL NCPOPT(NCVERBOS)
C
C     OPEN FILE FOR OUTPUT
C
      IF (.NOT. NEW) THEN
C
C        OPEN EXISTING NETCDF FILE
C
         CDFID = NCOPN(FNOUT,NCWRITE,IRET)
         IF (IRET .NE. 0) THEN
            WRITE(0,570) FNOUT
  570       FORMAT(/,' STSWM: FATAL ERROR IN ROUTINE OUTPTP:',/,
     $         ' CANNOT OPEN NETCDF FILE: ',A80,/,
     $         ' CHECK FOR CORRECT FILE (NAMELIST PARAMETER FNOUT) !')
            STOP
         ENDIF
C
C        GET DIMENSION IDS
C
         CXDIM = NCDID(CDFID,'complex',IRET)
         MMDIM = NCDID(CDFID,'m-wave+1',IRET)
         NNDIM = NCDID(CDFID,'n-wave+1',IRET)
         TDIM  = NCDID(CDFID,'timestep',IRET)
C
C        GET VARIABLE IDS
C
         CASEID = NCVID(CDFID,'testcase',IRET)
         ROTID  = NCVID(CDFID,'rot_angle',IRET)
         TIMEID = NCVID(CDFID,'time',IRET)
         MMID   = NCVID(CDFID,'mm',IRET)
         NNID   = NCVID(CDFID,'nn',IRET)
         KKID   = NCVID(CDFID,'kk',IRET)
         ZETAID = NCVID(CDFID,'vorticity',IRET)
         DIVID  = NCVID(CDFID,'divergence',IRET)
         PHIID  = NCVID(CDFID,'geopotential',IRET)
      ELSE
C
C        DEFINE NETCDF FILE
C
         CDFID = NCCRE(FNOUT,NCNOCLOB,IRET)
         IF (IRET .NE. 0) THEN
            WRITE(0,580) FNOUT
  580       FORMAT(/,' STSWM: FATAL ERROR IN ROUTINE OUTPTP:',/,
     $         ' CANNOT CREATE NETCDF FILE: ',A80,/,
     $         ' CHECK FOR CORRECT FILE (NAMELIST PARAMETER FNOUT) !')
            STOP
         ENDIF
C
C        DEFINE DIMENSIONS
C
         CXDIM = NCDDEF(CDFID,'complex',2,IRET)
         MMDIM = NCDDEF(CDFID,'m-wave+1',RTRUNC+1,IRET)
         NNDIM = NCDDEF(CDFID,'n-wave+1',RTRUNC+1,IRET)
         TDIM  = NCDDEF(CDFID,'timestep',NCUNLIM,IRET)
C
C        DEFINE VARIABLES
C
         CASEID = NCVDEF(CDFID,'testcase',NCLONG,0,0,IRET)
         ROTID  = NCVDEF(CDFID,'rot_angle',NCFLOAT,0,0,IRET)
         DIMS(1) = TDIM
         TIMEID = NCVDEF(CDFID,'time',NCFLOAT,1,DIMS,IRET)
         MMID   = NCVDEF(CDFID,'mm',NCLONG,0,0,IRET)
         NNID   = NCVDEF(CDFID,'nn',NCLONG,0,0,IRET)
         KKID   = NCVDEF(CDFID,'kk',NCLONG,0,0,IRET)
         DIMS(4) = TDIM
         DIMS(3) = NNDIM
         DIMS(2) = MMDIM
         DIMS(1) = CXDIM
         ZETAID = NCVDEF(CDFID,'vorticity',NCFLOAT,4,DIMS,IRET)
         DIVID  = NCVDEF(CDFID,'divergence',NCFLOAT,4,DIMS,IRET)
         PHIID  = NCVDEF(CDFID,'geopotential',NCFLOAT,4,DIMS,IRET)
C
C        ASSIGN ATTRIBUTES
C
         CALL NCAPTC(CDFID,ROTID, 'units',NCCHAR,6,'radian',IRET)
         CALL NCAPTC(CDFID,TIMEID,'units',NCCHAR,4,'hour',IRET)
         CALL NCAPTC(CDFID,ZETAID,'units',NCCHAR,8,'second-1',IRET)
         CALL NCAPTC(CDFID,DIVID, 'units',NCCHAR,8,'second-1',IRET)
         CALL NCAPTC(CDFID,PHIID, 'units',NCCHAR,15,
     $               'meter2 second-2',IRET)
C
C        LEAVE DEFINE MODE
C
         CALL NCENDF(CDFID,IRET)
         IF (IRET .NE. 0) THEN
            WRITE(0,590)
  590       FORMAT(/,' STSWM: FATAL ERROR IN ROUTINE OUTPTP:',/,
     $            ' CANNOT WRITE TIMESTEP VARIABLE TO NETCDF FILE')
            STOP
         ENDIF
C
C        SET TEST RUN PARAMETERS
C
         CALL NCVPT1(CDFID,CASEID,DIMS,ICOND,IRET)
         IF (IRET .NE. 0) THEN
            WRITE(0,600)
  600       FORMAT(/,' STSWM: FATAL ERROR IN ROUTINE OUTPTP:',/,
     $            ' CANNOT WRITE TESTCASE VARIABLE TO NETCDF FILE')
            STOP
         ENDIF
         CALL NCVPT1(CDFID,MMID,DIMS,RTRUNC,IRET)
         IF (IRET .NE. 0) THEN
            WRITE(0,610)
  610       FORMAT(/,' STSWM: FATAL ERROR IN ROUTINE OUTPTP:',/,
     $            ' CANNOT WRITE M-TRUNCATION VARIABLE TO NETCDF FILE')
            STOP
         ENDIF
         CALL NCVPT1(CDFID,NNID,DIMS,RTRUNC,IRET)
         IF (IRET .NE. 0) THEN
            WRITE(0,620)
  620       FORMAT(/,' STSWM: FATAL ERROR IN ROUTINE OUTPTP:',/,
     $            ' CANNOT WRITE N-TRUNCATION VARIABLE TO NETCDF FILE')
            STOP
         ENDIF
         CALL NCVPT1(CDFID,KKID,DIMS,RTRUNC,IRET)
         IF (IRET .NE. 0) THEN
            WRITE(0,630)
  630       FORMAT(/,' STSWM: FATAL ERROR IN ROUTINE OUTPTP:',/,
     $            ' CANNOT WRITE K-TRUNCATION VARIABLE TO NETCDF FILE')
            STOP
         ENDIF
         CALL NCVPT1(CDFID,ROTID,DIMS,ALPHA,IRET)
         IF (IRET .NE. 0) THEN
            WRITE(0,640)
  640       FORMAT(/,' STSWM: FATAL ERROR IN ROUTINE OUTPTP:',/,
     $            ' CANNOT WRITE ROTATION ANGLE VARIABLE TO',
     $            ' NETCDF FILE')
            STOP
         ENDIF
      ENDIF
C
C     CHECK TRUNCATION LIMITS FOR OUTPUT OF MODEL STATE
C
      IF ((RTRUNC .GT. MM) .OR. (RTRUNC .GT. NN)) THEN
         WRITE(0,650) RTRUNC,MM,NN
  650    FORMAT (/,' STSWM: WARNING IN SUBROUTINE OUTPTP:',/,
     $   ' SOME SPECTRAL COEFFICIENTS FOR MODEL STATE OUTPUT',/,
     $   ' WILL BE SET TO ZERO, SINCE TRUNCATION PARAMETER RTRUNC',/,
     $   ' IN FILE ''PARAMS.i'' GREATER THAN MM OR NN !',/,
     $   ' RTRUNC = ',I4,' MM = ',I4,' NN = ',I4)
      ELSE
         WRITE(0,660) RTRUNC
  660    FORMAT (/,' STSWM: NOTE IN ROUTINE OUTPTP:',/,
     $       ' SPECTRAL COEFFICIENTS ARE TRUNCATED',/,
     $       ' AT T-',I4,' AS SPECIFIED BY PARAMETER RTRUNC')
      ENDIF
C
C     SET TIMESTEP
C
      DIMS(1) = STEP
      CALL NCVPT1(CDFID,TIMEID,DIMS,TAU,IRET)
      IF (IRET .NE. 0) THEN
         WRITE(0,670)
  670    FORMAT(/,' STSWM: FATAL ERROR IN ROUTINE OUTPTP:',/,
     $         ' CANNOT WRITE TIMESTEP VARIABLE TO NETCDF FILE')
         STOP
      ENDIF
C
C     WRITE VORTICITY
C
      CALL PRNT(CDFID,ZETAID,STEP,ZETASC,LDIAG)                             
C                                                                               
C     WRITE DIVERGENCE
C                                                                               
      CALL PRNT(CDFID,DIVID,STEP,DIVSC,LDIAG)                             
C
C     ADD STEADY STATE GLOBAL AVERAGE FLOW PHIBAR
C     ASSOC. LEGENDRE POLYNOMIAL P(m=0,n=0) = 1/2*SQRT(2)
C
      PHISC(1) = PHISC(1) + PHIBAR*SQRT(2.0)
C                                                                               
C     WRITE GEOPOTENTIAL
C                                                                               
      CALL PRNT(CDFID,PHIID,STEP,PHISC,LDIAG)                             
C
C     CLOSE FILE
C
      CALL NCCLOS(CDFID,IRET)
      IF (IRET .NE. 0) THEN
         WRITE(0,680)
  680    FORMAT(/,' STSWM: FATAL ERROR IN ROUTINE OUTPTP:',/,
     $         ' CANNOT CLOSE NETCDF FILE')
         STOP
      ENDIF
C
      RETURN                                                                    
      END                                                                       

