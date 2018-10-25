      SUBROUTINE PLOTS(NFLAG,PHI,PANL,U,V,ZETA,DIV,L3D,LN)
C                                                                              
C THIS ROUTINE PRODUCES A VARIETY OF PLOTS FOR STANDARD STATE               
C VARIABLES IN THE ASSOCIATED SHALLOW WATER MODEL CODE.  IT RELIES          
C ON THE NCAR PLOT PACKAGE (I.E., IT'S NOT A GENERIC PLOTTING               
C ROUTINE).  OPTIONS ARE TURNED ON AND OFF WITH LOGICAL AND OTHER           
C REAL VARIABLES SPECIFIED IN COMMON BLOCK /COMPLT/.
C
C CALLED BY: ERRANL, STSWM
C CALLS: CPCNRC, CPMPXY, CPSETI, CPSETR, FRAME, GLAT, 
C        MAPDRW, MAPROJ, MAPRS, MAPSET, MAPSTC, MAPSTI,
C        PLCHMQ, SET, VELVCT
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
C Plotting Data
      INCLUDE 'complt.i'
C Workspace
      INCLUDE 'wrkspc.i'
C
C------------ Arguments ------------------------------------------------
C
C     Input
C
C NFLAG = 0 ==> PLOT COMPUTED SOLUTION
C       = 1 ==> PLOT ERRORS
C       = 2 ==> PLOT ANALYTIC SOLUTION
      INTEGER NFLAG
C HEIGHT FIELD
      REAL PHI(NLON+2,NLAT)
C ANALYTIC HEIGHT FIELD (CASE 1 OVERLAYED PLOTS)
      REAL PANL(NLON+2,NLAT)
C EASTWARD WIND FIELD
      REAL U(NLON+2,NLAT)
C NORTHWARD WIND FIELD
      REAL V(NLON+2,NLAT)
C VORTICITY FIELD
      REAL ZETA(NLON+2,NLAT,L3D)
C DIVERGENCE FIELD
      REAL DIV(NLON+2,NLAT,L3D)
C NUMBER OF TIMELEVELS 
C     (=1 IF CALLED FROM SUBROUTINE ERRANL FOR ANALYTIC SOL. AND ERROR,
C      =LVLS IF CALLED FROM MAIN PROGRAM FOR COMPUTED FIELDS)
      INTEGER L3D
C CURRENT INDEX OF PROGNOSTIC VARIABLES
      INTEGER LN
C
C------ Local Parameters -----------------------------------------------
C
      REAL PI
      PARAMETER (PI=3.141592653589793)
C
C     DASH PATTERNS
C     DASHFL = FULL LINES
C     DASHAS = ALL STIPPLED
C     DASHNS = NEGATIVE CONTOURS STIPPLED, NONNEGATIVE CONTOURS FULL LINES
C
      INTEGER DASHFL,DASHAS,DASHNS
      PARAMETER (DASHFL = 0, DASHAS = 585, DASHNS = -585)
C
C------ Local Variables ------------------------------------------------
C
      REAL UPLT(NLON+1,NLAT),VPLT(NLON+1,NLAT)
      EQUIVALENCE (UPLT(1,1), WS6(1,1)),
     $            (VPLT(1,1), WS7(1,1))
C
C     TEMPORARIES FOR POLE CLUTTER
C
      INTEGER I,J,ISFAC,ISKIP,ICNTR
      REAL RLAT
C
C     CHARACTER ARRAYS FOR PLOT TITLES
C
      CHARACTER*60 GPH1,GPH2,GPH3,GPH4,GPH5,GPH6,GPH7,GPHSPC
      CHARACTER*39 MODEL
C
C     CONTOUR INTERVAL FOR PLOTS
C
      REAL GINC, UINC, VINC, ZINC, DINC, VCTLO, VCTHI
      SAVE GINC, UINC, VINC, ZINC, DINC, VCTLO, VCTHI
C
C     PROJECTION BOUNDARIES
C     (RESTRICTED FOR POLAR STEREOGRAPHIC PROJECTION)
C
      REAL PRB1(2),PRB2(2),PRB3(2),PRB4(2)
      SAVE PRB1,PRB2,PRB3,PRB4
C
C----- External Functions ----------------------------------------------
C
C     LATITUDES OF GRID
C
      EXTERNAL GLAT
      REAL GLAT
C
C----- Initialized Variables -------------------------------------------
C
C     DEFAULT IS AUTOMATIC CHOICE OF CONTOUR INTERVAL FOR PLTS
C     (VALUE = 0 IMPLIES AUTOMATIC CHOICE OF BOUNDS)
C
      DATA GINC / 0.0 /
      DATA UINC / 0.0 /
      DATA VINC / 0.0 /
      DATA ZINC / 0.0 /
      DATA DINC / 0.0 /
      DATA VCTLO,VCTHI  / 2*0.0 /
      DATA PRB1,PRB2,PRB3,PRB4 / 8*0.0 /
C
C----- Executable Statements -------------------------------------------
C
C     CHECK IF PLOTTING REQUESTED
C
      IF (LGPHS) THEN                                                   
C                                                                               
C     MODEL CHARACTERIZATION
C
      WRITE(MODEL,600) ICOND,CHEXP,STRUNC,TAU,INT(DT)
  600 FORMAT('TEST ',I2,',EXP.',A4,',',A6,',T=',F5.1,',DT=',I4)
C
C     CHECK WHICH TYPE OF PLOTS FOR CORRECT TITLE
C
      IF (NFLAG .EQ. 0) THEN
C                                                                               
C        ENCODE TITLES FOR COMPUTED FIELD GRAPHS   
C                                                                               
         GPH1 = 'MODEL HEIGHT         ' // MODEL
         GPH2 = 'MODEL U VELOCITY     ' // MODEL
         GPH3 = 'MODEL V VELOCITY     ' // MODEL
         GPH4 = 'MODEL WIND VECTOR    ' // MODEL
         GPH5 = 'MODEL HEIGHT/WIND    ' // MODEL
         GPH6 = 'MODEL VORTICITY      ' // MODEL
         GPH7 = 'MODEL DIVERGENCE     ' // MODEL
C
      ELSE IF (NFLAG .EQ. 1) THEN
C                                                                               
C        ENCODE TITLES DENOTING ERRORS 
C                                                                               
         GPH1 = 'HEIGHT ERROR         ' // MODEL
         GPH2 = 'U VELOCITY ERROR     ' // MODEL
         GPH3 = 'V VELOCITY ERROR     ' // MODEL
         GPH4 = 'WIND VECTOR ERROR    ' // MODEL
         GPH5 = 'HEIGHT/WIND ERROR    ' // MODEL
         GPH6 = 'VORTICITY ERROR      ' // MODEL
         GPH7 = 'DIVERGENCE ERROR     ' // MODEL
C
      ELSE IF (NFLAG .EQ. 2) THEN
C                                                                               
C        ENCODE TITLES FOR ANALYTIC SOLUTION
C                                                                               
         GPH1 = 'ANALYTIC HEIGHT      ' // MODEL
         GPH2 = 'ANALYTIC U VELOCITY  ' // MODEL
         GPH3 = 'ANALYTIC V VELOCITY  ' // MODEL
         GPH4 = 'ANALYTIC WIND VECTOR ' // MODEL
         GPH5 = 'ANALYTIC HEIGHT/WIND ' // MODEL
         GPH6 = 'ANALYTIC VORTICITY   ' // MODEL
         GPH7 = 'ANALYTIC DIVERGENCE  ' // MODEL
         GPHSPC='MODEL & ANAL. HEIGHT ' // MODEL
C
      ELSE
C
         WRITE(0,830) NFLAG
  830    FORMAT(/,' STSWM: FATAL ERROR IN ROUTINE PLOTS:',/,
     $    ' ILLEGAL NFLAG ARGUMENT IN CALL',/,
     $    ' NFLAG = ',I2)
         STOP
C
      ENDIF
C
C     SET CONTOUR INTERVAL FOR CASE 1
C
      IF ((ICOND .EQ. 1) .AND. (NFLAG .NE. 1)) THEN
         GINC = 100.0
      ELSE
         GINC = 0.0
      ENDIF
C
C     DEFINE GAUSSIAN GRID MAPPING (SEE SUBROUTINES CPMPXY/CPMVXY)
C
      CALL CPSETI('MAP - MAPPING FLAG',3)
C
C     NUMBER OF GRIDPOINTS IN ARRAYS
C
      CALL CPSETR('XC1 - GAUSSIAN REMAPPING',1.)
      CALL CPSETR('XCM - GAUSSIAN REMAPPING',REAL(NLON+1))
      CALL CPSETR('YC1 - GAUSSIAN REMAPPING',1.)
      CALL CPSETR('YCN - GAUSSIAN REMAPPING',REAL(NLAT))
C
C     INVISIBLE(PROJECTION) AND DELETED(CROWDING) POINTS 
C
      CALL CPSETR('SPV - SPECIAL VALUE',SPVAL)
C
C     OUT OF RANGE VALUES
C
      CALL CPSETR('ORV - OUT OF RANGE VALUE',1.0E12)
C
C     SET FLAG FOR CONTINENTAL OUTLINES
C
      IF (LCONT) THEN                                                           
         CALL MAPSTC('OU','CO')
         CALL MAPSTI('DO - DOTTING OF OUTLINES',1)
      ELSE                                                                      
         CALL MAPSTC('OU','NO')
      ENDIF                                                                     
C
C     SET UP MAP PROJECTION
C
      IF (LPSP) THEN              
C
C        NORTH POLAR STEREOGRAPHIC PROJECTION
C
         CALL MAPROJ('ST',90.0-ALPHA*180.0/PI,-180.0,0.0)
         LNORTH = .TRUE.
C
C        COPY ROTATION ANGLE INTO POROT FOR USE IN ROUTINE CPMPXY
C
         POROT = ALPHA
C
C        HEMISPHERE BOUNDARIES
C
         CALL MAPSTI('EL',1)
         PRB1(1) = 90.0
         PRB2(1) = 90.0
         PRB3(1) = 90.0
         PRB4(1) = 90.0
         CALL MAPSET('AN',PRB1,PRB2,PRB3,PRB4)
C
      ELSEIF (LOP) THEN
C
C        ORTHOGRAPHIC PROJECTION
C
         CALL MAPROJ('OR',POLAT,POLNG,POROT)
         CALL MAPSET('MA',PRB1,PRB2,PRB3,PRB4)
         CALL MAPSTI('EL',1)
C
      ELSEIF (LCP) THEN
C
C        CYLINDRICAL PROJECTION
C
         CALL MAPROJ('CE',POLAT,POLNG,POROT)
         CALL MAPSET('MA',PRB1,PRB2,PRB3,PRB4)
      ENDIF
C
C     SHOW LAT/LONG GRID
C
      CALL MAPSTI('GR',30)
C
C     INITIALIZE MAPPING
C
      CALL MAPINT
C
C     DISABLE MAPPING-CHANGE BY CPCNRC
C
      CALL CPSETI('SET',0)
C
C     USE ARRAYS UPLT, VPLT FOR PLOTTING (SPCL VALUE MODS)          
C     PLOTTING ARRAYS HAVE AN EXTRA ZONAL POINT TO CLOSE LINES ON PLOTS,        
C     I.E., PHI(1,J)=PHI(NLON+1,J) AND RLOND(1)=RLOND(NLON+1)                     
C                                                                               
C     LOAD PLOTTING ARRAYS WITH U, V
C                                                                               
      IF (LG .OR. LU .OR. LV .OR. LVV .OR. LVVG 
     $       .OR. LZ .OR. LD) THEN                                              
      DO 215 J=1,NLAT                                                        
C
C        COPY WIND DATA FOR VECTOR PLOTS
C
         IF (LVV .OR. LVVG) THEN
C
C        ELIMINATE CLUTTER OF POINTS NEAR THE POLE WHEN PLOTTING
C        VELOCITY VECTORS IN STEREOPGRAPHIC OR ORTHOGRAPHIC PROJECTION
C
         IF (LPSP .OR. LOP) THEN
            RLAT = GLAT(J)
C
C           LATITUDE = RLAT = GLAT(J)
C
            ISFAC = 1.0/COS(RLAT)
            ISKIP = 1
            ICNTR = 0
C
  795       CONTINUE
            ISFAC = ISFAC/2
            IF (ISFAC .EQ. 0) GO TO 796
            ISKIP = ISKIP*2
            ICNTR = ICNTR+1
            IF (ICNTR .GE. 16) THEN
               WRITE(0,840) ICNTR
  840          FORMAT(/,' STSWM: FATAL ERROR IN ROUTINE PLOTS:',/,
     $              ' TOO MANY CLUTTERING POINTS NEAR POLE',/,
     $              ' ICNTR = ',I2)
               STOP
            ENDIF
            GO TO 795
C
  796       CONTINUE
         ELSE
C
C           NO CLUTTER IN CYLINDRICAL PROJECTION
C
            ISKIP = 1
         ENDIF
C
C        KEEP ONLY EACH (ISKIP)TH POINT
C
         DO 799 I=1,NLON
            IF (MOD(I,ISKIP) .NE. 0) THEN
               UPLT(I,J) = SPVAL
               VPLT(I,J) = SPVAL
            ELSE
               UPLT(I,J) = U(I,J)
               VPLT(I,J) = V(I,J)
            ENDIF
  799    CONTINUE
         ENDIF
C
C           EXTRA POINT FOR CONTINOUS CONTOUR PLOTS
C
            U(NLON+1,J) = U(1,J)
            V(NLON+1,J) = V(1,J)
            PHI(NLON+1,J)  = PHI(1,J)
            PANL(NLON+1,J) = PANL(1,J)
            ZETA(NLON+1,J,LN) = ZETA(1,J,LN)
            DIV(NLON+1,J,LN)  = DIV(1,J,LN)
  215    CONTINUE                                                               
      ENDIF                                                                     
C                                                                               
C     PLOT NUMBER 1 --- CONTOUR PLOT OF HEIGHT        
C                                                                               
      IF (LG) THEN                                                              
         CALL MAPDRW
         IF ((ICOND .EQ. 1) .AND. (NFLAG .EQ. 2)) THEN
C
C           SHOW MODEL HEIGHT (ALL DASHED)
C           DEFINE CONTOUR INTERVAL 100 m
C
            CALL CPCNRC(PHI,NLON+2,NLON+1,NLAT,GINC,DASHAS)
C
C           OVERLAP WITH ANALYTIC SOLUTION (FULL LINES)
C           DEFINE CONTOUR INTERVAL 100 m
C
            CALL CPSETC('ILT - INFORMATIONAL LABEL',' ')
            CALL CPCNRC(PANL,NLON+2,NLON+1,NLAT,GINC,DASHFL)
            CALL CPSETC('ILT - INFORMATIONAL LABEL',
     $          'CONTOUR FROM $CMN$ TO $CMX$ BY $CIU$')
            CALL LABTOP(GPHSPC,0.015)
         ELSEIF ((ICOND .EQ. 2) .AND. (NFLAG .EQ. 1)) THEN
               CALL CPCNRC(PHI,NLON+2,NLON+1,NLAT,GINC,DASHFL)
               CALL LABTOP(GPH1,0.015)
         ELSEIF (ICOND .EQ. 3) THEN
               CALL CPCNRC(PHI,NLON+2,NLON+1,NLAT,GINC,DASHFL)
               CALL LABTOP(GPH1,0.015)
         ELSE
            CALL CPCNRC(PHI,NLON+2,NLON+1,NLAT,GINC,DASHNS)
            CALL LABTOP(GPH1,0.015)
         ENDIF
         CALL FRAME                                                             
         IF (LPSP) THEN
            LNORTH = .FALSE.
            CALL MAPROJ('ST',-90.0+ALPHA*180.0/PI,0.0,0.0)
            CALL MAPDRW
            IF ((ICOND .EQ. 1) .AND. (NFLAG .EQ. 2)) THEN
C
C              SHOW MODEL HEIGHT (ALL DASHED)
C              DEFINE CONTOUR INTERVAL 100 m
C
               CALL CPCNRC(PHI,NLON+2,NLON+1,NLAT,GINC,DASHAS)
C
C              OVERLAP WITH ANALYTIC SOLUTION (FULL LINES)
C              DEFINE CONTOUR INTERVAL 100 m
C
               CALL CPSETC('ILT - INFORMATIONAL LABEL',' ')
               CALL CPCNRC(PANL,NLON+2,NLON+1,NLAT,GINC,DASHFL)
               CALL CPSETC('ILT - INFORMATIONAL LABEL',
     $             'CONTOUR FROM $CMN$ TO $CMX$ BY $CIU$')
               CALL LABTOP(GPHSPC,0.015)
            ELSEIF ((ICOND .EQ. 2) .AND. (NFLAG .EQ. 1)) THEN
               CALL CPCNRC(PHI,NLON+2,NLON+1,NLAT,GINC,DASHFL)
               CALL LABTOP(GPH1,0.015)
            ELSE
               CALL CPCNRC(PHI,NLON+2,NLON+1,NLAT,GINC,DASHNS)
               CALL LABTOP(GPH1,0.015)
            ENDIF
            CALL FRAME
            CALL MAPROJ('ST',90.0-ALPHA*180.0/PI,-180.0,0.0)
            LNORTH = .TRUE.
         ENDIF
      ENDIF                                                                     
C                                                                               
C     PLOT NUMBER 2 --- CONTOUR PLOT OF U VELOCITY          
C                                                                               
      IF (LU) THEN                                                              
         CALL MAPDRW
         CALL CPCNRC(U,NLON+2,NLON+1,NLAT,UINC,DASHNS)          
         CALL LABTOP(GPH2,0.015)
         CALL FRAME                                                             
         IF (LPSP) THEN
            LNORTH = .FALSE.
            CALL MAPROJ('ST',-90.0+ALPHA*180.0/PI,0.0,0.0)
            CALL MAPDRW
            CALL CPCNRC(U,NLON+2,NLON+1,NLAT,UINC,DASHNS)
            CALL LABTOP(GPH2,0.015)
            CALL FRAME
            CALL MAPROJ('ST',90.0-ALPHA*180.0/PI,-180.0,0.0)
            LNORTH = .TRUE.
         ENDIF
      ENDIF                                                                     
C                                                                               
C     PLOT NUMBER 3 --- CONTOUR PLOT OF V VELOCITY          
C                                                                               
      IF (LV) THEN                                                              
         CALL MAPDRW
         CALL CPCNRC(V,NLON+2,NLON+1,NLAT,VINC,DASHNS)          
         CALL LABTOP(GPH3,0.015)
         CALL FRAME                                                             
         IF (LPSP) THEN
            LNORTH = .FALSE.
            CALL MAPROJ('ST',-90.0+ALPHA*180.0/PI,0.0,0.0)
            CALL MAPDRW
            CALL CPCNRC(V,NLON+2,NLON+1,NLAT,VINC,DASHNS)
            CALL LABTOP(GPH3,0.015)
            CALL FRAME
            CALL MAPROJ('ST',90.0-ALPHA*180.0/PI,-180.0,0.0)
            LNORTH = .TRUE.
         ENDIF
      ENDIF                                                                     
C
C     PLOT NUMBER 4 --- VELOCITY VECTORS                    
C     BE CAREFUL WITH THIS PLOT AT HIGH RESOLUTION                              
C     FURTHER MODIFICATION OF THE DATA MAY BE REQUIRED BEYOND T-63              
C                                                                               
      IF (LVV) THEN                                                             
         CALL MAPDRW
         CALL VELVCT(UPLT,NLON+1,VPLT,NLON+1,NLON,NLAT,                     
     $                VCTLO,VCTHI,1,16,3,SPV)  
         CALL LABTOP(GPH4,0.015)
         CALL FRAME                                                             
         IF (LPSP) THEN
            LNORTH = .FALSE.
            CALL MAPROJ('ST',-90.0+ALPHA*180.0/PI,0.0,0.0)
            CALL MAPDRW
            CALL VELVCT(UPLT,NLON+1,VPLT,NLON+1,NLON,NLAT,
     $                VCTLO,VCTHI,1,16,3,SPV)
            CALL LABTOP(GPH4,0.015)
            CALL FRAME
            CALL MAPROJ('ST',90.0-ALPHA*180.0/PI,-180.0,0.0)
            LNORTH = .TRUE.
         ENDIF
      ENDIF                                                                     
C                                                                               
C     PLOT NUMBER 5 --- OVERLAY OF PLOTS 1 AND 4            
C                                                                               
      IF (LVVG) THEN                                                            
         CALL MAPDRW
         CALL CPSETC('ILT - INFORMATIONAL LABEL',' ')
         CALL CPCNRC(PHI,NLON+2,NLON+1,NLAT,GINC,DASHNS)
         CALL CPSETC('ILT - INFORMATIONAL LABEL',
     $    'CONTOUR FROM $CMN$ TO $CMX$ BY $CIU$')
         CALL VELVCT(UPLT,NLON+1,VPLT,NLON+1,NLON,NLAT,                     
     $                VCTLO,VCTHI,1,16,3,SPV)
         CALL LABTOP(GPH5,0.015)
         CALL FRAME                                                             
         IF (LPSP) THEN
            LNORTH = .FALSE.
            CALL MAPROJ('ST',-90.0+ALPHA*180.0/PI,0.0,0.0)
            CALL MAPDRW
            CALL CPSETC('ILT - INFORMATIONAL LABEL',' ')
            CALL CPCNRC(PHI,NLON+2,NLON+1,NLAT,GINC,DASHNS)
            CALL CPSETC('ILT - INFORMATIONAL LABEL',
     $    'CONTOUR FROM $CMN$ TO $CMX$ BY $CIU$')
            CALL VELVCT(UPLT,NLON+1,VPLT,NLON+1,NLON,NLAT,
     $                VCTLO,VCTHI,1,16,3,SPV)
            CALL LABTOP(GPH5,0.015)
            CALL FRAME
            CALL MAPROJ('ST',90.0-ALPHA*180.0/PI,-180.0,0.0)
            LNORTH = .TRUE.
         ENDIF
      ENDIF                                                                     
C
C     PLOT NUMBER 6 --- CONTOUR PLOT OF VORTICITY                   
C                                                                               
      IF (LZ) THEN                                                              
         CALL MAPDRW
         CALL CPCNRC(ZETA(1,1,LN),NLON+2,NLON+1,NLAT,ZINC,DASHNS)            
         CALL LABTOP(GPH6,0.015)       
         CALL FRAME                                                             
         IF (LPSP) THEN
            LNORTH = .FALSE.
            CALL MAPROJ('ST',-90.0+ALPHA*180.0/PI,0.0,0.0)
            CALL MAPDRW
            CALL CPCNRC(ZETA(1,1,LN),NLON+2,NLON+1,NLAT,ZINC,DASHNS)
            CALL LABTOP(GPH6,0.015)
            CALL FRAME
            CALL MAPROJ('ST',90.0-ALPHA*180.0/PI,-180.0,0.0)
            LNORTH = .TRUE.
         ENDIF
      ENDIF                                                                     
C                                                                               
C     PLOT NUMBER 7 --- CONTOUR PLOT OF DIVERGENCE                  
C                                                                               
      IF (LD) THEN                                                              
         CALL MAPDRW
         CALL CPCNRC(DIV(1,1,LN),NLON+2,NLON+1,NLAT,DINC,DASHNS)            
         CALL LABTOP(GPH7,0.015)       
         CALL FRAME                                                             
         IF (LPSP) THEN
            LNORTH = .FALSE.
            CALL MAPROJ('ST',-90.0+ALPHA*180.0/PI,0.0,0.0)
            CALL MAPDRW
            CALL CPCNRC(DIV(1,1,LN),NLON+2,NLON+1,NLAT,DINC,DASHNS)
            CALL LABTOP(GPH7,0.015)
            CALL FRAME
            CALL MAPROJ('ST',90.0-ALPHA*180.0/PI,-180.0,0.0)
            LNORTH = .TRUE.
         ENDIF
      ENDIF                                                                     
C
C     !!!!!!!!!!!!!!!!!!!!!!!! END GRAPHICS !!!!!!!!!!!!!!!!!!!!!!!!!!!!        
C                                                                               
      ENDIF
C
      RETURN                                                                    
      END                                                                       
