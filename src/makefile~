SHELL=/bin/sh
#
# UNIX make file for the Spectral Transform Shallow Water Model
# Version 2.0    August 1992
# Ruediger Jakob
#
# The command 
#
# make
#
# will compile the Fortran modules indicated by SRCS into the object
# modules indicated by OBJS and then use LDR to produce an executable
# named COMMAND.
#
# Please choose appropriate machine type with following flags:
#
# sun  = SparcStation
# cray = Cray Y-MP
# ibm  = IBM RS 6000
#
MACHINE=   sun
# MACHINE= cray
# MACHINE= ibm
#
# Link with NCAR Graphics library/ Dummy definitions
#
GRAPHICS= -lncarg -lncarg_gks -lncarg_loc
# GRAPHICS= ncarg.o
#
# Choose double precision version of graphics code (64 bit, for Sun)
#
# NCARG_PARAMETER_FILE= /local/isis/sun4/lib/ncarg_dp/ncarg.params
#
# Link with netCDF library/ Dummy definitions
#
NETCDF= -lnetcdf
# NETCDF= netcdf.o
#
# NetCDF include file name
#
CDFINC= /usr/local/include/netcdf.inc
#
# Link with NAG library / Fortran replacement
#
# D01AHE= -lnag
D01AHE= adquad.o
#
# Link with ECMWF FFT library / Fortran replacement
#
# FFT= -lecmfft
FFT= fft99f.o
#
# The parameters SRCS and OBJS should not need to be changed.  If, however,
# you need to add a new module add the name of the source module to the
# SRCS parameter and add the name of the resulting object file to the OBJS
# parameter.  The new modules are not limited to fortran, but may be C, YACC,
# LEX, or CAL.  An explicit rule will need to be added for PASCAL modules.
#

SRCS=   advect.f     agpwrt.f     anlytc.f     bf2.f        bubfnc.f \
	calp.f       ceps.f       comp1.f      cpcnrc.f     cpmpxy.f \
	cpmvxy.f     dbubf.f      dcalp.f      dceps.f      dftfax.f \
	dft991.f     drawcl.f     dtrns.f      duv.f        d2bubf.f \
	epslon.f     erranl.f     eval.f       force.f      ftrnde.f \
	ftrndi.f     ftrnex.f     ftrnim.f     ftrnpe.f     ftrnpi.f \
	ftrnve.f     fu.f         func2.f      glat.f       glats.f  \
	glon.f       infld.f      init.f       inptp.f      input.f  \
	labtop.f     nrgtcs.f     ordleg.f     ortho.f      outptp.f \
	plots.f      prnt.f       rotate.f     shader.f     shtrns.f \
	spcanl.f     step.f       stswm.f      us.f         uv.f     \
	velvct.f     weight.f     zd.f

OBJS=   advect.o     agpwrt.o     anlytc.o     bf2.o        bubfnc.o \
	calp.o       ceps.o       comp1.o      cpcnrc.o     cpmpxy.o \
	cpmvxy.o     dbubf.o      dcalp.o      dceps.o      dftfax.o \
	dft991.o     drawcl.o     dtrns.o      duv.o        d2bubf.o \
	epslon.o     erranl.o     eval.o       force.o      ftrnde.o \
	ftrndi.o     ftrnex.o     ftrnim.o     ftrnpe.o     ftrnpi.o \
	ftrnve.o     fu.o         func2.o      glat.o       glats.o  \
	glon.o       infld.o      init.o       inptp.o      input.o  \
	labtop.o     nrgtcs.o     ordleg.o     ortho.o      outptp.o \
	plots.o      prnt.o       rotate.o     shader.o     shtrns.o \
	spcanl.o     step.o       stswm.o      us.o         uv.o     \
	velvct.o     weight.o     zd.o

#
# Tunable parameters
#
# CF      	Name of the Cray Fortran compiler
# SUNFC         Name of Sun Fortran compiler 
# IBMFC         Name of IBM Fortran compiler
# CRAYFFLAGS	Flags to the Cray Fortran compiler
# SUNFFLAGS     Flags to the Sun Fortran compiler
# IBMFFLAGS     Flags to the IBM Fortran compiler
# CRAYLDR	Program to load the Cray objects into an executable
# SUNLDR        Program to load the Sun objects into an executable
# IBMLDR        Program to load the IBM objects into an executable
# SUNLIBS	List of Sun libraries to use 
# CRAYLIBS      List of Cray libraries to use
# IBMLIBS       List of IBM libraries to use
# COMMAND	Name of the executable
#
CF=		cf77
SUNFC= 		f77 
IBMFC=		xlf
CRAYFFLAGS= 	-Wf"-eisz" -Zv \
                -Wd"-Iglat,glon,weight -Sglat.f,glon.f,weight.f"
SUNFFLAGS= 	-g 
IBMFFLAGS=	
CRAYLDR=	segldr -f indef 
SUNLDR=         f77
IBMLDR=         xlf
SUNLIBS=        -L /lib,/usr/lib,/usr/local/lib \
		$(D01AHE) $(FFT) $(NETCDF) $(GRAPHICS)
CRAYLIBS= 	-L /lib,/usr/lib,/usr/local/lib \
		$(D01AHE) $(FFT) $(NETCDF) $(GRAPHICS) 
IBMLIBS=	-L /lib,/usr/lib,/usr/local/lib \
		$(D01AHE) $(FFT) $(NETCDF) $(GRAPHICS)
COMMAND=	stswm
#
# Lines from here on down should not need to be changed.  They are the
# actual rules which make uses to build stswm.
#
$(COMMAND): 	$(OBJS) adquad.o fft99f.o netcdf.o ncarg.o
	@echo "make: Linking object modules and libraries" ; \
	 if test $(MACHINE) = sun ; \
	 then \
	    $(SUNLDR) $(SUNFFLAGS) -o $(COMMAND) $(OBJS) $(SUNLIBS) ; \
	 elif test $(MACHINE) = cray ; \
	 then \
	    $(CRAYLDR) -o $(COMMAND) $(OBJS) $(CRAYLIBS) ; \
	 elif test $(MACHINE) = ibm ; \
	 then \
	    $(IBMLDR) $(IBMFFLAGS) -o $(COMMAND) $(OBJS) $(IBMLIBS) ; \
	 else \
	    echo "make: undefined MACHINE choice" ; \
	 fi
#
# Fortran Compiler only
#
.f.o:
	@echo "make: Compiling model source file" ; \
	 if test $(MACHINE) = sun ; \
	 then \
	    $(SUNFC) $(SUNFFLAGS) $(TARGET_ARCH) -c $< ; \
	 elif test $(MACHINE) = cray ; \
	 then \
	    $(CF) $(CRAYFFLAGS) -c $< ; \
	 elif test $(MACHINE) = ibm ; \
	 then \
	    $(IBMFC) $(IBMFFLAGS) -c $< ; \
	 else \
	    echo "make: undefined MACHINE choice" ; \
	 fi
#
# Library replacements: NAG integration routine
#
adquad.o: ../lib/adquad.f
	@if test $(D01AHE) = adquad.o ; \
	 then \
	   echo "make: Compiling NAG integration replacement file" ; \
	   if test $(MACHINE) = sun ; \
	   then \
	     $(SUNFC) $(SUNFFLAGS) $(TARGET_ARCH) -c ../lib/adquad.f ; \
	   elif test $(MACHINE) = cray ; \
	   then \
	     $(CF) $(CRAYFFLAGS) -c ../lib/adquad.f ; \
	   elif test $(MACHINE) = ibm ; \
	   then \
	     $(IBMFC) $(IBMFFLAGS) -c ../lib/adquad.f ; \
	   else \
	     echo "make: undefined MACHINE choice" ; \
	   fi \
	 fi
#
# Library replacements: ECMWF FFT routine
#
fft99f.o: ../lib/fft99f.f
	@if test $(FFT) = fft99f.o ; \
	 then \
	   echo "make: Compiling ECMWF FFT replacement file" ; \
	   if test $(MACHINE) = sun ; \
	   then \
	     $(SUNFC) $(SUNFFLAGS) $(TARGET_ARCH) -c ../lib/fft99f.f ; \
	   elif test $(MACHINE) = cray ; \
	   then \
	     $(CF) $(CRAYFFLAGS) -c ../lib/fft99f.f ; \
	   elif test $(MACHINE) = ibm ; \
	   then \
	     $(IBMFC) $(IBMFFLAGS) -c ../lib/fft99f.f ; \
	   else \
	     echo "make: undefined MACHINE choice" ; \
	   fi \
	 fi
#
# NetCDF Library: Dummy routine definitions
#
netcdf.o: ../lib/netcdf.f
	@if test $(NETCDF) = netcdf.o ; \
	 then \
	   echo "make: Compiling NetCDF Dummy Definitions" ; \
	   if test $(MACHINE) = sun ; \
	   then \
	     $(SUNFC) $(SUNFFLAGS) $(TARGET_ARCH) -c ../lib/netcdf.f ; \
	   elif test $(MACHINE) = cray ; \
	   then \
	     $(CF) $(CRAYFFLAGS) -c ../lib/netcdf.f ; \
	   elif test $(MACHINE) = ibm ; \
	   then \
	     $(IBMFC) $(IBMFFLAGS) -c ../lib/netcdf.f ; \
	   else \
	     echo "make: undefined MACHINE choice" ; \
	   fi \
	 fi
#
# NetCDF Library declarations
#
netcdf.inc: $(CDFINC)
	@cp $(CDFINC) .
#
# NCAR Graphics Library: Dummy routine definitions
#
ncarg.o: ../lib/ncarg.f
	@if test '$(GRAPHICS)' = ncarg.o ; \
	 then \
	   echo "make: Compiling NCAR Graphics Dummy Definitions" ; \
	   if test $(MACHINE) = sun ; \
	   then \
	     $(SUNFC) $(SUNFFLAGS) $(TARGET_ARCH) -c ../lib/ncarg.f ; \
	   elif test $(MACHINE) = cray ; \
	   then \
	     $(CF) $(CRAYFFLAGS) -c ../lib/ncarg.f ; \
	   elif test $(MACHINE) = ibm ; \
	   then \
	     $(IBMFC) $(IBMFFLAGS) -c ../lib/ncarg.f ; \
	   else \
	     echo "make: undefined MACHINE choice" ; \
	   fi \
	 fi
#
# Dependencies on include files
#
advect.o: params.i          trnsfm.i wrkspc.i         tdvars.i
agpwrt.o:
anlytc.o: params.i consts.i                   finit.i          complt.i
bf2.o:
bubfnc.o:                                     case4.i
calp.o:   params.i
ceps.o:   params.i
comp1.o:  params.i consts.i trnsfm.i wrkspc.i finit.i tdvars.i
cpcnrc.o:
cpmpxy.o: params.i                                             complt.i
cpmvxy.o: params.i                                             complt.i
dbubf.o:                                      case4.i
dcalp.o:  params.i
dceps.o:
dftfax.o:
dft991.o:
drawcl.o:
dtrns.o:  params.i
duv.o:    params.i
d2bubf.o:                                     case4.i
epslon.o: 
erranl.o: params.i consts.i                   finit.i
eval.o:   params.i
force.o:  params.i consts.i          wrkspc.i finit.i
ftrnde.o: params.i          trnsfm.i
ftrndi.o: params.i          trnsfm.i
ftrnex.o: params.i consts.i trnsfm.i wrkspc.i finit.i tdvars.i
ftrnim.o: params.i consts.i trnsfm.i wrkspc.i finit.i tdvars.i
ftrnpe.o: params.i          trnsfm.i
ftrnpi.o: params.i          trnsfm.i
ftrnve.o: params.i          trnsfm.i
fu.o:              consts.i
func2.o:           consts.i
glat.o:   params.i          trnsfm.i
glats.o:
glon.o:   params.i
infld.o:                    netcdf.inc
init.o:   params.i consts.i                   finit.i
inptp.o:           consts.i netcdf.inc
input.o:  params.i consts.i trnsfm.i                           complt.i
labtop.o:
nrgtcs.o: params.i consts.i
ordleg.o:
ortho.o:  params.i          trnsfm.i wrkspc.i
outptp.o: params.i consts.i netcdf.inc
plots.o:  params.i consts.i                                    complt.i
prnt.o:   params.i          netcdf.inc
rotate.o: 
shader.o:
shtrns.o: params.i          trnsfm.i wrkspc.i
spcanl.o: params.i consts.i trnsfm.i
step.o:   params.i consts.i                           tdvars.i
stswm.o:  params.i consts.i trnsfm.i wrkspc.i finit.i tdvars.i
us.o:              consts.i
uv.o:     params.i          trnsfm.i wrkspc.i
velvct.o:
weight.o: params.i          trnsfm.i
zd.o:     params.i          trnsfm.i wrkspc.i

