File: README


SPECTRAL TRANSFORM SHALLOW WATER MODEL (Version 2.0)
Copyright (C) 1992
University Corporation for Atmospheric Research
All Rights Reserved

       Ruediger Jakob
National Center for Atmospheric Research
    Boulder, CO 80307-3000

         August 1992


Contents
--------

1. Software Distribution Conditions
2. Description of Software
3. Directory of Files
4. Corrections and Changes


1. Software Distribution Conditions
-----------------------------------

This spectral transform shallow water model (STSWM) is made available
``as is'' without warranty, expressed or implied, as to its suitability
for general purpose applications. The source code has been made available 
for scientific research and educational purposes only, and should not
be construed to be in the public domain. UCAR does not indemnify any 
infringement of copyright, patent, or trademark through use or 
modification of this software. NCAR/CGD Division support for the software 
is limited to copies of the source code and the distribution of 
a technical note (i.e. no consulting support will be provided).

NCAR should be acknowledged as the source of the software used in any
resulting research or publications. If users make substantial
modifications to the model source code, any references in publications
of simulation results should clearly indicate that it is a derivative
of the NCAR STSWM. Users are requested to send reprints of their work 
with this model when available to the NCAR CGD Division Office.


2. Description of Software
--------------------------

This software is based on the NCAR Technical Note TN-343+STR
``Description of a Global Shallow Water Model Based on the Spectral
Transform Method'' by J.J. Hack and R. Jakob.
A detailed "Description of Software for the Spectral
Transform Shallow Water Model" can be found in the text file
'description.txt' in directory 'docu'. Please read this document 
carefully, as it explains how to compile and run the code, 
and what input and output files are used by the model. 

In addition to the source code, the following libraries and 
subroutines are used: NCAR Graphics library, Fast Fourier Transform 
and adaptive quadrature routines, and the NetCDF library for portable 
storage of reference solutions. The following paragraphs explain how 
to obtain these codes:

Fast Fourier Transform 

The model uses the Fast Fourier Transform routines SET99 and
FFT991 from the European Center for Medium Range Weather
Forecasting (ECMWF). A FORTRAN version of these routines can be
obtained via anonymous FTP from the machine 'ftp.ucar.edu',
file 'dsl/lib/ecmfft/fft99f.f'. A modified version of this
code has been copied into directory 'lib'.

Integration Routine

The model uses an adaptive quadrature routine
D01AHE from the Numerical Algorithms Group (NAG) library, Mark 14.
A FORTRAN version with a similar functionality can be
obtained via anonymous FTP from the machine 'ftp.ucar.edu',
file 'dsl/lib/ncarm/adquad.f'. A modified version which has the 
same name and interface as D01AHE has been copied into 
directory 'lib'.

NCAR Graphics Library

For obtaining a copy of the NCAR graphics package, please contact
the Scientific Computing Division at NCAR, phone: (303) 497-1201
or send email to scdinfo@ncar.ucar.edu. This is the only library
that is not available free of charge.

NetCDF Library

The NetCDF software is distributed via anonymous FTP by the Unidata
Program Center. For UNIX systems, a compressed tar file can be
accessed (in binary mode) from the file pub/netcdf.tar.Z in the
anonymous FTP directory of the machine unidata.ucar.edu.
VMS sites can get a backup saveset of the same software from the
anonymous FTP directory of machine laurel.ucar.edu.


3. Directory of Files
---------------------

This directory contains the following:

README    - the document you are reading
docu/     - subdirectory with descriptive documents
lib/      - subdirectory with source code of FFT and integration routines
netcdf/   - subdirectory with reference solutions in netCDF format
src/      - subdirectory with shallow water model source code 
test1/    - subdirectory for test case 1: Advection of Low
test2/    - subdirectory for test case 2: Solid body rotation
test3/    - subdirectory for test case 3: Mid-latitude Jet Stream
test4/    - subdirectory for test case 4: Forced Low in Jet Stream
test5/    - subdirectory for test case 5: Flow around mountain
test6/    - subdirectory for test case 6: Rossby-Haurwitz Wave
test7/    - subdirectory for test case 7: 500 mb height real data

Except for NetCDF files with filename extension .cdf, all files are 
plain text files. NetCDF files are in binary format: remember to
put FTP into binary mode for accessing these files.

4. Corrections and Changes
--------------------------

Software bugs, along with suggested fixes, should be reported 
to the email address

stswm@ncar.ucar.edu

Messages to this address will be forwarded to the authors. After
verification, any corrections or changes will be appended here.

Feb 25, 1992 (Ruediger Jakob)
Changes for Version 1.1:

- semi-implicit time stepping is now done with spectral coefficients
- vorticity/divergence forcing is no longer spectrally transformed
- corrections to the conservation analysis routines
- added experiment number on all plots
- polar stereographic plots of north and south pole
- cleaned up FFT and integration library replacement routines
- added shading of negative contours and better font type

Mar 10, 1992 (Ruediger Jakob)
Changes for Version 1.2:

- new code for out-of range error detection in graphics
- changed shading for negative parts of plots
- modified makefile to handle replacement files for FFT
  and NAG integration routine auomatically
- eliminated copying in plotting subroutine for non-vector plots
- improved font type for time series plots
- updated code description using R. Sato's suggestions

Apr 10, 1992 (Ruediger Jakob)
Changes for Version 1.3:

- changed name of integration subroutine D01AHF to D01AHE for
  new NAG library Mar 14 naming convention. Also changed name
  of replacement routine adquad.

May 22, 1992 (Ruediger Jakob)
Changes for Version 1.4:

- PHIBAR input variable is now computed automatically
- wind speed and height have changed for test case 5 (flow over
  mountain); the geopotential diffusion is now computed on the 
  height of the free surface, to prevent spurios diffusion caused
  by surface topography
- Expanded makefile to also handle IBM RS 6000
- Use new files netcdf.f and ncarg.f in /lib to solve problem
  of unresolved externals when linking
- Created subdirectory /netcdf for high resolution reference
  solutions (test cases 5,6 & 7)

Aug 18, 1992 (Ruediger Jakob)
Changes for Version 2.0:

- each subroutine is now in a separate file; filename is lower-case
  of subroutine name to conform with UNIX conventions
- a wrong parameter value for NLAT at T-42 in version 1.4 was corrected
- the makefile for Sun was changed back to single precision
- the Newton-iteration convergence test in GLATS was relaxed by a 
   factor of 2.0 since it did not converge on the Suns
- the argument list for the routines were reordered to have input
  arguments before output arguments
- the declarations for test case 4 were taken out of FINIT.i into
  a separate include file case4.i
- the routine HIRES was removed; the two calls to INPTP and EVAL
  were moved to routine ANLYTC
- the routine GSSLAT was removed; the copying for the symmetric other
  hemisphere is now done in GLATS
- a routine ORTHO was added to check the orthogonality and ortho-
  normality of the associated Legendre functions. It is enabled with
  the flag LORTHO in STSWM.
- routine ANLYTC now computes the analytic solution for the entire 
  grid, not only one latitude.
- the number of transforms for the semi-implicit case and explicit
  case was reduced using H. Ritchie's algorithm (see MWR 1988 paper)
- To untangle the complex COMP1 routine, it now calls three subroutines
  for the advection equation, and semi-implicit/explicit timestepping
  of the full equations: ADVECT, FTRNEX, FTRNIM
- A fatal error in SPCANL caused by kinetic/potential energy = 0.0
  for any N in the logarithmic plots was fixed
- High/Low labels in the contour plots can be enabled/disabled with the
  flag PLABEL in CPCNRC
- NetCDF files can now be read at any MULTIPLE of the time interval used 
  when writing the file
- the flag LDIFF in STSWM can be used for postprocessing of NetCDF
  spectral coefficient files (this facility has been requested widely
  but is still in the test phase). Use the NAMELIST parameters FNIN 
  and DT to specify file name and time interval of analysis.

Sep 22, 1992 (Ruediger Jakob)
Correction to Version 2.0:

- Added dummy routine definitions in /lib/netcdf.f and /lib/ncarg.f
  (NCPOPT, NCDINQ, and VELDAT) to prevent 'unresolved external names'
  when no NetCDF or NCAR graphics library is present.

Oct 22, 1992 (Ruediger Jakob)
Correction to Version 2.0:

- Updated reference to JCP paper in description.txt

-------- end of document ----------------------------------------------
