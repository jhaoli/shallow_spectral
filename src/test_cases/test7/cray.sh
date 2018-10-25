# QSUB -lm 2.0Mw
# QSUB -lM 2.0Mw
# QSUB -q prem
# QSUB -lt 1100
# QSUB -lT 1200
# QSUB -s /bin/sh -x

#
# experiment number 
#
export EXPERIMENT
EXPERIMENT=0007
#
# define Mass Store Directory for model output	
#
export MS_DIR
MS_DIR=/JAKOB/SHALLOW/EXP/EXP.${EXPERIMENT}
#
# define Input directory for experiment definition
#
export SRC_DIR
SRC_DIR=/crestone/u0/jakob/shallow
#
# define output filename for NCAR graphics
#
export NCARG_GKS_OUTPUT
NCARG_GKS_OUTPUT=gmeta
#
# use temporary directory on Cray
#
cd ${TMPDIR}
#
# copy executable into TMPDIR
#
cp ${SRC_DIR}/src/stswm . 
#
# copy input spectral coefficients into TMPDIR
#
msread VDGDATA.cdf /JAKOB/SHALLOW/EXP/EXP.0077/REFDATA.cdf
#
# copy input parameters into TMPDIR
#
cp ${SRC_DIR}/test7/exp.${EXPERIMENT} exp.data
#
# save job script on Mass Store
#
mswrite -t 4000 ${SRC_DIR}/test7/cray.sh ${MS_DIR}/cray.sh
#
# save experiment input file on Mass Store
#
mswrite -t 4000 exp.data    ${MS_DIR}/exp.data
#
# setup job accounting
#
ja >> stdout 2>> stderr
#
# execute model
#
stswm < exp.data >> stdout 2>> stderr
#
# summarize job costs
#
ja -s -t >> stdout 2>> stderr
#
# save reference solution on Mass Store
#
mswrite -t 4000 REF0007.cdf ${MS_DIR}/REFDATA.cdf
#
# save text and graphics output files on Mass Store
#
mswrite -t 4000 stdout      ${MS_DIR}/stdout
mswrite -t 4000 stderr      ${MS_DIR}/stderr
mswrite -t 4000 gmeta       ${MS_DIR}/gmeta
#
# copy output files to temporary data directory for interactive
# postprocessing on home machine
#
netng FLNM=stdout       flnm=/d1/jakob/stdout.${EXPERIMENT} \
                        host=cgdisis.cgd.ucar.edu
netng FLNM=stderr       flnm=/d1/jakob/stderr.${EXPERIMENT} \
                        host=cgdisis.cgd.ucar.edu
netng FLNM=gmeta  DF=bi flnm=/d1/jakob/gmeta.${EXPERIMENT}  \
                        host=cgdisis.cgd.ucar.edu
#
# end of file
#
