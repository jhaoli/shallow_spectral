#!/bin/sh -x
#
# define unique experiment number
#
export EXPERIMENT;
EXPERIMENT=0003;
#
# define output filename for NCAR graphics
#
export NCARG_GKS_OUTPUT;
NCARG_GKS_OUTPUT=gmeta.${EXPERIMENT};
#
# remove old files
#
rm gmeta.${EXPERIMENT};
rm stdout.${EXPERIMENT}; 
rm stderr.${EXPERIMENT};
#
# execute model
#
../src/stswm < exp.${EXPERIMENT} \
             >> stdout.${EXPERIMENT} \
             2>> stderr.${EXPERIMENT};
#
# end of file
#
