#!/bin/sh -x
#
# define unique experiment number
#
export EXPERIMENT;
EXPERIMENT=0001;
#
# remove old files
#
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
