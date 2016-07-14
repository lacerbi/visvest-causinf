#!/bin/sh
echo "Usage: makestats job# [file#]"

PROJECT="VestBMS"
cd /scratch/la67/${PROJECT}

module purge
#. /etc/profile.d/modules.sh

# Use Intel compiler
module load matlab
source /home/la67/MATLAB/setpath.sh
export MATLABPATH=${MATLABPATH}

FILEID=${1}
FILENAME="'joblist-${FILEID}.txt'"
echo "Input #: ${1}   Output file: ${FILENAME}"

#Number of running processors is second argument
if [[ ! -z "$2" ]]; then
        NPROCS=$2
else
	NPROCS=Inf
fi

NREPLICAS="1"

BASEDIR="stats${1}"
mkdir ${BASEDIR}
cd ${BASEDIR}
rm *.job
rm *.o*
rm *.e*
rm *.log

cat<<EOF | matlab -nodisplay
ModelWork_makeJobList('$PROJECT',[],${1},${NREPLICAS},${NPROCS})
EOF

cd ..
cd /home/la67/${PROJECT}/scripts
