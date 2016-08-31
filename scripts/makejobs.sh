#!/bin/sh
echo "Usage: makejobs job# [file#]"

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

#Replicas are second argument
if [[ ! -z "$2" ]]; then
        NREPLICAS=$2
else
        NREPLICAS="1:3"
fi

#Number of running processors is third argument
if [[ ! -z "$3" ]]; then
        NPROCS=$3
else
	NPROCS=Inf
fi

BASEDIR="run${1}"
mkdir ${BASEDIR}
cd ${BASEDIR}
rm *.job
rm *.o*
rm *.e*
rm *.log

echo "Creating jobs for ${BASEDIR}, replicas=${NREPLICAS}, nprocs=${NPROCS}"

cat<<EOF | matlab -nodisplay
ModelWork_makeJobList('$PROJECT',[],${1},${NREPLICAS},${NPROCS})
EOF

cd ..
cd /home/la67/${PROJECT}/scripts
