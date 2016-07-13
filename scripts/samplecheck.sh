#!/bin/sh
echo "Usage: samplecheck job#"

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
NREPLICAS="3"

BASEDIR="run${1}"
cd ${BASEDIR}
rm *.job
rm *.o*
rm *.e*
rm *.log

cat<<EOF | matlab -nodisplay
temp = load('../${PROJECT}_${RUN}.mat','mbag');
mbag = temp.mbag
exitflag = ModelWork_samplingCheck(mbag,'write',${1},${NREPLICAS})
EOF

cd ..
cd /home/la67/${PROJECT}/scripts
