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

#Running processors is second argument
REPLICAS="1:3"

BASEDIR="run${1}"
cd ${BASEDIR}
rm *.job
#rm *.o*
#rm *.e*

cat<<EOF | matlab -nodisplay
temp = load('../${PROJECT}_${1}.mat','mbag');
mbag = temp.mbag;
exitflag = ModelWork_samplingCheck(mbag,'write',${1},${REPLICAS})
EOF

cd ..
cd /home/la67/${PROJECT}/scripts
