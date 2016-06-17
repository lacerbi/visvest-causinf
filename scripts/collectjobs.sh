#!/bin/sh
echo "Usage: collectjobs job# [file#]"

PROJECT="VestBMS"
BASEDIR="run${1}"
cd /scratch/la67/${PROJECT}/${BASEDIR}

module purge
#. /etc/profile.d/modules.sh

# Use Intel compiler
module load matlab
source /home/la67/MATLAB/setpath.sh
export MATLABPATH=${MATLABPATH}

rm *.job
rm *.o*
rm *.e*
rm *.log

cat<<EOF | matlab -nodisplay
[mbag,modelsummary] = ModelWork_collectFits('$PROJECT','./*',[],[]);
modelsummary
save('${PROJECT}_${1}.mat','mbag','modelsummary');
EOF

cd ..
cd /home/la67/${PROJECT}/scripts
