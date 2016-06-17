#!/bin/sh
#PBS -o localhost:${PBS_O_WORKDIR}/
#PBS -e localhost:${PBS_O_WORKDIR}/
#PBS -M la67@nyu.edu

module purge
#. /etc/profile.d/modules.sh

# Use Intel compiler
module load matlab
export MATLABPATH=${MATLABPATH}:/home/la67/${NAME}/matlab:/home/la67/MATLAB
source /home/la67/MATLAB/setpath.sh

#Check if running as an array job
if [[ ! -z "$PBS_ARRAYID" ]]; then
	IID=${PBS_ARRAYID}
fi
#Check if running as an array job
if [[ ! -z "$SGE_TASK_ID" ]]; then
        IID=${SGE_TASK_ID}
fi

# Run the program
echo ${WORKDIR} ${PROJECT} ${IID}.job

cat<<EOF | matlab -nodisplay
addpath(genpath('/home/la67/MATLAB'));
addpath(genpath('/home/la67/${NAME}'));
cd('${WORKDIR}');
nsamples=${NSAMPLES}; % MCMC samples
continueflag=${CONTINUE}; % Continue flag
ModelWork_batchEval('${PROJECT}',[],'${IID}.job','procid',${IID},'nsamples',nsamples,'continueflag',continueflag);
EOF
