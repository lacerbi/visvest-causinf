#!/bin/sh
echo "Usage: collectjobs job1 [job2...]"
PROJECT="VestBMS"

module purge
#. /etc/profile.d/modules.sh

# Use Intel compiler
module load matlab
source ${HOME}/MATLAB/setpath.sh
export MATLABPATH=${MATLABPATH}

cat<<EOF | matlab -nodisplay
RUNS=str2num('$@')
for iRun=1:numel(RUNS)
	cd(['${SCRATCH}/${PROJECT}/run' num2str(RUNS(iRun))]);
	[mbag,modelsummary] = ModelWork_collectFits('$PROJECT','./*',[],[]);
	modelsummary
	filename=['${SCRATCH}/${PROJECT}/${PROJECT}_' num2str(RUNS(iRun)) '.mat'];
	save(filename,'mbag','modelsummary');
end
EOF

cd ${HOME}/${PROJECT}/scripts
