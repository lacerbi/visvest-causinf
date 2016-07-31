#!/bin/sh
echo "Usage: collectstats job1 [job2...]"
PROJECT="VestBMS"

module purge
#. /etc/profile.d/modules.sh

# Use Intel compiler
module load matlab
source /home/la67/MATLAB/setpath.sh
export MATLABPATH=${MATLABPATH}
mkdir /scratch/la67/${PROJECT}/finalstats

cat<<EOF | matlab -nodisplay
RUNS=str2num('$@')
for iRun=1:numel(RUNS)
	cd(['/scratch/la67/${PROJECT}/stats' num2str(RUNS(iRun))]);
	[mbag,modelsummary] = ModelWork_collectFits('$PROJECT','./*',[],[]);
	modelsummary
	filename=['/scratch/la67/${PROJECT}/finalstats/${PROJECT}_' num2str(RUNS(iRun)) '.mat'];
	save(filename,'mbag','modelsummary');
end
EOF

cd /home/la67/${PROJECT}/scripts
