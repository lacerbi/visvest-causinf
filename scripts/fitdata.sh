#!/bin/sh

#PBS -o localhost:${PBS_O_WORKDIR}/
#PBS -e localhost:${PBS_O_WORKDIR}/

module purge
#. /etc/profile.d/modules.sh

# Use Intel compiler
module load matlab
export MATLABPATH=${MATLABPATH}:${HOME}/${PROJECT}:${HOME}/MATLAB
source ${HOME}/MATLAB/setpath.sh

#Check if running as an array job
if [[ ! -z "$PBS_ARRAYID" ]]; then
	IID=${PBS_ARRAYID}
fi
if [[ ! -z "$SGE_TASK_ID" ]]; then
        IID=${SGE_TASK_ID}
fi
if [[ ! -z "$SLURM_ARRAY_TASK_ID" ]]; then
        IID=${SLURM_ARRAY_TASK_ID}
fi


# Run the program
echo ${WORKDIR} ${PROJECT} ${IID}.job

cat<<EOF | matlab -nodisplay
addpath(genpath('${HOME}/MATLAB'));
addpath(genpath('${HOME}/${PROJECT}'));
cd('${WORKDIR}');
nsamples=${NSAMPLES}; % MCMC samples
nburnin=${NBURNIN}; % MCMC burn-in
storedsamples=${STOREDSAMPLES}; % Stored MCMC samples
optfevals=${OPTFEVALS}; % Optimization function evaluations
continueflag=${CONTINUE}; % Continue flag
loadmbag=${LOADMBAG};
fitstep=[];
datafile=['${DATAFILENAME}'];
if ~isempty(datafile)
	temp=load(datafile); data=temp.data;
else
	data=[];
end
if ~isempty(optfevals) && optfevals==0 && nsamples==0
	recomputesamplingmetrics=1;
	computemarginallike=1;
	refit=1;
else
	recomputesamplingmetrics=0;
	computemarginallike=0;
	refit=0;
end
if loadmbag
	temp = load('../${PROJECT}_${RUN}.mat','mbag');
	mbag = temp.mbag
	newsampling=[${NEWSAMPLING}];
	if isempty(newsampling); newsampling = 0; end
	if newsampling
		display('Ignoring previous samples; sampling from scratch.');
		for i=1:numel(mbag.bag)
			mbag.bag{i}.sampling = [];
		end
		fitstep=1;
	end
else
	mbag = [];
end
ModelWork_batchEval('${PROJECT}',data,'${IID}.job','procid',${IID},'optfevals',optfevals,'nsamples',nsamples,'burnin',nburnin,'maxstoredsamples',storedsamples,'continueflag',continueflag,'mbag',mbag,'recomputesamplingmetrics',recomputesamplingmetrics,'computemarginallike',computemarginallike,'refit',refit,'fitstep',fitstep);
EOF
