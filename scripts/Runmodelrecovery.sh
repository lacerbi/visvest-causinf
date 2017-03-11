#!/bin/bash
PROJECT="VestBMS"
SHORTNAME="VB"
BASEDIR="${HOME}/${PROJECT}"
SOURCEDIR="${BASEDIR}"
JOBSCRIPT="${BASEDIR}/scripts/fitdata.sh"

#Job parameters
RUN=${1}
WORKDIR="${SCRATCH}/${PROJECT}/run${RUN}"
mkdir ${WORKDIR}
cd ${WORKDIR}
MAXID=$(ls -l *.job | wc -l)
MAXRT=NaN
VERBOSE=1
NSAMPLES=0
NBURNIN=0
STOREDSAMPLES=5000
OPTFEVALS="[]"
CONTINUE=1
LOADMBAG=0
NEWSAMPLING=0
#DATAFILENAME="VestBMS_fakedata_unity.mat"
DATAFILENAME="VestBMS_fakedata_joint.mat"

RUNTIME=24:00:00
NODES="1"
PPN="1"
MEM="3GB"
RESOURCES="nodes=${NODES}:ppn=${PPN},mem=${MEM},walltime=${RUNTIME}"


#Job list is second argument
if [[ ! -z "$2" ]]; then
        JOBLIST=$2
else
	#Get all files in directory
	JOBLIST=$(ls -l *.job | rev | cut -f 1 -d " " | rev | cut -f 1 -d ".")
	JOBLIST=$(echo $JOBLIST)
	JOBLIST=$(echo ${JOBLIST// /,})
fi


#if [[ -z ${1} ]]; then
#        JOBLIST="1-$MAXID"
#        NEWJOB=1
#else
#        JOB=${1}
#        NEWJOB=0
#        echo "JOB=${JOB}" >> ${BASEDIR}/reruns.log
#fi

#Convert from spaces to commas
JOBLIST=${JOBLIST// /,}
echo JOBS $JOBLIST

JOBNAME=${SHORTNAME}${RUN}

if [ ${CLUSTER} = "Prince" ]; then
        # running on Prince
        sbatch --verbose --array=${JOBLIST} --mail-type=FAIL --mail-user=${USER}@nyu.edu --mem=${MEM} --time=${RUNTIME} --nodes=${NODES} --ntasks-per-node=${PPN} --export=PROJECT=${PROJECT},RUN=${RUN},MAXID=$MAXID,WORKDIR=$WORKDIR,USER=$USER,VERBOSE=${VERBOSE},OPTFEVALS=${OPTFEVALS},NSAMPLES=${NSAMPLES},NBURNIN=${NBURNIN},STOREDSAMPLES=${STOREDSAMPLES},NEWSAMPLING=${NEWSAMPLING},CONTINUE=${CONTINUE},LOADMBAG=${LOADMBAG},DATAFILENAME=${DATAFILENAME} --job-name=${JOBNAME} ${JOBSCRIPT}
else
	qsub -t ${JOBLIST} -q normal -v PROJECT=${PROJECT},RUN=${RUN},MAXID=$MAXID,MAXRT=$MAXRT,WORKDIR=$WORKDIR,USER=$USER,VERBOSE=${VERBOSE},OPTFEVALS=${OPTFEVALS},NSAMPLES=${NSAMPLES},NBURNIN=${NBURNIN},STOREDSAMPLES=${STOREDSAMPLES},NEWSAMPLING=${NEWSAMPLING},CONTINUE=${CONTINUE},LOADMBAG=${LOADMBAG},DATAFILENAME=${DATAFILENAME} -l ${RESOURCES} -M ${USER}@nyu.edu -N ${JOBNAME} ${JOBSCRIPT}
fi
