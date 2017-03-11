#!/bin/bash
PROJECT=VestBMS
SHORTNAME=VB
BASEDIR="/home/la67/${PROJECT}"
SOURCEDIR="${BASEDIR}"
JOBSCRIPT="${BASEDIR}/scripts/fitdata.sh"

#Job parameters
RUN=${1}
WORKDIR="/scratch/la67/${PROJECT}/run${RUN}"
mkdir ${WORKDIR}
cd ${WORKDIR}
MAXID=$(ls -l *.job | wc -l)
RUNTIME=36:00:00
MAXRT=NaN
VERBOSE=1
NSAMPLES=0
NBURNIN=0
STOREDSAMPLES=5000
OPTFEVALS="[]"
CONTINUE=1
LOADMBAG=0
NEWSAMPLING=0

#Job list is second argument
if [[ ! -z "$2" ]]; then
        JOBLIST=$2
else
	#Get all files in directory
	JOBLIST=$(ls -l *.job | rev | cut -f 1 -d " " | rev | cut -f 1 -d ".")
	JOBLIST=$(echo $JOBLIST)
	JOBLIST=$(echo ${JOBLIST// /,})
fi

#RESOURCES="nodes=1:ppn=1,mem=4GB,walltime=${RUNTIME},feature=ivybridge_20p_64GB_3000MHz"
RESOURCES="nodes=1:ppn=1,mem=2GB,walltime=${RUNTIME}"

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
qsub -t ${JOBLIST} -v PROJECT=${PROJECT},RUN=${RUN},MAXID=$MAXID,MAXRT=$MAXRT,WORKDIR=$WORKDIR,VERBOSE=${VERBOSE},OPTFEVALS=${OPTFEVALS},NSAMPLES=${NSAMPLES},NBURNIN=${NBURNIN},STOREDSAMPLES=${STOREDSAMPLES},NEWSAMPLING=${NEWSAMPLING},CONTINUE=${CONTINUE},LOADMBAG=${LOADMBAG} -l ${RESOURCES} -N ${JOBNAME} ${JOBSCRIPT}
