#!/bin/bash
PROJECT="VestBMS"
STATSDIR=${1}
cd /scratch/la67/${PROJECT}/stats${STATSDIR}

#Get list of all log likelihoods (old and new)
LL_LIST="$(cat VB${STATSDIR}.o* | grep hood:)"
LL_LIST="$(echo "${LL_LIST}" | rev | cut -d" " -f1 | cut -c 2- | rev | tr "\n" " ")"

FIRST=1
for LL in ${LL_LIST}; do
	if [ "$FIRST" -eq "1" ]; then
		LL_OLD=${LL}
		FIRST=0
	else
		LL_DIFF=$(echo "${LL} - ${LL_OLD}" | bc)
		echo ${LL_DIFF}
		FIRST=1
	fi
done

cd /home/la67/${PROJECT}/scripts

