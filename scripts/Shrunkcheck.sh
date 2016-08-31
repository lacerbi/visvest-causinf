#!/bin/sh
echo "Usage: Shrunkcheck job#"

PROJECT="VestBMS"
cd /scratch/la67/${PROJECT}

module purge
#. /etc/profile.d/modules.sh

BASEDIR="run${1}"
cd ${BASEDIR}

FIRST=1
FILES=$(ls VB*.o*)
for i in $FILES; do
	OUT=$(cat ${i} | grep Shrunk | wc | awk '{print $1;}')

	if [ $OUT -gt 0 ]; then
		if [ ${FIRST} -gt 0 ]; then
			LIST="${i#*-}"
			FIRST=0
		else
			LIST="${LIST},${i#*-}"
		fi
	fi
done

echo ${LIST}


cd ..
cd /home/la67/${PROJECT}/scripts
