#!/bin/bash
PROJECT="VestBMS"
cd /scratch/la67/${PROJECT}

ALLDIRS=$(ls -d */)
for STATSDIR in ${ALLDIRS}; do

	echo "cd /scratch/la67/${PROJECT}/${STATSDIR}"
	cd /scratch/la67/${PROJECT}/${STATSDIR}

	rm *.log
	rm *.job
	rm VB*

	MAINDIR=$(ls -d */)
	cd ${MAINDIR}

	if [ ! -z "${MAINDIR}" ]; then

		SUBDIRS=$(ls -d */)
		for DIR in ${SUBDIRS}; do
			echo "${STATSDIR}${MAINDIR}${DIR}"
			cd ${DIR};
			mkdir tmp;
			ls *.tmp;
			mv *.tmp tmp/;
			rm -r tmp/;
			cd ..;
		done
	fi
done

cd /home/la67/${PROJECT}/scripts

