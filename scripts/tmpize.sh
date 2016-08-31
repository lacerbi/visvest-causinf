#!/bin/bash
PROJECT="VestBMS"
STATSDIR=${1}
cd /scratch/la67/${PROJECT}/run${STATSDIR}

MAINDIR=$(ls -d */)
cd ${MAINDIR}

SUBDIRS=$(ls -d */)
for DIR in ${SUBDIRS}; do
	echo "run${STATSDIR}/${MAINDIR}${DIR}"
	cd ${DIR};
	mkdir tmp;
	ls *.tmp;
	mv *.tmp tmp/;
	cd ..;
done

cd /home/la67/${PROJECT}/scripts

