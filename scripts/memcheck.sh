#!/bin/bash
USERID=${USER}

#Get list of memory usage in all running jobs (expand job arrays)
LIST="$(qstat -u ${USERID} -f -t | grep resources_used.mem)"
LIST="$(echo "${LIST}" | rev | cut -d" " -f1 | rev | tr -d "." | tr "\n" " ")"

#Find maximum memory usage
MAXMEM=0
for IMEM in ${LIST}; do
	if [ "${IMEM%??}" -gt "$MAXMEM" ]; then
		MAXMEM=${IMEM%??}
	fi
done
let MAXMEM="${MAXMEM} / 1024"

echo "Maximum memory usage across all running jobs for user ${USERID} is ${MAXMEM}MB."
