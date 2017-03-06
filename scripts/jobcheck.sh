#!/bin/bash

PROJECT="VestBMS"

# Initialize variables
LISTD=
LISTS=
LISTF=
LISTR=
FIRSTD=1
FIRSTS=1
FIRSTF=1
FIRSTR=1
NSUCC=0
NRUN=0
IDS=$(seq 1 $2)

# Start scan
cd ${SCRATCH}/${PROJECT}/run$1
for i in $IDS; do

        # Check if the job has started
        if [ ! -f diary-$i.log ]
        then
                if [ $FIRSTD -eq 1 ]
                then
                        LISTD=$i
                        FIRSTD=0
                else
                        LISTD=$LISTD,$i
                fi
        fi

        # Check if the job has failed
        if [ -f fail-$i.log ]
        then
                if [ $FIRSTF -eq 1 ]
                then
                        LISTF=$i
                        FIRSTF=0
                else
                        LISTF=$LISTF,$i
                fi
        fi

        # Check completed jobs
        if [ -f success-$i.log ]
        then
                if [ $FIRSTS -eq 1 ]
                then
                        LISTS=$i
                        FIRSTS=0
                else
                        LISTS=$LISTS,$i
                fi
		let NSUCC="$NSUCC + 1"
        fi

	# Check still running jobs
        if [ ! -f fail-$i.log ] && [ ! -f success-$i.log ]
        then
                if [ $FIRSTR -eq 1 ]
                then
                        LISTR=$i
                        FIRSTR=0
                else
                        LISTR=$LISTR,$i
                fi
		let NRUN="$NRUN + 1"
        fi



done

# Print report

echo
echo "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
echo
echo "Report for RUN${1}:"
echo
if [ $FIRSTD -eq 1 ]
then
        echo "All jobs have started."
else
        echo "Jobs that did NOT start:"
        echo $LISTD
fi
echo
if [ $FIRSTF -eq 1 ]
then
        echo "No jobs have failed (so far)."
else
        echo "Failed jobs:"
        echo $LISTF
fi
echo
if [ $FIRSTS -eq 1 ]
then
        echo "No jobs have completed successfully (so far?)."
else
        echo "Completed jobs ($NSUCC out of $2):"
        echo $LISTS
fi
echo
if [ $FIRSTR -eq 1 ]
then
        echo "No jobs are still running."
else
        echo "Unfinished jobs ($NRUN out of $2):"
        echo $LISTR
fi
echo
echo "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
echo

cd ..
