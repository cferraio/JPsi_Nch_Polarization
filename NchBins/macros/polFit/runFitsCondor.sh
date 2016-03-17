#!/bin/sh


#pt=1
#cpm=1
CONDOR_JOB="runfitscondor.jdl"

for pt in 1 2
do
	for cpm in 1 2 3 4 5 6
	do
		cp runcondorFits.jdl $CONDOR_JOB
		echo "condor-simple.py $pt $cpm" >> $CONDOR_JOB
		echo "Queue 1" >> $CONDOR_JOB	
		condor_submit $CONDOR_JOB
done
done