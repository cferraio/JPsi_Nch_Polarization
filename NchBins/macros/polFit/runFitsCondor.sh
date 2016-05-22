#!/bin/sh


#pt=1
#cpm=1
CONDOR_JOB="runfitscondor.jdl"
n=0
nfits=50

#lessthen=$nfits-1
while [[ $n -le $nfits-1 ]]

do
#for pt in 2
for pt in 1 2
do
	for cpm in 1 2 3 4 5 6 7 8 9 10
#	for cpm in 10	
	do
		cp runcondorFits.jdl $CONDOR_JOB
		echo "$n condor-simple.py $pt $cpm" >> $CONDOR_JOB
		echo "Queue 1" >> $CONDOR_JOB	
		condor_submit $CONDOR_JOB
done
done
n=$(( n+1 ))
done