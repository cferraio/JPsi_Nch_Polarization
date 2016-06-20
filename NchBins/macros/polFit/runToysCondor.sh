#!/bin/sh


#pt=1
#cpm=1
CONDOR_JOB="runtoyscondor.jdl"
n=0
ngen=50

lessthen=$ngen-1
while [[ $n -le $ngen-1 ]]
do
for polsig in 1 2 4 5
do
#for pt in 2
for pt in 1 2
do
	for cpm in 1 2 3 4 5 6 7 8 9 10
#	for cpm in 10	
	do
		sleep 1
		cp runcondorToys.jdl $CONDOR_JOB
		echo "$pt $cpm $n $polsig" >> $CONDOR_JOB 
		echo "Queue 1" >> $CONDOR_JOB	
		condor_submit $CONDOR_JOB
done
done
done 
n=$(( n+1 ))
done