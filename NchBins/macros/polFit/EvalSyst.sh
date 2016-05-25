#!/bin/sh

source /cvmfs/sft.cern.ch/lcg/app/releases/ROOT/5.34.05/x86_64-slc5-gcc43-opt/root/bin/thisroot.sh                                     
source /cvmfs/sft.cern.ch/lcg/app/releases/ROOT/5.34.05/x86_64-slc5-gcc43-opt/root/bin/setxrd.sh /cvmfs/sft.cern.ch/lcg/external/xrootd/3.2.4/x86_64-slc5-gcc46-opt/

homedir=$PWD
cd ../..
basedir=$PWD
cd macros/polFit
storagedir=`more storagedir`
#storagedir=$basedir/Psi/ToyMC

########## INPUTS ##########

for nState in 4;do


#this is the systematic error!!! 
#ifCorrectCentralResultsForBias: here is the central
for JobID1 in 19May16_MassUpdateFixedErrBars_FracL75; do

#this is the default!!!
#JobID2=Psi$[nState-3]S_fLSB0_18March2013 #default_Linlin_20March2013   
JobID2=19May16_MassUpdateFixedErrBars_FracL25

#define name of the directory
SystID=BackgroundModel_19May16

if [ $nState -eq 4 ]
then
cpmBinMin=1
cpmBinMax=10
ptBinMin=1
ptBinMax=2
rapBinMin=1
rapBinMax=1
fi
if [ $nState -eq 5 ]
then
cpmBinMin=1
cpmBinMax=5
ptBinMin=1
ptBinMax=5
rapBinMin=1
rapBinMax=3
fi

statErrConsideration=false
centralsPlusSyst=false #take one from centrals and one from systematics
differentErrors=false #'take central value from JobID2, take error from JobID1'
sqrt12=false #divide mean value by sqrt12
removeFirstPoint=false
TU=false #calculate mean = sqrt(error1^2 - error2^2)
########################################

cd ${homedir}

#create directory
touch EvaluateSyst.cc
make

mkdir Systematics
mkdir Systematics/${SystID}

SystDir=Systematics/${SystID}/${JobID1}_TO_${JobID2}

mkdir ${SystDir}
#copy files
#cp MakefileEvalSyst ${SystDir}/Makefile
#cp ../../interface/rootIncludes.inc ${SystDir}/rootIncludes.inc
#cp ../../interface/commonVar_Psi$[nState-3]S.h ${SystDir}/commonVar.h
#cp EvaluateSyst.cc ${SystDir}/EvaluateSyst.cc
#cp ../../interface/ToyMC_Psi$[nState-3]S.h ${SystDir}/ToyMC.h
#cd ${SystDir}
#touch EvaluateSyst.cc
#make

./EvaluateSyst JobID1=${JobID1} JobID2=${JobID2} SystDir=${SystDir} storagedir=${storagedir} basedir=${basedir} cpmBinMin=${cpmBinMin} cpmBinMax=${cpmBinMax} ptBinMin=${ptBinMin} ptBinMax=${ptBinMax} rapBinMin=${rapBinMin} rapBinMax=${rapBinMax} nState=${nState} statErrConsideration=${statErrConsideration} centralsPlusSyst=${centralsPlusSyst} differentErrors=${differentErrors} sqrt12=${sqrt12} removeFirstPoint=${removeFirstPoint} TU=${TU}
rm EvaluateSyst
done
done

