#!/bin/sh

export VO_CMS_SW_DIR=/sharesoft/cmssw #comment out for non-condor
. $VO_CMS_SW_DIR/cmsset_default.sh #comment out for non-condor
cd /home/ferraioc/PolNew/CMSSW_5_3_20/src/ChicPol/macros/polFit #comment out for non-condor
eval `scramv1 runtime -sh` #comment out for non-condor

source /cvmfs/sft.cern.ch/lcg/app/releases/ROOT/5.34.05/x86_64-slc5-gcc43-opt/root/bin/thisroot.sh                                     
source /cvmfs/sft.cern.ch/lcg/app/releases/ROOT/5.34.05/x86_64-slc5-gcc43-opt/root/bin/setxrd.sh /cvmfs/sft.cern.ch/lcg/external/xrootd/3.2.4/x86_64-slc5-gcc46-opt/

CONDOR_JOB="runfitscondor.jdl"
n=0

storagedir=/data/users/ferraioc/Polarization/JPsi/Data
homedir=/home/ferraioc/PolNew/CMSSW_5_3_20/src/JPsi_Nch_Polarization/NchInt/macros/polFit



########################
########################
########################
JobID=NonPrompt2011
DataID=_ctauScen0_FracLSB-1_NonPrompt2011
#JobID=2012Data_2011effs
#DataID=_ctauScen0_FracLSB-1_2012Run3


nState=4
nSample=50000
nfits=5
FidCuts=11
########################
########################
########################

cd ${homedir}
cd ..
cd ..
basedir=$PWD
cd macros/polFit

datadir=${basedir}/macros/DataFiles/SetOfCuts${FidCuts}${DataID}/Psi$[nState-3]S/tmpFiles

TreeID=Psi$[nState-3]S

cd ${homedir}

mkdir -p ${storagedir}/${JobID}

cp ../../interface/rootIncludes.inc ${storagedir}/${JobID}/rootIncludes.inc
cp ../../interface/commonVar_Psi$[nState-3]S.h ${storagedir}/${JobID}/commonVar.h
cp ../../interface/ToyMC_Psi$[nState-3]S.h ${storagedir}/${JobID}/ToyMC.h
cp ../../interface/effsAndCuts_Psi$[nState-3]S.h ${storagedir}/${JobID}/effsAndCuts.h

cp ${basedir}/macros/polFit/polGenRecFitPlot.cc ${storagedir}/${JobID}/polGenRecFitPlot.cc
cp ${basedir}/macros/polFit/polRapPtPlot.cc ${storagedir}/${JobID}/polRapPtPlot.cc
cp ${basedir}/macros/polFit/PlotFinalResults.cc ${storagedir}/${JobID}/PlotFinalResults.cc
cp ${basedir}/macros/polFit/Makefile ${storagedir}/${JobID}/Makefile
cp ${basedir}/macros/polFit/polGen.C ${storagedir}/${JobID}/polGen.C
cp ${basedir}/macros/polFit/polRec.C ${storagedir}/${JobID}/polRec.C
cp ${basedir}/macros/polFit/polFit.C ${storagedir}/${JobID}/polFit.C
cp ${basedir}/macros/polFit/polPlot.C ${storagedir}/${JobID}/polPlot.C



cd ${storagedir}/${JobID}
cp ${basedir}/macros/polFit/runDataFits.sh .

touch polGenRecFitPlot.cc
make

cd ${homedir}

while [[ $n -le $nfits-1 ]]
do
for rap in 1 2 #3
do
	for pt in 2 3 4 5 6 7 8 9 10 11 12
	do
		cp runcondorFits.jdl $CONDOR_JOB
		echo "$nSample $DataID $storagedir $homedir $n $rap $pt $nState $JobID  " >> $CONDOR_JOB
		echo "Queue 1" >> $CONDOR_JOB	
		condor_submit $CONDOR_JOB
done
done
n=$(( n+1 ))
done