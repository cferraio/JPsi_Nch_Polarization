#!/bin/sh
source /cvmfs/sft.cern.ch/lcg/app/releases/ROOT/5.34.05/x86_64-slc5-gcc43-opt/root/bin/thisroot.sh                                     
source /cvmfs/sft.cern.ch/lcg/app/releases/ROOT/5.34.05/x86_64-slc5-gcc43-opt/root/bin/setxrd.sh /cvmfs/sft.cern.ch/lcg/external/xrootd/3.2.4/x86_64-slc5-gcc46-opt/


homedir=$PWD
cd ..
cd ..
basedir=$PWD
cd macros/polFit
storagedir=`more storagedir` #please define the directory storagedir in the file macros/polFit/storagedir
#storagedir=${basedir}/Psi/Data

########## INPUTS ##########

for nState in 4 5;do

cp ../../interface/rootIncludes.inc               rootIncludes.inc
cp ../../interface/commonVar_Psi$[nState-3]S.h    commonVar.h
cp ../../interface/ToyMC_Psi$[nState-3]S.h        ToyMC.h

SystID=Framework

nSystematics=3

JobID1=FrameworkI
JobID2=FrameworkII
JobID3=FrameworkIII
JobID4=
JobID5=
JobID6=
JobID7=
JobID8=
JobID9=


#SystID=TotalSyst
#
#nSystematics=7
#
#JobID1=BestSyst_Bkg
#JobID2=BestSyst_FrameworkI
#JobID3=BestSyst_Param
#JobID4=BestSyst_Sig_NoUnpol
#JobID5=BestSyst_TnP
#JobID6=ConstSyst
#JobID7=BestSyst_28p_SQRT12
#JobID8=
#JobID9=

if [ $nState -eq 4 ]
then
ptBinMin=1
ptBinMax=12
fi 
if [ $nState -eq 5 ]
then
ptBinMin=1
ptBinMax=5
fi

########################################

cd ${homedir}

touch AverageSystematics.cc
make

mkdir Systematics
mkdir Systematics/${SystID}

SystDir=Systematics/${SystID}/AverageSyst

mkdir ${SystDir}

./AverageSystematics ${JobID1}=JobID1 ${JobID2}=JobID2 ${JobID3}=JobID3 ${JobID4}=JobID4 ${JobID5}=JobID5 ${JobID6}=JobID6 ${JobID7}=JobID7 ${JobID8}=JobID8 ${JobID9}=JobID9 ${SystID}=SystID ${storagedir}=storagedir ${basedir}=basedir ${ptBinMin}ptBinMin ${ptBinMax}ptBinMax ${nState}nState ${nSystematics}nSystematics

rm AverageSystematics
done


