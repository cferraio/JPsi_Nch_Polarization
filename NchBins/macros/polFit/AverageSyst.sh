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

for nState in 4;do

cp ../../interface/rootIncludes.inc               rootIncludes.inc
cp ../../interface/commonVar_Psi$[nState-3]S.h    commonVar.h
cp ../../interface/ToyMC_Psi$[nState-3]S.h        ToyMC.h

SystID=FrameworkII_19May2016

nSystematics=4

JobID1=Sig_frame3scen1_Bkg_frame1scen3
JobID2=Sig_frame3scen2_Bkg_frame1scen3
JobID3=Sig_frame3scen4_Bkg_frame1scen3
JobID4=Sig_frame3scen5_Bkg_frame1scen3
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
ptBinMax=2
cpmBinMin=1
cpmBinMax=10
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

./AverageSystematics ${JobID1}=JobID1 ${JobID2}=JobID2 ${JobID3}=JobID3 ${JobID4}=JobID4 ${JobID5}=JobID5 ${JobID6}=JobID6 ${JobID7}=JobID7 ${JobID8}=JobID8 ${JobID9}=JobID9 ${SystID}=SystID ${storagedir}=storagedir ${basedir}=basedir ${ptBinMin}ptBinMin ${ptBinMax}ptBinMax ${cpmBinMin}cpmBinMin ${cpmBinMax}cpmBinMax ${nState}nState ${nSystematics}nSystematics

rm AverageSystematics
done


