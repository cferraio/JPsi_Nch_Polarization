#!/bin/sh
source /afs/ihep.ac.cn/users/z/zhangll/fs/rootset.sh 

homedir=$PWD
cd ..
cd ..
basedir=$PWD
cd macros/polFit
#storagedir=`more storagedir`/Data #please define the directory storagedir in the file macros/polFit/storagedir
storagedir=$basedir/Psi/Data

########## INPUTS ##########

for nState in 4;do

cp ../../interface/rootIncludes.inc               rootIncludes.inc
cp ../../interface/commonVar_Psi$[nState-3]S.h    commonVar.h
cp ../../interface/ToyMC_Psi$[nState-3]S.h        ToyMC.h

#touch AlterPPD.cc
#make

JobID=AlteredPPD_BKGlinPLUSRestSquaredGauss_5nRand_Apr5

DefaultID=Psi$[nState-3]S_ctauScen0_FracLSB-1_19Mar2013_1FracBg_1BgModel_1Rho
ShiftID=BiasCorrection_1S2S3SAug12_1Sig
ShiftResults=0

nSystematics=2

#ProbDist 1...gauss, 2...linear

SystID1Base=BackgroundModel
SystID1Specify=fLSB100_TO_fLSB0 #in the cc file, this values are multiplied by sqrt12
SystID1Title=BKGmodel
SystID1ProbDist=2

SystID2Base=TotalSyst
SystID2Specify=SquaredSystApr4_ExceptBKGmodel
SystID2Title=TotalSyst_ExceptBKGmodel
SystID2ProbDist=1

if [ $nState -eq 4 ]
then 
rapBinMin=1
rapBinMax=2
ptBinMin=1
ptBinMax=12
fi 
if [ $nState -eq 5 ]
then 
rapBinMin=1
rapBinMax=3
ptBinMin=1
ptBinMax=5
fi

#SystID1Base=TheGreatRun_BKGmodel
#SystID1Specify=BestSyst_28p_SQRT12 #in the cc file, this values are multiplied by sqrt12
#SystID1Title=BKGmodel
#SystID1ProbDist=2
#
#SystID2Base=TotalSyst
#SystID2Specify=TotalSquaredSystMay21_ExceptBKGmodel
#SystID2Title=TotalSyst_ExceptBKGmodel
#SystID2ProbDist=1

#SystID2Base=TheGreatRun_Rho
#SystID2Specify=ConstSyst
#SystID2Title=RhoFactor
#SystID2ProbDist=1

#SystID3Base=TheGreatRun_Param
#SystID3Specify=BestSyst
#SystID3Title=Parametrization
#SystID3ProbDist=1
#
#SystID4Base=TheGreatRun_TnP
#SystID4Specify=BestSyst
#SystID4Title=TnP_model
#SystID4ProbDist=1
#
#SystID5Base=TheGreatRun_FrameworkII
#SystID5Specify=BestSyst_Bkg
#SystID5Title=FrameworkIII
#SystID5ProbDist=1
#
#SystID6Base=TheGreatRun_FrameworkII
#SystID6Specify=BestSyst_Sig_NoUnpol
#SystID6Title=FrameworkII
#SystID6ProbDist=1
#
#SystID7Base=TheGreatRun_FrameworkI
#SystID7Specify=BestSyst
#SystID7Title=FrameworkI
#SystID7ProbDist=1

########################################

cd ${homedir}

touch AlterPPD.cc
make

JobIDDir=${storagedir}/${DefaultID}
JobIDDirAlterPPD=${storagedir}/${DefaultID}_${JobID}

mkdir ${JobIDDirAlterPPD}

rap_=${rapBinMin}
while [ $rap_ -le ${rapBinMax} ]
do
pT_=${ptBinMin}
while [ $pT_ -le ${ptBinMax} ]
do

rm ${JobIDDirAlterPPD}/results_Psi$[nState-3]S_rap${rap_}_pT${pT_}.root
cp ${JobIDDir}/results_Psi$[nState-3]S_rap${rap_}_pT${pT_}.root ${JobIDDirAlterPPD}/results_Psi$[nState-3]S_rap${rap_}_pT${pT_}.root

./AlterPPD ${DefaultID}=DefaultID ${ShiftID}=ShiftID ${JobID}=JobID ${SystID1Base}=SystID1Base ${SystID1Specify}=SystID1Specify ${SystID1Title}=SystID1Title ${SystID2Base}=SystID2Base ${SystID2Specify}=SystID2Specify ${SystID2Title}=SystID2Title ${SystID3Base}=SystID3Base ${SystID3Specify}=SystID3Specify ${SystID3Title}=SystID3Title ${SystID4Base}=SystID4Base ${SystID4Specify}=SystID4Specify ${SystID4Title}=SystID4Title ${SystID5Base}=SystID5Base ${SystID5Specify}=SystID5Specify ${SystID5Title}=SystID5Title ${SystID6Base}=SystID6Base ${SystID6Specify}=SystID6Specify ${SystID6Title}=SystID6Title ${SystID7Base}=SystID7Base ${SystID7Specify}=SystID7Specify ${SystID7Title}=SystID7Title ${SystID8Base}=SystID8Base ${SystID8Specify}=SystID8Specify ${SystID8Title}=SystID8Title ${SystID1ProbDist}SystID1ProbDist ${SystID2ProbDist}SystID2ProbDist ${SystID3ProbDist}SystID3ProbDist ${SystID4ProbDist}SystID4ProbDist ${SystID5ProbDist}SystID5ProbDist ${SystID6ProbDist}SystID6ProbDist ${SystID7ProbDist}SystID7ProbDist ${SystID8ProbDist}SystID8ProbDist ${basedir}=basedir ${storagedir}=storagedir ${pT_}ptBinMin ${pT_}ptBinMax ${rap_}rapBinMin ${rap_}rapBinMax ${nSystematics}nSystematics ${nState}nState ShiftResults=${ShiftResults} 


pT_=$[pT_+1]
done
rap_=$[rap_+1]
done

done
