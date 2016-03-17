#!/bin/sh

########## INPUTS ##########

nState=1

JobID=Toy_TheGreatRun_date_1S_Sanity_TnP_TighterCuts

nGenerations=50
MergeFiles=1

rapBinMin=1
rapBinMax=1
ptBinMin=1
ptBinMax=3

polScenSig=3
polScenBkg=3
frameSig=1
frameBkg=1

MPValgo=3 		#1...mean,2...gauss,3...gauss-loop with chi2<2
nSigma=2


########################################

gen=false
rec=false
fit=false
plot=true

NewAccCalc=true
nEff=109
UseMCeff=true
nDileptonEff=1
UseMCDileptoneff=true
nRhoFactor=1

FidCuts=9

nSample=10000
nSkipGen=0

#GENERATION SETTINGS
ConstEvents=10000
UseConstEv=true

UseDifferingEff=false
nEffRec=104
UseMCReceff=false
nDileptonEffRec=201
UseMCDileptonReceff=true
nRecRhoFactor=322

TreeID=${nState}SUps

##############################################

cd ..
cd ..
basedir=$PWD
cd macros/polFit
storagedir=`more storagedir`/ToyMC #please define the directory storagedir in the file macros/polFit/storagedir



touch polGenRecFitPlot.cc
make

ScenDir=Sig_frame${frameSig}scen${polScenSig}_Bkg_frame${frameBkg}scen${polScenBkg}

mkdir ${storagedir}
mkdir ${storagedir}/${JobID}
mkdir ${storagedir}/${JobID}/${ScenDir}
cp ${basedir}/macros/polFit/polGenRecFitPlot ${storagedir}/${JobID}/${ScenDir}/polGenRecFitPlot
cp ${basedir}/macros/polFit/polGen.C ${storagedir}/${JobID}/${ScenDir}/polGen.C
cp ${basedir}/macros/polFit/polRec.C ${storagedir}/${JobID}/${ScenDir}/polRec.C
cp ${basedir}/macros/polFit/polFit.C ${storagedir}/${JobID}/${ScenDir}/polFit.C
cp ${basedir}/macros/polFit/polPlot.C ${storagedir}/${JobID}/${ScenDir}/polPlot.C

cd ${storagedir}/${JobID}/${ScenDir}

rap_=${rapBinMin}
while [ $rap_ -le ${rapBinMax} ]
do
pT_=${ptBinMin}
while [ $pT_ -le ${ptBinMax} ]
do

nGen_=${nSkipGen}
nGen_=$[nGen_+1]
while [ $nGen_ -le ${nGenerations} ]
do


resultfilename=resultsMerged_${nState}SUps_rap${rap_}_pT${pT_}.root


GenResultName=${storagedir}/${JobID}/${ScenDir}/rap${rap_}_pT${pT_}/Generation${nGen_}/results.root

if [ $MergeFiles -eq 1 ]
then

if [ $nGen_ -eq 1 ]
then
cp ${GenResultName} ${resultfilename}
fi

if [ $nGen_ -ge 2 ]
then
mv ${resultfilename} BUFFER_${resultfilename}
hadd -f ${resultfilename} BUFFER_${resultfilename} ${GenResultName}
rm BUFFER_${resultfilename}
fi

fi



nGen_=$[nGen_+1]
done

if [ $MergeFiles -eq 1 ]
then
mv ${resultfilename} ${storagedir}/${JobID}/results_${nState}SUps_rap${rap_}_pT${pT_}.root
fi

cp ${storagedir}/${JobID}/${ScenDir}/polGenRecFitPlot ${storagedir}/${JobID}/${ScenDir}/polGenRecFitPlot_rap${rap_}_pt${pT_}
./polGenRecFitPlot_rap${rap_}_pt${pT_} ${nGen_}ThisGen ${JobID}=JobID/${ScenDir} ${storagedir}=storagedir ${basedir}=basedir ${nGenerations}nGenerations ${polScenSig}polScenSig ${frameSig}frameSig ${polScenBkg}polScenBkg ${frameBkg}frameBkg ${rap_}rapBinMin ${rap_}rapBinMax ${pT_}ptBinMin ${pT_}ptBinMax ${nEff}nEff ${nDileptonEff}nDiEff ${FidCuts}FidCuts ${nSample}nSample ${ConstEvents}ConstEvents ${nSkipGen}nSkipGen UseConstEv=${UseConstEv} gen=false rec=false fit=false plot=${plot} ${TreeID}=TreeID ${datadir}=realdatadir UseMCeff=${UseMCeff} UseMCDileptoneff=${UseMCDileptoneff} ${nRhoFactor}nRhoFactor ${MPValgo}MPValgo scalePlots=true NewAccCalc=${NewAccCalc}
rm polGenRecFitPlot_rap${rap_}_pt${pT_}


pT_=$[pT_+1]
done
rap_=$[rap_+1]
done


