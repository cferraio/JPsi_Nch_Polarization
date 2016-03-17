#!/bin/sh

export VO_CMS_SW_DIR=/sharesoft/cmssw #comment out for non-condor
. $VO_CMS_SW_DIR/cmsset_default.sh #comment out for non-condor
cd /home/ferraioc/Polarization/CMSSW_4_4_2/src/JPsi/NchBins/macros/polFit #comment out for non-condor
eval `scramv1 runtime -sh` #comment out for non-condor


source /cvmfs/sft.cern.ch/lcg/app/releases/ROOT/5.34.05/x86_64-slc5-gcc43-opt/root/bin/thisroot.sh                                     
source /cvmfs/sft.cern.ch/lcg/app/releases/ROOT/5.34.05/x86_64-slc5-gcc43-opt/root/bin/setxrd.sh /cvmfs/sft.cern.ch/lcg/external/xrootd/3.2.4/x86_64-slc5-gcc46-opt/
#### following will be changed by "launchFit.sh" ####
# HOMEDIR
# NState 
# RapBinMin
# RapBinMax
# PtBinMin
# PtBinMax
# NSkipGen
######################################### 

storagedir=$3
HOMEDIR=$4
setptbin=$5
setcpmbin=$6
setbgframe=$7

echo "storage dir: $storagedir"
echo "home dir: $homedir"
echo "ptbin: $setptbin"
echo "cpmbin: $setcpmbin"


#homedir=HOMEDIR
cd ${homedir}
cd ..
cd ..
basedir=$PWD
cd macros/polFit
#storagedir=`more storagedir`/Data #please define the directory storagedir in the file macros/polFit/storagedir
storagedir=/data/users/ferraioc/Polarization/JPsi/NchBins
datadir_Start=${basedir}/macros/DataFiles

########## INPUTS ##########

#Batch submission system: 0/1
useBatch=1

fracL=50 #in percent #MC closure: 25 for data sigmas, 50 for MC sigmas
nSigma=3.00 #needed in 2 decimal accuracy (x.yz)

for nState in 4;do

JobID=FirstNchRun1S #Please define nSigma and fracL yourself in the JobID, if needed
#JobID=JOBID

rapBinMin=1
rapBinMax=1
ptBinMin=$setptbin
ptBinMax=$setptbin
cpmBinMin=$setcpmbin
cpmBinMax=$setcpmbin

FidCuts=11

nEff=1060				#1101 MCtruthFineEta, 1080 MCTnPparam      #1030=soft-1060=tight-1070=mixed-111=soft-112=tight
UseMCeff=false

nDileptonEff=1
UseMCDileptoneff=true

nRhoFactor=1
#nRhoFactor=325

useAmapApproach=false
nAmap=1                    #frame/state/sigma/ID ( ID= 2 digits )
nDenominatorAmap=1		    	#the number here corresponds to the same notation as nEff

#nSample=50000
nSample=50000

nFits=50
nSkipGen=NSkipGen

DataID=_FirstFullNch1S

MPValgo=3 		#1...mean,2...gauss,3...gauss-loop with chi2<2

########################################
#useCentralFracL=0

datadir=${datadir_Start}/SetOfCuts${FidCuts}${DataID}/Psi$[nState-3]S/tmpFiles

TreeID=Psi$[nState-3]S

cd ${homedir}


polScenSig=3
polScenBkg=3
frameSig=1
frameBkg=1
ConstEvents=15000
UseConstEv=true
nGenerations=${nFits}
gen=false
rec=false
fit=true
plot=false
NewAccCalc=false

ScenDir=Default_ScenDir
mkdir -p ${storagedir}/${JobID}

cp ${basedir}/macros/polFit/polGenRecFitPlot.cc ${storagedir}/${JobID}/polGenRecFitPlot.cc
cp ${basedir}/macros/polFit/polRapPtPlot.cc ${storagedir}/${JobID}/polRapPtPlot.cc
cp ${basedir}/macros/polFit/PlotFinalResults.cc ${storagedir}/${JobID}/PlotFinalResults.cc
cp ${basedir}/macros/polFit/Makefile ${storagedir}/${JobID}/Makefile
cp ${basedir}/macros/polFit/polGen.C ${storagedir}/${JobID}/polGen.C
cp ${basedir}/macros/polFit/polRec.C ${storagedir}/${JobID}/polRec.C
cp ${basedir}/macros/polFit/polFit.C ${storagedir}/${JobID}/polFit.C
cp ${basedir}/macros/polFit/polPlot.C ${storagedir}/${JobID}/polPlot.C

cp ../../interface/rootIncludes.inc ${storagedir}/${JobID}/rootIncludes.inc
cp ../../interface/commonVar_Psi$[nState-3]S.h ${storagedir}/${JobID}/commonVar.h
cp ../../interface/ToyMC_Psi$[nState-3]S.h ${storagedir}/${JobID}/ToyMC.h
cp ../../interface/effsAndCuts_Psi$[nState-3]S.h ${storagedir}/${JobID}/effsAndCuts.h

cd ${storagedir}/${JobID}
cp ${basedir}/macros/polFit/runDataFits.sh .

rap_=${rapBinMin}
while [ $rap_ -le ${rapBinMax} ]
do
pT_=${ptBinMin}
while [ $pT_ -le ${ptBinMax} ]
do
cpm_=${cpmBinMin}
while [ $cpm_ -le ${cpmBinMax} ]
do

#touch polGenRecFitPlot.cc
#make

nGen_=${nSkipGen}
nGen_=$[nGen_+1]
nMaxGen=$[nGenerations+nSkipGen]
while [ $nGen_ -le $nMaxGen ]
do

if [ $useBatch -eq 0 ]
then

resultfilename=resultsMerged_Psi$[nState-3]S_rap${rap_}_pT${pT_}_cpm${cpm_}.root
nActualGen=$[nGen_-nSkipGen]
if [ $nSkipGen -ge 0 ]
then
if [ $nActualGen -eq 1 ]
then
cp results_Psi$[nState-3]S_rap${rap_}_pT${pT_}_cpm${cpm_}.root ${resultfilename}
fi
fi

fi

cp ${storagedir}/${JobID}/polGenRecFitPlot ${storagedir}/${JobID}/polGenRecFitPlot_Psi$[nState-3]S_rap${rap_}_pt${pT_}_cpm${cpm_}_${nGen_}
./polGenRecFitPlot_Psi$[nState-3]S_rap${rap_}_pt${pT_}_cpm${cpm_}_${nGen_} ${nGen_}ThisGen ${JobID}=JobID ${storagedir}=storagedir ${basedir}=basedir ${nGenerations}nGenerations ${polScenSig}polScenSig ${frameSig}frameSig ${polScenBkg}polScenBkg ${frameBkg}frameBkg ${rap_}rapBinMin ${rap_}rapBinMax ${pT_}ptBinMin ${pT_}ptBinMax ${cpm_}cpmBinMin ${cpm_}cpmBinMax ${nEff}nEff ${nDileptonEff}nDiEff ${FidCuts}FidCuts ${nSample}nSample ${ConstEvents}ConstEvents ${nSkipGen}nSkipGen UseConstEv=${UseConstEv} gen=${gen} rec=${rec} fit=${fit} plot=${plot} ${TreeID}=TreeID ${datadir}=realdatadir UseMCeff=${UseMCeff} UseMCDileptoneff=${UseMCDileptoneff} ${nRhoFactor}nRhoFactor ${MPValgo}MPValgo NewAccCalc=${NewAccCalc} useAmapApproach=${useAmapApproach} ${nAmap}nAmap ${nDenominatorAmap}nDenominatorAmap useBatch=${useBatch}
rm polGenRecFitPlot_Psi$[nState-3]S_rap${rap_}_pt${pT_}_cpm${cpm_}_${nGen_}


if [ $useBatch -eq 0 ]
then

mv results_Psi$[nState-3]S_rap${rap_}_pT${pT_}_cpm${cpm_}.root results_Fit${nGen_}_Psi$[nState-3]S_rap${rap_}_pT${pT_}_cpm${cpm_}.root

if [ $nGen_ -eq 1 ]
then
cp results_Fit${nGen_}_Psi$[nState-3]S_rap${rap_}_pT${pT_}_cpm${cpm_}.root ${resultfilename}
fi

if [ $nGen_ -ge 2 ]
then
mv ${resultfilename} BUFFER_${resultfilename}
hadd -f ${resultfilename} BUFFER_${resultfilename} results_Fit${nGen_}_Psi$[nState-3]S_rap${rap_}_pT${pT_}_cpm${cpm_}.root
rm BUFFER_${resultfilename}
fi

#cp ${resultfilename} results_MergedUpToFit${nGen_}_$[nState-3]SUps_rap${rap_}_pT${pT_}.root

fi

nGen_=$[nGen_+1]
done

if [ $useBatch -eq 0 ]
then

mv ${resultfilename} results_Psi$[nState-3]S_rap${rap_}_pT${pT_}_cpm${cpm_}.root

cp ${storagedir}/${JobID}/polGenRecFitPlot ${storagedir}/${JobID}/polGenRecFitPlot_Psi$[nState-3]S_rap${rap_}_pt${pT_}_cpm${cpm_}
./polGenRecFitPlot_Psi$[nState-3]S_rap${rap_}_pt${pT_}_cpm${cpm_} ${nGen_}ThisGen ${JobID}=JobID ${storagedir}=storagedir ${basedir}=basedir ${nGenerations}nGenerations ${polScenSig}polScenSig ${frameSig}frameSig ${polScenBkg}polScenBkg ${frameBkg}frameBkg ${rap_}rapBinMin ${rap_}rapBinMax ${pT_}ptBinMin ${pT_}ptBinMax ${cpm_}cpmBinMin ${cpm_}cpmBinMax ${nEff}nEff ${nDileptonEff}nDiEff ${FidCuts}FidCuts ${nSample}nSample ${ConstEvents}ConstEvents ${nSkipGen}nSkipGen UseConstEv=${UseConstEv} gen=false rec=false fit=false plot=true ${TreeID}=TreeID ${datadir}=realdatadir UseMCeff=${UseMCeff} UseMCDileptoneff=${UseMCDileptoneff} ${nRhoFactor}nRhoFactor ${MPValgo}MPValgo scalePlots=true NewAccCalc=${NewAccCalc} useAmapApproach=${useAmapApproach} ${nAmap}nAmap ${nDenominatorAmap}nDenominatorAmap
rm polGenRecFitPlot_Psi$[nState-3]S_rap${rap_}_pt${pT_}_cpm${cpm_}

fi

cpm_=$[cpm_+1]
done
pT_=$[pT_+1]
done
rap_=$[rap_+1]
done

done

rm *.so
rm *.d

mkdir ../tmp
