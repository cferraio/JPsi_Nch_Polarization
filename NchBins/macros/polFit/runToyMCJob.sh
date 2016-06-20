#!/bin/sh
export VO_CMS_SW_DIR=/sharesoft/cmssw #comment out for non-condor
. $VO_CMS_SW_DIR/cmsset_default.sh #comment out for non-condor
cd /home/ferraioc/PolNew/CMSSW_5_3_20/src/JPsi_Nch_Polarization/NchBins/macros/polFit #comment out for non-condor
eval `scramv1 runtime -sh` #comment out for non-condor

source /cvmfs/sft.cern.ch/lcg/app/releases/ROOT/5.34.05/x86_64-slc5-gcc43-opt/root/bin/thisroot.sh                                     
source /cvmfs/sft.cern.ch/lcg/app/releases/ROOT/5.34.05/x86_64-slc5-gcc43-opt/root/bin/setxrd.sh /cvmfs/sft.cern.ch/lcg/external/xrootd/3.2.4/x86_64-slc5-gcc46-opt/

#### following will be changed by "launchFitToy.sh" ####
# HOMEDIR
# NState 
# RapBinMin
# RapBinMax
# PtBinMin
# PtBinMax
# NSkipGen
# PolScenSig
# PolScenBkg
# FrameSig
# FrameBkg
# USEConstEv
######################################### 

########## INPUTS ##########

storagedir=$3
homedir=$4
setptbin=$5
setcpmbin=$6
setgen=$7
sigpol=$8

echo "storage dir: $storagedir"
echo "home dir: $homedir"
echo "ptbin: $setptbin"
echo "cpmbin: $setcpmbin"
echo "sigpol: $sigpol"


nState=4


#JobID=ToyMC_Psi$[nState-3]S_Eff_pTshift_minus_8Jan2013
#JobID=ToyMC_Psi$[nState-3]S_13Dec2012_100K
JobID=FrameworkII_19May2016
#JobID=ToyMC_Psi$[nState-3]S_13Dec2012

nGenerations=1

rapBinMin=1
rapBinMax=1
ptBinMin=$setptbin
ptBinMax=$setptbin
cpmBinMin=$setcpmbin
cpmBinMax=$setcpmbin

polScenSig=$sigpol
polScenBkg=3
frameSig=3 #1=CS, 2=HX, 3=PX
frameBkg=1

nEff=1060 #MC-truth
#nEff=1062
UseMCeff=false
nDileptonEff=1
UseMCDileptoneff=false
nRhoFactor=1

useAmapApproach=false       #if false, the next two lines do not matter
nAmap=32104                 #frame/state/sigma/ID ( ID= 2 digits )
nDenominatorAmap=113 		    #the number here corresponds to the same notation as nEff
 
FidCuts=11

nSample=10000
nSkipGen=$setgen

#GENERATION SETTINGS
ConstEvents=200000
UseConstEv=true #if false, the number of events is taken from ToyMC.h
UseDifferingEff=false #if false, the next five lines do not matter
nEffRec=1060 #1101
UseMCReceff=false
nDileptonEffRec=1
UseMCDileptonReceff=true
nRecRhoFactor=1

gen=true
rec=true
fit=true
plot=false
deletePseudoData=true

MPValgo=3 		#1...mean,2...gauss,3...gauss-loop with chi2<2
nSigma=1
NewAccCalc=false


########################################

#homedir=homedir
cd ${homedir}
cd ..
cd ..
basedir=$PWD
cd macros/polFit
#storagedir=`more storagedir`/ToyMC #please define the directory storagedir in the file macros/polFit/storagedir
storagedir=${storagedir}/ToyMC

ScenDir=Sig_frame${frameSig}scen${polScenSig}_Bkg_frame${frameBkg}scen${polScenBkg}

mkdir -p ${storagedir}/${JobID}/${ScenDir}

cp ${basedir}/macros/polFit/polGenRecFitPlot.cc ${storagedir}/${JobID}/${ScenDir}/polGenRecFitPlot.cc
cp ${basedir}/macros/polFit/polRapPtPlot.cc ${storagedir}/${JobID}/${ScenDir}/polRapPtPlot.cc
cp ${basedir}/macros/polFit/PlotFinalResults.cc ${storagedir}/${JobID}/${ScenDir}/PlotFinalResults.cc
cp ${basedir}/macros/polFit/Makefile ${storagedir}/${JobID}/${ScenDir}/Makefile
cp ${basedir}/macros/polFit/polGen.C ${storagedir}/${JobID}/${ScenDir}/polGen.C
cp ${basedir}/macros/polFit/polRec.C ${storagedir}/${JobID}/${ScenDir}/polRec.C
cp ${basedir}/macros/polFit/polFit.C ${storagedir}/${JobID}/${ScenDir}/polFit.C
cp ${basedir}/macros/polFit/polPlot.C ${storagedir}/${JobID}/${ScenDir}/polPlot.C

cp ../../interface/rootIncludes.inc ${storagedir}/${JobID}/${ScenDir}/rootIncludes.inc
cp ../../interface/commonVar_Psi$[nState-3]S.h ${storagedir}/${JobID}/${ScenDir}/commonVar.h
cp ../../interface/ToyMC_Psi$[nState-3]S.h ${storagedir}/${JobID}/${ScenDir}/ToyMC.h
cp ../../interface/effsAndCuts_Psi$[nState-3]S.h ${storagedir}/${JobID}/${ScenDir}/effsAndCuts.h

cd ${storagedir}/${JobID}/${ScenDir}
cp ${basedir}/macros/polFit/runToyMC.sh .

touch polGenRecFitPlot.cc
make

rap_=${rapBinMin}
while [ $rap_ -le ${rapBinMax} ]
do
pT_=${ptBinMin}
while [ $pT_ -le ${ptBinMax} ]
do
cpm_=${cpmBinMin}
while [ $cpm_ -le ${cpmBinMax} ]
do


nGen_=${nSkipGen}
nGen_=$[nGen_+1]
nMaxGen=$[nGenerations+nSkipGen]
while [ $nGen_ -le $nMaxGen ]
do

plot=${plot}

if [ ${nGen_} -eq 1 ]
then
plot=true
fi


cp polGenRecFitPlot polGenRecFitPlot_rap${rap_}_pt${pT_}_cpm${cpm_}_Gen${nGen_}
./polGenRecFitPlot_rap${rap_}_pt${pT_}_cpm${cpm_}_Gen${nGen_} ${nGen_}ThisGen ${JobID}=JobID ${storagedir}=storagedir ${basedir}=basedir ${nGenerations}nGenerations ${polScenSig}polScenSig ${frameSig}frameSig ${polScenBkg}polScenBkg ${frameBkg}frameBkg ${rap_}rapBinMin ${rap_}rapBinMax ${pT_}ptBinMin ${pT_}ptBinMax ${cpm_}cpmBinMin ${cpm_}cpmBinMax ${nEff}nEff ${nDileptonEff}nDiEff ${nEffRec}nRecEff ${nDileptonEffRec}nRecDiEff ${FidCuts}FidCuts ${nSample}nSample ${ConstEvents}ConstEvents ${nSkipGen}nSkipGen UseConstEv=${UseConstEv} gen=${gen} rec=${rec} fit=${fit} plot=${plot} UseDifferingEff=${UseDifferingEff} UseMCeff=${UseMCeff} UseMCReceff=${UseMCReceff} UseMCDileptoneff=${UseMCDileptoneff} UseMCDileptonReceff=${UseMCDileptonReceff}  ${nRhoFactor}nRhoFactor ${nRecRhoFactor}nRecRhoFactor ${MPValgo}MPValgo ${nSigma}nSigma ${nState}nState NewAccCalc=${NewAccCalc} deletePseudoData=${deletePseudoData} useAmapApproach=${useAmapApproach} ${nAmap}nAmap ${nDenominatorAmap}nDenominatorAmap
rm polGenRecFitPlot_rap${rap_}_pt${pT_}_cpm${cpm_}_Gen${nGen_}

nGen_=$[nGen_+1]
done
cpm_=$[cpm_+1]
done
pT_=$[pT_+1]
done
rap_=$[rap_+1]
done


#rm polGen.C 
#rm polRec.C 
#rm polFit.C 
#rm polPlot.C 

