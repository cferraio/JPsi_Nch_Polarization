#!/bin/sh
source /cvmfs/sft.cern.ch/lcg/app/releases/ROOT/5.34.05/x86_64-slc5-gcc43-opt/root/bin/thisroot.sh                                     
source /cvmfs/sft.cern.ch/lcg/app/releases/ROOT/5.34.05/x86_64-slc5-gcc43-opt/root/bin/setxrd.sh /cvmfs/sft.cern.ch/lcg/external/xrootd/3.2.4/x86_64-slc5-gcc46-opt/


homedir=$PWD
cd ${homedir}
cd ..
cd ..
basedir=$PWD
cd macros/polFit
#storagedir=`more storagedir`/Data #please define the directory storagedir in the file macros/polFit/storagedir
storagedir=`more storagedir`/Data
datadir_Start=${basedir}/macros/DataFiles

########## INPUTS ##########

#Batch submission system: 0/1
useBatch=0

#fracL=50 #in percent #MC closure: 25 for data sigmas, 50 for MC sigmas
#nSigma=3.00 #needed in 2 decimal accuracy (x.yz)

for nState in 4;do

StatVarTotBGfraction=0     #apply statistical fluctuations on f_background
StatVarTotBGmodel=0        #apply statistical fluctuations on Bg model
StatVarRho=0               #apply statistical fluctuations on rho factor

#JobID=Psi$[nState-3]S_${nSigma}Sigma_11Dec2012
#JobID=Psi$[nState-3]S_${nSigma}Sigma_11Dec2012_noRhoFactor
#JobID=Psi$[nState-3]S_ctauScen0_FracLSB-1_noRhoFactor
#JobID=Psi$[nState-3]S_ctauScen0_FracLSB-1_rho_20Feb2013
#JobID=Psi$[nState-3]S_ctauScen0_FracLSB-1_rho_21Feb2013
#JobID=Psi$[nState-3]S_ctauScen0_FracLSB-1_rho_25Feb2013
#JobID=Psi$[nState-3]S_ctauScen0_FracLSB-1_rho_25Feb2013_massRange
#JobID=Psi$[nState-3]S_ctauScen0_FracLSB-1_rho_25Feb2013_massRange_Bin20_2_2
#JobID=Psi$[nState-3]S_ctauScen0_FracLSB-1_rho_25Feb2013_massRange_Bin5_8_8
#JobID=Psi$[nState-3]S_ctauScen0_FracLSB-1_rho_26Feb2013_BgNoRebin
#JobID=Psi$[nState-3]S_ctauScen0_FracLSB-1_16Mar2013
#JobID=Psi$[nState-3]S_ctauScen0_FracLSB-1_16Mar2013_scaleFracBg
#JobID=Psi$[nState-3]S_ctauScen0_FracLSB-1_19Mar2013_${StatVarTotBGfraction}FracBg_${StatVarTotBGmodel}BgModel_${StatVarRho}Rho
JobID=FirstTry

rapBinMin=1
rapBinMax=2
ptBinMin=1
ptBinMax=12

FidCuts=11

nEff=1103				#1101 MCtruthFineEta, 1080 MCTnPparam      #1030=soft-1060=tight-1070=mixed-111=soft-112=tight
UseMCeff=false

nDileptonEff=1
UseMCDileptoneff=true

#nRhoFactor=1
#nRhoFactor=325 ## old 
#nRhoFactor=326 ## newest
nRhoFactor=329 ## from ilse 17 sep 2015

useAmapApproach=false
nAmap=1                    #frame/state/sigma/ID ( ID= 2 digits )
nDenominatorAmap=1		     #the number here corresponds to the same notation as nEff

nSample=20000

nFits=1
nSkipGen=0

#DataID=_FrameworkTest_5Dec2012
#DataID=_ctauScen0_FracLSB-1
#DataID=_ctauScen0_FracLSB-1_21Feb2013 ##with correct PR region definition
#DataID=_ctauScen0_FracLSB-1_25Feb2013 ##with correct PR region definition
#DataID=_ctauScen0_FracLSB-1_25Feb2013_Bin20_2_2 ##with correct PR region definition
#DataID=_ctauScen0_FracLSB-1_25Feb2013_Bin5_8_8 ##with correct PR region definition
#DataID=_ctauScen0_FracLSB-1_26Feb2013_BgNoRebin ##with correct PR region definition
#DataID=_ctauScen0_FracLSB-1_newMLfit_4Mar2013
DataID=_ctauScen0_FracLSB-1_newMLfit_31Jul2015
#DataID=_ctauScen0_FracLSB-1_newMLfit_4Mar2013_scaleFracBg

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

touch polGenRecFitPlot.cc
make

rap_=${rapBinMin}
while [ $rap_ -le ${rapBinMax} ]
do
pT_=${ptBinMin}
while [ $pT_ -le ${ptBinMax} ]
do

nGen_=${nSkipGen}
nGen_=$[nGen_+1]
nMaxGen=$[nGenerations+nSkipGen]
while [ $nGen_ -le $nMaxGen ]
do

if [ $useBatch -eq 0 ]
then

resultfilename=resultsMerged_Psi$[nState-3]S_rap${rap_}_pT${pT_}.root
nActualGen=$[nGen_-nSkipGen]
if [ $nSkipGen -ge 0 ]
then
if [ $nActualGen -eq 1 ]
then
cp results_Psi$[nState-3]S_rap${rap_}_pT${pT_}.root ${resultfilename}
fi
fi

fi

cp ${storagedir}/${JobID}/polGenRecFitPlot ${storagedir}/${JobID}/polGenRecFitPlot_Psi$[nState-3]S_rap${rap_}_pt${pT_}_Gen${nGen_}
./polGenRecFitPlot_Psi$[nState-3]S_rap${rap_}_pt${pT_}_Gen${nGen_} ${nGen_}ThisGen ${JobID}=JobID ${storagedir}=storagedir ${basedir}=basedir ${nGenerations}nGenerations ${polScenSig}polScenSig ${frameSig}frameSig ${polScenBkg}polScenBkg ${frameBkg}frameBkg ${rap_}rapBinMin ${rap_}rapBinMax ${pT_}ptBinMin ${pT_}ptBinMax ${nEff}nEff ${nDileptonEff}nDiEff ${FidCuts}FidCuts ${nSample}nSample ${ConstEvents}ConstEvents ${nSkipGen}nSkipGen UseConstEv=${UseConstEv} gen=${gen} rec=${rec} fit=${fit} plot=${plot} ${TreeID}=TreeID ${datadir}=realdatadir UseMCeff=${UseMCeff} UseMCDileptoneff=${UseMCDileptoneff} ${nRhoFactor}nRhoFactor ${MPValgo}MPValgo NewAccCalc=${NewAccCalc} useAmapApproach=${useAmapApproach} ${nAmap}nAmap ${nDenominatorAmap}nDenominatorAmap useBatch=${useBatch} StatVarTotBGfraction=${StatVarTotBGfraction} StatVarTotBGmodel=${StatVarTotBGmodel} StatVarRho=${StatVarRho}
rm polGenRecFitPlot_Psi$[nState-3]S_rap${rap_}_pt${pT_}_Gen${nGen_}


if [ $useBatch -eq 0 ]
then

mv results_Psi$[nState-3]S_rap${rap_}_pT${pT_}.root results_Fit${nGen_}_Psi$[nState-3]S_rap${rap_}_pT${pT_}.root

if [ $nGen_ -eq 1 ]
then
cp results_Fit${nGen_}_Psi$[nState-3]S_rap${rap_}_pT${pT_}.root ${resultfilename}
fi

if [ $nGen_ -ge 2 ]
then
mv ${resultfilename} BUFFER_${resultfilename}
hadd -f ${resultfilename} BUFFER_${resultfilename} results_Fit${nGen_}_Psi$[nState-3]S_rap${rap_}_pT${pT_}.root
rm BUFFER_${resultfilename}
fi

#cp ${resultfilename} results_MergedUpToFit${nGen_}_$[nState-3]SUps_rap${rap_}_pT${pT_}.root

fi

nGen_=$[nGen_+1]
done

if [ $useBatch -eq 0 ]
then

mv ${resultfilename} results_Psi$[nState-3]S_rap${rap_}_pT${pT_}.root

cp ${storagedir}/${JobID}/polGenRecFitPlot ${storagedir}/${JobID}/polGenRecFitPlot_Psi$[nState-3]S_rap${rap_}_pt${pT_}
./polGenRecFitPlot_Psi$[nState-3]S_rap${rap_}_pt${pT_} ${nGen_}ThisGen ${JobID}=JobID ${storagedir}=storagedir ${basedir}=basedir ${nGenerations}nGenerations ${polScenSig}polScenSig ${frameSig}frameSig ${polScenBkg}polScenBkg ${frameBkg}frameBkg ${rap_}rapBinMin ${rap_}rapBinMax ${pT_}ptBinMin ${pT_}ptBinMax ${nEff}nEff ${nDileptonEff}nDiEff ${FidCuts}FidCuts ${nSample}nSample ${ConstEvents}ConstEvents ${nSkipGen}nSkipGen UseConstEv=${UseConstEv} gen=false rec=false fit=false plot=true ${TreeID}=TreeID ${datadir}=realdatadir UseMCeff=${UseMCeff} UseMCDileptoneff=${UseMCDileptoneff} ${nRhoFactor}nRhoFactor ${MPValgo}MPValgo scalePlots=true NewAccCalc=${NewAccCalc} useAmapApproach=${useAmapApproach} ${nAmap}nAmap ${nDenominatorAmap}nDenominatorAmap ${nState}nState
rm polGenRecFitPlot_Psi$[nState-3]S_rap${rap_}_pt${pT_}

fi

pT_=$[pT_+1]
done
rap_=$[rap_+1]
done

done

rm *.so
rm *.d

mkdir ../tmp
