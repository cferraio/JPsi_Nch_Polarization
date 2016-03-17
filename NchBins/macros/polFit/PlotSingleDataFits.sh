#!/bin/sh

homedir=$PWD
cd ..
cd ..
basedir=$PWD
cd macros/polFit
#storagedir=`more storagedir`/Data #please define the directory storagedir in the file macros/polFit/storagedir
storagedir=$basedir/Psi/Data
datadir_Start=${basedir}/macros/DataFiles

########## INPUTS ##########

#fracL=0 #in percent
nSigma=3.00 #needed in 2 decimal accuracy (x.yz)

for nState in 4;do

#JobID=Psi$[nState-3]S_${nSigma}Sigma_11Dec2012
JobID=Psi$[nState-3]S_${nSigma}Sigma_11Dec2012_noRhoFactor

rapBinMin=1
rapBinMax=2
ptBinMin=1
ptBinMax=4

#NOTE: take care of nFits!
nFits=10

########################################

FidCuts=11
nEff=1060
UseMCeff=false
nDileptonEff=1
UseMCDileptoneff=true
nRhoFactor=325
nSample=60000
nSkipGen=0
DataID=_FrameworkTest_5Dec2012
MPValgo=3 		#1...mean,2...gauss,3...gauss-loop with chi2<2

additionalName=MPV${MPValgo}

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
fit=false
plot=true
NewAccCalc=false

touch polGenRecFitPlot.cc
make

ScenDir=Default_ScenDir
mkdir ${storagedir}
mkdir ${storagedir}/${JobID}

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
cp ${basedir}/macros/polFit/PlotSingleDataFits.sh .

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

resultfilename=resultsMerged_Psi$[nState-3]S_rap${rap_}_pT${pT_}.root
nActualGen=$[nGen_-nSkipGen]

nGen_=$[nGen_+1]
done


cp ${storagedir}/${JobID}/polGenRecFitPlot ${storagedir}/${JobID}/polGenRecFitPlot_rap${rap_}_pt${pT_}
./polGenRecFitPlot_rap${rap_}_pt${pT_} ${nGen_}ThisGen ${JobID}=JobID ${storagedir}=storagedir ${basedir}=basedir ${nGenerations}nGenerations ${polScenSig}polScenSig ${frameSig}frameSig ${polScenBkg}polScenBkg ${frameBkg}frameBkg ${rap_}rapBinMin ${rap_}rapBinMax ${pT_}ptBinMin ${pT_}ptBinMax ${nEff}nEff ${nDileptonEff}nDiEff ${FidCuts}FidCuts ${nSample}nSample ${nState}nState ${ConstEvents}ConstEvents ${nSkipGen}nSkipGen UseConstEv=${UseConstEv} gen=false rec=false fit=false plot=${plot} ${TreeID}=TreeID ${datadir}=realdatadir UseMCeff=${UseMCeff} UseMCDileptoneff=${UseMCDileptoneff} ${nRhoFactor}nRhoFactor ${MPValgo}MPValgo scalePlots=true NewAccCalc=${NewAccCalc}
rm polGenRecFitPlot_rap${rap_}_pt${pT_}


cp ${basedir}/latex/IndividualFitResults.tex IndividualFitResults.tex
pdflatex "\newcommand\TreeBinID{${TreeID}_rap${rap_}_pT${pT_}}\input{IndividualFitResults.tex}"
mv IndividualFitResults.pdf ${basedir}/macros/polFit/FiguresData/${JobID}/${TreeID}/IndividualFitResults_rap${rap_}pt${pT_}_${additionalName}.pdf

pT_=$[pT_+1]
done
rap_=$[rap_+1]
done


done

rm *.so
rm *.d
#rm polGenRecFitPlot

mkdir ../tmp
#mv *.C ../tmp

