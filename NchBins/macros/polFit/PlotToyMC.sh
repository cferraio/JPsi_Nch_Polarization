#!/bin/sh
source /cvmfs/sft.cern.ch/lcg/app/releases/ROOT/5.34.05/x86_64-slc5-gcc43-opt/root/bin/thisroot.sh                                     
source /cvmfs/sft.cern.ch/lcg/app/releases/ROOT/5.34.05/x86_64-slc5-gcc43-opt/root/bin/setxrd.sh /cvmfs/sft.cern.ch/lcg/external/xrootd/3.2.4/x86_64-slc5-gcc46-opt/

homedir=$PWD
cd ..
cd ..
basedir=$PWD
cd macros/polFit
#storagedir=`more storagedir`/ToyMC #please define the directory storagedir in the file macros/polFit/storagedir
storagedir=/data/users/ferraioc/Polarization/JPsi/NchBins/ToyMC

########## INPUTS ##########

for nState in 4; do

cp ../../interface/rootIncludes.inc               rootIncludes.inc
cp ../../interface/commonVar_Psi$[nState-3]S.h    commonVar.h
cp ../../interface/ToyMC_Psi$[nState-3]S.h        ToyMC.h
cp ../../interface/effsAndCuts_Psi$[nState-3]S.h  effsAndCuts.h
touch polRapPtPlot.cc
make

for JobID in FrameworkI_19May2016; do
#for JobID in ToyMC_Psi$[nState-3]S_13Dec2012; do

echo ${JobID}


if [ $nState -eq 4 ]
then
ptBinMin=1
ptBinMax=2
cpmBinMin=1
cpmBinMax=10
fi
if [ $nState -eq 5 ]
then
ptBinMin=2
ptBinMax=6
fi

frameSig=1
for polScenSig in 3;do

frameBkg=1
for polScenBkg in 3;do

nGenerations=50

MPValgo=3 		#1...mean,2...gauss,3...gauss-loop with chi2<2
additionalName=MPV${MPValgo}

############################

TreeID=Psi$[nState-3]S

cd ${basedir}/macros/polFit

rapBinMin=1 #don't change
if [ $nState -eq 4 ]
then
rapBinMax=2 #don't change
fi
if [ $nState -eq 5 ]
then
rapBinMax=3 #don't change
fi

ScenDir=Sig_frame${frameSig}scen${polScenSig}_Bkg_frame${frameBkg}scen${polScenBkg}

mkdir ${basedir}/macros/polFit/FiguresToyMC
mkdir ${basedir}/macros/polFit/FiguresToyMC/${JobID}
mkdir ${basedir}/macros/polFit/FiguresToyMC/${JobID}/${ScenDir}

cd ${storagedir}/${JobID}
mkdir ${ScenDir}

cp ${basedir}/macros/polFit/polRapPtPlot .

./polRapPtPlot ${cpmBinMin}cpmBinMin ${cpmBinMax}cpmBinMax ${ptBinMin}ptBinMin ${ptBinMax}ptBinMax ${rapBinMin}rapBinMin ${rapBinMax}rapBinMax ${frameSig}frameSig ${polScenSig}polScen ${MPValgo}MPValgo ${nGenerations}nGenerations ${ScenDir}=dirstruct ${nState}nState

mv ${ScenDir}/TGraphResults_${TreeID}_temp.root ${ScenDir}/TGraphResults_${TreeID}.root 

cp ${basedir}/latex/PullSummaryResults_Psi$[nState-3]S.tex ${ScenDir}/PullSummaryResults_${ScenDir}.tex
cp ${basedir}/latex/ParameterSummaryResults_Psi$[nState-3]S.tex ${ScenDir}/ParameterSummaryResults_${ScenDir}.tex
cp ${basedir}/latex/ToyResults.tex ${ScenDir}/ToyResults_${ScenDir}.tex

pdflatex ToyNumericalResults_${ScenDir}.tex
mv ToyNumericalResults_${ScenDir}.pdf ${basedir}/macros/polFit/FiguresToyMC/${JobID}/${ScenDir}/ToyNumericalResults_${ScenDir}_${additionalName}.pdf
rm *.aux
rm *.log

cd ${ScenDir}
pdflatex PullSummaryResults_${ScenDir}.tex
pdflatex ParameterSummaryResults_${ScenDir}.tex

rap_=${rapBinMin}
while [ $rap_ -le ${rapBinMax} ]
do
pT_=${ptBinMin}
while [ $pT_ -le ${ptBinMax} ]
do
cpm_=${cpmBinMin}
while [ $cpm_ -le ${cpmBinMax} ]
do

pdflatex "\newcommand\rapptcpm{rap${rap_}pt${pT_}cpm${cpm_}}\input{ToyResults_${ScenDir}.tex}"
mv ToyResults_${ScenDir}.pdf ${basedir}/macros/polFit/FiguresToyMC/${JobID}/${ScenDir}/ToyResults_${ScenDir}_rap${rap_}pt${pT_}_cpm${cpm_}_${additionalName}.pdf

cpm_=$[cpm_+1]
done
pT_=$[pT_+1]
done
rap_=$[rap_+1]
done

mv PullSummaryResults_${ScenDir}.pdf ${basedir}/macros/polFit/FiguresToyMC/${JobID}/${ScenDir}/PullSummaryResults_${ScenDir}_${additionalName}.pdf
mv ParameterSummaryResults_${ScenDir}.pdf ${basedir}/macros/polFit/FiguresToyMC/${JobID}/${ScenDir}/ParameterSummaryResults_${ScenDir}_${additionalName}.pdf

rm *.aux
rm *.log
rm *.tex

cd ..
rm polRapPtPlot

done
done
done

cd ${basedir}/macros/polFit
rm polRapPtPlot
done
