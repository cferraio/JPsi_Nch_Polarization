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
NSigma=3.00 #needed in 2 decimal accuracy (x.yz)

for nState in 4;do

cp ../../interface/ToyMC_Psi$[nState-3]S.h       ToyMC.h
cp ../../interface/effsAndCuts_Psi$[nState-3]S.h effsAndCuts.h
cp ../../interface/rootIncludes.inc              rootIncludes.inc
cp ../../interface/commonVar_Psi$[nState-3]S.h   commonVar.h

#JobID=CentralsApr8_CentralsFromAlteredPPDApr5_1SigmaStatError
#JobID=CentralsApr16_CentralsFromAlteredPPDApr5_1SigmaStatError
#JobID=CentralsApr22_CentralsFromAlteredPPDApr5_1SigmaStatError_fLSB75_25
#JobID=CentralsApr23_CentralsFromAlteredPPDApr23_1SigmaStatError_fLSB75_25_noRhoPt35_Pt12_52
#JobID=CentralsMay3_CentralsFromAlteredPPDMay3_1SigmaStatError_correctfLSB_noRhoPt35
#JobID=CentralsMay9_CentralsFromAlteredPPDMay9_1SigmaStatError_ctauScen5
#JobID=CentralsMay9_CentralsFromAlteredPPDMay9_1SigmaStatError_ctauScen5_CL1sigma
JobID=CentralsMay28_CentralsFromAlteredPPDMay26_1SigmaStatError_ctauScen5_RhoAll

additionalName=_Psi$[nState-3]S

PlotMatt=0
PlotMattForICHEP=0
PlotCompare=0

PlotAsymm=0
PlotFinalData=1
PlotSystematics=0
PlotLegend=0
PlotBrazilian=1
FitGraph=0
DrawLatexStuff=1
DrawPreliminary=1
MultiPanelPlots=1
MPCentralsWithTotalSystID=CentralsMay28_CentralsFromAlteredPPDMay26_1SigmaStatError_ctauScen5_RhoAll
#CentralsMay9_CentralsFromAlteredPPDMay9_1SigmaStatError_ctauScen5
#MPCentralsWithTotalSystID=CentralsMay3_CentralsFromAlteredPPDMay3_1SigmaStatError_correctfLSB_noRhoPt35
#MPCentralsWithTotalSystID=CentralsApr23_CentralsFromAlteredPPDApr23_1SigmaStatError_fLSB75_25_noRhoPt35
#CentralsApr22_CentralsFromAlteredPPDApr5_1SigmaStatError_fLSB75_25 
#CentralsApr16_CentralsFromAlteredPPDApr5_1SigmaStatError
PlotAlteredPPDResults=1
PlotCL1sigma=0

#UncommentIFplotSystOverview
#PlotFinalData=0
#PlotSystematics=1
#PlotLegend=1
#PlotBrazilian=0
#FitGraph=0
#DrawLatexStuff=1
#MultiPanelPlots=0

#UncommentIFplotCDFComparison
#PlotMatt=1
#PlotCompare=1
#PlotAsymm=0
#PlotFinalData=1
#PlotSystematics=0
#PlotLegend=0
#PlotBrazilian=0
#FitGraph=0
#DrawLatexStuff=1
#MultiPanelPlots=1

DefaultID=ctauScen5_FracLSB-1_26May2013_1FracBg_1BgModel_1RhoAll_AlteredPPD_BKGlinPLUSRestSquaredGauss_5nRand_May28
#DefaultID=ctauScen5_FracLSB-1_7May2013_1FracBg_1BgModel_AlteredPPD_BKGlinPLUSRestSquaredGauss_5nRand_May9
#DefaultID=ctauScen0_FracLSB-1_30Apr2013_1FracBg_1BgModel_noRhoPt35_correctfLSB_AlteredPPD_BKGAndRholinPLUSRestSquaredGauss_5nRand_May3
#DefaultID=ctauScen0_FracLSB-1_23Apr2013_1FracBg_1BgModel_noRhoPt35_AlteredPPD_BKGAndRholinPLUSRestSquaredGauss_5nRand_Apr23
#DefaultID=ctauScen0_FracLSB-1_19Mar2013_1FracBg_1BgModel_1Rho_AlteredPPD_BKGlinPLUSRestSquaredGauss_5nRand_Apr22
#DefaultID=ctauScen0_FracLSB-1_19Mar2013_1FracBg_1BgModel_1Rho_AlteredPPD_BKGlinPLUSRestSquaredGauss_5nRand_Apr16
#DefaultID=ctauScen0_FracLSB-1_19Mar2013_1FracBg_1BgModel_1Rho_AlteredPPD_BKGlinPLUSRestSquaredGauss_5nRand_Apr5
CompareID1=Psi$[nState-3]S_${NSigma}Sigma_11Dec2012_noRhoFactor
CompareID2=Data_TheGreatRun_10B_May11_NewestCentrals_AlteredPPD_May17_BKGlinPLUSRestSquaredGauss_10nRand
CompareID3=BG0_Mar19_HighCtauSigCheck3p25
CompareID4=BG0_Mar19_HighCtauSigCheck3p5
nComp=0

LegendEntryDefID=with_RhoFactor
LegendEntryCompID1=no_RhoFactor
LegendEntryCompID2=xxx
LegendEntryCompID3=xxx
LegendEntryCompID4=xxx

nSystematics=0
#nSystematics=7

SystID1Base=TheGreatRun_BKGmodel
SystID1Specify=BestSyst_28p_SQRT12
SystID1Title=BKGmodel

SystID2Base=TheGreatRun_Rho
SystID2Specify=ConstSyst
SystID2Title=RhoFactor

SystID3Base=TheGreatRun_Param
SystID3Specify=BestSyst
SystID3Title=Parametrization

SystID4Base=TheGreatRun_TnP
SystID4Specify=BestSyst
SystID4Title=TnP_model

SystID5Base=TheGreatRun_FrameworkII
SystID5Specify=BestSyst_Bkg
SystID5Title=FrameworkIII

SystID6Base=TheGreatRun_FrameworkII
SystID6Specify=BestSyst_Sig_NoUnpol
SystID6Title=FrameworkII

SystID7Base=TheGreatRun_FrameworkI
SystID7Specify=BestSyst
SystID7Title=FrameworkI

if [ $nState -eq 4 ]
then
ptBinMin=3
ptBinMax=12
fi

if [ $nState -eq 5 ]
then
ptBinMin=2
ptBinMax=5
fi

########################################

cd ${homedir}

touch PlotFinalResults.cc
make

JobIDDir=FinalResults/${JobID}
mkdir -p  FinalResults/${JobID}/Psi$[nState-3]S
mkdir -p ${JobIDDir}

cp PlotFinalResults PlotFinalResults_Psi$[nState-3]S
./PlotFinalResults_Psi$[nState-3]S PlotAlteredPPDResults=${PlotAlteredPPDResults} MultiPanelPlots=${MultiPanelPlots} DrawLatexStuff=${DrawLatexStuff} ${MPCentralsWithTotalSystID}=MPCentralsWithTotalSystID ${DefaultID}=DefaultID ${CompareID1}=CompareID1 ${CompareID2}=CompareID2 ${CompareID3}=CompareID3 ${CompareID4}=CompareID4 ${JobID}=JobID ${SystID1Base}=SystID1Base ${SystID1Specify}=SystID1Specify ${SystID1Title}=SystID1Title ${SystID2Base}=SystID2Base ${SystID2Specify}=SystID2Specify ${SystID2Title}=SystID2Title ${basedir}=basedir ${storagedir}=storagedir ${ptBinMin}ptBinMin ${ptBinMax}ptBinMax ${nSystematics}nSystematics ${nComp}nComp ${nState}nState ${SystID3Base}=SystID3Base ${SystID3Specify}=SystID3Specify ${SystID3Title}=SystID3Title ${SystID4Base}=SystID4Base ${SystID4Specify}=SystID4Specify ${SystID4Title}=SystID4Title ${SystID5Base}=SystID5Base ${SystID5Specify}=SystID5Specify ${SystID5Title}=SystID5Title ${SystID6Base}=SystID6Base ${SystID6Specify}=SystID6Specify ${SystID6Title}=SystID6Title ${SystID7Base}=SystID7Base ${SystID7Specify}=SystID7Specify ${SystID7Title}=SystID7Title ${SystID8Base}=SystID8Base ${SystID8Specify}=SystID8Specify ${SystID8Title}=SystID8Title PlotMatt=${PlotMatt} PlotAsymm=${PlotAsymm} PlotCompare=${PlotCompare} PlotFinalData=${PlotFinalData} PlotSystematics=${PlotSystematics} PlotLegend=${PlotLegend} PlotBrazilian=${PlotBrazilian} FitGraph=${FitGraph} DrawPreliminary=${DrawPreliminary} PlotMattForICHEP=${PlotMattForICHEP} ${LegendEntryDefID}=LegendEntryDefID ${LegendEntryCompID1}=LegendEntryCompID1 PlotCL1sigma=${PlotCL1sigma}
rm PlotFinalResults_Psi$[nState-3]S
rm PlotFinalResults

cd ${homedir}/FinalResults/${JobID}/Psi$[nState-3]S
if [ ${PlotFinalData} -eq 1 ]
then
cp ${basedir}/latex/FinalDataResults_Psi$[nState-3]S.tex ./FinalDataResults.tex
pdflatex FinalDataResults.tex
mv FinalDataResults.pdf FinalDataResults${additionalName}.pdf
cd Figures
pdflatex FinalNumericalResults.tex
rm *.aux
rm *.log
cd ../
fi

if [ ${PlotSystematics} -eq 1 ]
then
cp ${basedir}/latex/Systematics.tex .
pdflatex Systematics.tex
mv Systematics.pdf Systematics${additionalName}.pdf
fi
rm *.aux
rm *.log
rm *.tex

#rm -r Figures${additionalName}
mkdir Figures${additionalName}
mv Figures/* Figures${additionalName}/

cd ${homedir}
rm ToyMC.h
rm effsAndCuts.h
rm rootIncludes.inc
rm commonVar.h

done

