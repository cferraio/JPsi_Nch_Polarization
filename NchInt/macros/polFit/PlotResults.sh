#!/bin/sh

homedir=$PWD
cd ..
cd ..
basedir=$PWD
cd macros/polFit
#storagedir=`more storagedir`/Data #please define the directory storagedir in the file macros/polFit/storagedir
#storagedir=`more storagedir`/ToyMC #please define the directory storagedir in the file macros/polFit/storagedir
storagedir=$basedir/Psi/Data

########## INPUTS ##########
NSigma=3.00 #needed in 2 decimal accuracy (x.yz)

for nState in 4 5;do

JobID=Psi$[nState-3]S_${NSigma}Sigma_11Dec2012
#JobID=Psi$[nState-3]S_${NSigma}Sigma_11Dec2012_noRhoFactor
#JobID=Psi$[nState-3]S_${NSigma}Sigma_11Dec2012_compareRhoFactor
additionalName=_Psi$[nState-3]S

PlotMatt=0
PlotCompare=0

PlotAsymm=0
PlotFinalData=1
PlotSystematics=0
PlotLegend=0
PlotBG0plots=0
DeltaTildeplots=0
SBmSigPlots=0
SteerIndividuals=0
BGratioFits=0
BGratioChi2Fits=0
rapBinComb=0
ExtendLegendInX=0
ShiftInX=0
PlotVsComp=0
DrawLatexStuff=1
DrawPreliminary=0

DefaultID=${JobID}
#DefaultID=Psi$[nState-3]S_${NSigma}Sigma_11Dec2012_noRhoFactor
CompareID1=Psi$[nState-3]S_${NSigma}Sigma_11Dec2012_noRhoFactor
CompareID2=MCclosure_July27_Ups1S_MCtruthFineEta_0toP1Sigma
CompareID3=MCclosure_July27_Ups1S_MCtruthFineEta_P1toP3Sigma
CompareID4=MCclosure_Sept9_Ups3S_3DataSig_GEN_pT_Eff_RECOdata
nComp=0

LegendEntryDefID=with_RhoFactor
LegendEntryCompID1=no_RhoFactor
LegendEntryCompID2=1S_0toP1Sigma
LegendEntryCompID3=1S_P1toP3Sigma
LegendEntryCompID4=3Sig_RECOdata_GENeff


nSystematics=0

if [ $nState -eq 4 ]
then
ptBinMin=1
ptBinMax=12
fi

if [ $nState -eq 5 ]
then
ptBinMin=1
ptBinMax=6
fi

### Background Polarization plots
#DefaultID=BG0_Mar19_HighCtauSigCheck3p0
#CompareID1=BG0_Mar19_LSB_LowCtauSigCheck3p0
#CompareID2=BG0_Mar19_RSB_LowCtauSigCheck3p0
#PlotBG0plots=1

### Background Polarization plots Split LSB
#DefaultID=Data_TheGreatARC_1rap_June7_BG0_MassSplit3_RSB4
#CompareID1=Data_TheGreatARC_1rap_June7_BG0_MassSplit3_RSB1
#CompareID2=Data_TheGreatARC_1rap_June7_BG0_MassSplit3_RSB2
#CompareID3=Data_TheGreatARC_1rap_June7_BG0_MassSplit3_RSB3
#PlotBG0plots=1
#rapBinComb=1



##### EtaDepEff_ParamSyst
#SystID1Base=EtaDepEff
#SystID1Specify=MCclosure_July25_Ups3S_MCTnPparam_TO_MCclosure_July25_Ups3S_MCtruthFineEta
#SystID1Title=3S_MC^{truth}_{fine_#eta}-MC^{TnP}_{param}
#
#SystID2Base=EtaDepEff
#SystID2Specify=Toy
#SystID2Title=Reco-MC^{truth}_{fine_#eta}_Extr-MC^{TnP}_{param}
#
#SystID3Base=TheGreatRun_Param
#SystID3Specify=BestSyst
#SystID3Title=ParamSyst


##### MassSigmaDep
#SystID1Base=MassSigmaDep
#SystID1Specify=MCclosure_July25_Ups3S_MCtruthFineEta_1Sig
#SystID1Title=MCclosure3S_1Sig
#
#SystID2Base=MassSigmaDep
#SystID2Specify=OldDataDiff1S
#SystID2Title=1S_OLDdata_1p0Sig-2p5Sig
#
#SystID3Base=MassSigmaDep
#SystID3Specify=OldDataDiff2S
#SystID3Title=2S_OLDdata_1p0Sig-2p5Sig
#
#SystID4Base=MassSigmaDep
#SystID4Specify=OldDataDiff3S
#SystID4Title=3S_OLDdata_1p0Sig-2p5Sig


#SystID1Base=MassSigmaDep
#SystID1Specify=Toy
#SystID1Title=Toy_1Sig


##### TheGreatRun_RunDependence
#SystID1Base=TheGreatRun_RunDependence
#SystID1Specify=SEC_Data_TheGreatCrisis_June14_2Sigma_Run2data_Run2eff_TO_Data_TheGreatCrisis_June14_2Sigma_Run1data_Run1eff
#SystID1Title=Run2_vs_Run1
#
#SystID2Base=TheGreatRun_RunDependence
#SystID2Specify=SEC_Data_TheGreatCrisis_June14_2Sigma_Run3data_Run3eff_TO_Data_TheGreatCrisis_June14_2Sigma_Run1data_Run1eff
#SystID2Title=Run3_vs_Run1
#
#SystID3Base=TheGreatRun_RunDependence
#SystID3Specify=SEC_Data_TheGreatCrisis_June14_2Sigma_Run3data_Run3eff_TO_Data_TheGreatCrisis_June14_2Sigma_Run2data_Run2eff
#SystID3Title=Run3_vs_Run2


##### TheGreatRun_FrameworkII
#SystID1Base=TheGreatRun_FrameworkII
#SystID1Specify=SigDataResults1S
#SystID1Title=1S
#
#SystID2Base=TheGreatRun_FrameworkII
#SystID2Specify=SigDataResults2S
#SystID2Title=2S
#
#SystID3Base=TheGreatRun_FrameworkII
#SystID3Specify=SigDataResults3S
#SystID3Title=3S

##### TheGreatRun_ParamMar19
#SystID1Base=TheGreatRun_ParamMar19
#SystID1Specify=BestSyst
#SystID1Title=SystParam.
#
#SystID2Base=TheGreatRun_TnP
#SystID2Specify=BestSyst
#SystID2Title=SystTnP.

##### DataRatioFitsNEWoverOLD
#SystID1Base=TheGreatRun_DataRatio
#SystID1Specify=June12_2DNEWoverOLDRatio
#SystID1Title=FitToRatioNEWoverOLD

##### BGratioFits
#BGratioFits
#SystID1Base=TheGreatRun_BKGratio
#SystID1Specify=1DfitsAllCells_Comb
#SystID1Title=1S
#
#SystID2Base=TheGreatRun_BKGratio
#SystID2Specify=1DfitsAllCells_Comb
#SystID2Title=2S
#
#SystID3Base=TheGreatRun_BKGratio
#SystID3Specify=1DfitsAllCells_Comb
#SystID3Title=3S


##### SBmSig plots
#SBmSig
#SystID1Base=TheGreatRun_BKGfits
#SystID1Specify=LSBmSig
#SystID1Title=LSB-Sig
#
#SystID2Base=TheGreatRun_BKGfits
#SystID2Specify=RSBmSig
#SystID2Title=RSB-Sig
#SBmSigPlots=1

##### Background model systematics
#BKGmodel
#SystID1Base=TheGreatRun_BKGmodel
#SystID1Specify=Data_TheGreatRun_10B_Apr26_BkgFracLminus28_TO_Data_TheGreatRun_10B_Apr18_NewestCentrals
#SystID1Title=f_{LSB}-0.28
#
#SystID2Base=TheGreatRun_BKGmodel
#SystID2Specify=Data_TheGreatRun_10B_Apr26_BkgFracLplus28_TO_Data_TheGreatRun_10B_Apr18_NewestCentrals
#SystID2Title=f_{LSB}+0.28

##### NEW Background model systematics
#BKGmodel
#SystID1Base=TheGreatRun_BKGmodel
#SystID1Specify=Data_TheGreatRun_10B_May20_NewestCentrals_BGmodel_MINUS_28_46_30_TO_Data_TheGreatRun_10B_May20_NewestCentrals_Original
#SystID1Title=f_{LSB}-0.30
#
#SystID2Base=TheGreatRun_BKGmodel
#SystID2Specify=Data_TheGreatRun_10B_May20_NewestCentrals_BGmodel_PLUS_28_46_30_TO_Data_TheGreatRun_10B_May20_NewestCentrals_Original
#SystID2Title=f_{LSB}+0.30

##### Background frameworkII systematics
#BKGframe
#SystID1Base=TheGreatRun_FrameworkII
#SystID1Specify=Bkg_Scen7
#SystID1Title=#lambda_{#theta}^{PX}_0_#lambda_{#phi}^{PX}_0.8
#
#SystID2Base=TheGreatRun_FrameworkII
#SystID2Specify=Bkg_Scen8
#SystID2Title=#lambda_{#theta}^{PX}_2_#lambda_{#phi}^{PX}_0.4
#
#SystID3Base=TheGreatRun_FrameworkII
#SystID3Specify=Bkg_Scen9
#SystID3Title=#lambda_{#theta}^{PX}_4_#lambda_{#phi}^{PX}_-0.4

##### Signal frameworkII systematics
#FrameworkIIsig
#SystID1Base=TheGreatRun_FrameworkII
#SystID1Specify=Sig_M05_10B
#SystID1Title=#lambda_{#theta}^{PX}_-0.5
#
#SystID2Base=TheGreatRun_FrameworkII
#SystID2Specify=Sig_P05_10B
#SystID2Title=#lambda_{#theta}^{PX}_+0.5

##### FrameworkI systematics
#SigFrameI_true
#SystID1Base=TheGreatRun_FrameworkI
#SystID1Specify=FrameworkI_3S
#SystID1Title=#Upsilon3S
#
#SystID2Base=TheGreatRun_FrameworkI
#SystID2Specify=FrameworkI_2S
#SystID2Title=#Upsilon2S
#
#SystID3Base=TheGreatRun_FrameworkI
#SystID3Specify=FrameworkI_1S
#SystID3Title=#Upsilon1S
#
#SystID4Base=TheGreatRun_FrameworkI
#SystID4Specify=FrameworkI_3S_Median
#SystID4Title=#Upsilon3SMedian
#
#SystID5Base=TheGreatRun_FrameworkI
#SystID5Specify=FrameworkI_2S_Median
#SystID5Title=#Upsilon2SMedian
#
#SystID6Base=TheGreatRun_FrameworkI
#SystID6Specify=FrameworkI_1S_Median
#SystID6Title=#Upsilon1SMedian

##### TnP
#TnP
#SystID1Base=TheGreatRun_TnP
#SystID1Specify=Toy_TheGreatRun_Mar19_1S_TnP
#SystID1Title=TnP_Study

##### Dilep
#DilepMC
#SystID1Base=TheGreatRun_Dilep
#SystID1Specify=DilepMC
#SystID1Title=#mu#mu-Vert.Eff.

##### DeltaTilde
#DeltaTilde
#SystID1Base=TheGreatRun_DeltaTilde
#SystID1Specify=CStoHX
#SystID1Title=#tilde{#lambda}_{CS}-#tilde{#lambda}_{HX}
#
#SystID2Base=TheGreatRun_DeltaTilde
#SystID2Specify=PXtoCS
#SystID2Title=#tilde{#lambda}_{PX}-#tilde{#lambda}_{CS}
#
#SystID3Base=TheGreatRun_DeltaTilde
#SystID3Specify=HXtoPX
#SystID3Title=#tilde{#lambda}_{HX}-#tilde{#lambda}_{PX}
#
#SystID4Base=TotalSyst
#SystID4Specify=TotalSquaredSystMay21
#SystID4Title=TotalSystematic
#
#SystID5Base=TotalSyst
#SystID5Specify=TotalSquaredSystMay21
#SystID5Title=TotalSystematic

##### Parametrization
#Parametrization_true
#SystID1Base=TheGreatRun_Param
#SystID1Specify=effshiftMINUS
#SystID1Title=effshiftMINUS
#
#SystID2Base=TheGreatRun_Param
#SystID2Specify=effshiftPLUS
#SystID2Title=effshiftPLUS
#
#SystID3Base=TheGreatRun_Param
#SystID3Specify=pTshiftMINUS
#SystID3Title=pTshiftMINUS
#
#SystID4Base=TheGreatRun_Param
#SystID4Specify=pTshiftPLUS
#SystID4Title=pTshiftPLUS
#
#SystID5Base=TheGreatRun_Param
#SystID5Specify=pTscaleMINUS
#SystID5Title=pTscaleMINUS
#
#SystID6Base=TheGreatRun_Param
#SystID6Specify=pTscalePLUS
#SystID6Title=pTscalePLUS


###### BGratio
#SystID1Base=TheGreatRun_BKGmodel
#SystID1Specify=BGratioFit
#SystID1Title=1S
#
#SystID2Base=TheGreatRun_BKGmodel
#SystID2Specify=BGratioFit
#SystID2Title=2S
#
#SystID3Base=TheGreatRun_BKGmodel
#SystID3Specify=BGratioFit
#SystID3Title=3S

###### BGratioChi2
#SystID1Base=TheGreatRun_BKGratio
#SystID1Specify=Chi2_SBratio_May24_2D_AbsCosthMax1_with1DPlots
#SystID1Title=asdf

###### NORMALS
#SystID1Base=TheGreatRun_ParamMar19
#SystID1Specify=BestSyst
#SystID1Title=Syst.Old
#
#SystID2Base=TheGreatRun_Param
#SystID2Specify=BestSyst
#SystID2Title=Syst.New
#
#SystID3Base=TheGreatRun_DeltaTilde
#SystID3Specify=HXtoPX
#SystID3Title=#tilde{#lambda}_{HX}-#tilde{#lambda}_{PX}

#SystID4Base=TotalSyst
#SystID4Specify=TotalSquaredSystApr27
#SystID4Title=TotalSystematic
#
#SystID5Base=TotalSyst
#SystID5Specify=TotalSquaredSystApr27
#SystID5Title=TotalSystematic
#
#SystID6Base=TheGreatRun_Param
#SystID6Specify=effshiftMINUS_10B
#SystID6Title=effshiftMINUS
#
#SystID7Base=TheGreatRun_BKGMassScan
#SystID7Specify=MassScan11
#SystID7Title=MassScan11
#
#SystID8Base=TheGreatRun_BKGMassScan
#SystID8Specify=MassScan12
#SystID8Title=MassScan12




########################################

cd ${homedir}

cp ../../interface/ToyMC_Psi$[nState-3]S.h ToyMC.h
cp ../../interface/effsAndCuts_Psi$[nState-3]S.h effsAndCuts.h
cp ../../interface/rootIncludes.inc rootIncludes.inc
cp ../../interface/commonVar_Psi$[nState-3]S.h commonVar.h

touch PlotFinalResults.cc
make

mkdir -p FinalResults/${JobID}/Psi$[nState-3]S

JobIDDir=FinalResults/${JobID}

mkdir -p ${JobIDDir}


cp PlotFinalResults PlotFinalResults_Psi$[nState-3]S
./PlotFinalResults_Psi$[nState-3]S ${DefaultID}=DefaultID ${CompareID1}=CompareID1 ${CompareID2}=CompareID2 ${CompareID3}=CompareID3 ${CompareID4}=CompareID4 ${JobID}=JobID ${SystID1Base}=SystID1Base ${SystID1Specify}=SystID1Specify ${SystID1Title}=SystID1Title ${SystID2Base}=SystID2Base ${SystID2Specify}=SystID2Specify ${SystID2Title}=SystID2Title ${basedir}=basedir ${storagedir}=storagedir ${ptBinMin}ptBinMin ${ptBinMax}ptBinMax ${nSystematics}nSystematics ${nComp}nComp ${nState}nState ${SystID3Base}=SystID3Base ${SystID3Specify}=SystID3Specify ${SystID3Title}=SystID3Title ${SystID4Base}=SystID4Base ${SystID4Specify}=SystID4Specify ${SystID4Title}=SystID4Title ${SystID5Base}=SystID5Base ${SystID5Specify}=SystID5Specify ${SystID5Title}=SystID5Title ${SystID6Base}=SystID6Base ${SystID6Specify}=SystID6Specify ${SystID6Title}=SystID6Title ${SystID7Base}=SystID7Base ${SystID7Specify}=SystID7Specify ${SystID7Title}=SystID7Title ${SystID8Base}=SystID8Base ${SystID8Specify}=SystID8Specify ${SystID8Title}=SystID8Title PlotMatt=${PlotMatt} PlotAsymm=${PlotAsymm} PlotCompare=${PlotCompare} PlotFinalData=${PlotFinalData} PlotSystematics=${PlotSystematics} PlotLegend=${PlotLegend} PlotBG0plots=${PlotBG0plots} DeltaTildeplots=${DeltaTildeplots} SBmSigPlots=${SBmSigPlots} CompareSyst=${CompareSyst} SteerIndividuals=${SteerIndividuals} BGratioFits=${BGratioFits} BGratioChi2Fits=${BGratioChi2Fits} rapBinComb=${rapBinComb} SetCompStyle=${SetCompStyle} ${LegendEntryDefID}=LegendEntryDefID ${LegendEntryCompID1}=LegendEntryCompID1 ${LegendEntryCompID2}=LegendEntryCompID2 ${LegendEntryCompID3}=LegendEntryCompID3 ${LegendEntryCompID4}=LegendEntryCompID4 ExtendLegendInX=${ExtendLegendInX} ShiftInX=${ShiftInX} PlotVsComp=${PlotVsComp} DrawLatexStuff=${DrawLatexStuff} DrawPreliminary=${DrawPreliminary}
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
cd ../
fi

if [ ${PlotSystematics} -eq 1 ]
then
cp ${basedir}/latex/Systematics_Psi$[nState-3]S.tex ./Systematics.tex
pdflatex Systematics.tex
mv Systematics.pdf Systematics${additionalName}.pdf
fi
rm *.aux
rm *.log
rm *.tex

#rm -r Figures${additionalName}
mkdir Figures${additionalName}
mv Figures/* Figures${additionalName}/

done

