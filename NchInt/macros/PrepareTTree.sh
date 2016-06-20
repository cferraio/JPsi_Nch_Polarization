#!/bin/sh

source /cvmfs/sft.cern.ch/lcg/app/releases/ROOT/5.34.05/x86_64-slc5-gcc43-opt/root/bin/thisroot.sh                                     
source /cvmfs/sft.cern.ch/lcg/app/releases/ROOT/5.34.05/x86_64-slc5-gcc43-opt/root/bin/setxrd.sh /cvmfs/sft.cern.ch/lcg/external/xrootd/3.2.4/x86_64-slc5-gcc46-opt/

#homedir=HOMEDIR
#cd ${homedir}
Cdir=$PWD

cd ..
basedir=$PWD
cd macros

# input arguments
for nState in 4;do    #1,2,3,Upsi(1S,2S,3S); 4,Jpsi 5,PsiPrime
for FidCuts in 11;do #defines the set of cuts to be used, see macros/polFit/effsAndCuts.h
cd $Cdir


rapMin=1     #takes bins, not actual values
rapMax=2     #if you only want to process 1 y bin, rapMax = rapMin
ptMin=1      #takes bins, not acutal values
ptMax=12      #if you only want to process 1 pt bin, ptMax = ptMin
Plotting=2   #plotting macro: 1 = plot all, 2 = plot mass, 3 = plot lifetime sidebands, 4 = plot lifetime singal region, 
	           # 5 = sidebands, separate pull and distribution, 6 = signal region, separate pull and distribution

rejectCowboys=true
RequestTrigger=true
MC=false
correctCtau=false   #correct pseudo-proper lifetime to l_new = l * MpsiPDG / Mpsi, with l = Lxy * Mpsi / pT
drawRapPt2D=false  #draw Rap-Pt 2D map of Psi

doCtauUncer=true
PolLSB=false       #measure polarization of the left sideband
PolRSB=false       #measure polarization of the right sideband
PolNP=false        #measure polarization of the non prompt events
forceBinning=true  #set binning of Psi1S consistently to non prompt binning and Psi2S consistently to background binning
folding=true       #folding is applied to all background histograms
normApproach=false #normalization 
ctauScen=0         #0:default(1s:2.5,2s:2.0), 1:(1s:3.5,2s:3.0), 2:(1s:1.5,2s:1.0), 3:100mm 1S and 2S, 
FracLSB=-1         #-1:defalut, 0, 100
scaleFracBg=false
fitMassPR=false
fitMassNP=false

DataID=Psi$[nState-3]S_ctauScen0_FracLSB-1_16Mar2013
polDataPath=${basedir}/Psi/Data/${DataID}

#Define JobID
JobID=ctauScen${ctauScen}_FracLSB${FracLSB}_2012higherptnoDoubleMu
#JobID=ctauScen${ctauScen}_FracLSB${FracLSB}_2012Run1withCowboys
#JobID=ctauScen${ctauScen}_FracLSB${FracLSB}_2011Data
#JobID=ctauScen${ctauScen}_FracLSB${FracLSB}_newMLfit_4Mar2013_NP
#JobID=ctauScen${ctauScen}_FracLSB${FracLSB}_newMLfit_4Mar2013
#JobID=ctauScen${ctauScen}_FracLSB${FracLSB}_newMLfit_4Mar2013_0fracBg
#JobID=ctauScen${ctauScen}_FracLSB${FracLSB}_newMLfit_4Mar2013_1sigMass
#JobID=ctauScen${ctauScen}_FracLSB${FracLSB}_newMLfit_correctCtau_11April2013

# input files
# In case of more input Files: define inputTreeX and adapt the line starting with inputTrees, at the moment up to 4 files implemented
if [ ${nState} -eq 4 ] 
then
#inputTree1=/data/users/ferraioc/TTree_Onia2MuMu_v30_PromptRecoAB_10May2012_Jpsi.root #2011 file
#inputTree1=/data/users/ferraioc/Polarization/2012ppOniaData/2012ABCDMuOnia_ADoubMu_MuPk_jpsi_v8.root
#inputTree1=/data/users/ferraioc/Polarization/2012ppOniaData/2012ABCDMuOnia_ADoubMu_jpsi_v8.root
#inputTree1=/data/users/ferraioc/Polarization/2012ppOniaData/r2012D_MuPk_jpsi_v8_1.root
#inputTree1=/data/users/ferraioc/Polarization/2012ppOniaData/r2012A_MuOnia_jpsi_v8.root
#inputTree2=/data/users/ferraioc/Polarization/2012ppOniaData/r2012B_MuOnia_jpsi_v8.root
#inputTree3=/data/users/ferraioc/Polarization/2012ppOniaData/r2012C_MuOnia_jpsi_v8v2.root
#inputTree4=/data/users/ferraioc/Polarization/2012ppOniaData/r2012D_MuOnia_jpsi_v8v2.root
inputTree1=/data/users/ferraioc/Polarization/2012ppOniaData/r2012A_DoubMu_jpsi_v8.root
inputTree2=/data/users/ferraioc/Polarization/2012ppOniaData/r2012C_MuPk_jpsi_v8_1.root
inputTree3=/data/users/ferraioc/Polarization/2012ppOniaData/r2012D_MuPk_jpsi_v8_1.root
inputTree4=/data/users/ferraioc/Polarization/2012ppOniaData/r2012B_MuPk_jpsi_v8.root
if [ ${MC} = 'true' ]
then
inputTree1=/scratch/ikratsch/Polarization/Jpsi/InputFiles/TTree_Psi1S_Gun_Pt9p5_70p5_19Dec2012.root
fi
fi

if [ ${nState} -eq 5 ]
then
#inputTree1=/data/users/ferraioc/Polarization/2012ppOniaData/r2012A_DoubMu_psi2s_v8.root
#inputTree1=/data/users/ferraioc/Polarization/2012ppOniaData/r2012A_MuOnia_psi2s_v8.root
#inputTree1=/data/users/ferraioc/Polarization/2012ppOniaData/r2012B_MuOnia_psi2s_v8.root
inputTree1=/data/users/ferraioc/TTree_Onia2MuMu_v30_PromptRecoAB_10May2012_Psi.root
if [ ${MC} = 'true' ]
then
inputTree1=/scratch/ikratsch/Polarization/Jpsi/InputFiles/TTree_Psi2S_Gun_Pt6p5_50p5_19Dec2012.root
fi
fi

################ EXECUTABLES #################

#following flags decide if the step is executed (1) or not (0):
#IMPORTANT: for MC set execute_runWorkspace, execute_MassFit and execute_runLifetimeFit to 0
execute_runData=1			           #independent of rapMin, rapMax, ptMin, ptMax
execute_runWorkspace=0			     #independent of rapMin, rapMax, ptMin, ptMax
execute_runMassFit=0			       #can be executed for different pt and y bins
execute_runLifetimeFit=0         #can be executed for different pt and y bins
execute_runPlotMassLifetime=0    #can be executed for different pt and y bins
execut_PlotFitPar=0              #independent of rapMin, rapMax, ptMin, ptMax
execute_runBkgHistos=0           #can be executed for different pt and y bins
execute_PlotCosThetaPhiBG=0 		 #This step only has to be executed once for each set of cuts (indep. of FracLSB and nSigma)
execute_PlotCosThetaPhiDistribution=0 #This step only has to be executed once for each set of cuts (indep. of FracLSB and nSigma)

#################################

# Make directories
CutDir=${Cdir}/DataFiles/SetOfCuts${FidCuts}_${JobID}

WorkDir=${CutDir}/Psi$[nState-3]S
mkdir -p ${CutDir}
mkdir -p ${WorkDir}
cp ../interface/commonVar_Psi$[nState-3]S.h ${WorkDir}/commonVar.h 

mkdir -p DataFiles
mkdir -p ${WorkDir}/tmpFiles/backupWorkSpace
mkdir -p ${WorkDir}/Figures
mkdir -p ${WorkDir}/PDF
mkdir -p ${WorkDir}/Fit

# Copy files to directory
cp Makefile ${WorkDir}/Makefile
cp ../interface/rootIncludes.inc ${WorkDir}/rootIncludes.inc

cp runData.cc ${WorkDir}/runData.cc
cp PolData.C ${WorkDir}/PolData.C
cp PolData.h ${WorkDir}/PolData.h
cp ../interface/effsAndCuts_Psi$[nState-3]S.h ${WorkDir}/effsAndCuts.h

cp runWorkspace.cc ${WorkDir}/runWorkspace.cc
cp createWorkspace.C ${WorkDir}/createWorkspace.C

cp runMassFit.cc ${WorkDir}/runMassFit.cc
cp massFit.cc ${WorkDir}/massFit.cc

cp runLifetimeFit.cc ${WorkDir}/runLifetimeFit.cc
cp lifetimeFit.cc ${WorkDir}/lifetimeFit.cc
cp ../interface/calculatePar.cc ${WorkDir}/calculatePar.cc
cp ../interface/RooUtils.h ${WorkDir}/RooUtils.h

cp runPlotMassLifetime.cc ${WorkDir}/runPlotMassLifetime.cc
cp PlotMassLifetime.cc ${WorkDir}/PlotMassLifetime.cc

cp PlotFitPar.cc ${WorkDir}/PlotFitPar.cc

cp runBkgHistos.cc ${WorkDir}/runBkgHistos.cc
cp bkgHistos.C ${WorkDir}/bkgHistos.C
cp calcPol.C ${WorkDir}/calcPol.C

cp PlotCosThetaPhiBG.cc ${WorkDir}/PlotCosThetaPhiBG.cc
cp PlotCosThetaPhiDistribution.cc ${WorkDir}/PlotCosThetaPhiDistribution.cc

cp ../latex/Mass_fitParameter.tex ${WorkDir}/Mass_fitParameter.tex
cp ../latex/Lifetime_fitParameter.tex ${WorkDir}/Lifetime_fitParameter.tex
cp ../latex/myStyle.tex ${WorkDir}/myStyle.tex
cp ../latex/evaluateCtau.tex ${WorkDir}/evaluateCtau.tex
cp ../latex/NumEvents.tex ${WorkDir}/NumEvents.tex

cp ../latex/cosThetaPhi_$[nState-3]S_BG.tex        ${WorkDir}/cosThetaPhi_$[nState-3]S_BG.tex
cp ../latex/cosThetaPhi_$[nState-3]S_BG_highct.tex ${WorkDir}/cosThetaPhi_$[nState-3]S_BG_highct.tex
cp ../latex/cosThetaPhi_$[nState-3]S_NPBG.tex      ${WorkDir}/cosThetaPhi_$[nState-3]S_NPBG.tex
cp ../latex/cosThetaPhi_$[nState-3]S_TBG.tex       ${WorkDir}/cosThetaPhi_$[nState-3]S_TBG.tex
cp ../latex/cosThetaPhi_$[nState-3]S.tex           ${WorkDir}/cosThetaPhi_$[nState-3]S.tex
cp ../latex/MassLifetime_Psi$[nState-3]S.tex       ${WorkDir}/MassLifetime_Psi$[nState-3]S.tex

cd ${WorkDir}

make

inputTrees="inputTree=${inputTree1} inputTree=${inputTree2} inputTree=${inputTree3} inputTree=${inputTree4} inputTree=${inputTree5} inputTree=${inputTree6} inputTree=${inputTree7} inputTree=${inputTree8}"
if [ ${execute_runData} -eq 1 ]
then
./runData ${inputTrees} rejectCowboys=${rejectCowboys} FidCuts=${FidCuts} nState=${nState} MC=${MC} RequestTrigger=${RequestTrigger}
fi

if [ ${execute_runWorkspace} -eq 1 ]
then
./runWorkspace nState=${nState} correctCtau=${correctCtau} drawRapPt2D=${drawRapPt2D}
fi

if [ ${execute_runMassFit} -eq 1 ]
then
rootfile=fit_Psi$[nState-3]S_rap${rapMin}_pt${ptMin}.root
cp tmpFiles/backupWorkSpace/fit_Psi$[nState-3]S* tmpFiles/
cp runMassFit runMassFit_$[nState-3]S_rap${rapMin}_pt${ptMin}
./runMassFit_$[nState-3]S_rap${rapMin}_pt${ptMin} rapMin=${rapMin} rapMax=${rapMax} ptMin=${ptMin} ptMax=${ptMax} nState=${nState} fitMassPR=${fitMassPR} fitMassNP=${fitMassNP}
rm runMassFit_$[nState-3]S_rap${rapMin}_pt${ptMin}
fi

if [ ${execute_runLifetimeFit} -eq 1 ]
then
cp runLifetimeFit runLifetimeFit_$[nState-3]S_rap${rapMin}_pt${ptMin}
./runLifetimeFit_$[nState-3]S_rap${rapMin}_pt${ptMin} rapMin=${rapMin} rapMax=${rapMax} ptMin=${ptMin} ptMax=${ptMax} nState=${nState}
rm runLifetimeFit_$[nState-3]S_rap${rapMin}_pt${ptMin}
fi

if [ ${execute_runPlotMassLifetime} -eq 1 ]
then
cp runPlotMassLifetime runPlotMassLifetime_$[nState-3]S_rap${rapMin}_pt${ptMin}
./runPlotMassLifetime_$[nState-3]S_rap${rapMin}_pt${ptMin} rapMin=${rapMin} rapMax=${rapMax} ptMin=${ptMin} ptMax=${ptMax} nState=${nState} Plotting=${Plotting}
rm runPlotMassLifetime_$[nState-3]S_rap${rapMin}_pt${ptMin}
#pdflatex MassLifetime_Psi$[nState-3]S.tex
#mv MassLifetime_Psi$[nState-3]S.pdf PDF/MassLifetime_Psi$[nState-3]S.pdf
fi

if [ ${execut_PlotFitPar} -eq 1 ]
then
./PlotFitPar nState=${nState} doCtauUncer=${doCtauUncer}
#pdflatex Lifetime_fitParameter.tex
pdflatex Mass_fitParameter.tex
#pdflatex evaluateCtau.tex
#pdflatex evaluateCtau.tex
#pdflatex NumEvents.tex
#pdflatex NumEvents.tex
#mv Lifetime_fitParameter.pdf PDF/Lifetime_fitParameter.pdf
mv Mass_fitParameter.pdf PDF/Mass_fitParameter.pdf
#mv evaluateCtau.pdf PDF/evaluateCtau.pdf
#mv NumEvents.pdf PDF/NumEvents.pdf
fi

if [ ${execute_runBkgHistos} -eq 1 ]
then
cp runBkgHistos runBkgHistos_$[nState-3]S_rap${rapMin}_pt${ptMin}
./runBkgHistos_$[nState-3]S_rap${rapMin}_pt${ptMin} rapMin=${rapMin} rapMax=${rapMax} ptMin=${ptMin} ptMax=${ptMax} nState=${nState} MC=${MC} doCtauUncer=${doCtauUncer} PolLSB=${PolLSB} PolRSB=${PolRSB} PolNP=${PolNP} ctauScen=${ctauScen} FracLSB=${FracLSB} forceBinning=${forceBinning} folding=${folding} normApproach=${normApproach} scaleFracBg=${scaleFracBg} ${polDataPath}=polDataPath
rm runBkgHistos_$[nState-3]S_rap${rapMin}_pt${ptMin}
fi

if [ ${execute_PlotCosThetaPhiBG} -eq 1 ]
then
./PlotCosThetaPhiBG nState=${nState}
pdflatex cosThetaPhi_$[nState-3]S_BG.tex
pdflatex cosThetaPhi_$[nState-3]S_BG_highct.tex
pdflatex cosThetaPhi_$[nState-3]S_NPBG.tex
pdflatex cosThetaPhi_$[nState-3]S_TBG.tex
mv cosThetaPhi_$[nState-3]S_BG.pdf PDF/cosThetaPhi_$[nState-3]S_BG.pdf
mv cosThetaPhi_$[nState-3]S_BG_highct.pdf PDF/cosThetaPhi_$[nState-3]S_BG_highct.pdf
mv cosThetaPhi_$[nState-3]S_NPBG.pdf PDF/cosThetaPhi_$[nState-3]S_NPBG.pdf
mv cosThetaPhi_$[nState-3]S_TBG.pdf PDF/cosThetaPhi_$[nState-3]S_TBG.pdf
fi

if [ ${execute_PlotCosThetaPhiDistribution} -eq 1 ]
then
./PlotCosThetaPhiDistribution ${nState}nState ${WorkDir}=DataPath
pdflatex cosThetaPhi_$[nState-3]S.tex
mv cosThetaPhi_$[nState-3]S.pdf PDF/cosThetaPhi_$[nState-3]S.pdf
fi

#rm runData
#rm runWorkspace
#rm runMassFit
#rm runLifetimeFit
#rm runPlotMassLifetime
#rm runBkgHistos
#rm PlotFitPar
#rm PlotCosThetaPhiBG
#rm PlotCosThetaPhiDistribution
#rm *.tex
rm *.aux
rm *.log
rm *.so
rm *.d
rm *.nav 
rm *.out 
rm *.snm 
rm *.toc 

done
done

