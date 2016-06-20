/*
 * PlotFinalResults.cc
 *
 *  Created on: Dec 3, 2011
 *      Author: valentinknuenz
 *
 *  Adapted to Psi(1S,2S): Dec 12 2012, LinlinZhang
 *  Upsilon(1S,2S,3S): nState = 1,2,3
 *  Psi(1S,2S):        nState = 4,5
 */

#include <iostream>
#include <string>
#include <sstream>
using namespace std;

#include "rootIncludes.inc"
#include "commonVar.h"
#include "ToyMC.h"

int main(int argc, char** argv) {

	Char_t *JobID = "JobID"; //Storage Directory

	Char_t *SystID1Base = "Default"; //Storage Directory
	Char_t *SystID1Specify = "Default";
	Char_t *SystID1Title = "Default";

	Char_t *SystID2Base = "Default"; //Storage Directory
	Char_t *SystID2Specify = "Default";
	Char_t *SystID2Title = "Default";

	Char_t *SystID3Base = "Default"; //Storage Directory
	Char_t *SystID3Specify = "Default";
	Char_t *SystID3Title = "Default";

	Char_t *SystID4Base = "Default"; //Storage Directory
	Char_t *SystID4Specify = "Default";
	Char_t *SystID4Title = "Default";

	Char_t *SystID5Base = "Default"; //Storage Directory
	Char_t *SystID5Specify = "Default";
	Char_t *SystID5Title = "Default";

	Char_t *SystID6Base = "Default"; //Storage Directory
	Char_t *SystID6Specify = "Default";
	Char_t *SystID6Title = "Default";

	Char_t *SystID7Base = "Default"; //Storage Directory
	Char_t *SystID7Specify = "Default";
	Char_t *SystID7Title = "Default";

	Char_t *SystID8Base = "Default"; //Storage Directory
	Char_t *SystID8Specify = "Default";
	Char_t *SystID8Title = "Default";

	Char_t *basedir = "Default";
	Char_t *storagedir = "Default";
	Char_t *DefaultID = "Default";
	Char_t *CompareID1 = "CompareID1";
	Char_t *CompareID2 = "CompareID2";
	Char_t *CompareID3 = "CompareID3";
	Char_t *CompareID4 = "CompareID4";

	Char_t *LegendEntryDefID  = "Default";
	Char_t *LegendEntryCompID1= "Default";
	Char_t *LegendEntryCompID2= "Default";
	Char_t *LegendEntryCompID3= "Default";
	Char_t *LegendEntryCompID4= "Default";

	Char_t *MPCentralsWithTotalSystID = "MPCentralsWithTotalSystID";

	int ptBinMin=1;
	int ptBinMax=1;
	int nSystematics=1;
	int nComp=1;
	int nState=1;

	bool PlotMatt(false);
	bool PlotAsymm(false);
	bool PlotCompare(false);
	bool PlotFinalData(false);
	bool PlotSystematics(false);
	bool PlotLegend(false);
	bool PlotBrazilian(false);
	bool FitGraph(false);
	bool DrawLatexStuff(false);
	bool MultiPanelPlots(false);
	bool PlotBG0plots(false);
	bool DeltaTildeplots(false);
	bool SBmSigPlots(false);
	bool PlotAlteredPPDResults(false);
	bool CompareSyst(false);
	bool SteerIndividuals(false);
	bool BGratioFits(false);
	bool BGratioChi2Fits(false);
	bool rapBinComb(false);
	bool SetCompStyle(false);
	bool DrawPreliminary(true);
	bool PlotMattForICHEP(false);
	bool ExtendLegendInX(false);
	bool ShiftInX(true);
	bool ShiftCompareInX(false);
	bool PlotVsComp(false);
	bool PlotSysSquare(false);
	bool PlotCL1sigma(false);

	for( int i=0;i < argc; ++i ) {

		if(std::string(argv[i]).find("JobID") != std::string::npos) {char* JobIDchar = argv[i]; char* JobIDchar2 = strtok (JobIDchar, "="); JobID = JobIDchar2; cout<<"JobID = "<<JobID<<endl;}
		if(std::string(argv[i]).find("DefaultID") != std::string::npos) {char* DefaultIDchar = argv[i]; char* DefaultIDchar2 = strtok (DefaultIDchar, "="); DefaultID = DefaultIDchar2; cout<<"DefaultID = "<<DefaultID<<endl;}
		if(std::string(argv[i]).find("basedir") != std::string::npos) {char* basedirchar = argv[i]; char* basedirchar2 = strtok (basedirchar, "="); basedir = basedirchar2; cout<<"basedir = "<<basedir<<endl;}
		if(std::string(argv[i]).find("storagedir") != std::string::npos) {char* storagedirchar = argv[i]; char* storagedirchar2 = strtok (storagedirchar, "="); storagedir = storagedirchar2; cout<<"storagedir = "<<storagedir<<endl;}
		if(std::string(argv[i]).find("SystID1Base") != std::string::npos) {char* SystID1Basechar = argv[i]; char* SystID1Basechar2 = strtok (SystID1Basechar, "="); SystID1Base = SystID1Basechar2; cout<<"SystID1Base = "<<SystID1Base<<endl;}
		if(std::string(argv[i]).find("SystID1Specify") != std::string::npos) {char* SystID1Specifychar = argv[i]; char* SystID1Specifychar2 = strtok (SystID1Specifychar, "="); SystID1Specify = SystID1Specifychar2; cout<<"SystID1Specify = "<<SystID1Specify<<endl;}
		if(std::string(argv[i]).find("SystID1Title") != std::string::npos) {char* SystID1Titlechar = argv[i]; char* SystID1Titlechar2 = strtok (SystID1Titlechar, "="); SystID1Title = SystID1Titlechar2; cout<<"SystID1Title = "<<SystID1Title<<endl;}
		if(std::string(argv[i]).find("SystID2Base") != std::string::npos) {char* SystID2Basechar = argv[i]; char* SystID2Basechar2 = strtok (SystID2Basechar, "="); SystID2Base = SystID2Basechar2; cout<<"SystID2Base = "<<SystID2Base<<endl;}
		if(std::string(argv[i]).find("SystID2Specify") != std::string::npos) {char* SystID2Specifychar = argv[i]; char* SystID2Specifychar2 = strtok (SystID2Specifychar, "="); SystID2Specify = SystID2Specifychar2; cout<<"SystID2Specify = "<<SystID2Specify<<endl;}
		if(std::string(argv[i]).find("SystID2Title") != std::string::npos) {char* SystID2Titlechar = argv[i]; char* SystID2Titlechar2 = strtok (SystID2Titlechar, "="); SystID2Title = SystID2Titlechar2; cout<<"SystID2Title = "<<SystID2Title<<endl;}
		if(std::string(argv[i]).find("SystID3Base") != std::string::npos) {char* SystID3Basechar = argv[i]; char* SystID3Basechar2 = strtok (SystID3Basechar, "="); SystID3Base = SystID3Basechar2; cout<<"SystID3Base = "<<SystID3Base<<endl;}
		if(std::string(argv[i]).find("SystID3Specify") != std::string::npos) {char* SystID3Specifychar = argv[i]; char* SystID3Specifychar2 = strtok (SystID3Specifychar, "="); SystID3Specify = SystID3Specifychar2; cout<<"SystID3Specify = "<<SystID3Specify<<endl;}
		if(std::string(argv[i]).find("SystID3Title") != std::string::npos) {char* SystID3Titlechar = argv[i]; char* SystID3Titlechar2 = strtok (SystID3Titlechar, "="); SystID3Title = SystID3Titlechar2; cout<<"SystID3Title = "<<SystID3Title<<endl;}
		if(std::string(argv[i]).find("SystID4Base") != std::string::npos) {char* SystID4Basechar = argv[i]; char* SystID4Basechar2 = strtok (SystID4Basechar, "="); SystID4Base = SystID4Basechar2; cout<<"SystID4Base = "<<SystID4Base<<endl;}
		if(std::string(argv[i]).find("SystID4Specify") != std::string::npos) {char* SystID4Specifychar = argv[i]; char* SystID4Specifychar2 = strtok (SystID4Specifychar, "="); SystID4Specify = SystID4Specifychar2; cout<<"SystID4Specify = "<<SystID4Specify<<endl;}
		if(std::string(argv[i]).find("SystID4Title") != std::string::npos) {char* SystID4Titlechar = argv[i]; char* SystID4Titlechar2 = strtok (SystID4Titlechar, "="); SystID4Title = SystID4Titlechar2; cout<<"SystID4Title = "<<SystID4Title<<endl;}
		if(std::string(argv[i]).find("SystID5Base") != std::string::npos) {char* SystID5Basechar = argv[i]; char* SystID5Basechar2 = strtok (SystID5Basechar, "="); SystID5Base = SystID5Basechar2; cout<<"SystID5Base = "<<SystID5Base<<endl;}
		if(std::string(argv[i]).find("SystID5Specify") != std::string::npos) {char* SystID5Specifychar = argv[i]; char* SystID5Specifychar2 = strtok (SystID5Specifychar, "="); SystID5Specify = SystID5Specifychar2; cout<<"SystID5Specify = "<<SystID5Specify<<endl;}
		if(std::string(argv[i]).find("SystID5Title") != std::string::npos) {char* SystID5Titlechar = argv[i]; char* SystID5Titlechar2 = strtok (SystID5Titlechar, "="); SystID5Title = SystID5Titlechar2; cout<<"SystID5Title = "<<SystID5Title<<endl;}
		if(std::string(argv[i]).find("SystID6Base") != std::string::npos) {char* SystID6Basechar = argv[i]; char* SystID6Basechar2 = strtok (SystID6Basechar, "="); SystID6Base = SystID6Basechar2; cout<<"SystID6Base = "<<SystID6Base<<endl;}
		if(std::string(argv[i]).find("SystID6Specify") != std::string::npos) {char* SystID6Specifychar = argv[i]; char* SystID6Specifychar2 = strtok (SystID6Specifychar, "="); SystID6Specify = SystID6Specifychar2; cout<<"SystID6Specify = "<<SystID6Specify<<endl;}
		if(std::string(argv[i]).find("SystID6Title") != std::string::npos) {char* SystID6Titlechar = argv[i]; char* SystID6Titlechar2 = strtok (SystID6Titlechar, "="); SystID6Title = SystID6Titlechar2; cout<<"SystID6Title = "<<SystID6Title<<endl;}
		if(std::string(argv[i]).find("SystID7Base") != std::string::npos) {char* SystID7Basechar = argv[i]; char* SystID7Basechar2 = strtok (SystID7Basechar, "="); SystID7Base = SystID7Basechar2; cout<<"SystID7Base = "<<SystID7Base<<endl;}
		if(std::string(argv[i]).find("SystID7Specify") != std::string::npos) {char* SystID7Specifychar = argv[i]; char* SystID7Specifychar2 = strtok (SystID7Specifychar, "="); SystID7Specify = SystID7Specifychar2; cout<<"SystID7Specify = "<<SystID7Specify<<endl;}
		if(std::string(argv[i]).find("SystID7Title") != std::string::npos) {char* SystID7Titlechar = argv[i]; char* SystID7Titlechar2 = strtok (SystID7Titlechar, "="); SystID7Title = SystID7Titlechar2; cout<<"SystID7Title = "<<SystID7Title<<endl;}
		if(std::string(argv[i]).find("SystID8Base") != std::string::npos) {char* SystID8Basechar = argv[i]; char* SystID8Basechar2 = strtok (SystID8Basechar, "="); SystID8Base = SystID8Basechar2; cout<<"SystID8Base = "<<SystID8Base<<endl;}
		if(std::string(argv[i]).find("SystID8Specify") != std::string::npos) {char* SystID8Specifychar = argv[i]; char* SystID8Specifychar2 = strtok (SystID8Specifychar, "="); SystID8Specify = SystID8Specifychar2; cout<<"SystID8Specify = "<<SystID8Specify<<endl;}
		if(std::string(argv[i]).find("SystID8Title") != std::string::npos) {char* SystID8Titlechar = argv[i]; char* SystID8Titlechar2 = strtok (SystID8Titlechar, "="); SystID8Title = SystID8Titlechar2; cout<<"SystID8Title = "<<SystID8Title<<endl;}

		if(std::string(argv[i]).find("CompareID1") != std::string::npos) {char* CompareID1char = argv[i]; char* CompareID1char2 = strtok (CompareID1char, "="); CompareID1 = CompareID1char2; cout<<"CompareID1 = "<<CompareID1<<endl;}
		if(std::string(argv[i]).find("CompareID2") != std::string::npos) {char* CompareID2char = argv[i]; char* CompareID2char2 = strtok (CompareID2char, "="); CompareID2 = CompareID2char2; cout<<"CompareID2 = "<<CompareID2<<endl;}
		if(std::string(argv[i]).find("CompareID3") != std::string::npos) {char* CompareID3char = argv[i]; char* CompareID3char2 = strtok (CompareID3char, "="); CompareID3 = CompareID3char2; cout<<"CompareID3 = "<<CompareID3<<endl;}
		if(std::string(argv[i]).find("CompareID4") != std::string::npos) {char* CompareID4char = argv[i]; char* CompareID4char2 = strtok (CompareID4char, "="); CompareID4 = CompareID4char2; cout<<"CompareID4 = "<<CompareID4<<endl;}
		if(std::string(argv[i]).find("MPCentralsWithTotalSystID") != std::string::npos) {char* MPCentralsWithTotalSystIDchar = argv[i]; char* MPCentralsWithTotalSystIDchar2 = strtok (MPCentralsWithTotalSystIDchar, "="); MPCentralsWithTotalSystID = MPCentralsWithTotalSystIDchar2; cout<<"MPCentralsWithTotalSystID = "<<MPCentralsWithTotalSystID<<endl;}

		if(std::string(argv[i]).find("LegendEntryDefID") != std::string::npos) {char* LegendEntryDefIDchar = argv[i]; char* LegendEntryDefIDchar2 = strtok (LegendEntryDefIDchar, "="); LegendEntryDefID = LegendEntryDefIDchar2; cout<<"LegendEntryDefID = "<<LegendEntryDefID<<endl;}
		if(std::string(argv[i]).find("LegendEntryCompID1") != std::string::npos) {char* LegendEntryCompID1char = argv[i]; char* LegendEntryCompID1char2 = strtok (LegendEntryCompID1char, "="); LegendEntryCompID1 = LegendEntryCompID1char2; cout<<"LegendEntryCompID1 = "<<LegendEntryCompID1<<endl;}
		if(std::string(argv[i]).find("LegendEntryCompID2") != std::string::npos) {char* LegendEntryCompID2char = argv[i]; char* LegendEntryCompID2char2 = strtok (LegendEntryCompID2char, "="); LegendEntryCompID2 = LegendEntryCompID2char2; cout<<"LegendEntryCompID2 = "<<LegendEntryCompID2<<endl;}
		if(std::string(argv[i]).find("LegendEntryCompID3") != std::string::npos) {char* LegendEntryCompID3char = argv[i]; char* LegendEntryCompID3char2 = strtok (LegendEntryCompID3char, "="); LegendEntryCompID3 = LegendEntryCompID3char2; cout<<"LegendEntryCompID3 = "<<LegendEntryCompID3<<endl;}
		if(std::string(argv[i]).find("LegendEntryCompID4") != std::string::npos) {char* LegendEntryCompID4char = argv[i]; char* LegendEntryCompID4char2 = strtok (LegendEntryCompID4char, "="); LegendEntryCompID4 = LegendEntryCompID4char2; cout<<"LegendEntryCompID4 = "<<LegendEntryCompID4<<endl;}

		if(std::string(argv[i]).find("ptBinMin") != std::string::npos) {char* ptBinMinchar = argv[i]; char* ptBinMinchar2 = strtok (ptBinMinchar, "p"); ptBinMin = atof(ptBinMinchar2); cout<<"ptBinMin = "<<ptBinMin<<endl;}
		if(std::string(argv[i]).find("ptBinMax") != std::string::npos) {char* ptBinMaxchar = argv[i]; char* ptBinMaxchar2 = strtok (ptBinMaxchar, "p"); ptBinMax = atof(ptBinMaxchar2); cout<<"ptBinMax = "<<ptBinMax<<endl;}
		if(std::string(argv[i]).find("nSystematics") != std::string::npos) {char* nSystematicschar = argv[i]; char* nSystematicschar2 = strtok (nSystematicschar, "p"); nSystematics = atof(nSystematicschar2); cout<<"nSystematics = "<<nSystematics<<endl;}
		if(std::string(argv[i]).find("nComp") != std::string::npos) {char* nCompchar = argv[i]; char* nCompchar2 = strtok (nCompchar, "p"); nComp = atof(nCompchar2); cout<<"nComp = "<<nComp<<endl;}
		if(std::string(argv[i]).find("nState") != std::string::npos) {char* nStatechar = argv[i]; char* nStatechar2 = strtok (nStatechar, "p"); nState = atof(nStatechar2); cout<<"nState = "<<nState<<endl;}

		if(std::string(argv[i]).find("PlotMatt=1") != std::string::npos) {PlotMatt=true; cout<<"Plot Matts results"<<endl;}
		if(std::string(argv[i]).find("PlotAsymm=1") != std::string::npos) {PlotAsymm=true; cout<<"Plot Asymms results"<<endl;}
		if(std::string(argv[i]).find("PlotCompare=1") != std::string::npos) {PlotCompare=true; cout<<"Plot Comparison results"<<endl;}
		if(std::string(argv[i]).find("PlotFinalData=1") != std::string::npos) {PlotFinalData=true; cout<<"Plot Data results"<<endl;}
		if(std::string(argv[i]).find("PlotSystematics=1") != std::string::npos) {PlotSystematics=true; cout<<"Plot Systematics results"<<endl;}
		if(std::string(argv[i]).find("PlotLegend=1") != std::string::npos) {PlotLegend=true; cout<<"Plot Legend"<<endl;}
		if(std::string(argv[i]).find("PlotBrazilian=1") != std::string::npos) {PlotBrazilian=true; cout<<"Plot Brazilian"<<endl;}
		if(std::string(argv[i]).find("FitGraph=1") != std::string::npos) {FitGraph=true; cout<<"FitGraph"<<endl;}
		if(std::string(argv[i]).find("DrawLatexStuff=1") != std::string::npos) {DrawLatexStuff=true; cout<<"DrawLatexStuff"<<endl;}
		if(std::string(argv[i]).find("MultiPanelPlots=1") != std::string::npos) {MultiPanelPlots=true; cout<<"MultiPanelPlots"<<endl;}
		if(std::string(argv[i]).find("PlotBG0plots=1") != std::string::npos) {PlotBG0plots=true; cout<<"PlotBG0plots"<<endl;}
		if(std::string(argv[i]).find("DeltaTildeplots=1") != std::string::npos) {DeltaTildeplots=true; cout<<"DeltaTildeplots"<<endl;}
		if(std::string(argv[i]).find("SBmSigPlots=1") != std::string::npos) {SBmSigPlots=true; cout<<"SBmSigPlots"<<endl;}
		if(std::string(argv[i]).find("PlotAlteredPPDResults=1") != std::string::npos) {PlotAlteredPPDResults=true; cout<<"PlotAlteredPPDResults"<<endl;}
		if(std::string(argv[i]).find("CompareSyst=1") != std::string::npos) {CompareSyst=true; cout<<"CompareSyst"<<endl;}
		if(std::string(argv[i]).find("SteerIndividuals=1") != std::string::npos) {SteerIndividuals=true; cout<<"SteerIndividuals"<<endl;}
		if(std::string(argv[i]).find("BGratioFits=1") != std::string::npos) {BGratioFits=true; cout<<"BGratioFits"<<endl;}
		if(std::string(argv[i]).find("BGratioChi2Fits=1") != std::string::npos) {BGratioChi2Fits=true; cout<<"BGratioChi2Fits"<<endl;}
		if(std::string(argv[i]).find("rapBinComb=1") != std::string::npos) {rapBinComb=true; cout<<"rapBinComb"<<endl;}
		if(std::string(argv[i]).find("SetCompStyle=1") != std::string::npos) {SetCompStyle=true; cout<<"SetCompStyle"<<endl;}
		if(std::string(argv[i]).find("DrawPreliminary=0") != std::string::npos) {DrawPreliminary=false; cout<<"DrawPreliminary"<<endl;}
		if(std::string(argv[i]).find("PlotMattForICHEP=1") != std::string::npos) {PlotMattForICHEP=true; cout<<"PlotMattForICHEP"<<endl;}
		if(std::string(argv[i]).find("ExtendLegendInX=1") != std::string::npos) {ExtendLegendInX=true; cout<<"ExtendLegendInX"<<endl;}
		if(std::string(argv[i]).find("ShiftInX=0") != std::string::npos) {ShiftInX=false; cout<<"ShiftInX false"<<endl;}
		if(std::string(argv[i]).find("ShiftCompareInX=1") != std::string::npos) {ShiftCompareInX=true; cout<<"ShiftCompareInX"<<endl;}
		if(std::string(argv[i]).find("PlotVsComp=1") != std::string::npos) {PlotVsComp=true; cout<<"PlotVsComp"<<endl;}
		if(std::string(argv[i]).find("PlotSysSquare=1") != std::string::npos) {PlotSysSquare=true; cout<<"PlotSysSquare"<<endl;}
		if(std::string(argv[i]).find("PlotCL1sigma=1") != std::string::npos) {PlotCL1sigma=true; cout<<"PlotCL1sigma"<<endl;}

	}

	gStyle->SetFrameBorderMode(0);
	double ColordBandWidth=1.;
	if(nState==5) ColordBandWidth = 1.2; //May2
	double DeltaXminOVERALL=0.;//0.9999;
	if(ShiftInX) DeltaXminOVERALL=0.9999;
	bool ShiftXminOVERALL=true;

	double DeltaXCompare=0.;
	if(ShiftCompareInX) DeltaXCompare=0.5;//0.9999

	double PlotpTMinInitial = 6., PlotpTMaxInitial = 72.;
	if(PlotFinalData) PlotpTMinInitial = 10.;
	double PlotpTMin = PlotpTMinInitial, 
				 PlotpTMax = PlotpTMaxInitial;

	int OneSigColor=416;
	int TwoSigColor=400;//858,898
	int ThreeSigColor=423;//423

	int CL1sigmaColor=kGreen;

	char filename[200];
	sprintf(filename,"%s/%s/TGraphResults_Psi%dS.root",storagedir,DefaultID,nState-3);
	if(CompareSyst)sprintf(filename,"%s/macros/polFit/Systematics/%s/%s/TGraphResults_Psi%dS.root",basedir,SystID1Base,SystID1Specify,nState-3);
	TFile *infileRes = new TFile(filename,"READ");

	sprintf(filename,"%s/%s/TGraphResults_Psi%dS_2sigma.root",storagedir,DefaultID,nState-3);
	TFile *infileRes2sigma = new TFile(filename,"READ");

	sprintf(filename,"%s/%s/TGraphResults_Psi%dS_3sigma.root",storagedir,DefaultID,nState-3);
	TFile *infileRes3sigma = new TFile(filename,"READ");

	if(!PlotBrazilian){
		infileRes2sigma=infileRes;
		infileRes3sigma=infileRes;
	}

	if(PlotMattForICHEP) PlotBrazilian=false;

	sprintf(filename,"/afs/ihep.ac.cn/users/z/zhangll/fs/work/polarization/PsiPol2011/macros/polFit/Systematics/TotalSyst/%s/TGraphResults_Psi%dS.root",MPCentralsWithTotalSystID,nState-3);
	TFile *infileStat = new TFile(filename,"READ");
	if(!PlotAlteredPPDResults) infileStat=infileRes;

	sprintf(filename,"%s/macros/polFit/Systematics/%s/%s/TGraphResults_Psi%dS.root",basedir,SystID1Base,SystID1Specify,nState-3);
	TFile *infileSyst1 = new TFile(filename,"READ");

	sprintf(filename,"%s/macros/polFit/Systematics/%s/%s/TGraphResults_Psi%dS.root",basedir,SystID2Base,SystID2Specify,nState-3);
	TFile *infileSyst2 = new TFile(filename,"READ");

	sprintf(filename,"%s/macros/polFit/Systematics/%s/%s/TGraphResults_Psi%dS.root",basedir,SystID3Base,SystID3Specify,nState-3);
	TFile *infileSyst3 = new TFile(filename,"READ");

	sprintf(filename,"%s/macros/polFit/Systematics/%s/%s/TGraphResults_Psi%dS.root",basedir,SystID4Base,SystID4Specify,nState-3);
	TFile *infileSyst4 = new TFile(filename,"READ");

	sprintf(filename,"%s/macros/polFit/Systematics/%s/%s/TGraphResults_Psi%dS.root",basedir,SystID5Base,SystID5Specify,nState-3);
	TFile *infileSyst5 = new TFile(filename,"READ");

	sprintf(filename,"%s/macros/polFit/Systematics/%s/%s/TGraphResults_Psi%dS.root",basedir,SystID6Base,SystID6Specify,nState-3);
	TFile *infileSyst6 = new TFile(filename,"READ");

	sprintf(filename,"%s/macros/polFit/Systematics/%s/%s/TGraphResults_Psi%dS.root",basedir,SystID7Base,SystID7Specify,nState-3);
	TFile *infileSyst7 = new TFile(filename,"READ");

	sprintf(filename,"%s/macros/polFit/Systematics/%s/%s/TGraphResults_Psi%dS.root",basedir,SystID8Base,SystID8Specify,nState-3);
	TFile *infileSyst8 = new TFile(filename,"READ");

	//sprintf(filename,"%s/macros/polFit/MattRes/TGraphResults_Psi%dS_MattResults.root",basedir,nState-3);
	//TFile *MattFileStat = new TFile(filename,"READ");
	//if(nState>1) MattFileStat=infileRes;

	//sprintf(filename,"%s/macros/polFit/MattRes/TGraphResults_Psi%dS_MattSyst.root",basedir,nState-3);
	//TFile *MattFileSyst = new TFile(filename,"READ");
	//if(nState>1) MattFileSyst=infileRes;

	sprintf(filename,"%s/macros/polFit/CDFRes/TGraphResults_Psi%dS_CDFStat.root",basedir,nState-3);
	TFile *MattFileStat = new TFile(filename,"READ");

	sprintf(filename,"%s/macros/polFit/CDFRes/TGraphResults_Psi%dS_CDFSyst.root",basedir,nState-3);
	TFile *MattFileSyst = new TFile(filename,"READ");

	sprintf(filename,"%s/macros/polFit/CDFRes/TGraphResults_Psi%dS_CDFTotal.root",basedir,nState-3);
	TFile *MattFileTotal = new TFile(filename,"READ");

	sprintf(filename,"%s/macros/polFit/CDFRes/TGraphResults_Psi1S_CDFStat.root",basedir);
	TFile *infileMP1SCDF_Stat = new TFile(filename,"READ");
	sprintf(filename,"%s/macros/polFit/CDFRes/TGraphResults_Psi2S_CDFStat.root",basedir);
	TFile *infileMP2SCDF_Stat = new TFile(filename,"READ");
	sprintf(filename,"%s/macros/polFit/CDFRes/TGraphResults_Psi3S_CDFStat.root",basedir);
	TFile *infileMP3SCDF_Stat = new TFile(filename,"READ");

	sprintf(filename,"%s/macros/polFit/CDFRes/TGraphResults_Psi1S_CDFTotal.root",basedir);
	TFile *infileMP1SCDF_Total = new TFile(filename,"READ");
	sprintf(filename,"%s/macros/polFit/CDFRes/TGraphResults_Psi2S_CDFTotal.root",basedir);
	TFile *infileMP2SCDF_Total = new TFile(filename,"READ");
	sprintf(filename,"%s/macros/polFit/CDFRes/TGraphResults_Psi3S_CDFTotal.root",basedir);
	TFile *infileMP3SCDF_Total = new TFile(filename,"READ");

	sprintf(filename,"%s/%s/TGraphResults_Psi%dS.root",storagedir,CompareID1,nState-3);
	if(CompareSyst)sprintf(filename,"%s/macros/polFit/Systematics/%s/%s/TGraphResults_Psi1S.root",basedir,SystID2Base,SystID2Specify,nState);
	TFile *CompareFile1 = new TFile(filename,"READ");
	if(nComp<1) CompareFile1=infileRes;

	sprintf(filename,"%s/%s/TGraphResults_Psi%dS.root",storagedir,CompareID2,nState-3);
	if(CompareSyst)sprintf(filename,"%s/macros/polFit/Systematics/%s/%s/TGraphResults_Psi2S.root",basedir,SystID3Base,SystID3Specify,nState);
	TFile *CompareFile2 = new TFile(filename,"READ");
	if(nComp<2) CompareFile2=infileRes;

	sprintf(filename,"%s/%s/TGraphResults_Psi%dS.root",storagedir,CompareID3,nState-3);
	if(CompareSyst)sprintf(filename,"%s/macros/polFit/Systematics/%s/%s/TGraphResults_Psi2S.root",basedir,SystID4Base,SystID4Specify,nState);
	TFile *CompareFile3 = new TFile(filename,"READ");
	if(nComp<3) CompareFile3=infileRes;

	sprintf(filename,"%s/%s/TGraphResults_Psi%dS.root",storagedir,CompareID4,nState-3);
	TFile *CompareFile4 = new TFile(filename,"READ");
	if(nComp<4) CompareFile4=infileRes;




	sprintf(filename,"/afs/ihep.ac.cn/users/z/zhangll/fs/work/polarization/PsiPol2011/macros/polFit/Systematics/TotalSyst/%s/TGraphResults_Psi1S.root",MPCentralsWithTotalSystID);
	TFile *infileMP1 = new TFile(filename,"READ");
	if(!MultiPanelPlots) infileMP1=infileRes;

	sprintf(filename,"/afs/ihep.ac.cn/users/z/zhangll/fs/work/polarization/PsiPol2011/macros/polFit/Systematics/TotalSyst/%s/TGraphResults_Psi2S.root",MPCentralsWithTotalSystID);
	TFile *infileMP2 = new TFile(filename,"READ");
	if(!MultiPanelPlots) infileMP2=infileRes;

	sprintf(filename,"/afs/ihep.ac.cn/users/z/zhangll/fs/work/polarization/PsiPol2011/macros/polFit/Systematics/TotalSyst/%s/TGraphResults_Psi3S.root",MPCentralsWithTotalSystID);
	TFile *infileMP3 = new TFile(filename,"READ");
	if(!MultiPanelPlots) infileMP3=infileRes;


	sprintf(filename,"%s/%s/TGraphResults_Psi1S_1sigma.root",storagedir,DefaultID);
	TFile *infileMP1_1sig = new TFile(filename,"READ");
	if(!MultiPanelPlots) infileMP1_1sig=infileRes;
	sprintf(filename,"%s/%s/TGraphResults_Psi1S_2sigma.root",storagedir,DefaultID);
	TFile *infileMP1_2sig = new TFile(filename,"READ");
	if(!MultiPanelPlots) infileMP1_2sig=infileRes;
	sprintf(filename,"%s/%s/TGraphResults_Psi1S_3sigma.root",storagedir,DefaultID);
	TFile *infileMP1_3sig = new TFile(filename,"READ");
	if(!MultiPanelPlots) infileMP1_3sig=infileRes;

	sprintf(filename,"%s/%s/TGraphResults_Psi2S_1sigma.root",storagedir,DefaultID);
	TFile *infileMP2_1sig = new TFile(filename,"READ");
	if(!MultiPanelPlots) infileMP2_1sig=infileRes;
	sprintf(filename,"%s/%s/TGraphResults_Psi2S_2sigma.root",storagedir,DefaultID);
	TFile *infileMP2_2sig = new TFile(filename,"READ");
	if(!MultiPanelPlots) infileMP2_2sig=infileRes;
	sprintf(filename,"%s/%s/TGraphResults_Psi2S_3sigma.root",storagedir,DefaultID);
	TFile *infileMP2_3sig = new TFile(filename,"READ");
	if(!MultiPanelPlots) infileMP2_3sig=infileRes;

	sprintf(filename,"%s/%s/TGraphResults_Psi3S_1sigma.root",storagedir,DefaultID);
	TFile *infileMP3_1sig = new TFile(filename,"READ");
	if(!MultiPanelPlots) infileMP3_1sig=infileRes;
	sprintf(filename,"%s/%s/TGraphResults_Psi3S_2sigma.root",storagedir,DefaultID);
	TFile *infileMP3_2sig = new TFile(filename,"READ");
	if(!MultiPanelPlots) infileMP3_2sig=infileRes;
	sprintf(filename,"%s/%s/TGraphResults_Psi3S_3sigma.root",storagedir,DefaultID);
	TFile *infileMP3_3sig = new TFile(filename,"READ");
	if(!MultiPanelPlots) infileMP3_3sig=infileRes;

	if(!PlotCompare) CompareFile1=infileRes;
	if(!PlotCompare||nComp<2) CompareFile2=infileRes;
	if(!PlotCompare||nComp<3) CompareFile3=infileRes;
	if(!PlotCompare||nComp<4) CompareFile4=infileRes;
	if(!PlotBrazilian) infileRes2sigma=infileRes;
	if(!PlotBrazilian) infileRes3sigma=infileRes;

	if(SBmSigPlots){
		CompareFile1=infileSyst1;
		CompareFile2=infileSyst2;
	}
	if(BGratioFits){
		sprintf(filename,"%s/macros/polFit/Systematics/%s/%s/TGraphResults_Psi1S.root",basedir,SystID1Base,SystID1Specify);
		CompareFile1 = new TFile(filename,"READ");
		sprintf(filename,"%s/macros/polFit/Systematics/%s/%s/TGraphResults_Psi2S.root",basedir,SystID2Base,SystID2Specify);
		CompareFile2 = new TFile(filename,"READ");
		sprintf(filename,"%s/macros/polFit/Systematics/%s/%s/TGraphResults_Psi3S.root",basedir,SystID3Base,SystID3Specify);
		CompareFile3 = new TFile(filename,"READ");
	}



	char JobDir[200],StateDir[200],FigDir[200],SuppDir[200];
	sprintf(JobDir,"%s/macros/polFit/FinalResults/%s",basedir,JobID);
	sprintf(StateDir,"%s/Psi%dS",JobDir,nState-3);
	sprintf(FigDir,"%s/Figures",StateDir);
	sprintf(SuppDir,"%s",JobDir);
	gSystem->mkdir(JobDir);
	gSystem->mkdir(StateDir);
	gSystem->mkdir(FigDir);

	char GraphName[200],axislabel[200], GraphNameOtherRap[200];
	double yMin,yMax;


	// Declare variables needed for table production

	const int tabPtBins=ptBinMax-ptBinMin+1;
	const int iParameters=18;

	int  nRapBins =2;
	if(nState==5) nRapBins=3;
	cout << "nRapBins: " << nRapBins << endl;

	double val_table[iParameters+1][nRapBins][tabPtBins];
	double errHigh_table[iParameters+1][nRapBins][tabPtBins];
	double errLow_table[iParameters+1][nRapBins][tabPtBins];
	double syst_table[iParameters+1][nRapBins][tabPtBins][nSystematics+1];
	double errHighTotal1_table[iParameters+1][nRapBins][tabPtBins];
	double errLowTotal1_table[iParameters+1][nRapBins][tabPtBins];
	double errHighTotal2_table[iParameters+1][nRapBins][tabPtBins];
	double errLowTotal2_table[iParameters+1][nRapBins][tabPtBins];
	double errHighTotal3_table[iParameters+1][nRapBins][tabPtBins];
	double errLowTotal3_table[iParameters+1][nRapBins][tabPtBins];
	double pTmean_table[iParameters+1][nRapBins][tabPtBins];


	TCanvas *MPcanvasCS;
	TCanvas *MPcanvasHX;
	TCanvas *MPcanvasPX;
	TCanvas *MPcanvasCS_New;
	TCanvas *MPcanvasHX_New;
	TCanvas *MPcanvasPX_New;
	TCanvas *MPcanvasCS_Psi;
	TCanvas *MPcanvasHX_Psi;
	TCanvas *MPcanvasPX_Psi;
	TCanvas *MPcanvasCS_rap1;
	TCanvas *MPcanvasHX_rap1;
	TCanvas *MPcanvasPX_rap1;
	TCanvas *MPcanvasCS_rap2;
	TCanvas *MPcanvasHX_rap2;
	TCanvas *MPcanvasPX_rap2;
	TCanvas *MPcanvasCS_rap3;
	TCanvas *MPcanvasHX_rap3;
	TCanvas *MPcanvasPX_rap3;
	TCanvas *MPcanvasTilde;
	TCanvas *MPcanvasTilde_rap1;
	TCanvas *MPcanvasTilde_rap2;
	TCanvas *MPcanvasTilde_New;
	TCanvas *MPcanvasTilde_Psi;
	TCanvas *MPcanvasCDF;

	//================================================================
	//==need to be changed for different rap bins(Psi1S:4, Psi2S:5)===
	//================================================================

	for(int iLam = 1; iLam<iParameters+1; iLam++){

		for(int rapBin = 1; rapBin < nRapBins+1; rapBin++){


			if(iLam==1)  sprintf(GraphName,"lth_CS_rap%d",rapBin);
			if(iLam==2)  sprintf(GraphName,"lph_CS_rap%d",rapBin);
			if(iLam==3)  sprintf(GraphName,"ltp_CS_rap%d",rapBin);
			if(iLam==4)  sprintf(GraphName,"lthstar_CS_rap%d",rapBin);
			if(iLam==5)  sprintf(GraphName,"lphstar_CS_rap%d",rapBin);
			if(iLam==6)  sprintf(GraphName,"ltilde_CS_rap%d",rapBin);

			if(iLam==7)  sprintf(GraphName,"lth_HX_rap%d",rapBin);
			if(iLam==8)  sprintf(GraphName,"lph_HX_rap%d",rapBin);
			if(iLam==9)  sprintf(GraphName,"ltp_HX_rap%d",rapBin);
			if(iLam==10) sprintf(GraphName,"lthstar_HX_rap%d",rapBin);
			if(iLam==11) sprintf(GraphName,"lphstar_HX_rap%d",rapBin);
			if(iLam==12) sprintf(GraphName,"ltilde_HX_rap%d",rapBin);

			if(iLam==13) sprintf(GraphName,"lth_PX_rap%d",rapBin);
			if(iLam==14) sprintf(GraphName,"lph_PX_rap%d",rapBin);
			if(iLam==15) sprintf(GraphName,"ltp_PX_rap%d",rapBin);
			if(iLam==16) sprintf(GraphName,"lthstar_PX_rap%d",rapBin);
			if(iLam==17) sprintf(GraphName,"lphstar_PX_rap%d",rapBin);
			if(iLam==18) sprintf(GraphName,"ltilde_PX_rap%d",rapBin);

			if(rapBin==1){
				if(iLam==1)  sprintf(GraphNameOtherRap,"lth_CS_rap2");
				if(iLam==2)  sprintf(GraphNameOtherRap,"lph_CS_rap2");
				if(iLam==3)  sprintf(GraphNameOtherRap,"ltp_CS_rap2");
				if(iLam==4)  sprintf(GraphNameOtherRap,"lthstar_CS_rap2");
				if(iLam==5)  sprintf(GraphNameOtherRap,"lphstar_CS_rap2");
				if(iLam==6)  sprintf(GraphNameOtherRap,"ltilde_CS_rap2");

				if(iLam==7)  sprintf(GraphNameOtherRap,"lth_HX_rap2");
				if(iLam==8)  sprintf(GraphNameOtherRap,"lph_HX_rap2");
				if(iLam==9)  sprintf(GraphNameOtherRap,"ltp_HX_rap2");
				if(iLam==10) sprintf(GraphNameOtherRap,"lthstar_HX_rap2");
				if(iLam==11) sprintf(GraphNameOtherRap,"lphstar_HX_rap2");
				if(iLam==12) sprintf(GraphNameOtherRap,"ltilde_HX_rap2");

				if(iLam==13) sprintf(GraphNameOtherRap,"lth_PX_rap2");
				if(iLam==14) sprintf(GraphNameOtherRap,"lph_PX_rap2");
				if(iLam==15) sprintf(GraphNameOtherRap,"ltp_PX_rap2");
				if(iLam==16) sprintf(GraphNameOtherRap,"lthstar_PX_rap2");
				if(iLam==17) sprintf(GraphNameOtherRap,"lphstar_PX_rap2");
				if(iLam==18) sprintf(GraphNameOtherRap,"ltilde_PX_rap2");
			}

			if(rapBin==2){
				if(iLam==1)  sprintf(GraphNameOtherRap,"lth_CS_rap1");
				if(iLam==2)  sprintf(GraphNameOtherRap,"lph_CS_rap1");
				if(iLam==3)  sprintf(GraphNameOtherRap,"ltp_CS_rap1");
				if(iLam==4)  sprintf(GraphNameOtherRap,"lthstar_CS_rap1");
				if(iLam==5)  sprintf(GraphNameOtherRap,"lphstar_CS_rap1");
				if(iLam==6)  sprintf(GraphNameOtherRap,"ltilde_CS_rap1");

				if(iLam==7)  sprintf(GraphNameOtherRap,"lth_HX_rap1");
				if(iLam==8)  sprintf(GraphNameOtherRap,"lph_HX_rap1");
				if(iLam==9)  sprintf(GraphNameOtherRap,"ltp_HX_rap1");
				if(iLam==10) sprintf(GraphNameOtherRap,"lthstar_HX_rap1");
				if(iLam==11) sprintf(GraphNameOtherRap,"lphstar_HX_rap1");
				if(iLam==12) sprintf(GraphNameOtherRap,"ltilde_HX_rap1");

				if(iLam==13) sprintf(GraphNameOtherRap,"lth_PX_rap1");
				if(iLam==14) sprintf(GraphNameOtherRap,"lph_PX_rap1");
				if(iLam==15) sprintf(GraphNameOtherRap,"ltp_PX_rap1");
				if(iLam==16) sprintf(GraphNameOtherRap,"lthstar_PX_rap1");
				if(iLam==17) sprintf(GraphNameOtherRap,"lphstar_PX_rap1");
				if(iLam==18) sprintf(GraphNameOtherRap,"ltilde_PX_rap1");
			}



			char beginLamLabel[200];
			if(PlotSystematics) sprintf(beginLamLabel,"#Delta");
			if(PlotFinalData) sprintf(beginLamLabel,"");
			if(PlotSystematics&&!PlotAsymm) sprintf(beginLamLabel,"#sigma");
			if(PlotSystematics&&PlotSysSquare) sprintf(beginLamLabel,"#sigma^{2}");
			char endLamLabel[200];
			sprintf(endLamLabel,"");


			if(iLam==1)  sprintf(axislabel,"%s#lambda^{CS}_{#vartheta}%s",beginLamLabel,endLamLabel);
			if(iLam==2)  sprintf(axislabel,"%s#lambda^{CS}_{#varphi}%s",beginLamLabel,endLamLabel);
			if(iLam==3)  sprintf(axislabel,"%s#lambda^{CS}_{#vartheta#varphi}%s",beginLamLabel,endLamLabel);
			if(iLam==4)  sprintf(axislabel,"%s#lambda^{*CS}_{#vartheta}%s",beginLamLabel,endLamLabel);
			if(iLam==5)  sprintf(axislabel,"%s#lambda^{*CS}_{#varphi}%s",beginLamLabel,endLamLabel);
			if(iLam==6)  sprintf(axislabel,"%s#tilde{#lambda}^{CS}%s",beginLamLabel,endLamLabel);

			if(iLam==7)  sprintf(axislabel,"%s#lambda^{HX}_{#vartheta}%s",beginLamLabel,endLamLabel);
			if(iLam==8)  sprintf(axislabel,"%s#lambda^{HX}_{#varphi}%s",beginLamLabel,endLamLabel);
			if(iLam==9)  sprintf(axislabel,"%s#lambda^{HX}_{#vartheta#varphi}%s",beginLamLabel,endLamLabel);
			if(iLam==10) sprintf(axislabel,"%s#lambda^{*HX}_{#vartheta}%s",beginLamLabel,endLamLabel);
			if(iLam==11) sprintf(axislabel,"%s#lambda^{*HX}_{#varphi}%s",beginLamLabel,endLamLabel);
			if(iLam==12) sprintf(axislabel,"%s#tilde{#lambda}^{HX}%s",beginLamLabel,endLamLabel);

			if(iLam==13) sprintf(axislabel,"%s#lambda^{PX}_{#vartheta}%s",beginLamLabel,endLamLabel);
			if(iLam==14) sprintf(axislabel,"%s#lambda^{PX}_{#varphi}%s",beginLamLabel,endLamLabel);
			if(iLam==15) sprintf(axislabel,"%s#lambda^{PX}_{#vartheta#varphi}%s",beginLamLabel,endLamLabel);
			if(iLam==16) sprintf(axislabel,"%s#lambda^{*PX}_{#vartheta}%s",beginLamLabel,endLamLabel);
			if(iLam==17) sprintf(axislabel,"%s#lambda^{*PX}_{#varphi}%s",beginLamLabel,endLamLabel);
			if(iLam==18) sprintf(axislabel,"%s#tilde{#lambda}^{PX}%s",beginLamLabel,endLamLabel);//IfLamTildeClosure no PX
			if(iLam==18&&DeltaTildeplots) sprintf(axislabel,"%s#tilde{#lambda}%s",beginLamLabel,endLamLabel);

			if(PlotSystematics&&PlotSysSquare){
				if(iLam==1)  sprintf(axislabel,"%s (#lambda^{CS}_{#vartheta}) %s",beginLamLabel,endLamLabel);
				if(iLam==2)  sprintf(axislabel,"%s (#lambda^{CS}_{#varphi}) %s",beginLamLabel,endLamLabel);
				if(iLam==3)  sprintf(axislabel,"%s (#lambda^{CS}_{#vartheta#varphi}) %s",beginLamLabel,endLamLabel);
				if(iLam==4)  sprintf(axislabel,"%s (#lambda^{*CS}_{#vartheta}) %s",beginLamLabel,endLamLabel);
				if(iLam==5)  sprintf(axislabel,"%s (#lambda^{*CS}_{#varphi}) %s",beginLamLabel,endLamLabel);
				if(iLam==6)  sprintf(axislabel,"%s (#tilde{#lambda}^{CS}) %s",beginLamLabel,endLamLabel);

				if(iLam==7)  sprintf(axislabel,"%s (#lambda^{HX}_{#vartheta}) %s",beginLamLabel,endLamLabel);
				if(iLam==8)  sprintf(axislabel,"%s (#lambda^{HX}_{#varphi}) %s",beginLamLabel,endLamLabel);
				if(iLam==9)  sprintf(axislabel,"%s (#lambda^{HX}_{#vartheta#varphi}) %s",beginLamLabel,endLamLabel);
				if(iLam==10) sprintf(axislabel,"%s (#lambda^{*HX}_{#vartheta}) %s",beginLamLabel,endLamLabel);
				if(iLam==11) sprintf(axislabel,"%s (#lambda^{*HX}_{#varphi}) %s",beginLamLabel,endLamLabel);
				if(iLam==12) sprintf(axislabel,"%s (#tilde{#lambda}^{HX}) %s",beginLamLabel,endLamLabel);

				if(iLam==13) sprintf(axislabel,"%s (#lambda^{PX}_{#vartheta}) %s",beginLamLabel,endLamLabel);
				if(iLam==14) sprintf(axislabel,"%s (#lambda^{PX}_{#varphi}) %s",beginLamLabel,endLamLabel);
				if(iLam==15) sprintf(axislabel,"%s (#lambda^{PX}_{#vartheta#varphi}) %s",beginLamLabel,endLamLabel);
				if(iLam==16) sprintf(axislabel,"%s (#lambda^{*PX}_{#vartheta}) %s",beginLamLabel,endLamLabel);
				if(iLam==17) sprintf(axislabel,"%s (#lambda^{*PX}_{#varphi}) %s",beginLamLabel,endLamLabel);
				if(iLam==18) sprintf(axislabel,"%s (#tilde{#lambda}^{PX}) %s",beginLamLabel,endLamLabel);//IfLamTildeClosure no PX
				if(iLam==18&&DeltaTildeplots) sprintf(axislabel,"%s (#tilde{#lambda}) %s",beginLamLabel,endLamLabel);
			}


			if(iLam==1)  sprintf(filename,"%s/FinalResults_CS_lth_rap%d.pdf",FigDir,rapBin);
			if(iLam==2)  sprintf(filename,"%s/FinalResults_CS_lph_rap%d.pdf",FigDir,rapBin);
			if(iLam==3)  sprintf(filename,"%s/FinalResults_CS_ltp_rap%d.pdf",FigDir,rapBin);
			if(iLam==4)  sprintf(filename,"%s/FinalResults_CS_lthstar_rap%d.pdf",FigDir,rapBin);
			if(iLam==5)  sprintf(filename,"%s/FinalResults_CS_lphstar_rap%d.pdf",FigDir,rapBin);
			if(iLam==6)  sprintf(filename,"%s/FinalResults_CS_ltilde_rap%d.pdf",FigDir,rapBin);

			if(iLam==7)  sprintf(filename,"%s/FinalResults_HX_lth_rap%d.pdf",FigDir,rapBin);
			if(iLam==8)  sprintf(filename,"%s/FinalResults_HX_lph_rap%d.pdf",FigDir,rapBin);
			if(iLam==9)  sprintf(filename,"%s/FinalResults_HX_ltp_rap%d.pdf",FigDir,rapBin);
			if(iLam==10) sprintf(filename,"%s/FinalResults_HX_lthstar_rap%d.pdf",FigDir,rapBin);
			if(iLam==11) sprintf(filename,"%s/FinalResults_HX_lphstar_rap%d.pdf",FigDir,rapBin);
			if(iLam==12) sprintf(filename,"%s/FinalResults_HX_ltilde_rap%d.pdf",FigDir,rapBin);

			if(iLam==13) sprintf(filename,"%s/FinalResults_PX_lth_rap%d.pdf",FigDir,rapBin);
			if(iLam==14) sprintf(filename,"%s/FinalResults_PX_lph_rap%d.pdf",FigDir,rapBin);
			if(iLam==15) sprintf(filename,"%s/FinalResults_PX_ltp_rap%d.pdf",FigDir,rapBin);
			if(iLam==16) sprintf(filename,"%s/FinalResults_PX_lthstar_rap%d.pdf",FigDir,rapBin);
			if(iLam==17) sprintf(filename,"%s/FinalResults_PX_lphstar_rap%d.pdf",FigDir,rapBin);
			if(iLam==18) sprintf(filename,"%s/FinalResults_PX_ltilde_rap%d.pdf",FigDir,rapBin);

			//if(iLam==18) sprintf(axislabel,"#tilde{#lambda}");

			if(BGratioChi2Fits) sprintf(axislabel,"#chi^{2} / ndf");

			int nFrame=0;

			if(iLam>0&&iLam<7) nFrame=1;
			if(iLam>6&&iLam<13) nFrame=2;
			if(iLam>12&&iLam<19) nFrame=3;

			yMin=-1.09;
			yMax=1.09;

			//yMin=-2.09;
			//yMax=2.09;

			//yMin=-0.1;//-0.209;
			//yMax=0.209;

			//yMin=-1.49; //plot noRho
			//yMax=1.49; //plot noRho

			if(iLam==2||iLam==8||iLam==14||iLam==3||iLam==9||iLam==15){
				yMin=-0.35;
				yMax=0.35;
				//if(iLam==2||iLam==8||iLam==14){
				//	yMin=-0.6;
				//	yMax=0.6;
				//}

				//yMin=-0.85;
				//yMax=0.85;

				//yMin=-0.05;//-0.109;
				//yMax=0.109;

				//yMin=-0.55; //plot noRho
				//yMax=0.55; //plot noRho
			}

			if(iLam==6||iLam==12||iLam==18){
				yMin=-1.09;
				yMax=1.09;

				//yMin=-2.09;
				//yMax=2.09;

				//yMin=-0.1;//-0.219;
				//yMax=0.219;

			}

			if(PlotBrazilian){
				yMin=-1.5;
				yMax=1.5;
				//yMin=-1.275;
				//yMax=1.275;

				if(iLam==2||iLam==8||iLam==14||iLam==3||iLam==9||iLam==15){
					//yMin=-0.55;
					//yMax=0.55;
					yMin=-0.6;
					yMax=0.6;
				}

				if(iLam==6||iLam==12||iLam==18){
					yMin=-1.5;
					yMax=2.1;
					//yMin=-1.;
					//yMax=1.75;
				}
			}
			//yMin=-5.5; //ifPull
			//yMax=5.5; //ifPull

			// Background Pol Fits:

			if(PlotBG0plots){
				yMin=-2;
				yMax=5;

				if(iLam==2||iLam==8||iLam==14||iLam==3||iLam==9||iLam==15){
					yMin=-1.5;
					yMax=1.5;
				}

				if(iLam==6||iLam==12||iLam==18){
					yMin=-7.5;
					yMax=20;
				}
			}

			if(SBmSigPlots){
				yMin=-2;
				yMax=2;

				if(iLam==2||iLam==8||iLam==14||iLam==3||iLam==9||iLam==15){
					yMin=-0.6;
					yMax=0.6;
				}

				if(iLam==6||iLam==12||iLam==18){
					yMin=-6;
					yMax=6;
				}

			}

			if(BGratioChi2Fits){
				yMin=0.1;
				yMax=30;

			}
			/*yMin=-4;
				yMax=4;

				if(iLam==2||iLam==8||iLam==14||iLam==3||iLam==9||iLam==15){
				yMin=-4;
				yMax=4;
				}

				if(iLam==6||iLam==12||iLam==18){
				yMin=-4;
				yMax=4;
				} */

			//if(CompareSyst||BGratioFits)
			if(BGratioFits){
				yMin=-2.5;
				yMax=2.5;

				if(iLam==2||iLam==8||iLam==14||iLam==3||iLam==9||iLam==15){
					yMin=-2.5;
					yMax=2.5;
				}

				if(iLam==6||iLam==12||iLam==18){
					yMin=-2.5;
					yMax=2.5;
				}
			}

			if(PlotMatt&&PlotMattForICHEP){
				yMin=-1.85;
				yMax=1.85;
				if(nState==1) yMin=-1.1;
				if(nState==1) yMax=1.1;

				if(iLam==2||iLam==8||iLam==14||iLam==3||iLam==9||iLam==15){
					yMin=-0.5;
					yMax=0.5;
				}

				if(iLam==6||iLam==12||iLam==18){
					yMin=-1.3;
					yMax=1.3;
				}
			}

			if(PlotSystematics&&PlotSysSquare){
				yMin=-0.06; //-0.15;
				yMax=0.06; //0.15;
				if(nState==5){
					yMin=-0.085; //-0.1;
					yMax=0.085; //0.1;
					if(iLam==6||iLam==12||iLam==18){
						yMin=-0.2;
						yMax=0.2;
					}
				}

				if(iLam==2||iLam==8||iLam==14||iLam==3||iLam==9||iLam==15){
					yMin=-0.02; //-0.025
					yMax=0.02;  //0.025
					if(iLam==2||iLam==8||iLam==14){
						yMin=-0.012;
						yMax=0.012; 
					}
					if(nState==4){
						yMin=-0.004; //-0.005;
						yMax=0.004; //0.005;
					}
				} 

				if(iLam==1||iLam==2||iLam==3){
					yMin=-0.02; //-0.05;
					yMax=0.02; //0.05;
					if(nState==4&&iLam==3){
						yMin=-0.005;
						yMax=0.005; 
					}
					if(nState==5){
						yMin=-0.035; //-0.05;
						yMax=0.035; //0.05;
						if(iLam==3){
							yMin=-0.02;
							yMax=0.02; 
						}
					}
				}

				//if(iLam==6||iLam==12||iLam==18){
				//	yMin=-0.25;
				//	yMax=0.25;
				//}

			}

			TGraphAsymmErrors* graphDefaultRes = (TGraphAsymmErrors*) infileRes->Get(GraphName);
			TGraphAsymmErrors* graphDefaultRes2sigma = (TGraphAsymmErrors*) infileRes2sigma->Get(GraphName);
			TGraphAsymmErrors* graphDefaultRes3sigma = (TGraphAsymmErrors*) infileRes3sigma->Get(GraphName);
			TGraphAsymmErrors* graphDefaultStat = (TGraphAsymmErrors*) infileStat->Get(GraphName);
			TGraphAsymmErrors* graphSyst1 = (TGraphAsymmErrors*) infileSyst1->Get(GraphName);
			TGraphAsymmErrors* graphSyst2 = (TGraphAsymmErrors*) infileSyst2->Get(GraphName);
			TGraphAsymmErrors* graphSyst3 = (TGraphAsymmErrors*) infileSyst3->Get(GraphName);
			TGraphAsymmErrors* graphSyst4 = (TGraphAsymmErrors*) infileSyst4->Get(GraphName);
			TGraphAsymmErrors* graphSyst5 = (TGraphAsymmErrors*) infileSyst5->Get(GraphName);
			TGraphAsymmErrors* graphSyst6 = (TGraphAsymmErrors*) infileSyst6->Get(GraphName);
			TGraphAsymmErrors* graphSyst7 = (TGraphAsymmErrors*) infileSyst7->Get(GraphName);
			TGraphAsymmErrors* graphSyst8 = (TGraphAsymmErrors*) infileSyst8->Get(GraphName);

			TGraphAsymmErrors* graphMattStat = (TGraphAsymmErrors*) MattFileStat->Get(GraphName);
			TGraphAsymmErrors* graphMattSyst = (TGraphAsymmErrors*) MattFileSyst->Get(GraphName);
			TGraphAsymmErrors* graphMattTotal = (TGraphAsymmErrors*) MattFileTotal->Get(GraphName);
			TGraphAsymmErrors* graphCompareFile1 = (TGraphAsymmErrors*) CompareFile1->Get(GraphName);
			TGraphAsymmErrors* graphCompareFile2 = (TGraphAsymmErrors*) CompareFile2->Get(GraphName);
			TGraphAsymmErrors* graphCompareFile3 = (TGraphAsymmErrors*) CompareFile3->Get(GraphName);
			TGraphAsymmErrors* graphCompareFile4 = (TGraphAsymmErrors*) CompareFile4->Get(GraphName);

			TGraphAsymmErrors* graphNLONRQCD = (TGraphAsymmErrors*) MattFileTotal->Get(GraphName);
			TGraphAsymmErrors* graphNNLO = (TGraphAsymmErrors*) MattFileTotal->Get(GraphName);


			int nBinspT=ptBinMax-ptBinMin+1;
			double ptCentre_[nBinspT];
			double ptCentreErr_low[nBinspT];
			double ptCentreErr_high[nBinspT];
			double lmean[nBinspT];
			double lmean_errlow[nBinspT];
			double lmean_errhigh[nBinspT];
			double lmean_errmean[nBinspT];
			double lmean_errmean_minus[nBinspT];
			double lmeanTotal1_errlow[nBinspT];
			double lmeanTotal1_errhigh[nBinspT];
			double lmeanTotal2_errlow[nBinspT];
			double lmeanTotal2_errhigh[nBinspT];
			double lmeanTotal3_errlow[nBinspT];
			double lmeanTotal3_errhigh[nBinspT];
			double ptCentre_ForTable[nBinspT];
			double lmeanBuff[nBinspT];

			double fit_lmean_errmean[nBinspT];
			double fit_lmean_errlow[nBinspT];
			double fit_lmean_errhigh[nBinspT];

			double SystError1[nBinspT];
			double SystError2[nBinspT];
			double SystError3[nBinspT];
			double SystError4[nBinspT];
			double SystError5[nBinspT];
			double SystError6[nBinspT];
			double SystError7[nBinspT];
			double SystError8[nBinspT];
			double SystError[nBinspT];
			double ErrSystError1[nBinspT];
			double ErrSystError2[nBinspT];
			double ErrSystError3[nBinspT];
			double ErrSystError4[nBinspT];
			double ErrSystError5[nBinspT];
			double ErrSystError6[nBinspT];
			double ErrSystError7[nBinspT];
			double ErrSystError8[nBinspT];

			double SystError12[nBinspT];
			double SystError123[nBinspT];
			double SystError1234[nBinspT];
			double SystError12345[nBinspT];
			double SystError123456[nBinspT];
			double SystError1234567[nBinspT];
			double SystError12345678[nBinspT];

			double Buffer[nBinspT];


			int pt=0;
			for(int ptBin = ptBinMin; ptBin < ptBinMax+1; ptBin++) {

				graphDefaultRes->GetPoint(ptBin-1,ptCentre_[pt],lmean[pt]);
				cout<<"debug: ptCentre["<<pt<<"]: "<<ptCentre_[pt]<<endl;
				cout<<"debug: lmean["<<pt<<"]: "<<lmean[pt]<<endl;
				ptCentreErr_high[pt]=graphDefaultRes->GetErrorXhigh(ptBin-1);
				ptCentreErr_low[pt]=graphDefaultRes->GetErrorXlow(ptBin-1);
				lmean_errhigh[pt]=graphDefaultRes->GetErrorYhigh(ptBin-1);
				lmean_errlow[pt]=graphDefaultRes->GetErrorYlow(ptBin-1);
				lmean_errmean[pt]=(lmean_errlow[pt]+lmean_errhigh[pt])/2.;
				lmean_errmean_minus[pt]=-(lmean_errlow[pt]+lmean_errhigh[pt])/2.;

				//cout<<"lmean_errlow: "<<lmean_errlow[pt]<<endl;
				//cout<<"lmean_errhigh: "<<lmean_errhigh[pt]<<endl;

				lmeanTotal1_errlow[pt]=lmean_errlow[pt];
				lmeanTotal1_errhigh[pt]=lmean_errhigh[pt];
				lmeanTotal2_errlow[pt]=graphDefaultRes2sigma->GetErrorYlow(ptBin-1);
				lmeanTotal2_errhigh[pt]=graphDefaultRes2sigma->GetErrorYhigh(ptBin-1);
				lmeanTotal3_errlow[pt]=graphDefaultRes3sigma->GetErrorYlow(ptBin-1);
				lmeanTotal3_errhigh[pt]=graphDefaultRes3sigma->GetErrorYhigh(ptBin-1);

				if(PlotAlteredPPDResults){
					graphDefaultRes->GetPoint(ptBin-1,ptCentre_[pt],lmeanBuff[pt]);
					lmean_errhigh[pt]=graphDefaultStat->GetErrorYhigh(ptBin-1);
					lmean_errlow[pt]=graphDefaultStat->GetErrorYlow(ptBin-1);
					lmean_errmean[pt]=(lmean_errlow[pt]+lmean_errhigh[pt])/2.;
					lmean_errmean_minus[pt]=-(lmean_errlow[pt]+lmean_errhigh[pt])/2.;
				}
				ptCentre_ForTable[pt]=ptCentre_[pt];

				if(PlotSystematics&&PlotSysSquare) {
					lmean_errmean[pt] = TMath::Power(lmean_errmean[pt],2);
					lmean_errmean_minus[pt] = -TMath::Power(lmean_errmean_minus[pt],2);
				}

				if(nSystematics>0) {graphSyst1->GetPoint(ptBin-1,Buffer[pt],SystError1[pt]);	ErrSystError1[pt]=graphSyst1->GetErrorY(pt);    if(!PlotAsymm) SystError1[pt]=TMath::Abs(SystError1[pt]); }
				if(nSystematics>1) {graphSyst2->GetPoint(ptBin-1,Buffer[pt],SystError2[pt]);	ErrSystError2[pt]=graphSyst2->GetErrorY(pt);    if(!PlotAsymm) SystError2[pt]=TMath::Abs(SystError2[pt]); }
				if(nSystematics>2) {graphSyst3->GetPoint(ptBin-1,Buffer[pt],SystError3[pt]);	ErrSystError3[pt]=graphSyst3->GetErrorY(pt);    if(!PlotAsymm) SystError3[pt]=TMath::Abs(SystError3[pt]); }
				if(nSystematics>3) {graphSyst4->GetPoint(ptBin-1,Buffer[pt],SystError4[pt]);	ErrSystError4[pt]=graphSyst4->GetErrorY(pt);    if(!PlotAsymm) SystError4[pt]=TMath::Abs(SystError4[pt]); }
				if(nSystematics>4) {graphSyst5->GetPoint(ptBin-1,Buffer[pt],SystError5[pt]);	ErrSystError5[pt]=graphSyst5->GetErrorY(pt);    if(!PlotAsymm) SystError5[pt]=TMath::Abs(SystError5[pt]); }
				if(nSystematics>5) {graphSyst6->GetPoint(ptBin-1,Buffer[pt],SystError6[pt]);	ErrSystError6[pt]=graphSyst6->GetErrorY(pt);    if(!PlotAsymm) SystError6[pt]=TMath::Abs(SystError6[pt]); }
				if(nSystematics>6) {graphSyst7->GetPoint(ptBin-1,Buffer[pt],SystError7[pt]);	ErrSystError7[pt]=graphSyst7->GetErrorY(pt);    if(!PlotAsymm) SystError7[pt]=TMath::Abs(SystError7[pt]); }
				if(nSystematics>7) {graphSyst8->GetPoint(ptBin-1,Buffer[pt],SystError8[pt]);	ErrSystError8[pt]=graphSyst8->GetErrorY(pt);    if(!PlotAsymm) SystError8[pt]=TMath::Abs(SystError8[pt]); }
				double SquaredSysts;
				if(nSystematics==1) SquaredSysts=SystError1[pt]*SystError1[pt];
				if(nSystematics==2) SquaredSysts=SystError1[pt]*SystError1[pt]+SystError2[pt]*SystError2[pt];
				if(nSystematics==3) SquaredSysts=SystError1[pt]*SystError1[pt]+SystError2[pt]*SystError2[pt]+SystError3[pt]*SystError3[pt];
				if(nSystematics==4) SquaredSysts=SystError1[pt]*SystError1[pt]+SystError2[pt]*SystError2[pt]+SystError3[pt]*SystError3[pt]+SystError4[pt]*SystError4[pt];
				if(nSystematics==5) SquaredSysts=SystError1[pt]*SystError1[pt]+SystError2[pt]*SystError2[pt]+SystError3[pt]*SystError3[pt]+SystError4[pt]*SystError4[pt]+SystError5[pt]*SystError5[pt];
				if(nSystematics==6) SquaredSysts=SystError1[pt]*SystError1[pt]+SystError2[pt]*SystError2[pt]+SystError3[pt]*SystError3[pt]+SystError4[pt]*SystError4[pt]+SystError5[pt]*SystError5[pt]+SystError6[pt]*SystError6[pt];
				if(nSystematics==7) SquaredSysts=SystError1[pt]*SystError1[pt]+SystError2[pt]*SystError2[pt]+SystError3[pt]*SystError3[pt]+SystError4[pt]*SystError4[pt]+SystError5[pt]*SystError5[pt]+SystError6[pt]*SystError6[pt]+SystError7[pt]*SystError7[pt];
				if(nSystematics==8) SquaredSysts=SystError1[pt]*SystError1[pt]+SystError2[pt]*SystError2[pt]+SystError3[pt]*SystError3[pt]+SystError4[pt]*SystError4[pt]+SystError5[pt]*SystError5[pt]+SystError6[pt]*SystError6[pt]+SystError7[pt]*SystError7[pt]+SystError8[pt]*SystError8[pt];
				SystError[pt]=TMath::Sqrt(SquaredSysts);

				SystError12[pt]=SystError1[pt]+SystError2[pt]; if(PlotAsymm) SystError12[pt]=SystError2[pt];
				SystError123[pt]=SystError1[pt]+SystError2[pt]+SystError3[pt]; if(PlotAsymm) SystError123[pt]=SystError3[pt];
				SystError1234[pt]=SystError1[pt]+SystError2[pt]+SystError3[pt]+SystError4[pt]; if(PlotAsymm) SystError1234[pt]=SystError4[pt];
				SystError12345[pt]=SystError1[pt]+SystError2[pt]+SystError3[pt]+SystError4[pt]+SystError5[pt]; if(PlotAsymm) SystError12345[pt]=SystError5[pt]; if(DeltaTildeplots) SystError12345[pt]=-SystError5[pt];//IfLamTildeClosure
				SystError123456[pt]=SystError1[pt]+SystError2[pt]+SystError3[pt]+SystError4[pt]+SystError5[pt]+SystError6[pt]; if(PlotAsymm) SystError123456[pt]=SystError6[pt];
				SystError1234567[pt]=SystError1[pt]+SystError2[pt]+SystError3[pt]+SystError4[pt]+SystError5[pt]+SystError6[pt]+SystError7[pt]; if(PlotAsymm) SystError1234567[pt]=SystError7[pt];
				SystError12345678[pt]=SystError1[pt]+SystError2[pt]+SystError3[pt]+SystError4[pt]+SystError5[pt]+SystError6[pt]+SystError7[pt]+SystError8[pt]; if(PlotAsymm) SystError12345678[pt]=SystError8[pt];

				fit_lmean_errmean[pt]=TMath::Sqrt(SystError[pt]*SystError[pt]+lmean_errmean[pt]*lmean_errmean[pt]);
				fit_lmean_errlow[pt]=TMath::Sqrt(SystError[pt]*SystError[pt]+lmean_errlow[pt]*lmean_errlow[pt]);
				fit_lmean_errhigh[pt]=TMath::Sqrt(SystError[pt]*SystError[pt]+lmean_errhigh[pt]*lmean_errhigh[pt]);

				// Fill table variables

				val_table[iLam][rapBin-1][pt]=lmean[pt];
				errHigh_table[iLam][rapBin-1][pt]=lmean_errhigh[pt];
				errLow_table[iLam][rapBin-1][pt]=lmean_errlow[pt];
				syst_table[iLam][rapBin-1][pt][0]=SystError[pt];
				syst_table[iLam][rapBin-1][pt][1]=TMath::Abs(SystError1[pt]);
				syst_table[iLam][rapBin-1][pt][2]=TMath::Abs(SystError2[pt]);
				syst_table[iLam][rapBin-1][pt][3]=TMath::Abs(SystError3[pt]);
				syst_table[iLam][rapBin-1][pt][4]=TMath::Abs(SystError4[pt]);
				syst_table[iLam][rapBin-1][pt][5]=TMath::Abs(SystError5[pt]);
				syst_table[iLam][rapBin-1][pt][6]=TMath::Abs(SystError6[pt]);
				syst_table[iLam][rapBin-1][pt][7]=TMath::Abs(SystError7[pt]);
				syst_table[iLam][rapBin-1][pt][8]=TMath::Abs(SystError8[pt]);

				errHighTotal1_table[iLam][rapBin-1][pt]=lmeanTotal1_errhigh[pt];
				errLowTotal1_table[iLam][rapBin-1][pt]=lmeanTotal1_errlow[pt];
				errHighTotal2_table[iLam][rapBin-1][pt]=lmeanTotal2_errhigh[pt];
				errLowTotal2_table[iLam][rapBin-1][pt]=lmeanTotal2_errlow[pt];
				errHighTotal3_table[iLam][rapBin-1][pt]=lmeanTotal3_errhigh[pt];
				errLowTotal3_table[iLam][rapBin-1][pt]=lmeanTotal3_errlow[pt];

				pTmean_table[iLam][rapBin-1][pt]=ptCentre_ForTable[pt];

				pt++;
			}

			gStyle->SetPalette(1,0);
			//gStyle->SetPadBottomMargin(0.12);
			//gStyle->SetPadLeftMargin(0.13);
			//gStyle->SetPadRightMargin(0.15);

			gStyle->SetPadBottomMargin(0.12); //0.11
			gStyle->SetPadLeftMargin(0.12); //0.08
			gStyle->SetPadRightMargin(0.02);
			gStyle->SetPadTopMargin(0.02);

			//gStyle->SetTitleFontSize(0.);
			gStyle->SetTitleSize(0.043,"X");
			gStyle->SetTitleSize(0.043,"Y");
			gStyle->SetLabelSize(0.045,"X");
			gStyle->SetLabelSize(0.045,"Y");

			gStyle->SetTickLength(-0.02, "xyz");
			gStyle->SetLabelOffset(0.02, "x");
			gStyle->SetLabelOffset(0.02, "y");
			gStyle->SetTitleOffset(1.3, "x");
			gStyle->SetTitleOffset(1.3, "y");
			gStyle->SetTitleFillColor(kWhite);

			if(PlotSystematics&&PlotSysSquare){
				gStyle->SetPadLeftMargin(0.15);
				gStyle->SetTitleOffset(1.7, "y");
			}

			TCanvas *plotCanvas = new TCanvas("plotCanvas","plotCanvas",1000,800);

			plotCanvas->SetFillColor(kWhite);
			//plotCanvas->SetGrid();
			plotCanvas->GetFrame()->SetFillColor(kWhite);
			plotCanvas->GetFrame()->SetBorderSize(0);
			//plotCanvas->SetRightMargin(0.05) ;


			TH1F *plotHisto = new TH1F;
			if(!PlotMatt) plotHisto = plotCanvas->DrawFrame(PlotpTMin-DeltaXminOVERALL,yMin,PlotpTMax,yMax); //to be consistant for Psi 1S and 2S
			if(PlotMatt) plotHisto = plotCanvas->DrawFrame(0.,yMin,onia::pTRange[rapBin][ptBinMax],yMax);
			if(PlotVsComp) plotHisto = plotCanvas->DrawFrame(10.1,yMin,10.7,yMax);
			plotHisto->SetXTitle("#it{p}_{T} [GeV]");
			if(PlotVsComp) plotHisto->SetXTitle("M_{#mu#mu} [GeV]");
			plotHisto->SetYTitle(axislabel);
			//plotHisto->GetYaxis()->SetTitleOffset(1.5);

			TLegend* plotcompLegend=new TLegend(0.12,0.12,0.62,0.3);
			plotcompLegend->SetFillColor(0);
			//plotcompLegend->SetTextFont(72);
			plotcompLegend->SetTextSize(0.04);
			plotcompLegend->SetBorderSize(1);
			char complegendentry[200];

			TGraphAsymmErrors *graphSyst = new TGraphAsymmErrors(nBinspT,ptCentre_,lmean,ptCentreErr_low,ptCentreErr_high,SystError,SystError);
			graphSyst->SetFillColor(kCyan-9);
			graphSyst->SetFillStyle(3001);

			char drawGraphStyle[200];
			sprintf(drawGraphStyle,"PE");



			///////////////////////////////////////////////////////////////
			//////////////Actual FinalDataResults Plotting ////////////////
			///////////////////////////////////////////////////////////////

			double BG0MarkerSize=1.;


			if(!PlotBrazilian&&!SBmSigPlots&&!BGratioFits&&!SteerIndividuals&&!PlotMatt&&!PlotVsComp) graphSyst->Draw("2");//Comment if PlotBG0plots Low
			graphDefaultRes->SetMarkerColor(ToyMC::MarkerColor[nFrame]);
			graphDefaultRes->SetLineColor(ToyMC::MarkerColor[nFrame]);
			graphDefaultRes->SetMarkerStyle(ToyMC::MarkerStyle[nState][rapBin]);
			graphDefaultRes->SetMarkerSize(ToyMC::MarkerSize[nState][rapBin]);
			if(PlotBG0plots){
				graphDefaultRes->SetMarkerStyle(20);
				graphDefaultRes->SetMarkerSize(BG0MarkerSize);

			}
			if(SetCompStyle){
				graphDefaultRes->SetMarkerStyle(24);
				graphDefaultRes->SetMarkerSize(BG0MarkerSize);

			}

			////plot noRho
			//graphDefaultRes->RemovePoint(0);
			//graphDefaultRes->RemovePoint(0);
			//graphDefaultRes->RemovePoint(0);
			//graphDefaultRes->RemovePoint(0);
			//graphDefaultRes->RemovePoint(0);
			//graphDefaultRes->RemovePoint(0);
			//graphDefaultRes->RemovePoint(0);
			//graphDefaultRes->RemovePoint(0);
			//graphDefaultRes->RemovePoint(0);

			if(!PlotBrazilian&&!SBmSigPlots&&!BGratioFits&&!SteerIndividuals&&!PlotMatt&&!PlotVsComp) graphDefaultRes->Draw(drawGraphStyle);//Comment if PlotBG0plots Low

			if(PlotBrazilian&&!SteerIndividuals&&!PlotVsComp){
				int ptOriginal;
				int nBinsOriginal;
				double ptCentre_Original[nBinspT];
				double lmean_Original[nBinspT];
				double lmean_errlow_Original[nBinspT];
				double lmean_errhigh_Original[nBinspT];
				double ptCentre_errlow_Original[nBinspT];
				double ptCentre_errhigh_Original[nBinspT];

				//nBinsOriginal=graphSyst->GetN();
				nBinsOriginal=nBinspT;
				ptOriginal=0;
				for(int ptBinOriginal=ptBinMin;ptBinOriginal<ptBinMax+1;ptBinOriginal++){
					graphSyst->GetPoint(ptBinOriginal-1,ptCentre_Original[ptOriginal],lmean_Original[ptOriginal]);
					lmean_errhigh_Original[ptOriginal]=graphSyst->GetErrorYhigh(ptBinOriginal-1);
					lmean_errlow_Original[ptOriginal]=graphSyst->GetErrorYlow(ptBinOriginal-1);
					ptCentre_errhigh_Original[ptOriginal]=graphSyst->GetErrorXhigh(ptBinOriginal-1);
					ptCentre_errlow_Original[ptOriginal]=graphSyst->GetErrorXlow(ptBinOriginal-1);
					///Alter TGraph
					ptCentre_errhigh_Original[ptOriginal]=0.;
					ptCentre_errlow_Original[ptOriginal]=0.;
					ptOriginal++;
				}
				graphSyst = new TGraphAsymmErrors(nBinsOriginal,ptCentre_Original,lmean_Original,ptCentre_errlow_Original,ptCentre_errhigh_Original,lmean_errlow_Original,lmean_errhigh_Original);

				graphSyst->SetLineColor(kBlack);
				graphSyst->SetMarkerColor(ToyMC::MarkerColor[nFrame]);
				graphSyst->SetLineColor(ToyMC::MarkerColor[nFrame]);
				graphSyst->SetMarkerStyle(ToyMC::MarkerStyle[nState][rapBin]);
				graphSyst->SetMarkerSize(ToyMC::MarkerSize[nState][rapBin]);

				ptOriginal=0;
				for(int ptBinOriginal=ptBinMin;ptBinOriginal<ptBinMax+1;ptBinOriginal++){
					graphDefaultStat->GetPoint(ptBinOriginal-1,ptCentre_Original[ptOriginal],lmean_Original[ptOriginal]);
					lmean_errhigh_Original[ptOriginal]=graphDefaultStat->GetErrorYhigh(ptBinOriginal-1);
					lmean_errlow_Original[ptOriginal]=graphDefaultStat->GetErrorYlow(ptBinOriginal-1);
					ptCentre_errhigh_Original[ptOriginal]=graphDefaultStat->GetErrorXhigh(ptBinOriginal-1);
					ptCentre_errlow_Original[ptOriginal]=graphDefaultStat->GetErrorXlow(ptBinOriginal-1);
					/// Alter TGraph
					ptCentre_errhigh_Original[ptOriginal]=0.;
					ptCentre_errlow_Original[ptOriginal]=0.;
					ptOriginal++;
				}
				graphDefaultStat = new TGraphAsymmErrors(nBinsOriginal,ptCentre_Original,lmean_Original,ptCentre_errlow_Original,ptCentre_errhigh_Original,lmean_errlow_Original,lmean_errhigh_Original);

				graphDefaultStat->SetLineColor(kBlack);
				graphDefaultStat->SetMarkerColor(ToyMC::MarkerColor[nFrame]);
				graphDefaultStat->SetLineColor(ToyMC::MarkerColor[nFrame]);
				graphDefaultStat->SetMarkerStyle(ToyMC::MarkerStyle[nState][rapBin]);
				graphDefaultStat->SetMarkerSize(ToyMC::MarkerSize[nState][rapBin]);

				ptOriginal=0;
				for(int ptBinOriginal=ptBinMin;ptBinOriginal<ptBinMax+1;ptBinOriginal++){
					graphDefaultRes->GetPoint(ptBinOriginal-1,ptCentre_Original[ptOriginal],lmean_Original[ptOriginal]);
					lmean_errhigh_Original[ptOriginal]=graphDefaultRes->GetErrorYhigh(ptBinOriginal-1);
					lmean_errlow_Original[ptOriginal]=graphDefaultRes->GetErrorYlow(ptBinOriginal-1);
					ptCentre_errhigh_Original[ptOriginal]=graphDefaultRes->GetErrorXhigh(ptBinOriginal-1);
					ptCentre_errlow_Original[ptOriginal]=graphDefaultRes->GetErrorXlow(ptBinOriginal-1);
					/// Alter TGraph
					ptCentre_errhigh_Original[ptOriginal]=ColordBandWidth;
					ptCentre_errlow_Original[ptOriginal]=ColordBandWidth;
					ptOriginal++;
				}
				graphDefaultRes = new TGraphAsymmErrors(nBinsOriginal,ptCentre_Original,lmean_Original,ptCentre_errlow_Original,ptCentre_errhigh_Original,lmean_errlow_Original,lmean_errhigh_Original);

				ptOriginal=0;
				for(int ptBinOriginal=ptBinMin;ptBinOriginal<ptBinMax+1;ptBinOriginal++){
					graphDefaultRes2sigma->GetPoint(ptBinOriginal-1,ptCentre_Original[ptOriginal],lmean_Original[ptOriginal]);
					lmean_errhigh_Original[ptOriginal]=graphDefaultRes2sigma->GetErrorYhigh(ptBinOriginal-1);
					lmean_errlow_Original[ptOriginal]=graphDefaultRes2sigma->GetErrorYlow(ptBinOriginal-1);
					ptCentre_errhigh_Original[ptOriginal]=graphDefaultRes2sigma->GetErrorXhigh(ptBinOriginal-1);
					ptCentre_errlow_Original[ptOriginal]=graphDefaultRes2sigma->GetErrorXlow(ptBinOriginal-1);
					/// Alter TGraph
					ptCentre_errhigh_Original[ptOriginal]=ColordBandWidth;
					ptCentre_errlow_Original[ptOriginal]=ColordBandWidth;
					ptOriginal++;
				}
				graphDefaultRes2sigma = new TGraphAsymmErrors(nBinsOriginal,ptCentre_Original,lmean_Original,ptCentre_errlow_Original,ptCentre_errhigh_Original,lmean_errlow_Original,lmean_errhigh_Original);

				ptOriginal=0;
				for(int ptBinOriginal=ptBinMin;ptBinOriginal<ptBinMax+1;ptBinOriginal++){
					graphDefaultRes3sigma->GetPoint(ptBinOriginal-1,ptCentre_Original[ptOriginal],lmean_Original[ptOriginal]);
					lmean_errhigh_Original[ptOriginal]=graphDefaultRes3sigma->GetErrorYhigh(ptBinOriginal-1);
					lmean_errlow_Original[ptOriginal]=graphDefaultRes3sigma->GetErrorYlow(ptBinOriginal-1);
					ptCentre_errhigh_Original[ptOriginal]=graphDefaultRes3sigma->GetErrorXhigh(ptBinOriginal-1);
					ptCentre_errlow_Original[ptOriginal]=graphDefaultRes3sigma->GetErrorXlow(ptBinOriginal-1);
					/// Alter TGraph
					ptCentre_errhigh_Original[ptOriginal]=ColordBandWidth;
					ptCentre_errlow_Original[ptOriginal]=ColordBandWidth;
					ptOriginal++;
				}
				graphDefaultRes3sigma = new TGraphAsymmErrors(nBinsOriginal,ptCentre_Original,lmean_Original,ptCentre_errlow_Original,ptCentre_errhigh_Original,lmean_errlow_Original,lmean_errhigh_Original);

				graphDefaultRes->SetFillColor(OneSigColor);
				graphDefaultRes->SetFillStyle(1001);
				graphDefaultRes2sigma->SetFillColor(TwoSigColor);
				graphDefaultRes2sigma->SetFillStyle(1001);
				graphDefaultRes3sigma->SetFillColor(ThreeSigColor);
				graphDefaultRes3sigma->SetFillStyle(1001);
				graphDefaultRes3sigma->Draw("2");
				graphDefaultRes2sigma->Draw("2");
				graphDefaultRes->Draw("2");
				if(!PlotAlteredPPDResults) graphSyst->Draw(drawGraphStyle);
				if(PlotAlteredPPDResults) graphDefaultStat->Draw(drawGraphStyle);
			} //PlotBrazilian

			if(PlotVsComp){
				//code lambda vs mass plots here
			}

			if(SteerIndividuals){ //not yet changed 12.12.2012 Linlin
				plotHisto->SetYTitle("#lambda");

				double ptCentre_SteerInd[nBinspT];
				double lmean_SteerInd[nBinspT];
				double lmean_errlow_SteerInd[nBinspT];
				double lmean_errhigh_SteerInd[nBinspT];
				double ptCentre_errlow_SteerInd[nBinspT];
				double ptCentre_errhigh_SteerInd[nBinspT];

				double SteerIndShiftOriginal=0.25;
				double SteerIndShift;

				const int nSteers=4;
				int MarkerColorSteerInd[nSteers+1]={0,632,616, 600, 432};
				int MarkerColorStyleInd[nSteers+1]={0,20,20,20,20};
				double MarkerColorSizeInd[nSteers+1]={0,1.75,1.75,1.75,1.75};

				TGraphAsymmErrors* SteerIndividualGraph1;
				TGraphAsymmErrors* SteerIndividualGraph2;
				TGraphAsymmErrors* SteerIndividualGraph3;
				TGraphAsymmErrors* SteerIndividualGraph4;

				if(rapBin==1){
					SteerIndividualGraph1 = (TGraphAsymmErrors*) infileRes->Get("ltp_PX_rap1");
					SteerIndividualGraph2 = (TGraphAsymmErrors*) infileRes->Get("ltp_PX_rap2");
					SteerIndividualGraph3 = (TGraphAsymmErrors*) infileRes->Get("ltp_CS_rap1");
					SteerIndividualGraph4 = (TGraphAsymmErrors*) infileRes->Get("ltp_CS_rap2");

					SteerIndShift=-1.5*SteerIndShiftOriginal;
					int ptSteerInd=0;
					for(int ptBinSteer=6;ptBinSteer<11;ptBinSteer++){
						SteerIndividualGraph1->GetPoint(ptBinSteer-1,ptCentre_SteerInd[ptSteerInd],lmean_SteerInd[ptSteerInd]);
						lmean_errhigh_SteerInd[ptSteerInd]=SteerIndividualGraph1->GetErrorYhigh(ptBinSteer-1);
						lmean_errlow_SteerInd[ptSteerInd]=SteerIndividualGraph1->GetErrorYlow(ptBinSteer-1);
						ptCentre_errhigh_SteerInd[ptSteerInd]=SteerIndividualGraph1->GetErrorXhigh(ptBinSteer-1);
						ptCentre_errlow_SteerInd[ptSteerInd]=SteerIndividualGraph1->GetErrorXlow(ptBinSteer-1);
						/// Alter TGraph
						ptCentre_SteerInd[ptSteerInd]=ptCentre_SteerInd[ptSteerInd]+SteerIndShift;
						ptCentre_errhigh_SteerInd[ptSteerInd]=ptCentre_errhigh_SteerInd[ptSteerInd]-SteerIndShift;
						ptCentre_errlow_SteerInd[ptSteerInd]=ptCentre_errlow_SteerInd[ptSteerInd]+SteerIndShift;

						ptSteerInd++;
					}
					SteerIndividualGraph1 = new TGraphAsymmErrors(nBinspT,ptCentre_SteerInd,lmean_SteerInd,ptCentre_errlow_SteerInd,ptCentre_errhigh_SteerInd,lmean_errlow_SteerInd,lmean_errhigh_SteerInd);

					SteerIndShift=-.5*SteerIndShiftOriginal;
					ptSteerInd=0;
					for(int ptBinSteer=6;ptBinSteer<11;ptBinSteer++){
						SteerIndividualGraph2->GetPoint(ptBinSteer-1,ptCentre_SteerInd[ptSteerInd],lmean_SteerInd[ptSteerInd]);
						lmean_errhigh_SteerInd[ptSteerInd]=SteerIndividualGraph2->GetErrorYhigh(ptBinSteer-1);
						lmean_errlow_SteerInd[ptSteerInd]=SteerIndividualGraph2->GetErrorYlow(ptBinSteer-1);
						ptCentre_errhigh_SteerInd[ptSteerInd]=SteerIndividualGraph2->GetErrorXhigh(ptBinSteer-1);
						ptCentre_errlow_SteerInd[ptSteerInd]=SteerIndividualGraph2->GetErrorXlow(ptBinSteer-1);
						/// Alter TGraph
						ptCentre_SteerInd[ptSteerInd]=ptCentre_SteerInd[ptSteerInd]+SteerIndShift;
						ptCentre_errhigh_SteerInd[ptSteerInd]=ptCentre_errhigh_SteerInd[ptSteerInd]-SteerIndShift;
						ptCentre_errlow_SteerInd[ptSteerInd]=ptCentre_errlow_SteerInd[ptSteerInd]+SteerIndShift;

						ptSteerInd++;
					}
					SteerIndividualGraph2 = new TGraphAsymmErrors(nBinspT,ptCentre_SteerInd,lmean_SteerInd,ptCentre_errlow_SteerInd,ptCentre_errhigh_SteerInd,lmean_errlow_SteerInd,lmean_errhigh_SteerInd);

					SteerIndShift=0.5*SteerIndShiftOriginal;
					ptSteerInd=0;
					for(int ptBinSteer=6;ptBinSteer<11;ptBinSteer++){
						SteerIndividualGraph3->GetPoint(ptBinSteer-1,ptCentre_SteerInd[ptSteerInd],lmean_SteerInd[ptSteerInd]);
						lmean_errhigh_SteerInd[ptSteerInd]=SteerIndividualGraph3->GetErrorYhigh(ptBinSteer-1);
						lmean_errlow_SteerInd[ptSteerInd]=SteerIndividualGraph3->GetErrorYlow(ptBinSteer-1);
						ptCentre_errhigh_SteerInd[ptSteerInd]=SteerIndividualGraph3->GetErrorXhigh(ptBinSteer-1);
						ptCentre_errlow_SteerInd[ptSteerInd]=SteerIndividualGraph3->GetErrorXlow(ptBinSteer-1);
						/// Alter TGraph
						lmean_SteerInd[ptSteerInd]=-lmean_SteerInd[ptSteerInd];
						lmean_errhigh_SteerInd[ptSteerInd]=SteerIndividualGraph3->GetErrorYlow(ptBinSteer-1);
						lmean_errlow_SteerInd[ptSteerInd]=SteerIndividualGraph3->GetErrorYhigh(ptBinSteer-1);
						ptCentre_SteerInd[ptSteerInd]=ptCentre_SteerInd[ptSteerInd]+SteerIndShift;
						ptCentre_errhigh_SteerInd[ptSteerInd]=ptCentre_errhigh_SteerInd[ptSteerInd]-SteerIndShift;
						ptCentre_errlow_SteerInd[ptSteerInd]=ptCentre_errlow_SteerInd[ptSteerInd]+SteerIndShift;

						ptSteerInd++;
					}
					SteerIndividualGraph3 = new TGraphAsymmErrors(nBinspT,ptCentre_SteerInd,lmean_SteerInd,ptCentre_errlow_SteerInd,ptCentre_errhigh_SteerInd,lmean_errlow_SteerInd,lmean_errhigh_SteerInd);

					SteerIndShift=1.5*SteerIndShiftOriginal;
					ptSteerInd=0;
					for(int ptBinSteer=6;ptBinSteer<11;ptBinSteer++){
						SteerIndividualGraph4->GetPoint(ptBinSteer-1,ptCentre_SteerInd[ptSteerInd],lmean_SteerInd[ptSteerInd]);
						lmean_errhigh_SteerInd[ptSteerInd]=SteerIndividualGraph4->GetErrorYhigh(ptBinSteer-1);
						lmean_errlow_SteerInd[ptSteerInd]=SteerIndividualGraph4->GetErrorYlow(ptBinSteer-1);
						ptCentre_errhigh_SteerInd[ptSteerInd]=SteerIndividualGraph4->GetErrorXhigh(ptBinSteer-1);
						ptCentre_errlow_SteerInd[ptSteerInd]=SteerIndividualGraph4->GetErrorXlow(ptBinSteer-1);
						/// Alter TGraph
						lmean_SteerInd[ptSteerInd]=-lmean_SteerInd[ptSteerInd];
						lmean_errhigh_SteerInd[ptSteerInd]=SteerIndividualGraph4->GetErrorYlow(ptBinSteer-1);
						lmean_errlow_SteerInd[ptSteerInd]=SteerIndividualGraph4->GetErrorYhigh(ptBinSteer-1);
						ptCentre_SteerInd[ptSteerInd]=ptCentre_SteerInd[ptSteerInd]+SteerIndShift;
						ptCentre_errhigh_SteerInd[ptSteerInd]=ptCentre_errhigh_SteerInd[ptSteerInd]-SteerIndShift;
						ptCentre_errlow_SteerInd[ptSteerInd]=ptCentre_errlow_SteerInd[ptSteerInd]+SteerIndShift;

						ptSteerInd++;
					}
					SteerIndividualGraph4 = new TGraphAsymmErrors(nBinspT,ptCentre_SteerInd,lmean_SteerInd,ptCentre_errlow_SteerInd,ptCentre_errhigh_SteerInd,lmean_errlow_SteerInd,lmean_errhigh_SteerInd);

				}

				if(rapBin==2){
					SteerIndividualGraph1 = (TGraphAsymmErrors*) infileRes->Get("lph_PX_rap1");
					SteerIndividualGraph2 = (TGraphAsymmErrors*) infileRes->Get("lph_PX_rap2");
					SteerIndividualGraph3 = (TGraphAsymmErrors*) infileRes->Get("lth_CS_rap1");
					SteerIndividualGraph4 = (TGraphAsymmErrors*) infileRes->Get("lth_CS_rap2");

					SteerIndShift=-1.5*SteerIndShiftOriginal;
					int ptSteerInd=0;
					for(int ptBinSteer=6;ptBinSteer<11;ptBinSteer++){
						SteerIndividualGraph1->GetPoint(ptBinSteer-1,ptCentre_SteerInd[ptSteerInd],lmean_SteerInd[ptSteerInd]);
						lmean_errhigh_SteerInd[ptSteerInd]=SteerIndividualGraph1->GetErrorYhigh(ptBinSteer-1);
						lmean_errlow_SteerInd[ptSteerInd]=SteerIndividualGraph1->GetErrorYlow(ptBinSteer-1);
						ptCentre_errhigh_SteerInd[ptSteerInd]=SteerIndividualGraph1->GetErrorXhigh(ptBinSteer-1);
						ptCentre_errlow_SteerInd[ptSteerInd]=SteerIndividualGraph1->GetErrorXlow(ptBinSteer-1);
						/// Alter TGraph
						ptCentre_SteerInd[ptSteerInd]=ptCentre_SteerInd[ptSteerInd]+SteerIndShift;
						ptCentre_errhigh_SteerInd[ptSteerInd]=ptCentre_errhigh_SteerInd[ptSteerInd]-SteerIndShift;
						ptCentre_errlow_SteerInd[ptSteerInd]=ptCentre_errlow_SteerInd[ptSteerInd]+SteerIndShift;

						ptSteerInd++;
					}
					SteerIndividualGraph1 = new TGraphAsymmErrors(nBinspT,ptCentre_SteerInd,lmean_SteerInd,ptCentre_errlow_SteerInd,ptCentre_errhigh_SteerInd,lmean_errlow_SteerInd,lmean_errhigh_SteerInd);

					SteerIndShift=-.5*SteerIndShiftOriginal;
					ptSteerInd=0;
					for(int ptBinSteer=6;ptBinSteer<11;ptBinSteer++){
						SteerIndividualGraph2->GetPoint(ptBinSteer-1,ptCentre_SteerInd[ptSteerInd],lmean_SteerInd[ptSteerInd]);
						lmean_errhigh_SteerInd[ptSteerInd]=SteerIndividualGraph2->GetErrorYhigh(ptBinSteer-1);
						lmean_errlow_SteerInd[ptSteerInd]=SteerIndividualGraph2->GetErrorYlow(ptBinSteer-1);
						ptCentre_errhigh_SteerInd[ptSteerInd]=SteerIndividualGraph2->GetErrorXhigh(ptBinSteer-1);
						ptCentre_errlow_SteerInd[ptSteerInd]=SteerIndividualGraph2->GetErrorXlow(ptBinSteer-1);
						/// Alter TGraph
						ptCentre_SteerInd[ptSteerInd]=ptCentre_SteerInd[ptSteerInd]+SteerIndShift;
						ptCentre_errhigh_SteerInd[ptSteerInd]=ptCentre_errhigh_SteerInd[ptSteerInd]-SteerIndShift;
						ptCentre_errlow_SteerInd[ptSteerInd]=ptCentre_errlow_SteerInd[ptSteerInd]+SteerIndShift;

						ptSteerInd++;
					}
					SteerIndividualGraph2 = new TGraphAsymmErrors(nBinspT,ptCentre_SteerInd,lmean_SteerInd,ptCentre_errlow_SteerInd,ptCentre_errhigh_SteerInd,lmean_errlow_SteerInd,lmean_errhigh_SteerInd);

					SteerIndShift=.5*SteerIndShiftOriginal;
					ptSteerInd=0;
					for(int ptBinSteer=6;ptBinSteer<11;ptBinSteer++){
						SteerIndividualGraph3->GetPoint(ptBinSteer-1,ptCentre_SteerInd[ptSteerInd],lmean_SteerInd[ptSteerInd]);
						lmean_errhigh_SteerInd[ptSteerInd]=SteerIndividualGraph3->GetErrorYhigh(ptBinSteer-1);
						lmean_errlow_SteerInd[ptSteerInd]=SteerIndividualGraph3->GetErrorYlow(ptBinSteer-1);
						ptCentre_errhigh_SteerInd[ptSteerInd]=SteerIndividualGraph3->GetErrorXhigh(ptBinSteer-1);
						ptCentre_errlow_SteerInd[ptSteerInd]=SteerIndividualGraph3->GetErrorXlow(ptBinSteer-1);
						/// Alter TGraph
						lmean_SteerInd[ptSteerInd]=-lmean_SteerInd[ptSteerInd];
						lmean_errhigh_SteerInd[ptSteerInd]=SteerIndividualGraph4->GetErrorYlow(ptBinSteer-1);
						lmean_errlow_SteerInd[ptSteerInd]=SteerIndividualGraph4->GetErrorYhigh(ptBinSteer-1);
						ptCentre_SteerInd[ptSteerInd]=ptCentre_SteerInd[ptSteerInd]+SteerIndShift;
						ptCentre_errhigh_SteerInd[ptSteerInd]=ptCentre_errhigh_SteerInd[ptSteerInd]-SteerIndShift;
						ptCentre_errlow_SteerInd[ptSteerInd]=ptCentre_errlow_SteerInd[ptSteerInd]+SteerIndShift;

						ptSteerInd++;
					}
					SteerIndividualGraph3 = new TGraphAsymmErrors(nBinspT,ptCentre_SteerInd,lmean_SteerInd,ptCentre_errlow_SteerInd,ptCentre_errhigh_SteerInd,lmean_errlow_SteerInd,lmean_errhigh_SteerInd);

					SteerIndShift=1.5*SteerIndShiftOriginal;
					ptSteerInd=0;
					for(int ptBinSteer=6;ptBinSteer<11;ptBinSteer++){
						SteerIndividualGraph4->GetPoint(ptBinSteer-1,ptCentre_SteerInd[ptSteerInd],lmean_SteerInd[ptSteerInd]);
						lmean_errhigh_SteerInd[ptSteerInd]=SteerIndividualGraph4->GetErrorYhigh(ptBinSteer-1);
						lmean_errlow_SteerInd[ptSteerInd]=SteerIndividualGraph4->GetErrorYlow(ptBinSteer-1);
						ptCentre_errhigh_SteerInd[ptSteerInd]=SteerIndividualGraph4->GetErrorXhigh(ptBinSteer-1);
						ptCentre_errlow_SteerInd[ptSteerInd]=SteerIndividualGraph4->GetErrorXlow(ptBinSteer-1);
						/// Alter TGraph
						lmean_SteerInd[ptSteerInd]=-lmean_SteerInd[ptSteerInd];
						lmean_errhigh_SteerInd[ptSteerInd]=SteerIndividualGraph4->GetErrorYlow(ptBinSteer-1);
						lmean_errlow_SteerInd[ptSteerInd]=SteerIndividualGraph4->GetErrorYhigh(ptBinSteer-1);
						ptCentre_SteerInd[ptSteerInd]=ptCentre_SteerInd[ptSteerInd]+SteerIndShift;
						ptCentre_errhigh_SteerInd[ptSteerInd]=ptCentre_errhigh_SteerInd[ptSteerInd]-SteerIndShift;
						ptCentre_errlow_SteerInd[ptSteerInd]=ptCentre_errlow_SteerInd[ptSteerInd]+SteerIndShift;

						ptSteerInd++;
					}
					SteerIndividualGraph4 = new TGraphAsymmErrors(nBinspT,ptCentre_SteerInd,lmean_SteerInd,ptCentre_errlow_SteerInd,ptCentre_errhigh_SteerInd,lmean_errlow_SteerInd,lmean_errhigh_SteerInd);


				}

				SteerIndividualGraph1->SetMarkerColor(MarkerColorSteerInd[1]);
				SteerIndividualGraph1->SetLineColor(MarkerColorSteerInd[1]);
				SteerIndividualGraph1->SetMarkerStyle(MarkerColorStyleInd[1]);
				SteerIndividualGraph1->SetMarkerSize(MarkerColorSizeInd[1]);

				SteerIndividualGraph2->SetMarkerColor(MarkerColorSteerInd[2]);
				SteerIndividualGraph2->SetLineColor(MarkerColorSteerInd[2]);
				SteerIndividualGraph2->SetMarkerStyle(MarkerColorStyleInd[2]);
				SteerIndividualGraph2->SetMarkerSize(MarkerColorSizeInd[2]);

				SteerIndividualGraph3->SetMarkerColor(MarkerColorSteerInd[3]);
				SteerIndividualGraph3->SetLineColor(MarkerColorSteerInd[3]);
				SteerIndividualGraph3->SetMarkerStyle(MarkerColorStyleInd[3]);
				SteerIndividualGraph3->SetMarkerSize(MarkerColorSizeInd[3]);

				SteerIndividualGraph4->SetMarkerColor(MarkerColorSteerInd[4]);
				SteerIndividualGraph4->SetLineColor(MarkerColorSteerInd[4]);
				SteerIndividualGraph4->SetMarkerStyle(MarkerColorStyleInd[4]);
				SteerIndividualGraph4->SetMarkerSize(MarkerColorSizeInd[4]);

				SteerIndividualGraph1->Draw(drawGraphStyle);
				SteerIndividualGraph2->Draw(drawGraphStyle);
				SteerIndividualGraph3->Draw(drawGraphStyle);
				SteerIndividualGraph4->Draw(drawGraphStyle);

				TLegend* plotSteerIndLegend=new TLegend(0.55,0.7,0.95,0.9);
				plotSteerIndLegend->SetFillColor(0);
				//plotSteerIndLegend->SetTextFont(72);
				plotSteerIndLegend->SetTextSize(0.04);
				plotSteerIndLegend->SetBorderSize(1);
				char SteerIndlegendentry[200];

				if(rapBin==1){
					sprintf(SteerIndlegendentry,"#lambda_{#vartheta#varphi}, |y|<0.6, PX");
					plotSteerIndLegend->AddEntry(SteerIndividualGraph1,SteerIndlegendentry,"elp");
					sprintf(SteerIndlegendentry,"#lambda_{#vartheta#varphi}, 0.6<|y|<1.2, PX");
					plotSteerIndLegend->AddEntry(SteerIndividualGraph2,SteerIndlegendentry,"elp");
					sprintf(SteerIndlegendentry,"-#lambda_{#vartheta#varphi}, |y|<0.6, CS");
					plotSteerIndLegend->AddEntry(SteerIndividualGraph3,SteerIndlegendentry,"elp");
					sprintf(SteerIndlegendentry,"-#lambda_{#vartheta#varphi}, 0.6<|y|<1.2, CS");
					plotSteerIndLegend->AddEntry(SteerIndividualGraph4,SteerIndlegendentry,"elp");
					plotSteerIndLegend->Draw();
				}

				if(rapBin==2){
					sprintf(SteerIndlegendentry,"#lambda_{#varphi}, |y|<0.6, PX");
					plotSteerIndLegend->AddEntry(SteerIndividualGraph1,SteerIndlegendentry,"elp");
					sprintf(SteerIndlegendentry,"#lambda_{#varphi}, 0.6<|y|<1.2, PX");
					plotSteerIndLegend->AddEntry(SteerIndividualGraph2,SteerIndlegendentry,"elp");
					sprintf(SteerIndlegendentry,"-#lambda_{#vartheta}, |y|<0.6, CS");
					plotSteerIndLegend->AddEntry(SteerIndividualGraph3,SteerIndlegendentry,"elp");
					sprintf(SteerIndlegendentry,"-#lambda_{#vartheta}, 0.6<|y|<1.2, CS");
					plotSteerIndLegend->AddEntry(SteerIndividualGraph4,SteerIndlegendentry,"elp");
					plotSteerIndLegend->Draw();
				}


				TLine* extreme0 = new TLine( PlotpTMin-DeltaXminOVERALL, 0, PlotpTMax ,0);
				extreme0->SetLineWidth( 2 );
				extreme0->SetLineStyle( 1 );
				extreme0->SetLineColor( kBlack );
				extreme0->Draw( "same" );
			}// SteerIndividuals

			if(PlotMatt&&!PlotMattForICHEP){ //not yet changed 12.12.2012 Linlin
				graphMattTotal->SetLineColor(kCyan);
				graphMattTotal->SetFillColor(kCyan);
				graphMattTotal->SetFillStyle(3002);//3002
				graphMattTotal->SetLineWidth(0.);

				graphMattStat->SetMarkerColor(kRed);
				graphMattStat->SetLineColor(kRed);
				graphMattStat->SetMarkerStyle(25);

				graphDefaultStat->SetLineColor(kBlack);
				graphDefaultStat->SetMarkerColor(ToyMC::MarkerColor[nFrame]);
				graphDefaultStat->SetLineColor(ToyMC::MarkerColor[nFrame]);
				graphDefaultStat->SetMarkerStyle(ToyMC::MarkerStyle[nState][rapBin]);
				graphDefaultStat->SetMarkerSize(ToyMC::MarkerSize[nState][rapBin]);

				graphDefaultRes->SetFillColor(kGreen);
				graphDefaultRes->SetLineColor(kGreen);
				graphDefaultRes->SetFillStyle(3002);

				graphMattTotal->Draw("2");
				graphDefaultRes->Draw("2");

				graphMattStat->Draw(drawGraphStyle);
				if(PlotAlteredPPDResults) graphDefaultStat->Draw(drawGraphStyle);

				TLegend* plotMattLegend=new TLegend(0.4,0.7,0.95,0.9);
				plotMattLegend->SetFillColor(0);
				//plotMattLegend->SetTextFont(72);
				plotMattLegend->SetTextSize(0.04);
				plotMattLegend->SetBorderSize(1);
				char Mattlegendentry[200];
				sprintf(Mattlegendentry,"CMS tot.uncert., 68.3%% CL");
				plotMattLegend->AddEntry(graphDefaultRes,Mattlegendentry,"f");
				sprintf(Mattlegendentry,"CMS stat.uncert., 68.3%% CL");
				plotMattLegend->AddEntry(graphDefaultStat,Mattlegendentry,"elp");
				sprintf(Mattlegendentry,"CDF tot.uncert., 68.3%% CL");
				plotMattLegend->AddEntry(graphMattTotal,Mattlegendentry,"f");
				sprintf(Mattlegendentry,"CDF stat.uncert., 68.3%% CL");
				plotMattLegend->AddEntry(graphMattStat,Mattlegendentry,"elp");
				plotMattLegend->Draw();
			} // PlotMatt&&!PlotMattForICHEP

			if(PlotMatt&&PlotMattForICHEP){ //not yet changed 12.12.2012 Linlin
				plotCanvas->SetRightMargin(0.05);
				plotCanvas->SetTopMargin(0.05);

				TLine* extreme0 = new TLine( 0, 0, PlotpTMax ,0);
				extreme0->SetLineWidth( 2 );
				extreme0->SetLineStyle( 2 );
				extreme0->SetLineColor( kBlack );
				extreme0->Draw( "same" );

				graphMattTotal->SetMarkerColor(kGreen+2);
				graphMattTotal->SetLineColor(kGreen+2);
				graphMattTotal->SetMarkerStyle(20);
				graphMattTotal->SetMarkerSize(1.5);

				graphDefaultRes->SetLineColor(kBlue);
				graphDefaultRes->SetMarkerColor(kBlue);
				graphDefaultRes->SetMarkerStyle(20);
				graphDefaultRes->SetMarkerSize(1.5);

				graphNLONRQCD->SetLineColor(kRed);
				graphNNLO->SetLineColor(616);
				graphNNLO->SetLineStyle(2);

				graphMattTotal->Draw(drawGraphStyle);
				graphDefaultRes->Draw(drawGraphStyle);

				double ICHEPlegendYMin=0.15;
				double ICHEPlegendYMax=0.275;
				if(nFrame==3) ICHEPlegendYMax=ICHEPlegendYMin+(ICHEPlegendYMax-ICHEPlegendYMin)*1/3;
				if(nState==3) ICHEPlegendYMax=ICHEPlegendYMin+(ICHEPlegendYMax-ICHEPlegendYMin)*4/3;
				if(nState==1) ICHEPlegendYMax=ICHEPlegendYMin+(ICHEPlegendYMax-ICHEPlegendYMin)*2/3;

				TLegend* plotICHEPLegend=new TLegend(0.1375,ICHEPlegendYMin-0.025,0.75,ICHEPlegendYMax);
				plotICHEPLegend->SetFillColor(0);
				//plotICHEPLegend->SetTextFont(72);
				plotICHEPLegend->SetTextSize(0.039);
				plotICHEPLegend->SetBorderSize(0);
				//plotICHEPLegend->SetHeader("Data error bars: tot. uncert., 68.3% CL");
				char ICHEPlegendentry[200];
				sprintf(ICHEPlegendentry,"CMS, tot. uncert., 68.3%% CL");
				plotICHEPLegend->AddEntry(graphDefaultRes,ICHEPlegendentry,"elp");
				sprintf(ICHEPlegendentry,"CDF PRL 108, 151802 (2012), tot. uncert., 68.3%% CL");
				if(nFrame!=3) plotICHEPLegend->AddEntry(graphMattTotal,ICHEPlegendentry,"elp");
				sprintf(ICHEPlegendentry,"NLO NRQCD at  #sqrt{s} = 1.96 TeV,  PRD83, 114021 (2011)");
				if(nState==3) plotICHEPLegend->AddEntry(graphNLONRQCD,ICHEPlegendentry,"l");
				sprintf(ICHEPlegendentry,"NNLO* CSM at  #sqrt{s} = 1.8 TeV, PRL101, 152001 (2008)");
				if(nState==3) plotICHEPLegend->AddEntry(graphNNLO,ICHEPlegendentry,"l");
				plotICHEPLegend->SetMargin(0.135);
				plotICHEPLegend->Draw();

				double ICHEPFontSize=0.04;

				char ICHEPtext[200];
				sprintf(ICHEPtext,"#psi(%dS)",nState-3);
				TLatex *texICHEP1 = new TLatex(40.5,yMin+(yMax-yMin)*0.8125,ICHEPtext);
				texICHEP1->SetTextSize(ICHEPFontSize*1.75);
				texICHEP1->Draw( "same" );

				sprintf(ICHEPtext,"CMS  pp  #sqrt{s} = 7 TeV  L = 4.9 fb^{-1}");
				TLatex *texICHEP2 = new TLatex(2,yMin+(yMax-yMin)*0.925,ICHEPtext);
				texICHEP2->SetTextSize(ICHEPFontSize);
				texICHEP2->Draw( "same" );

				sprintf(ICHEPtext,"CDF  p#bar{p}  #sqrt{s} = 1.96 TeV");
				TLatex *texICHEP3 = new TLatex(2,yMin+(yMax-yMin)*0.85,ICHEPtext);
				texICHEP3->SetTextSize(ICHEPFontSize);
				if(nFrame!=3) texICHEP3->Draw( "same" );

				double xICHEPrap=33.5;
				if(rapBin==2) xICHEPrap=28.;

				if(nFrame==1&&rapBin==1) sprintf(ICHEPtext,"CS frame, |#it{y}| < 0.6");
				if(nFrame==2&&rapBin==1) sprintf(ICHEPtext,"HX frame, |#it{y}| < 0.6");
				if(nFrame==3&&rapBin==1) sprintf(ICHEPtext,"PX frame, |#it{y}| < 0.6");
				if(nFrame==1&&rapBin==2) sprintf(ICHEPtext,"CS frame, 0.6 < |#it{y}| < 1.2");
				if(nFrame==2&&rapBin==2) sprintf(ICHEPtext,"HX frame, 0.6 < |#it{y}| < 1.2");
				if(nFrame==3&&rapBin==2) sprintf(ICHEPtext,"PX frame, 0.6 < |#it{y}| < 1.2");
				if(nFrame==1&&rapBin==3) sprintf(ICHEPtext,"CS frame, 1.2 < |#it{y}| < 1.5");
				if(nFrame==2&&rapBin==3) sprintf(ICHEPtext,"HX frame, 1.2 < |#it{y}| < 1.5");
				if(nFrame==3&&rapBin==3) sprintf(ICHEPtext,"PX frame, 1.2 < |#it{y}| < 1.5");
				TLatex *texICHEP4 = new TLatex(xICHEPrap,yMin+(yMax-yMin)*0.925,ICHEPtext);
				texICHEP4->SetTextSize(ICHEPFontSize);
				texICHEP4->Draw( "same" );

				plotHisto->GetYaxis()->SetTitleOffset(15);

				//sprintf(ICHEPtext,"%s",axislabel);
				sprintf(ICHEPtext,"#lambda_{#vartheta}");
				TLatex *texICHEP5 = new TLatex(-7.65,yMin+(yMax-yMin)*0.485,ICHEPtext);
				texICHEP5->SetTextSize(ICHEPFontSize*1.75);
				texICHEP5->Draw( "same" );
			} //PlotMatt&&PlotMattForICHEP

			if(CompareSyst){
				TLine* extreme0 = new TLine( PlotpTMin-DeltaXminOVERALL, 0, PlotpTMax ,0);
				extreme0->SetLineWidth( 2 );
				extreme0->SetLineStyle( 1 );
				extreme0->SetLineColor( kBlack );
				extreme0->Draw( "same" );
			}//CompareSyst

			sprintf(complegendentry,"%s",LegendEntryDefID);//FindLegend
			if(!SBmSigPlots&&!BGratioFits) plotcompLegend->AddEntry(graphDefaultRes,complegendentry,"lp");

			if(FitGraph){
				TGraphAsymmErrors *fitGraph = new TGraphAsymmErrors(nBinspT,ptCentre_,lmean,0,0,fit_lmean_errlow,fit_lmean_errhigh);
				//TGraphErrors *fitGraph = new TGraphErrors(nBinspT,ptCentre_,lmean,0,fit_lmean_errmean);

				double pTmaxPlot=72.; //50;
				double pTminPlot=6.;  //10;

				double FontSize=0.0215;
				double xText=15;

				char FitOptions[200];
				sprintf(FitOptions,"EFNR");

				TF1* fLin = new TF1("fLin","pol1",pTminPlot,pTmaxPlot);
				fitGraph->Fit("fLin",FitOptions);
				fLin->SetLineWidth(0.4);
				fLin->SetLineColor(kRed);
				fLin->SetLineStyle(1);
				fLin->Draw("same");

				double Lin_p0 = fLin->GetParameter(0);
				double err_Lin_p0 = fLin->GetParError(0);
				double Lin_p1 = fLin->GetParameter(1);
				double err_Lin_p1 = fLin->GetParError(1);
				double Lin_chi2=fLin->GetChisquare();
				int Lin_NDF=fLin->GetNDF();
				double Lin_BIC=Lin_chi2+2*TMath::Log(5.);

				double highest=yMax-(yMax-yMin)*0.05;

				char text[200];
				sprintf(text,"#color[2]{Fitting pol1:} p0 = %1.4f #pm %1.4f, p1 = %1.4f #pm %1.4f, #chi^{2}/ndf = %1.4f / %d, BIC = %1.3f", Lin_p0, err_Lin_p0, Lin_p1, err_Lin_p1,Lin_chi2,Lin_NDF,Lin_BIC);
				TLatex *textLin1 = new TLatex(pTminPlot*1.2,yMin+(yMax-yMin)*0.95,text);
				textLin1->SetTextSize(FontSize)                                                                                                                                                                                                                                             ;
				textLin1->Draw( "same" )                                                                                                                                                                                                                                                 ;

				TF1* fConst = new TF1("fConst","pol0",pTminPlot,pTmaxPlot);
				fitGraph->Fit("fConst",FitOptions);
				fConst->SetLineWidth(0.4);
				fConst->SetLineStyle(1);
				fConst->SetLineColor(kBlue);
				fConst->Draw("same");

				double Const_p0 = fConst->GetParameter(0);
				double err_Const_p0 = fConst->GetParError(0);
				double Const_chi2=fConst->GetChisquare();
				int Const_NDF=fConst->GetNDF();
				double Const_BIC=Const_chi2+1*TMath::Log(5.);


				sprintf(text,"#color[4]{Fitting constant:} const. = %1.4f #pm %1.4f, #chi^{2}/ndf = %1.4f / %d, BIC = %1.3f", Const_p0, err_Const_p0, Const_chi2, Const_NDF,Const_BIC);
				TLatex *textConst1 = new TLatex(pTminPlot*1.2,yMin+(yMax-yMin)*0.9,text);
				textConst1->SetTextSize(FontSize)                                                                                                                                                                                                                                             ;
				textConst1->Draw( "same" )                                                                                                                                                                                                                                                 ;

				sprintf(text,"H_{0} (unpolarized) test: %1.1f#sigma", Const_p0/err_Const_p0);
				TLatex *textConst2 = new TLatex(pTminPlot*1.2,yMin+(yMax-yMin)*0.8,text);
				textConst2->SetTextSize(1.75*FontSize)                                                                                                                                                                                                                                             ;
				textConst2->Draw( "same" )                                                                                                                                                                                                                                                 ;

				sprintf(text,"H_{1} (constant) test: %1.1f#sigma", Lin_p1/err_Lin_p1);
				TLatex *textConst2p5 = new TLatex(pTminPlot*1.2,yMin+(yMax-yMin)*0.725,text);
				textConst2p5->SetTextSize(1.75*FontSize)                                                                                                                                                                                                                                             ;
				textConst2p5->Draw( "same" )                                                                                                                                                                                                                                                 ;

				if(Const_BIC<Lin_BIC) sprintf(text,"BIC prefers #color[4]{constant} fit");
				else sprintf(text,"BIC prefers #color[2]{linear} fit");
				TLatex *textConst3 = new TLatex(pTminPlot*1.2,yMin+(yMax-yMin)*0.65,text);
				textConst3->SetTextSize(1.75*FontSize)                                                                                                                                                                                                                                             ;
				textConst3->Draw( "same" )                                                                                                                                                                                                                                                 ;
			} //FitGraph


			if(SBmSigPlots){ //not yet changed 12.12.2012 Linlin
				TGraphAsymmErrors *LSBgraph;
				TGraphAsymmErrors *RSBgraph;

				bool CombineRapsForSB=true;

				double FreePoints=5.;
				if(CombineRapsForSB) FreePoints=10.;

				if(!CombineRapsForSB){
					LSBgraph = graphCompareFile1;
					RSBgraph = graphCompareFile2;
				} 


				if(CombineRapsForSB){

					TGraphAsymmErrors* graphCompareFile1_OtherRap = (TGraphAsymmErrors*) CompareFile1->Get(GraphNameOtherRap);
					TGraphAsymmErrors* graphCompareFile2_OtherRap = (TGraphAsymmErrors*) CompareFile2->Get(GraphNameOtherRap);

					const int nBinspTSB=14;

					double ptCentre_SB_LSB[nBinspTSB];
					double lmean_SB_LSB[nBinspTSB];
					double fit_lmean_errlow_SB_LSB[nBinspTSB];
					double fit_lmean_errhigh_SB_LSB[nBinspTSB];

					double ptCentre_SB_RSB[nBinspTSB];
					double lmean_SB_RSB[nBinspTSB];
					double fit_lmean_errlow_SB_RSB[nBinspTSB];
					double fit_lmean_errhigh_SB_RSB[nBinspTSB];

					int runningSB=0;
					for(int rapBinSB=1;rapBinSB<3;rapBinSB++){
						for(int ptBinSB=6;ptBinSB<13;ptBinSB++){

							if(rapBinSB==1){
								graphCompareFile1_OtherRap->GetPoint(ptBinSB-1,ptCentre_SB_LSB[runningSB],lmean_SB_LSB[runningSB]);
								fit_lmean_errhigh_SB_LSB[runningSB]=graphCompareFile1_OtherRap->GetErrorYhigh(ptBinSB-1);
								fit_lmean_errlow_SB_LSB[runningSB]=graphCompareFile1_OtherRap->GetErrorYlow(ptBinSB-1);
							}
							if(rapBinSB==2){
								graphCompareFile1->GetPoint(ptBinSB-1,ptCentre_SB_LSB[runningSB],lmean_SB_LSB[runningSB]);
								fit_lmean_errhigh_SB_LSB[runningSB]=graphCompareFile1->GetErrorYhigh(ptBinSB-1);
								fit_lmean_errlow_SB_LSB[runningSB]=graphCompareFile1->GetErrorYlow(ptBinSB-1);
							}

							if(rapBinSB==1){
								graphCompareFile2_OtherRap->GetPoint(ptBinSB-1,ptCentre_SB_RSB[runningSB],lmean_SB_RSB[runningSB]);
								fit_lmean_errhigh_SB_RSB[runningSB]=graphCompareFile2_OtherRap->GetErrorYhigh(ptBinSB-1);
								fit_lmean_errlow_SB_RSB[runningSB]=graphCompareFile2_OtherRap->GetErrorYlow(ptBinSB-1);
							}
							if(rapBinSB==2){
								graphCompareFile2->GetPoint(ptBinSB-1,ptCentre_SB_RSB[runningSB],lmean_SB_RSB[runningSB]);
								fit_lmean_errhigh_SB_RSB[runningSB]=graphCompareFile2->GetErrorYhigh(ptBinSB-1);
								fit_lmean_errlow_SB_RSB[runningSB]=graphCompareFile2->GetErrorYlow(ptBinSB-1);
							}

							runningSB++;
						}
					}


					LSBgraph = new TGraphAsymmErrors(nBinspTSB,ptCentre_SB_LSB,lmean_SB_LSB,0,0,fit_lmean_errlow_SB_LSB,fit_lmean_errhigh_SB_LSB);
					RSBgraph = new TGraphAsymmErrors(nBinspTSB,ptCentre_SB_RSB,lmean_SB_RSB,0,0,fit_lmean_errlow_SB_RSB,fit_lmean_errhigh_SB_RSB);

					char savequick[200];
					sprintf(savequick,"tmp/LSBgraph_rap%d.root",rapBin);
					LSBgraph->SaveAs(savequick);

				}


				double pTmaxPlot=72.; //50;
				double pTminPlot=6.;  //10;

				double FontSize=0.0215;
				double xText=15;
				char text[200];

				char FitOptions[200];
				sprintf(FitOptions,"EFNR");


				bool DrawSBLatex=false;
				bool MonotonicityConst=true;
				bool MonotonicityLin=false;


				if(MonotonicityConst){


					TF1* fConstLSB = new TF1("fConstLSB","pol0",pTminPlot,pTmaxPlot);
					LSBgraph->Fit("fConstLSB",FitOptions);
					fConstLSB->SetLineWidth(0.4);
					fConstLSB->SetLineStyle(1);
					fConstLSB->SetLineColor(kRed);
					fConstLSB->Draw("same");

					double ConstLSB_p0 = fConstLSB->GetParameter(0);
					double err_ConstLSB_p0 = fConstLSB->GetParError(0);
					double ConstLSB_chi2=fConstLSB->GetChisquare();
					int ConstLSB_NDF=fConstLSB->GetNDF();
					double ConstLSB_BIC=ConstLSB_chi2+1*TMath::Log(FreePoints);

					TF1* fConstLSB_m = new TF1("fConstLSB_m","pol0",pTminPlot,pTmaxPlot);
					fConstLSB_m->SetParameter(0,ConstLSB_p0+err_ConstLSB_p0);
					fConstLSB_m->SetLineWidth(0.2);
					fConstLSB_m->SetLineColor(kRed);
					fConstLSB_m->SetLineStyle(2);
					fConstLSB_m->Draw("same");
					TF1* fConstLSB_p = new TF1("fConstLSB_p","pol0",pTminPlot,pTmaxPlot);
					fConstLSB_p->SetParameter(0,ConstLSB_p0-err_ConstLSB_p0);
					fConstLSB_p->SetLineWidth(0.2);
					fConstLSB_p->SetLineColor(kRed);
					fConstLSB_p->SetLineStyle(2);
					fConstLSB_p->Draw("same");

					if(DrawSBLatex){
						sprintf(text,"#color[2]{Fitting constant, LSB:} const. = %1.4f #pm %1.4f, #chi^{2}/ndf = %1.4f / %d, BIC = %1.3f", ConstLSB_p0, err_ConstLSB_p0, ConstLSB_chi2, ConstLSB_NDF,ConstLSB_BIC);
						TLatex *textConstLSB = new TLatex(pTminPlot*1.2,yMin+(yMax-yMin)*0.95,text);
						textConstLSB->SetTextSize(FontSize)                                                                                                                                                                                                                                             ;
						textConstLSB->Draw( "same" )                                                                                                                                                                                                                                                 ;
					}


					TF1* fConstRSB = new TF1("fConstRSB","pol0",pTminPlot,pTmaxPlot);
					RSBgraph->Fit("fConstRSB",FitOptions);
					fConstRSB->SetLineWidth(0.4);
					fConstRSB->SetLineStyle(1);
					fConstRSB->SetLineColor(kBlue);
					fConstRSB->Draw("same");

					double ConstRSB_p0 = fConstRSB->GetParameter(0);
					double err_ConstRSB_p0 = fConstRSB->GetParError(0);
					double ConstRSB_chi2=fConstRSB->GetChisquare();
					int ConstRSB_NDF=fConstRSB->GetNDF();
					double ConstRSB_BIC=ConstRSB_chi2+1*TMath::Log(FreePoints);

					TF1* fConstRSB_m = new TF1("fConstRSB_m","pol0",pTminPlot,pTmaxPlot);
					fConstRSB_m->SetParameter(0,ConstRSB_p0+err_ConstRSB_p0);
					fConstRSB_m->SetLineWidth(0.2);
					fConstRSB_m->SetLineColor(kBlue);
					fConstRSB_m->SetLineStyle(2);
					fConstRSB_m->Draw("same");
					TF1* fConstRSB_p = new TF1("fConstRSB_p","pol0",pTminPlot,pTmaxPlot);
					fConstRSB_p->SetParameter(0,ConstRSB_p0-err_ConstRSB_p0);
					fConstRSB_p->SetLineWidth(0.2);
					fConstRSB_p->SetLineColor(kBlue);
					fConstRSB_p->SetLineStyle(2);
					fConstRSB_p->Draw("same");

					if(DrawSBLatex){
						sprintf(text,"#color[4]{Fitting constant, RSB:} const. = %1.4f #pm %1.4f, #chi^{2}/ndf = %1.4f / %d, BIC = %1.3f", ConstRSB_p0, err_ConstRSB_p0, ConstRSB_chi2, ConstRSB_NDF,ConstRSB_BIC);
						TLatex *textConstRSB = new TLatex(pTminPlot*1.2,yMin+(yMax-yMin)*0.9,text);
						textConstRSB->SetTextSize(FontSize)                                                                                                                                                                                                                                             ;
						textConstRSB->Draw( "same" )                                                                                                                                                                                                                                                 ;
					}



					bool LSBaboveSignal=false;
					bool LSBaboveSignalSignificantly=false;
					bool LSBbelowSignal=false;
					bool LSBbelowSignalSignificantly=false;

					if(ConstLSB_p0>0){
						LSBaboveSignal=true;
						LSBbelowSignal=false;
						if(ConstLSB_p0-err_ConstLSB_p0>0) LSBaboveSignalSignificantly=true;
						if(ConstLSB_p0-err_ConstLSB_p0<0) LSBaboveSignalSignificantly=false;
					}
					if(ConstLSB_p0<0){
						LSBaboveSignal=false;
						LSBbelowSignal=true;
						if(ConstLSB_p0+err_ConstLSB_p0<0) LSBbelowSignalSignificantly=true;
						if(ConstLSB_p0+err_ConstLSB_p0>0) LSBbelowSignalSignificantly=false;
					}

					bool RSBaboveSignal=false;
					bool RSBaboveSignalSignificantly=false;
					bool RSBbelowSignal=false;
					bool RSBbelowSignalSignificantly=false;

					if(ConstRSB_p0>0){
						RSBaboveSignal=true;
						RSBbelowSignal=false;
						if(ConstRSB_p0-err_ConstRSB_p0>0) RSBaboveSignalSignificantly=true;
						if(ConstRSB_p0-err_ConstRSB_p0<0) RSBaboveSignalSignificantly=false;
					}
					if(ConstRSB_p0<0){
						RSBaboveSignal=false;
						RSBbelowSignal=true;
						if(ConstRSB_p0+err_ConstRSB_p0<0) RSBbelowSignalSignificantly=true;
						if(ConstRSB_p0+err_ConstRSB_p0>0) RSBbelowSignalSignificantly=false;
					}

					bool MonotonicityInconclusive=false;
					bool MonotonicityConclusivePositive=false;
					bool MonotonicityConclusiveNegative=false;

					if(!RSBbelowSignalSignificantly||!LSBbelowSignalSignificantly){
						MonotonicityInconclusive=true;
						MonotonicityConclusivePositive=false;
						MonotonicityConclusiveNegative=false;
					}
					if(RSBbelowSignalSignificantly&&LSBaboveSignalSignificantly) MonotonicityConclusivePositive=true;
					if(RSBaboveSignalSignificantly&&LSBbelowSignalSignificantly) MonotonicityConclusivePositive=true;

					if(RSBbelowSignalSignificantly&&LSBbelowSignalSignificantly) MonotonicityConclusiveNegative=true;
					if(RSBaboveSignalSignificantly&&LSBaboveSignalSignificantly) MonotonicityConclusiveNegative=true;


					if(MonotonicityInconclusive) sprintf(text,"Inconclusive");
					if(MonotonicityConclusivePositive) sprintf(text,"Monotonic behaviour");
					if(MonotonicityConclusiveNegative) sprintf(text,"Non-monotonicity behaviour");
					TLatex *textMonotonicityResult = new TLatex(pTminPlot*1.2,yMin+(yMax-yMin)*0.9,text);
					textMonotonicityResult->SetTextSize(1.75*FontSize)                                                                                                                                                                                                                                             ;
					textMonotonicityResult->Draw( "same" )                                                                                                                                                                                                                                                 ;


				}

				///////////////////////////////////////////////////

				if(MonotonicityLin){


					TF1* fLinLSB = new TF1("fLinLSB","pol1",pTminPlot,pTmaxPlot);
					LSBgraph->Fit("fLinLSB",FitOptions);
					fLinLSB->SetLineWidth(0.4);
					fLinLSB->SetLineColor(kRed);
					fLinLSB->SetLineStyle(1);
					fLinLSB->Draw("same");

					double LinLSB_p0 = fLinLSB->GetParameter(0);
					double err_LinLSB_p0 = fLinLSB->GetParError(0);
					double LinLSB_p1 = fLinLSB->GetParameter(1);
					double err_LinLSB_p1 = fLinLSB->GetParError(1);
					double LinLSB_chi2=fLinLSB->GetChisquare();
					int LinLSB_NDF=fLinLSB->GetNDF();
					double LinLSB_BIC=LinLSB_chi2+2*TMath::Log(FreePoints);

					TF1* fLinLSB_m = new TF1("fLinLSB_m","pol1",pTminPlot,pTmaxPlot);
					fLinLSB_m->SetParameter(0,LinLSB_p0+err_LinLSB_p0);
					fLinLSB_m->SetParameter(1,LinLSB_p1-err_LinLSB_p1);
					fLinLSB_m->SetLineWidth(0.2);
					fLinLSB_m->SetLineColor(kRed);
					fLinLSB_m->SetLineStyle(2);
					fLinLSB_m->Draw("same");
					TF1* fLinLSB_p = new TF1("fLinLSB_p","pol1",pTminPlot,pTmaxPlot);
					fLinLSB_p->SetParameter(0,LinLSB_p0-err_LinLSB_p0);
					fLinLSB_p->SetParameter(1,LinLSB_p1+err_LinLSB_p1);
					fLinLSB_p->SetLineWidth(0.2);
					fLinLSB_p->SetLineColor(kRed);
					fLinLSB_p->SetLineStyle(2);
					fLinLSB_p->Draw("same");

					if(DrawSBLatex){
						sprintf(text,"#color[2]{Fitting pol1, LSB:} p0 = %1.4f #pm %1.4f, p1 = %1.4f #pm %1.4f, #chi^{2}/ndf = %1.4f / %d, BIC = %1.3f", LinLSB_p0, err_LinLSB_p0, LinLSB_p1, err_LinLSB_p1,LinLSB_chi2,LinLSB_NDF,LinLSB_BIC);
						TLatex *textLinLSB = new TLatex(pTminPlot*1.2,yMin+(yMax-yMin)*0.95,text);
						textLinLSB->SetTextSize(FontSize)                                                                                                                                                                                                                                             ;
						textLinLSB->Draw( "same" )                                                                                                                                                                                                                                                 ;
					}



					TF1* fLinRSB = new TF1("fLinRSB","pol1",pTminPlot,pTmaxPlot);
					RSBgraph->Fit("fLinRSB",FitOptions);
					fLinRSB->SetLineWidth(0.4);
					fLinRSB->SetLineColor(kBlue);
					fLinRSB->SetLineStyle(1);
					fLinRSB->Draw("same");

					double LinRSB_p0 = fLinRSB->GetParameter(0);
					double err_LinRSB_p0 = fLinRSB->GetParError(0);
					double LinRSB_p1 = fLinRSB->GetParameter(1);
					double err_LinRSB_p1 = fLinRSB->GetParError(1);
					double LinRSB_chi2=fLinRSB->GetChisquare();
					int LinRSB_NDF=fLinRSB->GetNDF();
					double LinRSB_BIC=LinRSB_chi2+2*TMath::Log(FreePoints);

					TF1* fLinRSB_m = new TF1("fLinRSB_m","pol1",pTminPlot,pTmaxPlot);
					fLinRSB_m->SetParameter(0,LinRSB_p0+err_LinRSB_p0);
					fLinRSB_m->SetParameter(1,LinRSB_p1-err_LinRSB_p1);
					fLinRSB_m->SetLineWidth(0.2);
					fLinRSB_m->SetLineColor(kBlue);
					fLinRSB_m->SetLineStyle(2);
					fLinRSB_m->Draw("same");
					TF1* fLinRSB_p = new TF1("fLinRSB_p","pol1",pTminPlot,pTmaxPlot);
					fLinRSB_p->SetParameter(0,LinRSB_p0-err_LinRSB_p0);
					fLinRSB_p->SetParameter(1,LinRSB_p1+err_LinRSB_p1);
					fLinRSB_p->SetLineWidth(0.2);
					fLinRSB_p->SetLineColor(kBlue);
					fLinRSB_p->SetLineStyle(2);
					fLinRSB_p->Draw("same");

					if(DrawSBLatex){
						sprintf(text,"#color[4]{Fitting pol1, RSB:} p0 = %1.4f #pm %1.4f, p1 = %1.4f #pm %1.4f, #chi^{2}/ndf = %1.4f / %d, BIC = %1.3f", LinRSB_p0, err_LinRSB_p0, LinRSB_p1, err_LinRSB_p1,LinRSB_chi2,LinRSB_NDF,LinRSB_BIC);
						TLatex *textLinRSB = new TLatex(pTminPlot*1.2,yMin+(yMax-yMin)*0.9,text);
						textLinRSB->SetTextSize(FontSize)                                                                                                                                                                                                                                             ;
						textLinRSB->Draw( "same" )                                                                                                                                                                                                                                                 ;
					}
				}

				/*	      if(Const_BIC<Lin_BIC) sprintf(text,"BIC prefers #color[4]{constant} fit");
									else sprintf(text,"BIC prefers #color[2]{linear} fit");
									TLatex *textConst3 = new TLatex(pTminPlot*1.2,yMin+(yMax-yMin)*0.8,text);
									textConst3->SetTextSize(1.75*FontSize)                                                                                                                                                                                                                                             ;
									textConst3->Draw( "same" )                                                                                                                                                                                                                                                 ;
									*/


				TLine* extreme0 = new TLine( PlotpTMin-DeltaXminOVERALL, 0, PlotpTMax ,0);
				extreme0->SetLineWidth( 2 );
				extreme0->SetLineStyle( 1 );
				extreme0->SetLineColor( kBlack );
				extreme0->Draw( "same" );
			} //SBmSigPlots



			if(BGratioFits){ //not yet changed 12.12.2012 Linlin
				TGraphAsymmErrors *Ups1SBGratiograph;
				TGraphAsymmErrors *Ups2SBGratiograph;
				TGraphAsymmErrors *Ups3SBGratiograph;

				bool CombineRapsForSB=true;

				double FreePoints=5.;
				if(CombineRapsForSB) FreePoints=10.;

				if(!CombineRapsForSB){
					Ups1SBGratiograph = graphCompareFile1;
					Ups2SBGratiograph = graphCompareFile2;
					Ups3SBGratiograph = graphCompareFile3;
				}


				if(CombineRapsForSB){

					TGraphAsymmErrors* graphCompareFile1_OtherRap = (TGraphAsymmErrors*) CompareFile1->Get(GraphNameOtherRap);
					TGraphAsymmErrors* graphCompareFile2_OtherRap = (TGraphAsymmErrors*) CompareFile2->Get(GraphNameOtherRap);
					TGraphAsymmErrors* graphCompareFile3_OtherRap = (TGraphAsymmErrors*) CompareFile3->Get(GraphNameOtherRap);

					const int nBinspTSB=14;

					double ptCentre_SB_Ups1SBGratio[nBinspTSB];
					double lmean_SB_Ups1SBGratio[nBinspTSB];
					double fit_lmean_errlow_SB_Ups1SBGratio[nBinspTSB];
					double fit_lmean_errhigh_SB_Ups1SBGratio[nBinspTSB];

					double ptCentre_SB_Ups2SBGratio[nBinspTSB];
					double lmean_SB_Ups2SBGratio[nBinspTSB];
					double fit_lmean_errlow_SB_Ups2SBGratio[nBinspTSB];
					double fit_lmean_errhigh_SB_Ups2SBGratio[nBinspTSB];

					double ptCentre_SB_Ups3SBGratio[nBinspTSB];
					double lmean_SB_Ups3SBGratio[nBinspTSB];
					double fit_lmean_errlow_SB_Ups3SBGratio[nBinspTSB];
					double fit_lmean_errhigh_SB_Ups3SBGratio[nBinspTSB];

					int runningSB=0;
					for(int rapBinSB=1;rapBinSB<3;rapBinSB++){
						for(int ptBinSB=6;ptBinSB<11;ptBinSB++){

							if(rapBinSB==1){
								graphCompareFile1_OtherRap->GetPoint(ptBinSB-1,ptCentre_SB_Ups1SBGratio[runningSB],lmean_SB_Ups1SBGratio[runningSB]);
								fit_lmean_errhigh_SB_Ups1SBGratio[runningSB]=graphCompareFile1_OtherRap->GetErrorYhigh(ptBinSB-1);
								fit_lmean_errlow_SB_Ups1SBGratio[runningSB]=graphCompareFile1_OtherRap->GetErrorYlow(ptBinSB-1);
							}
							if(rapBinSB==2){
								graphCompareFile1->GetPoint(ptBinSB-1,ptCentre_SB_Ups1SBGratio[runningSB],lmean_SB_Ups1SBGratio[runningSB]);
								fit_lmean_errhigh_SB_Ups1SBGratio[runningSB]=graphCompareFile1->GetErrorYhigh(ptBinSB-1);
								fit_lmean_errlow_SB_Ups1SBGratio[runningSB]=graphCompareFile1->GetErrorYlow(ptBinSB-1);
							}

							if(rapBinSB==1){
								graphCompareFile2_OtherRap->GetPoint(ptBinSB-1,ptCentre_SB_Ups2SBGratio[runningSB],lmean_SB_Ups2SBGratio[runningSB]);
								fit_lmean_errhigh_SB_Ups2SBGratio[runningSB]=graphCompareFile2_OtherRap->GetErrorYhigh(ptBinSB-1);
								fit_lmean_errlow_SB_Ups2SBGratio[runningSB]=graphCompareFile2_OtherRap->GetErrorYlow(ptBinSB-1);
							}
							if(rapBinSB==2){
								graphCompareFile2->GetPoint(ptBinSB-1,ptCentre_SB_Ups2SBGratio[runningSB],lmean_SB_Ups2SBGratio[runningSB]);
								fit_lmean_errhigh_SB_Ups2SBGratio[runningSB]=graphCompareFile2->GetErrorYhigh(ptBinSB-1);
								fit_lmean_errlow_SB_Ups2SBGratio[runningSB]=graphCompareFile2->GetErrorYlow(ptBinSB-1);
							}

							if(rapBinSB==1){
								graphCompareFile3_OtherRap->GetPoint(ptBinSB-1,ptCentre_SB_Ups3SBGratio[runningSB],lmean_SB_Ups3SBGratio[runningSB]);
								fit_lmean_errhigh_SB_Ups3SBGratio[runningSB]=graphCompareFile3_OtherRap->GetErrorYhigh(ptBinSB-1);
								fit_lmean_errlow_SB_Ups3SBGratio[runningSB]=graphCompareFile3_OtherRap->GetErrorYlow(ptBinSB-1);
							}
							if(rapBinSB==2){
								graphCompareFile3->GetPoint(ptBinSB-1,ptCentre_SB_Ups3SBGratio[runningSB],lmean_SB_Ups3SBGratio[runningSB]);
								fit_lmean_errhigh_SB_Ups3SBGratio[runningSB]=graphCompareFile3->GetErrorYhigh(ptBinSB-1);
								fit_lmean_errlow_SB_Ups3SBGratio[runningSB]=graphCompareFile3->GetErrorYlow(ptBinSB-1);
							}

							runningSB++;
						}
					}


					Ups1SBGratiograph = new TGraphAsymmErrors(nBinspTSB,ptCentre_SB_Ups1SBGratio,lmean_SB_Ups1SBGratio,0,0,fit_lmean_errlow_SB_Ups1SBGratio,fit_lmean_errhigh_SB_Ups1SBGratio);
					Ups2SBGratiograph = new TGraphAsymmErrors(nBinspTSB,ptCentre_SB_Ups2SBGratio,lmean_SB_Ups2SBGratio,0,0,fit_lmean_errlow_SB_Ups2SBGratio,fit_lmean_errhigh_SB_Ups2SBGratio);
					Ups3SBGratiograph = new TGraphAsymmErrors(nBinspTSB,ptCentre_SB_Ups3SBGratio,lmean_SB_Ups3SBGratio,0,0,fit_lmean_errlow_SB_Ups3SBGratio,fit_lmean_errhigh_SB_Ups3SBGratio);

					//			char savequick[200];
					//			sprintf(savequick,"tmp/Ups1SBGratiograph_rap%d.root",rapBin);
					//			Ups1SBGratiograph->SaveAs(savequick);

				}

				TLine* extreme0 = new TLine( PlotpTMin-DeltaXminOVERALL, 0, PlotpTMax ,0);
				extreme0->SetLineWidth( 2 );
				extreme0->SetLineStyle( 1 );
				extreme0->SetLineColor( kBlack );
				extreme0->Draw( "same" );

				double pTmaxPlot=72.; //50;
				double pTminPlot=6.;  //10;

				double FontSize=0.0215;
				double xText=15;
				char text[200];

				char FitOptions[200];
				sprintf(FitOptions,"EFNR");


				bool DrawSBLatex=false;
				bool MonotonicityConst=true;
				bool MonotonicityLin=false;


				if(MonotonicityConst){


					TF1* fConstUps1SBGratio = new TF1("fConstUps1SBGratio","pol0",pTminPlot,pTmaxPlot);
					Ups1SBGratiograph->Fit("fConstUps1SBGratio",FitOptions);
					fConstUps1SBGratio->SetLineWidth(0.4);
					fConstUps1SBGratio->SetLineStyle(1);
					fConstUps1SBGratio->SetLineColor(kRed);
					fConstUps1SBGratio->Draw("same");

					double ConstUps1SBGratio_p0 = fConstUps1SBGratio->GetParameter(0);
					double err_ConstUps1SBGratio_p0 = fConstUps1SBGratio->GetParError(0);
					double ConstUps1SBGratio_chi2=fConstUps1SBGratio->GetChisquare();
					int ConstUps1SBGratio_NDF=fConstUps1SBGratio->GetNDF();
					double ConstUps1SBGratio_BIC=ConstUps1SBGratio_chi2+1*TMath::Log(FreePoints);

					TF1* fConstUps1SBGratio_m = new TF1("fConstUps1SBGratio_m","pol0",pTminPlot,pTmaxPlot);
					fConstUps1SBGratio_m->SetParameter(0,ConstUps1SBGratio_p0+err_ConstUps1SBGratio_p0);
					fConstUps1SBGratio_m->SetLineWidth(0.2);
					fConstUps1SBGratio_m->SetLineColor(kRed);
					fConstUps1SBGratio_m->SetLineStyle(2);
					fConstUps1SBGratio_m->Draw("same");
					TF1* fConstUps1SBGratio_p = new TF1("fConstUps1SBGratio_p","pol0",pTminPlot,pTmaxPlot);
					fConstUps1SBGratio_p->SetParameter(0,ConstUps1SBGratio_p0-err_ConstUps1SBGratio_p0);
					fConstUps1SBGratio_p->SetLineWidth(0.2);
					fConstUps1SBGratio_p->SetLineColor(kRed);
					fConstUps1SBGratio_p->SetLineStyle(2);
					fConstUps1SBGratio_p->Draw("same");

					if(DrawSBLatex){
						sprintf(text,"#color[2]{Fitting constant, Ups1SBGratio:} const. = %1.4f #pm %1.4f, #chi^{2}/ndf = %1.4f / %d, BIC = %1.3f", ConstUps1SBGratio_p0, err_ConstUps1SBGratio_p0, ConstUps1SBGratio_chi2, ConstUps1SBGratio_NDF,ConstUps1SBGratio_BIC);
						TLatex *textConstUps1SBGratio = new TLatex(pTminPlot*1.2,yMin+(yMax-yMin)*0.95,text);
						textConstUps1SBGratio->SetTextSize(FontSize)                                                                                                                                                                                                                                             ;
						textConstUps1SBGratio->Draw( "same" )                                                                                                                                                                                                                                                 ;
					}


					TF1* fConstUps2SBGratio = new TF1("fConstUps2SBGratio","pol0",pTminPlot,pTmaxPlot);
					Ups2SBGratiograph->Fit("fConstUps2SBGratio",FitOptions);
					fConstUps2SBGratio->SetLineWidth(0.4);
					fConstUps2SBGratio->SetLineStyle(1);
					fConstUps2SBGratio->SetLineColor(kBlue);
					if(nComp>2) fConstUps2SBGratio->Draw("same");

					double ConstUps2SBGratio_p0 = fConstUps2SBGratio->GetParameter(0);
					double err_ConstUps2SBGratio_p0 = fConstUps2SBGratio->GetParError(0);
					double ConstUps2SBGratio_chi2=fConstUps2SBGratio->GetChisquare();
					int ConstUps2SBGratio_NDF=fConstUps2SBGratio->GetNDF();
					double ConstUps2SBGratio_BIC=ConstUps2SBGratio_chi2+1*TMath::Log(FreePoints);

					TF1* fConstUps2SBGratio_m = new TF1("fConstUps2SBGratio_m","pol0",pTminPlot,pTmaxPlot);
					fConstUps2SBGratio_m->SetParameter(0,ConstUps2SBGratio_p0+err_ConstUps2SBGratio_p0);
					fConstUps2SBGratio_m->SetLineWidth(0.2);
					fConstUps2SBGratio_m->SetLineColor(kBlue);
					fConstUps2SBGratio_m->SetLineStyle(2);
					if(nComp>2) fConstUps2SBGratio_m->Draw("same");
					TF1* fConstUps2SBGratio_p = new TF1("fConstUps2SBGratio_p","pol0",pTminPlot,pTmaxPlot);
					fConstUps2SBGratio_p->SetParameter(0,ConstUps2SBGratio_p0-err_ConstUps2SBGratio_p0);
					fConstUps2SBGratio_p->SetLineWidth(0.2);
					fConstUps2SBGratio_p->SetLineColor(kBlue);
					fConstUps2SBGratio_p->SetLineStyle(2);
					if(nComp>2) fConstUps2SBGratio_p->Draw("same");

					if(DrawSBLatex){
						sprintf(text,"#color[4]{Fitting constant, Ups2SBGratio:} const. = %1.4f #pm %1.4f, #chi^{2}/ndf = %1.4f / %d, BIC = %1.3f", ConstUps2SBGratio_p0, err_ConstUps2SBGratio_p0, ConstUps2SBGratio_chi2, ConstUps2SBGratio_NDF,ConstUps2SBGratio_BIC);
						TLatex *textConstUps2SBGratio = new TLatex(pTminPlot*1.2,yMin+(yMax-yMin)*0.9,text);
						textConstUps2SBGratio->SetTextSize(FontSize)                                                                                                                                                                                                                                             ;
						if(nComp>2) textConstUps2SBGratio->Draw( "same" )                                                                                                                                                                                                                                                 ;
					}


					TF1* fConstUps3SBGratio = new TF1("fConstUps3SBGratio","pol0",pTminPlot,pTmaxPlot);
					Ups3SBGratiograph->Fit("fConstUps3SBGratio",FitOptions);
					fConstUps3SBGratio->SetLineWidth(0.4);
					fConstUps3SBGratio->SetLineStyle(1);
					fConstUps3SBGratio->SetLineColor(kGreen-2);
					if(nComp>2) fConstUps3SBGratio->Draw("same");

					double ConstUps3SBGratio_p0 = fConstUps3SBGratio->GetParameter(0);
					double err_ConstUps3SBGratio_p0 = fConstUps3SBGratio->GetParError(0);
					double ConstUps3SBGratio_chi2=fConstUps3SBGratio->GetChisquare();
					int ConstUps3SBGratio_NDF=fConstUps3SBGratio->GetNDF();
					double ConstUps3SBGratio_BIC=ConstUps3SBGratio_chi2+1*TMath::Log(FreePoints);

					TF1* fConstUps3SBGratio_m = new TF1("fConstUps3SBGratio_m","pol0",pTminPlot,pTmaxPlot);
					fConstUps3SBGratio_m->SetParameter(0,ConstUps3SBGratio_p0+err_ConstUps3SBGratio_p0);
					fConstUps3SBGratio_m->SetLineWidth(0.2);
					fConstUps3SBGratio_m->SetLineColor(kGreen-2);
					fConstUps3SBGratio_m->SetLineStyle(2);
					if(nComp>2) fConstUps3SBGratio_m->Draw("same");
					TF1* fConstUps3SBGratio_p = new TF1("fConstUps3SBGratio_p","pol0",pTminPlot,pTmaxPlot);
					fConstUps3SBGratio_p->SetParameter(0,ConstUps3SBGratio_p0-err_ConstUps3SBGratio_p0);
					fConstUps3SBGratio_p->SetLineWidth(0.2);
					fConstUps3SBGratio_p->SetLineColor(kGreen-2);
					fConstUps3SBGratio_p->SetLineStyle(2);
					if(nComp>2) fConstUps3SBGratio_p->Draw("same");

					if(DrawSBLatex){
						sprintf(text,"#color[4]{Fitting constant, Ups3SBGratio:} const. = %1.4f #pm %1.4f, #chi^{2}/ndf = %1.4f / %d, BIC = %1.3f", ConstUps3SBGratio_p0, err_ConstUps3SBGratio_p0, ConstUps3SBGratio_chi2, ConstUps3SBGratio_NDF,ConstUps3SBGratio_BIC);
						TLatex *textConstUps3SBGratio = new TLatex(pTminPlot*1.2,yMin+(yMax-yMin)*0.85,text);
						textConstUps3SBGratio->SetTextSize(FontSize)                                                                                                                                                                                                                                             ;
						if(nComp>2) textConstUps3SBGratio->Draw( "same" )                                                                                                                                                                                                                                                 ;
					}

					bool Ups1SBGratioCompZero=false;
					bool Ups2SBGratioCompZero=false;
					bool Ups3SBGratioCompZero=false;

					if(TMath::Abs(ConstUps1SBGratio_p0/err_ConstUps1SBGratio_p0)<1) Ups1SBGratioCompZero=true;
					if(TMath::Abs(ConstUps2SBGratio_p0/err_ConstUps2SBGratio_p0)<1) Ups2SBGratioCompZero=true;
					if(TMath::Abs(ConstUps3SBGratio_p0/err_ConstUps3SBGratio_p0)<1) Ups2SBGratioCompZero=true;

					if(nComp>2)  sprintf(text,"Polarization of BG ratios compatible with zero:");
					if(nComp<2&&Ups1SBGratioCompZero)  sprintf(text,"Polarization of BG ratios compatible with zero");
					if(nComp<2&&!Ups1SBGratioCompZero)  sprintf(text,"Polarization of BG ratios not compatible with zero");
					TLatex *textMonotonicityResult = new TLatex(pTminPlot*1.2,yMin+(yMax-yMin)*0.9,text);
					textMonotonicityResult->SetTextSize(1.75*FontSize)                                                                                                                                                                                                                                             ;
					textMonotonicityResult->Draw( "same" )                                                                                                                                                                                                                                                 ;

					if(Ups1SBGratioCompZero&&Ups2SBGratioCompZero&&Ups3SBGratioCompZero) sprintf(text,"#psi(1S), #psi(2S) and #psi(3S) mass regions");
					if(Ups1SBGratioCompZero&&Ups2SBGratioCompZero&&!Ups3SBGratioCompZero) sprintf(text,"#psi(1S) and #psi(2S) mass regions");
					if(Ups1SBGratioCompZero&&!Ups2SBGratioCompZero&&Ups3SBGratioCompZero) sprintf(text,"#psi(1S) and #psi(3S) mass regions");
					if(!Ups1SBGratioCompZero&&Ups2SBGratioCompZero&&Ups3SBGratioCompZero) sprintf(text,"#psi(2S) and #psi(3S) mass regions");
					if(Ups1SBGratioCompZero&&!Ups2SBGratioCompZero&&!Ups3SBGratioCompZero) sprintf(text,"#psi(1S) mass region");
					if(!Ups1SBGratioCompZero&&Ups2SBGratioCompZero&&!Ups3SBGratioCompZero) sprintf(text,"#psi(2S) mass region");
					if(!Ups1SBGratioCompZero&&!Ups2SBGratioCompZero&&Ups3SBGratioCompZero) sprintf(text,"#psi(3S) mass region");
					if(!Ups1SBGratioCompZero&&!Ups2SBGratioCompZero&&!Ups3SBGratioCompZero) sprintf(text,"none");

					TLatex *textMonotonicityResult2 = new TLatex(pTminPlot*1.2,yMin+(yMax-yMin)*0.825,text);
					textMonotonicityResult2->SetTextSize(1.75*FontSize)                                                                                                                                                                                                                                             ;
					if(nComp>2) textMonotonicityResult2->Draw( "same" )                                                                                                                                                                                                                                                 ;

				}

				///////////////////////////////////////////////////

				if(MonotonicityLin){


					TF1* fLinUps1SBGratio = new TF1("fLinUps1SBGratio","pol1",pTminPlot,pTmaxPlot);
					Ups1SBGratiograph->Fit("fLinUps1SBGratio",FitOptions);
					fLinUps1SBGratio->SetLineWidth(0.4);
					fLinUps1SBGratio->SetLineColor(kRed);
					fLinUps1SBGratio->SetLineStyle(1);
					fLinUps1SBGratio->Draw("same");

					double LinUps1SBGratio_p0 = fLinUps1SBGratio->GetParameter(0);
					double err_LinUps1SBGratio_p0 = fLinUps1SBGratio->GetParError(0);
					double LinUps1SBGratio_p1 = fLinUps1SBGratio->GetParameter(1);
					double err_LinUps1SBGratio_p1 = fLinUps1SBGratio->GetParError(1);
					double LinUps1SBGratio_chi2=fLinUps1SBGratio->GetChisquare();
					int LinUps1SBGratio_NDF=fLinUps1SBGratio->GetNDF();
					double LinUps1SBGratio_BIC=LinUps1SBGratio_chi2+2*TMath::Log(FreePoints);

					TF1* fLinUps1SBGratio_m = new TF1("fLinUps1SBGratio_m","pol1",pTminPlot,pTmaxPlot);
					fLinUps1SBGratio_m->SetParameter(0,LinUps1SBGratio_p0+err_LinUps1SBGratio_p0);
					fLinUps1SBGratio_m->SetParameter(1,LinUps1SBGratio_p1-err_LinUps1SBGratio_p1);
					fLinUps1SBGratio_m->SetLineWidth(0.2);
					fLinUps1SBGratio_m->SetLineColor(kRed);
					fLinUps1SBGratio_m->SetLineStyle(2);
					fLinUps1SBGratio_m->Draw("same");
					TF1* fLinUps1SBGratio_p = new TF1("fLinUps1SBGratio_p","pol1",pTminPlot,pTmaxPlot);
					fLinUps1SBGratio_p->SetParameter(0,LinUps1SBGratio_p0-err_LinUps1SBGratio_p0);
					fLinUps1SBGratio_p->SetParameter(1,LinUps1SBGratio_p1+err_LinUps1SBGratio_p1);
					fLinUps1SBGratio_p->SetLineWidth(0.2);
					fLinUps1SBGratio_p->SetLineColor(kRed);
					fLinUps1SBGratio_p->SetLineStyle(2);
					fLinUps1SBGratio_p->Draw("same");

					if(DrawSBLatex){
						sprintf(text,"#color[2]{Fitting pol1, Ups1SBGratio:} p0 = %1.4f #pm %1.4f, p1 = %1.4f #pm %1.4f, #chi^{2}/ndf = %1.4f / %d, BIC = %1.3f", LinUps1SBGratio_p0, err_LinUps1SBGratio_p0, LinUps1SBGratio_p1, err_LinUps1SBGratio_p1,LinUps1SBGratio_chi2,LinUps1SBGratio_NDF,LinUps1SBGratio_BIC);
						TLatex *textLinUps1SBGratio = new TLatex(pTminPlot*1.2,yMin+(yMax-yMin)*0.95,text);
						textLinUps1SBGratio->SetTextSize(FontSize)                                                                                                                                                                                                                                             ;
						textLinUps1SBGratio->Draw( "same" )                                                                                                                                                                                                                                                 ;
					}



					TF1* fLinUps2SBGratio = new TF1("fLinUps2SBGratio","pol1",pTminPlot,pTmaxPlot);
					Ups2SBGratiograph->Fit("fLinUps2SBGratio",FitOptions);
					fLinUps2SBGratio->SetLineWidth(0.4);
					fLinUps2SBGratio->SetLineColor(kBlue);
					fLinUps2SBGratio->SetLineStyle(1);
					fLinUps2SBGratio->Draw("same");

					double LinUps2SBGratio_p0 = fLinUps2SBGratio->GetParameter(0);
					double err_LinUps2SBGratio_p0 = fLinUps2SBGratio->GetParError(0);
					double LinUps2SBGratio_p1 = fLinUps2SBGratio->GetParameter(1);
					double err_LinUps2SBGratio_p1 = fLinUps2SBGratio->GetParError(1);
					double LinUps2SBGratio_chi2=fLinUps2SBGratio->GetChisquare();
					int LinUps2SBGratio_NDF=fLinUps2SBGratio->GetNDF();
					double LinUps2SBGratio_BIC=LinUps2SBGratio_chi2+2*TMath::Log(FreePoints);

					TF1* fLinUps2SBGratio_m = new TF1("fLinUps2SBGratio_m","pol1",pTminPlot,pTmaxPlot);
					fLinUps2SBGratio_m->SetParameter(0,LinUps2SBGratio_p0+err_LinUps2SBGratio_p0);
					fLinUps2SBGratio_m->SetParameter(1,LinUps2SBGratio_p1-err_LinUps2SBGratio_p1);
					fLinUps2SBGratio_m->SetLineWidth(0.2);
					fLinUps2SBGratio_m->SetLineColor(kBlue);
					fLinUps2SBGratio_m->SetLineStyle(2);
					fLinUps2SBGratio_m->Draw("same");
					TF1* fLinUps2SBGratio_p = new TF1("fLinUps2SBGratio_p","pol1",pTminPlot,pTmaxPlot);
					fLinUps2SBGratio_p->SetParameter(0,LinUps2SBGratio_p0-err_LinUps2SBGratio_p0);
					fLinUps2SBGratio_p->SetParameter(1,LinUps2SBGratio_p1+err_LinUps2SBGratio_p1);
					fLinUps2SBGratio_p->SetLineWidth(0.2);
					fLinUps2SBGratio_p->SetLineColor(kBlue);
					fLinUps2SBGratio_p->SetLineStyle(2);
					fLinUps2SBGratio_p->Draw("same");

					if(DrawSBLatex){
						sprintf(text,"#color[4]{Fitting pol1, Ups2SBGratio:} p0 = %1.4f #pm %1.4f, p1 = %1.4f #pm %1.4f, #chi^{2}/ndf = %1.4f / %d, BIC = %1.3f", LinUps2SBGratio_p0, err_LinUps2SBGratio_p0, LinUps2SBGratio_p1, err_LinUps2SBGratio_p1,LinUps2SBGratio_chi2,LinUps2SBGratio_NDF,LinUps2SBGratio_BIC);
						TLatex *textLinUps2SBGratio = new TLatex(pTminPlot*1.2,yMin+(yMax-yMin)*0.9,text);
						textLinUps2SBGratio->SetTextSize(FontSize)                                                                                                                                                                                                                                             ;
						textLinUps2SBGratio->Draw( "same" )                                                                                                                                                                                                                                                 ;
					}
				}

			} //BGratioFits


			if(PlotCompare){ //not yet changed 12.12.2012 Linlin
				////plot noRho
				//graphCompareFile1->RemovePoint(0);
				//graphCompareFile1->RemovePoint(0);
				//graphCompareFile1->RemovePoint(0);
				//graphCompareFile1->RemovePoint(0);
				//graphCompareFile1->RemovePoint(0);
				//graphCompareFile1->RemovePoint(0);
				//graphCompareFile1->RemovePoint(0);
				//graphCompareFile1->RemovePoint(0);
				//graphCompareFile1->RemovePoint(0);

				graphCompareFile1->SetMarkerColor(kRed);
				graphCompareFile1->SetLineColor(kRed);
				graphCompareFile1->SetMarkerStyle(28);
				//graphCompareFile1->SetLineWidth(3);
				graphCompareFile2->SetMarkerColor(kBlue);
				graphCompareFile2->SetLineColor(kBlue);
				graphCompareFile2->SetMarkerStyle(30);
				graphCompareFile3->SetMarkerColor(kGreen-2);
				graphCompareFile3->SetLineColor(kGreen-2);
				graphCompareFile3->SetMarkerStyle(3);
				graphCompareFile4->SetMarkerColor(kOrange);
				graphCompareFile4->SetLineColor(kOrange);
				graphCompareFile4->SetMarkerStyle(5);
				if(PlotBG0plots||SetCompStyle){
					graphCompareFile1->SetMarkerStyle(24);
					graphCompareFile2->SetMarkerStyle(24);
					graphCompareFile3->SetMarkerStyle(24);
					graphCompareFile4->SetMarkerStyle(24);
					graphCompareFile1->SetMarkerSize(BG0MarkerSize);
					graphCompareFile2->SetMarkerSize(BG0MarkerSize);
					graphCompareFile3->SetMarkerSize(BG0MarkerSize);
					graphCompareFile4->SetMarkerSize(BG0MarkerSize);

				}

				TLine* extreme0 = new TLine( PlotpTMin-DeltaXminOVERALL, 0, PlotpTMax ,0);
				extreme0->SetLineWidth( 1 );
				extreme0->SetLineStyle( 2 );
				extreme0->SetLineColor( kBlack );
				extreme0->Draw( "same" );


				int nBinspT=ptBinMax-ptBinMin+1;
				double ptCentre[nBinspT];
				double ptCentreErr_low[nBinspT];
				double ptCentreErr_high[nBinspT];
				double lmean[nBinspT];
				double lmean_errlow[nBinspT];
				double lmean_errhigh[nBinspT];

				int pt=0;

				if(nComp>0){//FindLegend
					pt=0;
					for(int ptBin = ptBinMin; ptBin < ptBinMax+1; ptBin++) {
						graphCompareFile1->GetPoint(ptBin-1,ptCentre[pt],lmean[pt]);
						ptCentreErr_high[pt]=graphCompareFile1->GetErrorXhigh(ptBin-1);
						ptCentreErr_low[pt]=graphCompareFile1->GetErrorXlow(ptBin-1);
						lmean_errhigh[pt]=graphCompareFile1->GetErrorYhigh(ptBin-1);
						lmean_errlow[pt]=graphCompareFile1->GetErrorYlow(ptBin-1);
						if(ShiftCompareInX){
							ptCentre[pt] = ptCentre[pt] - DeltaXCompare;
							ptCentreErr_high[pt] = ptCentreErr_high[pt]+DeltaXCompare; //0.
							ptCentreErr_low[pt] = ptCentreErr_low[pt]-DeltaXCompare; //0.
						}
						pt++;
					}
					graphCompareFile1 = new TGraphAsymmErrors(nBinspT,ptCentre,lmean,ptCentreErr_low,ptCentreErr_high,lmean_errlow,lmean_errhigh);
					graphCompareFile1->SetMarkerColor(kRed);
					graphCompareFile1->SetLineColor(kRed);
					graphCompareFile1->SetMarkerStyle(28);
					graphCompareFile1->SetMarkerSize(2.75);

					graphCompareFile1->Draw(drawGraphStyle);
					sprintf(complegendentry,"#lambda(LSB)-#lambda(sig. region)");
					sprintf(complegendentry,"%s",LegendEntryCompID1);
					plotcompLegend->AddEntry(graphCompareFile1,complegendentry,"lp");
				}
				if(nComp>1){
					pt=0;
					for(int ptBin = ptBinMin; ptBin < ptBinMax+1; ptBin++) {
						graphCompareFile2->GetPoint(ptBin-1,ptCentre[pt],lmean[pt]);
						ptCentreErr_high[pt]=graphCompareFile2->GetErrorXhigh(ptBin-1);
						ptCentreErr_low[pt]=graphCompareFile2->GetErrorXlow(ptBin-1);
						lmean_errhigh[pt]=graphCompareFile2->GetErrorYhigh(ptBin-1);
						lmean_errlow[pt]=graphCompareFile2->GetErrorYlow(ptBin-1);
						if(ShiftCompareInX){
							ptCentre[pt] = ptCentre[pt] + DeltaXCompare;
							ptCentreErr_high[pt] = ptCentreErr_high[pt]-DeltaXCompare; //0.
							ptCentreErr_low[pt] = ptCentreErr_low[pt]+DeltaXCompare; //0.
						} 
						pt++;
					}
					graphCompareFile2 = new TGraphAsymmErrors(nBinspT,ptCentre,lmean,ptCentreErr_low,ptCentreErr_high,lmean_errlow,lmean_errhigh);
					graphCompareFile2->SetMarkerColor(kBlue);
					graphCompareFile2->SetLineColor(kBlue);
					graphCompareFile2->SetMarkerStyle(30);
					graphCompareFile2->SetMarkerSize(2.75);

					graphCompareFile2->Draw(drawGraphStyle);
					sprintf(complegendentry,"#lambda(RSB)-#lambda(sig. region)");
					sprintf(complegendentry,"%s",LegendEntryCompID2);
					plotcompLegend->AddEntry(graphCompareFile2,complegendentry,"lp");
				}
				if(nComp>2){
					pt=0;
					for(int ptBin = ptBinMin; ptBin < ptBinMax+1; ptBin++) {
						graphCompareFile3->GetPoint(ptBin-1,ptCentre[pt],lmean[pt]);
						ptCentreErr_high[pt]=graphCompareFile3->GetErrorXhigh(ptBin-1);
						ptCentreErr_low[pt]=graphCompareFile3->GetErrorXlow(ptBin-1);
						lmean_errhigh[pt]=graphCompareFile3->GetErrorYhigh(ptBin-1);
						lmean_errlow[pt]=graphCompareFile3->GetErrorYlow(ptBin-1);
						if(ShiftCompareInX){
							ptCentre[pt] = ptCentre[pt] - 2*DeltaXCompare;
							ptCentreErr_high[pt] = ptCentreErr_high[pt]+2*DeltaXCompare; //0.
							ptCentreErr_low[pt] = ptCentreErr_low[pt]-2*DeltaXCompare; //0.
						}
						pt++;
					}
					graphCompareFile3 = new TGraphAsymmErrors(nBinspT,ptCentre,lmean,ptCentreErr_low,ptCentreErr_high,lmean_errlow,lmean_errhigh);
					graphCompareFile3->SetMarkerColor(kGreen-2);
					graphCompareFile3->SetLineColor(kGreen-2);
					graphCompareFile3->SetMarkerStyle(3);
					graphCompareFile3->SetMarkerSize(2.75);

					graphCompareFile3->Draw(drawGraphStyle);
					sprintf(complegendentry,"Right sided 1.5 sigma");
					sprintf(complegendentry,"%s",LegendEntryCompID3);
					plotcompLegend->AddEntry(graphCompareFile3,complegendentry,"lp");
				}
				if(nComp>3){
					pt=0;
					for(int ptBin = ptBinMin; ptBin < ptBinMax+1; ptBin++) {
						graphCompareFile4->GetPoint(ptBin-1,ptCentre[pt],lmean[pt]);
						ptCentreErr_high[pt]=graphCompareFile4->GetErrorXhigh(ptBin-1);
						ptCentreErr_low[pt]=graphCompareFile4->GetErrorXlow(ptBin-1);
						lmean_errhigh[pt]=graphCompareFile4->GetErrorYhigh(ptBin-1);
						lmean_errlow[pt]=graphCompareFile4->GetErrorYlow(ptBin-1);
						if(ShiftCompareInX){
							ptCentre[pt] = ptCentre[pt] + 2*DeltaXCompare;
							ptCentreErr_high[pt] = ptCentreErr_high[pt]-2*DeltaXCompare; //0.
							ptCentreErr_low[pt] = ptCentreErr_low[pt]+2*DeltaXCompare; //0.
						} 
						pt++;
					}
					graphCompareFile4 = new TGraphAsymmErrors(nBinspT,ptCentre,lmean,ptCentreErr_low,ptCentreErr_high,lmean_errlow,lmean_errhigh);
					graphCompareFile4->SetMarkerColor(kOrange);
					graphCompareFile4->SetLineColor(kOrange);
					graphCompareFile4->SetMarkerStyle(5);
					graphCompareFile4->SetMarkerSize(2.75);

					graphCompareFile4->Draw(drawGraphStyle);
					sprintf(complegendentry,"%s",LegendEntryCompID4);
					plotcompLegend->AddEntry(graphCompareFile4,complegendentry,"lp");
				}


			} //PlotCompare

			char texTex[200];
			if(rapBin==1) sprintf(texTex,"      |#it{y}| < 0.6");
			if(rapBin==2) sprintf(texTex,"0.6 < |#it{y}| < 1.2");
			if(rapBin==3) sprintf(texTex,"1.2 < |#it{y}| < 1.5");
			if(rapBinComb){ 
				sprintf(texTex,"      |#it{y}| < 1.2");
				if(nState>3) sprintf(texTex,"      |#it{y}| < 1.5");
			}
			TLatex *text = new TLatex(PlotpTMax*0.75,yMin+(yMax-yMin)*0.066,texTex);
			text->SetTextSize(0.045); //0.035
			if(!SteerIndividuals&&!PlotMattForICHEP&&!PlotVsComp) text->Draw( "same" );

			if(PlotBG0plots){
				char texTex2[200];
				sprintf(texTex2,"      |c#tau/#sigma_{c#tau}| < 2");
				TLatex *text2 = new TLatex(PlotpTMax*0.7,yMin+(yMax-yMin)*0.122,texTex2);
				text2->SetTextSize(0.05);
				text2->Draw( "same" );
			} 



			if(PlotFinalData&&!PlotMattForICHEP){
				cout<<"DRAW CMS preliminary Latex"<<endl;
				double CentralsFontSize=0.035;
				char text[200];
				sprintf(text,"CMS preliminary");
				cout<<text<<endl;
				TLatex *CentralsText1 = new TLatex(PlotpTMin+(PlotpTMax-PlotpTMin)*0.02-DeltaXminOVERALL,yMin+(yMax-yMin)*0.76,text);
				CentralsText1->SetTextSize(CentralsFontSize*1.25);
				if(DrawPreliminary) CentralsText1->Draw( "same" );
				sprintf(text,"L_{int} = 4.9 fb^{-1}");
				cout<<text<<endl;
				TLatex *CentralsText2 = new TLatex(PlotpTMin+(PlotpTMax-PlotpTMin)*0.02-DeltaXminOVERALL,yMin+(yMax-yMin)*0.935,text);
				CentralsText2->SetTextSize(CentralsFontSize);
				if(DrawLatexStuff) CentralsText2->Draw( "same" );
				sprintf(text,"pp    #sqrt{s} = 7 TeV");
				cout<<text<<endl;
				TLatex *CentralsText3 = new TLatex(PlotpTMin+(PlotpTMax-PlotpTMin)*0.02-DeltaXminOVERALL,yMin+(yMax-yMin)*0.875,text);
				CentralsText3->SetTextSize(CentralsFontSize);
				if(DrawLatexStuff) CentralsText3->Draw( "same" );


				TLine* extreme0 = new TLine( PlotpTMin-DeltaXminOVERALL, 0, PlotpTMax ,0);
				extreme0->SetLineWidth( 1 );
				extreme0->SetLineStyle( 2 );
				extreme0->SetLineColor( kBlack );
				extreme0->Draw( "same" );

				TLine* extreme1 = new TLine( PlotpTMin-DeltaXminOVERALL, 1, PlotpTMax , 1);
				extreme1->SetLineWidth( 1 );
				extreme1->SetLineStyle( 2 );
				extreme1->SetLineColor( kBlack );
				if(iLam==1||iLam==7||iLam==13) extreme1->Draw( "same" );

				TLine* extreme2 = new TLine( PlotpTMin-DeltaXminOVERALL, -1, PlotpTMax ,-1);
				extreme2->SetLineWidth( 1 );
				extreme2->SetLineStyle( 2 );
				extreme2->SetLineColor( kBlack );
				if(iLam==1||iLam==7||iLam==13) extreme2->Draw( "same" );
				if(iLam==6||iLam==12||iLam==18) extreme2->Draw( "same" );

			} //PlotFinalData DrawLatexStuff !PlotMattForICHEP

			if(PlotLegend) plotcompLegend->Draw();

			//plotCanvas->SetLogx(true);
			if(PlotFinalData) plotCanvas->SaveAs(filename);
			plotCanvas->Close();

			delete plotCanvas;

			///////////////////////////////////////////////////////////////
			////////////// End FinalDataResults Plotting //////////////////
			///////////////////////////////////////////////////////////////

			if(MultiPanelPlots&&iLam!=4&&iLam!=5&&iLam!=10&&iLam!=11&&iLam!=16&&iLam!=17){

				bool NEW_MPplots=false;
				bool Psi_MPplots=true;
				bool Psi_MPplots_Old = false;
				bool Psi_MPplots_New = true;

				cout<<"Drawing MultiPanel"<<endl;

				PlotpTMin = 10.;
				PlotpTMax = 75.;//PlotpTMaxInitial;
				if(nState==5){ //May2
					PlotpTMin = 12.;
					PlotpTMax = 52.;
				}

				int MPframe;
				int iPanel;

				// Pad Definitions
				float Top_margin   = 0.;//0
				float Left_margin  = 0.15;//0.025
				float Right_margin = 0.01;//0.015
				const int nPanels=3;
				double lowestBottomMargin=0.2;//0.2
				double PadCoordXMax=0.99;
				double PadCoordYMax=0.985;//0.985
				double deltaCoordY=PadCoordYMax/(double(nPanels-1)+1./(1-lowestBottomMargin));
				double startValCoordY=deltaCoordY/(1-lowestBottomMargin);
				double PadCoordY[nPanels+1]={0.,startValCoordY,startValCoordY+deltaCoordY,PadCoordYMax}; //[0.,0.379,0.682,0.985]
				double PadCoordX[3]={0.1,0.5,0.9};

				// Canvas Definitions
				int MPcanvasXpixelInitial=3000;
				int MPcanvasYpixelInitial=3000;
				int MPcanvasXpixel=MPcanvasXpixelInitial;
				int MPcanvasYpixel=MPcanvasYpixelInitial;
				if(nState==5) MPcanvasXpixel = MPcanvasXpixelInitial * 1.5;

				// Axis Definitions
				double yMinMP=yMin+0.01;
				double yMaxMP=yMax-0.01;
				double LabelSize=0.065;
				double TitleSize=0.085;
				double titleoffset=-0.65;

				double ticksize=0.015;
				int AxisDivisions=510; //510
				double deltaTrickAxisMin=-0.001;
				double deltaTrickAxisMax=-0.001;
				if(rapBin==2) deltaTrickAxisMax=+0.001;

				// Latex definitions
				double whereTexteInPlotX;
				double whereTexteInPlotY;
				double labelOffsetX=0.02;
				double YaxistitleLatexSize=0.12;
				double MPlatexX=11.5;
				if(nState==5) MPlatexX = 12.5;
				double MPlatexYmax=(yMax-yMin)*0.8775+yMin;//0.35;(yMax-yMin)*0.85+yMin
				double MPlatexDeltaYmax=0.09*(yMax-yMin);
				double CentralsFontSizeMP=0.0875;
				// inner legend definitions
				double textSizeRap=0.07825;
				double xRapText;
				double xRapTextTilde;
				double yRapText=0.06;
				if(rapBin==1) xRapText=PlotpTMax * 0.8; //0.825
				if(rapBin==2 || rapBin==3) xRapText=PlotpTMax * 0.69; //0.725
				if(rapBin==1) xRapTextTilde=PlotpTMax * 0.71;
				if(rapBin==2 || rapBin==3) xRapTextTilde=PlotpTMax * 0.60;
				if(nState==5){
					if(rapBin==1) xRapTextTilde=PlotpTMax * 0.66;
					if(rapBin==2 || rapBin==3) xRapTextTilde=PlotpTMax * 0.55;
				}
				double xabcdefText=PlotpTMax * 0.11; //0.225

				double XaxislabelLatexSize=0.0245;
				double YtitleAngle=0.;
				double XtitlePosition=0.; //4.
				double XtitlePositionYshift=0.025;
				if(nState==5) XtitlePosition = 4.5; //May2

				// marker definitions
				double MarkerSizeMP[4]={2.75,2.75,2.75,4.15};// for each frame
				int MarkerColorMP[4] = {1,1,632,600};//{0,600,632,418}; // for each frame
				int MarkerStyleMP[4] = {20,24,25,27}; // for each frame

				// legend
				double errorLegendX1=0.265; //0.165
				double errorLegendX2=0.665; //0.565
				double errorLegendY1=0.7; //0.655
				double errorLegendY2=0.95;
				double errorLegendFontSize=0.07;//0.08


				double xRapText_MPnew=44;
				//======================================================================================

				cout<<"begin Frame dependent plots"<<endl;

				if(Psi_MPplots){
					MPcanvasXpixel=MPcanvasXpixelInitial;
					if(nState==5) MPcanvasXpixel = MPcanvasXpixelInitial * (2.+1./(1-Left_margin)) / (1.+1./(1-Left_margin)) ; // * 1.46
					MPcanvasYpixel = MPcanvasYpixelInitial;
					cout<<"MPcanvasXpixel: "<<MPcanvasXpixel<<" MPcanvasYpixel: "<<MPcanvasYpixel<<endl;

					if(iLam!=6&&iLam!=12&&iLam!=18){

						if(iLam>0&&iLam<7) MPframe=1;
						if(iLam>6&&iLam<13) MPframe=2;
						if(iLam>12&&iLam<19) MPframe=3;

						if(iLam==1||iLam==7||iLam==13) iPanel=1;
						if(iLam==2||iLam==8||iLam==14) iPanel=2;
						if(iLam==3||iLam==9||iLam==15) iPanel=3;

						cout<<"MultiPanel canvas"<<endl;

						if(iLam==1&&rapBin==1){
							MPcanvasCS = new TCanvas("MPcanvasCS", "MPcanvasCS",MPcanvasXpixel,MPcanvasYpixel);
							MPcanvasCS->SetFillColor(kWhite);
							MPcanvasCS->GetFrame()->SetFillColor(kWhite);
							MPcanvasCS->GetFrame()->SetBorderSize(0);
						}
						if(iLam==7&&rapBin==1){
							MPcanvasHX = new TCanvas("MPcanvasHX", "MPcanvasHX",MPcanvasXpixel,MPcanvasYpixel);
							MPcanvasHX->SetFillColor(kWhite);
							MPcanvasHX->GetFrame()->SetFillColor(kWhite);
							MPcanvasHX->GetFrame()->SetBorderSize(0);
						}
						if(iLam==13&&rapBin==1){
							MPcanvasPX = new TCanvas("MPcanvasPX", "MPcanvasPX",MPcanvasXpixel,MPcanvasYpixel);
							MPcanvasPX->SetFillColor(kWhite);
							MPcanvasPX->GetFrame()->SetFillColor(kWhite);
							MPcanvasPX->GetFrame()->SetBorderSize(0);
						}

						if(MPframe==1) MPcanvasCS->cd();
						if(MPframe==2) MPcanvasHX->cd();
						if(MPframe==3) MPcanvasPX->cd();

						cout<<"MultiPanel pad"<<endl;
						TPad *MPpad;
						double X1panel = 0.373;
						double deltaCoordX = 0.;
						if(nState==5){
							deltaCoordX = PadCoordXMax/(2. + 1./(1-Left_margin));
							X1panel = deltaCoordX/(1.-Left_margin);
							if(rapBin==1) MPpad = new TPad("MPpad","MPpad",0.,PadCoordY[nPanels-iPanel],X1panel,PadCoordY[nPanels-iPanel+1]);
							if(rapBin==2) MPpad = new TPad("MPpad","MPpad",X1panel,PadCoordY[nPanels-iPanel],
									(PadCoordXMax+X1panel)/2.,PadCoordY[nPanels-iPanel+1]);
							if(rapBin==3) MPpad = new TPad("MPpad","MPpad",(PadCoordXMax+X1panel)/2.,PadCoordY[nPanels-iPanel],
									PadCoordXMax,PadCoordY[nPanels-iPanel+1]);
						}
						else{
							deltaCoordX = PadCoordXMax/(1. + 1./(1-Left_margin));
							X1panel = deltaCoordX/(1.-Left_margin);
							if(rapBin==1) MPpad = new TPad("MPpad","MPpad",0.,PadCoordY[nPanels-iPanel],X1panel,PadCoordY[nPanels-iPanel+1]);
							if(rapBin==2) MPpad = new TPad("MPpad","MPpad",X1panel,PadCoordY[nPanels-iPanel],
									PadCoordXMax,PadCoordY[nPanels-iPanel+1]);
							//if(rapBin==1) MPpad = new TPad("MPpad","MPpad",0.,PadCoordY[nPanels-iPanel],0.53,PadCoordY[nPanels-iPanel+1]);
							//if(rapBin==2) MPpad = new TPad("MPpad","MPpad",0.53,PadCoordY[nPanels-iPanel],1.,PadCoordY[nPanels-iPanel+1]);
						}

						MPpad->Draw();
						MPpad->cd();
						MPpad->SetFillColor(kWhite);
						MPpad->SetFrameFillColor(kWhite);
						MPpad->SetBorderSize(0);
						MPpad->SetLeftMargin(0.);
						if(rapBin==1) MPpad->SetLeftMargin(Left_margin);
						MPpad->SetRightMargin(0.);
						if(nState<=4 && rapBin==2) MPpad->SetRightMargin(Right_margin);
						if(nState==5 && rapBin==3) MPpad->SetRightMargin(Right_margin);
						MPpad->SetTopMargin(Top_margin+0.0025);
						MPpad->SetBottomMargin(0.0);
						if(iPanel==nPanels) MPpad->SetBottomMargin(lowestBottomMargin);


						cout<<"MultiPanel hist"<<endl;
						TH1F *MPhist = new TH1F;
						if(MPframe==1)MPhist = MPcanvasCS->DrawFrame(PlotpTMin-DeltaXminOVERALL,yMinMP,PlotpTMax,yMaxMP);
						if(MPframe==2)MPhist = MPcanvasHX->DrawFrame(PlotpTMin-DeltaXminOVERALL,yMinMP,PlotpTMax,yMaxMP);
						if(MPframe==3)MPhist = MPcanvasPX->DrawFrame(PlotpTMin-DeltaXminOVERALL,yMinMP,PlotpTMax,yMaxMP);

						MPhist->SetXTitle("#it{p}_{T} [GeV]");
						MPhist->GetXaxis()->SetTitleOffset(-1.35);

						MPhist->SetYTitle(axislabel);
						MPhist->GetYaxis()->SetTitleOffset(titleoffset);
						MPhist->GetYaxis()->SetTitleSize(0.);
						if(iPanel==nPanels) MPhist->GetYaxis()->SetTitleOffset(titleoffset*1.35);
						if(iPanel==nPanels) MPhist->GetYaxis()->SetTitleSize(0.*(1-lowestBottomMargin));

						MPhist->GetYaxis()->SetLabelSize(LabelSize*1.25);
						MPhist->GetXaxis()->SetLabelSize(0.);
						MPhist->GetYaxis()->SetLabelOffset(-0.015);
						MPhist->GetXaxis()->SetLabelOffset(-0.06);

						if(iPanel==nPanels) MPhist->GetYaxis()->SetLabelSize(LabelSize*1.25*(1-lowestBottomMargin));
						MPhist->GetXaxis()->SetTitleSize(TitleSize*0.85);
						MPhist->GetXaxis()->SetAxisColor(kWhite);
						MPhist->GetYaxis()->SetAxisColor(kWhite);
						MPhist->GetXaxis()->SetTicks("-");
						MPhist->GetYaxis()->SetTicks("+");


						TLegend* MPframedepLegend;
						//June3
						MPframedepLegend=new TLegend(errorLegendX1-0.08,errorLegendY1,errorLegendX2-0.33,errorLegendY2);
						if(nState>3&&MPframe==1) 
							MPframedepLegend=new TLegend(errorLegendX1-0.08,errorLegendY1+0.05,errorLegendX2-0.33,errorLegendY2+0.02);

						//MPframedepLegend=new TLegend(errorLegendX1+0.27,errorLegendY1,errorLegendX2+0.1,errorLegendY2);
						//if(nState==5) MPframedepLegend=new TLegend(errorLegendX1+0.25,errorLegendY1-0.05,errorLegendX2,errorLegendY2); //May2

						MPframedepLegend->SetFillColor(0);
						//MPframedepLegend->SetTextFont(72);
						MPframedepLegend->SetTextSize(errorLegendFontSize);

						if(nState>3&&MPframe==1) //June3
							MPframedepLegend->SetTextSize(errorLegendFontSize*0.75);

						if((nState==2||nState==3)&&MPframe==1) MPframedepLegend->SetTextSize(errorLegendFontSize*(1-lowestBottomMargin));
						MPframedepLegend->SetBorderSize(0);

						char MPframedepLegendEntry[200];

						graphSyst->SetMarkerSize(MarkerSizeMP[0]);
						graphSyst->SetMarkerStyle(MarkerStyleMP[0]);
						graphSyst->SetMarkerColor(MarkerColorMP[0]);

						graphDefaultStat->SetMarkerSize(MarkerSizeMP[0]);
						graphDefaultStat->SetMarkerStyle(MarkerStyleMP[0]);
						graphDefaultStat->SetMarkerColor(MarkerColorMP[0]);

						if(!PlotCL1sigma){
							graphDefaultRes3sigma->Draw("2");
							graphDefaultRes2sigma->Draw("2");
						}
						graphDefaultRes->Draw("2");
						if(!PlotAlteredPPDResults) graphSyst->Draw(drawGraphStyle);
						if(PlotAlteredPPDResults) graphDefaultStat->Draw(drawGraphStyle);


						graphDefaultRes->SetLineColor(kGreen);
						graphDefaultRes2sigma->SetLineColor(kYellow);
						graphDefaultRes3sigma->SetLineColor(kCyan-9);

						if(PlotAlteredPPDResults) sprintf(MPframedepLegendEntry,"Stat. uncert., 68.3 %% CL");
						if(!PlotAlteredPPDResults) sprintf(MPframedepLegendEntry,"Tot. sys. uncert.");
						MPframedepLegend->AddEntry(graphSyst,MPframedepLegendEntry,"ple");
						if(PlotAlteredPPDResults) sprintf(MPframedepLegendEntry,"Tot.  uncert., 68.3 %% CL");
						if(!PlotAlteredPPDResults) sprintf(MPframedepLegendEntry,"68.3 %% CL");
						MPframedepLegend->AddEntry(graphDefaultRes,MPframedepLegendEntry,"f");
						if(!PlotCL1sigma){
							if(PlotAlteredPPDResults) sprintf(MPframedepLegendEntry,"Tot.  uncert., 95.5 %% CL");
							if(!PlotAlteredPPDResults) sprintf(MPframedepLegendEntry,"95.5 %% CL");
							MPframedepLegend->AddEntry(graphDefaultRes2sigma,MPframedepLegendEntry,"f");
							if(PlotAlteredPPDResults) sprintf(MPframedepLegendEntry,"Tot.  uncert., 99.7 %% CL");
							if(!PlotAlteredPPDResults) sprintf(MPframedepLegendEntry,"99.7 %% CL");
							MPframedepLegend->AddEntry(graphDefaultRes3sigma,MPframedepLegendEntry,"f");
						}

						if(nState<=4 && rapBin==2) deltaTrickAxisMax=-0.001;
						if(nState==5 && rapBin==3) deltaTrickAxisMax=-0.001;

						TGaxis *axisMPY1 = new TGaxis(PlotpTMin-DeltaXminOVERALL,yMinMP,PlotpTMin-DeltaXminOVERALL,yMaxMP,
								yMinMP,yMaxMP,AxisDivisions,"-US");
						axisMPY1->SetTickSize(ticksize);
						if(iPanel==nPanels) axisMPY1->SetTickSize(ticksize/(1-lowestBottomMargin));
						axisMPY1->Draw("same");

						TGaxis *axisMPY2 = new TGaxis(PlotpTMax,yMinMP,PlotpTMax,yMaxMP,yMinMP,yMaxMP,AxisDivisions,"+US");
						axisMPY2->SetTickSize(ticksize);
						if(iPanel==nPanels) axisMPY2->SetTickSize(ticksize/(1-lowestBottomMargin));
						axisMPY2->Draw("same");


						TGaxis *axisMPX1 = new TGaxis(PlotpTMin-DeltaXminOVERALL,yMinMP,PlotpTMax,yMinMP,
								PlotpTMin-DeltaXminOVERALL+deltaTrickAxisMin,PlotpTMax+deltaTrickAxisMax,AxisDivisions,"+S");
						axisMPX1->SetTickSize(ticksize*2);
						if(iPanel==nPanels) axisMPX1->SetLabelSize(LabelSize);
						if(iPanel<nPanels) axisMPX1->SetLabelSize(0);
						axisMPX1->SetLabelOffset(labelOffsetX);
						if(iPanel==nPanels) axisMPX1->SetTickSize(ticksize/(1-lowestBottomMargin));
						axisMPX1->Draw("same");

						TGaxis *axisMPX2 = new TGaxis(PlotpTMin-DeltaXminOVERALL,yMaxMP,PlotpTMax,yMaxMP,
								PlotpTMin-DeltaXminOVERALL+deltaTrickAxisMin,PlotpTMax+deltaTrickAxisMax,AxisDivisions,"-US");
						axisMPX2->SetTickSize(ticksize*2);
						if(iPanel==nPanels) axisMPX2->SetTickSize(ticksize/(1-lowestBottomMargin));
						axisMPX2->Draw("same");

						XtitlePositionYshift=0.06;
						if(nState==5) XtitlePositionYshift = 0.075; //May2
						whereTexteInPlotX=XtitlePosition;
						whereTexteInPlotY=(yMaxMP+yMinMP)/2.-(yMaxMP-yMinMP)*XtitlePositionYshift;
						if(nState==5) whereTexteInPlotY=(yMaxMP-yMinMP)*XtitlePositionYshift; //May2

						XtitlePositionYshift=0.025;

						char axistitleMPdep[200];
						if(iLam==1||iLam==7||iLam==13)  sprintf(axistitleMPdep,"#lambda_{#vartheta}");
						if(iLam==2||iLam==8||iLam==14)  sprintf(axistitleMPdep,"#lambda_{#varphi}");
						if(iLam==3||iLam==9||iLam==15)  sprintf(axistitleMPdep,"#lambda_{#vartheta#varphi}");

						if(iPanel==nPanels) YaxistitleLatexSize=YaxistitleLatexSize*(1-lowestBottomMargin);
						TLatex *MPYtitletext = new TLatex(whereTexteInPlotX,whereTexteInPlotY ,axistitleMPdep);
						MPYtitletext->SetTextSize(YaxistitleLatexSize);
						if(iPanel==nPanels) MPYtitletext->SetTextSize(YaxistitleLatexSize*(1-lowestBottomMargin));
						MPYtitletext->SetTextColor(kBlack);
						MPYtitletext->PaintLatex(whereTexteInPlotX,whereTexteInPlotY, YtitleAngle, YaxistitleLatexSize, axislabel);
						MPYtitletext->Draw( "same" );

						char frameMPtex[200];
						if(MPframe==1) sprintf(frameMPtex,"CS");
						if(MPframe==2) sprintf(frameMPtex,"HX");
						if(MPframe==3) sprintf(frameMPtex,"PX");
						char texTexMP[200];
						if(rapBin==1) sprintf(texTexMP,"|#it{y}| < 0.6", nState-3, frameMPtex);
						if(rapBin==2) sprintf(texTexMP,"0.6 < |#it{y}| < 1.2", nState-3, frameMPtex);
						if(rapBin==3) sprintf(texTexMP,"1.2 < |#it{y}| < 1.5", nState-3, frameMPtex);
						TLatex *textMP = new TLatex(xRapText,yMin+(yMax-yMin)*yRapText,texTexMP);
						textMP->SetTextSize(textSizeRap);
						if(iPanel==nPanels) textMP->SetTextSize(textSizeRap*(1-lowestBottomMargin));
						textMP->Draw( "same" );

						char abcdef[200];
						if(rapBin==1&&iPanel==1) sprintf(abcdef,"a)");
						if(rapBin==1&&iPanel==2) sprintf(abcdef,"b)");
						if(rapBin==1&&iPanel==3) sprintf(abcdef,"c)");
						if(rapBin==2&&iPanel==1) sprintf(abcdef,"d)");
						if(rapBin==2&&iPanel==2) sprintf(abcdef,"e)");
						if(rapBin==2&&iPanel==3) sprintf(abcdef,"f)");
						if(rapBin==3&&iPanel==1) sprintf(abcdef,"g)");
						if(rapBin==3&&iPanel==2) sprintf(abcdef,"h)");
						if(rapBin==3&&iPanel==3) sprintf(abcdef,"i)");
						cout<<abcdef<<endl;
						TLatex *tex_abcdef = new TLatex(xabcdefText,yMin+(yMax-yMin)*yRapText*0.92,abcdef);
						tex_abcdef->SetTextSize(textSizeRap);
						if(iPanel==nPanels) tex_abcdef->SetTextSize(textSizeRap*(1-lowestBottomMargin));
						//tex_abcdef->Draw( "same" );

						if(PlotFinalData&&DrawLatexStuff){

							TLine* extreme0MP = new TLine( PlotpTMin-DeltaXminOVERALL, 0, PlotpTMax ,0);
							extreme0MP->SetLineWidth( 1 );
							extreme0MP->SetLineStyle( 2 );
							extreme0MP->SetLineColor( kBlack );
							extreme0MP->Draw( "same" );

							TLine* extreme1MP = new TLine( PlotpTMin-DeltaXminOVERALL, 1, PlotpTMax, 1);
							extreme1MP->SetLineWidth( 1 );
							extreme1MP->SetLineStyle( 2 );
							extreme1MP->SetLineColor( kBlack );
							if(iLam==1||iLam==7||iLam==13) extreme1MP->Draw( "same" );

							TLine* extreme2MP = new TLine( PlotpTMin-DeltaXminOVERALL, -1, PlotpTMax ,-1);
							extreme2MP->SetLineWidth( 1 );
							extreme2MP->SetLineStyle( 2 );
							extreme2MP->SetLineColor( kBlack );
							if(iLam==1||iLam==7||iLam==13) extreme2MP->Draw( "same" );
							if(iLam==6||iLam==12||iLam==18) extreme2MP->Draw( "same" );
						}


						if(rapBin==1&&iPanel==1){
							cout<<"DRAW CMS preliminary Latex"<<endl;
							char text[200];
							sprintf(text,"CMS  pp  #sqrt{s} = 7 TeV   L = 4.9 fb^{-1}");
							TLatex *CentralsText1MP = new TLatex(MPlatexX,MPlatexYmax,text);
							CentralsText1MP->SetTextSize(CentralsFontSizeMP);
							CentralsText1MP->Draw( "same" );

							sprintf(text,"Preliminary");
							TLatex *CentralsText2MP = new TLatex(PlotpTMax-22.,MPlatexYmax-1.5*MPlatexDeltaYmax,text);
							if(nState==5) CentralsText2MP = new TLatex(PlotpTMax-13.,MPlatexYmax-1.5*MPlatexDeltaYmax,text); //May2
							CentralsText2MP->SetTextSize(CentralsFontSizeMP);
							CentralsText2MP->Draw( "same" );

							//sprintf(text,"L = 4.9 fb^{-1}");
							//TLatex *CentralsText2MP = new TLatex(MPlatexX,MPlatexYmax-2*MPlatexDeltaYmax,text);
							//CentralsText2MP->SetTextSize(CentralsFontSizeMP);
							//CentralsText2MP->Draw( "same" );
							//sprintf(text,"pp    #sqrt{s} = 7 TeV");
							//TLatex *CentralsText3MP = new TLatex(MPlatexX,MPlatexYmax-MPlatexDeltaYmax,text);
							//CentralsText3MP->SetTextSize(CentralsFontSizeMP);
							//CentralsText3MP->Draw( "same" );
						}

						if(rapBin==1&&iPanel==3&&nState>3&&MPframe==1){ //June3
							MPframedepLegend->Draw("same");
						}

						if(rapBin==1&&iPanel==2&&nState>3&&MPframe>1){ //June3
							MPframedepLegend->Draw("same");
						}

						if(nState==2||nState==3){
							if(rapBin==1&&iPanel==3&&MPframe==1){
								MPframedepLegend->Draw("same");
							}
							if(rapBin==1&&iPanel==2&&MPframe!=1){
								MPframedepLegend->Draw("same");
							}
						}
						if(rapBin==2&&iPanel==1){

							char frameMPtex[200];
							if(MPframe==1) sprintf(frameMPtex,"CS frame");
							if(MPframe==2) sprintf(frameMPtex,"HX frame");
							if(MPframe==3) sprintf(frameMPtex,"PX frame");
							char textStateFrame[200];
							sprintf(textStateFrame,"#psi(%dS), %s", nState-3, frameMPtex);
							if(nState==4) sprintf(textStateFrame,"J/#psi, %s", frameMPtex);
							TLatex *TexStateFrame = new TLatex(MPlatexX,MPlatexYmax,textStateFrame);
							TexStateFrame->SetTextSize(CentralsFontSizeMP);
							TexStateFrame->Draw( "same" );

						}

						if(MPframe==1) MPcanvasCS->cd();
						if(MPframe==2) MPcanvasHX->cd();
						if(MPframe==3) MPcanvasPX->cd();

						whereTexteInPlotX=0.488;
						whereTexteInPlotY=startValCoordY-deltaCoordY-1.375*labelOffsetX;

						TLatex *MPXlabeltext = new TLatex(whereTexteInPlotX,whereTexteInPlotY ,"1");
						MPXlabeltext->SetTextSize(XaxislabelLatexSize);
						MPXlabeltext->SetTextColor(kBlack);
						//if(iPanel==nPanels) MPXlabeltext->Draw( "same" );

						if(iLam==3&&((nState<=4&&rapBin==2)||(nState==5&&rapBin==3))){  
							sprintf(filename,"%s/FinalResultsCS_Psi%dS.pdf",FigDir,nState-3);
							if(PlotFinalData) MPcanvasCS->SaveAs(filename);
							MPcanvasCS->Close();
						}
						if(iLam==9&&((nState<=4&&rapBin==2)||(nState==5&&rapBin==3))){  
							sprintf(filename,"%s/FinalResultsHX_Psi%dS.pdf",FigDir,nState-3);
							if(PlotFinalData) MPcanvasHX->SaveAs(filename);
							MPcanvasHX->Close();
						}
						if(iLam==15&&((nState<=4&&rapBin==2)||(nState==5&&rapBin==3))){  
							sprintf(filename,"%s/FinalResultsPX_Psi%dS.pdf",FigDir,nState-3);
							if(PlotFinalData) MPcanvasPX->SaveAs(filename);
							MPcanvasPX->Close();
						}

					}//end Frame dependent plots
				}


				cout<<"begin Frame independent plots"<<endl;
				cout<<"if(iLam==6||iLam==12||iLam==18)"<<endl;

				if(Psi_MPplots){
					//lowestBottomMargin =  lowestBottomMargin - 0.08;
					lowestBottomMargin =  lowestBottomMargin;

					if(nState>3) {
						if(nState==5) MPcanvasXpixel = MPcanvasXpixelInitial * (2.+1./(1-Left_margin)) / (1.+1./(1-Left_margin) ); //* 1.46
						else MPcanvasXpixel = MPcanvasXpixelInitial;
						MPcanvasYpixel = MPcanvasYpixelInitial  / ( (2.+1./(1-lowestBottomMargin)) * (1.-lowestBottomMargin) ); // *0.384 
					}
					cout<<"MPcanvasXpixel: "<<MPcanvasXpixel<<" MPcanvasYpixel: "<<MPcanvasYpixel<<endl;

					if(iLam==6||iLam==12||iLam==18){

						int mainframe;
						if(iLam==6) mainframe=1;
						if(iLam==12) mainframe=2;
						if(iLam==18) mainframe=3;

						cout<<"iLam = "<<iLam<<endl;

						if((iLam==6||iLam==12||iLam==18)&&rapBin==1){
							MPcanvasTilde = new TCanvas("MPcanvasTilde", "MPcanvasTilde",MPcanvasXpixel,MPcanvasYpixel);
							MPcanvasTilde->SetFillColor(kWhite);
							MPcanvasTilde->GetFrame()->SetFillColor(kWhite);
							MPcanvasTilde->GetFrame()->SetBorderSize(0);
						}

						cout<<"MultiPanel canvas"<<endl;

						MPcanvasTilde->cd();

						cout<<"MultiPanel pad"<<endl;
						TPad *MPpad;
						double X1panel = 0.373;
						double deltaCoordX = 0.;
						if(nState==5){
							deltaCoordX = PadCoordXMax/(2. + 1./(1-Left_margin));
							X1panel = deltaCoordX/(1.-Left_margin);
							if(rapBin==1) MPpad = new TPad("MPpad","MPpad",0.,   0.,X1panel,1.);
							if(rapBin==2) MPpad = new TPad("MPpad","MPpad",X1panel,0.,(PadCoordXMax+X1panel)/2.,1.);
							if(rapBin==3) MPpad = new TPad("MPpad","MPpad",(PadCoordXMax+X1panel)/2.,0.,PadCoordXMax,   1.);
						}
						else{
							//deltaCoordX = PadCoordXMax/(1.+1./(1-Left_margin));
							deltaCoordX = PadCoordXMax/(1./(1-Left_margin) + 1./(1-Right_margin));
							X1panel = deltaCoordX/(1.-Left_margin);
							if(rapBin==1) MPpad = new TPad("MPpad","MPpad",0., 0.,X1panel,1.);
							if(rapBin==2) MPpad = new TPad("MPpad","MPpad",X1panel,0.,PadCoordXMax, 1.);
						}

						MPpad->Draw();
						MPpad->cd();
						MPpad->SetFillColor(kWhite);
						MPpad->SetFrameFillColor(kWhite);
						MPpad->SetBorderSize(0);
						MPpad->SetLeftMargin(0.);
						if(rapBin==1) MPpad->SetLeftMargin(Left_margin);
						MPpad->SetRightMargin(0.);
						if(nState<=4 && rapBin==2) MPpad->SetRightMargin(Right_margin); //+0.05
						if(nState==5 && rapBin==3) MPpad->SetRightMargin(Right_margin); //+0.05
						MPpad->SetTopMargin(Top_margin+0.025);
						MPpad->SetBottomMargin(0.0);
						MPpad->SetBottomMargin(lowestBottomMargin);

						cout<<"MultiPanel hist"<<endl;
						TH1F *MPhist = new TH1F;
						MPhist = MPcanvasTilde->DrawFrame(PlotpTMin-DeltaXminOVERALL,yMinMP,PlotpTMax,yMaxMP);

						MPhist->SetXTitle("#it{p}_{T} [GeV]");
						MPhist->GetXaxis()->SetTitleOffset(-1.35); //-1.2
						//if(rapBin==1) MPhist->GetXaxis()->SetTitleOffset(-1.35);

						MPhist->SetYTitle(axislabel);
						MPhist->GetYaxis()->SetTitleOffset(titleoffset*0.2);
						MPhist->GetYaxis()->SetTitleSize(0.);
						if(rapBin==1) MPhist->GetYaxis()->SetTitleOffset(titleoffset*1.35);
						if(rapBin==1) MPhist->GetYaxis()->SetTitleSize(0.*(1-lowestBottomMargin));

						MPhist->GetYaxis()->SetLabelSize(LabelSize*1.25); //0.5
						MPhist->GetXaxis()->SetLabelSize(0.);
						MPhist->GetYaxis()->SetLabelOffset(-0.015);
						MPhist->GetXaxis()->SetLabelOffset(-0.08);
						//if(rapBin==1) MPhist->GetXaxis()->SetLabelOffset(-0.08*(1-Left_margin));

						if(rapBin==1) MPhist->GetYaxis()->SetLabelSize(LabelSize*1.25*(1-lowestBottomMargin));
						MPhist->GetXaxis()->SetTitleSize(TitleSize*0.85); //0.55
						MPhist->GetXaxis()->SetAxisColor(kWhite);
						MPhist->GetYaxis()->SetAxisColor(kWhite);
						MPhist->GetXaxis()->SetTicks("-");
						MPhist->GetYaxis()->SetTicks("+");

						TLegend* MPframedepLegendError;
						MPframedepLegendError=new TLegend(errorLegendX1-0.25,errorLegendY1+0.1,errorLegendX2,errorLegendY2);
						if(nState==5) MPframedepLegendError=new TLegend(errorLegendX1-0.25,errorLegendY1+0.05,errorLegendX2-0.2,errorLegendY2); //May2
						MPframedepLegendError->SetFillColor(0);
						//MPframedepLegendError->SetTextFont(72);
						MPframedepLegendError->SetTextSize(errorLegendFontSize*0.7);
						MPframedepLegendError->SetBorderSize(0);

						char MPframedepLegendEntry[200];


						TGraphAsymmErrors* graphMP1;
						TGraphAsymmErrors* graphMP2;

						TGraphAsymmErrors* graphMP1_1sig;
						TGraphAsymmErrors* graphMP1_2sig;
						TGraphAsymmErrors* graphMP1_3sig;
						TGraphAsymmErrors* graphMP2_1sig;
						TGraphAsymmErrors* graphMP2_2sig;
						TGraphAsymmErrors* graphMP2_3sig;

						TLegend* MPtildeLegend;
						MPtildeLegend=new TLegend(0.7,0.75,0.9,0.95);
						MPtildeLegend->SetFillColor(0);
						//MPtildeLegend->SetTextFont(72);
						MPtildeLegend->SetTextSize(0.05);
						MPtildeLegend->SetBorderSize(0);
						char MPtildeLegendEntry[200];

						for(int iFrameMP=1;iFrameMP<4;iFrameMP++){

							char GraphNameMP[200];

							if(iFrameMP==1){
								sprintf(GraphNameMP,"ltilde_CS_rap%d",rapBin);
							}
							if(iFrameMP==2){
								sprintf(GraphNameMP,"ltilde_HX_rap%d",rapBin);
							}
							if(iFrameMP==3){
								sprintf(GraphNameMP,"ltilde_PX_rap%d",rapBin);
							}

							if(rapBin<3) graphMP1 = (TGraphAsymmErrors*) infileMP1->Get(GraphNameMP);
							graphMP2 = (TGraphAsymmErrors*) infileMP2->Get(GraphNameMP);

							int MarkerDefinitionForThisBin[4][4]={{0,0,0,0},{0,1,2,3},{0,2,1,3},{0,2,3,1}};

							int ptBinMinMP1=3, ptBinMaxMP1=12;
							int ptBinMinMP2=2, ptBinMaxMP2=5;
							int nBinspTMP1 = ptBinMaxMP1-ptBinMinMP1+1, 
									nBinspTMP2 = ptBinMaxMP2-ptBinMinMP2+1;

							double ptCentreMP1[nBinspTMP1];
							double ptCentreErr_lowMP1[nBinspTMP1];
							double ptCentreErr_highMP1[nBinspTMP1];
							double lmeanMP1[nBinspTMP1];
							double lmean_errlowMP1[nBinspTMP1];
							double lmean_errhighMP1[nBinspTMP1];

							double ptCentreMP2[nBinspTMP2];
							double ptCentreErr_lowMP2[nBinspTMP2];
							double ptCentreErr_highMP2[nBinspTMP2];
							double lmeanMP2[nBinspTMP2];
							double lmean_errlowMP2[nBinspTMP2];
							double lmean_errhighMP2[nBinspTMP2];

							double ShiftTildePlot;
							double ShiftTildePlotZero=ColordBandWidth/2.+0.2; //0.5;
							if(nState==4) ShiftTildePlotZero=ColordBandWidth/2.;
							//if(nState==5) ShiftTildePlotZero = 0.7;
							bool RemoveHorizontalErrorBar=true;

							if(MarkerDefinitionForThisBin[mainframe][iFrameMP]==1) ShiftTildePlot=0.;
							if(MarkerDefinitionForThisBin[mainframe][iFrameMP]==2) ShiftTildePlot=ShiftTildePlotZero;
							if(MarkerDefinitionForThisBin[mainframe][iFrameMP]==3) ShiftTildePlot=-ShiftTildePlotZero;

							int pt=0;
							if(rapBin<3){
								for(int ptBin = ptBinMinMP1; ptBin < ptBinMaxMP1+1; ptBin++) {

									graphMP1->GetPoint(ptBin-1,ptCentreMP1[pt],lmeanMP1[pt]);
									ptCentreErr_highMP1[pt]=graphMP1->GetErrorXhigh(ptBin-1);
									ptCentreErr_lowMP1[pt]=graphMP1->GetErrorXlow(ptBin-1);
									lmean_errhighMP1[pt]=graphMP1->GetErrorYhigh(ptBin-1);
									lmean_errlowMP1[pt]=graphMP1->GetErrorYlow(ptBin-1);

									ptCentreMP1[pt]=ptCentreMP1[pt]+ShiftTildePlot;
									ptCentreErr_highMP1[pt]=ptCentreErr_highMP1[pt]-ShiftTildePlot;
									ptCentreErr_lowMP1[pt]=ptCentreErr_lowMP1[pt]-ShiftTildePlot;
									if(RemoveHorizontalErrorBar){ ptCentreErr_highMP1[pt]=0.; ptCentreErr_lowMP1[pt]=0;}

									pt++;
								}
								graphMP1 = new TGraphAsymmErrors(nBinspTMP1,ptCentreMP1,lmeanMP1,ptCentreErr_lowMP1,ptCentreErr_highMP1,lmean_errlowMP1,lmean_errhighMP1);
							}
							pt=0;
							for(int ptBin = ptBinMinMP2; ptBin < ptBinMaxMP2+1; ptBin++) {

								graphMP2->GetPoint(ptBin-1,ptCentreMP2[pt],lmeanMP2[pt]);
								ptCentreErr_highMP2[pt]=graphMP2->GetErrorXhigh(ptBin-1);
								ptCentreErr_lowMP2[pt]=graphMP2->GetErrorXlow(ptBin-1);
								lmean_errhighMP2[pt]=graphMP2->GetErrorYhigh(ptBin-1);
								lmean_errlowMP2[pt]=graphMP2->GetErrorYlow(ptBin-1);

								ptCentreMP2[pt]=ptCentreMP2[pt]+ShiftTildePlot;
								ptCentreErr_highMP2[pt]=ptCentreErr_highMP2[pt]-ShiftTildePlot;
								ptCentreErr_lowMP2[pt]=ptCentreErr_lowMP2[pt]-ShiftTildePlot;
								if(RemoveHorizontalErrorBar){ ptCentreErr_highMP2[pt]=0.; ptCentreErr_lowMP2[pt]=0;}

								pt++;
							}
							graphMP2 = new TGraphAsymmErrors(nBinspTMP2,ptCentreMP2,lmeanMP2,ptCentreErr_lowMP2,ptCentreErr_highMP2,lmean_errlowMP2,lmean_errhighMP2);

							if(rapBin<3){
								graphMP1->SetMarkerColor(MarkerColorMP[MarkerDefinitionForThisBin[mainframe][iFrameMP]]);
								graphMP1->SetLineColor(MarkerColorMP[MarkerDefinitionForThisBin[mainframe][iFrameMP]]);
								graphMP1->SetMarkerStyle(MarkerStyleMP[MarkerDefinitionForThisBin[mainframe][iFrameMP]]);
								graphMP1->SetMarkerSize(MarkerSizeMP[MarkerDefinitionForThisBin[mainframe][iFrameMP]]);
							}

							graphMP2->SetMarkerColor(MarkerColorMP[MarkerDefinitionForThisBin[mainframe][iFrameMP]]);
							graphMP2->SetLineColor(MarkerColorMP[MarkerDefinitionForThisBin[mainframe][iFrameMP]]);
							graphMP2->SetMarkerStyle(MarkerStyleMP[MarkerDefinitionForThisBin[mainframe][iFrameMP]]);
							graphMP2->SetMarkerSize(MarkerSizeMP[MarkerDefinitionForThisBin[mainframe][iFrameMP]]);


							if(mainframe==1&&iFrameMP==1){ 
								sprintf(MPtildeLegendEntry,"CS"); MPtildeLegend->AddEntry(graphMP1,MPtildeLegendEntry,"p");
								//sprintf(MPtildeLegendEntry,"CS"); MPtildeLegend->AddEntry(graphMP1,MPtildeLegendEntry,"lp");
							}
							if(mainframe!=1&&iFrameMP==1){ 
								sprintf(MPtildeLegendEntry,"CS"); MPtildeLegend->AddEntry(graphMP1,MPtildeLegendEntry,"p"); 
							}

							if(mainframe==2&&iFrameMP==2){ 
								sprintf(MPtildeLegendEntry,"HX"); MPtildeLegend->AddEntry(graphMP1,MPtildeLegendEntry,"p");
								//sprintf(MPtildeLegendEntry,"HX"); MPtildeLegend->AddEntry(graphMP1,MPtildeLegendEntry,"lp"); 
							}
							if(mainframe!=2&&iFrameMP==2){ 
								sprintf(MPtildeLegendEntry,"HX"); MPtildeLegend->AddEntry(graphMP1,MPtildeLegendEntry,"p");
							}

							if(mainframe==3&&iFrameMP==3){ 
								sprintf(MPtildeLegendEntry,"PX"); MPtildeLegend->AddEntry(graphMP1,MPtildeLegendEntry,"p");
								//sprintf(MPtildeLegendEntry,"PX"); MPtildeLegend->AddEntry(graphMP1,MPtildeLegendEntry,"lp"); 
							}

							if(mainframe!=3&&iFrameMP==3){ 
								sprintf(MPtildeLegendEntry,"PX"); MPtildeLegend->AddEntry(graphMP1,MPtildeLegendEntry,"p");
							}



							if(mainframe==1){
								sprintf(GraphNameMP,"ltilde_CS_rap%d",rapBin);
							}
							if(mainframe==2){
								sprintf(GraphNameMP,"ltilde_HX_rap%d",rapBin);
							}
							if(mainframe==3){
								sprintf(GraphNameMP,"ltilde_PX_rap%d",rapBin);
							}

							if(rapBin<3){
								graphMP1_1sig = (TGraphAsymmErrors*) infileMP1_1sig->Get(GraphNameMP);
								graphMP1_2sig = (TGraphAsymmErrors*) infileMP1_2sig->Get(GraphNameMP);
								graphMP1_3sig = (TGraphAsymmErrors*) infileMP1_3sig->Get(GraphNameMP);
							}
							graphMP2_1sig = (TGraphAsymmErrors*) infileMP2_1sig->Get(GraphNameMP);
							graphMP2_2sig = (TGraphAsymmErrors*) infileMP2_2sig->Get(GraphNameMP);
							graphMP2_3sig = (TGraphAsymmErrors*) infileMP2_3sig->Get(GraphNameMP);

							//remove points from the TGraph
							for(int ipoint=1; ipoint<ptBinMinMP1; ipoint++){
								if(rapBin<3){
									graphMP1_1sig->RemovePoint(0);
									graphMP1_2sig->RemovePoint(0);
									graphMP1_3sig->RemovePoint(0);
								}
							}
							for(int ipoint=1; ipoint<ptBinMinMP2; ipoint++){
								graphMP2_1sig->RemovePoint(0);
								graphMP2_2sig->RemovePoint(0);
								graphMP2_3sig->RemovePoint(0);
							}
							//end--remove points from the TGraph

							//remove X bin error to ColordBandWidth
							for(int ipoint=0; ipoint<nBinspTMP1; ipoint++){
								if(rapBin<3){
									graphMP1_1sig->SetPointEXlow(ipoint, ColordBandWidth); graphMP1_1sig->SetPointEXhigh(ipoint, ColordBandWidth);
									graphMP1_2sig->SetPointEXlow(ipoint, ColordBandWidth); graphMP1_2sig->SetPointEXhigh(ipoint, ColordBandWidth);
									graphMP1_3sig->SetPointEXlow(ipoint, ColordBandWidth); graphMP1_3sig->SetPointEXhigh(ipoint, ColordBandWidth);
								}
							}
							for(int ipoint=0; ipoint<nBinspTMP2; ipoint++){
								graphMP2_1sig->SetPointEXlow(ipoint, ColordBandWidth); graphMP2_1sig->SetPointEXhigh(ipoint, ColordBandWidth);
								graphMP2_2sig->SetPointEXlow(ipoint, ColordBandWidth); graphMP2_2sig->SetPointEXhigh(ipoint, ColordBandWidth);
								graphMP2_3sig->SetPointEXlow(ipoint, ColordBandWidth); graphMP2_3sig->SetPointEXhigh(ipoint, ColordBandWidth);
							}
							//end--remove X bin error 

							graphMP1_1sig->SetFillColor(kGreen);
							graphMP1_1sig->SetFillStyle(1001);
							if(PlotCL1sigma) {
								graphMP1_1sig->SetLineColor(CL1sigmaColor);
								graphMP1_1sig->SetFillColor(CL1sigmaColor);
							}
							graphMP1_2sig->SetFillColor(kYellow);
							graphMP1_2sig->SetFillStyle(1001);
							graphMP1_3sig->SetFillColor(kCyan-9);
							graphMP1_3sig->SetFillStyle(1001);

							graphMP2_1sig->SetFillColor(kGreen);
							graphMP2_1sig->SetFillStyle(1001);
							if(PlotCL1sigma) {
								graphMP2_1sig->SetLineColor(CL1sigmaColor);
								graphMP2_1sig->SetFillColor(CL1sigmaColor);
							}
							graphMP2_2sig->SetFillColor(kYellow);
							graphMP2_2sig->SetFillStyle(1001);
							graphMP2_3sig->SetFillColor(kCyan-9);
							graphMP2_3sig->SetFillStyle(1001);

							if(nState==4){

								if(iFrameMP==1){
									if(!PlotCL1sigma){
										graphMP1_3sig->Draw("2");
										graphMP1_2sig->Draw("2");
									}
									graphMP1_1sig->Draw("2");
									if(mainframe==1) graphMP1->Draw(drawGraphStyle);
								}
								if(iFrameMP==2){
									if(mainframe==2) graphMP1->Draw(drawGraphStyle);
								}
								if(iFrameMP==3){
									if(mainframe==3) graphMP1->Draw(drawGraphStyle);
								}

								if(iFrameMP==1){
									if(mainframe!=1) graphMP1->Draw("PX");
								}
								if(iFrameMP==2){
									if(mainframe!=2) graphMP1->Draw("PX");
								}
								if(iFrameMP==3){
									if(mainframe!=3) graphMP1->Draw("PX");
								}

							}
							if(nState==5){

								if(iFrameMP==1){
									if(!PlotCL1sigma){
										graphMP2_3sig->Draw("2");
										graphMP2_2sig->Draw("2");
									}
									graphMP2_1sig->Draw("2");
									if(mainframe==1) graphMP2->Draw(drawGraphStyle);
								}
								if(iFrameMP==2){
									if(mainframe==2) graphMP2->Draw(drawGraphStyle);
								}
								if(iFrameMP==3){
									if(mainframe==3) graphMP2->Draw(drawGraphStyle);
								}

								if(iFrameMP==1){
									if(mainframe!=1) graphMP2->Draw("PX");
								}
								if(iFrameMP==2){
									if(mainframe!=2) graphMP2->Draw("PX");
								}
								if(iFrameMP==3){
									if(mainframe!=3) graphMP2->Draw("PX");
								}

							}
						}

						TGaxis *axisMPY1 = new TGaxis(PlotpTMin-DeltaXminOVERALL,yMinMP,PlotpTMin-DeltaXminOVERALL,yMaxMP,
								yMinMP,yMaxMP,AxisDivisions,"-US");
						axisMPY1->SetTickSize(ticksize);
						axisMPY1->SetTickSize(ticksize/(1-lowestBottomMargin));
						axisMPY1->Draw("same");

						TGaxis *axisMPY2 = new TGaxis(PlotpTMax,yMinMP,PlotpTMax,yMaxMP,yMinMP,yMaxMP,AxisDivisions,"+US");
						axisMPY2->SetTickSize(ticksize);
						axisMPY2->SetTickSize(ticksize/(1-lowestBottomMargin));
						axisMPY2->Draw("same");

						TGaxis *axisM3S1 = new TGaxis(PlotpTMin-DeltaXminOVERALL,yMinMP,PlotpTMax,yMinMP,
								PlotpTMin-DeltaXminOVERALL+deltaTrickAxisMin,PlotpTMax+deltaTrickAxisMax,AxisDivisions,"+S");
						axisM3S1->SetTickSize(ticksize*2);
						axisM3S1->SetLabelSize(LabelSize); //*0.62
						//if(rapBin==1) axisM3S1->SetLabelSize(LabelSize*(1-Left_margin)); //*0.72
						axisM3S1->SetLabelOffset(labelOffsetX);//*0.7
						if(rapBin==1) axisM3S1->SetLabelOffset(labelOffsetX/(1-Left_margin)); //*0.8
						axisM3S1->SetTickSize(ticksize*0.7/(1-lowestBottomMargin));
						axisM3S1->Draw("same");

						TGaxis *axisM3S2 = new TGaxis(PlotpTMin-DeltaXminOVERALL,yMaxMP,PlotpTMax,yMaxMP,
								PlotpTMin-DeltaXminOVERALL+deltaTrickAxisMin,PlotpTMax+deltaTrickAxisMax,AxisDivisions,"-US");
						axisM3S2->SetTickSize(ticksize*2);
						axisM3S2->SetTickSize(ticksize/(1-lowestBottomMargin));
						axisM3S2->Draw("same");

						whereTexteInPlotX=XtitlePosition;
						whereTexteInPlotY=(yMaxMP+yMinMP)/2.-(yMaxMP-yMinMP)*XtitlePositionYshift;

						char axistitleMPtilde[200];
						sprintf(axistitleMPtilde,"#tilde{#lambda}");
						YaxistitleLatexSize=YaxistitleLatexSize*(1-lowestBottomMargin);//*0.7
						TLatex *MPYtitletext = new TLatex(whereTexteInPlotX,whereTexteInPlotY ,axistitleMPtilde);
						MPYtitletext->SetTextSize(YaxistitleLatexSize);
						MPYtitletext->SetTextSize(YaxistitleLatexSize*(1-lowestBottomMargin));
						MPYtitletext->SetTextColor(kBlack);
						MPYtitletext->PaintLatex(whereTexteInPlotX,whereTexteInPlotY, YtitleAngle, YaxistitleLatexSize, axislabel);
						MPYtitletext->Draw( "same" );

						char texTexMP[200];
						if(rapBin==1) sprintf(texTexMP,"#psi(%dS), |#it{y}| < 0.6", nState-3);
						if(rapBin==2) sprintf(texTexMP,"#psi(%dS), 0.6 < |#it{y}| < 1.2", nState-3);
						if(rapBin==3) sprintf(texTexMP,"#psi(%dS), 0.6 < |#it{y}| < 1.5", nState-3);
						if(nState==4){
							if(rapBin==1) sprintf(texTexMP,"J/#psi, |#it{y}| < 0.6");
							if(rapBin==2) sprintf(texTexMP,"J/#psi, 0.6 < |#it{y}| < 1.2");
						}
						TLatex *textMP = new TLatex(xRapTextTilde,yMin+(yMax-yMin)*yRapText,texTexMP);//*0.92
						textMP->SetTextSize(textSizeRap*0.8);//*0.7
						//if(rapBin==1) textMP->SetTextSize(textSizeRap*(1-Left_margin));//*0.75
						textMP->Draw( "same" );

						char abcdef[200];
						if(rapBin==1) sprintf(abcdef,"a)");
						if(rapBin==2) sprintf(abcdef,"b)");
						if(rapBin==3) sprintf(abcdef,"c)");
						cout<<abcdef<<endl;
						TLatex *tex_abcdef = new TLatex(xabcdefText,yMin+(yMax-yMin)*yRapText*0.92,abcdef);
						tex_abcdef->SetTextSize(textSizeRap);
						tex_abcdef->SetTextSize(textSizeRap*(1-lowestBottomMargin));
						//tex_abcdef->Draw( "same" );

						if(PlotFinalData&&DrawLatexStuff){

							TLine* extreme0MP = new TLine( PlotpTMin-DeltaXminOVERALL, 0, PlotpTMax ,0);
							extreme0MP->SetLineWidth( 1 );
							extreme0MP->SetLineStyle( 2 );
							extreme0MP->SetLineColor( kBlack );
							extreme0MP->Draw( "same" );

							TLine* extreme1MP = new TLine( PlotpTMin-DeltaXminOVERALL, 1, PlotpTMax , 1);
							extreme1MP->SetLineWidth( 1 );
							extreme1MP->SetLineStyle( 2 );
							extreme1MP->SetLineColor( kBlack );
							if(iLam==1||iLam==7||iLam==13) extreme1MP->Draw( "same" );

							TLine* extreme2MP = new TLine( PlotpTMin-DeltaXminOVERALL, -1, PlotpTMax ,-1);
							extreme2MP->SetLineWidth( 1 );
							extreme2MP->SetLineStyle( 2 );
							extreme2MP->SetLineColor( kBlack );
							if(iLam==1||iLam==7||iLam==13) extreme2MP->Draw( "same" );
							if(iLam==6||iLam==12||iLam==18) extreme2MP->Draw( "same" );

						}


						if(rapBin==1){
							cout<<"DRAW CMS preliminary Latex"<<endl;
							char text[200];
							sprintf(text,"CMS  pp  #sqrt{s} = 7 TeV   L = 4.9 fb^{-1}");
							TLatex *CentralsText1MP = new TLatex(MPlatexX,MPlatexYmax,text);
							CentralsText1MP->SetTextSize(CentralsFontSizeMP*0.75);
							CentralsText1MP->Draw( "same" );

							sprintf(text,"Preliminary");
							TLatex *CentralsText2MP = new TLatex(PlotpTMax-22,MPlatexYmax-1.5*MPlatexDeltaYmax,text);
							if(nState==5) CentralsText2MP = new TLatex(PlotpTMax-13.,MPlatexYmax-1.5*MPlatexDeltaYmax,text); //May2
							CentralsText2MP->SetTextSize(CentralsFontSizeMP*0.75);
							CentralsText2MP->Draw( "same" );

							//sprintf(text,"L = 4.9 fb^{-1}");
							//TLatex *CentralsText2MP = new TLatex(MPlatexX,MPlatexYmax-2*MPlatexDeltaYmax,text);
							//CentralsText2MP->SetTextSize(CentralsFontSizeMP);
							//CentralsText2MP->Draw( "same" );
							//sprintf(text,"pp    #sqrt{s} = 7 TeV");
							//TLatex *CentralsText3MP = new TLatex(MPlatexX,MPlatexYmax-MPlatexDeltaYmax,text);
							//CentralsText3MP->SetTextSize(CentralsFontSizeMP);
							//CentralsText3MP->Draw( "same" );

						}

						if(rapBin==2){
							graphMP1_1sig->SetLineColor(kGreen);
							graphMP1_2sig->SetLineColor(kYellow);
							graphMP1_3sig->SetLineColor(kCyan-9);

							TGraphAsymmErrors *legendPhantom = (TGraphAsymmErrors*) infileMP2_3sig->Get("ltilde_CS_rap1");

							legendPhantom->SetMarkerColor(MarkerColorMP[0]);
							legendPhantom->SetLineColor(MarkerColorMP[0]);
							legendPhantom->SetMarkerStyle(MarkerStyleMP[1]);
							legendPhantom->SetMarkerSize(MarkerSizeMP[0]);


							if(PlotAlteredPPDResults) sprintf(MPframedepLegendEntry,"Stat. uncert., 68.3 %% CL");
							if(!PlotAlteredPPDResults) sprintf(MPframedepLegendEntry,"Tot. sys. uncert.");
							MPframedepLegendError->AddEntry(legendPhantom,MPframedepLegendEntry,"ple");
							if(PlotAlteredPPDResults) sprintf(MPframedepLegendEntry,"Tot.  uncert., 68.3 %% CL");
							if(!PlotAlteredPPDResults) sprintf(MPframedepLegendEntry,"68.3 %% CL");
							MPframedepLegendError->AddEntry(graphMP1_1sig,MPframedepLegendEntry,"f");
							if(!PlotCL1sigma){
								if(PlotAlteredPPDResults) sprintf(MPframedepLegendEntry,"Tot.  uncert., 95.5 %% CL");
								if(!PlotAlteredPPDResults) sprintf(MPframedepLegendEntry,"95.5 %% CL");
								MPframedepLegendError->AddEntry(graphMP1_2sig,MPframedepLegendEntry,"f");
								if(PlotAlteredPPDResults) sprintf(MPframedepLegendEntry,"Tot.  uncert., 99.7 %% CL");
								if(!PlotAlteredPPDResults) sprintf(MPframedepLegendEntry,"99.7 %% CL");
								MPframedepLegendError->AddEntry(graphMP1_3sig,MPframedepLegendEntry,"f");
							}

							MPtildeLegend->Draw("same");
							MPframedepLegendError->Draw("same");

						}

						//if(rapBin==2){
						//	char frameMPtex[200];
						//	if(MPframe==1) sprintf(frameMPtex,"CS frame");
						//	if(MPframe==2) sprintf(frameMPtex,"HX frame");
						//	if(MPframe==3) sprintf(frameMPtex,"PX frame");
						//	char textStateFrame[200];
						//	sprintf(textStateFrame,"%s", frameMPtex);
						//	TLatex *TexStateFrame = new TLatex(MPlatexX,MPlatexYmax,textStateFrame);
						//	TexStateFrame->SetTextSize(CentralsFontSizeMP);
						//	TexStateFrame->Draw( "same" );
						//}

						MPcanvasTilde->cd();

						whereTexteInPlotX=0.488;
						whereTexteInPlotY=startValCoordY-deltaCoordY-1.425*labelOffsetX;

						TLatex *M3Slabeltext = new TLatex(whereTexteInPlotX,whereTexteInPlotY ,"1");
						M3Slabeltext->SetTextSize(XaxislabelLatexSize);
						M3Slabeltext->SetTextColor(kBlack);
						//M3Slabeltext->Draw( "same" );

						if((iLam==6||iLam==12||iLam==18)&&((nState<=4&&rapBin==2)||(nState==5&&rapBin==3))){
							if(mainframe==1) sprintf(filename,"%s/FinalResultsTildeCS_Psi%dS.pdf",FigDir,nState-3);
							if(mainframe==2) sprintf(filename,"%s/FinalResultsTildeHX_Psi%dS.pdf",FigDir,nState-3);
							if(mainframe==3) sprintf(filename,"%s/FinalResultsTildePX_Psi%dS.pdf",FigDir,nState-3);
							if(PlotFinalData) MPcanvasTilde->SaveAs(filename);
							MPcanvasTilde->Close();
						}
					}//end Frame independent plots
				}

				//======================================================================================


				cout<<"begin Frame dependent plots -- NEW"<<endl;

				if(Psi_MPplots_New){
					MPcanvasXpixel=MPcanvasXpixelInitial;
					MPcanvasYpixel = MPcanvasYpixelInitial;
					cout<<"MPcanvasXpixel: "<<MPcanvasXpixel<<" MPcanvasYpixel: "<<MPcanvasYpixel<<endl;

					if(iLam!=6&&iLam!=12&&iLam!=18&&rapBin==1){

						//yMinMP = -0.99; yMaxMP = 0.99; 	
						//if(iLam==2||iLam==3||iLam==8||iLam==9||iLam==14||iLam==15){
						//	yMinMP = -0.33; yMaxMP = 0.33;
						//}

						if(iLam>0&&iLam<7) MPframe=1;
						if(iLam>6&&iLam<13) MPframe=2;
						if(iLam>12&&iLam<19) MPframe=3;

						if(iLam==1||iLam==7||iLam==13) iPanel=1;
						if(iLam==2||iLam==8||iLam==14) iPanel=2;
						if(iLam==3||iLam==9||iLam==15) iPanel=3;

						cout<<"MultiPanel canvas"<<endl;

						if(iLam==1&&rapBin==1){
							MPcanvasCS_New = new TCanvas("MPcanvasCS_New", "MPcanvasCS_New",MPcanvasXpixel,MPcanvasYpixel);
							MPcanvasCS_New->SetFillColor(kWhite);
							MPcanvasCS_New->GetFrame()->SetFillColor(kWhite);
							MPcanvasCS_New->GetFrame()->SetBorderSize(0);
						}
						if(iLam==7&&rapBin==1){
							MPcanvasHX_New = new TCanvas("MPcanvasHX_New", "MPcanvasHX_New",MPcanvasXpixel,MPcanvasYpixel);
							MPcanvasHX_New->SetFillColor(kWhite);
							MPcanvasHX_New->GetFrame()->SetFillColor(kWhite);
							MPcanvasHX_New->GetFrame()->SetBorderSize(0);
						}
						if(iLam==13&&rapBin==1){
							MPcanvasPX_New = new TCanvas("MPcanvasPX_New", "MPcanvasPX_New",MPcanvasXpixel,MPcanvasYpixel);
							MPcanvasPX_New->SetFillColor(kWhite);
							MPcanvasPX_New->GetFrame()->SetFillColor(kWhite);
							MPcanvasPX_New->GetFrame()->SetBorderSize(0);
						}

						int totalState = 2;
						for(int iStateMP=1;iStateMP<totalState+1;iStateMP++){ //

							PlotpTMin = 12.; PlotpTMax = 72.;
							//if(iStateMP==1){ PlotpTMin = 10.; PlotpTMax = 75.; }
							//if(iStateMP==2){ PlotpTMin = 12.; PlotpTMax = 52.; }

							if(MPframe==1) MPcanvasCS_New->cd();
							if(MPframe==2) MPcanvasHX_New->cd();
							if(MPframe==3) MPcanvasPX_New->cd();

							cout<<"MultiPanel pad"<<endl;
							cout<<"MultiPanel pad"<<endl;
							cout<<"iLam "<<iLam<<endl;
							cout<<"MPframe "<<MPframe<<endl;
							cout<<"iStateMP "<<iStateMP<<endl;
							cout<<"iPanel "<<iPanel<<endl;

							TPad *MPpad;
							double X1panel = 0.373;
							double deltaCoordX = 0.;
							deltaCoordX = PadCoordXMax/(1. + 1./(1-Left_margin));
							X1panel = deltaCoordX/(1.-Left_margin);
							if(iStateMP==1) MPpad = new TPad("MPpad","MPpad",0.,PadCoordY[nPanels-iPanel],X1panel,PadCoordY[nPanels-iPanel+1]);
							if(iStateMP==2) MPpad = new TPad("MPpad","MPpad",X1panel,PadCoordY[nPanels-iPanel], PadCoordXMax,PadCoordY[nPanels-iPanel+1]);

							MPpad->Draw();
							MPpad->cd();
							MPpad->SetFillColor(kWhite);
							MPpad->SetFrameFillColor(kWhite);
							MPpad->SetBorderSize(0);
							MPpad->SetLeftMargin(0.);
							if(iStateMP==1) MPpad->SetLeftMargin(Left_margin);
							MPpad->SetRightMargin(0.);
							if(iStateMP==2) MPpad->SetRightMargin(Right_margin);
							MPpad->SetTopMargin(Top_margin+0.0025);
							MPpad->SetBottomMargin(0.0);
							if(iPanel==nPanels) MPpad->SetBottomMargin(lowestBottomMargin);


							cout<<"MultiPanel hist"<<endl;
							TH1F *MPhist = new TH1F;
							if(MPframe==1)MPhist = MPcanvasCS_New->DrawFrame(PlotpTMin-DeltaXminOVERALL,yMinMP,PlotpTMax,yMaxMP);
							if(MPframe==2)MPhist = MPcanvasHX_New->DrawFrame(PlotpTMin-DeltaXminOVERALL,yMinMP,PlotpTMax,yMaxMP);
							if(MPframe==3)MPhist = MPcanvasPX_New->DrawFrame(PlotpTMin-DeltaXminOVERALL,yMinMP,PlotpTMax,yMaxMP);

							MPhist->SetXTitle("#it{p}_{T} [GeV]");
							MPhist->GetXaxis()->SetTitleOffset(-1.35);

							MPhist->SetYTitle(axislabel);
							MPhist->GetYaxis()->SetTitleOffset(titleoffset);
							MPhist->GetYaxis()->SetTitleSize(0.);
							if(iPanel==nPanels) MPhist->GetYaxis()->SetTitleOffset(titleoffset*1.35);
							if(iPanel==nPanels) MPhist->GetYaxis()->SetTitleSize(0.*(1-lowestBottomMargin));

							MPhist->GetYaxis()->SetLabelSize(LabelSize*1.25);
							MPhist->GetXaxis()->SetLabelSize(0.);
							MPhist->GetYaxis()->SetLabelOffset(-0.015);
							MPhist->GetXaxis()->SetLabelOffset(-0.06);

							if(iPanel==nPanels) MPhist->GetYaxis()->SetLabelSize(LabelSize*1.25*(1-lowestBottomMargin));
							MPhist->GetXaxis()->SetTitleSize(TitleSize*0.85);
							MPhist->GetXaxis()->SetAxisColor(kWhite);
							MPhist->GetYaxis()->SetAxisColor(kWhite);
							MPhist->GetXaxis()->SetTicks("-");
							MPhist->GetYaxis()->SetTicks("+");


							TLegend* MPframedepLegend;
							MPframedepLegend=new TLegend(errorLegendX1-0.18,errorLegendY1+0.05,errorLegendX2-0.33,errorLegendY2+0.02);

							MPframedepLegend->SetFillColor(0);
							//MPframedepLegend->SetTextFont(72);
							MPframedepLegend->SetTextSize(errorLegendFontSize);

							if(MPframe==1) //June3
								MPframedepLegend->SetTextSize(errorLegendFontSize*0.75);

							MPframedepLegend->SetBorderSize(0);

							char MPframedepLegendEntry[200];

							TGraphAsymmErrors* graphMP_1sig_rap1;
							TGraphAsymmErrors* graphMP_1sig_rap2;
							TGraphAsymmErrors* graphMP_1sig_rap3;

							char GraphNameRap1[100], GraphNameRap2[100], GraphNameRap3[100];

							//rap1
							if(iLam==1)  sprintf(GraphNameRap1,"lth_CS_rap1");
							if(iLam==2)  sprintf(GraphNameRap1,"lph_CS_rap1");
							if(iLam==3)  sprintf(GraphNameRap1,"ltp_CS_rap1");
							if(iLam==4)  sprintf(GraphNameRap1,"lthstar_CS_rap1");
							if(iLam==5)  sprintf(GraphNameRap1,"lphstar_CS_rap1");
							if(iLam==6)  sprintf(GraphNameRap1,"ltilde_CS_rap1");

							if(iLam==7)  sprintf(GraphNameRap1,"lth_HX_rap1");
							if(iLam==8)  sprintf(GraphNameRap1,"lph_HX_rap1");
							if(iLam==9)  sprintf(GraphNameRap1,"ltp_HX_rap1");
							if(iLam==10) sprintf(GraphNameRap1,"lthstar_HX_rap1");
							if(iLam==11) sprintf(GraphNameRap1,"lphstar_HX_rap1");
							if(iLam==12) sprintf(GraphNameRap1,"ltilde_HX_rap1");

							if(iLam==13) sprintf(GraphNameRap1,"lth_PX_rap1");
							if(iLam==14) sprintf(GraphNameRap1,"lph_PX_rap1");
							if(iLam==15) sprintf(GraphNameRap1,"ltp_PX_rap1");
							if(iLam==16) sprintf(GraphNameRap1,"lthstar_PX_rap1");
							if(iLam==17) sprintf(GraphNameRap1,"lphstar_PX_rap1");
							if(iLam==18) sprintf(GraphNameRap1,"ltilde_PX_rap1");

							//rap2
							if(iLam==1)  sprintf(GraphNameRap2,"lth_CS_rap2");
							if(iLam==2)  sprintf(GraphNameRap2,"lph_CS_rap2");
							if(iLam==3)  sprintf(GraphNameRap2,"ltp_CS_rap2");
							if(iLam==4)  sprintf(GraphNameRap2,"lthstar_CS_rap2");
							if(iLam==5)  sprintf(GraphNameRap2,"lphstar_CS_rap2");
							if(iLam==6)  sprintf(GraphNameRap2,"ltilde_CS_rap2");

							if(iLam==7)  sprintf(GraphNameRap2,"lth_HX_rap2");
							if(iLam==8)  sprintf(GraphNameRap2,"lph_HX_rap2");
							if(iLam==9)  sprintf(GraphNameRap2,"ltp_HX_rap2");
							if(iLam==10) sprintf(GraphNameRap2,"lthstar_HX_rap2");
							if(iLam==11) sprintf(GraphNameRap2,"lphstar_HX_rap2");
							if(iLam==12) sprintf(GraphNameRap2,"ltilde_HX_rap2");

							if(iLam==13) sprintf(GraphNameRap2,"lth_PX_rap2");
							if(iLam==14) sprintf(GraphNameRap2,"lph_PX_rap2");
							if(iLam==15) sprintf(GraphNameRap2,"ltp_PX_rap2");
							if(iLam==16) sprintf(GraphNameRap2,"lthstar_PX_rap2");
							if(iLam==17) sprintf(GraphNameRap2,"lphstar_PX_rap2");
							if(iLam==18) sprintf(GraphNameRap2,"ltilde_PX_rap2");

							//rap3
							if(iLam==1)  sprintf(GraphNameRap3,"lth_CS_rap3");
							if(iLam==2)  sprintf(GraphNameRap3,"lph_CS_rap3");
							if(iLam==3)  sprintf(GraphNameRap3,"ltp_CS_rap3");
							if(iLam==4)  sprintf(GraphNameRap3,"lthstar_CS_rap3");
							if(iLam==5)  sprintf(GraphNameRap3,"lphstar_CS_rap3");
							if(iLam==6)  sprintf(GraphNameRap3,"ltilde_CS_rap3");

							if(iLam==7)  sprintf(GraphNameRap3,"lth_HX_rap3");
							if(iLam==8)  sprintf(GraphNameRap3,"lph_HX_rap3");
							if(iLam==9)  sprintf(GraphNameRap3,"ltp_HX_rap3");
							if(iLam==10) sprintf(GraphNameRap3,"lthstar_HX_rap3");
							if(iLam==11) sprintf(GraphNameRap3,"lphstar_HX_rap3");
							if(iLam==12) sprintf(GraphNameRap3,"ltilde_HX_rap3");

							if(iLam==13) sprintf(GraphNameRap3,"lth_PX_rap3");
							if(iLam==14) sprintf(GraphNameRap3,"lph_PX_rap3");
							if(iLam==15) sprintf(GraphNameRap3,"ltp_PX_rap3");
							if(iLam==16) sprintf(GraphNameRap3,"lthstar_PX_rap3");
							if(iLam==17) sprintf(GraphNameRap3,"lphstar_PX_rap3");
							if(iLam==18) sprintf(GraphNameRap3,"ltilde_PX_rap3");

							if(iStateMP==1) {
								graphMP_1sig_rap1 = (TGraphAsymmErrors*) infileMP1_1sig -> Get(GraphNameRap1);
								graphMP_1sig_rap2 = (TGraphAsymmErrors*) infileMP1_1sig -> Get(GraphNameRap2);
							}
							if(iStateMP==2) {
								graphMP_1sig_rap1 = (TGraphAsymmErrors*) infileMP2_1sig -> Get(GraphNameRap1);
								graphMP_1sig_rap2 = (TGraphAsymmErrors*) infileMP2_1sig -> Get(GraphNameRap2);
								graphMP_1sig_rap3 = (TGraphAsymmErrors*) infileMP2_1sig -> Get(GraphNameRap3);
							}

							int ptBinMinMP=3, ptBinMaxMP=12;
							if(iStateMP==2) { ptBinMinMP=2; ptBinMaxMP=5; }
							int nBinspTMP = ptBinMaxMP-ptBinMinMP+1;

							double ptCentreMP[nBinspTMP];
							double ptCentreErr_lowMP[nBinspTMP];
							double ptCentreErr_highMP[nBinspTMP];
							double lmeanMP[nBinspTMP];
							double lmean_errlowMP[nBinspTMP];
							double lmean_errhighMP[nBinspTMP];

							double ShiftPlot=0.5;
							bool RemoveHorizontalErrorBar=true;

							int pt=0;

							for(int ptBin = ptBinMinMP; ptBin < ptBinMaxMP+1; ptBin++) {

								graphMP_1sig_rap1->GetPoint(ptBin-1,ptCentreMP[pt],lmeanMP[pt]);
								ptCentreErr_highMP[pt]=graphMP_1sig_rap1->GetErrorXhigh(ptBin-1);
								ptCentreErr_lowMP[pt]=graphMP_1sig_rap1->GetErrorXlow(ptBin-1);
								lmean_errhighMP[pt]=graphMP_1sig_rap1->GetErrorYhigh(ptBin-1);
								lmean_errlowMP[pt]=graphMP_1sig_rap1->GetErrorYlow(ptBin-1);

								ptCentreMP[pt]=ptCentreMP[pt] -  ShiftPlot;
								ptCentreErr_highMP[pt]=ptCentreErr_highMP[pt] - ShiftPlot;
								ptCentreErr_lowMP[pt]=ptCentreErr_lowMP[pt] - ShiftPlot;
								if(RemoveHorizontalErrorBar){ ptCentreErr_highMP[pt]=0.; ptCentreErr_lowMP[pt]=0;}

								pt++;
							}
							graphMP_1sig_rap1 = new TGraphAsymmErrors(nBinspTMP,ptCentreMP,lmeanMP, ptCentreErr_lowMP,ptCentreErr_highMP,lmean_errlowMP,lmean_errhighMP);
							pt=0;
							for(int ptBin = ptBinMinMP; ptBin < ptBinMaxMP+1; ptBin++) {

								graphMP_1sig_rap2->GetPoint(ptBin-1,ptCentreMP[pt],lmeanMP[pt]);
								ptCentreErr_highMP[pt]=graphMP_1sig_rap2->GetErrorXhigh(ptBin-1);
								ptCentreErr_lowMP[pt]=graphMP_1sig_rap2->GetErrorXlow(ptBin-1);
								lmean_errhighMP[pt]=graphMP_1sig_rap2->GetErrorYhigh(ptBin-1);
								lmean_errlowMP[pt]=graphMP_1sig_rap2->GetErrorYlow(ptBin-1);

								ptCentreMP[pt]=ptCentreMP[pt];
								ptCentreErr_highMP[pt]=ptCentreErr_highMP[pt];
								ptCentreErr_lowMP[pt]=ptCentreErr_lowMP[pt];
								if(RemoveHorizontalErrorBar){ ptCentreErr_highMP[pt]=0.; ptCentreErr_lowMP[pt]=0;}

								pt++;
							}
							graphMP_1sig_rap2 = new TGraphAsymmErrors(nBinspTMP,ptCentreMP,lmeanMP, ptCentreErr_lowMP,ptCentreErr_highMP,lmean_errlowMP,lmean_errhighMP);

							if(iStateMP==2){
								pt=0;
								for(int ptBin = ptBinMinMP; ptBin < ptBinMaxMP+1; ptBin++) {

									graphMP_1sig_rap3->GetPoint(ptBin-1,ptCentreMP[pt],lmeanMP[pt]);
									ptCentreErr_highMP[pt]=graphMP_1sig_rap3->GetErrorXhigh(ptBin-1);
									ptCentreErr_lowMP[pt]=graphMP_1sig_rap3->GetErrorXlow(ptBin-1);
									lmean_errhighMP[pt]=graphMP_1sig_rap3->GetErrorYhigh(ptBin-1);
									lmean_errlowMP[pt]=graphMP_1sig_rap3->GetErrorYlow(ptBin-1);

									ptCentreMP[pt]=ptCentreMP[pt] + ShiftPlot;
									ptCentreErr_highMP[pt]=ptCentreErr_highMP[pt] + ShiftPlot;
									ptCentreErr_lowMP[pt]=ptCentreErr_lowMP[pt] + ShiftPlot;
									if(RemoveHorizontalErrorBar){ ptCentreErr_highMP[pt]=0.; ptCentreErr_lowMP[pt]=0;}

									pt++;
								}
								graphMP_1sig_rap3 = new TGraphAsymmErrors(nBinspTMP,ptCentreMP,lmeanMP, ptCentreErr_lowMP,ptCentreErr_highMP,lmean_errlowMP,lmean_errhighMP);
							}

							graphMP_1sig_rap1 -> SetLineColor(kBlack);
							graphMP_1sig_rap1 -> SetMarkerColor(kBlack);
							graphMP_1sig_rap1 -> SetMarkerStyle(20);
							graphMP_1sig_rap1 -> SetMarkerSize(2.75);

							graphMP_1sig_rap2 -> SetLineColor(kBlue); 
							graphMP_1sig_rap2 -> SetMarkerColor(kBlue);
							graphMP_1sig_rap2 -> SetMarkerStyle(25);
							graphMP_1sig_rap2 -> SetMarkerSize(2.75);

							if(iStateMP==2){
								graphMP_1sig_rap3 -> SetLineColor(kRed);
								graphMP_1sig_rap3 -> SetMarkerColor(kRed);
								graphMP_1sig_rap3 -> SetMarkerStyle(30);
								graphMP_1sig_rap3 -> SetMarkerSize(2.75);
							}

							graphMP_1sig_rap1 -> Draw(drawGraphStyle);
							graphMP_1sig_rap2 -> Draw("psame");
							if(iStateMP==2)
								graphMP_1sig_rap3 -> Draw("psame");

							sprintf(MPframedepLegendEntry,"|#it{y}| < 0.6");
							MPframedepLegend->AddEntry(graphMP_1sig_rap1,MPframedepLegendEntry,"pl");
							sprintf(MPframedepLegendEntry,"0.6 < |#it{y}| < 1.2");
							MPframedepLegend->AddEntry(graphMP_1sig_rap2,MPframedepLegendEntry,"pl");
							if(iStateMP==2) {
								sprintf(MPframedepLegendEntry,"1.2 < |#it{y}| < 1.5");
								MPframedepLegend->AddEntry(graphMP_1sig_rap3,MPframedepLegendEntry,"pl");
							}

							if(iStateMP==2) deltaTrickAxisMax=-0.001;

							TGaxis *axisMPY1 = new TGaxis(PlotpTMin-DeltaXminOVERALL,yMinMP,PlotpTMin-DeltaXminOVERALL,yMaxMP,
									yMinMP,yMaxMP,AxisDivisions,"-US");
							//yMinMP,yMaxMP,308,"-US");
							axisMPY1->SetTickSize(ticksize);
							//axisMPY1->SetNdivisions(406);
							if(iPanel==nPanels) axisMPY1->SetTickSize(ticksize/(1-lowestBottomMargin));
							axisMPY1->Draw("same");

							TGaxis *axisMPY2 = new TGaxis(PlotpTMax,yMinMP,PlotpTMax,yMaxMP,yMinMP,yMaxMP,AxisDivisions,"+US");
							axisMPY2->SetTickSize(ticksize);
							if(iPanel==nPanels) axisMPY2->SetTickSize(ticksize/(1-lowestBottomMargin));
							axisMPY2->Draw("same");


							TGaxis *axisMPX1 = new TGaxis(PlotpTMin-DeltaXminOVERALL,yMinMP,PlotpTMax,yMinMP,
									PlotpTMin-DeltaXminOVERALL+deltaTrickAxisMin,PlotpTMax+deltaTrickAxisMax,AxisDivisions,"+S");
							axisMPX1->SetTickSize(ticksize*2);
							if(iPanel==nPanels) axisMPX1->SetLabelSize(LabelSize);
							if(iPanel<nPanels) axisMPX1->SetLabelSize(0);
							axisMPX1->SetLabelOffset(labelOffsetX);
							if(iPanel==nPanels) axisMPX1->SetTickSize(ticksize/(1-lowestBottomMargin));
							axisMPX1->Draw("same");

							TGaxis *axisMPX2 = new TGaxis(PlotpTMin-DeltaXminOVERALL,yMaxMP,PlotpTMax,yMaxMP,
									PlotpTMin-DeltaXminOVERALL+deltaTrickAxisMin,PlotpTMax+deltaTrickAxisMax,AxisDivisions,"-US");
							axisMPX2->SetTickSize(ticksize*2);
							if(iPanel==nPanels) axisMPX2->SetTickSize(ticksize/(1-lowestBottomMargin));
							axisMPX2->Draw("same");

							XtitlePositionYshift=0.06;
							whereTexteInPlotX=2; 
							whereTexteInPlotY=(yMaxMP+yMinMP)/2.-(yMaxMP-yMinMP)*XtitlePositionYshift;

							XtitlePositionYshift=0.025;

							char axistitleMPdep[200];
							if(iLam==1||iLam==7||iLam==13)  sprintf(axistitleMPdep,"#lambda_{#vartheta}");
							if(iLam==2||iLam==8||iLam==14)  sprintf(axistitleMPdep,"#lambda_{#varphi}");
							if(iLam==3||iLam==9||iLam==15)  sprintf(axistitleMPdep,"#lambda_{#vartheta#varphi}");

							TLatex *MPYtitletext = new TLatex(whereTexteInPlotX,whereTexteInPlotY ,axistitleMPdep);
							MPYtitletext->SetTextSize(YaxistitleLatexSize);
							MPYtitletext->SetTextColor(kBlack);
							MPYtitletext->PaintLatex(whereTexteInPlotX,whereTexteInPlotY, YtitleAngle, YaxistitleLatexSize, axislabel);
							MPYtitletext->Draw( "same" );

							double increaseSize=1.45;
							double SpecialShiftUpsilonLabel=0.9;
							double SpecialShiftUpsilonLabelx=-3;
							SpecialShiftUpsilonLabel=1.; //14;

							char texTexMP[200];
							if(iStateMP==1) sprintf(texTexMP,"J/#psi");
							if(iStateMP==2) sprintf(texTexMP,"#psi(2S)");
							xRapText_MPnew = 44 + 20;
							if(iStateMP==2) xRapText_MPnew = 44 + 16;

							TLatex *textMP = new TLatex(xRapText_MPnew+SpecialShiftUpsilonLabelx,yMin+(yMax-yMin)*yRapText*SpecialShiftUpsilonLabel,texTexMP);
							textMP->SetTextSize(textSizeRap*increaseSize);
							if(iPanel==nPanels) textMP->SetTextSize(textSizeRap*(1-lowestBottomMargin)*increaseSize);
							textMP->Draw( "same" );
							xRapText_MPnew = 44;

							char abcdef[200];
							if(iStateMP==1&&iPanel==1) sprintf(abcdef,"a)");
							if(iStateMP==1&&iPanel==2) sprintf(abcdef,"b)");
							if(iStateMP==1&&iPanel==3) sprintf(abcdef,"c)");
							if(iStateMP==2&&iPanel==1) sprintf(abcdef,"d)");
							if(iStateMP==2&&iPanel==2) sprintf(abcdef,"e)");
							if(iStateMP==2&&iPanel==3) sprintf(abcdef,"f)");
							if(iStateMP==3&&iPanel==1) sprintf(abcdef,"g)");
							if(iStateMP==3&&iPanel==2) sprintf(abcdef,"h)");
							if(iStateMP==3&&iPanel==3) sprintf(abcdef,"i)");
							cout<<abcdef<<endl;
							TLatex *tex_abcdef = new TLatex(xabcdefText,yMin+(yMax-yMin)*yRapText*0.92,abcdef);
							tex_abcdef->SetTextSize(textSizeRap);
							if(iPanel==nPanels) tex_abcdef->SetTextSize(textSizeRap*(1-lowestBottomMargin));
							//tex_abcdef->Draw( "same" );

							if(PlotFinalData&&DrawLatexStuff){

								TLine* extreme0MP = new TLine( PlotpTMin-DeltaXminOVERALL, 0, PlotpTMax ,0);
								extreme0MP->SetLineWidth( 1 );
								extreme0MP->SetLineStyle( 2 );
								extreme0MP->SetLineColor( kBlack );
								extreme0MP->Draw( "same" );

								TLine* extreme1MP = new TLine( PlotpTMin-DeltaXminOVERALL, 1, PlotpTMax, 1);
								extreme1MP->SetLineWidth( 1 );
								extreme1MP->SetLineStyle( 2 );
								extreme1MP->SetLineColor( kBlack );
								if(iLam==1||iLam==7||iLam==13) extreme1MP->Draw( "same" );

								TLine* extreme2MP = new TLine( PlotpTMin-DeltaXminOVERALL, -1, PlotpTMax ,-1);
								extreme2MP->SetLineWidth( 1 );
								extreme2MP->SetLineStyle( 2 );
								extreme2MP->SetLineColor( kBlack );
								if(iLam==1||iLam==7||iLam==13) extreme2MP->Draw( "same" );
								if(iLam==6||iLam==12||iLam==18) extreme2MP->Draw( "same" );
							}


							if(iStateMP==1&&iPanel==1){
								cout<<"DRAW CMS preliminary Latex"<<endl;
								char text[200];
								sprintf(text,"CMS  pp  #sqrt{s} = 7 TeV   L = 4.9 fb^{-1}");
								double MPlatexX_buffer = 11.5;
								TLatex *CentralsText1MP = new TLatex(MPlatexX_buffer,MPlatexYmax,text);
								CentralsText1MP->SetTextSize(CentralsFontSizeMP);
								CentralsText1MP->Draw( "same" );

								sprintf(text,"Preliminary");
								TLatex *CentralsText2MP = new TLatex(PlotpTMax-22.,MPlatexYmax-1.5*MPlatexDeltaYmax,text);
								CentralsText2MP->SetTextSize(CentralsFontSizeMP);
								CentralsText2MP->Draw( "same" );

							}

							//if(iStateMP==2&&iPanel==3&&MPframe==1){
							//	MPframedepLegend->Draw("same");
							//}

							if(iStateMP==2&&iPanel==2&&MPframe>=1){
								MPframedepLegend->Draw("same");
							}

							if(iStateMP==2&&iPanel==1){

								char frameMPtex[200];
								if(MPframe==1) sprintf(frameMPtex,"CS frame");
								if(MPframe==2) sprintf(frameMPtex,"HX frame");
								if(MPframe==3) sprintf(frameMPtex,"PX frame");
								char textStateFrame[200];
								sprintf(textStateFrame,"#psi(%dS), %s", nState-3, frameMPtex);
								if(nState==4) sprintf(textStateFrame,"J/#psi, %s", frameMPtex);
								TLatex *TexStateFrame = new TLatex(MPlatexX,MPlatexYmax,textStateFrame);
								TexStateFrame->SetTextSize(CentralsFontSizeMP);
								//TexStateFrame->Draw( "same" );

							}

						}//end iState Loop

						if(MPframe==1) MPcanvasCS_New->cd();
						if(MPframe==2) MPcanvasHX_New->cd();
						if(MPframe==3) MPcanvasPX_New->cd();

						whereTexteInPlotX=0.488;
						whereTexteInPlotY=startValCoordY-deltaCoordY-1.375*labelOffsetX;

						TLatex *MPXlabeltext = new TLatex(whereTexteInPlotX,whereTexteInPlotY ,"1");
						MPXlabeltext->SetTextSize(XaxislabelLatexSize);
						MPXlabeltext->SetTextColor(kBlack);
						//if(iPanel==nPanels) MPXlabeltext->Draw( "same" );

						if(iLam==3){  
							sprintf(filename,"%s/FinalResultsCS.pdf",FigDir);
							if(PlotFinalData) MPcanvasCS_New->SaveAs(filename);
							sprintf(filename,"%s/FinalResultsCS.C",FigDir);
							if(PlotFinalData) MPcanvasCS_New->SaveAs(filename);
							MPcanvasCS_New->Close();
						}
						if(iLam==9){
							sprintf(filename,"%s/FinalResultsHX.pdf",FigDir);
							if(PlotFinalData) MPcanvasHX_New->SaveAs(filename);
							sprintf(filename,"%s/FinalResultsHX.C",FigDir);
							if(PlotFinalData) MPcanvasHX_New->SaveAs(filename);
							MPcanvasHX_New->Close();
						}
						if(iLam==15){
							sprintf(filename,"%s/FinalResultsPX.pdf",FigDir);
							if(PlotFinalData) MPcanvasPX_New->SaveAs(filename);
							sprintf(filename,"%s/FinalResultsPX.C",FigDir);
							if(PlotFinalData) MPcanvasPX_New->SaveAs(filename);
							MPcanvasPX_New->Close();
						}

					}//end Frame dependent plots
				}


				cout<<"begin Frame independent plots -- NEW"<<endl;
				cout<<"if(iLam==6||iLam==12||iLam==18)"<<endl;

				if(Psi_MPplots_New){
					lowestBottomMargin =  lowestBottomMargin;

					MPcanvasXpixel = MPcanvasXpixelInitial;
					MPcanvasYpixel = MPcanvasYpixelInitial  / ( (2.+1./(1-lowestBottomMargin)) * (1.-lowestBottomMargin) ); // *0.384 
					cout<<"MPcanvasXpixel: "<<MPcanvasXpixel<<" MPcanvasYpixel: "<<MPcanvasYpixel<<endl;

					if( (iLam==6||iLam==12||iLam==18 ) && rapBin<3 ){

						int mainframe;
						if(iLam==6) mainframe=1;
						if(iLam==12) mainframe=2;
						if(iLam==18) mainframe=3;

						cout<<"iLam = "<<iLam<<endl;

						if((iLam==6||iLam==12||iLam==18)&&rapBin==1){
							MPcanvasTilde_rap1 = new TCanvas("MPcanvasTilde_rap1", "MPcanvasTilde_rap1",MPcanvasXpixel,MPcanvasYpixel);
							MPcanvasTilde_rap1->SetFillColor(kWhite);
							MPcanvasTilde_rap1->GetFrame()->SetFillColor(kWhite);
							MPcanvasTilde_rap1->GetFrame()->SetBorderSize(0);
						}
						if((iLam==6||iLam==12||iLam==18)&&rapBin==2){
							MPcanvasTilde_rap2 = new TCanvas("MPcanvasTilde_rap2", "MPcanvasTilde_rap2",MPcanvasXpixel,MPcanvasYpixel);
							MPcanvasTilde_rap2->SetFillColor(kWhite);
							MPcanvasTilde_rap2->GetFrame()->SetFillColor(kWhite);
							MPcanvasTilde_rap2->GetFrame()->SetBorderSize(0);
						}

						cout<<"MultiPanel canvas"<<endl;

						int totalState = 2;
						for(int iStateMP=1;iStateMP<totalState+1;iStateMP++){ //

							PlotpTMin = 12.; PlotpTMax = 72.;

							if((iLam==6||iLam==12||iLam==18)&&rapBin==1) MPcanvasTilde_rap1 -> cd() ;
							if((iLam==6||iLam==12||iLam==18)&&rapBin==2) MPcanvasTilde_rap2 -> cd() ;

							cout<<"MultiPanel pad"<<endl;
							cout<<"iLam "<<iLam<<endl;
							cout<<"mainframe "<<mainframe<<endl;
							cout<<"iStateMP "<<iStateMP<<endl;

							TPad *MPpad;
							double X1panel = 0.373;
							double deltaCoordX = 0.;
							deltaCoordX = PadCoordXMax/(1./(1-Left_margin) + 1./(1-Right_margin));
							X1panel = deltaCoordX/(1.-Left_margin);
							if(iStateMP==1) MPpad = new TPad("MPpad","MPpad",0., 0.,X1panel,1.);
							if(iStateMP==2) MPpad = new TPad("MPpad","MPpad",X1panel,0.,PadCoordXMax, 1.);

							MPpad->Draw();
							MPpad->cd();
							MPpad->SetFillColor(kWhite);
							MPpad->SetFrameFillColor(kWhite);
							MPpad->SetBorderSize(0);
							MPpad->SetLeftMargin(0.);
							if(iStateMP==1) MPpad->SetLeftMargin(Left_margin);
							MPpad->SetRightMargin(0.);
							if(iStateMP==2) MPpad->SetRightMargin(Right_margin); //+0.05
							MPpad->SetTopMargin(Top_margin+0.025);
							MPpad->SetBottomMargin(0.0);
							MPpad->SetBottomMargin(lowestBottomMargin);

							cout<<"MultiPanel hist"<<endl;
							TH1F *MPhist = new TH1F;
							MPhist = MPcanvasTilde_rap1->DrawFrame(PlotpTMin-DeltaXminOVERALL,yMinMP,PlotpTMax,yMaxMP);

							MPhist->SetXTitle("#it{p}_{T} [GeV]");
							MPhist->GetXaxis()->SetTitleOffset(-1.35); //-1.2

							MPhist->SetYTitle(axislabel);
							MPhist->GetYaxis()->SetTitleOffset(titleoffset*0.2);
							MPhist->GetYaxis()->SetTitleSize(0.);
							if(iStateMP==1) MPhist->GetYaxis()->SetTitleOffset(titleoffset*1.35);
							if(iStateMP==1) MPhist->GetYaxis()->SetTitleSize(0.*(1-lowestBottomMargin));

							MPhist->GetYaxis()->SetLabelSize(LabelSize*1.25); //0.5
							MPhist->GetXaxis()->SetLabelSize(0.);
							MPhist->GetYaxis()->SetLabelOffset(-0.015);
							MPhist->GetXaxis()->SetLabelOffset(-0.08);
							//if(rapBin==1) MPhist->GetXaxis()->SetLabelOffset(-0.08*(1-Left_margin));

							if(iStateMP==1) MPhist->GetYaxis()->SetLabelSize(LabelSize*1.25*(1-lowestBottomMargin));
							MPhist->GetXaxis()->SetTitleSize(TitleSize*0.85); //0.55
							MPhist->GetXaxis()->SetAxisColor(kWhite);
							MPhist->GetYaxis()->SetAxisColor(kWhite);
							MPhist->GetXaxis()->SetTicks("-");
							MPhist->GetYaxis()->SetTicks("+");

							TLegend* MPframedepLegendError;
							MPframedepLegendError=new TLegend(errorLegendX1-0.25,errorLegendY1+0.1,errorLegendX2,errorLegendY2);
							//if(nState==5) MPframedepLegendError=new TLegend(errorLegendX1-0.25,errorLegendY1+0.05,errorLegendX2-0.2,errorLegendY2); //May2
							MPframedepLegendError->SetFillColor(0);
							//MPframedepLegendError->SetTextFont(72);
							MPframedepLegendError->SetTextSize(errorLegendFontSize*0.7);
							MPframedepLegendError->SetBorderSize(0);

							char MPframedepLegendEntry[200];


							TGraphAsymmErrors* graphMP_1sig;

							TLegend* MPtildeLegend;
							MPtildeLegend=new TLegend(0.7,0.75,0.9,0.95);
							MPtildeLegend->SetFillColor(0);
							//MPtildeLegend->SetTextFont(72);
							MPtildeLegend->SetTextSize(0.05);
							MPtildeLegend->SetBorderSize(0);
							char MPtildeLegendEntry[200];

							for(int iFrameMP=1;iFrameMP<4;iFrameMP++){

								char GraphNameMP[200];

								if(iFrameMP==1){
									sprintf(GraphNameMP,"ltilde_CS_rap%d",rapBin);
								}
								if(iFrameMP==2){
									sprintf(GraphNameMP,"ltilde_HX_rap%d",rapBin);
								}
								if(iFrameMP==3){
									sprintf(GraphNameMP,"ltilde_PX_rap%d",rapBin);
								}

								if(iStateMP==1) graphMP_1sig = (TGraphAsymmErrors*) infileMP1_1sig->Get(GraphNameMP);
								if(iStateMP==2) graphMP_1sig = (TGraphAsymmErrors*) infileMP2_1sig->Get(GraphNameMP);

								int MarkerDefinitionForThisBin[4][4]={{0,0,0,0},{0,1,2,3},{0,2,1,3},{0,2,3,1}};

								int ptBinMinMP=3, ptBinMaxMP=12;
								if(iStateMP==2) { ptBinMinMP=2; ptBinMaxMP=5; }
								int nBinspTMP = ptBinMaxMP-ptBinMinMP+1;

								double ptCentreMP[nBinspTMP];
								double ptCentreErr_lowMP[nBinspTMP];
								double ptCentreErr_highMP[nBinspTMP];
								double lmeanMP[nBinspTMP];
								double lmean_errlowMP[nBinspTMP];
								double lmean_errhighMP[nBinspTMP];

								double ShiftTildePlot;
								double ShiftTildePlotZero=ColordBandWidth/2.+0.2; //0.5;
								if(iStateMP==1) ShiftTildePlotZero=ColordBandWidth/2.;
								//if(iStateMP==2) ShiftTildePlotZero = 0.7;
								bool RemoveHorizontalErrorBar=true;

								if(MarkerDefinitionForThisBin[mainframe][iFrameMP]==1) ShiftTildePlot=0.;
								if(MarkerDefinitionForThisBin[mainframe][iFrameMP]==2) ShiftTildePlot=ShiftTildePlotZero;
								if(MarkerDefinitionForThisBin[mainframe][iFrameMP]==3) ShiftTildePlot=-ShiftTildePlotZero;

								int pt=0;
								for(int ptBin = ptBinMinMP; ptBin < ptBinMaxMP+1; ptBin++) {

									graphMP_1sig->GetPoint(ptBin-1,ptCentreMP[pt],lmeanMP[pt]);
									ptCentreErr_highMP[pt]=graphMP_1sig->GetErrorXhigh(ptBin-1);
									ptCentreErr_lowMP[pt]=graphMP_1sig->GetErrorXlow(ptBin-1);
									lmean_errhighMP[pt]=graphMP_1sig->GetErrorYhigh(ptBin-1);
									lmean_errlowMP[pt]=graphMP_1sig->GetErrorYlow(ptBin-1);

									ptCentreMP[pt]=ptCentreMP[pt]+ShiftTildePlot;
									ptCentreErr_highMP[pt]=ptCentreErr_highMP[pt]-ShiftTildePlot;
									ptCentreErr_lowMP[pt]=ptCentreErr_lowMP[pt]-ShiftTildePlot;
									if(RemoveHorizontalErrorBar){ ptCentreErr_highMP[pt]=0.; ptCentreErr_lowMP[pt]=0;}

									pt++;
								}
								graphMP_1sig = new TGraphAsymmErrors(nBinspTMP,ptCentreMP,lmeanMP,ptCentreErr_lowMP,ptCentreErr_highMP,lmean_errlowMP,lmean_errhighMP);

								graphMP_1sig->SetMarkerColor(MarkerColorMP[MarkerDefinitionForThisBin[mainframe][iFrameMP]]);
								graphMP_1sig->SetLineColor(MarkerColorMP[MarkerDefinitionForThisBin[mainframe][iFrameMP]]);
								graphMP_1sig->SetMarkerStyle(MarkerStyleMP[MarkerDefinitionForThisBin[mainframe][iFrameMP]]);
								graphMP_1sig->SetMarkerSize(MarkerSizeMP[MarkerDefinitionForThisBin[mainframe][iFrameMP]]);

								if(mainframe==1&&iFrameMP==1){ 
									sprintf(MPtildeLegendEntry,"CS"); MPtildeLegend->AddEntry(graphMP_1sig,MPtildeLegendEntry,"p");
									//sprintf(MPtildeLegendEntry,"CS"); MPtildeLegend->AddEntry(graphMP_1sig,MPtildeLegendEntry,"lp");
								}
								if(mainframe!=1&&iFrameMP==1){ 
									sprintf(MPtildeLegendEntry,"CS"); MPtildeLegend->AddEntry(graphMP_1sig,MPtildeLegendEntry,"p"); 
								}

								if(mainframe==2&&iFrameMP==2){ 
									sprintf(MPtildeLegendEntry,"HX"); MPtildeLegend->AddEntry(graphMP_1sig,MPtildeLegendEntry,"p");
									//sprintf(MPtildeLegendEntry,"HX"); MPtildeLegend->AddEntry(graphMP_1sig,MPtildeLegendEntry,"lp"); 
								}
								if(mainframe!=2&&iFrameMP==2){ 
									sprintf(MPtildeLegendEntry,"HX"); MPtildeLegend->AddEntry(graphMP_1sig,MPtildeLegendEntry,"p");
								}

								if(mainframe==3&&iFrameMP==3){ 
									sprintf(MPtildeLegendEntry,"PX"); MPtildeLegend->AddEntry(graphMP_1sig,MPtildeLegendEntry,"p");
									//sprintf(MPtildeLegendEntry,"PX"); MPtildeLegend->AddEntry(graphMP_1sig,MPtildeLegendEntry,"lp"); 
								}

								if(mainframe!=3&&iFrameMP==3){ 
									sprintf(MPtildeLegendEntry,"PX"); MPtildeLegend->AddEntry(graphMP_1sig,MPtildeLegendEntry,"p");
								}

								if(iFrameMP==1){ if(mainframe==1) graphMP_1sig->Draw(drawGraphStyle); }
								if(iFrameMP==2){ if(mainframe==2) graphMP_1sig->Draw(drawGraphStyle); }
								if(iFrameMP==3){ if(mainframe==3) graphMP_1sig->Draw(drawGraphStyle); }

								if(iFrameMP==1){ if(mainframe!=1) graphMP_1sig->Draw("PX"); }
								if(iFrameMP==2){ if(mainframe!=2) graphMP_1sig->Draw("PX"); }
								if(iFrameMP==3){ if(mainframe!=3) graphMP_1sig->Draw("PX"); }

							}

							TGaxis *axisMPY1 = new TGaxis(PlotpTMin-DeltaXminOVERALL,yMinMP,PlotpTMin-DeltaXminOVERALL,yMaxMP,
									yMinMP,yMaxMP,AxisDivisions,"-US");
							axisMPY1->SetTickSize(ticksize);
							axisMPY1->SetTickSize(ticksize/(1-lowestBottomMargin));
							axisMPY1->Draw("same");

							TGaxis *axisMPY2 = new TGaxis(PlotpTMax,yMinMP,PlotpTMax,yMaxMP,yMinMP,yMaxMP,AxisDivisions,"+US");
							axisMPY2->SetTickSize(ticksize);
							axisMPY2->SetTickSize(ticksize/(1-lowestBottomMargin));
							axisMPY2->Draw("same");

							TGaxis *axisM3S1 = new TGaxis(PlotpTMin-DeltaXminOVERALL,yMinMP,PlotpTMax,yMinMP,
									PlotpTMin-DeltaXminOVERALL+deltaTrickAxisMin,PlotpTMax+deltaTrickAxisMax,AxisDivisions,"+S");
							axisM3S1->SetTickSize(ticksize*2);
							axisM3S1->SetLabelSize(LabelSize); //*0.62
							//if(rapBin==1) axisM3S1->SetLabelSize(LabelSize*(1-Left_margin)); //*0.72
							axisM3S1->SetLabelOffset(labelOffsetX);//*0.7
							if(iStateMP==1) axisM3S1->SetLabelOffset(labelOffsetX/(1-Left_margin)); //*0.8
							axisM3S1->SetTickSize(ticksize*0.7/(1-lowestBottomMargin));
							axisM3S1->Draw("same");

							TGaxis *axisM3S2 = new TGaxis(PlotpTMin-DeltaXminOVERALL,yMaxMP,PlotpTMax,yMaxMP,
									PlotpTMin-DeltaXminOVERALL+deltaTrickAxisMin,PlotpTMax+deltaTrickAxisMax,AxisDivisions,"-US");
							axisM3S2->SetTickSize(ticksize*2);
							axisM3S2->SetTickSize(ticksize/(1-lowestBottomMargin));
							axisM3S2->Draw("same");

							whereTexteInPlotX=2;//XtitlePosition;
							whereTexteInPlotY=(yMaxMP+yMinMP)/2.-(yMaxMP-yMinMP)*XtitlePositionYshift;

							char axistitleMPtilde[200];
							sprintf(axistitleMPtilde,"#tilde{#lambda}");
							YaxistitleLatexSize=YaxistitleLatexSize*(1-lowestBottomMargin);//*0.7
							TLatex *MPYtitletext = new TLatex(whereTexteInPlotX,whereTexteInPlotY ,axistitleMPtilde);
							MPYtitletext->SetTextSize(YaxistitleLatexSize);
							MPYtitletext->SetTextSize(YaxistitleLatexSize*(1-lowestBottomMargin));
							MPYtitletext->SetTextColor(kBlack);
							MPYtitletext->PaintLatex(whereTexteInPlotX,whereTexteInPlotY, YtitleAngle, YaxistitleLatexSize, axislabel);
							MPYtitletext->Draw( "same" );

							char texTexMP[200];
							if(iStateMP==1&&rapBin==1) sprintf(texTexMP,"J/#psi, |#it{y}| < 0.6", nState-3);
							if(iStateMP==1&&rapBin==2) sprintf(texTexMP,"J/#psi, 0.6 < |#it{y}| < 1.2", nState-3);
							if(iStateMP==2&&rapBin==1) sprintf(texTexMP,"#psi(2S), |#it{y}| < 0.6", nState-3);
							if(iStateMP==2&&rapBin==2) sprintf(texTexMP,"#psi(2S), 0.6 < |#it{y}| < 1.2", nState-3);

							xRapTextTilde = PlotpTMax * 0.66 ;
							if(iStateMP==1) xRapTextTilde = PlotpTMax * 0.71;
							if(rapBin==2) {
								xRapTextTilde = PlotpTMax * 0.55;
								if(iStateMP==1) xRapTextTilde = PlotpTMax * 0.6;
							}

							TLatex *textMP = new TLatex(xRapTextTilde,yMin+(yMax-yMin)*yRapText,texTexMP);//*0.92
							textMP->SetTextSize(textSizeRap*0.8);//*0.7
							//if(rapBin==1) textMP->SetTextSize(textSizeRap*(1-Left_margin));//*0.75
							textMP->Draw( "same" );

							char abcdef[200];
							if(iStateMP==1) sprintf(abcdef,"a)");
							if(iStateMP==2) sprintf(abcdef,"b)");
							if(iStateMP==3) sprintf(abcdef,"c)");
							cout<<abcdef<<endl;
							TLatex *tex_abcdef = new TLatex(xabcdefText,yMin+(yMax-yMin)*yRapText*0.92,abcdef);
							tex_abcdef->SetTextSize(textSizeRap);
							tex_abcdef->SetTextSize(textSizeRap*(1-lowestBottomMargin));
							//tex_abcdef->Draw( "same" );

							if(PlotFinalData&&DrawLatexStuff){

								TLine* extreme0MP = new TLine( PlotpTMin-DeltaXminOVERALL, 0, PlotpTMax ,0);
								extreme0MP->SetLineWidth( 1 );
								extreme0MP->SetLineStyle( 2 );
								extreme0MP->SetLineColor( kBlack );
								extreme0MP->Draw( "same" );

								TLine* extreme1MP = new TLine( PlotpTMin-DeltaXminOVERALL, 1, PlotpTMax , 1);
								extreme1MP->SetLineWidth( 1 );
								extreme1MP->SetLineStyle( 2 );
								extreme1MP->SetLineColor( kBlack );
								if(iLam==1||iLam==7||iLam==13) extreme1MP->Draw( "same" );

								TLine* extreme2MP = new TLine( PlotpTMin-DeltaXminOVERALL, -1, PlotpTMax ,-1);
								extreme2MP->SetLineWidth( 1 );
								extreme2MP->SetLineStyle( 2 );
								extreme2MP->SetLineColor( kBlack );
								if(iLam==1||iLam==7||iLam==13) extreme2MP->Draw( "same" );
								//if(iLam==6||iLam==12||iLam==18) extreme2MP->Draw( "same" );

							}


							if(iStateMP==1){
								cout<<"DRAW CMS preliminary Latex"<<endl;
								char text[200];
								sprintf(text,"CMS  pp  #sqrt{s} = 7 TeV   L = 4.9 fb^{-1}");
								TLatex *CentralsText1MP = new TLatex(MPlatexX,MPlatexYmax,text);
								CentralsText1MP->SetTextSize(CentralsFontSizeMP*0.75);
								CentralsText1MP->Draw( "same" );

								sprintf(text,"Preliminary");
								TLatex *CentralsText2MP = new TLatex(PlotpTMax-22,MPlatexYmax-1.5*MPlatexDeltaYmax,text);
								if(nState==5) CentralsText2MP = new TLatex(PlotpTMax-13.,MPlatexYmax-1.5*MPlatexDeltaYmax,text); //May2
								CentralsText2MP->SetTextSize(CentralsFontSizeMP*0.75);
								CentralsText2MP->Draw( "same" );

							}

							if(iStateMP==2) MPtildeLegend->Draw("same");

						}//end iState Loop

						if(rapBin==1) MPcanvasTilde_rap1->cd();
						if(rapBin==2) MPcanvasTilde_rap2->cd();

						whereTexteInPlotX=0.488;
						whereTexteInPlotY=startValCoordY-deltaCoordY-1.425*labelOffsetX;

						TLatex *M3Slabeltext = new TLatex(whereTexteInPlotX,whereTexteInPlotY ,"1");
						M3Slabeltext->SetTextSize(XaxislabelLatexSize);
						M3Slabeltext->SetTextColor(kBlack);
						//M3Slabeltext->Draw( "same" );

						if(iLam==6||iLam==12||iLam==18){
							if(mainframe==1&&rapBin==1) sprintf(filename,"%s/FinalResultsTildeCS_rap1.pdf",FigDir);
							if(mainframe==2&&rapBin==1) sprintf(filename,"%s/FinalResultsTildeHX_rap1.pdf",FigDir);
							if(mainframe==3&&rapBin==1) sprintf(filename,"%s/FinalResultsTildePX_rap1.pdf",FigDir);
							if(mainframe==1&&rapBin==2) sprintf(filename,"%s/FinalResultsTildeCS_rap2.pdf",FigDir);
							if(mainframe==2&&rapBin==2) sprintf(filename,"%s/FinalResultsTildeHX_rap2.pdf",FigDir);
							if(mainframe==3&&rapBin==2) sprintf(filename,"%s/FinalResultsTildePX_rap2.pdf",FigDir);
							if(PlotFinalData) { 
								if(rapBin==1) MPcanvasTilde_rap1->SaveAs(filename); 
								if(rapBin==2) MPcanvasTilde_rap2->SaveAs(filename); 
							}

							if(mainframe==1&&rapBin==1) sprintf(filename,"%s/FinalResultsTildeCS_rap1.C",FigDir);
							if(mainframe==2&&rapBin==1) sprintf(filename,"%s/FinalResultsTildeHX_rap1.C",FigDir);
							if(mainframe==3&&rapBin==1) sprintf(filename,"%s/FinalResultsTildePX_rap1.C",FigDir);
							if(mainframe==1&&rapBin==2) sprintf(filename,"%s/FinalResultsTildeCS_rap2.C",FigDir);
							if(mainframe==2&&rapBin==2) sprintf(filename,"%s/FinalResultsTildeHX_rap2.C",FigDir);
							if(mainframe==3&&rapBin==2) sprintf(filename,"%s/FinalResultsTildePX_rap2.C",FigDir);
							if(PlotFinalData) { 
								if(rapBin==1) MPcanvasTilde_rap1->SaveAs(filename); 
								if(rapBin==2) MPcanvasTilde_rap2->SaveAs(filename); 
							}

							if(rapBin==1)MPcanvasTilde_rap1->Close();
							if(rapBin==2)MPcanvasTilde_rap2->Close();
						}


					}//end Frame independent plots
				}

				//======================================================================================



				////////////////////////// OLD PLOT /////////////////////////////////////////
				cout<<"begin Frame dependent plots"<<endl;

				//begin Frame dependent plots
				if(Psi_MPplots_Old){
					if(nState>3) {
						MPcanvasXpixel = MPcanvasXpixelInitial * 1.5;
						if(nState==4)
							MPcanvasYpixel = MPcanvasYpixelInitial * 0.65;
					}
					cout<<"MPcanvasXpixel: "<<MPcanvasXpixel<<" MPcanvasYpixel: "<<MPcanvasYpixel<<endl;

					if(iLam!=6&&iLam!=12&&iLam!=18){

						if(iLam>0&&iLam<7) MPframe=1;
						if(iLam>6&&iLam<13) MPframe=2;
						if(iLam>12&&iLam<19) MPframe=3;

						if(iLam==1||iLam==7||iLam==13) iPanel=1;
						if(iLam==2||iLam==8||iLam==14) iPanel=2;
						if(iLam==3||iLam==9||iLam==15) iPanel=3;

						cout<<"MultiPanel canvas"<<endl;

						if(iLam==1&&rapBin==1){
							MPcanvasCS_Psi = new TCanvas("MPcanvasCS_Psi", "MPcanvasCS_Psi",MPcanvasXpixel,MPcanvasYpixel);
							MPcanvasCS_Psi->SetFillColor(kWhite);
							MPcanvasCS_Psi->GetFrame()->SetFillColor(kWhite);
							MPcanvasCS_Psi->GetFrame()->SetBorderSize(0);
						}
						if(iLam==7&&rapBin==1){
							MPcanvasHX_Psi = new TCanvas("MPcanvasHX_Psi", "MPcanvasHX_Psi",MPcanvasXpixel,MPcanvasYpixel);
							MPcanvasHX_Psi->SetFillColor(kWhite);
							MPcanvasHX_Psi->GetFrame()->SetFillColor(kWhite);
							MPcanvasHX_Psi->GetFrame()->SetBorderSize(0);
						}
						if(iLam==13&&rapBin==1){
							MPcanvasPX_Psi = new TCanvas("MPcanvasPX_Psi", "MPcanvasPX_Psi",MPcanvasXpixel,MPcanvasYpixel);
							MPcanvasPX_Psi->SetFillColor(kWhite);
							MPcanvasPX_Psi->GetFrame()->SetFillColor(kWhite);
							MPcanvasPX_Psi->GetFrame()->SetBorderSize(0);
						}

						if(MPframe==1) MPcanvasCS_Psi->cd();
						if(MPframe==2) MPcanvasHX_Psi->cd();
						if(MPframe==3) MPcanvasPX_Psi->cd();

						cout<<"MultiPanel pad"<<endl;
						TPad *MPpad_PsiDep;
						if(nState==5){
							if(rapBin==1) MPpad_PsiDep = new TPad("MPpad_PsiDep","MPpad_PsiDep",PadCoordY[iPanel-1],PadCoordY[2],PadCoordY[iPanel],PadCoordY[3]);
							if(rapBin==2) MPpad_PsiDep = new TPad("MPpad_PsiDep","MPpad_PsiDep",PadCoordY[iPanel-1],PadCoordY[1],PadCoordY[iPanel],PadCoordY[2]);
							if(rapBin==3) MPpad_PsiDep = new TPad("MPpad_PsiDep","MPpad_PsiDep",PadCoordY[iPanel-1],PadCoordY[0],PadCoordY[iPanel],PadCoordY[1]);
						}
						else{
							if(rapBin==1) MPpad_PsiDep = new TPad("MPpad_PsiDep","MPpad_PsiDep",PadCoordY[iPanel-1] ,0.55, PadCoordY[iPanel] ,1.);
							if(rapBin==2) MPpad_PsiDep = new TPad("MPpad_PsiDep","MPpad_PsiDep",PadCoordY[iPanel-1] ,0.0, PadCoordY[iPanel] ,0.55);
						}

						MPpad_PsiDep->Draw();
						MPpad_PsiDep->cd();

						MPpad_PsiDep->SetFillColor(kWhite);
						MPpad_PsiDep->SetFrameFillColor(kWhite);
						MPpad_PsiDep->SetBorderSize(0);
						MPpad_PsiDep->SetLeftMargin(0.);
						MPpad_PsiDep->SetRightMargin(0.);
						if(iPanel==1) MPpad_PsiDep->SetLeftMargin(Left_margin);
						if(iPanel==3) MPpad_PsiDep->SetRightMargin(Right_margin);
						MPpad_PsiDep->SetTopMargin(Top_margin+0.0025);
						if(rapBin==1) MPpad_PsiDep->SetTopMargin(Top_margin+0.035);
						MPpad_PsiDep->SetBottomMargin(0.0);
						if((nState==4&&rapBin==2)||(nState==5&&rapBin==3)) MPpad_PsiDep->SetBottomMargin(lowestBottomMargin);

						cout<<"MultiPanel hist"<<endl;
						TH1F *MPhist_PsiDep = new TH1F;
						if(MPframe==1)MPhist_PsiDep = MPcanvasCS_Psi->DrawFrame(PlotpTMin-DeltaXminOVERALL,yMinMP,PlotpTMax,yMaxMP);
						if(MPframe==2)MPhist_PsiDep = MPcanvasHX_Psi->DrawFrame(PlotpTMin-DeltaXminOVERALL,yMinMP,PlotpTMax,yMaxMP);
						if(MPframe==3)MPhist_PsiDep = MPcanvasPX_Psi->DrawFrame(PlotpTMin-DeltaXminOVERALL,yMinMP,PlotpTMax,yMaxMP);

						MPhist_PsiDep->SetXTitle("#it{p}_{T} [GeV]");
						MPhist_PsiDep->GetXaxis()->SetTitleOffset(-1.35);

						MPhist_PsiDep->SetYTitle(axislabel);
						MPhist_PsiDep->GetYaxis()->SetTitleOffset(titleoffset);
						MPhist_PsiDep->GetYaxis()->SetTitleSize(0.);

						if((nState==4&&rapBin==2)||(nState==5&&rapBin==3)) 
							MPhist_PsiDep->GetYaxis()->SetTitleOffset(titleoffset*1.35);
						if((nState==4&&rapBin==2)||(nState==5&&rapBin==3)) 
							MPhist_PsiDep->GetYaxis()->SetTitleSize(0.*(1-lowestBottomMargin));

						MPhist_PsiDep->GetYaxis()->SetLabelSize(LabelSize*1.25);
						MPhist_PsiDep->GetXaxis()->SetLabelSize(0.);
						MPhist_PsiDep->GetYaxis()->SetLabelOffset(-0.015);
						if(iPanel>1) MPhist_PsiDep->GetYaxis()->SetLabelOffset(0.08);
						MPhist_PsiDep->GetXaxis()->SetLabelOffset(-0.06);

						if((nState==4&&rapBin==2)||(nState==5&&rapBin==3)) 
							MPhist_PsiDep->GetYaxis()->SetLabelSize(LabelSize*1.25*(1-lowestBottomMargin));
						MPhist_PsiDep->GetXaxis()->SetTitleSize(TitleSize*0.85);
						MPhist_PsiDep->GetXaxis()->SetAxisColor(kWhite);
						MPhist_PsiDep->GetYaxis()->SetAxisColor(kWhite);
						MPhist_PsiDep->GetXaxis()->SetTicks("-");
						MPhist_PsiDep->GetYaxis()->SetTicks("+");


						TLegend* MPframedepLegend;
						MPframedepLegend=new TLegend(errorLegendX1,errorLegendY1,errorLegendX2,errorLegendY2);
						if((nState==2||nState==3)&&MPframe==1) 
							MPframedepLegend=new TLegend(errorLegendX1, errorLegendY2-(errorLegendY2-errorLegendY1)*(1-lowestBottomMargin),
									errorLegendX2,errorLegendY2);

						MPframedepLegend->SetFillColor(0);
						//MPframedepLegend->SetTextFont(72);
						MPframedepLegend->SetTextSize(errorLegendFontSize);
						if((nState==2||nState==3)&&MPframe==1) MPframedepLegend->SetTextSize(errorLegendFontSize*(1-lowestBottomMargin));
						MPframedepLegend->SetBorderSize(0);

						char MPframedepLegendEntry[200];

						graphSyst->SetMarkerSize(MarkerSizeMP[0]);
						graphSyst->SetMarkerStyle(MarkerStyleMP[0]);
						graphSyst->SetMarkerColor(MarkerColorMP[0]);

						graphDefaultStat->SetMarkerSize(MarkerSizeMP[0]);
						graphDefaultStat->SetMarkerStyle(MarkerStyleMP[0]);
						graphDefaultStat->SetMarkerColor(MarkerColorMP[0]);

						graphDefaultRes3sigma->Draw("2");
						graphDefaultRes2sigma->Draw("2");
						graphDefaultRes->Draw("2");
						if(!PlotAlteredPPDResults) graphSyst->Draw(drawGraphStyle);
						if(PlotAlteredPPDResults) graphDefaultStat->Draw(drawGraphStyle);

						graphDefaultRes->SetLineColor(kGreen);
						graphDefaultRes2sigma->SetLineColor(kYellow);
						graphDefaultRes3sigma->SetLineColor(kCyan-9);

						if(PlotAlteredPPDResults) sprintf(MPframedepLegendEntry,"Stat. uncert., 68.3 %% CL");
						if(!PlotAlteredPPDResults) sprintf(MPframedepLegendEntry,"Tot. sys. uncert.");
						MPframedepLegend->AddEntry(graphSyst,MPframedepLegendEntry,"ple");
						if(PlotAlteredPPDResults) sprintf(MPframedepLegendEntry,"Tot. uncert., 68.3 %% CL");
						if(!PlotAlteredPPDResults) sprintf(MPframedepLegendEntry,"68.3 %% CL");
						MPframedepLegend->AddEntry(graphDefaultRes,MPframedepLegendEntry,"f");
						if(PlotAlteredPPDResults) sprintf(MPframedepLegendEntry,"Tot. uncert., 95.5 %% CL");
						if(!PlotAlteredPPDResults) sprintf(MPframedepLegendEntry,"95.5 %% CL");
						MPframedepLegend->AddEntry(graphDefaultRes2sigma,MPframedepLegendEntry,"f");
						if(PlotAlteredPPDResults) sprintf(MPframedepLegendEntry,"Tot. uncert., 99.7 %% CL");
						if(!PlotAlteredPPDResults) sprintf(MPframedepLegendEntry,"99.7 %% CL");
						MPframedepLegend->AddEntry(graphDefaultRes3sigma,MPframedepLegendEntry,"f");

						if(nState<=4 && rapBin==2) deltaTrickAxisMax=-0.001;
						if(nState==5 && rapBin==3) deltaTrickAxisMax=-0.001;

						TGaxis *axisMPY1 = new TGaxis(PlotpTMin-DeltaXminOVERALL,yMinMP,PlotpTMin-DeltaXminOVERALL,yMaxMP,
								yMinMP,yMaxMP,AxisDivisions,"-US");
						axisMPY1->SetTickSize(ticksize);
						if((nState==4&&rapBin==2)||(nState==5&&rapBin==3)) axisMPY1->SetTickSize(ticksize/(1-lowestBottomMargin));
						axisMPY1->Draw("same");

						TGaxis *axisMPY2 = new TGaxis(PlotpTMax,yMinMP,PlotpTMax,yMaxMP,yMinMP,yMaxMP,AxisDivisions,"+US");
						axisMPY2->SetTickSize(ticksize);
						if((nState==4&&rapBin==2)||(nState==5&&rapBin==3)) axisMPY2->SetTickSize(ticksize/(1-lowestBottomMargin));
						////axisMPY2->Draw("same");

						TGaxis *axisMPX1 = new TGaxis(PlotpTMin-DeltaXminOVERALL,yMinMP,PlotpTMax,yMinMP,
								PlotpTMin-DeltaXminOVERALL+deltaTrickAxisMin,PlotpTMax+deltaTrickAxisMax,AxisDivisions,"+S");
						axisMPX1->SetTickSize(ticksize*2);
						if((nState==4&&rapBin==2)||(nState==5&&rapBin==3)) axisMPX1->SetLabelSize(LabelSize);
						else axisMPX1->SetLabelSize(0);
						axisMPX1->SetLabelOffset(labelOffsetX);
						if((nState==4&&rapBin==2)||(nState==5&&rapBin==3)) axisMPX1->SetTickSize(ticksize/(1-lowestBottomMargin));
						axisMPX1->Draw("same");

						TGaxis *axisMPX2 = new TGaxis(PlotpTMin-DeltaXminOVERALL,yMaxMP,PlotpTMax,yMaxMP,
								PlotpTMin-DeltaXminOVERALL+deltaTrickAxisMin,PlotpTMax+deltaTrickAxisMax,AxisDivisions,"-US");
						axisMPX2->SetTickSize(ticksize*2);
						if((nState==4&&rapBin==2)||(nState==5&&rapBin==3)) axisMPX2->SetTickSize(ticksize/(1-lowestBottomMargin));
						axisMPX2->Draw("same");

						whereTexteInPlotX=XtitlePosition;
						whereTexteInPlotY=(yMaxMP+yMinMP)/2.-(yMaxMP-yMinMP)*XtitlePositionYshift;

						char axistitleMPdep[200];
						if(iLam==1||iLam==7||iLam==13)  sprintf(axistitleMPdep,"#lambda_{#vartheta}");
						if(iLam==2||iLam==8||iLam==14)  sprintf(axistitleMPdep,"#lambda_{#varphi}");
						if(iLam==3||iLam==9||iLam==15)  sprintf(axistitleMPdep,"#lambda_{#vartheta#varphi}");

						if((nState==4&&rapBin==2)||(nState==5&&rapBin==3)) 
							YaxistitleLatexSize=YaxistitleLatexSize*(1-lowestBottomMargin);
						//TLatex *MPYtitletext = new TLatex(whereTexteInPlotX,whereTexteInPlotY ,axistitleMPdep);
						//MPYtitletext->SetTextSize(YaxistitleLatexSize);
						//if((nState==4&&rapBin==2)||(nState==5&&rapBin==3)) 
						//	MPYtitletext->SetTextSize(YaxistitleLatexSize*(1-lowestBottomMargin));
						//MPYtitletext->SetTextColor(kBlack);
						//MPYtitletext->PaintLatex(whereTexteInPlotX,whereTexteInPlotY, YtitleAngle, YaxistitleLatexSize, axislabel);
						//MPYtitletext->Draw( "same" );

						TLatex *MPYtitletext = new TLatex();
						MPYtitletext->SetTextSize(YaxistitleLatexSize);
						if((nState==4&&rapBin==2)||(nState==5&&rapBin==3)) MPYtitletext->SetTextSize(YaxistitleLatexSize);
						MPYtitletext->SetTextColor(kBlack);
						if(iPanel==1){
							whereTexteInPlotX=0.; whereTexteInPlotY=-0.05;
							MPYtitletext->DrawLatex(whereTexteInPlotX,whereTexteInPlotY,axistitleMPdep);

							sprintf(axistitleMPdep,"#lambda_{#varphi}");
							whereTexteInPlotX=66.; whereTexteInPlotY=-0.05;
							MPYtitletext->DrawLatex(whereTexteInPlotX,whereTexteInPlotY,axistitleMPdep);
						}
						if(iPanel==2){
							sprintf(axistitleMPdep,"#lambda_{#vartheta#varphi}");
							whereTexteInPlotX=64.; whereTexteInPlotY=-0.05;
							MPYtitletext->DrawLatex(whereTexteInPlotX,whereTexteInPlotY,axistitleMPdep);
						}


						char frameMPtex[200];
						if(MPframe==1) sprintf(frameMPtex,"CS");
						if(MPframe==2) sprintf(frameMPtex,"HX");
						if(MPframe==3) sprintf(frameMPtex,"PX");
						char texTexMP[200];
						if(rapBin==1) sprintf(texTexMP,"|#it{y}| < 0.6", nState-3, frameMPtex);
						if(rapBin==2) sprintf(texTexMP,"0.6 < |#it{y}| < 1.2", nState-3, frameMPtex);
						if(rapBin==3) sprintf(texTexMP,"1.2 < |#it{y}| < 1.5", nState-3, frameMPtex);
						TLatex *textMP = new TLatex(xRapText,yMin+(yMax-yMin)*yRapText,texTexMP);
						textMP->SetTextSize(textSizeRap);
						if((nState==4&&rapBin==2)||(nState==5&&rapBin==3)) textMP->SetTextSize(textSizeRap*(1-lowestBottomMargin));
						textMP->Draw( "same" );

						char abcdef[200];
						if(rapBin==1&&iPanel==1) sprintf(abcdef,"a)");
						if(rapBin==1&&iPanel==2) sprintf(abcdef,"b)");
						if(rapBin==1&&iPanel==3) sprintf(abcdef,"c)");
						if(rapBin==2&&iPanel==1) sprintf(abcdef,"d)");
						if(rapBin==2&&iPanel==2) sprintf(abcdef,"e)");
						if(rapBin==2&&iPanel==3) sprintf(abcdef,"f)");
						if(rapBin==3&&iPanel==1) sprintf(abcdef,"g)");
						if(rapBin==3&&iPanel==2) sprintf(abcdef,"h)");
						if(rapBin==3&&iPanel==3) sprintf(abcdef,"i)");
						cout<<abcdef<<endl;
						TLatex *tex_abcdef = new TLatex(xabcdefText,yMin+(yMax-yMin)*yRapText*0.92,abcdef);
						tex_abcdef->SetTextSize(textSizeRap);
						if((nState==4&&rapBin==2)||(nState==5&&rapBin==3)) tex_abcdef->SetTextSize(textSizeRap*(1-lowestBottomMargin));
						//tex_abcdef->Draw( "same" );

						if(PlotFinalData&&DrawLatexStuff){

							TLine* extreme0MP = new TLine( PlotpTMin-DeltaXminOVERALL, 0, PlotpTMax ,0);
							extreme0MP->SetLineWidth( 1 );
							extreme0MP->SetLineStyle( 2 );
							extreme0MP->SetLineColor( kBlack );
							extreme0MP->Draw( "same" );

							//TLine* extreme1MP = new TLine( PlotpTMin-DeltaXminOVERALL, 1, PlotpTMax, 1);
							//extreme1MP->SetLineWidth( 1 );
							//extreme1MP->SetLineStyle( 2 );
							//extreme1MP->SetLineColor( kBlack );
							//if(iLam==1||iLam==7||iLam==13) extreme1MP->Draw( "same" );

							//TLine* extreme2MP = new TLine( PlotpTMin-DeltaXminOVERALL, -1, PlotpTMax ,-1);
							//extreme2MP->SetLineWidth( 1 );
							//extreme2MP->SetLineStyle( 2 );
							//extreme2MP->SetLineColor( kBlack );
							//if(iLam==1||iLam==7||iLam==13) extreme2MP->Draw( "same" );
							//if(iLam==6||iLam==12||iLam==18) extreme2MP->Draw( "same" );
						}


						if(rapBin==1&&iPanel==1){
							//cout<<"DRAW CMS preliminary Latex"<<endl;
							//char text[200];
							//sprintf(text,"CMS preliminary   pp   #sqrt{s} = 7 TeV   L = 4.9 fb^{-1}");
							//TLatex *CentralsText1MP = new TLatex(MPlatexX,MPlatexYmax,text);
							//CentralsText1MP->SetTextSize(CentralsFontSizeMP);
							//CentralsText1MP->Draw( "same" );

							//sprintf(text,"L = 4.9 fb^{-1}");
							//TLatex *CentralsText2MP = new TLatex(MPlatexX,MPlatexYmax-2*MPlatexDeltaYmax,text);
							//CentralsText2MP->SetTextSize(CentralsFontSizeMP);
							//CentralsText2MP->Draw( "same" );
							//sprintf(text,"pp    #sqrt{s} = 7 TeV");
							//TLatex *CentralsText3MP = new TLatex(MPlatexX,MPlatexYmax-MPlatexDeltaYmax,text);
							//CentralsText3MP->SetTextSize(CentralsFontSizeMP);
							//CentralsText3MP->Draw( "same" );
						}

						if(rapBin==2&&iPanel==1&&(nState==1||nState==4)){
							MPframedepLegend->Draw("same");
						}
						if(nState==2||nState==3){
							if(rapBin==1&&iPanel==3&&MPframe==1){
								MPframedepLegend->Draw("same");
							}
							if(rapBin==1&&iPanel==2&&MPframe!=1){
								MPframedepLegend->Draw("same");
							}
						}
						if(rapBin==1&&iPanel==2){

							char frameMPtex[200];
							if(MPframe==1) sprintf(frameMPtex,"CS frame");
							if(MPframe==2) sprintf(frameMPtex,"HX frame");
							if(MPframe==3) sprintf(frameMPtex,"PX frame");
							char textStateFrame[200];
							sprintf(textStateFrame,"#psi(%dS), %s", nState-3, frameMPtex);
							TLatex *TexStateFrame = new TLatex(MPlatexX+10,MPlatexYmax,textStateFrame);
							TexStateFrame->SetTextSize(CentralsFontSizeMP);
							TexStateFrame->Draw( "same" );

						}

						if(MPframe==1) MPcanvasCS_Psi->cd();
						if(MPframe==2) MPcanvasHX_Psi->cd();
						if(MPframe==3) MPcanvasPX_Psi->cd();

						whereTexteInPlotX=0.488;
						whereTexteInPlotY=startValCoordY-deltaCoordY-1.375*labelOffsetX;

						TLatex *MPXlabeltext = new TLatex(whereTexteInPlotX,whereTexteInPlotY ,"1");
						MPXlabeltext->SetTextSize(XaxislabelLatexSize);
						MPXlabeltext->SetTextColor(kBlack);
						//if(iPanel==nPanels) MPXlabeltext->Draw( "same" );

						if(iLam==3&&((nState<=4&&rapBin==2)||(nState==5&&rapBin==3))){  
							sprintf(filename,"%s/FinalResultsCS_Psi%dS_OLD.pdf",FigDir,nState-3);
							if(PlotFinalData) MPcanvasCS_Psi->SaveAs(filename);
							MPcanvasCS_Psi->Close();
						}
						if(iLam==9&&((nState<=4&&rapBin==2)||(nState==5&&rapBin==3))){  
							sprintf(filename,"%s/FinalResultsHX_Psi%dS_OLD.pdf",FigDir,nState-3);
							if(PlotFinalData) MPcanvasHX_Psi->SaveAs(filename);
							MPcanvasHX_Psi->Close();
						}
						if(iLam==15&&((nState<=4&&rapBin==2)||(nState==5&&rapBin==3))){  
							sprintf(filename,"%s/FinalResultsPX_Psi%dS_OLD.pdf",FigDir,nState-3);
							if(PlotFinalData) MPcanvasPX_Psi->SaveAs(filename);
							MPcanvasPX_Psi->Close();
						}

					}//end Frame dependent plots
				}

				//===============================================================================================





				///////////////////////////////////////////////////////
				///////////////////////////////////////////////////////
				////////////// NEW MP PLOTS ///////////////////////////
				///////////////////////////////////////////////////////
				///////////////////////////////////////////////////////

				//// Pad Definitions
				//float Top_margin_MPnew   = 0.;//0
				//float Left_margin_MPnew  = 0.15;//0.025
				//float Right_margin_MPnew = 0.15;//0.005
				//const int nPanels_MPnew=3;
				//double lowestBottomMargin_MPnew=0.3;
				//double PadCoordYMax_MPnew=0.95;
				//double deltaCoordY_MPnew=PadCoordYMax_MPnew/(double(nPanels_MPnew-1)+1./(1-lowestBottomMargin_MPnew));
				//double startValCoordY_MPnew=deltaCoordY_MPnew/(1-lowestBottomMargin_MPnew);
				////double PadCoordY[nPanels+1]={0.,0.3,0.5,0.7,0.9};
				//double PadCoordY_MPnew[nPanels_MPnew+1]={0.,startValCoordY_MPnew,
				//startValCoordY_MPnew+deltaCoordY_MPnew,PadCoordYMax_MPnew};
				//double PadCoordX_MPnew[3]={0.1,0.5,0.9};

				// Canvas Definitions
				int MPcanvasXpixel_MPnew=MPcanvasXpixelInitial*1.5;
				int MPcanvasYpixel_MPnew=MPcanvasYpixelInitial;

				const int nXPanels_newMP=3;
				double x_1=1./(1.+(1-Left_margin)/(1-Right_margin)+1-Left_margin);//(1.-x_tilde)/2.;

				double x_tilde=x_1*(1-Left_margin);//1./(1.+2./(1-Left_margin));
				double PadCoordX_newMP[nXPanels_newMP+1]={0.,x_1,x_1+x_tilde,1.};
				xRapText_MPnew=44;

				cout<<"begin NEW Frame dependent plots"<<endl;
				//begin Frame dependent plots
				if(NEW_MPplots){
					if(iLam!=6&&iLam!=12&&iLam!=18){

						if(iLam>0&&iLam<7) MPframe=1;
						if(iLam>6&&iLam<13) MPframe=2;
						if(iLam>12&&iLam<19) MPframe=3;

						if(iLam==1||iLam==7||iLam==13) iPanel=1;
						if(iLam==2||iLam==8||iLam==14) iPanel=2;
						if(iLam==3||iLam==9||iLam==15) iPanel=3;

						cout<<"MultiPanel canvas"<<endl;

						if(iLam==1&&rapBin==1){
							MPcanvasCS_rap1 = new TCanvas("MPcanvasCS_rap1", "MPcanvasCS_rap1",MPcanvasXpixel_MPnew,MPcanvasYpixel_MPnew);
							MPcanvasCS_rap1->SetFillColor(kWhite);
							MPcanvasCS_rap1->GetFrame()->SetFillColor(kWhite);
							MPcanvasCS_rap1->GetFrame()->SetBorderSize(0);
						}
						if(iLam==7&&rapBin==1){
							MPcanvasHX_rap1 = new TCanvas("MPcanvasHX_rap1", "MPcanvasHX_rap1",MPcanvasXpixel_MPnew,MPcanvasYpixel_MPnew);
							MPcanvasHX_rap1->SetFillColor(kWhite);
							MPcanvasHX_rap1->GetFrame()->SetFillColor(kWhite);
							MPcanvasHX_rap1->GetFrame()->SetBorderSize(0);
						}
						if(iLam==13&&rapBin==1){
							MPcanvasPX_rap1 = new TCanvas("MPcanvasPX_rap1", "MPcanvasPX_rap1",MPcanvasXpixel_MPnew,MPcanvasYpixel_MPnew);
							MPcanvasPX_rap1->SetFillColor(kWhite);
							MPcanvasPX_rap1->GetFrame()->SetFillColor(kWhite);
							MPcanvasPX_rap1->GetFrame()->SetBorderSize(0);
						}
						if(iLam==1&&rapBin==2){
							MPcanvasCS_rap2 = new TCanvas("MPcanvasCS_rap2", "MPcanvasCS_rap2",MPcanvasXpixel_MPnew,MPcanvasYpixel_MPnew);
							MPcanvasCS_rap2->SetFillColor(kWhite);
							MPcanvasCS_rap2->GetFrame()->SetFillColor(kWhite);
							MPcanvasCS_rap2->GetFrame()->SetBorderSize(0);
						}
						if(iLam==7&&rapBin==2){
							MPcanvasHX_rap2 = new TCanvas("MPcanvasHX_rap2", "MPcanvasHX_rap2",MPcanvasXpixel_MPnew,MPcanvasYpixel_MPnew);
							MPcanvasHX_rap2->SetFillColor(kWhite);
							MPcanvasHX_rap2->GetFrame()->SetFillColor(kWhite);
							MPcanvasHX_rap2->GetFrame()->SetBorderSize(0);
						}
						if(iLam==13&&rapBin==2){
							MPcanvasPX_rap2 = new TCanvas("MPcanvasPX_rap2", "MPcanvasPX_rap2",MPcanvasXpixel_MPnew,MPcanvasYpixel_MPnew);
							MPcanvasPX_rap2->SetFillColor(kWhite);
							MPcanvasPX_rap2->GetFrame()->SetFillColor(kWhite);
							MPcanvasPX_rap2->GetFrame()->SetBorderSize(0);
						}
						if(iLam==1&&rapBin==3){
							MPcanvasCS_rap3 = new TCanvas("MPcanvasCS_rap3", "MPcanvasCS_rap3",MPcanvasXpixel_MPnew,MPcanvasYpixel_MPnew);
							MPcanvasCS_rap3->SetFillColor(kWhite);
							MPcanvasCS_rap3->GetFrame()->SetFillColor(kWhite);
							MPcanvasCS_rap3->GetFrame()->SetBorderSize(0);
						}
						if(iLam==7&&rapBin==3){
							MPcanvasHX_rap3 = new TCanvas("MPcanvasHX_rap3", "MPcanvasHX_rap3",MPcanvasXpixel_MPnew,MPcanvasYpixel_MPnew);
							MPcanvasHX_rap3->SetFillColor(kWhite);
							MPcanvasHX_rap3->GetFrame()->SetFillColor(kWhite);
							MPcanvasHX_rap3->GetFrame()->SetBorderSize(0);
						}
						if(iLam==13&&rapBin==3){
							MPcanvasPX_rap3 = new TCanvas("MPcanvasPX_rap3", "MPcanvasPX_rap3",MPcanvasXpixel_MPnew,MPcanvasYpixel_MPnew);
							MPcanvasPX_rap3->SetFillColor(kWhite);
							MPcanvasPX_rap3->GetFrame()->SetFillColor(kWhite);
							MPcanvasPX_rap3->GetFrame()->SetBorderSize(0);
						}


						int totalState = 3;
						if(nState>3) totalState = 2;
						for(int iStateMP=1;iStateMP<totalState+1;iStateMP++){

							if(MPframe==1&&rapBin==1) MPcanvasCS_rap1->cd();
							if(MPframe==2&&rapBin==1) MPcanvasHX_rap1->cd();
							if(MPframe==3&&rapBin==1) MPcanvasPX_rap1->cd();

							if(MPframe==1&&rapBin==2) MPcanvasCS_rap2->cd();
							if(MPframe==2&&rapBin==2) MPcanvasHX_rap2->cd();
							if(MPframe==3&&rapBin==2) MPcanvasPX_rap2->cd();

							if(MPframe==1&&rapBin==3) MPcanvasCS_rap3->cd();
							if(MPframe==2&&rapBin==3) MPcanvasHX_rap3->cd();
							if(MPframe==3&&rapBin==3) MPcanvasPX_rap3->cd();


							cout<<"MultiPanel pad"<<endl;
							cout<<"iLam "<<iLam<<endl;
							cout<<"MPframe "<<MPframe<<endl;
							cout<<"iStateMP "<<iStateMP<<endl;
							cout<<"iPanel "<<iPanel<<endl;
							cout<<"PadYlow "<<PadCoordY[nPanels-iPanel]<<endl;
							cout<<"PadYlow "<<PadCoordY[nPanels-iPanel+1]<<endl;

							TPad *MPpad;
							if(iStateMP==1) MPpad = new TPad("MPpad","MPpad",PadCoordX_newMP[iStateMP-1],PadCoordY[nPanels-iPanel],PadCoordX_newMP[iStateMP],PadCoordY[nPanels-iPanel+1]);
							if(iStateMP==2) MPpad = new TPad("MPpad","MPpad",PadCoordX_newMP[iStateMP-1],PadCoordY[nPanels-iPanel],PadCoordX_newMP[iStateMP],PadCoordY[nPanels-iPanel+1]);
							if(iStateMP==3) MPpad = new TPad("MPpad","MPpad",PadCoordX_newMP[iStateMP-1],PadCoordY[nPanels-iPanel],PadCoordX_newMP[iStateMP],PadCoordY[nPanels-iPanel+1]);
							MPpad->Draw();
							MPpad->cd();
							MPpad->SetFillColor(kWhite);
							MPpad->SetFrameFillColor(kWhite);
							MPpad->SetBorderSize(0);
							MPpad->SetLeftMargin(0.);
							if(iStateMP==1) MPpad->SetLeftMargin(Left_margin);
							MPpad->SetRightMargin(0.);
							if(iStateMP==3) MPpad->SetRightMargin(Right_margin);
							MPpad->SetTopMargin(Top_margin+0.0025);
							MPpad->SetBottomMargin(0.0);
							if(iPanel==nPanels) MPpad->SetBottomMargin(lowestBottomMargin);

							double deltaXaxisMin=0.;

							cout<<"MultiPanel hist"<<endl;
							TH1F *MPhist = new TH1F;
							if(MPframe==1&&rapBin==1)MPhist = MPcanvasCS_rap1->DrawFrame(onia::pTRange[rapBin][ptBinMin-1]-DeltaXminOVERALL-deltaXaxisMin,yMinMP,onia::pTRange[rapBin][ptBinMax],yMaxMP);
							if(MPframe==2&&rapBin==1)MPhist = MPcanvasHX_rap1->DrawFrame(onia::pTRange[rapBin][ptBinMin-1]-DeltaXminOVERALL-deltaXaxisMin,yMinMP,onia::pTRange[rapBin][ptBinMax],yMaxMP);
							if(MPframe==3&&rapBin==1)MPhist = MPcanvasPX_rap1->DrawFrame(onia::pTRange[rapBin][ptBinMin-1]-DeltaXminOVERALL-deltaXaxisMin,yMinMP,onia::pTRange[rapBin][ptBinMax],yMaxMP);
							if(MPframe==1&&rapBin==2)MPhist = MPcanvasCS_rap2->DrawFrame(onia::pTRange[rapBin][ptBinMin-1]-DeltaXminOVERALL-deltaXaxisMin,yMinMP,onia::pTRange[rapBin][ptBinMax],yMaxMP);
							if(MPframe==2&&rapBin==2)MPhist = MPcanvasHX_rap2->DrawFrame(onia::pTRange[rapBin][ptBinMin-1]-DeltaXminOVERALL-deltaXaxisMin,yMinMP,onia::pTRange[rapBin][ptBinMax],yMaxMP);
							if(MPframe==3&&rapBin==2)MPhist = MPcanvasPX_rap2->DrawFrame(onia::pTRange[rapBin][ptBinMin-1]-DeltaXminOVERALL-deltaXaxisMin,yMinMP,onia::pTRange[rapBin][ptBinMax],yMaxMP);

							MPhist->SetXTitle("#it{p}_{T} [GeV]");
							MPhist->GetXaxis()->SetTitleOffset(-1.35);

							MPhist->SetYTitle(axislabel);
							MPhist->GetYaxis()->SetTitleOffset(titleoffset);
							MPhist->GetYaxis()->SetTitleSize(0.);
							if(iPanel==nPanels) MPhist->GetYaxis()->SetTitleOffset(titleoffset*1.35);
							if(iPanel==nPanels) MPhist->GetYaxis()->SetTitleSize(0.*(1-lowestBottomMargin));

							MPhist->GetYaxis()->SetLabelSize(LabelSize*1.25);
							MPhist->GetXaxis()->SetLabelSize(0.);
							MPhist->GetYaxis()->SetLabelOffset(-0.015);
							MPhist->GetXaxis()->SetLabelOffset(-0.06);

							if(iPanel==nPanels) MPhist->GetYaxis()->SetLabelSize(LabelSize*1.25*(1-lowestBottomMargin));
							MPhist->GetXaxis()->SetTitleSize(TitleSize*0.85);
							MPhist->GetXaxis()->SetAxisColor(kWhite);
							MPhist->GetYaxis()->SetAxisColor(kWhite);
							MPhist->GetXaxis()->SetTicks("-");
							MPhist->GetYaxis()->SetTicks("+");


							double SpecialShift=0.01;
							TLegend* MPframedepLegend;
							MPframedepLegend=new TLegend(errorLegendX1,errorLegendY1,errorLegendX2,errorLegendY2);
							if(MPframe==1) MPframedepLegend=new TLegend(errorLegendX1,
									errorLegendY2-(errorLegendY2-errorLegendY1)*(1-lowestBottomMargin)+SpecialShift,
									errorLegendX2,errorLegendY2+SpecialShift);
							MPframedepLegend->SetFillColor(0);
							//MPframedepLegend->SetTextFont(72);
							MPframedepLegend->SetTextSize(errorLegendFontSize);
							if(MPframe==1) MPframedepLegend->SetTextSize(errorLegendFontSize*(1-lowestBottomMargin));
							MPframedepLegend->SetBorderSize(0);

							char MPframedepLegendEntry[200];

							TGraphAsymmErrors* graphMP_3sig_MPnew;
							TGraphAsymmErrors* graphMP_2sig_MPnew;
							TGraphAsymmErrors* graphMP_1sig_MPnew;
							TGraphAsymmErrors* graphMP_MPnew;

							if(iStateMP==1) graphMP_MPnew = (TGraphAsymmErrors*) infileMP1->Get(GraphName);
							if(iStateMP==2) graphMP_MPnew = (TGraphAsymmErrors*) infileMP2->Get(GraphName);
							if(iStateMP==3) graphMP_MPnew = (TGraphAsymmErrors*) infileMP3->Get(GraphName);

							if(iStateMP==1) graphMP_1sig_MPnew = (TGraphAsymmErrors*) infileMP1_1sig->Get(GraphName);
							if(iStateMP==2) graphMP_1sig_MPnew = (TGraphAsymmErrors*) infileMP2_1sig->Get(GraphName);
							if(iStateMP==3) graphMP_1sig_MPnew = (TGraphAsymmErrors*) infileMP3_1sig->Get(GraphName);
							if(iStateMP==1) graphMP_2sig_MPnew = (TGraphAsymmErrors*) infileMP1_2sig->Get(GraphName);
							if(iStateMP==2) graphMP_2sig_MPnew = (TGraphAsymmErrors*) infileMP2_2sig->Get(GraphName);
							if(iStateMP==3) graphMP_2sig_MPnew = (TGraphAsymmErrors*) infileMP3_2sig->Get(GraphName);
							if(iStateMP==1) graphMP_3sig_MPnew = (TGraphAsymmErrors*) infileMP1_3sig->Get(GraphName);
							if(iStateMP==2) graphMP_3sig_MPnew = (TGraphAsymmErrors*) infileMP2_3sig->Get(GraphName);
							if(iStateMP==3) graphMP_3sig_MPnew = (TGraphAsymmErrors*) infileMP3_3sig->Get(GraphName);

							int ptMP;
							int nBinsMP;
							double ptCentre_MP[nBinspT];
							double lmean_MP[nBinspT];
							double lmean_errlow_MP[nBinspT];
							double lmean_errhigh_MP[nBinspT];
							double ptCentre_errlow_MP[nBinspT];
							double ptCentre_errhigh_MP[nBinspT];

							nBinsMP=graphMP_1sig_MPnew->GetN();
							ptMP=0;
							for(int ptBinMP=1;ptBinMP<nBinsMP+1;ptBinMP++){
								graphMP_1sig_MPnew->GetPoint(ptBinMP-1,ptCentre_MP[ptMP],lmean_MP[ptMP]);
								lmean_errhigh_MP[ptMP]=graphMP_1sig_MPnew->GetErrorYhigh(ptBinMP-1);
								lmean_errlow_MP[ptMP]=graphMP_1sig_MPnew->GetErrorYlow(ptBinMP-1);
								ptCentre_errhigh_MP[ptMP]=graphMP_1sig_MPnew->GetErrorXhigh(ptBinMP-1);
								ptCentre_errlow_MP[ptMP]=graphMP_1sig_MPnew->GetErrorXlow(ptBinMP-1);
								/// Alter TGraph
								ptCentre_errhigh_MP[ptMP]=ColordBandWidth;
								ptCentre_errlow_MP[ptMP]=ColordBandWidth;
								ptMP++;
							}
							graphMP_1sig_MPnew = new TGraphAsymmErrors(nBinsMP,ptCentre_MP,lmean_MP,ptCentre_errlow_MP,ptCentre_errhigh_MP,lmean_errlow_MP,lmean_errhigh_MP);

							nBinsMP=graphMP_2sig_MPnew->GetN();
							ptMP=0;
							for(int ptBinMP=1;ptBinMP<nBinsMP+1;ptBinMP++){
								graphMP_2sig_MPnew->GetPoint(ptBinMP-1,ptCentre_MP[ptMP],lmean_MP[ptMP]);
								lmean_errhigh_MP[ptMP]=graphMP_2sig_MPnew->GetErrorYhigh(ptBinMP-1);
								lmean_errlow_MP[ptMP]=graphMP_2sig_MPnew->GetErrorYlow(ptBinMP-1);
								ptCentre_errhigh_MP[ptMP]=graphMP_2sig_MPnew->GetErrorXhigh(ptBinMP-1);
								ptCentre_errlow_MP[ptMP]=graphMP_2sig_MPnew->GetErrorXlow(ptBinMP-1);
								/// Alter TGraph
								ptCentre_errhigh_MP[ptMP]=ColordBandWidth;
								ptCentre_errlow_MP[ptMP]=ColordBandWidth;
								ptMP++;
							}
							graphMP_2sig_MPnew = new TGraphAsymmErrors(nBinsMP,ptCentre_MP,lmean_MP,ptCentre_errlow_MP,ptCentre_errhigh_MP,lmean_errlow_MP,lmean_errhigh_MP);

							nBinsMP=graphMP_3sig_MPnew->GetN();
							ptMP=0;
							for(int ptBinMP=1;ptBinMP<nBinsMP+1;ptBinMP++){
								graphMP_3sig_MPnew->GetPoint(ptBinMP-1,ptCentre_MP[ptMP],lmean_MP[ptMP]);
								lmean_errhigh_MP[ptMP]=graphMP_3sig_MPnew->GetErrorYhigh(ptBinMP-1);
								lmean_errlow_MP[ptMP]=graphMP_3sig_MPnew->GetErrorYlow(ptBinMP-1);
								ptCentre_errhigh_MP[ptMP]=graphMP_3sig_MPnew->GetErrorXhigh(ptBinMP-1);
								ptCentre_errlow_MP[ptMP]=graphMP_3sig_MPnew->GetErrorXlow(ptBinMP-1);
								/// Alter TGraph
								ptCentre_errhigh_MP[ptMP]=ColordBandWidth;
								ptCentre_errlow_MP[ptMP]=ColordBandWidth;
								ptMP++;
							}
							graphMP_3sig_MPnew = new TGraphAsymmErrors(nBinsMP,ptCentre_MP,lmean_MP,ptCentre_errlow_MP,ptCentre_errhigh_MP,lmean_errlow_MP,lmean_errhigh_MP);

							nBinsMP=graphMP_MPnew->GetN();
							ptMP=0;
							for(int ptBinMP=1;ptBinMP<nBinsMP+1;ptBinMP++){
								graphMP_MPnew->GetPoint(ptBinMP-1,ptCentre_MP[ptMP],lmean_MP[ptMP]);
								lmean_errhigh_MP[ptMP]=graphMP_MPnew->GetErrorYhigh(ptBinMP-1);
								lmean_errlow_MP[ptMP]=graphMP_MPnew->GetErrorYlow(ptBinMP-1);
								ptCentre_errhigh_MP[ptMP]=graphMP_MPnew->GetErrorXhigh(ptBinMP-1);
								ptCentre_errlow_MP[ptMP]=graphMP_MPnew->GetErrorXlow(ptBinMP-1);
								/// Alter TGraph
								ptCentre_errhigh_MP[ptMP]=0.;
								ptCentre_errlow_MP[ptMP]=0.;
								ptMP++;
							}

							bool RemoveHorizontalErrorBar=true;
							if(RemoveHorizontalErrorBar) graphMP_MPnew = new TGraphAsymmErrors(nBinsMP,ptCentre_MP,lmean_MP,ptCentre_errlow_MP,ptCentre_errhigh_MP,lmean_errlow_MP,lmean_errhigh_MP);

							graphMP_MPnew->SetMarkerSize(MarkerSizeMP[0]);
							graphMP_MPnew->SetMarkerStyle(MarkerStyleMP[0]);
							graphMP_MPnew->SetMarkerColor(MarkerColorMP[0]);

							graphMP_1sig_MPnew->SetFillColor(OneSigColor);
							graphMP_1sig_MPnew->SetFillStyle(1001);
							graphMP_2sig_MPnew->SetFillColor(TwoSigColor);
							graphMP_2sig_MPnew->SetFillStyle(1001);
							graphMP_3sig_MPnew->SetFillColor(ThreeSigColor);
							graphMP_3sig_MPnew->SetFillStyle(1001);


							graphMP_3sig_MPnew->Draw("2");
							graphMP_2sig_MPnew->Draw("2");
							graphMP_1sig_MPnew->Draw("2");
							if(PlotAlteredPPDResults) {
								graphMP_MPnew->Draw("[]");
								graphMP_MPnew->Draw(drawGraphStyle);
							}


							graphMP_1sig_MPnew->SetLineColor(OneSigColor);
							graphMP_2sig_MPnew->SetLineColor(TwoSigColor);
							graphMP_3sig_MPnew->SetLineColor(ThreeSigColor);

							if(PlotAlteredPPDResults) sprintf(MPframedepLegendEntry,"Stat. uncert., 68.3 %% CL");
							if(!PlotAlteredPPDResults) sprintf(MPframedepLegendEntry,"Tot. sys. uncert.");
							MPframedepLegend->AddEntry(graphMP_MPnew,MPframedepLegendEntry,"p");
							if(PlotAlteredPPDResults) sprintf(MPframedepLegendEntry,"Tot. uncert., 68.3 %% CL");
							if(!PlotAlteredPPDResults) sprintf(MPframedepLegendEntry,"68.3 %% CL");
							MPframedepLegend->AddEntry(graphMP_1sig_MPnew,MPframedepLegendEntry,"f");
							if(PlotAlteredPPDResults) sprintf(MPframedepLegendEntry,"Tot. uncert., 95.5 %% CL");
							if(!PlotAlteredPPDResults) sprintf(MPframedepLegendEntry,"95.5 %% CL");
							MPframedepLegend->AddEntry(graphMP_2sig_MPnew,MPframedepLegendEntry,"f");
							if(PlotAlteredPPDResults) sprintf(MPframedepLegendEntry,"Tot. uncert., 99.7 %% CL");
							if(!PlotAlteredPPDResults) sprintf(MPframedepLegendEntry,"99.7 %% CL");
							MPframedepLegend->AddEntry(graphMP_3sig_MPnew,MPframedepLegendEntry,"f");

							TGaxis *axisMPY1 = new TGaxis(onia::pTRange[rapBin][ptBinMin-1]-DeltaXminOVERALL,yMinMP,onia::pTRange[rapBin][ptBinMin-1]-DeltaXminOVERALL,yMaxMP,yMinMP,yMaxMP,AxisDivisions,"-US");
							axisMPY1->SetTickSize(ticksize);
							if(iPanel==nPanels) axisMPY1->SetTickSize(ticksize/(1-lowestBottomMargin));
							axisMPY1->Draw("same");

							TGaxis *axisMPY2 = new TGaxis(onia::pTRange[rapBin][ptBinMax],yMinMP,onia::pTRange[rapBin][ptBinMax],yMaxMP,yMinMP,yMaxMP,AxisDivisions,"+US");
							axisMPY2->SetTickSize(ticksize);
							if(iPanel==nPanels) axisMPY2->SetTickSize(ticksize/(1-lowestBottomMargin));
							axisMPY2->Draw("same");


							double deltaTrickAxisMax_MPnew;
							deltaTrickAxisMax_MPnew=-0.001;
							//if(iStateMP==3) deltaTrickAxisMax_MPnew=+0.001;

							TGaxis *axisMPX1 = new TGaxis(onia::pTRange[rapBin][ptBinMin-1]-DeltaXminOVERALL,yMinMP,onia::pTRange[rapBin][ptBinMax],yMinMP,onia::pTRange[rapBin][ptBinMin-1]-DeltaXminOVERALL+deltaTrickAxisMin,onia::pTRange[rapBin][ptBinMax]+deltaTrickAxisMax_MPnew,AxisDivisions,"+S");
							axisMPX1->SetTickSize(ticksize*2);
							if(iPanel==nPanels) axisMPX1->SetLabelSize(LabelSize);
							if(iPanel<nPanels) axisMPX1->SetLabelSize(0);
							axisMPX1->SetLabelOffset(labelOffsetX);
							if(iPanel==nPanels) axisMPX1->SetTickSize(ticksize/(1-lowestBottomMargin));
							axisMPX1->Draw("same");

							TGaxis *axisMPX2 = new TGaxis(onia::pTRange[rapBin][ptBinMin-1]-DeltaXminOVERALL,yMaxMP,onia::pTRange[rapBin][ptBinMax],yMaxMP,onia::pTRange[rapBin][ptBinMin-1]-DeltaXminOVERALL+deltaTrickAxisMin,onia::pTRange[rapBin][ptBinMax]+deltaTrickAxisMax_MPnew,AxisDivisions,"-US");
							axisMPX2->SetTickSize(ticksize*2);
							if(iPanel==nPanels) axisMPX2->SetTickSize(ticksize/(1-lowestBottomMargin));
							axisMPX2->Draw("same");

							whereTexteInPlotX=XtitlePosition;
							whereTexteInPlotY=(yMaxMP+yMinMP)/2.-(yMaxMP-yMinMP)*XtitlePositionYshift;

							char axistitleMPdep[200];
							if(iLam==1||iLam==7||iLam==13)  sprintf(axistitleMPdep,"#lambda_{#vartheta}");
							if(iLam==2||iLam==8||iLam==14)  sprintf(axistitleMPdep,"#lambda_{#varphi}");
							if(iLam==3||iLam==9||iLam==15)  sprintf(axistitleMPdep,"#lambda_{#vartheta#varphi}");

							if(iPanel==nPanels) YaxistitleLatexSize=YaxistitleLatexSize*(1-lowestBottomMargin);
							TLatex *MPYtitletext = new TLatex(whereTexteInPlotX,whereTexteInPlotY ,axistitleMPdep);
							MPYtitletext->SetTextSize(YaxistitleLatexSize);
							if(iPanel==nPanels) MPYtitletext->SetTextSize(YaxistitleLatexSize*(1-lowestBottomMargin));
							MPYtitletext->SetTextColor(kBlack);
							MPYtitletext->PaintLatex(whereTexteInPlotX-DeltaXminOVERALL,whereTexteInPlotY, YtitleAngle, YaxistitleLatexSize, axislabel);
							MPYtitletext->Draw( "same" );



							double increaseSize=1.45;
							double SpecialShiftUpsilonLabel=0.9;
							double SpecialShiftUpsilonLabelx=-3;
							SpecialShiftUpsilonLabel=14;//1.25;//FRchange

							char frameMPtex[200];
							if(MPframe==1) sprintf(frameMPtex,"CS");
							if(MPframe==2) sprintf(frameMPtex,"HX");
							if(MPframe==3) sprintf(frameMPtex,"PX");
							char texTexMP[200];
							sprintf(texTexMP,"#psi(%dS)", iStateMP);
							TLatex *textMP = new TLatex(xRapText_MPnew+SpecialShiftUpsilonLabelx,yMin+(yMax-yMin)*yRapText*SpecialShiftUpsilonLabel,texTexMP);
							textMP->SetTextSize(textSizeRap*increaseSize);
							if(iPanel==nPanels) textMP->SetTextSize(textSizeRap*(1-lowestBottomMargin)*increaseSize);
							textMP->Draw( "same" );

							char abcdef[200];
							if(iStateMP==1&&iPanel==1) sprintf(abcdef,"a)");
							if(iStateMP==1&&iPanel==2) sprintf(abcdef,"b)");
							if(iStateMP==1&&iPanel==3) sprintf(abcdef,"c)");
							if(iStateMP==2&&iPanel==1) sprintf(abcdef,"d)");
							if(iStateMP==2&&iPanel==2) sprintf(abcdef,"e)");
							if(iStateMP==2&&iPanel==3) sprintf(abcdef,"f)");
							if(iStateMP==3&&iPanel==1) sprintf(abcdef,"g)");
							if(iStateMP==3&&iPanel==2) sprintf(abcdef,"h)");
							if(iStateMP==3&&iPanel==3) sprintf(abcdef,"i)");

							cout<<abcdef<<endl;
							TLatex *tex_abcdef = new TLatex(xabcdefText,yMin+(yMax-yMin)*yRapText*0.92,abcdef);
							tex_abcdef->SetTextSize(textSizeRap);
							if(iPanel==nPanels) tex_abcdef->SetTextSize(textSizeRap*(1-lowestBottomMargin));
							//tex_abcdef->Draw( "same" );

							if(PlotFinalData&&DrawLatexStuff){

								TLine* extreme0MP = new TLine( PlotpTMin-DeltaXminOVERALL, 0, PlotpTMax ,0);
								extreme0MP->SetLineWidth( 1 );
								extreme0MP->SetLineStyle( 2 );
								extreme0MP->SetLineColor( kBlack );
								extreme0MP->Draw( "same" );

								TLine* extreme1MP = new TLine( PlotpTMin-DeltaXminOVERALL, 1, PlotpTMax , 1);
								extreme1MP->SetLineWidth( 1 );
								extreme1MP->SetLineStyle( 2 );
								extreme1MP->SetLineColor( kBlack );
								//if(iLam==1||iLam==7||iLam==13) extreme1MP->Draw( "same" );

								TLine* extreme2MP = new TLine( PlotpTMin-DeltaXminOVERALL, -1, PlotpTMax ,-1);
								extreme2MP->SetLineWidth( 1 );
								extreme2MP->SetLineStyle( 2 );
								extreme2MP->SetLineColor( kBlack );
								//if(iLam==1||iLam==7||iLam==13) extreme2MP->Draw( "same" );
								//if(iLam==6||iLam==12||iLam==18) extreme2MP->Draw( "same" );

							}

							MPlatexYmax=(yMax-yMin)*0.08+yMin;//FRchange

							if(iStateMP==1&&iPanel==1){
								cout<<"DRAW CMS preliminary Latex"<<endl;
								char text[200];
								sprintf(text,"CMS     pp      #sqrt{s} = 7 TeV     L = 4.9 fb^{-1}");
								TLatex *CentralsText1MP = new TLatex(MPlatexX-DeltaXminOVERALL,MPlatexYmax,text);
								CentralsText1MP->SetTextSize(CentralsFontSizeMP);
								CentralsText1MP->Draw( "same" );
								sprintf(text,"preliminary");
								CentralsText1MP = new TLatex(MPlatexX-DeltaXminOVERALL,(yMax-yMin)*0.2+yMin,text);
								CentralsText1MP->SetTextSize(CentralsFontSizeMP);
								if(DrawPreliminary) CentralsText1MP->Draw( "same" );
								//sprintf(text,"L = 4.9 fb^{-1}");
								//TLatex *CentralsText2MP = new TLatex(MPlatexX,MPlatexYmax-2*MPlatexDeltaYmax,text);
								//CentralsText2MP->SetTextSize(CentralsFontSizeMP);
								//CentralsText2MP->Draw( "same" );
								//sprintf(text,"pp    #sqrt{s} = 7 TeV");
								//TLatex *CentralsText3MP = new TLatex(MPlatexX,MPlatexYmax-MPlatexDeltaYmax,text);
								//CentralsText3MP->SetTextSize(CentralsFontSizeMP);
								//CentralsText3MP->Draw( "same" );
							}

							if(iStateMP==1&&iPanel==2&&MPframe!=1){
								MPframedepLegend->Draw("same");
								double xStatErrorLine=13.06-DeltaXminOVERALL+DeltaXminOVERALL*0.075;
								double StatErrorLineShift=0.1575;
								double StatErrorLineLength=(errorLegendY2-errorLegendY1)*0.25;
								double yMeanStatErrorLine=yMin+(yMax-yMin)*errorLegendY2-(errorLegendY2-errorLegendY1)*StatErrorLineShift;


								TLine* StatErrorLine = new TLine( xStatErrorLine, yMeanStatErrorLine-StatErrorLineLength/2., 
										xStatErrorLine ,yMeanStatErrorLine+StatErrorLineLength/2.);
								StatErrorLine->SetLineWidth( 1 );
								StatErrorLine->SetLineStyle( 1 );
								StatErrorLine->SetLineColor( kBlack );
								StatErrorLine->Draw( "same" );

							}
							if(iStateMP==1&&iPanel==3&&MPframe==1){
								MPframedepLegend->Draw("same");
								double xStatErrorLine=13.06-DeltaXminOVERALL+DeltaXminOVERALL*0.075;
								double StatErrorLineShift=0.1575;
								double StatErrorLineLength=(errorLegendY2-errorLegendY1)*0.25;
								double yMeanStatErrorLine=yMin+(yMax-yMin)*errorLegendY2-(errorLegendY2-errorLegendY1)*StatErrorLineShift;


								TLine* StatErrorLine = new TLine( xStatErrorLine, yMeanStatErrorLine-StatErrorLineLength/2., 
										xStatErrorLine ,yMeanStatErrorLine+StatErrorLineLength/2.);
								StatErrorLine->SetLineWidth( 1 );
								StatErrorLine->SetLineStyle( 1 );
								StatErrorLine->SetLineColor( kBlack );
								StatErrorLine->Draw( "same" );

							}

							if(iStateMP==2&&iPanel==1){//FRchange

								double DeltaXRap;
								if(rapBin==1) DeltaXRap=19.5-DeltaXminOVERALL;
								if(rapBin==2) DeltaXRap=14-DeltaXminOVERALL;

								MPlatexYmax=(yMax-yMin)*0.08+yMin;//FRchange
								//DeltaXRap=-DeltaXminOVERALL;//FRchange

								char frameMPtex[200];
								if(MPframe==1) sprintf(frameMPtex,"CS frame");
								if(MPframe==2) sprintf(frameMPtex,"HX frame");
								if(MPframe==3) sprintf(frameMPtex,"PX frame");
								char textStateFrame[200];
								if(rapBin==1) sprintf(textStateFrame,"%s, |#it{y}| < 0.6", frameMPtex);
								if(rapBin==2) sprintf(textStateFrame,"%s, 0.6 < |#it{y}| < 1.2", frameMPtex);
								TLatex *TexStateFrame = new TLatex(MPlatexX+DeltaXRap,MPlatexYmax,textStateFrame);
								TexStateFrame->SetTextSize(CentralsFontSizeMP);
								TexStateFrame->Draw( "same" );

							}


						}//end iState Loop

						if(MPframe==1&&rapBin==1) MPcanvasCS_rap1->cd();
						if(MPframe==2&&rapBin==1) MPcanvasHX_rap1->cd();
						if(MPframe==3&&rapBin==1) MPcanvasPX_rap1->cd();
						if(MPframe==1&&rapBin==2) MPcanvasCS_rap2->cd();
						if(MPframe==2&&rapBin==2) MPcanvasHX_rap2->cd();
						if(MPframe==3&&rapBin==2) MPcanvasPX_rap2->cd();
						if(MPframe==1&&rapBin==3) MPcanvasCS_rap3->cd();
						if(MPframe==2&&rapBin==3) MPcanvasHX_rap3->cd();
						if(MPframe==3&&rapBin==3) MPcanvasPX_rap3->cd();

						whereTexteInPlotX=0.36;
						whereTexteInPlotY=startValCoordY-deltaCoordY-1.35*labelOffsetX;

						TLatex *MPXlabeltext = new TLatex(whereTexteInPlotX,whereTexteInPlotY ,"1");
						MPXlabeltext->SetTextSize(XaxislabelLatexSize);
						MPXlabeltext->SetTextColor(kBlack);
						if(iPanel==nPanels&&!ShiftXminOVERALL) MPXlabeltext->Draw( "same" );

						whereTexteInPlotX+=x_tilde;
						MPXlabeltext = new TLatex(whereTexteInPlotX,whereTexteInPlotY ,"1");
						MPXlabeltext->SetTextSize(XaxislabelLatexSize);
						MPXlabeltext->SetTextColor(kBlack);
						if(iPanel==nPanels&&!ShiftXminOVERALL) MPXlabeltext->Draw( "same" );

						if(iLam==3&&rapBin==1){  sprintf(filename,"%s/FinalResultsCS_rap%d.pdf",FigDir,rapBin);
							if(PlotFinalData) MPcanvasCS_rap1->SaveAs(filename);
							sprintf(filename,"%s/FinalResultsCS_rap%d.C",FigDir,rapBin);
							if(PlotFinalData) MPcanvasCS_rap1->SaveAs(filename);
							MPcanvasCS_rap1->Close();
						}
						if(iLam==9&&rapBin==1){  sprintf(filename,"%s/FinalResultsHX_rap%d.pdf",FigDir,rapBin);
							if(PlotFinalData) MPcanvasHX_rap1->SaveAs(filename);
							sprintf(filename,"%s/FinalResultsHX_rap%d.C",FigDir,rapBin);
							if(PlotFinalData) MPcanvasHX_rap1->SaveAs(filename);
							MPcanvasHX_rap1->Close();
						}
						if(iLam==15&&rapBin==1){  sprintf(filename,"%s/FinalResultsPX_rap%d.pdf",FigDir,rapBin);
							if(PlotFinalData) MPcanvasPX_rap1->SaveAs(filename);
							sprintf(filename,"%s/FinalResultsPX_rap%d.C",FigDir,rapBin);
							if(PlotFinalData) MPcanvasPX_rap1->SaveAs(filename);
							MPcanvasPX_rap1->Close();
						}

						if(iLam==3&&rapBin==2){  sprintf(filename,"%s/FinalResultsCS_rap%d.pdf",FigDir,rapBin);
							if(PlotFinalData) MPcanvasCS_rap2->SaveAs(filename);
							sprintf(filename,"%s/FinalResultsCS_rap%d.C",FigDir,rapBin);
							if(PlotFinalData) MPcanvasCS_rap2->SaveAs(filename);
							MPcanvasCS_rap2->Close();
						}
						if(iLam==9&&rapBin==2){  sprintf(filename,"%s/FinalResultsHX_rap%d.pdf",FigDir,rapBin);
							if(PlotFinalData) MPcanvasHX_rap2->SaveAs(filename);
							sprintf(filename,"%s/FinalResultsHX_rap%d.C",FigDir,rapBin);
							if(PlotFinalData) MPcanvasHX_rap2->SaveAs(filename);
							MPcanvasHX_rap2->Close();
						}
						if(iLam==15&&rapBin==2){  sprintf(filename,"%s/FinalResultsPX_rap%d.pdf",FigDir,rapBin);
							if(PlotFinalData) MPcanvasPX_rap2->SaveAs(filename);
							sprintf(filename,"%s/FinalResultsPX_rap%d.C",FigDir,rapBin);
							if(PlotFinalData) MPcanvasPX_rap2->SaveAs(filename);
							MPcanvasPX_rap2->Close();
						}
					}//end Frame dependent plots

				}


				double TwoTOthreePanelScaleFactor=deltaCoordY*2.+deltaCoordY*lowestBottomMargin/(1.-lowestBottomMargin)+(1-PadCoordYMax);
				//0.83125/3.*2.+0.16875;
				int MPcanvasXpixel_MPnewTilde=MPcanvasXpixelInitial*1.5;
				int MPcanvasYpixel_MPnewTilde=MPcanvasYpixelInitial*TwoTOthreePanelScaleFactor;

				//errorLegendFontSize*=TwoTOthreePanelScaleFactor;

				const int nPanels_MPnew=2;

				double PadCoordYMax_MPnew=PadCoordYMax;
				double deltaCoordY_MPnew=PadCoordYMax_MPnew/(double(nPanels_MPnew-1)+1./(1-lowestBottomMargin));
				double startValCoordY_MPnew=deltaCoordY_MPnew/(1-lowestBottomMargin);
				//double PadCoordY_MPnew[nPanels_MPnew+1]={0.,startValCoordY_MPnew,
				// 	startValCoordY_MPnew+deltaCoordY_MPnew,PadCoordYMax_MPnew};
				double PadCoordY_MPnew[nPanels_MPnew+1]={0.,startValCoordY_MPnew,PadCoordYMax_MPnew};


				cout<<"begin NEW Frame independent plots"<<endl;
				cout<<"if(iLam==6||iLam==12||iLam==18)"<<endl;
				//begin Frame independent plots
				if(NEW_MPplots){
					if(iLam==6||iLam==12||iLam==18){

						int mainframe;
						if(iLam==6) mainframe=1;
						if(iLam==12) mainframe=2;
						if(iLam==18) mainframe=3;

						cout<<"iLam = "<<iLam<<endl;

						if((iLam==6||iLam==12||iLam==18)&&rapBin==1){
							MPcanvasTilde_Psi = new TCanvas("MPcanvasTilde_Psi", "MPcanvasTilde_Psi",MPcanvasXpixel_MPnewTilde,MPcanvasYpixel_MPnewTilde);
							MPcanvasTilde_Psi->SetFillColor(kWhite);
							MPcanvasTilde_Psi->GetFrame()->SetFillColor(kWhite);
							MPcanvasTilde_Psi->GetFrame()->SetBorderSize(0);
						}

						for(int iStateMP=1;iStateMP<4;iStateMP++){
							if(rapBin==1) iPanel=1;
							if(rapBin==2) iPanel=2;

							cout<<"MultiPanel canvas"<<endl;

							MPcanvasTilde_Psi->cd();

							cout<<"MultiPanel pad"<<endl;
							TPad *MPpad;
							MPpad = new TPad("MPpad","MPpad",PadCoordX_newMP[iStateMP-1],PadCoordY_MPnew[nPanels_MPnew-iPanel],PadCoordX_newMP[iStateMP],PadCoordY_MPnew[nPanels_MPnew-iPanel+1]);
							MPpad->Draw();
							MPpad->cd();
							MPpad->SetFillColor(kWhite);
							MPpad->SetFrameFillColor(kWhite);
							MPpad->SetBorderSize(0);
							MPpad->SetLeftMargin(0.);
							if(iStateMP==1) MPpad->SetLeftMargin(Left_margin);
							MPpad->SetRightMargin(0.);
							if(iStateMP==3) MPpad->SetRightMargin(Right_margin);
							MPpad->SetTopMargin(Top_margin+0.0025);
							MPpad->SetBottomMargin(0.0);
							if(iPanel==nPanels_MPnew) MPpad->SetBottomMargin(lowestBottomMargin);


							cout<<"MultiPanel hist"<<endl;
							TH1F *MPhist = new TH1F;
							MPhist = MPcanvasTilde_Psi->DrawFrame(onia::pTRange[rapBin][ptBinMin-1]-DeltaXminOVERALL,yMinMP,onia::pTRange[rapBin][ptBinMax],yMaxMP);

							MPhist->SetXTitle("#it{p}_{T} [GeV]");
							MPhist->GetXaxis()->SetTitleOffset(-1.35);

							MPhist->SetYTitle(axislabel);
							MPhist->GetYaxis()->SetTitleOffset(titleoffset);
							MPhist->GetYaxis()->SetTitleSize(0.);
							if(iPanel==nPanels_MPnew) MPhist->GetYaxis()->SetTitleOffset(titleoffset*1.35);
							if(iPanel==nPanels_MPnew) MPhist->GetYaxis()->SetTitleSize(0.*(1-lowestBottomMargin));

							MPhist->GetYaxis()->SetLabelSize(LabelSize*1.25);
							MPhist->GetXaxis()->SetLabelSize(0.);
							MPhist->GetYaxis()->SetLabelOffset(-0.015);
							MPhist->GetXaxis()->SetLabelOffset(-0.06);

							if(iPanel==nPanels_MPnew) MPhist->GetYaxis()->SetLabelSize(LabelSize*1.25*(1-lowestBottomMargin));
							MPhist->GetXaxis()->SetTitleSize(TitleSize*0.85);
							MPhist->GetXaxis()->SetAxisColor(kWhite);
							MPhist->GetYaxis()->SetAxisColor(kWhite);
							MPhist->GetXaxis()->SetTicks("-");
							MPhist->GetYaxis()->SetTicks("+");

							TLegend* MPframedepLegendError;
							MPframedepLegendError=new TLegend(errorLegendX1,errorLegendY2-(errorLegendY2-errorLegendY1)*(1-lowestBottomMargin),errorLegendX2,errorLegendY2);
							MPframedepLegendError->SetFillColor(0);
							//MPframedepLegendError->SetTextFont(72);
							MPframedepLegendError->SetTextSize(errorLegendFontSize*(1-lowestBottomMargin));
							MPframedepLegendError->SetBorderSize(0);

							char MPframedepLegendEntry[200];


							TGraphAsymmErrors* graphMP1;
							TGraphAsymmErrors* graphMP2;
							TGraphAsymmErrors* graphMP3;

							TGraphAsymmErrors* graphMP1_1sig;
							TGraphAsymmErrors* graphMP1_2sig;
							TGraphAsymmErrors* graphMP1_3sig;
							TGraphAsymmErrors* graphMP2_1sig;
							TGraphAsymmErrors* graphMP2_2sig;
							TGraphAsymmErrors* graphMP2_3sig;
							TGraphAsymmErrors* graphMP3_1sig;
							TGraphAsymmErrors* graphMP3_2sig;
							TGraphAsymmErrors* graphMP3_3sig;

							TLegend* MPtildeLegend;
							MPtildeLegend=new TLegend(0.8,0.75,1.,0.95);
							MPtildeLegend->SetFillColor(0);
							//			MPtildeLegend->SetTextFont(72);
							MPtildeLegend->SetTextSize(0.07);
							MPtildeLegend->SetBorderSize(0);
							char MPtildeLegendEntry[200];

							for(int iFrameMP=1;iFrameMP<4;iFrameMP++){

								char GraphNameMP[200];


								if(iFrameMP==1){
									sprintf(GraphNameMP,"ltilde_CS_rap%d",rapBin);
								}
								if(iFrameMP==2){
									sprintf(GraphNameMP,"ltilde_HX_rap%d",rapBin);
								}
								if(iFrameMP==3){
									sprintf(GraphNameMP,"ltilde_PX_rap%d",rapBin);
								}

								graphMP1 = (TGraphAsymmErrors*) infileMP1->Get(GraphNameMP);
								graphMP2 = (TGraphAsymmErrors*) infileMP2->Get(GraphNameMP);
								graphMP3 = (TGraphAsymmErrors*) infileMP3->Get(GraphNameMP);

								int MarkerDefinitionForThisBin[4][4]={{0,0,0,0},{0,1,2,3},{0,2,1,3},{0,2,3,1}};


								double ptCentreMP[nBinspT];
								double ptCentreErr_lowMP[nBinspT];
								double ptCentreErr_highMP[nBinspT];
								double lmeanMP[nBinspT];
								double lmean_errlowMP[nBinspT];
								double lmean_errhighMP[nBinspT];

								double ShiftTildePlot;
								double ShiftTildePlotZero=0.75;

								if(MarkerDefinitionForThisBin[mainframe][iFrameMP]==1) ShiftTildePlot=0.;
								if(MarkerDefinitionForThisBin[mainframe][iFrameMP]==2) ShiftTildePlot=ShiftTildePlotZero;
								if(MarkerDefinitionForThisBin[mainframe][iFrameMP]==3) ShiftTildePlot=-ShiftTildePlotZero;

								bool RemoveHorizontalErrorBar=true;

								int pt=0;
								for(int ptBin = ptBinMin; ptBin < ptBinMax+1; ptBin++) {

									graphMP1->GetPoint(ptBin-1,ptCentreMP[pt],lmeanMP[pt]);
									ptCentreErr_highMP[pt]=graphMP1->GetErrorXhigh(ptBin-1);
									ptCentreErr_lowMP[pt]=graphMP1->GetErrorXlow(ptBin-1);
									lmean_errhighMP[pt]=graphMP1->GetErrorYhigh(ptBin-1);
									lmean_errlowMP[pt]=graphMP1->GetErrorYlow(ptBin-1);

									ptCentreMP[pt]=ptCentreMP[pt]+ShiftTildePlot;
									ptCentreErr_highMP[pt]=ptCentreErr_highMP[pt]-ShiftTildePlot;
									ptCentreErr_lowMP[pt]=ptCentreErr_lowMP[pt]-ShiftTildePlot;
									if(RemoveHorizontalErrorBar) ptCentreErr_highMP[pt]=0;
									if(RemoveHorizontalErrorBar) ptCentreErr_lowMP[pt]=0;

									pt++;
								}

								graphMP1 = new TGraphAsymmErrors(nBinspT,ptCentreMP,lmeanMP,ptCentreErr_lowMP,ptCentreErr_highMP,lmean_errlowMP,lmean_errhighMP);

								pt=0;
								for(int ptBin = ptBinMin; ptBin < ptBinMax+1; ptBin++) {

									graphMP2->GetPoint(ptBin-1,ptCentreMP[pt],lmeanMP[pt]);
									ptCentreErr_highMP[pt]=graphMP2->GetErrorXhigh(ptBin-1);
									ptCentreErr_lowMP[pt]=graphMP2->GetErrorXlow(ptBin-1);
									lmean_errhighMP[pt]=graphMP2->GetErrorYhigh(ptBin-1);
									lmean_errlowMP[pt]=graphMP2->GetErrorYlow(ptBin-1);

									ptCentreMP[pt]=ptCentreMP[pt]+ShiftTildePlot;
									ptCentreErr_highMP[pt]=ptCentreErr_highMP[pt]-ShiftTildePlot;
									ptCentreErr_lowMP[pt]=ptCentreErr_lowMP[pt]-ShiftTildePlot;
									if(RemoveHorizontalErrorBar) ptCentreErr_highMP[pt]=0;
									if(RemoveHorizontalErrorBar) ptCentreErr_lowMP[pt]=0;

									pt++;
								}

								graphMP2 = new TGraphAsymmErrors(nBinspT,ptCentreMP,lmeanMP,ptCentreErr_lowMP,ptCentreErr_highMP,lmean_errlowMP,lmean_errhighMP);

								pt=0;
								for(int ptBin = ptBinMin; ptBin < ptBinMax+1; ptBin++) {

									graphMP3->GetPoint(ptBin-1,ptCentreMP[pt],lmeanMP[pt]);
									ptCentreErr_highMP[pt]=graphMP3->GetErrorXhigh(ptBin-1);
									ptCentreErr_lowMP[pt]=graphMP3->GetErrorXlow(ptBin-1);
									lmean_errhighMP[pt]=graphMP3->GetErrorYhigh(ptBin-1);
									lmean_errlowMP[pt]=graphMP3->GetErrorYlow(ptBin-1);

									ptCentreMP[pt]=ptCentreMP[pt]+ShiftTildePlot;
									ptCentreErr_highMP[pt]=ptCentreErr_highMP[pt]-ShiftTildePlot;
									ptCentreErr_lowMP[pt]=ptCentreErr_lowMP[pt]-ShiftTildePlot;
									if(RemoveHorizontalErrorBar) ptCentreErr_highMP[pt]=0;
									if(RemoveHorizontalErrorBar) ptCentreErr_lowMP[pt]=0;

									pt++;
								}

								graphMP3 = new TGraphAsymmErrors(nBinspT,ptCentreMP,lmeanMP,ptCentreErr_lowMP,ptCentreErr_highMP,lmean_errlowMP,lmean_errhighMP);

								graphMP1->SetMarkerColor(MarkerColorMP[MarkerDefinitionForThisBin[mainframe][iFrameMP]]);
								graphMP1->SetLineColor(MarkerColorMP[MarkerDefinitionForThisBin[mainframe][iFrameMP]]);
								graphMP1->SetMarkerStyle(MarkerStyleMP[MarkerDefinitionForThisBin[mainframe][iFrameMP]]);
								graphMP1->SetMarkerSize(MarkerSizeMP[MarkerDefinitionForThisBin[mainframe][iFrameMP]]);

								graphMP2->SetMarkerColor(MarkerColorMP[MarkerDefinitionForThisBin[mainframe][iFrameMP]]);
								graphMP2->SetLineColor(MarkerColorMP[MarkerDefinitionForThisBin[mainframe][iFrameMP]]);
								graphMP2->SetMarkerStyle(MarkerStyleMP[MarkerDefinitionForThisBin[mainframe][iFrameMP]]);
								graphMP2->SetMarkerSize(MarkerSizeMP[MarkerDefinitionForThisBin[mainframe][iFrameMP]]);

								graphMP3->SetMarkerColor(MarkerColorMP[MarkerDefinitionForThisBin[mainframe][iFrameMP]]);
								graphMP3->SetLineColor(MarkerColorMP[MarkerDefinitionForThisBin[mainframe][iFrameMP]]);
								graphMP3->SetMarkerStyle(MarkerStyleMP[MarkerDefinitionForThisBin[mainframe][iFrameMP]]);
								graphMP3->SetMarkerSize(MarkerSizeMP[MarkerDefinitionForThisBin[mainframe][iFrameMP]]);


								if(mainframe==1&&iFrameMP==1){ sprintf(MPtildeLegendEntry,"CS"); 
									MPtildeLegend->AddEntry(graphMP1,MPtildeLegendEntry,"p"); }
									if(mainframe!=1&&iFrameMP==1){ sprintf(MPtildeLegendEntry,"CS");
										MPtildeLegend->AddEntry(graphMP1,MPtildeLegendEntry,"p"); }

										if(mainframe==2&&iFrameMP==2){ sprintf(MPtildeLegendEntry,"HX"); 
											MPtildeLegend->AddEntry(graphMP1,MPtildeLegendEntry,"p"); }
											if(mainframe!=2&&iFrameMP==2){ sprintf(MPtildeLegendEntry,"HX"); 
												MPtildeLegend->AddEntry(graphMP1,MPtildeLegendEntry,"p"); }

												if(mainframe==3&&iFrameMP==3){ sprintf(MPtildeLegendEntry,"PX"); 
													MPtildeLegend->AddEntry(graphMP1,MPtildeLegendEntry,"p"); }
													if(mainframe!=3&&iFrameMP==3){ sprintf(MPtildeLegendEntry,"PX"); 
														MPtildeLegend->AddEntry(graphMP1,MPtildeLegendEntry,"p"); }

														if(mainframe==1){
															sprintf(GraphNameMP,"ltilde_CS_rap%d",rapBin);
														}
														if(mainframe==2){
															sprintf(GraphNameMP,"ltilde_HX_rap%d",rapBin);
														}
														if(mainframe==3){
															sprintf(GraphNameMP,"ltilde_PX_rap%d",rapBin);
														}

														graphMP1_1sig = (TGraphAsymmErrors*) infileMP1_1sig->Get(GraphNameMP);
														graphMP1_2sig = (TGraphAsymmErrors*) infileMP1_2sig->Get(GraphNameMP);
														graphMP1_3sig = (TGraphAsymmErrors*) infileMP1_3sig->Get(GraphNameMP);
														graphMP2_1sig = (TGraphAsymmErrors*) infileMP2_1sig->Get(GraphNameMP);
														graphMP2_2sig = (TGraphAsymmErrors*) infileMP2_2sig->Get(GraphNameMP);
														graphMP2_3sig = (TGraphAsymmErrors*) infileMP2_3sig->Get(GraphNameMP);
														graphMP3_1sig = (TGraphAsymmErrors*) infileMP3_1sig->Get(GraphNameMP);
														graphMP3_2sig = (TGraphAsymmErrors*) infileMP3_2sig->Get(GraphNameMP);
														graphMP3_3sig = (TGraphAsymmErrors*) infileMP3_3sig->Get(GraphNameMP);

														int ptMP;
														int nBinsMP;
														double ptCentre_MP[nBinspT];
														double lmean_MP[nBinspT];
														double lmean_errlow_MP[nBinspT];
														double lmean_errhigh_MP[nBinspT];
														double ptCentre_errlow_MP[nBinspT];
														double ptCentre_errhigh_MP[nBinspT];

														nBinsMP=graphMP1_1sig->GetN();
														ptMP=0;
														for(int ptBinMP=1;ptBinMP<nBinsMP+1;ptBinMP++){
															graphMP1_1sig->GetPoint(ptBinMP-1,ptCentre_MP[ptMP],lmean_MP[ptMP]);
															lmean_errhigh_MP[ptMP]=graphMP1_1sig->GetErrorYhigh(ptBinMP-1);
															lmean_errlow_MP[ptMP]=graphMP1_1sig->GetErrorYlow(ptBinMP-1);
															ptCentre_errhigh_MP[ptMP]=graphMP1_1sig->GetErrorXhigh(ptBinMP-1);
															ptCentre_errlow_MP[ptMP]=graphMP1_1sig->GetErrorXlow(ptBinMP-1);
															/// Alter TGraph
															ptCentre_errhigh_MP[ptMP]=ColordBandWidth;
															ptCentre_errlow_MP[ptMP]=ColordBandWidth;
															ptMP++;
														}
														graphMP1_1sig = new TGraphAsymmErrors(nBinsMP,ptCentre_MP,lmean_MP,ptCentre_errlow_MP,ptCentre_errhigh_MP,lmean_errlow_MP,lmean_errhigh_MP);

														nBinsMP=graphMP1_2sig->GetN();
														ptMP=0;
														for(int ptBinMP=1;ptBinMP<nBinsMP+1;ptBinMP++){
															graphMP1_2sig->GetPoint(ptBinMP-1,ptCentre_MP[ptMP],lmean_MP[ptMP]);
															lmean_errhigh_MP[ptMP]=graphMP1_2sig->GetErrorYhigh(ptBinMP-1);
															lmean_errlow_MP[ptMP]=graphMP1_2sig->GetErrorYlow(ptBinMP-1);
															ptCentre_errhigh_MP[ptMP]=graphMP1_2sig->GetErrorXhigh(ptBinMP-1);
															ptCentre_errlow_MP[ptMP]=graphMP1_2sig->GetErrorXlow(ptBinMP-1);
															/// Alter TGraph
															ptCentre_errhigh_MP[ptMP]=ColordBandWidth;
															ptCentre_errlow_MP[ptMP]=ColordBandWidth;
															ptMP++;
														}
														graphMP1_2sig = new TGraphAsymmErrors(nBinsMP,ptCentre_MP,lmean_MP,ptCentre_errlow_MP,ptCentre_errhigh_MP,lmean_errlow_MP,lmean_errhigh_MP);

														nBinsMP=graphMP1_3sig->GetN();
														ptMP=0;
														for(int ptBinMP=1;ptBinMP<nBinsMP+1;ptBinMP++){
															graphMP1_3sig->GetPoint(ptBinMP-1,ptCentre_MP[ptMP],lmean_MP[ptMP]);
															lmean_errhigh_MP[ptMP]=graphMP1_3sig->GetErrorYhigh(ptBinMP-1);
															lmean_errlow_MP[ptMP]=graphMP1_3sig->GetErrorYlow(ptBinMP-1);
															ptCentre_errhigh_MP[ptMP]=graphMP1_3sig->GetErrorXhigh(ptBinMP-1);
															ptCentre_errlow_MP[ptMP]=graphMP1_3sig->GetErrorXlow(ptBinMP-1);
															/// Alter TGraph
															ptCentre_errhigh_MP[ptMP]=ColordBandWidth;
															ptCentre_errlow_MP[ptMP]=ColordBandWidth;
															ptMP++;
														}
														graphMP1_3sig = new TGraphAsymmErrors(nBinsMP,ptCentre_MP,lmean_MP,ptCentre_errlow_MP,ptCentre_errhigh_MP,lmean_errlow_MP,lmean_errhigh_MP);


														nBinsMP=graphMP2_1sig->GetN();
														ptMP=0;
														for(int ptBinMP=1;ptBinMP<nBinsMP+1;ptBinMP++){
															graphMP2_1sig->GetPoint(ptBinMP-1,ptCentre_MP[ptMP],lmean_MP[ptMP]);
															lmean_errhigh_MP[ptMP]=graphMP2_1sig->GetErrorYhigh(ptBinMP-1);
															lmean_errlow_MP[ptMP]=graphMP2_1sig->GetErrorYlow(ptBinMP-1);
															ptCentre_errhigh_MP[ptMP]=graphMP2_1sig->GetErrorXhigh(ptBinMP-1);
															ptCentre_errlow_MP[ptMP]=graphMP2_1sig->GetErrorXlow(ptBinMP-1);
															/// Alter TGraph
															ptCentre_errhigh_MP[ptMP]=ColordBandWidth;
															ptCentre_errlow_MP[ptMP]=ColordBandWidth;
															ptMP++;
														}
														graphMP2_1sig = new TGraphAsymmErrors(nBinsMP,ptCentre_MP,lmean_MP,ptCentre_errlow_MP,ptCentre_errhigh_MP,lmean_errlow_MP,lmean_errhigh_MP);

														nBinsMP=graphMP2_2sig->GetN();
														ptMP=0;
														for(int ptBinMP=1;ptBinMP<nBinsMP+1;ptBinMP++){
															graphMP2_2sig->GetPoint(ptBinMP-1,ptCentre_MP[ptMP],lmean_MP[ptMP]);
															lmean_errhigh_MP[ptMP]=graphMP2_2sig->GetErrorYhigh(ptBinMP-1);
															lmean_errlow_MP[ptMP]=graphMP2_2sig->GetErrorYlow(ptBinMP-1);
															ptCentre_errhigh_MP[ptMP]=graphMP2_2sig->GetErrorXhigh(ptBinMP-1);
															ptCentre_errlow_MP[ptMP]=graphMP2_2sig->GetErrorXlow(ptBinMP-1);
															/// Alter TGraph
															ptCentre_errhigh_MP[ptMP]=ColordBandWidth;
															ptCentre_errlow_MP[ptMP]=ColordBandWidth;
															ptMP++;
														}
														graphMP2_2sig = new TGraphAsymmErrors(nBinsMP,ptCentre_MP,lmean_MP,ptCentre_errlow_MP,ptCentre_errhigh_MP,lmean_errlow_MP,lmean_errhigh_MP);

														nBinsMP=graphMP2_3sig->GetN();
														ptMP=0;
														for(int ptBinMP=1;ptBinMP<nBinsMP+1;ptBinMP++){
															graphMP2_3sig->GetPoint(ptBinMP-1,ptCentre_MP[ptMP],lmean_MP[ptMP]);
															lmean_errhigh_MP[ptMP]=graphMP2_3sig->GetErrorYhigh(ptBinMP-1);
															lmean_errlow_MP[ptMP]=graphMP2_3sig->GetErrorYlow(ptBinMP-1);
															ptCentre_errhigh_MP[ptMP]=graphMP2_3sig->GetErrorXhigh(ptBinMP-1);
															ptCentre_errlow_MP[ptMP]=graphMP2_3sig->GetErrorXlow(ptBinMP-1);
															/// Alter TGraph
															ptCentre_errhigh_MP[ptMP]=ColordBandWidth;
															ptCentre_errlow_MP[ptMP]=ColordBandWidth;
															ptMP++;
														}
														graphMP2_3sig = new TGraphAsymmErrors(nBinsMP,ptCentre_MP,lmean_MP,ptCentre_errlow_MP,ptCentre_errhigh_MP,lmean_errlow_MP,lmean_errhigh_MP);


														nBinsMP=graphMP3_1sig->GetN();
														ptMP=0;
														for(int ptBinMP=1;ptBinMP<nBinsMP+1;ptBinMP++){
															graphMP3_1sig->GetPoint(ptBinMP-1,ptCentre_MP[ptMP],lmean_MP[ptMP]);
															lmean_errhigh_MP[ptMP]=graphMP3_1sig->GetErrorYhigh(ptBinMP-1);
															lmean_errlow_MP[ptMP]=graphMP3_1sig->GetErrorYlow(ptBinMP-1);
															ptCentre_errhigh_MP[ptMP]=graphMP3_1sig->GetErrorXhigh(ptBinMP-1);
															ptCentre_errlow_MP[ptMP]=graphMP3_1sig->GetErrorXlow(ptBinMP-1);
															/// Alter TGraph
															ptCentre_errhigh_MP[ptMP]=ColordBandWidth;
															ptCentre_errlow_MP[ptMP]=ColordBandWidth;
															ptMP++;
														}
														graphMP3_1sig = new TGraphAsymmErrors(nBinsMP,ptCentre_MP,lmean_MP,ptCentre_errlow_MP,ptCentre_errhigh_MP,lmean_errlow_MP,lmean_errhigh_MP);

														nBinsMP=graphMP3_2sig->GetN();
														ptMP=0;
														for(int ptBinMP=1;ptBinMP<nBinsMP+1;ptBinMP++){
															graphMP3_2sig->GetPoint(ptBinMP-1,ptCentre_MP[ptMP],lmean_MP[ptMP]);
															lmean_errhigh_MP[ptMP]=graphMP3_2sig->GetErrorYhigh(ptBinMP-1);
															lmean_errlow_MP[ptMP]=graphMP3_2sig->GetErrorYlow(ptBinMP-1);
															ptCentre_errhigh_MP[ptMP]=graphMP3_2sig->GetErrorXhigh(ptBinMP-1);
															ptCentre_errlow_MP[ptMP]=graphMP3_2sig->GetErrorXlow(ptBinMP-1);
															/// Alter TGraph
															ptCentre_errhigh_MP[ptMP]=ColordBandWidth;
															ptCentre_errlow_MP[ptMP]=ColordBandWidth;
															ptMP++;
														}
														graphMP3_2sig = new TGraphAsymmErrors(nBinsMP,ptCentre_MP,lmean_MP,ptCentre_errlow_MP,ptCentre_errhigh_MP,lmean_errlow_MP,lmean_errhigh_MP);

														nBinsMP=graphMP3_3sig->GetN();
														ptMP=0;
														for(int ptBinMP=1;ptBinMP<nBinsMP+1;ptBinMP++){
															graphMP3_3sig->GetPoint(ptBinMP-1,ptCentre_MP[ptMP],lmean_MP[ptMP]);
															lmean_errhigh_MP[ptMP]=graphMP3_3sig->GetErrorYhigh(ptBinMP-1);
															lmean_errlow_MP[ptMP]=graphMP3_3sig->GetErrorYlow(ptBinMP-1);
															ptCentre_errhigh_MP[ptMP]=graphMP3_3sig->GetErrorXhigh(ptBinMP-1);
															ptCentre_errlow_MP[ptMP]=graphMP3_3sig->GetErrorXlow(ptBinMP-1);
															/// Alter TGraph
															ptCentre_errhigh_MP[ptMP]=ColordBandWidth;
															ptCentre_errlow_MP[ptMP]=ColordBandWidth;
															ptMP++;
														}
														graphMP3_3sig = new TGraphAsymmErrors(nBinsMP,ptCentre_MP,lmean_MP,ptCentre_errlow_MP,ptCentre_errhigh_MP,lmean_errlow_MP,lmean_errhigh_MP);


														graphMP1_1sig->SetFillColor(OneSigColor);
														graphMP1_1sig->SetFillStyle(1001);
														graphMP1_2sig->SetFillColor(TwoSigColor);
														graphMP1_2sig->SetFillStyle(1001);
														graphMP1_3sig->SetFillColor(ThreeSigColor);
														graphMP1_3sig->SetFillStyle(1001);

														graphMP2_1sig->SetFillColor(OneSigColor);
														graphMP2_1sig->SetFillStyle(1001);
														graphMP2_2sig->SetFillColor(TwoSigColor);
														graphMP2_2sig->SetFillStyle(1001);
														graphMP2_3sig->SetFillColor(ThreeSigColor);
														graphMP2_3sig->SetFillStyle(1001);

														graphMP3_1sig->SetFillColor(OneSigColor);
														graphMP3_1sig->SetFillStyle(1001);
														graphMP3_2sig->SetFillColor(TwoSigColor);
														graphMP3_2sig->SetFillStyle(1001);
														graphMP3_3sig->SetFillColor(ThreeSigColor);
														graphMP3_3sig->SetFillStyle(1001);

														if(iStateMP==1){

															if(iFrameMP==1){
																graphMP1_3sig->Draw("2");
																graphMP1_2sig->Draw("2");
																graphMP1_1sig->Draw("2");
																if(mainframe==1) graphMP1->Draw("[]");
																if(mainframe==1) graphMP1->Draw(drawGraphStyle);
															}
															if(iFrameMP==2){
																if(mainframe==2) graphMP1->Draw("[]");
																if(mainframe==2) graphMP1->Draw(drawGraphStyle);
															}
															if(iFrameMP==3){
																if(mainframe==3) graphMP1->Draw("[]");
																if(mainframe==3) graphMP1->Draw(drawGraphStyle);
															}

															if(iFrameMP==1){
																if(mainframe!=1) graphMP1->Draw("PX");
															}
															if(iFrameMP==2){
																if(mainframe!=2) graphMP1->Draw("PX");
															}
															if(iFrameMP==3){
																if(mainframe!=3) graphMP1->Draw("PX");
															}

														}
														if(iStateMP==2){

															if(iFrameMP==1){
																graphMP2_3sig->Draw("2");
																graphMP2_2sig->Draw("2");
																graphMP2_1sig->Draw("2");
																if(mainframe==1) graphMP2->Draw("[]");
																if(mainframe==1) graphMP2->Draw(drawGraphStyle);
															}
															if(iFrameMP==2){
																if(mainframe==2) graphMP2->Draw("[]");
																if(mainframe==2) graphMP2->Draw(drawGraphStyle);
															}
															if(iFrameMP==3){
																if(mainframe==3) graphMP2->Draw("[]");
																if(mainframe==3) graphMP2->Draw(drawGraphStyle);
															}

															if(iFrameMP==1){
																if(mainframe!=1) graphMP2->Draw("PX");
															}
															if(iFrameMP==2){
																if(mainframe!=2) graphMP2->Draw("PX");
															}
															if(iFrameMP==3){
																if(mainframe!=3) graphMP2->Draw("PX");
															}

														}
														if(iStateMP==3){

															if(iFrameMP==1){
																graphMP3_3sig->Draw("2");
																graphMP3_2sig->Draw("2");
																graphMP3_1sig->Draw("2");
																if(mainframe==1) graphMP3->Draw("[]");
																if(mainframe==1) graphMP3->Draw(drawGraphStyle);
															}
															if(iFrameMP==2){
																if(mainframe==2) graphMP3->Draw("[]");
																if(mainframe==2) graphMP3->Draw(drawGraphStyle);
															}
															if(iFrameMP==3){
																if(mainframe==3) graphMP3->Draw("[]");
																if(mainframe==3) graphMP3->Draw(drawGraphStyle);
															}

															if(iFrameMP==1){
																if(mainframe!=1) graphMP3->Draw("PX");
															}
															if(iFrameMP==2){
																if(mainframe!=2) graphMP3->Draw("PX");
															}
															if(iFrameMP==3){
																if(mainframe!=3) graphMP3->Draw("PX");
															}


														}
							}


							TGaxis *axisMPY1 = new TGaxis(onia::pTRange[rapBin][ptBinMin-1]-DeltaXminOVERALL,yMinMP,onia::pTRange[rapBin][ptBinMin-1]-DeltaXminOVERALL,yMaxMP,yMinMP,yMaxMP,AxisDivisions,"-US");
							axisMPY1->SetTickSize(ticksize);
							if(iPanel==nPanels_MPnew) axisMPY1->SetTickSize(ticksize/(1-lowestBottomMargin));
							axisMPY1->Draw("same");

							TGaxis *axisMPY2 = new TGaxis(onia::pTRange[rapBin][ptBinMax],yMinMP,onia::pTRange[rapBin][ptBinMax],yMaxMP,yMinMP,yMaxMP,AxisDivisions,"+US");
							axisMPY2->SetTickSize(ticksize);
							if(iPanel==nPanels_MPnew) axisMPY2->SetTickSize(ticksize/(1-lowestBottomMargin));
							axisMPY2->Draw("same");


							double deltaTrickAxisMax_MPnew;
							deltaTrickAxisMax_MPnew=-0.001;
							//if(iStateMP==3) deltaTrickAxisMax_MPnew=+0.001;

							TGaxis *axisM3S1 = new TGaxis(onia::pTRange[rapBin][ptBinMin-1]-DeltaXminOVERALL,yMinMP,onia::pTRange[rapBin][ptBinMax],yMinMP,onia::pTRange[rapBin][ptBinMin-1]-DeltaXminOVERALL+deltaTrickAxisMin,onia::pTRange[rapBin][ptBinMax]+deltaTrickAxisMax_MPnew,AxisDivisions,"+S");
							axisM3S1->SetTickSize(ticksize*2);
							if(iPanel==nPanels_MPnew) axisM3S1->SetLabelSize(LabelSize);
							if(iPanel<nPanels_MPnew) axisM3S1->SetLabelSize(0);
							axisM3S1->SetLabelOffset(labelOffsetX);
							if(iPanel==nPanels_MPnew) axisM3S1->SetTickSize(ticksize/(1-lowestBottomMargin));
							axisM3S1->Draw("same");

							TGaxis *axisM3S2 = new TGaxis(onia::pTRange[rapBin][ptBinMin-1]-DeltaXminOVERALL,yMaxMP,onia::pTRange[rapBin][ptBinMax],yMaxMP,onia::pTRange[rapBin][ptBinMin-1]-DeltaXminOVERALL+deltaTrickAxisMin,onia::pTRange[rapBin][ptBinMax]+deltaTrickAxisMax_MPnew,AxisDivisions,"-US");
							axisM3S2->SetTickSize(ticksize*2);
							if(iPanel==nPanels_MPnew) axisM3S2->SetTickSize(ticksize/(1-lowestBottomMargin));
							axisM3S2->Draw("same");

							whereTexteInPlotX=XtitlePosition;
							whereTexteInPlotY=(yMaxMP+yMinMP)/2.-(yMaxMP-yMinMP)*XtitlePositionYshift;

							char axistitleMPtilde[200];
							sprintf(axistitleMPtilde,"#tilde{#lambda}");
							if(iPanel==nPanels_MPnew) YaxistitleLatexSize=YaxistitleLatexSize*(1-lowestBottomMargin);
							TLatex *MPYtitletext = new TLatex(whereTexteInPlotX,whereTexteInPlotY ,axistitleMPtilde);
							MPYtitletext->SetTextSize(YaxistitleLatexSize);
							if(iPanel==nPanels_MPnew) MPYtitletext->SetTextSize(YaxistitleLatexSize*(1-lowestBottomMargin));
							MPYtitletext->SetTextColor(kBlack);
							MPYtitletext->PaintLatex(whereTexteInPlotX-DeltaXminOVERALL,whereTexteInPlotY, YtitleAngle, YaxistitleLatexSize, axislabel);
							MPYtitletext->Draw( "same" );



							double increaseSize=1.25;
							double SpecialShiftUpsilonLabel=0.9;
							double SpecialShiftUpsilonLabelx=-3;

							char texTexMP[200];
							if(rapBin==1) sprintf(texTexMP,"#psi(%dS), |#it{y}| < 0.6", iStateMP);
							if(rapBin==2) sprintf(texTexMP,"#psi(%dS), 0.6 < |#it{y}| < 1.2", iStateMP);
							TLatex *textMP = new TLatex(xRapTextTilde-DeltaXminOVERALL+SpecialShiftUpsilonLabelx,yMin+(yMax-yMin)*yRapText*0.92*SpecialShiftUpsilonLabel,texTexMP);
							textMP->SetTextSize(textSizeRap*increaseSize);
							if(iPanel==nPanels_MPnew) textMP->SetTextSize(textSizeRap*(1-lowestBottomMargin)*increaseSize);
							textMP->Draw( "same" );

							char abcdef[200];
							if(iStateMP==1&&iPanel==1) sprintf(abcdef,"a)");
							if(iStateMP==1&&iPanel==2) sprintf(abcdef,"b)");
							if(iStateMP==2&&iPanel==1) sprintf(abcdef,"c)");
							if(iStateMP==2&&iPanel==2) sprintf(abcdef,"d)");
							if(iStateMP==3&&iPanel==1) sprintf(abcdef,"e)");
							if(iStateMP==3&&iPanel==2) sprintf(abcdef,"f)");
							cout<<abcdef<<endl;
							TLatex *tex_abcdef = new TLatex(xabcdefText-DeltaXminOVERALL,yMin+(yMax-yMin)*yRapText*0.92,abcdef);
							tex_abcdef->SetTextSize(textSizeRap);
							if(iPanel==nPanels_MPnew) tex_abcdef->SetTextSize(textSizeRap*(1-lowestBottomMargin));
							//tex_abcdef->Draw( "same" );

							if(PlotFinalData&&DrawLatexStuff){

								TLine* extreme0MP = new TLine( PlotpTMin-DeltaXminOVERALL, 0, PlotpTMax ,0);
								extreme0MP->SetLineWidth( 1 );
								extreme0MP->SetLineStyle( 2 );
								extreme0MP->SetLineColor( kBlack );
								extreme0MP->Draw( "same" );

								TLine* extreme1MP = new TLine( PlotpTMin-DeltaXminOVERALL, 1, PlotpTMax , 1);
								extreme1MP->SetLineWidth( 1 );
								extreme1MP->SetLineStyle( 2 );
								extreme1MP->SetLineColor( kBlack );
								if(iLam==1||iLam==7||iLam==13) extreme1MP->Draw( "same" );

								TLine* extreme2MP = new TLine( PlotpTMin-DeltaXminOVERALL, -1, PlotpTMax ,-1);
								extreme2MP->SetLineWidth( 1 );
								extreme2MP->SetLineStyle( 2 );
								extreme2MP->SetLineColor( kBlack );
								if(iLam==1||iLam==7||iLam==13) extreme2MP->Draw( "same" );
								//if(iLam==6||iLam==12||iLam==18) extreme2MP->Draw( "same" );

							}


							if(iStateMP==1&&iPanel==1){
								cout<<"DRAW CMS preliminary Latex"<<endl;
								char text[200];
								sprintf(text,"CMS     pp      #sqrt{s} = 7 TeV     L = 4.9 fb^{-1}");
								TLatex *CentralsText1MP = new TLatex(MPlatexX-DeltaXminOVERALL,MPlatexYmax,text);
								CentralsText1MP->SetTextSize(CentralsFontSizeMP);
								CentralsText1MP->Draw( "same" );
								sprintf(text,"preliminary");
								CentralsText1MP = new TLatex(MPlatexX-DeltaXminOVERALL,(yMax-yMin)*0.75+yMin,text);
								CentralsText1MP->SetTextSize(CentralsFontSizeMP);
								if(DrawPreliminary) CentralsText1MP->Draw( "same" );
								//sprintf(text,"L = 4.9 fb^{-1}");
								//TLatex *CentralsText2MP = new TLatex(MPlatexX,MPlatexYmax-2*MPlatexDeltaYmax,text);
								//CentralsText2MP->SetTextSize(CentralsFontSizeMP);
								//CentralsText2MP->Draw( "same" );
								//sprintf(text,"pp    #sqrt{s} = 7 TeV");
								//TLatex *CentralsText3MP = new TLatex(MPlatexX,MPlatexYmax-MPlatexDeltaYmax,text);
								//CentralsText3MP->SetTextSize(CentralsFontSizeMP);
								//CentralsText3MP->Draw( "same" );

							}

							if(iStateMP==1&&iPanel==2){
								graphMP1_1sig->SetLineColor(OneSigColor);
								graphMP1_2sig->SetLineColor(TwoSigColor);
								graphMP1_3sig->SetLineColor(ThreeSigColor);

								TGraphAsymmErrors *legendPhantom = (TGraphAsymmErrors*) infileMP3_3sig->Get("ltilde_CS_rap1");

								legendPhantom->SetMarkerColor(MarkerColorMP[0]);
								legendPhantom->SetLineColor(MarkerColorMP[0]);
								legendPhantom->SetMarkerStyle(MarkerStyleMP[1]);
								legendPhantom->SetMarkerSize(MarkerSizeMP[0]);


								if(PlotAlteredPPDResults) sprintf(MPframedepLegendEntry,"Stat. uncert., 68.3 %% CL");
								if(!PlotAlteredPPDResults) sprintf(MPframedepLegendEntry,"Tot. sys. uncert.");
								MPframedepLegendError->AddEntry(legendPhantom,MPframedepLegendEntry,"p");
								if(PlotAlteredPPDResults) sprintf(MPframedepLegendEntry,"Tot. uncert., 68.3 %% CL");
								if(!PlotAlteredPPDResults) sprintf(MPframedepLegendEntry,"68.3 %% CL");
								MPframedepLegendError->AddEntry(graphMP1_1sig,MPframedepLegendEntry,"f");
								if(PlotAlteredPPDResults) sprintf(MPframedepLegendEntry,"Tot. uncert., 95.5 %% CL");
								if(!PlotAlteredPPDResults) sprintf(MPframedepLegendEntry,"95.5 %% CL");
								MPframedepLegendError->AddEntry(graphMP1_2sig,MPframedepLegendEntry,"f");
								if(PlotAlteredPPDResults) sprintf(MPframedepLegendEntry,"Tot. uncert., 99.7 %% CL");
								if(!PlotAlteredPPDResults) sprintf(MPframedepLegendEntry,"99.7 %% CL");
								MPframedepLegendError->AddEntry(graphMP1_3sig,MPframedepLegendEntry,"f");

								MPtildeLegend->Draw("same");
								MPframedepLegendError->Draw("same");

								double xStatErrorLine=13.06-DeltaXminOVERALL+DeltaXminOVERALL*0.075;
								double StatErrorLineShift=0.75;
								double errorLegendY1Tilde=errorLegendY2-(errorLegendY2-errorLegendY1)*(1-lowestBottomMargin)+0.05;
								double StatErrorLineLength=(errorLegendY2-errorLegendY1Tilde)*0.95;
								double yMeanStatErrorLine=yMin+(yMax-yMin)*errorLegendY2-(errorLegendY2-errorLegendY1Tilde)*StatErrorLineShift;


								TLine* StatErrorLine = new TLine( xStatErrorLine, yMeanStatErrorLine-StatErrorLineLength/2., xStatErrorLine ,yMeanStatErrorLine+StatErrorLineLength/2.);
								StatErrorLine->SetLineWidth( 1 );
								StatErrorLine->SetLineStyle( 1 );
								StatErrorLine->SetLineColor( kBlack );
								StatErrorLine->Draw( "same" );

							}

							//if(rapBin==2&&iPanel==1){
							//	char frameMPtex[200];
							//	if(MPframe==1) sprintf(frameMPtex,"CS frame");
							//	if(MPframe==2) sprintf(frameMPtex,"HX frame");
							//	if(MPframe==3) sprintf(frameMPtex,"PX frame");
							//	char textStateFrame[200];
							//	sprintf(textStateFrame,"%s", frameMPtex);
							//	TLatex *TexStateFrame = new TLatex(MPlatexX,MPlatexYmax,textStateFrame);
							//	TexStateFrame->SetTextSize(CentralsFontSizeMP);
							//	TexStateFrame->Draw( "same" );
							//}

						}//end iStateMP loop

						MPcanvasTilde_Psi->cd();

						//whereTexteInPlotX=0.488;
						//whereTexteInPlotY=startValCoordY-deltaCoordY-1.425*labelOffsetX;
						//TLatex *M3Slabeltext = new TLatex(whereTexteInPlotX,whereTexteInPlotY ,"1");
						//M3Slabeltext->SetTextSize(XaxislabelLatexSize);
						//M3Slabeltext->SetTextColor(kBlack);
						//if(iPanel==nPanels_MPnew) M3Slabeltext->Draw( "same" );

						whereTexteInPlotX=0.36035;
						whereTexteInPlotY=startValCoordY-deltaCoordY-0.006;

						TLatex *MPXlabeltext = new TLatex(whereTexteInPlotX,whereTexteInPlotY ,"1");
						MPXlabeltext->SetTextSize(XaxislabelLatexSize/TwoTOthreePanelScaleFactor);
						MPXlabeltext->SetTextColor(kBlack);
						if(iPanel==nPanels_MPnew&&!ShiftXminOVERALL) MPXlabeltext->Draw( "same" );

						whereTexteInPlotX+=x_tilde;
						MPXlabeltext = new TLatex(whereTexteInPlotX,whereTexteInPlotY ,"1");
						MPXlabeltext->SetTextSize(XaxislabelLatexSize/TwoTOthreePanelScaleFactor);
						MPXlabeltext->SetTextColor(kBlack);
						if(iPanel==nPanels_MPnew&&!ShiftXminOVERALL) MPXlabeltext->Draw( "same" );


						if((iLam==6||iLam==12||iLam==18)&&rapBin==2){
							if(mainframe==1) sprintf(filename,"%s/FinalResultsTildeCS.pdf",FigDir);
							if(mainframe==2) sprintf(filename,"%s/FinalResultsTildeHX.pdf",FigDir);
							if(mainframe==3) sprintf(filename,"%s/FinalResultsTildePX.pdf",FigDir);
							if(PlotFinalData) MPcanvasTilde_Psi->SaveAs(filename);
							if(mainframe==1) sprintf(filename,"%s/FinalResultsTildeCS.C",FigDir);
							if(mainframe==2) sprintf(filename,"%s/FinalResultsTildeHX.C",FigDir);
							if(mainframe==3) sprintf(filename,"%s/FinalResultsTildePX.C",FigDir);
							if(PlotFinalData) MPcanvasTilde_Psi->SaveAs(filename);
							MPcanvasTilde_Psi->Close();
						}
					}//end Frame independent plots NEW
				}


				double TwoTOthreePanelScaleFactorCDF=deltaCoordY*1+deltaCoordY*lowestBottomMargin/(1-lowestBottomMargin)+(1-PadCoordYMax);//0.83125/3.*1.+0.16875;
				int MPcanvasXpixel_MPnewTildeCDF=MPcanvasXpixelInitial*1.5;
				int MPcanvasYpixel_MPnewTildeCDF=MPcanvasYpixelInitial*TwoTOthreePanelScaleFactorCDF;

				//			errorLegendFontSize*=TwoTOthreePanelScaleFactor;

				const int nPanels_MPnewCDF=1;

				double PadCoordYMax_MPnewCDF=PadCoordYMax;
				double deltaCoordY_MPnewCDF=PadCoordYMax_MPnewCDF/(double(nPanels_MPnewCDF-1)+1./(1-lowestBottomMargin));
				double startValCoordY_MPnewCDF=deltaCoordY_MPnewCDF/(1-lowestBottomMargin);
				double PadCoordY_MPnewCDF[nPanels_MPnewCDF+1]={0.,PadCoordYMax_MPnewCDF};

				double PlotMattpTMin=0.;

				cout<<"begin NEW CDF MP plots"<<endl;
				cout<<"if(iLam==1||iLam==7){"<<endl;
				//begin Frame independent plots
				if(NEW_MPplots&&PlotMatt){
					if(iLam<1000&&rapBin==1){

						/*yMin=-1;
							yMax=1;

							if(iLam==2||iLam==8||iLam==14||iLam==3||iLam==9||iLam==15){
							yMin=-0.5;
							yMax=0.5;
							}

							if(iLam==6||iLam==12||iLam==18){
							yMin=-1;
							yMax=1;
							}

							yMinMP=yMin+0.01;
							yMaxMP=yMax-0.01;
							*/
						if(iLam==6||iLam==12||iLam==18||iLam==1||iLam==7||iLam==13){
							yMin=-1;
							yMax=1;
							yMinMP=yMin+0.01;
							yMaxMP=yMax-0.01;
						}

						int mainframe;

						if(iLam>0&&iLam<7) mainframe=1;
						if(iLam>6&&iLam<13) mainframe=2;
						if(iLam>12&&iLam<19) mainframe=3;

						cout<<"iLam = "<<iLam<<endl;
						cout<<"axislabel = "<<axislabel<<endl;

						MPcanvasCDF = new TCanvas("MPcanvasCDF", "MPcanvasCDF",MPcanvasXpixel_MPnewTildeCDF,MPcanvasYpixel_MPnewTildeCDF);
						MPcanvasCDF->SetFillColor(kWhite);
						MPcanvasCDF->GetFrame()->SetFillColor(kWhite);
						MPcanvasCDF->GetFrame()->SetBorderSize(0);


						for(int iStateMP=1;iStateMP<4;iStateMP++){
							iPanel=1;

							cout<<"MultiPanel canvas"<<endl;

							MPcanvasCDF->cd();

							cout<<"MultiPanel pad"<<endl;
							TPad *MPpad;
							MPpad = new TPad("MPpad","MPpad",PadCoordX_newMP[iStateMP-1],PadCoordY_MPnewCDF[nPanels_MPnewCDF-iPanel],PadCoordX_newMP[iStateMP],PadCoordY_MPnewCDF[nPanels_MPnewCDF-iPanel+1]);
							MPpad->Draw();
							MPpad->cd();
							MPpad->SetFillColor(kWhite);
							MPpad->SetFrameFillColor(kWhite);
							MPpad->SetBorderSize(0);
							MPpad->SetLeftMargin(0.);
							if(iStateMP==1) MPpad->SetLeftMargin(Left_margin);
							MPpad->SetRightMargin(0.);
							if(iStateMP==3) MPpad->SetRightMargin(Right_margin);
							MPpad->SetTopMargin(Top_margin+0.0025);
							MPpad->SetBottomMargin(0.0);
							if(iPanel==nPanels_MPnewCDF) MPpad->SetBottomMargin(lowestBottomMargin);


							cout<<"MultiPanel hist"<<endl;
							TH1F *MPhist = new TH1F;
							MPhist = MPcanvasCDF->DrawFrame(PlotMattpTMin,yMinMP,onia::pTRange[rapBin][ptBinMax],yMaxMP);

							MPhist->SetXTitle("#it{p}_{T} [GeV]");
							MPhist->GetXaxis()->SetTitleOffset(-1.35);

							MPhist->SetYTitle(axislabel);
							MPhist->GetYaxis()->SetTitleOffset(titleoffset);
							MPhist->GetYaxis()->SetTitleSize(0.);
							if(iPanel==nPanels_MPnewCDF) MPhist->GetYaxis()->SetTitleOffset(titleoffset*1.35);
							if(iPanel==nPanels_MPnewCDF) MPhist->GetYaxis()->SetTitleSize(0.*(1-lowestBottomMargin));

							MPhist->GetYaxis()->SetLabelSize(LabelSize*1.25);
							MPhist->GetXaxis()->SetLabelSize(0.);
							MPhist->GetYaxis()->SetLabelOffset(-0.015);
							MPhist->GetXaxis()->SetLabelOffset(-0.06);

							if(iPanel==nPanels_MPnewCDF) MPhist->GetYaxis()->SetLabelSize(LabelSize*1.25*(1-lowestBottomMargin));
							MPhist->GetXaxis()->SetTitleSize(TitleSize*0.85);
							MPhist->GetXaxis()->SetAxisColor(kWhite);
							MPhist->GetYaxis()->SetAxisColor(kWhite);
							MPhist->GetXaxis()->SetTicks("-");
							MPhist->GetYaxis()->SetTicks("+");
							if(iLam==6||iLam==12||iLam==18||iLam==1||iLam==7||iLam==13) MPhist->GetYaxis()->SetNdivisions(205);

							TLegend* MPframedepLegendError;
							MPframedepLegendError=new TLegend(errorLegendX1,errorLegendY2-(errorLegendY2-errorLegendY1)*(1-lowestBottomMargin),errorLegendX2,errorLegendY2);
							MPframedepLegendError->SetFillColor(0);
							//MPframedepLegendError->SetTextFont(72);
							MPframedepLegendError->SetTextSize(errorLegendFontSize*(1-lowestBottomMargin));
							MPframedepLegendError->SetBorderSize(0);

							char MPframedepLegendEntry[200];


							TGraphAsymmErrors* graphCMS_Stat;
							TGraphAsymmErrors* graphCDF_Stat;
							TGraphAsymmErrors* graphCMS_Total;
							TGraphAsymmErrors* graphCDF_Total;

							if(iStateMP==1) graphCDF_Stat = (TGraphAsymmErrors*) infileMP1SCDF_Stat->Get(GraphName);
							if(iStateMP==2) graphCDF_Stat = (TGraphAsymmErrors*) infileMP2SCDF_Stat->Get(GraphName);
							if(iStateMP==3) graphCDF_Stat = (TGraphAsymmErrors*) infileMP3SCDF_Stat->Get(GraphName);
							if(iStateMP==1) graphCDF_Total = (TGraphAsymmErrors*) infileMP1SCDF_Total->Get(GraphName);
							if(iStateMP==2) graphCDF_Total = (TGraphAsymmErrors*) infileMP2SCDF_Total->Get(GraphName);
							if(iStateMP==3) graphCDF_Total = (TGraphAsymmErrors*) infileMP3SCDF_Total->Get(GraphName);

							if(iStateMP==1) graphCMS_Stat = (TGraphAsymmErrors*) infileMP1->Get(GraphName);
							if(iStateMP==2) graphCMS_Stat = (TGraphAsymmErrors*) infileMP2->Get(GraphName);
							if(iStateMP==3) graphCMS_Stat = (TGraphAsymmErrors*) infileMP3->Get(GraphName);
							if(iStateMP==1) graphCMS_Total = (TGraphAsymmErrors*) infileMP1_1sig->Get(GraphName);
							if(iStateMP==2) graphCMS_Total = (TGraphAsymmErrors*) infileMP2_1sig->Get(GraphName);
							if(iStateMP==3) graphCMS_Total = (TGraphAsymmErrors*) infileMP3_1sig->Get(GraphName);


							double ptCentre_CDF[nBinspT];
							double lmean_CDF[nBinspT];
							double lmean_errlow_CDF[nBinspT];
							double lmean_errhigh_CDF[nBinspT];
							double ptCentre_errlow_CDF[nBinspT];
							double ptCentre_errhigh_CDF[nBinspT];

							int nBinsCDF=graphCDF_Total->GetN();
							int ptCDF=0;
							for(int ptBinCDF=1;ptBinCDF<nBinsCDF+1;ptBinCDF++){
								graphCDF_Total->GetPoint(ptBinCDF-1,ptCentre_CDF[ptCDF],lmean_CDF[ptCDF]);
								lmean_errhigh_CDF[ptCDF]=graphCDF_Total->GetErrorYhigh(ptBinCDF-1);
								lmean_errlow_CDF[ptCDF]=graphCDF_Total->GetErrorYlow(ptBinCDF-1);
								ptCentre_errhigh_CDF[ptCDF]=graphCDF_Total->GetErrorXhigh(ptBinCDF-1);
								ptCentre_errlow_CDF[ptCDF]=graphCDF_Total->GetErrorXlow(ptBinCDF-1);
								/// Alter TGraph
								ptCentre_errhigh_CDF[ptCDF]=ColordBandWidth;
								ptCentre_errlow_CDF[ptCDF]=ColordBandWidth;

								ptCDF++;
							}
							graphCDF_Total = new TGraphAsymmErrors(nBinsCDF,ptCentre_CDF,lmean_CDF,ptCentre_errlow_CDF,ptCentre_errhigh_CDF,lmean_errlow_CDF,lmean_errhigh_CDF);

							nBinsCDF=5;
							ptCDF=0;
							for(int ptBinCDF=6;ptBinCDF<11;ptBinCDF++){
								graphCMS_Total->GetPoint(ptBinCDF-1,ptCentre_CDF[ptCDF],lmean_CDF[ptCDF]);
								lmean_errhigh_CDF[ptCDF]=graphCMS_Total->GetErrorYhigh(ptBinCDF-1);
								lmean_errlow_CDF[ptCDF]=graphCMS_Total->GetErrorYlow(ptBinCDF-1);
								ptCentre_errhigh_CDF[ptCDF]=graphCMS_Total->GetErrorXhigh(ptBinCDF-1);
								ptCentre_errlow_CDF[ptCDF]=graphCMS_Total->GetErrorXlow(ptBinCDF-1);
								/// Alter TGraph
								ptCentre_errhigh_CDF[ptCDF]=ColordBandWidth;
								ptCentre_errlow_CDF[ptCDF]=ColordBandWidth;

								ptCDF++;
							}
							graphCMS_Total = new TGraphAsymmErrors(nBinsCDF,ptCentre_CDF,lmean_CDF,ptCentre_errlow_CDF,ptCentre_errhigh_CDF,lmean_errlow_CDF,lmean_errhigh_CDF);

							TLegend* MPtildeLegend;
							MPtildeLegend=new TLegend(0.8,0.75,1.,0.95);
							MPtildeLegend->SetFillColor(0);
							//MPtildeLegend->SetTextFont(72);
							MPtildeLegend->SetTextSize(0.07);
							MPtildeLegend->SetBorderSize(0);
							char MPtildeLegendEntry[200];


							//infileMP1SCDF_Total

							int CDFFillStyle=1001;//3002

							///// DRAW HERE
							graphCDF_Total->SetLineColor(kCyan-9);
							graphCDF_Total->SetFillColor(kCyan-9);
							graphCDF_Total->SetFillStyle(CDFFillStyle);//3002
							graphCDF_Total->SetLineWidth(0.);

							graphCDF_Stat->SetMarkerColor(kRed);
							graphCDF_Stat->SetLineColor(kRed);
							graphCDF_Stat->SetMarkerStyle(25);
							graphCDF_Stat->SetMarkerSize(2.5);

							graphCMS_Stat->SetLineColor(kBlack);
							graphCMS_Stat->SetMarkerColor(kBlack);
							graphCMS_Stat->SetMarkerStyle(20);
							graphCMS_Stat->SetMarkerSize(2.5);

							graphCMS_Total->SetFillColor(kGreen-7);
							graphCMS_Total->SetLineColor(kGreen-7);
							graphCMS_Total->SetFillStyle(CDFFillStyle);

							graphCDF_Total->Draw("2");
							graphCMS_Total->Draw("2");

							graphCDF_Stat->Draw("P[]");
							graphCDF_Stat->Draw("PE");
							if(PlotAlteredPPDResults){
								graphCMS_Stat->Draw("P[]");//drawGraphStyle
								graphCMS_Stat->Draw("PE");//drawGraphStyle
							}

							TGraphAsymmErrors* graphCMS_OtherTilde;
							TGraphAsymmErrors* graphCDF_OtherTilde;
							if(iLam==12){
								if(iStateMP==1) graphCMS_OtherTilde = (TGraphAsymmErrors*) infileMP1->Get("ltilde_CS_rap1");
								if(iStateMP==2) graphCMS_OtherTilde = (TGraphAsymmErrors*) infileMP2->Get("ltilde_CS_rap1");
								if(iStateMP==3) graphCMS_OtherTilde = (TGraphAsymmErrors*) infileMP3->Get("ltilde_CS_rap1");
								if(iStateMP==1) graphCDF_OtherTilde = (TGraphAsymmErrors*) infileMP1SCDF_Stat->Get("ltilde_CS_rap1");
								if(iStateMP==2) graphCDF_OtherTilde = (TGraphAsymmErrors*) infileMP2SCDF_Stat->Get("ltilde_CS_rap1");
								if(iStateMP==3) graphCDF_OtherTilde = (TGraphAsymmErrors*) infileMP3SCDF_Stat->Get("ltilde_CS_rap1");


								double ptCentre_CDF[nBinspT];
								double lmean_CDF[nBinspT];
								double lmean_errlow_CDF[nBinspT];
								double lmean_errhigh_CDF[nBinspT];
								double ptCentre_errlow_CDF[nBinspT];
								double ptCentre_errhigh_CDF[nBinspT];

								double CDFshift=0.5;
								int nBinsCDF=graphCDF_OtherTilde->GetN();
								int ptCDF=0;
								for(int ptBinCDF=1;ptBinCDF<nBinsCDF+1;ptBinCDF++){
									graphCDF_OtherTilde->GetPoint(ptBinCDF-1,ptCentre_CDF[ptCDF],lmean_CDF[ptCDF]);
									lmean_errhigh_CDF[ptCDF]=graphCDF_OtherTilde->GetErrorYhigh(ptBinCDF-1);
									lmean_errlow_CDF[ptCDF]=graphCDF_OtherTilde->GetErrorYlow(ptBinCDF-1);
									ptCentre_errhigh_CDF[ptCDF]=graphCDF_OtherTilde->GetErrorXhigh(ptBinCDF-1);
									ptCentre_errlow_CDF[ptCDF]=graphCDF_OtherTilde->GetErrorXlow(ptBinCDF-1);
									/// Alter TGraph
									ptCentre_CDF[ptCDF]+=CDFshift;
									ptCDF++;
								}
								graphCDF_OtherTilde = new TGraphAsymmErrors(nBinsCDF,ptCentre_CDF,lmean_CDF,ptCentre_errlow_CDF,ptCentre_errhigh_CDF,lmean_errlow_CDF,lmean_errhigh_CDF);

								nBinsCDF=5;
								ptCDF=0;
								for(int ptBinCDF=6;ptBinCDF<11;ptBinCDF++){
									graphCMS_OtherTilde->GetPoint(ptBinCDF-1,ptCentre_CDF[ptCDF],lmean_CDF[ptCDF]);
									lmean_errhigh_CDF[ptCDF]=graphCMS_OtherTilde->GetErrorYhigh(ptBinCDF-1);
									lmean_errlow_CDF[ptCDF]=graphCMS_OtherTilde->GetErrorYlow(ptBinCDF-1);
									ptCentre_errhigh_CDF[ptCDF]=graphCMS_OtherTilde->GetErrorXhigh(ptBinCDF-1);
									ptCentre_errlow_CDF[ptCDF]=graphCMS_OtherTilde->GetErrorXlow(ptBinCDF-1);
									/// Alter TGraph
									ptCentre_CDF[ptCDF]-=CDFshift;
									ptCDF++;
								}
								graphCMS_OtherTilde = new TGraphAsymmErrors(nBinsCDF,ptCentre_CDF,lmean_CDF,ptCentre_errlow_CDF,ptCentre_errhigh_CDF,lmean_errlow_CDF,lmean_errhigh_CDF);

								graphCMS_OtherTilde->SetMarkerColor(kBlue);
								graphCMS_OtherTilde->SetMarkerStyle(24);
								graphCMS_OtherTilde->SetMarkerSize(2.5);

								graphCDF_OtherTilde->SetMarkerColor(kMagenta);
								graphCDF_OtherTilde->SetMarkerStyle(21);
								graphCDF_OtherTilde->SetMarkerSize(2.5);

								graphCMS_OtherTilde->Draw("PX");
								graphCDF_OtherTilde->Draw("PX");

							}

							TGaxis *axisMPY1 = new TGaxis(PlotMattpTMin,yMinMP,PlotMattpTMin,yMaxMP,yMinMP,yMaxMP,AxisDivisions,"-US");
							axisMPY1->SetTickSize(ticksize);
							if(iPanel==nPanels_MPnewCDF) axisMPY1->SetTickSize(ticksize/(1-lowestBottomMargin));
							axisMPY1->Draw("same");

							TGaxis *axisMPY2 = new TGaxis(onia::pTRange[rapBin][ptBinMax],yMinMP,onia::pTRange[rapBin][ptBinMax],yMaxMP,yMinMP,yMaxMP,AxisDivisions,"+US");
							axisMPY2->SetTickSize(ticksize);
							if(iPanel==nPanels_MPnewCDF) axisMPY2->SetTickSize(ticksize/(1-lowestBottomMargin));
							axisMPY2->Draw("same");


							double deltaTrickAxisMax_MPnew;
							double deltaTrickAxisMin_MPnew;
							deltaTrickAxisMax_MPnew=-0.001;
							deltaTrickAxisMin_MPnew=-0.001;
							if(iStateMP>1) deltaTrickAxisMin_MPnew=+0.001;
							//if(iStateMP==3) deltaTrickAxisMax_MPnew=+0.001;

							TGaxis *axisM3S1 = new TGaxis(PlotMattpTMin,yMinMP,onia::pTRange[rapBin][ptBinMax],yMinMP,PlotMattpTMin+deltaTrickAxisMin_MPnew,onia::pTRange[rapBin][ptBinMax]+deltaTrickAxisMax_MPnew,AxisDivisions,"+S");
							axisM3S1->SetTickSize(ticksize*2);
							if(iPanel==nPanels_MPnewCDF) axisM3S1->SetLabelSize(LabelSize);
							if(iPanel<nPanels_MPnewCDF) axisM3S1->SetLabelSize(0);
							axisM3S1->SetLabelOffset(labelOffsetX);
							if(iPanel==nPanels_MPnewCDF) axisM3S1->SetTickSize(ticksize/(1-lowestBottomMargin));
							axisM3S1->Draw("same");

							TGaxis *axisM3S2 = new TGaxis(PlotMattpTMin,yMaxMP,onia::pTRange[rapBin][ptBinMax],yMaxMP,PlotMattpTMin+deltaTrickAxisMin_MPnew,onia::pTRange[rapBin][ptBinMax]+deltaTrickAxisMax_MPnew,AxisDivisions,"-US");
							axisM3S2->SetTickSize(ticksize*2);
							if(iPanel==nPanels_MPnewCDF) axisM3S2->SetTickSize(ticksize/(1-lowestBottomMargin));
							axisM3S2->Draw("same");

							whereTexteInPlotX=XtitlePosition-12.;
							whereTexteInPlotY=(yMaxMP+yMinMP)/2.-(yMaxMP-yMinMP)*XtitlePositionYshift;

							if(iLam==6||iLam==12||iLam==18||iLam==1||iLam==7||iLam==13) whereTexteInPlotY=0.65;

							char axistitleMPtilde[200];
							sprintf(axistitleMPtilde,"%s",axislabel);
							if(iPanel==nPanels_MPnewCDF) YaxistitleLatexSize=YaxistitleLatexSize*(1-lowestBottomMargin);
							TLatex *MPYtitletext = new TLatex(whereTexteInPlotX,whereTexteInPlotY ,axistitleMPtilde);
							MPYtitletext->SetTextSize(YaxistitleLatexSize);
							if(iPanel==nPanels_MPnewCDF) MPYtitletext->SetTextSize(YaxistitleLatexSize*(1-lowestBottomMargin));
							MPYtitletext->SetTextColor(kBlack);
							if(iLam==3||iLam==9) MPYtitletext->PaintLatex(whereTexteInPlotX,whereTexteInPlotY, YtitleAngle, YaxistitleLatexSize*2., axislabel);
							else MPYtitletext->PaintLatex(whereTexteInPlotX,whereTexteInPlotY, YtitleAngle, YaxistitleLatexSize, axislabel);//axislabel
							MPYtitletext->Draw( "same" );




							char texTexMP[200];
							if(rapBin==1) sprintf(texTexMP,"#psi(%dS)", iStateMP);
							TLatex *textMP = new TLatex(xRapTextTilde+6,yMin+(yMax-yMin)*yRapText*0.92,texTexMP);
							textMP->SetTextSize(textSizeRap);
							if(iPanel==nPanels_MPnewCDF) textMP->SetTextSize(textSizeRap*(1-lowestBottomMargin));
							textMP->Draw( "same" );

							char abcdef[200];
							if(iStateMP==1&&iPanel==1) sprintf(abcdef,"a)");
							if(iStateMP==1&&iPanel==2) sprintf(abcdef,"b)");
							if(iStateMP==2&&iPanel==1) sprintf(abcdef,"c)");
							if(iStateMP==2&&iPanel==2) sprintf(abcdef,"d)");
							if(iStateMP==3&&iPanel==1) sprintf(abcdef,"e)");
							if(iStateMP==3&&iPanel==2) sprintf(abcdef,"f)");
							cout<<abcdef<<endl;
							TLatex *tex_abcdef = new TLatex(xabcdefText,yMin+(yMax-yMin)*yRapText*0.92,abcdef);
							tex_abcdef->SetTextSize(textSizeRap);
							if(iPanel==nPanels_MPnewCDF) tex_abcdef->SetTextSize(textSizeRap*(1-lowestBottomMargin));
							//tex_abcdef->Draw( "same" );

							if(PlotFinalData&&DrawLatexStuff){

								TLine* extreme0MP = new TLine( PlotMattpTMin, 0, onia::pTRange[rapBin][ptBinMax] ,0);
								extreme0MP->SetLineWidth( 1 );
								extreme0MP->SetLineStyle( 2 );
								extreme0MP->SetLineColor( kBlack );
								extreme0MP->Draw( "same" );

								TLine* extreme1MP = new TLine( PlotMattpTMin, 1, onia::pTRange[rapBin][ptBinMax] , 1);
								extreme1MP->SetLineWidth( 1 );
								extreme1MP->SetLineStyle( 2 );
								extreme1MP->SetLineColor( kBlack );
								//if(iLam==1||iLam==7||iLam==13) extreme1MP->Draw( "same" );

								TLine* extreme2MP = new TLine( PlotMattpTMin, -1, onia::pTRange[rapBin][ptBinMax] ,-1);
								extreme2MP->SetLineWidth( 1 );
								extreme2MP->SetLineStyle( 2 );
								extreme2MP->SetLineColor( kBlack );
								//if(iLam==1||iLam==7||iLam==13) extreme2MP->Draw( "same" );
								//if(iLam==6||iLam==12||iLam==18) extreme2MP->Draw( "same" );

							}


							if(iStateMP==2&&iPanel==1){
								cout<<"DRAW CMS preliminary Latex"<<endl;
								char text[200];
								sprintf(text,"CMS     pp      #sqrt{s} = 7 TeV     L = 4.9 fb^{-1}");
								TLatex *CentralsText1MP = new TLatex(MPlatexX-10,MPlatexYmax,text);
								CentralsText1MP->SetTextSize(CentralsFontSizeMP*TwoTOthreePanelScaleFactor*0.925);
								//CentralsText1MP->Draw( "same" );
								sprintf(text,"preliminary");
								CentralsText1MP = new TLatex(MPlatexX-10,(yMax-yMin)*0.72+yMin,text);
								CentralsText1MP->SetTextSize(CentralsFontSizeMP*TwoTOthreePanelScaleFactor*0.925);
								//CentralsText1MP->Draw( "same" );

								//sprintf(text,"L = 4.9 fb^{-1}");
								//TLatex *CentralsText2MP = new TLatex(MPlatexX,MPlatexYmax-2*MPlatexDeltaYmax,text);
								//CentralsText2MP->SetTextSize(CentralsFontSizeMP);
								//CentralsText2MP->Draw( "same" );
								//sprintf(text,"pp    #sqrt{s} = 7 TeV");
								//TLatex *CentralsText3MP = new TLatex(MPlatexX,MPlatexYmax-MPlatexDeltaYmax,text);
								//CentralsText3MP->SetTextSize(CentralsFontSizeMP);
								//CentralsText3MP->Draw( "same" );

							}

							if(iStateMP==1&&iPanel==1){

								TLegend* plotCDFLegend=new TLegend(0.165,0.8,0.5,0.95);
								plotCDFLegend->SetFillColor(0);
								//plotCDFLegend->SetTextFont(72);
								plotCDFLegend->SetTextSize(0.05);
								plotCDFLegend->SetBorderSize(0);
								char Mattlegendentry[200];

								//sprintf(Mattlegendentry,"CMS tot. uncert., 68.3%% CL");
								//plotCDFLegend->AddEntry(graphCMS_Total,Mattlegendentry,"f");
								//sprintf(Mattlegendentry,"CMS stat. uncert., 68.3%% CL");
								//plotCDFLegend->AddEntry(graphCMS_Stat,Mattlegendentry,"elp");
								//sprintf(Mattlegendentry,"CDF tot. uncert., 68.3%% CL");
								//plotCDFLegend->AddEntry(graphCDF_Total,Mattlegendentry,"f");
								//sprintf(Mattlegendentry,"CDF stat. uncert., 68.3%% CL");
								//plotCDFLegend->AddEntry(graphCDF_Stat,Mattlegendentry,"elp");

								sprintf(Mattlegendentry,"CMS");
								plotCDFLegend->AddEntry(graphCMS_Stat,Mattlegendentry,"elp");
								sprintf(Mattlegendentry,"CDF");
								plotCDFLegend->AddEntry(graphCDF_Stat,Mattlegendentry,"elp");

								plotCDFLegend->Draw("same");

								TLegend* plotCDFLegend2=new TLegend(0.8,0.8,0.95,0.95);
								plotCDFLegend2->SetFillColor(0);
								//plotCDFLegend2->SetTextFont(72);
								plotCDFLegend2->SetTextSize(0.05);
								plotCDFLegend2->SetBorderSize(0);

								sprintf(Mattlegendentry,"CMS #tilde{#lambda}^{CS}");
								plotCDFLegend2->AddEntry(graphCMS_OtherTilde,Mattlegendentry,"p");
								sprintf(Mattlegendentry,"CDF #tilde{#lambda}^{CS}");
								plotCDFLegend2->AddEntry(graphCDF_OtherTilde,Mattlegendentry,"p");

								if(iLam==12) plotCDFLegend2->Draw("same");
							}

							//if(rapBin==2&&iPanel==1){
							//	char frameMPtex[200];
							//	if(MPframe==1) sprintf(frameMPtex,"CS frame");
							//	if(MPframe==2) sprintf(frameMPtex,"HX frame");
							//	if(MPframe==3) sprintf(frameMPtex,"PX frame");
							//	char textStateFrame[200];
							//	sprintf(textStateFrame,"%s", frameMPtex);
							//	TLatex *TexStateFrame = new TLatex(MPlatexX,MPlatexYmax,textStateFrame);
							//	TexStateFrame->SetTextSize(CentralsFontSizeMP);
							//	TexStateFrame->Draw( "same" );
							//}

						}//end iStateMP loop

						MPcanvasCDF->cd();

						//whereTexteInPlotX=0.488;
						//whereTexteInPlotY=startValCoordY-deltaCoordY-1.425*labelOffsetX;
						//TLatex *M3Slabeltext = new TLatex(whereTexteInPlotX,whereTexteInPlotY ,"1");
						//M3Slabeltext->SetTextSize(XaxislabelLatexSize);
						//M3Slabeltext->SetTextColor(kBlack);
						//if(iPanel==nPanels_MPnewCDF) M3Slabeltext->Draw( "same" );

						whereTexteInPlotX=0.365;
						whereTexteInPlotY=startValCoordY-deltaCoordY+0.05;

						TLatex *MPXlabeltext = new TLatex(whereTexteInPlotX,whereTexteInPlotY ,"0");
						MPXlabeltext->SetTextSize(XaxislabelLatexSize/TwoTOthreePanelScaleFactorCDF*1.05);
						MPXlabeltext->SetTextColor(kBlack);
						if(iPanel==nPanels_MPnewCDF) MPXlabeltext->Draw( "same" );

						whereTexteInPlotX+=x_tilde;
						MPXlabeltext = new TLatex(whereTexteInPlotX,whereTexteInPlotY ,"0");
						MPXlabeltext->SetTextSize(XaxislabelLatexSize/TwoTOthreePanelScaleFactorCDF*1.05);
						MPXlabeltext->SetTextColor(kBlack);
						if(iPanel==nPanels_MPnewCDF) MPXlabeltext->Draw( "same" );


						if((iLam==1||iLam==7)){
							if(mainframe==1) sprintf(filename,"%s/FinalResultsCDF_lamthCS.pdf",FigDir);
							if(mainframe==2) sprintf(filename,"%s/FinalResultsCDF_lamthHX.pdf",FigDir);
							//MPcanvasCDF->SetRightMargin(0.05);
							//MPcanvasCDF->Modified();

							if(PlotFinalData) MPcanvasCDF->SaveAs(filename);
							MPcanvasCDF->Close();
						}
						if((iLam==2||iLam==8)){
							if(mainframe==1) sprintf(filename,"%s/FinalResultsCDF_lamphCS.pdf",FigDir);
							if(mainframe==2) sprintf(filename,"%s/FinalResultsCDF_lamphHX.pdf",FigDir);
							if(PlotFinalData) MPcanvasCDF->SaveAs(filename);
							MPcanvasCDF->Close();
						}
						if((iLam==3||iLam==9)){
							if(mainframe==1) sprintf(filename,"%s/FinalResultsCDF_lamtpCS.pdf",FigDir);
							if(mainframe==2) sprintf(filename,"%s/FinalResultsCDF_lamtpHX.pdf",FigDir);
							if(PlotFinalData) MPcanvasCDF->SaveAs(filename);
							MPcanvasCDF->Close();
						}
						if((iLam==6||iLam==12)){
							if(mainframe==1) sprintf(filename,"%s/FinalResultsCDF_lamtildeCS.pdf",FigDir);
							if(mainframe==2) sprintf(filename,"%s/FinalResultsCDF_lamtildeHX.pdf",FigDir);
							if(PlotFinalData) MPcanvasCDF->SaveAs(filename);
							MPcanvasCDF->Close();
						}

					}//end CDF MP plots NEW
				}



			} //MultiPanelPlots





			double PlotpTMin = PlotpTMinInitial, 
						 PlotpTMax = PlotpTMaxInitial;

			// Systematic uncertainties

			if(iLam==1)  sprintf(filename,"%s/Systematics_CS_lth_rap%d.pdf",FigDir,rapBin);
			if(iLam==2)  sprintf(filename,"%s/Systematics_CS_lph_rap%d.pdf",FigDir,rapBin);
			if(iLam==3)  sprintf(filename,"%s/Systematics_CS_ltp_rap%d.pdf",FigDir,rapBin);
			if(iLam==4)  sprintf(filename,"%s/Systematics_CS_lthstar_rap%d.pdf",FigDir,rapBin);
			if(iLam==5)  sprintf(filename,"%s/Systematics_CS_lphstar_rap%d.pdf",FigDir,rapBin);
			if(iLam==6)  sprintf(filename,"%s/Systematics_CS_ltilde_rap%d.pdf",FigDir,rapBin);

			if(iLam==7)  sprintf(filename,"%s/Systematics_HX_lth_rap%d.pdf",FigDir,rapBin);
			if(iLam==8)  sprintf(filename,"%s/Systematics_HX_lph_rap%d.pdf",FigDir,rapBin);
			if(iLam==9)  sprintf(filename,"%s/Systematics_HX_ltp_rap%d.pdf",FigDir,rapBin);
			if(iLam==10) sprintf(filename,"%s/Systematics_HX_lthstar_rap%d.pdf",FigDir,rapBin);
			if(iLam==11) sprintf(filename,"%s/Systematics_HX_lphstar_rap%d.pdf",FigDir,rapBin);
			if(iLam==12) sprintf(filename,"%s/Systematics_HX_ltilde_rap%d.pdf",FigDir,rapBin);

			if(iLam==13) sprintf(filename,"%s/Systematics_PX_lth_rap%d.pdf",FigDir,rapBin);
			if(iLam==14) sprintf(filename,"%s/Systematics_PX_lph_rap%d.pdf",FigDir,rapBin);
			if(iLam==15) sprintf(filename,"%s/Systematics_PX_ltp_rap%d.pdf",FigDir,rapBin);
			if(iLam==16) sprintf(filename,"%s/Systematics_PX_lthstar_rap%d.pdf",FigDir,rapBin);
			if(iLam==17) sprintf(filename,"%s/Systematics_PX_lphstar_rap%d.pdf",FigDir,rapBin);
			if(iLam==18) sprintf(filename,"%s/Systematics_PX_ltilde_rap%d.pdf",FigDir,rapBin);

			TCanvas *SystCanvas = new TCanvas("SystCanvas","SystCanvas",1000,800);

			SystCanvas->SetFillColor(kWhite);
			SystCanvas->SetGrid();
			if(PlotSystematics&&PlotSysSquare) {SystCanvas->SetGridx(0);SystCanvas->SetGridy(0);}
			SystCanvas->GetFrame()->SetFillColor(kWhite);
			SystCanvas->GetFrame()->SetBorderSize(0);
			//SystCanvas->SetRightMargin(0.05) ;

			double ParametrizedFontSize[8]={0.05,0.05,0.05,0.05,0.04,0.04,0.03,0.03};
			double LegendYmin[8]={0.8,0.75,0.7,0.65,0.65,0.65,0.6,0.6}; //{0.8,0.75,0.7,0.65,0.65,0.6,0.6,0.6};

			double LegendXmin=0.6; // for FrameworkIII: 0.5
			if(ExtendLegendInX) LegendXmin=0.25;

			TLegend* plotLegend; 
			plotLegend=new TLegend(LegendXmin,LegendYmin[nSystematics-1],0.98,0.98);
			//plotLegend=new TLegend(0.12,0.12,0.45,0.45);
			if(PlotSystematics&&PlotSysSquare) {
				plotLegend=new TLegend(LegendXmin+0.12,LegendYmin[nSystematics-1],0.98,0.98);
				if(nState==4) plotLegend=new TLegend(0.15,0.12,0.45,0.45);
			}
			plotLegend->SetFillColor(kWhite);
			//plotLegend->SetTextFont(72);
			plotLegend->SetTextSize(ParametrizedFontSize[nSystematics-1]);
			plotLegend->SetBorderSize(1);
			char legendentry[200];


			double lineWidth=3;
			sprintf(drawGraphStyle,"LX");
			//sprintf(drawGraphStyle,"PE");

			TH1F *SystHisto = new TH1F;
			SystHisto = SystCanvas->DrawFrame(PlotpTMin,yMin,PlotpTMax,yMax); //to be consistant for Psi 1S and 2S
			SystHisto->SetXTitle("#it{p}_{T} [GeV]");
			SystHisto->SetYTitle(axislabel);
			//SystHisto->GetYaxis()->SetTitleOffset(1.5);

			TGraphAsymmErrors *graphStat = new TGraphAsymmErrors(nBinspT,ptCentre_,lmean_errmean_minus,ptCentreErr_low,ptCentreErr_high,0,lmean_errmean);
			graphStat->SetFillColor(kGreen-2);
			graphStat->SetFillStyle(1001);
			graphStat->SetLineColor(kGreen-2);
			graphStat->SetLineWidth(lineWidth);
			if(!PlotAsymm) graphStat->Draw("2");
			sprintf(legendentry,"Statistics");

			if(nSystematics>7){
				TGraphAsymmErrors *graphSyst12345678 = new TGraphAsymmErrors(nBinspT,ptCentre_,SystError12345678,ptCentreErr_low,ptCentreErr_high,SystError12345678,0);
				graphSyst12345678->SetFillColor(41);
				graphSyst12345678->SetFillStyle(1001);
				graphSyst12345678->SetLineColor(41);
				graphSyst12345678->SetLineWidth(lineWidth);
				if(!PlotAsymm) graphSyst12345678->Draw("2");
				else graphSyst12345678->Draw(drawGraphStyle);
				if(!PlotAsymm) plotLegend->AddEntry(graphSyst12345678,SystID8Title,"f");
				else plotLegend->AddEntry(graphSyst12345678,SystID8Title,"l");
			}
			if(nSystematics>6){
				TGraphAsymmErrors *graphSyst1234567 = new TGraphAsymmErrors(nBinspT,ptCentre_,SystError1234567,ptCentreErr_low,ptCentreErr_high,SystError1234567,0);
				graphSyst1234567->SetFillColor(46);
				graphSyst1234567->SetFillStyle(1001);
				graphSyst1234567->SetLineColor(46);
				graphSyst1234567->SetLineWidth(lineWidth);
				if(!PlotAsymm) graphSyst1234567->Draw("2");
				else graphSyst1234567->Draw(drawGraphStyle);
				if(!PlotAsymm) plotLegend->AddEntry(graphSyst1234567,SystID7Title,"f");
				else plotLegend->AddEntry(graphSyst1234567,SystID7Title,"l");
			}
			if(nSystematics>5){
				TGraphAsymmErrors *graphSyst123456 = new TGraphAsymmErrors(nBinspT,ptCentre_,SystError123456,ptCentreErr_low,ptCentreErr_high,SystError123456,0);
				graphSyst123456->SetFillColor(kMagenta);//kYellow
				graphSyst123456->SetFillStyle(1001);
				graphSyst123456->SetLineColor(kMagenta);//kYellow
				graphSyst123456->SetLineWidth(lineWidth);
				if(!PlotAsymm) graphSyst123456->Draw("2");
				else graphSyst123456->Draw(drawGraphStyle);
				if(!PlotAsymm) plotLegend->AddEntry(graphSyst123456,SystID6Title,"f");
				else plotLegend->AddEntry(graphSyst123456,SystID6Title,"l");
			}
			if(nSystematics>4){
				TGraphAsymmErrors *graphSyst12345 = new TGraphAsymmErrors(nBinspT,ptCentre_,SystError12345,ptCentreErr_low,ptCentreErr_high,SystError12345,0);
				graphSyst12345->SetFillColor(9);//9 //IfLamTildeClosure: kOrange
				if(DeltaTildeplots) graphSyst12345->SetFillColor(kOrange);
				graphSyst12345->SetFillStyle(1001);
				graphSyst12345->SetLineColor(9);//9 //IfLamTildeClosure: kOrange
				if(DeltaTildeplots) graphSyst12345->SetLineColor(kOrange);
				graphSyst12345->SetLineWidth(lineWidth);
				if(!PlotAsymm||DeltaTildeplots) graphSyst12345->Draw("2");
				else graphSyst12345->Draw(drawGraphStyle);
				if(!PlotAsymm&&!DeltaTildeplots) plotLegend->AddEntry(graphSyst12345,SystID5Title,"f");
				if(PlotAsymm&&!DeltaTildeplots) plotLegend->AddEntry(graphSyst12345,SystID5Title,"l");

				//		graphSyst12345->Draw("2");//IfLamTildeClosure

			}
			if(nSystematics>3){
				TGraphAsymmErrors *graphSyst1234 = new TGraphAsymmErrors(nBinspT,ptCentre_,SystError1234,ptCentreErr_low,ptCentreErr_high,SystError1234,0);
				graphSyst1234->SetFillColor(kOrange);
				graphSyst1234->SetFillStyle(1001);
				graphSyst1234->SetLineColor(kOrange);
				graphSyst1234->SetLineWidth(lineWidth);
				if(!PlotAsymm||DeltaTildeplots) graphSyst1234->Draw("2");
				else graphSyst1234->Draw(drawGraphStyle);
				if(!PlotAsymm||DeltaTildeplots) plotLegend->AddEntry(graphSyst1234,SystID4Title,"f");
				else plotLegend->AddEntry(graphSyst1234,SystID4Title,"l");

				//graphSyst1234->Draw("2");//IfLamTildeClosure
				//plotLegend->AddEntry(graphSyst1234,SystID4Title,"f");//IfLamTildeClosure

			}
			if(nSystematics>2){
				TGraphAsymmErrors *graphSyst123 = new TGraphAsymmErrors(nBinspT,ptCentre_,SystError123,ptCentreErr_low,ptCentreErr_high,SystError123,0);
				graphSyst123->SetFillColor(8);
				graphSyst123->SetFillStyle(1001);
				graphSyst123->SetLineColor(8);
				graphSyst123->SetLineWidth(lineWidth);
				if(!PlotAsymm) graphSyst123->Draw("2"); //if!Pull
				else graphSyst123->Draw(drawGraphStyle); //if!Pull

				/*graphSyst123->SetMarkerStyle(20); //ifPull
					graphSyst123->SetMarkerSize(2); //ifPull
					graphSyst123->SetMarkerColor(kGreen+2); //ifPull
					graphSyst123->Draw("PX"); //ifPull
					*/
				if(!PlotAsymm) plotLegend->AddEntry(graphSyst123,SystID3Title,"f");
				else plotLegend->AddEntry(graphSyst123,SystID3Title,"l");
			}
			if(nSystematics>1){
				TGraphAsymmErrors *graphSyst12 = new TGraphAsymmErrors(nBinspT,ptCentre_,SystError12,ptCentreErr_low,ptCentreErr_high,SystError12,0);
				graphSyst12->SetFillColor(kRed);
				graphSyst12->SetFillStyle(1001);
				graphSyst12->SetLineColor(kRed);
				graphSyst12->SetLineWidth(lineWidth);
				if(!PlotAsymm) graphSyst12->Draw("2"); //if!Pull
				else graphSyst12->Draw(drawGraphStyle); //if!Pull

				if(!PlotAsymm) plotLegend->AddEntry(graphSyst12,SystID2Title,"f");
				else plotLegend->AddEntry(graphSyst12,SystID2Title,"l");
			}

			TGraphAsymmErrors *graphSyst1_ = new TGraphAsymmErrors(nBinspT,ptCentre_,SystError1,ptCentreErr_low,ptCentreErr_high,SystError1,0);
			graphSyst1_->SetFillColor(kBlue);
			graphSyst1_->SetFillStyle(1001);
			graphSyst1_->SetLineColor(kBlue);
			graphSyst1_->SetLineWidth(lineWidth);
			if(!PlotAsymm) graphSyst1_->Draw("2"); //if!Pull
			else graphSyst1_->Draw(drawGraphStyle); //if!Pull
			if(!PlotAsymm) plotLegend->AddEntry(graphSyst1_,SystID1Title,"f");
			else plotLegend->AddEntry(graphSyst1_,SystID1Title,"l");

			/*		graphSyst1_->SetMarkerStyle(20); //ifPull
						graphSyst1_->SetMarkerSize(2); //ifPull
						graphSyst1_->SetMarkerColor(kBlue); //ifPull
						graphSyst1_->Draw("PX"); //ifPull
						*/

			if(!PlotAsymm) plotLegend->AddEntry(graphStat,legendentry,"f");



			if(rapBin==1) sprintf(texTex,"      |#it{y}| < 0.6");
			if(rapBin==2) sprintf(texTex,"0.6 < |#it{y}| < 1.2");
			if(rapBin==3) sprintf(texTex,"1.2 < |#it{y}| < 1.5");
			//TLatex *Systtext = new TLatex(PlotpTMax*0.75,yMin+(yMax-yMin)*0.1,texTex);
			TLatex *Systtext = new TLatex(PlotpTMax*0.75,yMin+(yMax-yMin)*0.03,texTex);
			Systtext->SetTextSize(0.035);
			if(DrawLatexStuff) Systtext->Draw( "same" );


			if(PlotLegend) plotLegend->Draw();

			if(BGratioChi2Fits) SystCanvas->SetLogy(true);
			if(PlotSystematics) SystCanvas->SaveAs(filename);
			SystCanvas->Close();

			delete SystCanvas;


		} // rapBin


	} // iLam







	FILE *NumFile;

	bool SaveTables=true;
	if(SaveTables){

		char framerap[200];

		char NumFileName[200];
		sprintf(NumFileName,"%s/FinalNumericalResults.tex",FigDir);
		NumFile = fopen(NumFileName,"w");

		fprintf(NumFile, "\n");
		fprintf(NumFile,"\\documentclass{article}\n\\usepackage[applemac]{inputenc}\n\\usepackage{amsmath}\n\\usepackage{textcomp}\n\\pagestyle{plain}\n\\usepackage{graphicx}\n\\usepackage{multicol}\n\\usepackage{geometry}\n\\usepackage{subfigure}\n\\usepackage{booktabs}\n\\usepackage{setspace}\n\n\n\n\\begin{document}\n");

		fprintf(NumFile, "\n\n\n\n");

		double lth_tab;
		double ltherr_tab;
		double ltherr_high_tab;
		double lph_tab;
		double lpherr_tab;
		double lpherr_high_tab;
		double ltp_tab;
		double ltperr_tab;
		double ltperr_high_tab;
		double ltilde_tab;
		double ltildeerr_tab;
		double ltildeerr_high_tab;

		double ltherrTotal1_tab;
		double ltherrTotal1_high_tab;
		double lpherrTotal1_tab;
		double lpherrTotal1_high_tab;
		double ltperrTotal1_tab;
		double ltperrTotal1_high_tab;
		double ltildeerrTotal1_tab;
		double ltildeerrTotal1_high_tab;

		double ltherrTotal2_tab;
		double ltherrTotal2_high_tab;
		double lpherrTotal2_tab;
		double lpherrTotal2_high_tab;
		double ltperrTotal2_tab;
		double ltperrTotal2_high_tab;
		double ltildeerrTotal2_tab;
		double ltildeerrTotal2_high_tab;

		double ltherrTotal3_tab;
		double ltherrTotal3_high_tab;
		double lpherrTotal3_tab;
		double lpherrTotal3_high_tab;
		double ltperrTotal3_tab;
		double ltperrTotal3_high_tab;
		double ltildeerrTotal3_tab;
		double ltildeerrTotal3_high_tab;

		double pTmean_tab;

		int nTables=2;
		for(int iTab=1; iTab<nTables+1;iTab++){

			fprintf(NumFile, "\n\n\n\n");

			if(iTab==1){
				if(nState>3)
					fprintf(NumFile, "\\begin{table}[!H]\n\\centering\n \\caption{Results of polarization parameters of the $\\Psi(%dS)$ analysis}\n \\begin{tabular}{|c|cccc|}\n\\hline\n",nState-3);
				else 
					fprintf(NumFile, "\\begin{table}[!H]\n\\centering\n \\caption{Results of polarization parameters of the $\\Upsilon(%dS)$ analysis}\n \\begin{tabular}{|c|cccc|}\n\\hline\n",nState);
				fprintf(NumFile, "$p_{T}$ [GeV] & $\\lambda_{\\vartheta}$ & $\\lambda_{\\varphi}$ &  $\\lambda_{\\vartheta \\varphi}$ & $\\tilde{\\lambda}$ \\\\\n");
			}
			if(iTab==2){
				if(nState>3)
					fprintf(NumFile, "\\begin{table}[!H]\n\\centering\n \\caption{Total systematic uncertainty on the polarization parameters of the $\\Psi(%dS)$ analysis}\n \\begin{tabular}{|c|cccc|}\n\\hline\n",nState-3);
				else
					fprintf(NumFile, "\\begin{table}[!H]\n\\centering\n \\caption{Total systematic uncertainty on the polarization parameters of the $\\Upsilon(%dS)$ analysis}\n \\begin{tabular}{|c|cccc|}\n\\hline\n",nState);
				fprintf(NumFile, "$p_{T}$ [GeV] & $\\sigma^{syst}(\\lambda_{\\vartheta})$ & $\\sigma^{syst}(\\lambda_{\\varphi})$ &  $\\sigma^{syst}(\\lambda_{\\vartheta \\varphi})$ & $\\sigma^{syst}(\\tilde{\\lambda})$ \\\\\n");
			}

			for(int iFrame=1; iFrame<4; iFrame++){

				int rap=0;
				for(int rapBin = 1; rapBin < nRapBins+1; rapBin++) {

					if(iFrame==1) {sprintf(framerap,"\\hline \\multicolumn{5}{|c|}{CS frame, $%1.1f < |y| < %1.1f$}\\\\ \\hline \\rule{0pt}{4mm}\n",onia::rapForPTRange[rapBin-1],onia::rapForPTRange[rapBin]); fprintf(NumFile,framerap);}
					if(iFrame==2) {sprintf(framerap,"\\hline \\multicolumn{5}{|c|}{HX frame, $%1.1f < |y| < %1.1f$}\\\\ \\hline \\rule{0pt}{4mm}\n",onia::rapForPTRange[rapBin-1],onia::rapForPTRange[rapBin]); fprintf(NumFile,framerap);}
					if(iFrame==3) {sprintf(framerap,"\\hline \\multicolumn{5}{|c|}{PX frame, $%1.1f < |y| < %1.1f$}\\\\ \\hline \\rule{0pt}{4mm}\n",onia::rapForPTRange[rapBin-1],onia::rapForPTRange[rapBin]); fprintf(NumFile,framerap);}

					int pt=0;
					for(int ptBin = ptBinMin; ptBin < ptBinMax+1; ptBin++) {


						/*
							 syst_table[iLam-1][rapBin-1][pt][0]=SystError[pt];
							 syst_table[iLam-1][rapBin-1][pt][1]=TMath::Abs(SystError1[pt]);
							 syst_table[iLam-1][rapBin-1][pt][2]=TMath::Abs(SystError2[pt]);
							 syst_table[iLam-1][rapBin-1][pt][3]=TMath::Abs(SystError3[pt]);
							 syst_table[iLam-1][rapBin-1][pt][4]=TMath::Abs(SystError4[pt]);
							 syst_table[iLam-1][rapBin-1][pt][5]=TMath::Abs(SystError5[pt]);
							 syst_table[iLam-1][rapBin-1][pt][6]=TMath::Abs(SystError6[pt]);
							 syst_table[iLam-1][rapBin-1][pt][7]=TMath::Abs(SystError7[pt]);
							 syst_table[iLam-1][rapBin-1][pt][8]=TMath::Abs(SystError8[pt]);
							 */
						if(iTab==1){
							if(iFrame==1){
								lth_tab=val_table[1][rap][pt]			;   ltherr_tab=errLow_table[1][rap][pt];            ltherr_high_tab=errHigh_table[1][rap][pt];
								lph_tab=val_table[2][rap][pt]			;   lpherr_tab=errLow_table[2][rap][pt];            lpherr_high_tab=errHigh_table[2][rap][pt];
								ltp_tab=val_table[3][rap][pt]			;   ltperr_tab=errLow_table[3][rap][pt];            ltperr_high_tab=errHigh_table[3][rap][pt];
								ltilde_tab=val_table[6][rap][pt]		;   ltildeerr_tab=errLow_table[6][rap][pt];      	ltildeerr_high_tab=errHigh_table[6][rap][pt];
							}
							if(iFrame==2){
								lth_tab=val_table[7][rap][pt]			;   ltherr_tab=errLow_table[7][rap][pt];            ltherr_high_tab=errHigh_table[7][rap][pt];
								lph_tab=val_table[8][rap][pt]			;   lpherr_tab=errLow_table[8][rap][pt];            lpherr_high_tab=errHigh_table[8][rap][pt];
								ltp_tab=val_table[9][rap][pt]			;   ltperr_tab=errLow_table[9][rap][pt];            ltperr_high_tab=errHigh_table[9][rap][pt];
								ltilde_tab=val_table[12][rap][pt]		;   ltildeerr_tab=errLow_table[12][rap][pt];      	ltildeerr_high_tab=errHigh_table[12][rap][pt];
							}
							if(iFrame==3){
								lth_tab=val_table[13][rap][pt]			;   ltherr_tab=errLow_table[13][rap][pt];            ltherr_high_tab=errHigh_table[13][rap][pt];
								lph_tab=val_table[14][rap][pt]			;   lpherr_tab=errLow_table[14][rap][pt];            lpherr_high_tab=errHigh_table[14][rap][pt];
								ltp_tab=val_table[15][rap][pt]			;   ltperr_tab=errLow_table[15][rap][pt];            ltperr_high_tab=errHigh_table[15][rap][pt];
								ltilde_tab=val_table[18][rap][pt]		;   ltildeerr_tab=errLow_table[18][rap][pt];      	ltildeerr_high_tab=errHigh_table[18][rap][pt];
							}
							fprintf(NumFile, "%1.0f--%1.0f   &  $%1.3f _{-%1.3f}^{+%1.3f} $  & $%1.3f _{-%1.3f}^{+%1.3f}$  &  $%1.3f _{-%1.3f}^{+%1.3f}$ &  $%1.3f _{-%1.3f}^{+%1.3f}$ \\\\\n",onia::pTRange[rapBin][ptBin-1],onia::pTRange[rapBin][ptBin],lth_tab,ltherr_tab,ltherr_high_tab   ,lph_tab,lpherr_tab,lpherr_high_tab   ,ltp_tab,ltperr_tab,ltperr_high_tab   ,ltilde_tab,ltildeerr_tab,ltildeerr_high_tab );
						}

						if(iTab==2){
							if(iFrame==1){
								lth_tab=syst_table[1][rap][pt][0]		;
								lph_tab=syst_table[2][rap][pt][0]		;
								ltp_tab=syst_table[3][rap][pt][0]		;
								ltilde_tab=syst_table[6][rap][pt][0]	;
							}                                       ;
							if(iFrame==2){                          ;
								lth_tab=syst_table[7][rap][pt][0]		;
								lph_tab=syst_table[8][rap][pt][0]		;
								ltp_tab=syst_table[9][rap][pt][0]		;
								ltilde_tab=syst_table[12][rap][pt][0]	;
							}                                       ;
							if(iFrame==3){                          ;
								lth_tab=syst_table[13][rap][pt][0]		;
								lph_tab=syst_table[14][rap][pt][0]		;
								ltp_tab=syst_table[15][rap][pt][0]		;
								ltilde_tab=syst_table[18][rap][pt][0]	;
							}
							fprintf(NumFile, "%1.0f--%1.0f   &  $%1.3f  $  & $%1.3f $  &  $%1.3f $ &  $%1.3f $ \\\\\n",onia::pTRange[rapBin][ptBin-1],onia::pTRange[rapBin][ptBin],lth_tab,lph_tab,ltp_tab,ltilde_tab );
						}


						pt++;
					}



					rap++;

				}//end rapBin

			}//end iFrame



			fprintf(NumFile, "\\hline\n");
			fprintf(NumFile, "\\end{tabular}\n");
			if(iTab==1) fprintf(NumFile, "\\label{tab:results_Psi%d}\n",nState-3,nState-3);
			if(iTab==2) fprintf(NumFile, "\\label{tab:syst_Psi%d}\n",nState-3,nState-3);
			fprintf(NumFile, "\\end{table}\n");
			fprintf(NumFile, "\n");

		}//end iTab





		fprintf(NumFile, "\\end{document}");

		fclose(NumFile);



		FILE *SuppFile;


		char SuppFileName[200];
		sprintf(SuppFileName,"%s/SupplementalMaterial.txt",SuppDir);
		if(nState>3)
			sprintf(SuppFileName,"%s/SupplementalMaterial_Psi%dS.txt",SuppDir,nState-3);
		if(nState==1 || nState ==4 || nState == 5) SuppFile = fopen(SuppFileName,"w");
		else SuppFile = fopen(SuppFileName,"a");

		if(nState==1){

			fprintf(SuppFile,"                                                       MEASUREMENT OF THE Y(1S), Y(2S) AND Y(3S) POLARIZATIONS\n                                                                 IN PP COLLISIONS AT SQRT(S) = 7 TeV\n\n                                                                        SUPPLEMENTAL MATERIAL\n\n\n");

			fprintf(SuppFile,"These tables list the results of the angular anisotropy parameters Lambda-Theta, Lambda-Phi, Lambda-Theta-Phi and the frame-invariant parameter Lambda-Tilde, for the Y(nS)\nstates in the Collins-Soper (CS), helicity (HX), and perpendicular helicity (PX) frames, along with their total uncertainties (TU) (68.3%%, 95.5%%, and 99.7%% CL) and statistical\nuncertainties (SU) only (68.3%% CL) for different bins of Y transverse momentum pT and rapidity y.\n\n");

			fprintf(SuppFile,"Column headings are:\n");
			fprintf(SuppFile,"      pT-min: lower border of Y transverse momentum (GeV)\n");
			fprintf(SuppFile,"      pT-max: upper border of Y transverse momentum (GeV)\n");
			fprintf(SuppFile,"     |y|-min: lower border of Y rapidity\n");
			fprintf(SuppFile,"     |y|-max: upper border of Y rapidity\n");
			fprintf(SuppFile,"    Lambda-*: central value of the result\n");
			fprintf(SuppFile,"-T.U.68.3%%CL: negative total uncertainty, 68.3%% CL\n");
			fprintf(SuppFile,"+T.U.68.3%%CL: positive total uncertainty, 68.3%% CL\n");
			fprintf(SuppFile,"-T.U.95.5%%CL: negative total uncertainty, 95.5%% CL\n");
			fprintf(SuppFile,"+T.U.95.5%%CL: positive total uncertainty, 95.5%% CL\n");
			fprintf(SuppFile,"-T.U.99.7%%CL: negative total uncertainty, 99.7%% CL\n");
			fprintf(SuppFile,"+T.U.99.7%%CL: positive total uncertainty, 99.7%% CL\n");
			fprintf(SuppFile,"-S.U.68.3%%CL: negative statistical uncertainty, 68.3%% CL\n");
			fprintf(SuppFile,"+S.U.68.3%%CL: positive statistical uncertainty, 68.3%% CL\n\n\n");

		}

		if(nState==4 || nState==5){

			fprintf(SuppFile,"                                                        MEASUREMENT OF THE Psi(1S) and Psi(2S) POLARIZATIONS\n                                                                 IN PP COLLISIONS AT SQRT(S) = 7 TeV\n\n                                                                        SUPPLEMENTAL MATERIAL\n\n\n");

			fprintf(SuppFile,"These tables list the results of the angular anisotropy parameters Lambda-Theta, Lambda-Phi, Lambda-Theta-Phi and the frame-invariant parameter Lambda-Tilde, for the Psi(nS)\nstates in the Collins-Soper (CS), helicity (HX), and perpendicular helicity (PX) frames, along with their total uncertainties (TU) (68.3%%, 95.5%%, and 99.7%% CL) and statistical\nuncertainties (SU) only (68.3%% CL) for different bins of Psi transverse momentum pT and rapidity y.\n\n");

			fprintf(SuppFile,"Column headings are:\n");
			fprintf(SuppFile,"      pT-min: lower border of Psi transverse momentum (GeV)\n");
			fprintf(SuppFile,"      pT-max: upper border of Psi transverse momentum (GeV)\n");
			fprintf(SuppFile,"     |y|-min: lower border of Psi rapidity\n");
			fprintf(SuppFile,"     |y|-max: upper border of Psi rapidity\n");
			fprintf(SuppFile,"    Lambda-*: central value of the result\n");
			fprintf(SuppFile,"-T.U.68.3%%CL: negative total uncertainty, 68.3%% CL\n");
			fprintf(SuppFile,"+T.U.68.3%%CL: positive total uncertainty, 68.3%% CL\n");
			fprintf(SuppFile,"-T.U.95.5%%CL: negative total uncertainty, 95.5%% CL\n");
			fprintf(SuppFile,"+T.U.95.5%%CL: positive total uncertainty, 95.5%% CL\n");
			fprintf(SuppFile,"-T.U.99.7%%CL: negative total uncertainty, 99.7%% CL\n");
			fprintf(SuppFile,"+T.U.99.7%%CL: positive total uncertainty, 99.7%% CL\n");
			fprintf(SuppFile,"-S.U.68.3%%CL: negative statistical uncertainty, 68.3%% CL\n");
			fprintf(SuppFile,"+S.U.68.3%%CL: positive statistical uncertainty, 68.3%% CL\n\n\n");

		}


		double lambda_tab;
		double lambdaerr_tab;
		double lambdaerr_high_tab;

		double lambdaerrTotal1_tab;
		double lambdaerrTotal1_high_tab;

		double lambdaerrTotal2_tab;
		double lambdaerrTotal2_high_tab;

		double lambdaerrTotal3_tab;
		double lambdaerrTotal3_high_tab;

		char TabParChar[200];
		char TabParCharSpace[200];
		char TabRapChar[200];
		char TabFrameChar[200];


		for(int iLam=1;iLam<19;iLam++){
			if(iLam==4||iLam==5||iLam==10||iLam==11||iLam==16||iLam==17) continue;

			int Tabframe;
			if(iLam>0&&iLam<7) Tabframe=1;
			if(iLam>6&&iLam<13) Tabframe=2;
			if(iLam>12&&iLam<19) Tabframe=3;

			if(Tabframe==1) sprintf(TabFrameChar,"CS frame");
			if(Tabframe==2) sprintf(TabFrameChar,"HX frame");
			if(Tabframe==3) sprintf(TabFrameChar,"PX frame");

			if(iLam==1||iLam==7||iLam==13) sprintf(TabParChar,"Lambda-Theta");
			if(iLam==2||iLam==8||iLam==14) sprintf(TabParChar,"Lambda-Phi");
			if(iLam==3||iLam==9||iLam==15) sprintf(TabParChar,"Lambda-Theta-Phi");
			if(iLam==6||iLam==12||iLam==18) sprintf(TabParChar,"Lambda-Tilde");

			if(iLam==1||iLam==7||iLam==13)  sprintf(TabParCharSpace,"  Lambda-Theta  ");
			if(iLam==2||iLam==8||iLam==14)  sprintf(TabParCharSpace,"   Lambda-Phi   ");
			if(iLam==3||iLam==9||iLam==15)  sprintf(TabParCharSpace,"Lambda-Theta-Phi");
			if(iLam==6||iLam==12||iLam==18) sprintf(TabParCharSpace,"  Lambda-Tilde  ");

			int rap=0;
			for(int rapBin = 1; rapBin < nRapBins+1; rapBin++) {

				if(rapBin==1) sprintf(TabRapChar,"|y| < 0.6");
				if(rapBin==2) sprintf(TabRapChar,"0.6 < |y| < 1.2");
				if(rapBin==3) sprintf(TabRapChar,"1.2 < |y| < 1.5");

				if(rapBin==1){
					if(nState<4) // Y
						fprintf(SuppFile,"Table %s, %s, Y(%dS):\n",TabParChar, TabFrameChar,nState);
					else // Psi
						fprintf(SuppFile,"Table %s, %s, Psi(%dS):\n",TabParChar, TabFrameChar,nState-3);
				}
				if(rapBin==1)fprintf(SuppFile,"pT-min   pT-max   |y|-min   |y|-max   %s   -T.U.68.3%%CL    +T.U.68.3%%CL   -T.U.95.5%%CL   +T.U.95.5%%CL   -T.U.99.7%%CL   +T.U.99.7%%CL   -S.U.68.3%%CL   +S.U.68.3%%CL\n",TabParCharSpace);

				int pt=0;
				for(int ptBin = ptBinMin; ptBin < ptBinMax+1; ptBin++) {

					lambda_tab=val_table[iLam][rap][pt]			;   lambdaerr_tab=errLow_table[iLam][rap][pt];            lambdaerr_high_tab=errHigh_table[iLam][rap][pt];          lambdaerrTotal1_high_tab=errHighTotal1_table[iLam][rap][pt];        lambdaerrTotal1_tab=errLowTotal1_table[iLam][rap][pt];          lambdaerrTotal2_high_tab=errHighTotal2_table[iLam][rap][pt];        lambdaerrTotal2_tab=errLowTotal2_table[iLam][rap][pt];         lambdaerrTotal3_high_tab=errHighTotal3_table[iLam][rap][pt];        lambdaerrTotal3_tab=errLowTotal3_table[iLam][rap][pt];

					pTmean_tab=pTmean_table[1][rap][pt];

					if(onia::pTRange[rapBin][ptBin-1]<10)
						fprintf(SuppFile, "   %1.0f       %1.0f       %1.1f       %1.1f          % 1.3f           -%1.3f          +%1.3f         -%1.3f         +%1.3f         -%1.3f         +%1.3f         -%1.3f         +%1.3f\n",onia::pTRange[rapBin][ptBin-1],onia::pTRange[rapBin][ptBin],onia::rapForPTRange[rapBin-1],onia::rapForPTRange[rapBin],lambda_tab,lambdaerrTotal1_tab, lambdaerrTotal1_high_tab,lambdaerrTotal2_tab, lambdaerrTotal2_high_tab,lambdaerrTotal3_tab, lambdaerrTotal3_high_tab,lambdaerr_tab, lambdaerr_high_tab );
					else
						fprintf(SuppFile, "  %1.0f       %1.0f       %1.1f       %1.1f          % 1.3f           -%1.3f          +%1.3f         -%1.3f         +%1.3f         -%1.3f         +%1.3f         -%1.3f         +%1.3f\n",onia::pTRange[rapBin][ptBin-1],onia::pTRange[rapBin][ptBin],onia::rapForPTRange[rapBin-1],onia::rapForPTRange[rapBin],lambda_tab,lambdaerrTotal1_tab, lambdaerrTotal1_high_tab,lambdaerrTotal2_tab, lambdaerrTotal2_high_tab,lambdaerrTotal3_tab, lambdaerrTotal3_high_tab,lambdaerr_tab, lambdaerr_high_tab );

					pt++;


				}



				rap++;

				if(nState < 5 && rapBin==2) fprintf(SuppFile,"\n");
				if(nState ==5 && rapBin==3) fprintf(SuppFile,"\n");

			}//end rapBin

		}//end iFrame


		fclose(SuppFile);

	}


	return 0;
}
