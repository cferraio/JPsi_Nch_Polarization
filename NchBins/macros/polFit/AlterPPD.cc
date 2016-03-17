/*
 * AverageSystematics.cc
 *
 *  Created on: Dec 5, 2011
 *      Author: valentinknuenz
 */

#include <iostream>
#include <string>
#include <sstream>
using namespace std;

#include "rootIncludes.inc"
#include "commonVar.h"
#include "ToyMC.h"
#include "TRandom3.h"

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
	Char_t *ShiftID = "Default";

	int ptBinMin=1;
	int ptBinMax=1;
	int rapBinMin=1;
	int rapBinMax=1;
	int nSystematics=1;
	int nState=1;
	bool ShiftResults=false;

	int SystID1ProbDist;
	int SystID2ProbDist;
	int SystID3ProbDist;
	int SystID4ProbDist;
	int SystID5ProbDist;
	int SystID6ProbDist;
	int SystID7ProbDist;
	int SystID8ProbDist;

	for( int i=0;i < argc; ++i ) {

		if(std::string(argv[i]).find("JobID") != std::string::npos) {char* JobIDchar = argv[i]; char* JobIDchar2 = strtok (JobIDchar, "="); JobID = JobIDchar2; cout<<"JobID = "<<JobID<<endl;}
		if(std::string(argv[i]).find("DefaultID") != std::string::npos) {char* DefaultIDchar = argv[i]; char* DefaultIDchar2 = strtok (DefaultIDchar, "="); DefaultID = DefaultIDchar2; cout<<"DefaultID = "<<DefaultID<<endl;}
		if(std::string(argv[i]).find("ShiftID") != std::string::npos) {char* ShiftIDchar = argv[i]; char* ShiftIDchar2 = strtok (ShiftIDchar, "="); ShiftID = ShiftIDchar2; cout<<"ShiftID = "<<ShiftID<<endl;}
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


		if(std::string(argv[i]).find("ptBinMin") != std::string::npos) {char* ptBinMinchar = argv[i]; char* ptBinMinchar2 = strtok (ptBinMinchar, "p"); ptBinMin = atof(ptBinMinchar2); cout<<"ptBinMin = "<<ptBinMin<<endl;}
		if(std::string(argv[i]).find("ptBinMax") != std::string::npos) {char* ptBinMaxchar = argv[i]; char* ptBinMaxchar2 = strtok (ptBinMaxchar, "p"); ptBinMax = atof(ptBinMaxchar2); cout<<"ptBinMax = "<<ptBinMax<<endl;}
		if(std::string(argv[i]).find("rapBinMin") != std::string::npos) {char* rapBinMinchar = argv[i]; char* rapBinMinchar2 = strtok (rapBinMinchar, "p"); rapBinMin = atof(rapBinMinchar2); cout<<"rapBinMin = "<<rapBinMin<<endl;}
		if(std::string(argv[i]).find("rapBinMax") != std::string::npos) {char* rapBinMaxchar = argv[i]; char* rapBinMaxchar2 = strtok (rapBinMaxchar, "p"); rapBinMax = atof(rapBinMaxchar2); cout<<"rapBinMax = "<<rapBinMax<<endl;}
		if(std::string(argv[i]).find("nSystematics") != std::string::npos) {char* nSystematicschar = argv[i]; char* nSystematicschar2 = strtok (nSystematicschar, "p"); nSystematics = atof(nSystematicschar2); cout<<"nSystematics = "<<nSystematics<<endl;}
		if(std::string(argv[i]).find("nState") != std::string::npos) {char* nStatechar = argv[i]; char* nStatechar2 = strtok (nStatechar, "p"); nState = atof(nStatechar2); cout<<"nState = "<<nState<<endl;}

		if(std::string(argv[i]).find("SystID1ProbDist") != std::string::npos) {char* SystID1ProbDistchar = argv[i]; char* SystID1ProbDistchar2 = strtok (SystID1ProbDistchar, "p"); SystID1ProbDist = atof(SystID1ProbDistchar2); cout<<"SystID1ProbDist = "<<SystID1ProbDist<<endl;}
		if(std::string(argv[i]).find("SystID2ProbDist") != std::string::npos) {char* SystID2ProbDistchar = argv[i]; char* SystID2ProbDistchar2 = strtok (SystID2ProbDistchar, "p"); SystID2ProbDist = atof(SystID2ProbDistchar2); cout<<"SystID2ProbDist = "<<SystID2ProbDist<<endl;}
		if(std::string(argv[i]).find("SystID3ProbDist") != std::string::npos) {char* SystID3ProbDistchar = argv[i]; char* SystID3ProbDistchar2 = strtok (SystID3ProbDistchar, "p"); SystID3ProbDist = atof(SystID3ProbDistchar2); cout<<"SystID3ProbDist = "<<SystID3ProbDist<<endl;}
		if(std::string(argv[i]).find("SystID4ProbDist") != std::string::npos) {char* SystID4ProbDistchar = argv[i]; char* SystID4ProbDistchar2 = strtok (SystID4ProbDistchar, "p"); SystID4ProbDist = atof(SystID4ProbDistchar2); cout<<"SystID4ProbDist = "<<SystID4ProbDist<<endl;}
		if(std::string(argv[i]).find("SystID5ProbDist") != std::string::npos) {char* SystID5ProbDistchar = argv[i]; char* SystID5ProbDistchar2 = strtok (SystID5ProbDistchar, "p"); SystID5ProbDist = atof(SystID5ProbDistchar2); cout<<"SystID5ProbDist = "<<SystID5ProbDist<<endl;}
		if(std::string(argv[i]).find("SystID6ProbDist") != std::string::npos) {char* SystID6ProbDistchar = argv[i]; char* SystID6ProbDistchar2 = strtok (SystID6ProbDistchar, "p"); SystID6ProbDist = atof(SystID6ProbDistchar2); cout<<"SystID6ProbDist = "<<SystID6ProbDist<<endl;}
		if(std::string(argv[i]).find("SystID7ProbDist") != std::string::npos) {char* SystID7ProbDistchar = argv[i]; char* SystID7ProbDistchar2 = strtok (SystID7ProbDistchar, "p"); SystID7ProbDist = atof(SystID7ProbDistchar2); cout<<"SystID7ProbDist = "<<SystID7ProbDist<<endl;}
		if(std::string(argv[i]).find("SystID8ProbDist") != std::string::npos) {char* SystID8ProbDistchar = argv[i]; char* SystID8ProbDistchar2 = strtok (SystID8ProbDistchar, "p"); SystID8ProbDist = atof(SystID8ProbDistchar2); cout<<"SystID8ProbDist = "<<SystID8ProbDist<<endl;}
		if(std::string(argv[i]).find("ShiftResults=1") != std::string::npos) {ShiftResults=true; cout<<"ShiftResults"<<endl;}


	}


	int ptBin=ptBinMin;
	int rapBin=rapBinMin;


	char filename[200];

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

	sprintf(filename,"%s/%s/TGraphResults_Psi%dS.root",storagedir,ShiftID,nState-3);
	TFile *infileShift = new TFile(filename,"READ");

	if(!ShiftResults) infileShift=infileSyst1;

	char GraphName[200];
	const int nFrames=3;
	const int nParameters=4;

	double BufferDouble;

	double SystVariation[nSystematics+1][nFrames+1][nParameters+1];
	double Shift[nFrames+1][nParameters+1];
	int ProbDist[nSystematics+1];

	TGraphAsymmErrors* graphSyst;
	TGraphAsymmErrors* graphShift;


	//// Get Syst-values and Probability distributions  ///////////////

	for(int iLam = 1; iLam<nParameters+1; iLam++){
		for(int iFrame = 1; iFrame<nFrames+1; iFrame++){
			for(int iSyst = 1; iSyst<nSystematics+1; iSyst++){

				if(iLam==1&&iFrame==1)  sprintf(GraphName,"lth_CS_rap%d",rapBin);
				if(iLam==2&&iFrame==1)  sprintf(GraphName,"lph_CS_rap%d",rapBin);
				if(iLam==3&&iFrame==1)  sprintf(GraphName,"ltp_CS_rap%d",rapBin);
				if(iLam==4&&iFrame==1)  sprintf(GraphName,"ltilde_CS_rap%d",rapBin);

				if(iLam==1&&iFrame==2)  sprintf(GraphName,"lth_HX_rap%d",rapBin);
				if(iLam==2&&iFrame==2)  sprintf(GraphName,"lph_HX_rap%d",rapBin);
				if(iLam==3&&iFrame==2)  sprintf(GraphName,"ltp_HX_rap%d",rapBin);
				if(iLam==4&&iFrame==2)  sprintf(GraphName,"ltilde_HX_rap%d",rapBin);

				if(iLam==1&&iFrame==3)  sprintf(GraphName,"lth_PX_rap%d",rapBin);
				if(iLam==2&&iFrame==3)  sprintf(GraphName,"lph_PX_rap%d",rapBin);
				if(iLam==3&&iFrame==3)  sprintf(GraphName,"ltp_PX_rap%d",rapBin);
				if(iLam==4&&iFrame==3)  sprintf(GraphName,"ltilde_PX_rap%d",rapBin);

				if(iSyst==1) {graphSyst = (TGraphAsymmErrors*) infileSyst1->Get(GraphName); ProbDist[iSyst]=SystID1ProbDist; }
				if(iSyst==2) {graphSyst = (TGraphAsymmErrors*) infileSyst2->Get(GraphName); ProbDist[iSyst]=SystID2ProbDist; }
				if(iSyst==3) {graphSyst = (TGraphAsymmErrors*) infileSyst3->Get(GraphName); ProbDist[iSyst]=SystID3ProbDist; }
				if(iSyst==4) {graphSyst = (TGraphAsymmErrors*) infileSyst4->Get(GraphName); ProbDist[iSyst]=SystID4ProbDist; }
				if(iSyst==5) {graphSyst = (TGraphAsymmErrors*) infileSyst5->Get(GraphName); ProbDist[iSyst]=SystID5ProbDist; }
				if(iSyst==6) {graphSyst = (TGraphAsymmErrors*) infileSyst6->Get(GraphName); ProbDist[iSyst]=SystID6ProbDist; }
				if(iSyst==7) {graphSyst = (TGraphAsymmErrors*) infileSyst7->Get(GraphName); ProbDist[iSyst]=SystID7ProbDist; }
				if(iSyst==8) {graphSyst = (TGraphAsymmErrors*) infileSyst8->Get(GraphName); ProbDist[iSyst]=SystID8ProbDist; }

				///////////////////////////////////////////
				if(nState == 4 && iSyst==3 && ptBin>9) //Rho factor, Jpsi, pT > 35: SystID2ProbDist = 2
				 	ProbDist[iSyst]=2; 
				///////////////////////////////////////////

				graphSyst->GetPoint(ptBin-1,BufferDouble,SystVariation[iSyst][iFrame][iLam]);

				SystVariation[iSyst][iFrame][iLam]=TMath::Abs(SystVariation[iSyst][iFrame][iLam]);

				if(iSyst==1) SystVariation[iSyst][iFrame][iLam]=SystVariation[iSyst][iFrame][iLam]*TMath::Sqrt(12.);

				cout<<"Systematic variation of Syst"<<iSyst<<" for iPar"<<iLam<<" in iFrame"<<iFrame<<": "<<SystVariation[iSyst][iFrame][iLam]<<", Probability Distribution: "<<ProbDist[iSyst]<<endl;


			}
			if(ShiftResults){
				graphShift = (TGraphAsymmErrors*) infileShift->Get(GraphName);
				graphShift->GetPoint(ptBin-1,BufferDouble,Shift[iFrame][iLam]);
				cout<<"Shift for iPar"<<iLam<<" in iFrame"<<iFrame<<": "<<Shift[iFrame][iLam]<<endl;
			}
			if(!ShiftResults) Shift[iFrame][iLam]=0.;
		}
	}


	//// Start variation of PPD ///////////////

	sprintf(filename,"%s/%s_%s/results_Psi%dS_rap%d_pT%d.root",storagedir,DefaultID,JobID,nState-3,rapBin,ptBin);
	TFile *results = new TFile(filename,"UPDATE");


	TTree* lambdaCS = (TTree*)results->Get("lambdaCS");

	double lth_CS;        lambdaCS->SetBranchAddress("lth",         &lth_CS        );
	double lph_CS;        lambdaCS->SetBranchAddress("lph",         &lph_CS        );
	double ltp_CS;        lambdaCS->SetBranchAddress("ltp",         &ltp_CS        );
	double lthstar_CS;    lambdaCS->SetBranchAddress("lthstar",     &lthstar_CS    );
	double lphstar_CS;    lambdaCS->SetBranchAddress("lphstar",     &lphstar_CS    );
	double ltilde_CS;     lambdaCS->SetBranchAddress("ltilde",      &ltilde_CS     );
	int    positivity_CS; lambdaCS->SetBranchAddress("positivity",  &positivity_CS );

	TTree* lambdaHX = (TTree*)results->Get("lambdaHX");

	double lth_HX;        lambdaHX->SetBranchAddress("lth",         &lth_HX        );
	double lph_HX;        lambdaHX->SetBranchAddress("lph",         &lph_HX        );
	double ltp_HX;        lambdaHX->SetBranchAddress("ltp",         &ltp_HX        );
	double lthstar_HX;    lambdaHX->SetBranchAddress("lthstar",     &lthstar_HX    );
	double lphstar_HX;    lambdaHX->SetBranchAddress("lphstar",     &lphstar_HX    );
	double ltilde_HX;     lambdaHX->SetBranchAddress("ltilde",      &ltilde_HX     );
	int    positivity_HX; lambdaHX->SetBranchAddress("positivity",  &positivity_HX );

	TTree* lambdaPX = (TTree*)results->Get("lambdaPX");

	double lth_PX;        lambdaPX->SetBranchAddress("lth",         &lth_PX        );
	double lph_PX;        lambdaPX->SetBranchAddress("lph",         &lph_PX        );
	double ltp_PX;        lambdaPX->SetBranchAddress("ltp",         &ltp_PX        );
	double lthstar_PX;    lambdaPX->SetBranchAddress("lthstar",     &lthstar_PX    );
	double lphstar_PX;    lambdaPX->SetBranchAddress("lphstar",     &lphstar_PX    );
	double ltilde_PX;     lambdaPX->SetBranchAddress("ltilde",      &ltilde_PX     );
	int    positivity_PX; lambdaPX->SetBranchAddress("positivity",  &positivity_PX );

	TTree* lambdaCS_out = (TTree*)lambdaCS->CloneTree(0);
	TTree* lambdaHX_out = (TTree*)lambdaHX->CloneTree(0);
	TTree* lambdaPX_out = (TTree*)lambdaPX->CloneTree(0);

	TRandom3 *gRandom = new TRandom3(0);

	int n_step;
	int n_entries;
	char TreeOutName[200];
	const int nRandPerEntry=5;
	double randomNumber[nRandPerEntry+1][nParameters+1];

	///// CS variations /////////

	n_entries = int( lambdaCS->GetEntries() );

	cout << endl;
	cout << "Reading distribution of CS parameters and alter entries of PPD (" << n_entries << " entries)"<< endl;
	cout << "-------------------------------------------------------------" << endl;
	cout << "Progress: ";

	n_step = n_entries/50;

	double lth_CS_Buff;
	double lph_CS_Buff;
	double ltp_CS_Buff;
	double ltilde_CS_Buff;

	for ( int i_entry = 0; i_entry < n_entries; i_entry++ ) {

		if ( i_entry%n_step == 0 ) cout << "X";

		lambdaCS->GetEvent( i_entry );

		int iFrame=1;
		lth_CS_Buff=lth_CS+Shift[iFrame][1];
		lph_CS_Buff=lph_CS+Shift[iFrame][2];
		ltp_CS_Buff=ltp_CS+Shift[iFrame][3];
		ltilde_CS_Buff=ltilde_CS+Shift[iFrame][4];

		for(int iRandPerEntry = 1; iRandPerEntry<nRandPerEntry+1; iRandPerEntry++){

			for(int iLam = 1; iLam<nParameters+1; iLam++){
				randomNumber[iRandPerEntry][iLam]=0;
				for(int iSyst = 1; iSyst<nSystematics+1; iSyst++){

					if(ProbDist[iSyst]==1)
					 	randomNumber[iRandPerEntry][iLam]+= gRandom->Gaus( 0., SystVariation[iSyst][iFrame][iLam] );
					
					if(ProbDist[iSyst]==2)
					 	randomNumber[iRandPerEntry][iLam]+= gRandom->Uniform(-SystVariation[iSyst][iFrame][iLam]/2.,  
								SystVariation[iSyst][iFrame][iLam]/2.);

				}//iSyst
			}//iLam

			lth_CS=lth_CS_Buff+randomNumber[iRandPerEntry][1];
			lph_CS=lph_CS_Buff+randomNumber[iRandPerEntry][2];
			ltp_CS=ltp_CS_Buff+randomNumber[iRandPerEntry][3];
			ltilde_CS=ltilde_CS_Buff+randomNumber[iRandPerEntry][4];

			lambdaCS_out->Fill();

		}//iRandPerEntry

	}

	//sprintf(TreeOutName,"lambdaCS_%s",JobID);
	//lambdaCS_out->SetName(TreeOutName);

	delete lambdaCS;

	cout<<endl;

	///// HX variations /////////

	n_entries = int( lambdaHX->GetEntries() );

	cout << endl;
	cout << "Reading distribution of HX parameters and alter entries of PPD (" << n_entries << " entries)"<< endl;
	cout << "-------------------------------------------------------------" << endl;
	cout << "Progress: ";

	n_step = n_entries/50;

	double lth_HX_Buff;
	double lph_HX_Buff;
	double ltp_HX_Buff;
	double ltilde_HX_Buff;

	for ( int i_entry = 0; i_entry < n_entries; i_entry++ ) {

		if ( i_entry%n_step == 0 ) cout << "X";

		lambdaHX->GetEvent( i_entry );

		int iFrame=2;
		lth_HX_Buff=lth_HX+Shift[iFrame][1];
		lph_HX_Buff=lph_HX+Shift[iFrame][2];
		ltp_HX_Buff=ltp_HX+Shift[iFrame][3];
		ltilde_HX_Buff=ltilde_HX+Shift[iFrame][4];

		for(int iRandPerEntry = 1; iRandPerEntry<nRandPerEntry+1; iRandPerEntry++){
			for(int iLam = 1; iLam<nParameters+1; iLam++){
				randomNumber[iRandPerEntry][iLam]=0;
				for(int iSyst = 1; iSyst<nSystematics+1; iSyst++){

					if(ProbDist[iSyst]==1)
						randomNumber[iRandPerEntry][iLam]+= gRandom->Gaus( 0., SystVariation[iSyst][iFrame][iLam] );
					if(ProbDist[iSyst]==2)
						randomNumber[iRandPerEntry][iLam]+= gRandom->Uniform(-SystVariation[iSyst][iFrame][iLam]/2.,  
								SystVariation[iSyst][iFrame][iLam]/2.);

				}//iSyst
			}//iLam

			lth_HX=lth_HX_Buff+randomNumber[iRandPerEntry][1];
			lph_HX=lph_HX_Buff+randomNumber[iRandPerEntry][2];
			ltp_HX=ltp_HX_Buff+randomNumber[iRandPerEntry][3];
			ltilde_HX=ltilde_HX_Buff+randomNumber[iRandPerEntry][4];

			lambdaHX_out->Fill();

		}//iRandPerEntry

	}

	//sprintf(TreeOutName,"lambdaHX");
	//lambdaHX_out->SetName(TreeOutName);

	delete lambdaHX;

	cout<<endl;

	///// PX variations /////////

	n_entries = int( lambdaPX->GetEntries() );

	cout << endl;
	cout << "Reading distribution of PX parameters and alter entries of PPD (" << n_entries << " entries)"<< endl;
	cout << "-------------------------------------------------------------" << endl;
	cout << "Progress: ";

	n_step = n_entries/50;

	double lth_PX_Buff;
	double lph_PX_Buff;
	double ltp_PX_Buff;
	double ltilde_PX_Buff;

	for ( int i_entry = 0; i_entry < n_entries; i_entry++ ) {

		if ( i_entry%n_step == 0 ) cout << "X";

		lambdaPX->GetEvent( i_entry );

		int iFrame=3;
		lth_PX_Buff=lth_PX+Shift[iFrame][1];
		lph_PX_Buff=lph_PX+Shift[iFrame][2];
		ltp_PX_Buff=ltp_PX+Shift[iFrame][3];
		ltilde_PX_Buff=ltilde_PX+Shift[iFrame][4];

		for(int iRandPerEntry = 1; iRandPerEntry<nRandPerEntry+1; iRandPerEntry++){
			for(int iLam = 1; iLam<nParameters+1; iLam++){
				randomNumber[iRandPerEntry][iLam]=0;
				for(int iSyst = 1; iSyst<nSystematics+1; iSyst++){

					if(ProbDist[iSyst]==1)
					 	randomNumber[iRandPerEntry][iLam]+= gRandom->Gaus( 0., SystVariation[iSyst][iFrame][iLam] );
					if(ProbDist[iSyst]==2)
						randomNumber[iRandPerEntry][iLam]+= gRandom->Uniform(-SystVariation[iSyst][iFrame][iLam]/2.,  
								SystVariation[iSyst][iFrame][iLam]/2.);

				}//iSyst
			}//iLam

			lth_PX=lth_PX_Buff+randomNumber[iRandPerEntry][1];
			lph_PX=lph_PX_Buff+randomNumber[iRandPerEntry][2];
			ltp_PX=ltp_PX_Buff+randomNumber[iRandPerEntry][3];
			ltilde_PX=ltilde_PX_Buff+randomNumber[iRandPerEntry][4];

			lambdaPX_out->Fill();

		}//iRandPerEntry


	}

	//sprintf(TreeOutName,"lambdaPX_%s",JobID);
	//lambdaPX_out->SetName(TreeOutName);

	delete lambdaPX;

	cout<<endl;

	// Save Results File ///


	lambdaCS_out->Write("lambdaCS",6);
	lambdaHX_out->Write("lambdaHX",6);
	lambdaPX_out->Write("lambdaPX",6);

	delete lambdaCS_out;
	delete lambdaHX_out;
	delete lambdaPX_out;

	results->Write();
	results->Close();

	return 0;
}
