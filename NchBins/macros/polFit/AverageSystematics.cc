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

int main(int argc, char** argv) {
	Char_t *storagedir = "Default"; //Storage Directory
	Char_t *basedir = "Default"; //Storage Directory
	Char_t *JobID1 = "Default";
	Char_t *JobID2 = "Default";
	Char_t *JobID3 = "Default";
	Char_t *JobID4 = "Default";
	Char_t *JobID5 = "Default";
	Char_t *JobID6 = "Default";
	Char_t *JobID7 = "Default";
	Char_t *JobID8 = "Default";
	Char_t *JobID9 = "Default";
	Char_t *SystID = "Default";

	int ptBinMin=1;
	int ptBinMax=1;
	int nState=1;
	int nSystematics=1;

	for( int i=0;i < argc; ++i ) {

		if(std::string(argv[i]).find("JobID1") != std::string::npos) {char* JobID1char = argv[i]; char* JobID1char2 = strtok (JobID1char, "="); JobID1 = JobID1char2; cout<<"JobID1 = "<<JobID1<<endl;}
		if(std::string(argv[i]).find("JobID2") != std::string::npos) {char* JobID2char = argv[i]; char* JobID2char2 = strtok (JobID2char, "="); JobID2 = JobID2char2; cout<<"JobID2 = "<<JobID2<<endl;}
		if(std::string(argv[i]).find("JobID3") != std::string::npos) {char* JobID3char = argv[i]; char* JobID3char2 = strtok (JobID3char, "="); JobID3 = JobID3char2; cout<<"JobID3 = "<<JobID3<<endl;}
		if(std::string(argv[i]).find("JobID4") != std::string::npos) {char* JobID4char = argv[i]; char* JobID4char2 = strtok (JobID4char, "="); JobID4 = JobID4char2; cout<<"JobID4 = "<<JobID4<<endl;}
		if(std::string(argv[i]).find("JobID5") != std::string::npos) {char* JobID5char = argv[i]; char* JobID5char2 = strtok (JobID5char, "="); JobID5 = JobID5char2; cout<<"JobID5 = "<<JobID5<<endl;}
		if(std::string(argv[i]).find("JobID6") != std::string::npos) {char* JobID6char = argv[i]; char* JobID6char2 = strtok (JobID6char, "="); JobID6 = JobID6char2; cout<<"JobID6 = "<<JobID6<<endl;}
		if(std::string(argv[i]).find("JobID7") != std::string::npos) {char* JobID7char = argv[i]; char* JobID7char2 = strtok (JobID7char, "="); JobID7 = JobID7char2; cout<<"JobID7 = "<<JobID7<<endl;}
		if(std::string(argv[i]).find("JobID8") != std::string::npos) {char* JobID8char = argv[i]; char* JobID8char2 = strtok (JobID8char, "="); JobID8 = JobID8char2; cout<<"JobID8 = "<<JobID8<<endl;}
		if(std::string(argv[i]).find("JobID9") != std::string::npos) {char* JobID9char = argv[i]; char* JobID9char2 = strtok (JobID9char, "="); JobID9 = JobID9char2; cout<<"JobID9 = "<<JobID9<<endl;}

		if(std::string(argv[i]).find("SystID") != std::string::npos) {char* SystIDchar = argv[i]; char* SystIDchar2 = strtok (SystIDchar, "="); SystID = SystIDchar2; cout<<"SystID = "<<SystID<<endl;}
		if(std::string(argv[i]).find("storagedir") != std::string::npos) {char* storagedirchar = argv[i]; char* storagedirchar2 = strtok (storagedirchar, "="); storagedir = storagedirchar2; cout<<"storagedir = "<<storagedir<<endl;}
		if(std::string(argv[i]).find("basedir") != std::string::npos) {char* basedirchar = argv[i]; char* basedirchar2 = strtok (basedirchar, "="); basedir = basedirchar2; cout<<"basedir = "<<basedir<<endl;}

		if(std::string(argv[i]).find("ptBinMin") != std::string::npos) {char* ptBinMinchar = argv[i]; char* ptBinMinchar2 = strtok (ptBinMinchar, "p"); ptBinMin = atof(ptBinMinchar2); cout<<"ptBinMin = "<<ptBinMin<<endl;}
		if(std::string(argv[i]).find("ptBinMax") != std::string::npos) {char* ptBinMaxchar = argv[i]; char* ptBinMaxchar2 = strtok (ptBinMaxchar, "p"); ptBinMax = atof(ptBinMaxchar2); cout<<"ptBinMax = "<<ptBinMax<<endl;}
		if(std::string(argv[i]).find("nState") != std::string::npos) {char* nStatechar = argv[i]; char* nStatechar2 = strtok (nStatechar, "p"); nState = atof(nStatechar2); cout<<"nState = "<<nState<<endl;}
		if(std::string(argv[i]).find("nSystematics") != std::string::npos) {char* nSystematicschar = argv[i]; char* nSystematicschar2 = strtok (nSystematicschar, "p"); nSystematics = atof(nSystematicschar2); cout<<"nSystematics = "<<nSystematics<<endl;}

	}

	char tmpfilename[200];
	sprintf(tmpfilename,"%s/macros/polFit/Systematics/%s/AverageSyst/TGraphResults_Psi%dS.root",basedir,SystID,nState-3);
	gSystem->Unlink(tmpfilename);

	char filename[200];
	sprintf(filename,"%s/macros/polFit/Systematics/%s/%s/TGraphResults_Psi%dS.root",basedir,SystID,JobID1,nState-3);
	TFile *infile1 = new TFile(filename,"READ");
	sprintf(filename,"%s/macros/polFit/Systematics/%s/%s/TGraphResults_Psi%dS.root",basedir,SystID,JobID2,nState-3);
	TFile *infile2 = new TFile(filename,"READ");
	sprintf(filename,"%s/macros/polFit/Systematics/%s/%s/TGraphResults_Psi%dS.root",basedir,SystID,JobID3,nState-3);
	TFile *infile3 = new TFile(filename,"READ");
	sprintf(filename,"%s/macros/polFit/Systematics/%s/%s/TGraphResults_Psi%dS.root",basedir,SystID,JobID4,nState-3);
	TFile *infile4 = new TFile(filename,"READ");
	sprintf(filename,"%s/macros/polFit/Systematics/%s/%s/TGraphResults_Psi%dS.root",basedir,SystID,JobID5,nState-3);
	TFile *infile5 = new TFile(filename,"READ");
	sprintf(filename,"%s/macros/polFit/Systematics/%s/%s/TGraphResults_Psi%dS.root",basedir,SystID,JobID6,nState-3);
	TFile *infile6 = new TFile(filename,"READ");
	sprintf(filename,"%s/macros/polFit/Systematics/%s/%s/TGraphResults_Psi%dS.root",basedir,SystID,JobID7,nState-3);
	TFile *infile7 = new TFile(filename,"READ");
	sprintf(filename,"%s/macros/polFit/Systematics/%s/%s/TGraphResults_Psi%dS.root",basedir,SystID,JobID8,nState-3);
	TFile *infile8 = new TFile(filename,"READ");
	sprintf(filename,"%s/macros/polFit/Systematics/%s/%s/TGraphResults_Psi%dS.root",basedir,SystID,JobID9,nState-3);
	TFile *infile9 = new TFile(filename,"READ");

	char GraphName[200];

	int nRapBins = 2;
	if(nState==5) nRapBins =  3;
	cout<<"nRapBins: "<<nRapBins<<endl;

	for(int iLam = 1; iLam<19; iLam++){

		for(int rapBin = 1; rapBin <= nRapBins; rapBin++){

			sprintf(filename,"%s/macros/polFit/Systematics/%s/AverageSyst/TGraphResults_Psi%dS.root",basedir,SystID,nState-3);
			//ifChange  sprintf(filename,"%s/macros/polFit/Systematics/%s/%s_10B/TGraphResults_%dSUps.root",basedir,SystID,JobID1,nState);
			TFile *outfile = new TFile(filename,"UPDATE");


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
			cout<<"GraphName: "<<GraphName<<endl;

			TGraphAsymmErrors* graph1 = (TGraphAsymmErrors*) infile1->Get(GraphName);
			TGraphAsymmErrors* graph2 = (TGraphAsymmErrors*) infile2->Get(GraphName);
			TGraphAsymmErrors* graph3 = (TGraphAsymmErrors*) infile3->Get(GraphName);
			TGraphAsymmErrors* graph4 = (TGraphAsymmErrors*) infile4->Get(GraphName);
			TGraphAsymmErrors* graph5 = (TGraphAsymmErrors*) infile5->Get(GraphName);
			TGraphAsymmErrors* graph6 = (TGraphAsymmErrors*) infile6->Get(GraphName);
			TGraphAsymmErrors* graph7 = (TGraphAsymmErrors*) infile7->Get(GraphName);
			TGraphAsymmErrors* graph8 = (TGraphAsymmErrors*) infile8->Get(GraphName);
			TGraphAsymmErrors* graph9 = (TGraphAsymmErrors*) infile9->Get(GraphName);

			int nBinspT=ptBinMax-ptBinMin+1;
			double ptCentre_[nBinspT];
			double ptCentreErr_low[nBinspT];
			double ptCentreErr_high[nBinspT];
			double lmean[nBinspT];
			double errlmean[nBinspT];

			double fit_lmean[nSystematics];
			double fit_errlmean[nSystematics];
			double fit_X[nSystematics];

			double ptCentre__[nSystematics][nBinspT];
			double ptCentreErr_low_[nSystematics][nBinspT];
			double ptCentreErr_high_[nSystematics][nBinspT];

			double lmean_[nSystematics][nBinspT];


			TGraphAsymmErrors* graph_;

			int pt=0;
			for(int ptBin = ptBinMin; ptBin < ptBinMax+1; ptBin++) {

				double lmean_Buffer=0;
				double ptCentre_Buffer=0;
				double ptCentreErr_low_Buffer=0;
				double ptCentreErr_high_Buffer=0;

				for(int iSyst=0;iSyst<nSystematics;iSyst++){

					if(iSyst==0)graph_=graph1;
					if(iSyst==1)graph_=graph2;
					if(iSyst==2)graph_=graph3;
					if(iSyst==3)graph_=graph4;
					if(iSyst==4)graph_=graph5;
					if(iSyst==5)graph_=graph6;
					if(iSyst==6)graph_=graph7;
					if(iSyst==7)graph_=graph8;
					if(iSyst==8)graph_=graph9;

					graph_->GetPoint(pt,ptCentre__[iSyst][pt],lmean_[iSyst][pt]);
					ptCentreErr_high_[iSyst][pt]=graph_->GetErrorXhigh(pt);
					ptCentreErr_low_[iSyst][pt]=graph_->GetErrorXlow(pt);

					fit_errlmean[iSyst]=graph_->GetErrorYhigh(pt);
					fit_lmean[iSyst]=lmean_[iSyst][pt];
					fit_X[iSyst]=0.;

				}

				double pTcentreReal=0;
				double pTcentreReallow=0;
				double pTcentreRealhigh=0;

				double lmeanHighest = 0.;
				for(int iSyst=0;iSyst<nSystematics;iSyst++){
					//lmean_Buffer=lmean_Buffer+TMath::Abs(lmean_[iSyst][pt]); //ifMean
					//lmean_Buffer=lmean_Buffer+lmean_[iSyst][pt]; //ifChange
					//lmean_[iSyst][pt] = lmean_[iSyst][pt] / 2. ;  //fracBg0_TO_default_half
					lmean_Buffer=lmean_Buffer+lmean_[iSyst][pt]*lmean_[iSyst][pt]; //ifSquare
					//if(ptBin<10 && iSyst==0) lmean_Buffer = lmean_Buffer + lmean_[iSyst][pt]*lmean_[iSyst][pt];
					//if(ptBin>=10 && iSyst==1)  lmean_Buffer = lmean_Buffer + lmean_[iSyst][pt]*lmean_[iSyst][pt];

					cout<<"lmean["<<iSyst<<"]["<<pt<<"]: "<<lmean_[iSyst][pt]<<endl;
					//if(lmean_[iSyst][pt]*lmean_[iSyst][pt]>lmeanHighest){
					//	lmeanHighest = lmean_[iSyst][pt]*lmean_[iSyst][pt] ;
					//	lmean_Buffer = lmeanHighest ;
					//	cout<<"lmeanHighest: "<<sqrt(lmeanHighest)<<endl;
					//}

					ptCentre_Buffer=ptCentre_Buffer+ptCentre__[iSyst][pt];
					ptCentreErr_low_Buffer=ptCentreErr_low_Buffer+ptCentreErr_low_[iSyst][pt];
					ptCentreErr_high_Buffer=ptCentreErr_high_Buffer+ptCentreErr_high_[iSyst][pt];

					if(iSyst==0) {//delete loop
						pTcentreReal=ptCentre__[iSyst][pt];
						pTcentreReallow=ptCentreErr_low_[iSyst][pt];
						pTcentreRealhigh=ptCentreErr_high_[iSyst][pt];
					}

				}

				ptCentre_[pt]=ptCentre_Buffer/nSystematics;
				ptCentreErr_low[pt]=ptCentreErr_low_Buffer/nSystematics;
				ptCentreErr_high[pt]=ptCentreErr_high_Buffer/nSystematics;

				ptCentre_[pt]=pTcentreReal;//delete
				ptCentreErr_low[pt]=pTcentreReallow;//delete
				ptCentreErr_high[pt]=pTcentreRealhigh;//delete

				//lmean[pt]=lmean_Buffer/nSystematics; //ifMean
				//if(pt>3) {lmean[pt]=0; std::cout << "point 4 set to 0" << std::endl;} // ifrho
				lmean[pt]=TMath::Sqrt(lmean_Buffer); //ifSquare
				//lmean[pt]=lmean_Buffer; 

				std::cout << pt << ": pT = " << ptCentre_[pt] << ", lambda = " << lmean[pt] << std::endl;

				// IF FIT THE NSYSTEMATIC VALUES INSTEAD OF USING THE MEAN:::

				//if(pt==9) {
				//      fit_X[2]=-999.;
				//}

				//TGraphAsymmErrors *fitGraph = new TGraphAsymmErrors(nSystematics,fit_X,fit_lmean,0,0,fit_errlmean,fit_errlmean);
				//TF1* fConst = new TF1("fConst","pol0",-1,1);
				//fitGraph->Fit("fConst","EFNR");
				//
				//lmean[pt]=fConst->GetParameter(0);
				//errlmean[pt]=fConst->GetParError(0);

				// END Fit

				pt++;
			}

			//////////////// Change TGraphs
			/*
				 double ptCentre_Change[nBinspT-2];
				 double ptCentreErr_low_Change[nBinspT-2];
				 double ptCentreErr_highChange[nBinspT-2];
				 double lmeanChange[nBinspT-2];

				 pt=0;
				 for(int ptBin = ptBinMin; ptBin < ptBinMax+1-2; ptBin++) {
				 ptCentre_Change[pt]=ToyMC::ptCentre[rapBin-1][pt];
				 ptCentreErr_low_Change[pt]=TMath::Abs(ToyMC::ptCentre[rapBin-1][pt]-onia::pTRange[rapBin][pt]);
				 ptCentreErr_highChange[pt]=TMath::Abs(ToyMC::ptCentre[rapBin-1][pt]-onia::pTRange[rapBin][pt+1]);

				 if(pt<5)lmeanChange[pt]=lmean_[0][pt];
				 if(pt==5)lmeanChange[pt]=(lmean_[0][5]+lmean_[0][6])/2.;
				 if(pt==6)lmeanChange[pt]=(lmean_[0][7]+lmean_[0][8])/2.;
				 if(pt==7)lmeanChange[pt]=lmean_[0][9];
				 if(pt==8)lmeanChange[pt]=lmean_[0][10];
				 if(pt==9)lmeanChange[pt]=lmean_[0][11];

				 pt++;
				 }


				 TGraphAsymmErrors *graphSyst = new TGraphAsymmErrors(nBinspT-2,ptCentre_Change,lmeanChange,ptCentreErr_low_Change,ptCentreErr_highChange,0,0);
				 graphSyst->SetMarkerColor(ToyMC::MarkerColor[rapBin]);
				 graphSyst->SetLineColor(ToyMC::MarkerColor[rapBin]);
			//              graphSyst->SetMarkerStyle(ToyMC::MarkerStyle[rapBin]);
			graphSyst->SetMarkerSize(2);
			graphSyst->SetName(GraphName);
			*/
			//////////////// END Change TGraphs

			TGraphAsymmErrors *graphSyst = new TGraphAsymmErrors(nBinspT,ptCentre_,lmean,ptCentreErr_low,ptCentreErr_high,0,0);//Original
			//TGraphAsymmErrors *graphSyst = new TGraphAsymmErrors(nBinspT,ptCentre_,lmean,ptCentreErr_low,ptCentreErr_high,errlmean,errlmean);//If Fit the nSyst values with constant
			graphSyst->SetMarkerColor(ToyMC::MarkerColor[rapBin]);
			graphSyst->SetLineColor(ToyMC::MarkerColor[rapBin]);
			//graphSyst->SetMarkerStyle(ToyMC::MarkerStyle[rapBin]);
			graphSyst->SetMarkerSize(2);
			graphSyst->SetName(GraphName);

			outfile->cd();
			graphSyst->Draw("P");
			graphSyst->Write();

			outfile->Write();
			outfile->Close();
			delete outfile;
			outfile = NULL;


		}


	}





	return 0;
}
