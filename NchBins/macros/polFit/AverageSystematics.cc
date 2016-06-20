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
	int cpmBinMin=1;
	int cpmBinMax=1;
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
		if(std::string(argv[i]).find("cpmBinMin") != std::string::npos) {char* cpmBinMinchar = argv[i]; char* cpmBinMinchar2 = strtok (cpmBinMinchar, "p"); cpmBinMin = atof(cpmBinMinchar2); cout<<"cpmBinMin = "<<cpmBinMin<<endl;}
		if(std::string(argv[i]).find("cpmBinMax") != std::string::npos) {char* cpmBinMaxchar = argv[i]; char* cpmBinMaxchar2 = strtok (cpmBinMaxchar, "p"); cpmBinMax = atof(cpmBinMaxchar2); cout<<"cpmBinMax = "<<cpmBinMax<<endl;}
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

	int nRapBins = 1;
	int nPtBins = 2;
	if(nState==5) nRapBins =  3;
	cout<<"nRapBins: "<<nRapBins<<endl;

	for(int iLam = 1; iLam<19; iLam++){

		for(int rapBin = 1; rapBin <= nRapBins; rapBin++){
		  for(int ptBin = 1; ptBin <= nPtBins; ptBin++){

			sprintf(filename,"%s/macros/polFit/Systematics/%s/AverageSyst/TGraphResults_Psi%dS.root",basedir,SystID,nState-3);
			//ifChange  sprintf(filename,"%s/macros/polFit/Systematics/%s/%s_10B/TGraphResults_%dSUps.root",basedir,SystID,JobID1,nState);
			TFile *outfile = new TFile(filename,"UPDATE");


			if(iLam==1)  sprintf(GraphName,"lth_CS_rap%d_pt%d",rapBin,ptBin);
			if(iLam==2)  sprintf(GraphName,"lph_CS_rap%d_pt%d",rapBin,ptBin);
			if(iLam==3)  sprintf(GraphName,"ltp_CS_rap%d_pt%d",rapBin,ptBin);
			if(iLam==4)  sprintf(GraphName,"lthstar_CS_rap%d_pt%d",rapBin,ptBin);
			if(iLam==5)  sprintf(GraphName,"lphstar_CS_rap%d_pt%d",rapBin,ptBin);
			if(iLam==6)  sprintf(GraphName,"ltilde_CS_rap%d_pt%d",rapBin,ptBin);

			if(iLam==7)  sprintf(GraphName,"lth_HX_rap%d_pt%d",rapBin,ptBin);
			if(iLam==8)  sprintf(GraphName,"lph_HX_rap%d_pt%d",rapBin,ptBin);
			if(iLam==9)  sprintf(GraphName,"ltp_HX_rap%d_pt%d",rapBin,ptBin);
			if(iLam==10) sprintf(GraphName,"lthstar_HX_rap%d_pt%d",rapBin,ptBin);
			if(iLam==11) sprintf(GraphName,"lphstar_HX_rap%d_pt%d",rapBin,ptBin);
			if(iLam==12) sprintf(GraphName,"ltilde_HX_rap%d_pt%d",rapBin,ptBin);

			if(iLam==13) sprintf(GraphName,"lth_PX_rap%d_pt%d",rapBin,ptBin);
			if(iLam==14) sprintf(GraphName,"lph_PX_rap%d_pt%d",rapBin,ptBin);
			if(iLam==15) sprintf(GraphName,"ltp_PX_rap%d_pt%d",rapBin,ptBin);
			if(iLam==16) sprintf(GraphName,"lthstar_PX_rap%d_pt%d",rapBin,ptBin);
			if(iLam==17) sprintf(GraphName,"lphstar_PX_rap%d_pt%d",rapBin,ptBin);
			if(iLam==18) sprintf(GraphName,"ltilde_PX_rap%d_pt%d",rapBin,ptBin);
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

			int nBinscpm=cpmBinMax-cpmBinMin+1;
			double cpmCentre_[nBinscpm];
			double cpmCentreErr_low[nBinscpm];
			double cpmCentreErr_high[nBinscpm];
			double lmean[nBinscpm];
			double errlmean[nBinscpm];

			double fit_lmean[nSystematics];
			double fit_errlmean[nSystematics];
			double fit_X[nSystematics];

			double cpmCentre__[nSystematics][nBinscpm];
			double cpmCentreErr_low_[nSystematics][nBinscpm];
			double cpmCentreErr_high_[nSystematics][nBinscpm];

			double lmean_[nSystematics][nBinscpm];


			TGraphAsymmErrors* graph_;

			int cpm=0;
			for(int cpmBin = cpmBinMin; cpmBin < cpmBinMax+1; cpmBin++) {

				double lmean_Buffer=0;
				double cpmCentre_Buffer=0;
				double cpmCentreErr_low_Buffer=0;
				double cpmCentreErr_high_Buffer=0;

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

					graph_->GetPoint(cpm,cpmCentre__[iSyst][cpm],lmean_[iSyst][cpm]);
					cpmCentreErr_high_[iSyst][cpm]=graph_->GetErrorXhigh(cpm);
					cpmCentreErr_low_[iSyst][cpm]=graph_->GetErrorXlow(cpm);

					fit_errlmean[iSyst]=graph_->GetErrorYhigh(cpm);
					fit_lmean[iSyst]=lmean_[iSyst][cpm];
					fit_X[iSyst]=0.;

				}

				double cpmcentreReal=0;
				double cpmcentreReallow=0;
				double cpmcentreRealhigh=0;

				double lmeanHighest = 0.;
				for(int iSyst=0;iSyst<nSystematics;iSyst++){
					lmean_Buffer=lmean_Buffer+TMath::Abs(lmean_[iSyst][cpm]); //ifMean
					//lmean_Buffer=lmean_Buffer+lmean_[iSyst][cpm]; //ifChange
					//lmean_[iSyst][cpm] = lmean_[iSyst][cpm] / 2. ;  //fracBg0_TO_default_half
					//lmean_Buffer=lmean_Buffer+lmean_[iSyst][cpm]*lmean_[iSyst][cpm]; //ifSquare
					//if(cpmBin<10 && iSyst==0) lmean_Buffer = lmean_Buffer + lmean_[iSyst][cpm]*lmean_[iSyst][cpm];
					//if(cpmBin>=10 && iSyst==1)  lmean_Buffer = lmean_Buffer + lmean_[iSyst][cpm]*lmean_[iSyst][cpm];

					cout<<"lmean["<<iSyst<<"]["<<cpm<<"]: "<<lmean_[iSyst][cpm]<<endl;
					//if(lmean_[iSyst][cpm]*lmean_[iSyst][cpm]>lmeanHighest){
					//	lmeanHighest = lmean_[iSyst][cpm]*lmean_[iSyst][cpm] ;
					//	lmean_Buffer = lmeanHighest ;
					//	cout<<"lmeanHighest: "<<sqrt(lmeanHighest)<<endl;
					//}

					cpmCentre_Buffer=cpmCentre_Buffer+cpmCentre__[iSyst][cpm];
					cpmCentreErr_low_Buffer=cpmCentreErr_low_Buffer+cpmCentreErr_low_[iSyst][cpm];
					cpmCentreErr_high_Buffer=cpmCentreErr_high_Buffer+cpmCentreErr_high_[iSyst][cpm];

					if(iSyst==0) {//delete loop
						cpmcentreReal=cpmCentre__[iSyst][cpm];
						cpmcentreReallow=cpmCentreErr_low_[iSyst][cpm];
						cpmcentreRealhigh=cpmCentreErr_high_[iSyst][cpm];
					}

				}

				cpmCentre_[cpm]=cpmCentre_Buffer/nSystematics;
				cpmCentreErr_low[cpm]=cpmCentreErr_low_Buffer/nSystematics;
				cpmCentreErr_high[cpm]=cpmCentreErr_high_Buffer/nSystematics;

				cpmCentre_[cpm]=cpmcentreReal;//delete
				cpmCentreErr_low[cpm]=cpmcentreReallow;//delete
				cpmCentreErr_high[cpm]=cpmcentreRealhigh;//delete

				lmean[cpm]=lmean_Buffer/nSystematics; //ifMean
				//if(cpm>3) {lmean[cpm]=0; std::cout << "point 4 set to 0" << std::endl;} // ifrho
				//lmean[cpm]=TMath::Sqrt(lmean_Buffer); //ifSquare
				//lmean[cpm]=lmean_Buffer; 

				std::cout << cpm << ": cpm = " << cpmCentre_[cpm] << ", lambda = " << lmean[cpm] << std::endl;

				// IF FIT THE NSYSTEMATIC VALUES INSTEAD OF USING THE MEAN:::

				//if(cpm==9) {
				//      fit_X[2]=-999.;
				//}

				//TGraphAsymmErrors *fitGraph = new TGraphAsymmErrors(nSystematics,fit_X,fit_lmean,0,0,fit_errlmean,fit_errlmean);
				//TF1* fConst = new TF1("fConst","pol0",-1,1);
				//fitGraph->Fit("fConst","EFNR");
				//
				//lmean[cpm]=fConst->GetParameter(0);
				//errlmean[cpm]=fConst->GetParError(0);

				// END Fit

				cpm++;
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

			TGraphAsymmErrors *graphSyst = new TGraphAsymmErrors(nBinscpm,cpmCentre_,lmean,cpmCentreErr_low,cpmCentreErr_high,0,0);//Original
			//TGraphAsymmErrors *graphSyst = new TGraphAsymmErrors(nBinscpm,cpmCentre_,lmean,cpmCentreErr_low,cpmCentreErr_high,errlmean,errlmean);//If Fit the nSyst values with constant
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


	}





	return 0;
}
