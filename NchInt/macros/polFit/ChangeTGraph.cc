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
	Char_t *SystID = "Default";
	Char_t *JobID1 = "Default";

	int ptBinMin=1;
	int ptBinMax=1;
	int nState=1;

	for( int i=0;i < argc; ++i ) {

		if(std::string(argv[i]).find("JobID1") != std::string::npos) {char* JobID1char = argv[i]; char* JobID1char2 = strtok (JobID1char, "="); JobID1 = JobID1char2; cout<<"JobID1 = "<<JobID1<<endl;}

		if(std::string(argv[i]).find("SystID") != std::string::npos) {char* SystIDchar = argv[i]; char* SystIDchar2 = strtok (SystIDchar, "="); SystID = SystIDchar2; cout<<"SystID = "<<SystID<<endl;}
		if(std::string(argv[i]).find("storagedir") != std::string::npos) {char* storagedirchar = argv[i]; char* storagedirchar2 = strtok (storagedirchar, "="); storagedir = storagedirchar2; cout<<"storagedir = "<<storagedir<<endl;}
		if(std::string(argv[i]).find("basedir") != std::string::npos) {char* basedirchar = argv[i]; char* basedirchar2 = strtok (basedirchar, "="); basedir = basedirchar2; cout<<"basedir = "<<basedir<<endl;}

		if(std::string(argv[i]).find("ptBinMin") != std::string::npos) {char* ptBinMinchar = argv[i]; char* ptBinMinchar2 = strtok (ptBinMinchar, "p"); ptBinMin = atof(ptBinMinchar2); cout<<"ptBinMin = "<<ptBinMin<<endl;}
		if(std::string(argv[i]).find("ptBinMax") != std::string::npos) {char* ptBinMaxchar = argv[i]; char* ptBinMaxchar2 = strtok (ptBinMaxchar, "p"); ptBinMax = atof(ptBinMaxchar2); cout<<"ptBinMax = "<<ptBinMax<<endl;}
		if(std::string(argv[i]).find("nState") != std::string::npos) {char* nStatechar = argv[i]; char* nStatechar2 = strtok (nStatechar, "p"); nState = atof(nStatechar2); cout<<"nState = "<<nState<<endl;}

	}

	char tmpfilename[200];
	sprintf(tmpfilename,"%s/macros/polFit/Systematics/%s/ChangedTGraph/TGraphResults_Psi%dS.root",basedir,SystID,nState-3);
	gSystem->Unlink(tmpfilename);

	char filename[200];
	sprintf(filename,"%s/TGraphResults_Psi%dS.root",JobID1,nState-3);
	TFile *infile1 = new TFile(filename,"READ");

	char GraphName[200];

	for(int iFrame = 1; iFrame<4; iFrame++){

		int nRapBins = 2;
		if(nState==5) nRapBins = 3;

		for(int rapBin = 1; rapBin < nRapBins+1; rapBin++){

			cout << "rapBin: " << rapBin << endl;
			TGraphAsymmErrors* graph_lth;
			TGraphAsymmErrors* graph_lph;
			TGraphAsymmErrors* graph_ltp;
			TGraphAsymmErrors* graph_lthstar;
			TGraphAsymmErrors* graph_lphstar;
			TGraphAsymmErrors* graph_ltilde;

			TGraphAsymmErrors* graph_DltildeCS;
			TGraphAsymmErrors* graph_DltildeHX;
			TGraphAsymmErrors* graph_DltildePX;

			if(iFrame==1)  {
				sprintf(GraphName,"lth_CS_rap%d",rapBin);
				graph_lth = (TGraphAsymmErrors*) infile1->Get(GraphName);
				sprintf(GraphName,"lph_CS_rap%d",rapBin);
				graph_lph = (TGraphAsymmErrors*) infile1->Get(GraphName);
				sprintf(GraphName,"ltp_CS_rap%d",rapBin);
				graph_ltp = (TGraphAsymmErrors*) infile1->Get(GraphName);
				sprintf(GraphName,"lthstar_CS_rap%d",rapBin);
				graph_lthstar = (TGraphAsymmErrors*) infile1->Get(GraphName);
				sprintf(GraphName,"lphstar_CS_rap%d",rapBin);
				graph_lphstar = (TGraphAsymmErrors*) infile1->Get(GraphName);
				sprintf(GraphName,"ltilde_CS_rap%d",rapBin);
				graph_ltilde = (TGraphAsymmErrors*) infile1->Get(GraphName);

				sprintf(GraphName,"ltilde_CS_rap%d",rapBin);
				graph_DltildeCS = (TGraphAsymmErrors*) infile1->Get(GraphName);
				sprintf(GraphName,"ltilde_HX_rap%d",rapBin);
				graph_DltildeHX = (TGraphAsymmErrors*) infile1->Get(GraphName);
				sprintf(GraphName,"ltilde_PX_rap%d",rapBin);
				graph_DltildePX = (TGraphAsymmErrors*) infile1->Get(GraphName);
			}

			if(iFrame==2)  {
				sprintf(GraphName,"lth_HX_rap%d",rapBin);
				graph_lth = (TGraphAsymmErrors*) infile1->Get(GraphName);
				sprintf(GraphName,"lph_HX_rap%d",rapBin);
				graph_lph = (TGraphAsymmErrors*) infile1->Get(GraphName);
				sprintf(GraphName,"ltp_HX_rap%d",rapBin);
				graph_ltp = (TGraphAsymmErrors*) infile1->Get(GraphName);
				sprintf(GraphName,"lthstar_HX_rap%d",rapBin);
				graph_lthstar = (TGraphAsymmErrors*) infile1->Get(GraphName);
				sprintf(GraphName,"lphstar_HX_rap%d",rapBin);
				graph_lphstar = (TGraphAsymmErrors*) infile1->Get(GraphName);
				sprintf(GraphName,"ltilde_HX_rap%d",rapBin);
				graph_ltilde = (TGraphAsymmErrors*) infile1->Get(GraphName);

				sprintf(GraphName,"ltilde_CS_rap%d",rapBin);
				graph_DltildeCS = (TGraphAsymmErrors*) infile1->Get(GraphName);
				sprintf(GraphName,"ltilde_HX_rap%d",rapBin);
				graph_DltildeHX = (TGraphAsymmErrors*) infile1->Get(GraphName);
				sprintf(GraphName,"ltilde_PX_rap%d",rapBin);
				graph_DltildePX = (TGraphAsymmErrors*) infile1->Get(GraphName);
			}

			if(iFrame==3)  {
				sprintf(GraphName,"lth_PX_rap%d",rapBin);
				graph_lth = (TGraphAsymmErrors*) infile1->Get(GraphName);
				sprintf(GraphName,"lph_PX_rap%d",rapBin);
				graph_lph = (TGraphAsymmErrors*) infile1->Get(GraphName);
				sprintf(GraphName,"ltp_PX_rap%d",rapBin);
				graph_ltp = (TGraphAsymmErrors*) infile1->Get(GraphName);
				sprintf(GraphName,"lthstar_PX_rap%d",rapBin);
				graph_lthstar = (TGraphAsymmErrors*) infile1->Get(GraphName);
				sprintf(GraphName,"lphstar_PX_rap%d",rapBin);
				graph_lphstar = (TGraphAsymmErrors*) infile1->Get(GraphName);
				sprintf(GraphName,"ltilde_PX_rap%d",rapBin);
				graph_ltilde = (TGraphAsymmErrors*) infile1->Get(GraphName);

				sprintf(GraphName,"ltilde_CS_rap%d",rapBin);
				graph_DltildeCS = (TGraphAsymmErrors*) infile1->Get(GraphName);
				sprintf(GraphName,"ltilde_HX_rap%d",rapBin);
				graph_DltildeHX = (TGraphAsymmErrors*) infile1->Get(GraphName);
				sprintf(GraphName,"ltilde_PX_rap%d",rapBin);
				graph_DltildePX = (TGraphAsymmErrors*) infile1->Get(GraphName);
			}

			cout<<"TGraphs of all parameters loaded for frame "<<iFrame<<endl;



			/*		//Ups1:
						double lamtildeCS[2][10]={
						{-0.268025, -0.423266, -0.804402, -0.775923, -0.374039, -0.261726, -0.407957, -0.244021, -0.210692, 0.010428},{ -0.72549, -0.232061, -0.422995, -0.507062, -0.641677, -0.285836, -0.233993, -0.112918, -0.0927753, 0.11018}
						};

						double lamtildeHX[2][10]={
						{0.324909, 0.309998, -0.645745, -0.628765, -0.303808, -0.259237, -0.250806, -0.181788, -0.150838, 0.0462099},{-0.693424, -0.114897, -0.365809, -0.470233, -0.593041, -0.262521, -0.192101, -0.114675, -0.0880647, 0.128067}
						};

						double lamtildePX[2][10]={
						{ 0.253969, -0.204504, -0.68701, -0.696149, -0.299076, -0.245953, -0.396063, -0.23093, -0.202291, 0.03004},{-0.684952, -0.209087, -0.400518, -0.484155, -0.612468, -0.276295, -0.227099, -0.104696, -0.0858445, 0.129378}
						};
						*/
			/*		//Ups2:
						double lamtildeCS[2][10]={
						{0.579805, 0.117124, -0.845207, -0.676587, -0.177738, -0.0878507, -0.139639, -0.0322774, -0.160865, 0.111783},{-0.36563, -0.299825, -0.269512, -0.263333, -0.582337, -0.110133, -0.229866, -0.102893, 0.284703, 0.267933}
						};

						double lamtildeHX[2][10]={
						{0.928091, 0.120682, -0.272608, -0.641832, -0.0902441, -0.147928, -0.111298, 0.0997691, -0.115126, 0.149243},{ -0.122689, -0.153951, -0.196745, -0.190463, -0.467136, -0.0634411, -0.211757, -0.0884084, 0.298757, 0.293643}
						};

						double lamtildePX[2][10]={
						{1.24984, 0.40587, -0.645656, -0.540806, -0.0248343, -0.0360626, -0.11936, -0.0110759, -0.144741, 0.14429},{ -0.336529, -0.267787, -0.229459, -0.219332, -0.518478, -0.0854088, -0.210192, -0.0806282, 0.299751, 0.304441}
						};
						*/
			//Ups3:
			double lamtildeCS[2][10]={
				{ -1.48319, -1.16822, 0.229813, -0.491375, 0.00964503, -0.178115, -0.100937, 0.00287869, -0.00585711, 0.152644},{-0.193191, 0.0248918, 0.0133266, -0.377318, 0.0195074, -0.111712, 0.0765794, 0.369197, -0.118172, 0.730991}
			};

			double lamtildeHX[2][10]={
				{-0.730572, -0.417162, 0.340792, -0.41835, 0.339563, -0.0582737, -0.0689351, 0.0448716, 0.0882122, 0.199121},{0.018715, 0.197731, 0.14685, -0.297466, 0.095253, -0.101151, 0.132685, 0.387782, -0.0860426, 0.798841}
			};

			double lamtildePX[2][10]={
				{0.038359, -0.495607, 0.464968, -0.317869, 0.240574, -0.100853, -0.0657933, 0.0389629, 0.0162931, 0.201466},{ -0.202751, 0.0605227, 0.0659666, -0.320916, 0.106622, -0.0713573, 0.0980278, 0.402986, -0.0910948, 0.806114}
			};


			int nBinspT=ptBinMax-ptBinMin+1;
			double ptCentre[nBinspT];
			double ptCentreErr_low[nBinspT];
			double ptCentreErr_high[nBinspT];

			double lth_lmean[nBinspT];
			double lth_lmeanErr_low[nBinspT];
			double lth_lmeanErr_high[nBinspT];
			double lth_lmeanErr[nBinspT];

			double lph_lmean[nBinspT];
			double lph_lmeanErr_low[nBinspT];
			double lph_lmeanErr_high[nBinspT];
			double lph_lmeanErr[nBinspT];

			double ltp_lmean[nBinspT];
			double ltp_lmeanErr_low[nBinspT];
			double ltp_lmeanErr_high[nBinspT];
			double ltp_lmeanErr[nBinspT];

			double lthstar_lmean[nBinspT];
			double lthstar_lmeanErr_low[nBinspT];
			double lthstar_lmeanErr_high[nBinspT];
			double lthstar_lmeanErr[nBinspT];

			double lphstar_lmean[nBinspT];
			double lphstar_lmeanErr_low[nBinspT];
			double lphstar_lmeanErr_high[nBinspT];
			double lphstar_lmeanErr[nBinspT];

			double ltilde_lmean[nBinspT];
			double ltilde_lmeanErr_low[nBinspT];
			double ltilde_lmeanErr_high[nBinspT];
			double ltilde_lmeanErr[nBinspT];

			double DltildeCS_lmean[nBinspT];
			double DltildeHX_lmean[nBinspT];
			double DltildePX_lmean[nBinspT];

			int pt=0;
			for(int ptBin = ptBinMin; ptBin < ptBinMax+1; ptBin++) {

				graph_lth->GetPoint(pt,ptCentre[pt],lth_lmean[pt]);
				ptCentreErr_high[pt]=graph_lth->GetErrorXhigh(pt);
				ptCentreErr_low[pt]=graph_lth->GetErrorXlow(pt);

				graph_lth->GetPoint(pt,ptCentre[pt],lth_lmean[pt]);
				lth_lmeanErr_high[pt]=graph_lth->GetErrorYhigh(pt);
				lth_lmeanErr_low[pt]=graph_lth->GetErrorYlow(pt);
				lth_lmeanErr[pt]=(lth_lmeanErr_high[pt]+lth_lmeanErr_low[pt])/2.;

				graph_lph->GetPoint(pt,ptCentre[pt],lph_lmean[pt]);
				lph_lmeanErr_high[pt]=graph_lph->GetErrorYhigh(pt);
				lph_lmeanErr_low[pt]=graph_lph->GetErrorYlow(pt);
				lph_lmeanErr[pt]=(lph_lmeanErr_high[pt]+lph_lmeanErr_low[pt])/2.;

				graph_ltp->GetPoint(pt,ptCentre[pt],ltp_lmean[pt]);
				ltp_lmeanErr_high[pt]=graph_ltp->GetErrorYhigh(pt);
				ltp_lmeanErr_low[pt]=graph_ltp->GetErrorYlow(pt);
				ltp_lmeanErr[pt]=(ltp_lmeanErr_high[pt]+ltp_lmeanErr_low[pt])/2.;

				graph_lthstar->GetPoint(pt,ptCentre[pt],lthstar_lmean[pt]);
				lthstar_lmeanErr_high[pt]=graph_lthstar->GetErrorYhigh(pt);
				lthstar_lmeanErr_low[pt]=graph_lthstar->GetErrorYlow(pt);
				lthstar_lmeanErr[pt]=(lthstar_lmeanErr_high[pt]+lthstar_lmeanErr_low[pt])/2.;

				graph_lphstar->GetPoint(pt,ptCentre[pt],lphstar_lmean[pt]);
				lphstar_lmeanErr_high[pt]=graph_lphstar->GetErrorYhigh(pt);
				lphstar_lmeanErr_low[pt]=graph_lphstar->GetErrorYlow(pt);
				lphstar_lmeanErr[pt]=(lphstar_lmeanErr_high[pt]+lphstar_lmeanErr_low[pt])/2.;

				graph_ltilde->GetPoint(pt,ptCentre[pt],ltilde_lmean[pt]);
				ltilde_lmeanErr_high[pt]=graph_ltilde->GetErrorYhigh(pt);
				ltilde_lmeanErr_low[pt]=graph_ltilde->GetErrorYlow(pt);
				ltilde_lmeanErr[pt]=(ltilde_lmeanErr_high[pt]+ltilde_lmeanErr_low[pt])/2.;

				graph_DltildeCS->GetPoint(pt,ptCentre[pt],DltildeCS_lmean[pt]);
				graph_DltildeHX->GetPoint(pt,ptCentre[pt],DltildeHX_lmean[pt]);
				graph_DltildePX->GetPoint(pt,ptCentre[pt],DltildePX_lmean[pt]);

				//			cout<<", "<<ltilde_lmean[pt];
				pt++;
			}

			cout<<"Values of all parameters loaded for frame "<<iFrame<<endl;

			for(int iParam=1;iParam<7;iParam++){

				if(iParam==1&&iFrame==1) sprintf(GraphName,"lth_CS_rap%d",rapBin);
				if(iParam==2&&iFrame==1) sprintf(GraphName,"lph_CS_rap%d",rapBin);
				if(iParam==3&&iFrame==1) sprintf(GraphName,"ltp_CS_rap%d",rapBin);
				if(iParam==4&&iFrame==1) sprintf(GraphName,"lthstar_CS_rap%d",rapBin);
				if(iParam==5&&iFrame==1) sprintf(GraphName,"lphstar_CS_rap%d",rapBin);
				if(iParam==6&&iFrame==1) sprintf(GraphName,"ltilde_CS_rap%d",rapBin);
				if(iParam==1&&iFrame==2) sprintf(GraphName,"lth_HX_rap%d",rapBin);
				if(iParam==2&&iFrame==2) sprintf(GraphName,"lph_HX_rap%d",rapBin);
				if(iParam==3&&iFrame==2) sprintf(GraphName,"ltp_HX_rap%d",rapBin);
				if(iParam==4&&iFrame==2) sprintf(GraphName,"lthstar_HX_rap%d",rapBin);
				if(iParam==5&&iFrame==2) sprintf(GraphName,"lphstar_HX_rap%d",rapBin);
				if(iParam==6&&iFrame==2) sprintf(GraphName,"ltilde_HX_rap%d",rapBin);
				if(iParam==1&&iFrame==3) sprintf(GraphName,"lth_PX_rap%d",rapBin);
				if(iParam==2&&iFrame==3) sprintf(GraphName,"lph_PX_rap%d",rapBin);
				if(iParam==3&&iFrame==3) sprintf(GraphName,"ltp_PX_rap%d",rapBin);
				if(iParam==4&&iFrame==3) sprintf(GraphName,"lthstar_PX_rap%d",rapBin);
				if(iParam==5&&iFrame==3) sprintf(GraphName,"lphstar_PX_rap%d",rapBin);
				if(iParam==6&&iFrame==3) sprintf(GraphName,"ltilde_PX_rap%d",rapBin);


				double lmean[nBinspT];
				double lmeanErr_low[nBinspT];
				double lmeanErr_high[nBinspT];

				pt=0;
				for(int ptBin = ptBinMin; ptBin < ptBinMax+1; ptBin++) {

					if(iParam==1){
						lmean[pt]=lth_lmean[pt]/2.;//0.05;//1./(1-1./3.*ToyMC::fracSignal[rapBin-1][ptBin-1]*(3+lth_lmean[pt]))*lth_lmean[pt];
						lmeanErr_low[pt]=0.;//(1.-ToyMC::fracSignal[rapBin-1][ptBin-1])/TMath::Power((1-1./3.*ToyMC::fracSignal[rapBin-1][ptBin-1]*(3+lth_lmean[pt])),2)*lth_lmeanErr[pt];
						lmeanErr_high[pt]=0.;//lmeanErr_low[pt];
					}

					if(iParam==2){
						lmean[pt]=lph_lmean[pt]/2.;//0.01;//1./(1-1./3.*ToyMC::fracSignal[rapBin-1][ptBin-1]*(3+lth_lmean[pt]))*lph_lmean[pt];
						lmeanErr_low[pt]=0.;//TMath::Sqrt( (1.)/TMath::Power((1-1./3.*ToyMC::fracSignal[rapBin-1][ptBin-1]*(3+lth_lmean[pt])),2)*TMath::Power(lph_lmeanErr[pt],2)        +       (1./9.)*(TMath::Power(ToyMC::fracSignal[rapBin-1][ptBin-1],2)*TMath::Power(lph_lmean[pt],2))/TMath::Power((1-1./3.*ToyMC::fracSignal[rapBin-1][ptBin-1]*(3+lth_lmean[pt])),4)*TMath::Power(lth_lmeanErr[pt],2));
						lmeanErr_high[pt]=0.;//lmeanErr_low[pt];
					}

					if(iParam==3){
						lmean[pt]=ltp_lmean[pt]/2.;//0.01;//1./(1-1./3.*ToyMC::fracSignal[rapBin-1][ptBin-1]*(3+lth_lmean[pt]))*ltp_lmean[pt];
						lmeanErr_low[pt]=0.;//TMath::Sqrt( (1.)/TMath::Power((1-1./3.*ToyMC::fracSignal[rapBin-1][ptBin-1]*(3+lth_lmean[pt])),2)*TMath::Power(ltp_lmeanErr[pt],2)        +       (1./9.)*(TMath::Power(ToyMC::fracSignal[rapBin-1][ptBin-1],2)*TMath::Power(ltp_lmean[pt],2))/TMath::Power((1-1./3.*ToyMC::fracSignal[rapBin-1][ptBin-1]*(3+lth_lmean[pt])),4)*TMath::Power(lth_lmeanErr[pt],2));
						lmeanErr_high[pt]=0.;//lmeanErr_low[pt];
					}

					if(iParam==4){
						lmean[pt]=0;
						lmeanErr_low[pt]=0;
						lmeanErr_high[pt]=0;
					}

					if(iParam==5){
						lmean[pt]=0;
						lmeanErr_low[pt]=0;
						lmeanErr_high[pt]=0;
					}

					if(iParam==6){
						lmean[pt]=ltilde_lmean[pt]/2.;//0.05;//0.05;
						lmeanErr_low[pt]=0.;
						lmeanErr_high[pt]=0.;

						//lmean[pt]=lamtildeCS[rapBin-1][pt]-lamtildeHX[rapBin-1][pt];
						//lmean[pt]=lamtildeHX[rapBin-1][pt]-lamtildePX[rapBin-1][pt];
						//lmean[pt]=lamtildePX[rapBin-1][pt]-lamtildeCS[rapBin-1][pt];

						lmean[pt]=DltildeCS_lmean[pt]-DltildeHX_lmean[pt];
						//lmean[pt]=DltildeHX_lmean[pt]-DltildePX_lmean[pt];
						//lmean[pt]=DltildePX_lmean[pt]-DltildeCS_lmean[pt];


					}


					pt++;
				}

				cout<<"Values of parameter "<<iParam<<" changed for frame "<<iFrame<<endl;

				TGraphAsymmErrors *ChangedGraph = new TGraphAsymmErrors(nBinspT,ptCentre,lmean,ptCentreErr_low,ptCentreErr_high,lmeanErr_low,lmeanErr_high);
				ChangedGraph->SetMarkerColor(ToyMC::MarkerColor[rapBin]);
				ChangedGraph->SetLineColor(ToyMC::MarkerColor[rapBin]);
				ChangedGraph->SetMarkerSize(2);
				ChangedGraph->SetName(GraphName);

				sprintf(filename,"%s/macros/polFit/Systematics/%s/ChangedTGraph/TGraphResults_Psi%dS.root",basedir,SystID,nState-3);
				TFile *outfile = new TFile(filename,"UPDATE");

				outfile->cd();
				ChangedGraph->Draw("P");
				ChangedGraph->Write();

				outfile->Write();
				outfile->Close();
				delete outfile;
				outfile = NULL;


			}

			cout<<"Values of all parameters changed for frame "<<iFrame<<endl;

		} //rapBin

		cout<<"Switching Rapidity"<<endl;

	} //iFrame
	cout<<"Switching frame"<<endl;



	return 0;
}
