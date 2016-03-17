#include <iostream>
#include <sstream>
#include <iomanip>


//ROOT includes
#include "TH2F.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TDirectory.h"
#include "TSystem.h"
#include "TFormula.h"
#include "TMath.h"
#include "TROOT.h"
#include "TLatex.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TAxis.h"
#include "TPad.h"
#include "TRandom3.h"
#include "TF1.h"
#include "TGraphAsymmErrors.h"
#include "TH1.h"
#include "TH1F.h"
#include "TStyle.h"
#include "TGraphErrors.h"
#include "TLine.h"
#include "TPolyLine.h"


void smoothLowPtRho(){

	gROOT->SetBatch();

	bool divideByError=false;

	//char DataID[200];
	char directory[200];
	char Fig_directory[200];
	sprintf(Fig_directory,"FigBuffer/smoothLowPtRho");
	gSystem->mkdir(Fig_directory,true);

	int pTmin=2;
	int pTmax=5;
	int rapmin=1;
	int rapmax=1;
	int nStateMin=4;
	int nStateMax=5;

	//if(nState==4) { pTmin=3; pTmax=9; rapmin=1; rapmax=2;}
	//if(nState==5) { pTmin=2; pTmax=5; rapmin=1; rapmax=3;}
	//cout << "nState: " << nState << "\n"
	//	<< "pTmin: " << pTmin << "\n" 
	//	<< "pTmax: " << pTmax << "\n"
	//	<< "rapmin: " << rapmin << "\n"
	//	<< "rapmax: " << rapmax << endl;

	//char filename[500];
	//sprintf(DataID,"RhoFactor/ToyMCtest_22May");

	//if(divideByError)
	// 	sprintf(Fig_directory,"%s/%s_divideByError/Psi%dS",Fig_directory,DataID,nState-3);
	//else
	// 	sprintf(Fig_directory,"%s/%s/Psi%dS",Fig_directory,DataID,nState-3);
	//gSystem->mkdir(Fig_directory,true);

	//sprintf(filename,"Systematics/%s/TGraphResults_Psi%dS.root",DataID,nState-3);
	//cout << "filename: " << filename << endl;
	//TFile *infile = new TFile(filename,"R");
	//if(!infile) { cout<<"failed to open file.."<<endl; return; }

	TH1F *histMean[18];

	char frameName[100];
	char GraphName[100];
	char beginLamLabel[200];
	char endLamLabel[200];
	char axislabel[200];
	sprintf(beginLamLabel,"#Delta");
	sprintf(endLamLabel,"");


	for(int ipar=1; ipar<=18; ipar++){

		if(ipar==1)  sprintf(axislabel,"%s#lambda^{CS}_{#vartheta}%s",beginLamLabel,endLamLabel);
		if(ipar==2)  sprintf(axislabel,"%s#lambda^{CS}_{#varphi}%s",beginLamLabel,endLamLabel);
		if(ipar==3)  sprintf(axislabel,"%s#lambda^{CS}_{#vartheta#varphi}%s",beginLamLabel,endLamLabel);
		if(ipar==4)  sprintf(axislabel,"%s#lambda^{*CS}_{#vartheta}%s",beginLamLabel,endLamLabel);
		if(ipar==5)  sprintf(axislabel,"%s#lambda^{*CS}_{#varphi}%s",beginLamLabel,endLamLabel);
		if(ipar==6)  sprintf(axislabel,"%s#tilde{#lambda}^{CS}%s",beginLamLabel,endLamLabel);

		if(ipar==7)  sprintf(axislabel,"%s#lambda^{HX}_{#vartheta}%s",beginLamLabel,endLamLabel);
		if(ipar==8)  sprintf(axislabel,"%s#lambda^{HX}_{#varphi}%s",beginLamLabel,endLamLabel);
		if(ipar==9)  sprintf(axislabel,"%s#lambda^{HX}_{#vartheta#varphi}%s",beginLamLabel,endLamLabel);
		if(ipar==10) sprintf(axislabel,"%s#lambda^{*HX}_{#vartheta}%s",beginLamLabel,endLamLabel);
		if(ipar==11) sprintf(axislabel,"%s#lambda^{*HX}_{#varphi}%s",beginLamLabel,endLamLabel);
		if(ipar==12) sprintf(axislabel,"%s#tilde{#lambda}^{HX}%s",beginLamLabel,endLamLabel);

		if(ipar==13) sprintf(axislabel,"%s#lambda^{PX}_{#vartheta}%s",beginLamLabel,endLamLabel);
		if(ipar==14) sprintf(axislabel,"%s#lambda^{PX}_{#varphi}%s",beginLamLabel,endLamLabel);
		if(ipar==15) sprintf(axislabel,"%s#lambda^{PX}_{#vartheta#varphi}%s",beginLamLabel,endLamLabel);
		if(ipar==16) sprintf(axislabel,"%s#lambda^{*PX}_{#vartheta}%s",beginLamLabel,endLamLabel);
		if(ipar==17) sprintf(axislabel,"%s#lambda^{*PX}_{#varphi}%s",beginLamLabel,endLamLabel);
		if(ipar==18) sprintf(axislabel,"%s#tilde{#lambda}^{PX}%s",beginLamLabel,endLamLabel);

		if(divideByError) sprintf(axislabel, "%s / error", axislabel);


		int nBins=70;
		double xMin = -0.15, xMax = 0.2;
		if(divideByError) { nBins=200; xMin=-5.; xMax=5.; }

		if(ipar==2 || ipar==3 || ipar==8 || ipar==9  || ipar==14 || ipar==15) { xMin=-0.2; xMax=0.2; nBins=(xMax-xMin)/0.01; }
		if(ipar==1 || ipar==6 || ipar==7 || ipar==12 || ipar==13 || ipar==18) { xMin=-0.5; xMax=0.5; nBins=(xMax-xMin)/0.02; }
		cout << "nBins " << nBins << "\n" 
			<< "xMin " << xMin << "\n"
			<< "xMax " << xMax << endl;

		histMean[ipar-1] = new TH1F(Form("histMean_%d",ipar), "", nBins, xMin, xMax);
		histMean[ipar-1] -> GetXaxis() -> SetTitle(axislabel);
		histMean[ipar-1] -> GetXaxis() -> CenterTitle(true);
		histMean[ipar-1] -> GetXaxis() -> SetTitleSize(0.08);

		for(int nState=nStateMin; nState<=nStateMax; nState++){

			if(nState==4) { pTmin=3; pTmax=9; rapmin=1; rapmax=2;}
			if(nState==5) { pTmin=2; pTmax=5; rapmin=1; rapmax=3;}
			cout << "nState: " << nState << "\n"
				<< "pTmin: " << pTmin << "\n"
				<< "pTmax: " << pTmax << "\n"
				<< "rapmin: " << rapmin << "\n"
				<< "rapmax: " << rapmax << endl;


			/////////
			char DataID[200];
			sprintf(DataID,"RhoFactor/ToyMCtest_22May");
			if(nStateMin != nStateMax && nState==4){
				if(ipar==1) sprintf(Fig_directory,"%s/%s/PsiNS",Fig_directory,DataID);
			}

			if(nStateMin==nStateMax){
				if(divideByError){
					if(ipar==1) sprintf(Fig_directory,"%s/%s_divideByError/Psi%dS",Fig_directory,DataID,nState-3);
				}
				else{
					if(ipar==1) sprintf(Fig_directory,"%s/%s/Psi%dS",Fig_directory,DataID,nState-3);
				}
			}
			gSystem->mkdir(Fig_directory,true);

			char filename[500];
			sprintf(filename,"Systematics/%s/TGraphResults_Psi%dS.root",DataID,nState-3);
			cout << "filename: " << filename << endl;
			TFile *infile = new TFile(filename,"R");
			if(!infile) { cout<<"failed to open file.."<<endl; return; }
			///////////////


			for(int irap=rapmin; irap<=rapmax; irap++){

				if(ipar==1)  sprintf(GraphName,"lth_CS_rap%d",irap);
				if(ipar==2)  sprintf(GraphName,"lph_CS_rap%d",irap);
				if(ipar==3)  sprintf(GraphName,"ltp_CS_rap%d",irap);
				if(ipar==4)  sprintf(GraphName,"lthstar_CS_rap%d",irap);
				if(ipar==5)  sprintf(GraphName,"lphstar_CS_rap%d",irap);
				if(ipar==6)  sprintf(GraphName,"ltilde_CS_rap%d",irap);

				if(ipar==7)  sprintf(GraphName,"lth_HX_rap%d",irap);
				if(ipar==8)  sprintf(GraphName,"lph_HX_rap%d",irap);
				if(ipar==9)  sprintf(GraphName,"ltp_HX_rap%d",irap);
				if(ipar==10) sprintf(GraphName,"lthstar_HX_rap%d",irap);
				if(ipar==11) sprintf(GraphName,"lphstar_HX_rap%d",irap);
				if(ipar==12) sprintf(GraphName,"ltilde_HX_rap%d",irap);

				if(ipar==13) sprintf(GraphName,"lth_PX_rap%d",irap);
				if(ipar==14) sprintf(GraphName,"lph_PX_rap%d",irap);
				if(ipar==15) sprintf(GraphName,"ltp_PX_rap%d",irap);
				if(ipar==16) sprintf(GraphName,"lthstar_PX_rap%d",irap);
				if(ipar==17) sprintf(GraphName,"lphstar_PX_rap%d",irap);
				if(ipar==18) sprintf(GraphName,"ltilde_PX_rap%d",irap);


				cout << "GraphName: " << GraphName << endl;

				TGraphAsymmErrors* graph = (TGraphAsymmErrors*) infile->Get(GraphName);

				double ptCentre, lmean, lmeanerr;
				for(int ipt=pTmin; ipt<=pTmax; ipt++){
					graph->GetPoint(ipt-1, ptCentre, lmean);
					lmeanerr = graph->GetErrorY(ipt-1);

					if(lmean<-0.15 || lmean>0.2)
					 	cout << "lmean: " << lmean << endl;
					//cout << "lmean / error: " << lmean / lmeanerr << endl;

					if(divideByError)
						histMean[ipar-1]->Fill(lmean / lmeanerr);
					else
						histMean[ipar-1]->Fill(lmean);

				}//ipt

			}//irap


			//infile->Close();
		}//nState

	}//ipar

	gStyle->SetPadBottomMargin(0.21);
	gStyle->SetPadLeftMargin(0.08);
	gStyle->SetPadRightMargin(0.02);
	gStyle->SetPadTopMargin(0.05);
	gStyle->SetTitleFillColor(10);
	gStyle->SetTitleBorderSize(1);
	gStyle->SetOptFit(1);
	gStyle->SetOptStat(1);
	gStyle->SetTitleFont(22);
	gStyle->SetStatFont(22);
	gStyle->SetStatColor(10);
	gStyle->SetStatBorderSize(1);
	gStyle->SetLabelFont(22,"X");
	gStyle->SetLabelFont(22,"Y");
	gStyle->SetTitleXOffset(1.2);
	gStyle->SetTitleYOffset(1.2);gStyle->SetHistLineWidth(2);
	gStyle->SetStatX(0.9);gStyle->SetStatY(0.9);
	gStyle->SetTitleX(0.15);
	gStyle->SetTitleY(0.96);

	//TCanvas *plotCanvas = new TCanvas("plotCanvas","plotCanvas",1200,800);
	TCanvas *plotCanvas = new TCanvas("plotCanvas","plotCanvas",1000,800);
	gPad->SetFillColor(kWhite);

	char savename[200];
	for(int ipar=1; ipar<=18; ipar++){

		if(ipar==1)  sprintf(savename,"lth_CS");
		if(ipar==2)  sprintf(savename,"lph_CS");
		if(ipar==3)  sprintf(savename,"ltp_CS");
		if(ipar==4)  sprintf(savename,"lthstar_CS");
		if(ipar==5)  sprintf(savename,"lphstar_CS");
		if(ipar==6)  sprintf(savename,"ltilde_CS");

		if(ipar==7)  sprintf(savename,"lth_HX");
		if(ipar==8)  sprintf(savename,"lph_HX");
		if(ipar==9)  sprintf(savename,"ltp_HX");
		if(ipar==10) sprintf(savename,"lthstar_HX");
		if(ipar==11) sprintf(savename,"lphstar_HX");
		if(ipar==12) sprintf(savename,"ltilde_HX");

		if(ipar==13) sprintf(savename,"lth_PX");
		if(ipar==14) sprintf(savename,"lph_PX");
		if(ipar==15) sprintf(savename,"ltp_PX");
		if(ipar==16) sprintf(savename,"lthstar_PX");
		if(ipar==17) sprintf(savename,"lphstar_PX");
		if(ipar==18) sprintf(savename,"ltilde_PX");

		histMean[ipar-1]->Draw();
		plotCanvas->SaveAs(Form("%s/%s.pdf",Fig_directory,savename));

	}

	for(int ipar=1; ipar<=18; ipar++){

		if(ipar==1)  sprintf(savename,"lth_CS");
		if(ipar==2)  sprintf(savename,"lph_CS");
		if(ipar==3)  sprintf(savename,"ltp_CS");
		if(ipar==4)  sprintf(savename,"lthstar_CS");
		if(ipar==5)  sprintf(savename,"lphstar_CS");
		if(ipar==6)  sprintf(savename,"ltilde_CS");

		if(ipar==7)  sprintf(savename,"lth_HX");
		if(ipar==8)  sprintf(savename,"lph_HX");
		if(ipar==9)  sprintf(savename,"ltp_HX");
		if(ipar==10) sprintf(savename,"lthstar_HX");
		if(ipar==11) sprintf(savename,"lphstar_HX");
		if(ipar==12) sprintf(savename,"ltilde_HX");

		if(ipar==13) sprintf(savename,"lth_PX");
		if(ipar==14) sprintf(savename,"lph_PX");
		if(ipar==15) sprintf(savename,"ltp_PX");
		if(ipar==16) sprintf(savename,"lthstar_PX");
		if(ipar==17) sprintf(savename,"lphstar_PX");
		if(ipar==18) sprintf(savename,"ltilde_PX");

		cout << "hist RMS of " << savename << "   " << histMean[ipar-1]->GetRMS() << endl;
	}

}

