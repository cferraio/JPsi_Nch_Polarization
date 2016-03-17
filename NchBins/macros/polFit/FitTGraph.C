#include <iostream>
#include <string>
#include <sstream>
using namespace std;

#include "rootIncludes.inc"
#include "commonVar.h"
#include "ToyMC.h"

#include "TSystem.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TFile.h"
#include "TFrame.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TROOT.h"

#include "TH2F.h"
#include "TTree.h"
#include "TChain.h"
#include "TDirectory.h"
#include "TLegend.h"
#include "Riostream.h"
#include "TSystem.h"
#include "TString.h"
#include "TFile.h"
#include "TROOT.h"
#include "TTree.h"
#include "TMath.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "TRandom3.h"
#include "TF1.h"
#include "TF2.h"
#include "TH1.h"
#include "TH2.h"



void FitTGraph() {


	int ptBinMin=6;
	int ptBinMax=10;
	int nState=1;



	    char dirstruct[200];
		sprintf(dirstruct,"/afs/hephy.at/scratch/k/knuenz/CMSSW_4_2_4_patch2/src/UpsilonPol/macros/polFit/Systematics/TotalSyst/Apr25CentralWithStatSystSquaredSymmetrical");

		char filename[200];
		sprintf(filename,"%s/TGraphResults_%dSUps.root",dirstruct,nState);
		TFile *infile1 = new TFile(filename,"READ");

		char GraphName[200];

		for(int iLam = 1; iLam<19; iLam++){

		for(int rapBin = 1; rapBin < 3; rapBin++){


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

		int nFrame=0;

		if(iLam>0&&iLam<7) nFrame=1;
		if(iLam>6&&iLam<13) nFrame=2;
		if(iLam>12&&iLam<19) nFrame=3;

		TGraphAsymmErrors* graph1 = (TGraphAsymmErrors*) infile1->Get(GraphName);

		graph1->RemovePoint(0);
		graph1->RemovePoint(0);
		graph1->RemovePoint(0);
		graph1->RemovePoint(0);
		graph1->RemovePoint(0);


		double yMin;
		double yMax;

		yMin=-2.1;
		yMax=2.1;

		if(iLam==2||iLam==8||iLam==14||iLam==3||iLam==9||iLam==15){
			yMin=-0.55;
			yMax=0.55;
		}

		if(iLam==6||iLam==12||iLam==18){
			yMin=-2.1;
			yMax=2.1;

			double pTmin=10;
			double pTmax=50;

	      double FontSize=0.0215;
	      double xText=pTmax+(pTmin-plopTmintQmin)*0.025;

	      char FitOptions[200];
	      sprintf(FitOptions,"EFNR");

	    TCanvas* ScaleCanvas = new TCanvas("ScaleCanvas","ScaleCanvas",1000, 800);
	    ScaleCanvas->SetFillColor(kWhite);
	    ScaleCanvas->cd(1);
	    ScaleCanvas->SetLeftMargin(0.2);
	    gPad->SetFillColor(kWhite);

		  /////////////////////
		  bool DrawExt=true;

			TH1F *plotHisto = new TH1F;
			plotHisto = ScaleCanvas->DrawFrame(plotQmin,plotPESmin,plotQmax,plotPESmax);
			plotHisto->SetXTitle("Q [GeV]");
			plotHisto->SetYTitle("Photon energy scale");
			plotHisto->GetYaxis()->SetTitleOffset(2);


		  PESgraph->SetMarkerColor(kGreen-2);
		  PESgraph->SetMarkerStyle(21);
		  PESgraph->SetTitle(0);
		  PESgraph->Draw("P");
		  PESgraph->SaveAs("ScalePlots/PESgraph.root");


	      TF1* f1PES = new TF1("f1PES","pol1",plotQmin,plotQmax);
	      PESgraph->Fit("f1PES",FitOptions);
	      f1PES->SetLineWidth(0.4);
	      f1PES->SetLineColor(kRed);
	      f1PES->SetLineStyle(1);
//	      f1PES->Draw("same");

	      double PES_p0 = f1PES->GetParameter(0);
	      double err_PES_p0 = f1PES->GetParError(0);
	      double PES_p1 = f1PES->GetParameter(1);
	      double err_PES_p1 = f1PES->GetParError(1);
	      double PES_chi2=f1PES->GetChisquare();
	      double PES_NDF=f1PES->GetNDF();
	      double PES_BIC=PES_chi2+2*TMath::Log(4);

	      double highest=plotPESmax-(plotPESmax-plotPESmin)*0.05;

	      sprintf(text,"#color[2]{Fitting pol1:} p0 = %1.4f #pm %1.4f, p1 = %1.4f #pm %1.4f, #chi^{2}/ndf = %1.4f / %d, BIC = %1.3f", PES_p0, err_PES_p0, PES_p1, err_PES_p1,PES_chi2,PES_NDF,PES_BIC);
	      TLatex textPES1 = TLatex(xText,highest,text);
	      textPES1.SetTextSize(FontSize)                                                                                                                                                                                                                                             ;
	      textPES1.Draw( "same" )                                                                                                                                                                                                                                                 ;

	      TF1* f1PESconst = new TF1("f1PESconst","pol0",plotQmin,plotQmax);
	      PESgraph->Fit("f1PESconst",FitOptions);
	      f1PESconst->SetLineWidth(0.4);
	      f1PESconst->SetLineStyle(1);
	      f1PESconst->SetLineColor(kBlue);
	      f1PESconst->Draw("same");

	      double PES_const_p0 = f1PESconst->GetParameter(0);
	      double err_PES_const_p0 = f1PESconst->GetParError(0);
	      double PES_const_chi2=f1PESconst->GetChisquare();
	      double PES_const_NDF=f1PESconst->GetNDF();
	      double PES_const_BIC=PES_const_chi2+1*TMath::Log(4);

	      TF1* f1PES_m = new TF1("f1PES_m","pol0",plotQmin,plotQmax);
	      f1PES_m->SetParameter(0,PES_const_p0+err_PES_const_p0);
	      f1PES_m->SetLineWidth(0.2);
	      f1PES_m->SetLineColor(kBlue);
	      f1PES_m->SetLineStyle(2);
	      if(DrawExt) f1PES_m->Draw("same");
	      TF1* f1PES_p = new TF1("f1PES_p","pol0",plotQmin,plotQmax);
	      f1PES_p->SetParameter(0,PES_const_p0-err_PES_const_p0);
	      f1PES_p->SetLineWidth(0.2);
	      f1PES_p->SetLineColor(kBlue);
	      f1PES_p->SetLineStyle(2);
	      if(DrawExt) f1PES_p->Draw("same");


	      highest=plotPESmax-(plotPESmax-plotPESmin)*0.1;

	      sprintf(text,"#color[4]{Fitting constant:} PES = %1.4f #pm %1.4f, #chi^{2}/ndf = %1.4f / %d, BIC = %1.3f", PES_const_p0, err_PES_const_p0, PES_const_chi2, PES_const_NDF,PES_const_BIC);
	      TLatex textPES2 = TLatex(xText,highest,text);
	      textPES2.SetTextSize(FontSize)                                                                                                                                                                                                                                             ;
	      textPES2.Draw( "same" )                                                                                                                                                                                                                                                 ;



	      highest=plotPESmax-(plotPESmax-plotPESmin)*0.95; TLatex textPESchi1 = TLatex(Q[0]-0.05,highest,"#chi_{c1}"); textPESchi1.SetTextSize(FontSize); textPESchi1.Draw( "same" );
	      highest=plotPESmax-(plotPESmax-plotPESmin)*0.9; TLatex textPESchi2 = TLatex(Q[1]-0.025,highest,"#chi_{b}(1P)"); textPESchi2.SetTextSize(FontSize); textPESchi2.Draw( "same" );
	      highest=plotPESmax-(plotPESmax-plotPESmin)*0.95; TLatex textPESchi3 = TLatex(Q[2]+0.025,highest,"#chi_{c2}"); textPESchi3.SetTextSize(FontSize); textPESchi3.Draw( "same" );
	      highest=plotPESmax-(plotPESmax-plotPESmin)*0.9; TLatex textPESchi4 = TLatex(Q[3]-0.025,highest,"#chi_{b}(2P)"); textPESchi4.SetTextSize(FontSize); textPESchi4.Draw( "same" );

		  PESgraph->Draw("P");
	      ScaleCanvas->SaveAs("ScalePlots/PESfit.pdf");





		}


		}





}
