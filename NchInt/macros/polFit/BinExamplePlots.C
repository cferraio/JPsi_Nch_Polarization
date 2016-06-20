/*
 * BinExamplePlots.C
 *
 *  Created on: Apr 22, 2012
 *      Author: valentinknuenz
 */

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
#include "TStyle.h"
#include "TGraphErrors.h"
#include "TLine.h"
#include "TPolyLine.h"


void plotSummary( TTree* Results, char parameter[200], char parameterTitle[200], int iState , int irapBin, int ipTbin, char DataID[200] ,double nSigma) {

	int nBins=100;
	TH1D* h_param   = new TH1D( "h_param", "h_param", nBins, -1,1);
	char DrawChar[200];
	sprintf(DrawChar,"%s>>h_param",parameter);
	Results->Draw(DrawChar);
	TH1D* h_param_   = new TH1D( "h_param_", "h_param_", nBins, h_param->GetMean()-nSigma*h_param->GetRMS(), h_param->GetMean()+nSigma*h_param->GetRMS());
	sprintf(DrawChar,"%s>>h_param_",parameter);
	Results->Draw(DrawChar);

	TCanvas *c2 = new TCanvas("c2","c2",1200,800);
	h_param_->SetXTitle(parameterTitle);
	gPad->SetFillColor(kWhite);
	h_param_->SetTitle(0);
	h_param_->Draw();

	char savename[200];
	if(iState<=3) sprintf(savename,"FigBuffer/BinExamplePlots/%s/Ups%dS/%s_rap%d_pT%d.pdf",DataID,iState,parameter,irapBin,ipTbin);
	if(iState>3) sprintf(savename,"FigBuffer/BinExamplePlots/%s/Psi%dS/%s_rap%d_pT%d.pdf",DataID,iState-3,parameter,irapBin,ipTbin);
	c2->SaveAs(savename);

}


void BinExamplePlots(){

	gROOT->SetBatch();

	char DataID[200];
	//sprintf(DataID,"Psi2S_3.00Sigma_11Dec2012");
	//sprintf(DataID,"Psi2S_ctauScen0_FracLSB-1_rho_20Feb2013");

	bool PlotAcceptance=true;

	char directory[200];
	char Fig_directory[200];
	sprintf(Fig_directory,"FigBuffer/BinExamplePlots");
	gSystem->mkdir(Fig_directory,true);
	//sprintf(Fig_directory,"FigBuffer/BinExamplePlots/%s",DataID);
	//gSystem->mkdir(Fig_directory,true);

	bool PlotBGnormCorr=false;
	double normCorr1S[2][12], normCorr1SErr[2][12];
	double normCorr2S[3][5],  normCorr2SErr[3][5];
	TGraphErrors *BGnormCorr1S[2];
	TGraphErrors *BGnormCorr2S[3];
	double pT1S[2][12]={
		{11.0099, 12.9394, 14.9263, 16.9274, 18.9324, 20.9341, 23.3586, 27.1481, 32.1847, 37.2121, 44.0308, 56.8306},
		{10.9643, 12.9246, 14.9218, 16.9219, 18.9241, 20.9311, 23.3568, 27.1473, 32.1828, 37.2213, 44.0189, 56.8043}};
	double pT1SErr[2][12]={
		{0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.},
		{0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.}};
	double pT2S[3][5]= {
		{11.8943, 15.7294, 19.7371, 25.1272, 36.1102},
		{11.7868, 15.7024, 19.7270, 25.1119, 36.1893},
		{11.7788, 15.7210, 19.7275, 25.1509, 36.1740}};
	double pT2SErr[3][5]= {
		{0., 0., 0., 0., 0.},
		{0., 0., 0., 0., 0.},
		{0., 0., 0., 0., 0.}};

	int pTmin=6;
	int pTmax=10;
	int rapmin=1;
	int rapmax=2;
	int nStatemin=4;
	int nStatemax=5;

	for(int iState=nStatemin;iState<nStatemax+1;iState++){

		if(iState==4) {pTmin=3; pTmax=12; rapmin=1; rapmax=2;}
		if(iState==5) {pTmin=2; pTmax=5; rapmin=1; rapmax=3;}

		//sprintf(DataID,Form("Psi%dS_3.00Sigma_11Dec2012",iState-3));
		//sprintf(DataID,Form("Psi%dS_ctauScen0_FracLSB-1_rho_25Feb2013",iState-3));
		//sprintf(DataID,Form("Psi%dS_ctauScen0_FracLSB-1_rho_25Feb2013_massRange",iState-3));
		//sprintf(DataID,Form("Psi%dS_ctauScen0_FracLSB-1_rho_25Feb2013_massRange_Bin20_2_2",iState-3));
		//sprintf(DataID,Form("Psi%dS_ctauScen0_FracLSB-1_rho_25Feb2013_massRange_Bin5_8_8",iState-3));
		//sprintf(DataID,Form("Psi%dS_ctauScen0_FracLSB-1_rho_26Feb2013_BgNoRebin",iState-3));
		//sprintf(DataID,Form("Psi%dS_ctauScen0_FracLSB-1_16Mar2013",iState-3));
		//sprintf(DataID,Form("Psi%dS_ctauScen0_FracLSB-1_16Mar2013_2",iState-3));
		//sprintf(DataID,Form("Psi%dS_ctauScen0_FracLSB-1_16Mar2013_scaleFracBg",iState-3));
		//sprintf(DataID,Form("Psi%dS_ctauScen0_FracLSB-1_19Mar2013_0FracBg_0BgModel_1Rho",iState-3));
		//sprintf(DataID,Form("Psi%dS_ctauScen0_FracLSB-1_19Mar2013_0FracBg_1BgModel_0Rho",iState-3));
		//sprintf(DataID,Form("Psi%dS_ctauScen0_FracLSB-1_19Mar2013_1FracBg_0BgModel_0Rho",iState-3));
		//sprintf(DataID,Form("Psi%dS_ctauScen0_FracLSB-1_19Mar2013_1FracBg_1BgModel_1Rho",iState-3));
		//sprintf(DataID,Form("Psi%dS_ctauScen0_FracLSB-1_19Mar2013_1FracBg_1BgModel_1Rho_AlteredPPD_BKGlinPLUSRestSquaredGauss_5nRand_Apr5",iState-3));
		//sprintf(DataID,Form("Psi%dS_ctauScen3_FracLSB-1_29Apr2013_1FracBg_1BgModel_noRhoPt35",iState-3));
		///sprintf(DataID,Form("Psi%dS_ctauScen1_FracLSB-1_29Apr2013_1FracBg_1BgModel_noRhoPt35",iState-3));
		//sprintf(DataID,Form("Psi%dS_ctauScen5_FracLSB-1_7May2013_1FracBg_1BgModel",iState-3));
		//sprintf(DataID,Form("Psi%dS_ctauScen0_FracLSB75_19Mar2013_1FracBg_1BgModel_1Rho",iState-3));
		//sprintf(DataID,Form("Psi%dS_ctauScen5_FracLSB25_7May2013_1FracBg_1BgModel",iState-3));
		//sprintf(DataID,Form("Psi%dS_ctauScen5_FracLSB-1_22May2013_1FracBg_1BgModel_newEff",iState-3));
		sprintf(DataID,Form("Psi%dS_ctauScen5_FracLSB-1_26May2013_1FracBg_1BgModel_1RhoAll",iState-3));

		sprintf(Fig_directory,"FigBuffer/BinExamplePlots/%s",DataID);
		gSystem->mkdir(Fig_directory,true);

		if(iState<=3) sprintf(directory,"%s/Ups%dS",Fig_directory,iState);
		if(iState>3) sprintf(directory,"%s/Psi%dS",Fig_directory,iState-3);
		gSystem->mkdir(directory,true);

		for(int irapBin=rapmin;irapBin<rapmax+1;irapBin++){
			for(int ipTbin=pTmin;ipTbin<pTmax+1;ipTbin++){


				char Filename[200];
				if(iState<=3)
					sprintf(Filename,"/afs/ihep.ac.cn/users/z/zhangll/fs/work/polarization/PsiPol2011/Psi/Data/%s/results_%dSUps_rap%d_pT%d.root",DataID,iState,irapBin,ipTbin);
				if(iState>3)
					sprintf(Filename,"/afs/ihep.ac.cn/users/z/zhangll/fs/work/polarization/PsiPol2011/Psi/Data/%s/results_Psi%dS_rap%d_pT%d.root",DataID,iState-3,irapBin,ipTbin);
				TFile *resultFile = new TFile(Filename,"READ");

				TTree* Results=(TTree*)resultFile->Get("Results");
				TH1D* SubtractedBG_test=(TH1D*)resultFile->Get("SubtractedBG_test");
				TH1D* MetropolisHastingsAcceptanceCS=(TH1D*)resultFile->Get("MetropolisHastingsAcceptanceCS");
				TH1D* MetropolisHastingsAcceptanceHX=(TH1D*)resultFile->Get("MetropolisHastingsAcceptanceHX");
				TH1D* MetropolisHastingsAcceptancePX=(TH1D*)resultFile->Get("MetropolisHastingsAcceptancePX");

				double nSigma=4;

				for(int iLam=1;iLam<13;iLam++){

					char parameter[200];
					char parameterTitle[200];

					if(iLam==1) {sprintf(parameter,"lthCS"); sprintf(parameterTitle,"#lambda_{#theta}^{CS}");}
					if(iLam==2) {sprintf(parameter,"lphCS"); sprintf(parameterTitle,"#lambda_{#phi}^{CS}");}
					if(iLam==3) {sprintf(parameter,"ltpCS"); sprintf(parameterTitle,"#lambda_{#theta#phi}^{CS}");}
					if(iLam==4) {sprintf(parameter,"ltildeCS"); sprintf(parameterTitle,"#tilde{#lambda}^{CS}");}

					if(iLam==5) {sprintf(parameter,"lthHX"); sprintf(parameterTitle,"#lambda_{#theta}^{HX}");}
					if(iLam==6) {sprintf(parameter,"lphHX"); sprintf(parameterTitle,"#lambda_{#phi}^{HX}");}
					if(iLam==7) {sprintf(parameter,"ltpHX"); sprintf(parameterTitle,"#lambda_{#theta#phi}^{HX}");}
					if(iLam==8) {sprintf(parameter,"ltildeHX"); sprintf(parameterTitle,"#tilde{#lambda}^{HX}");}

					if(iLam==9) {sprintf(parameter,"lthPX"); sprintf(parameterTitle,"#lambda_{#theta}^{PX}");}
					if(iLam==10) {sprintf(parameter,"lphPX"); sprintf(parameterTitle,"#lambda_{#phi}^{PX}");}
					if(iLam==11) {sprintf(parameter,"ltpPX"); sprintf(parameterTitle,"#lambda_{#theta#phi}^{PX}");}
					if(iLam==12) {sprintf(parameter,"ltildePX"); sprintf(parameterTitle,"#tilde{#lambda}^{PX}");}

					plotSummary( Results, parameter, parameterTitle, iState , irapBin, ipTbin, DataID , nSigma);

				}



				double lth_min=-1.1;
				double lth_max=1.1;
				double lph_min=-1.1;
				double lph_max=1.1;

				TCanvas *c1 = new TCanvas("c1", "c1", 10, 28, 580,571);
				c1->Range(-237.541,-66.47556,187.377,434.8609);
				c1->SetFillColor(0);
				c1->SetBorderMode(0);
				c1->SetBorderSize(0);
				c1->SetLeftMargin(0.1354167);
				c1->SetRightMargin(0.01736111);
				c1->SetTopMargin(0.01841621);
				c1->SetBottomMargin(0.1325967);
				c1->SetFrameBorderMode(0);

				TH2D* h_lph_vs_lth_CS   = new TH2D( "h_lph_vs_lth_CS", "h_lph_vs_lth_CS", 100, lth_min,lth_max,100,lph_min,lph_max);
				Results->Draw("lphCS:lthCS>>h_lph_vs_lth_CS");
				TH2D* h_lph_vs_lth_HX   = new TH2D( "h_lph_vs_lth_HX", "h_lph_vs_lth_HX", 100,lth_min,lth_max,100,lph_min,lph_max);
				Results->Draw("lphHX:lthHX>> h_lph_vs_lth_HX");
				TH2D* h_lph_vs_lth_PX   = new TH2D( "h_lph_vs_lth_PX", "h_lph_vs_lth_PX", 100, lth_min,lth_max,100,lph_min,lph_max);
				Results->Draw("lphPX:lthPX>> h_lph_vs_lth_PX");

				h_lph_vs_lth_CS->SetXTitle("#lambda_{#theta}");
				h_lph_vs_lth_CS->SetYTitle("#lambda_{#phi}");


				h_lph_vs_lth_CS->SetTitle(0);
				h_lph_vs_lth_HX->SetTitle(0);
				h_lph_vs_lth_PX->SetTitle(0);
				h_lph_vs_lth_CS->SetStats(0);
				h_lph_vs_lth_HX->SetStats(0);
				h_lph_vs_lth_PX->SetStats(0);

				h_lph_vs_lth_CS->SetMarkerStyle(20);
				h_lph_vs_lth_HX->SetMarkerStyle(20);
				h_lph_vs_lth_PX->SetMarkerStyle(20);
				h_lph_vs_lth_CS->SetMarkerSize(0.5);
				h_lph_vs_lth_HX->SetMarkerSize(0.5);
				h_lph_vs_lth_PX->SetMarkerSize(0.5);
				h_lph_vs_lth_CS->SetMarkerColor(kBlue);
				h_lph_vs_lth_HX->SetMarkerColor(kRed);
				h_lph_vs_lth_PX->SetMarkerColor(kGreen);

				TH2D* h_lph_vs_lth = (TH2D*)h_lph_vs_lth_CS->Clone();
				h_lph_vs_lth->SetName("h_lph_vs_lth");
				h_lph_vs_lth->Reset();

				h_lph_vs_lth->GetXaxis()->SetTitle("#lambda_{#vartheta}");
				h_lph_vs_lth->GetXaxis()->SetLabelOffset(0.028);
				h_lph_vs_lth->GetXaxis()->SetTitleSize(0.05);
				h_lph_vs_lth->GetXaxis()->SetTickLength(-0.03);
				h_lph_vs_lth->GetXaxis()->SetTitleOffset(1.15);
				h_lph_vs_lth->GetYaxis()->SetTitle("#lambda_{#varphi}");
				h_lph_vs_lth->GetYaxis()->SetLabelOffset(0.032);
				h_lph_vs_lth->GetYaxis()->SetTitleSize(0.05);
				h_lph_vs_lth->GetYaxis()->SetTickLength(-0.03);
				h_lph_vs_lth->GetYaxis()->SetTitleOffset(1.3);
				h_lph_vs_lth->Draw("");

				double x_background_ph_vs_th[4] = { lth_min+0.01, lth_min+0.01, lth_max-0.01, lth_max-0.01 };
				double y_background_ph_vs_th[4] = { lph_min+0.01, lph_max-0.01, lph_max-0.01, lph_min+0.01 };
				TPolyLine *background_ph_vs_th = new TPolyLine( 4, x_background_ph_vs_th, y_background_ph_vs_th );
				background_ph_vs_th->SetFillColor(kGray);
				background_ph_vs_th->SetLineStyle(0);
				background_ph_vs_th->Draw("f same");

				double x_triangle_ph_vs_th[3] = {-1.,  1., 1.};
				double y_triangle_ph_vs_th[3] = { 0., -1., 1.};
				TPolyLine *triangle_ph_vs_th = new TPolyLine( 3, x_triangle_ph_vs_th, y_triangle_ph_vs_th );
				triangle_ph_vs_th->SetFillColor(kWhite);
				triangle_ph_vs_th->SetLineStyle(0);
				triangle_ph_vs_th->Draw("f same");

				TLine* lh_p1 = new TLine( lth_min, 1., lth_max, 1. );
				lh_p1->SetLineWidth( 1 );
				lh_p1->SetLineStyle( 2 );
				lh_p1->SetLineColor( kGray+2 );
				lh_p1->Draw( "same" );

				TLine* lh_m1 = new TLine( lth_min, -1., lth_max, -1. );
				lh_m1->SetLineWidth( 1 );
				lh_m1->SetLineStyle( 2 );
				lh_m1->SetLineColor( kGray+2 );
				lh_m1->Draw( "same" );

				TLine* lh_0 = new TLine( lth_min, 0., lth_max, 0. );
				lh_0->SetLineWidth( 1 );
				lh_0->SetLineStyle( 2 );
				lh_0->SetLineColor( kGray+2 );
				lh_0->Draw( "same" );

				TLine* lv_p1 = new TLine( 1., lph_min, 1., lph_max );
				lv_p1->SetLineWidth( 1 );
				lv_p1->SetLineStyle( 2 );
				lv_p1->SetLineColor( kGray+2 );
				lv_p1->Draw( "same" );

				TLine* lv_m1 = new TLine( -1., lph_min, -1., lph_max );
				lv_m1->SetLineWidth( 1 );
				lv_m1->SetLineStyle( 2 );
				lv_m1->SetLineColor( kGray+2 );
				lv_m1->Draw( "same" );

				TLine* lv_0 = new TLine( 0., lph_min, 0., lph_max );
				lv_0->SetLineWidth( 1 );
				lv_0->SetLineStyle( 2 );
				lv_0->SetLineColor( kGray+2 );
				lv_0->Draw( "same" );



				// CS frame

				/*      setContourHistogram ( h_lph_vs_lth_CS );
								h_lph_vs_lth_CS->SetLineColor( kBlue );
								h_lph_vs_lth_CS->SetLineWidth( 2 );
								h_lph_vs_lth_CS->SetLineStyle( 1 );
								h_lph_vs_lth_CS->Draw( "same" );

				// HX frame

				setContourHistogram ( h_lph_vs_lth_HX );
				h_lph_vs_lth_HX->SetLineColor( kRed );
				h_lph_vs_lth_HX->SetLineWidth( 2 );
				h_lph_vs_lth_HX->SetLineStyle( 1 );
				h_lph_vs_lth_HX->Draw( "cont2, same" );

				// PX frame

				setContourHistogram ( h_lph_vs_lth_PX );
				h_lph_vs_lth_PX->SetLineColor( kGreen+2 );
				h_lph_vs_lth_PX->SetLineWidth( 2 );
				h_lph_vs_lth_PX->SetLineStyle( 1 );
				h_lph_vs_lth_PX->Draw( "cont2, same" );
				*/

				h_lph_vs_lth_CS->Draw("same");
				h_lph_vs_lth_HX->Draw("same");
				h_lph_vs_lth_PX->Draw("same");

				TLegend* plotLegend=new TLegend(0.15,0.8,0.25,0.95);
				plotLegend->SetFillColor(kWhite);
				//    	 plotLegend->SetTextFont(72);
				plotLegend->SetTextSize(0.035);
				plotLegend->SetBorderSize(1);
				char legendentry[200];
				sprintf(legendentry,"CS");
				plotLegend->AddEntry(h_lph_vs_lth_CS,legendentry,"p");
				sprintf(legendentry,"HX");
				plotLegend->AddEntry(h_lph_vs_lth_HX,legendentry,"p");
				sprintf(legendentry,"PX");
				plotLegend->AddEntry(h_lph_vs_lth_PX,legendentry,"p");
				plotLegend->Draw(); plotLegend->Draw();

				char savename[200];
				if(iState<=3) sprintf(savename,"FigBuffer/BinExamplePlots/%s/Ups%dS/lth_vs_lph_rap%d_pT%d.pdf",DataID,iState,irapBin,ipTbin);
				if(iState>3)  sprintf(savename,"FigBuffer/BinExamplePlots/%s/Psi%dS/lth_vs_lph_rap%d_pT%d.pdf",DataID,iState-3,irapBin,ipTbin);
				c1->SaveAs(savename);



				/*      TCanvas *c2 = new TCanvas("c2","c2",800,800);
								gPad->SetFillColor(kWhite);
								h_lthlphCS->Draw();
								h_lthlphHX->Draw("same");
								h_lthlphPX->Draw("same");

*/

				TCanvas *c2 = new TCanvas("c2","c2",1200,800);
				SubtractedBG_test->SetXTitle("N_{BG,subtracted} / N_{BG,actual}");
				gPad->SetFillColor(kWhite);
				SubtractedBG_test->SetTitle(0);
				SubtractedBG_test->Draw();

				if(iState<=3) sprintf(savename,"FigBuffer/BinExamplePlots/%s/Ups%dS/SubtractedBackground_rap%d_pT%d.pdf",DataID,iState,irapBin,ipTbin);
				if(iState>3) sprintf(savename,"FigBuffer/BinExamplePlots/%s/Psi%dS/SubtractedBackground_rap%d_pT%d.pdf",DataID,iState-3,irapBin,ipTbin);
				c2->SaveAs(savename);


				///////////////////////
				if(PlotBGnormCorr){
					if(iState==4){
						normCorr1S   [irapBin-1][ipTbin-1] = SubtractedBG_test -> GetMean();
						normCorr1SErr[irapBin-1][ipTbin-1] = SubtractedBG_test -> GetRMS();
					}
					if(iState==5){
						normCorr2S   [irapBin-1][ipTbin-1] = SubtractedBG_test -> GetMean();
						normCorr2SErr[irapBin-1][ipTbin-1] = SubtractedBG_test -> GetRMS();
					}
				}
				///////////////////////

				if(PlotAcceptance){
					c2 = new TCanvas("c2","c2",1200,800);
					MetropolisHastingsAcceptanceCS->SetXTitle("Metr.Hast. acceptance, CS");
					gPad->SetFillColor(kWhite);
					MetropolisHastingsAcceptanceCS->SetTitle(0);
					MetropolisHastingsAcceptanceCS->Draw();

					if(iState<=3) sprintf(savename,"FigBuffer/BinExamplePlots/%s/Ups%dS/MetropolisHastingsAcceptanceCS_rap%d_pT%d.pdf",DataID,iState,irapBin,ipTbin);
					if(iState>3) sprintf(savename,"FigBuffer/BinExamplePlots/%s/Psi%dS/MetropolisHastingsAcceptanceCS_rap%d_pT%d.pdf",DataID,iState-3,irapBin,ipTbin);
					c2->SaveAs(savename);


					c2 = new TCanvas("c2","c2",1200,800);
					MetropolisHastingsAcceptanceHX->SetXTitle("Metr.Hast. acceptance, HX");
					gPad->SetFillColor(kWhite);
					MetropolisHastingsAcceptanceHX->SetTitle(0);
					MetropolisHastingsAcceptanceHX->Draw();

					if(iState<=3) sprintf(savename,"FigBuffer/BinExamplePlots/%s/Ups%dS/MetropolisHastingsAcceptanceHX_rap%d_pT%d.pdf",DataID,iState,irapBin,ipTbin);
					if(iState>3) sprintf(savename,"FigBuffer/BinExamplePlots/%s/Psi%dS/MetropolisHastingsAcceptanceHX_rap%d_pT%d.pdf",DataID,iState-3,irapBin,ipTbin);
					c2->SaveAs(savename);


					c2 = new TCanvas("c2","c2",1200,800);
					MetropolisHastingsAcceptancePX->SetXTitle("Metr.Hast. acceptance, PX");
					gPad->SetFillColor(kWhite);
					MetropolisHastingsAcceptancePX->SetTitle(0);
					MetropolisHastingsAcceptancePX->Draw();

					if(iState<=3) sprintf(savename,"FigBuffer/BinExamplePlots/%s/Ups%dS/MetropolisHastingsAcceptancePX_rap%d_pT%d.pdf",DataID,iState,irapBin,ipTbin);
					if(iState>3) sprintf(savename,"FigBuffer/BinExamplePlots/%s/Psi%dS/MetropolisHastingsAcceptancePX_rap%d_pT%d.pdf",DataID,iState-3,irapBin,ipTbin);
					c2->SaveAs(savename);
				}
				resultFile->Close();

			}  //ipTbin
		} //irapBin

	}  //iState

	if(PlotBGnormCorr){

		gStyle->SetPadBottomMargin(0.11);
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
		gStyle->SetTitleYOffset(1.2);
		gStyle->SetHistLineWidth(2);
		gStyle->SetStatX(0.9);
		gStyle->SetStatY(0.9);
		gStyle->SetTitleX(0.15);
		gStyle->SetTitleY(0.96);

		TCanvas *c3 = new TCanvas("c3","c3",1200,800);
		gPad->SetFillColor(kWhite);

		BGnormCorr1S[0] = new TGraphErrors(12,pT1S[0],normCorr1S[0],pT1SErr[0],normCorr1SErr[0]);
		BGnormCorr1S[1] = new TGraphErrors(12,pT1S[1],normCorr1S[1],pT1SErr[1],normCorr1SErr[1]);
		BGnormCorr2S[0] = new TGraphErrors(5,pT2S[0],normCorr2S[0],pT2SErr[0],normCorr2SErr[0]);
		BGnormCorr2S[1] = new TGraphErrors(5,pT2S[1],normCorr2S[1],pT2SErr[1],normCorr2SErr[1]);
		BGnormCorr2S[2] = new TGraphErrors(5,pT2S[2],normCorr2S[2],pT2SErr[2],normCorr2SErr[2]);

		BGnormCorr1S[0]->SetTitle("");
		BGnormCorr1S[0]->GetXaxis()->SetTitle("p_{T} (GeV)");
		BGnormCorr1S[0]->GetYaxis()->SetTitle("BGnormCorr");
		BGnormCorr1S[0]->GetYaxis()->SetRangeUser(0.6, 1.2);
		BGnormCorr1S[0]->GetXaxis()->SetLimits(6., 72.);

		BGnormCorr1S[0]->SetMarkerStyle(20);
		BGnormCorr1S[0]->SetMarkerColor(4);
		BGnormCorr1S[0]->SetLineColor(4);
		BGnormCorr1S[1]->SetMarkerStyle(24);
		BGnormCorr1S[1]->SetMarkerSize(1.2);
		BGnormCorr1S[1]->SetMarkerColor(4);
		BGnormCorr1S[1]->SetLineColor(4);

		BGnormCorr2S[0]->SetMarkerStyle(21);
		BGnormCorr2S[0]->SetMarkerColor(2);
		BGnormCorr2S[0]->SetLineColor(2);
		BGnormCorr2S[1]->SetMarkerStyle(25);
		BGnormCorr2S[1]->SetMarkerSize(1.2);
		BGnormCorr2S[1]->SetMarkerColor(2);
		BGnormCorr2S[1]->SetLineColor(2);
		BGnormCorr2S[2]->SetMarkerStyle(22);
		BGnormCorr2S[2]->SetMarkerSize(1.2);
		BGnormCorr2S[2]->SetMarkerColor(kMagenta+1);
		BGnormCorr2S[2]->SetLineColor(kMagenta+1);

		double blX = 0.7, blY = 0.7, trX = 0.95, trY = 0.93;
		TLegend* legend=new TLegend(blX,blY,trX,trY);
		legend->SetFillColor(kWhite);
		legend->SetTextFont(42);
		legend->SetTextSize(0.035);
		legend->SetBorderSize(0.);
		legend->AddEntry(BGnormCorr1S[0],"J/#psi |y| < 0.6","lp");
		legend->AddEntry(BGnormCorr1S[1],"J/#psi 0.6 < |y| < 1.2","lp");
		legend->AddEntry(BGnormCorr2S[0],"#psi(2S) |y| < 0.6","lp");
		legend->AddEntry(BGnormCorr2S[1],"#psi(2S) 0.6 < |y| < 1.2","lp");
		legend->AddEntry(BGnormCorr2S[2],"#psi(2S) 1.2 < |y| < 1.5","lp");

		TLine *line = new TLine(6., 1., 72., 1.);
		line->SetLineWidth(1.5); line->SetLineColor(kBlack); line->SetLineStyle(7);


		BGnormCorr1S[0]->Draw("AP");
		BGnormCorr1S[1]->Draw("P");
		BGnormCorr2S[0]->Draw("P");
		BGnormCorr2S[1]->Draw("P");
		BGnormCorr2S[2]->Draw("P");
		line->Draw("same");
		legend->Draw("same");

		c3 -> SaveAs(Form("FigBuffer/BinExamplePlots/%s/BGnormCorr.pdf",DataID));

	} //PlotBGnormCorr

	return;
}
