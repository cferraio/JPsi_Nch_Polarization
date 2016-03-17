#include <iostream>
#include <string>
#include <sstream>
using namespace std;

#include "rootIncludes.inc"
#include "commonVar.h"
#include "ToyMC.h"

// calculation of the frame-invariants lambdatheta^star and lambdaphi^star:

void calcLambdastar( double& lthstar, double& lphstar,
		double& lth,     double& lph,     double& ltp ) {

	double LamPlus = 0.25 * ( lth - lph + TMath::Sqrt(TMath::Power( lth - lph, 2. ) + 4. * ltp*ltp ) );
	double LamMnus = 0.25 * ( lth - lph - TMath::Sqrt(TMath::Power( lth - lph, 2. ) + 4. * ltp*ltp ) );

	lthstar = ( lth - 3.*LamPlus ) / ( 1. + LamPlus );
	lphstar = ( lph + LamPlus )    / ( 1. + LamPlus );

	double lphstarMnus = ( lph + LamMnus ) / ( 1. + LamMnus );

	if ( TMath::Abs(lphstarMnus) < TMath::Abs(lphstar) ) {
		lthstar = ( lth - 3.*LamMnus ) / ( 1. + LamMnus );
		lphstar = lphstarMnus;
	}
}


const int NrapBins_=2;
const int NptBins_=2;
const int NcpmBins_=2;

double Gauss(double* x, double* par)
{
	double a=2*TMath::Pi();
	double gaussval=par[2]/(par[1]*TMath::Sqrt(a))*TMath::Exp(-pow(x[0]-par[0],2)/(2*pow(par[1],2)));
	return gaussval;
}



void PlotObject(TH1D* object_histo, char Xaxistitle[200], double l_min_object, double l_max_object, double l_step_1D_object, double object_histo_mean, double object_histo_meanerr, double object_histo_sigma, double object_histo_sigmaerr, char filename [200], int nGenerations, bool pull, double injvalue){

	TCanvas* c1 = new TCanvas("c1", "c1", 10, 28, 588,563);
	c1->Range(-1.370833,-114.012,1.029167,745.8285);
	c1->SetFillColor(0);
	c1->SetBorderMode(0);
	c1->SetBorderSize(0);
	c1->SetLeftMargin(0.10215278);
	c1->SetRightMargin(0.10215278);
	c1->SetTopMargin(0.01841621);
	c1->SetBottomMargin(0.1325967);
	c1->SetFrameBorderMode(0);

	double plotborder = 1.1;
	double plotMax = plotborder * object_histo->GetMaximum();
	object_histo->SetMinimum(0.);
	object_histo->SetMaximum(plotMax);

	object_histo->GetXaxis()->SetTitle(Xaxistitle);
	object_histo->GetXaxis()->SetLabelOffset(0.028);
	object_histo->GetXaxis()->SetTitleSize(0.05);
	object_histo->GetXaxis()->SetTickLength(-0.03);
	object_histo->GetXaxis()->SetTitleOffset(1.20);
	//	 object_histo->GetYaxis()->SetTitle("Pull Distribution");
	object_histo->GetYaxis()->SetLabelOffset(0.032);
	object_histo->GetYaxis()->SetTitleSize(0.05);
	object_histo->GetYaxis()->SetTickLength(-0.03);
	object_histo->GetYaxis()->SetTitleOffset(1.55);
	object_histo->Draw("");
	object_histo->Draw( "same" );

	TF1* gauss = new TF1("gauss",Gauss,l_min_object,l_max_object,3);
	gauss->SetParameter(0,object_histo_mean);
	gauss->SetParameter(1,object_histo_sigma);
	gauss->SetParameter(2,nGenerations*l_step_1D_object);
	gauss->SetLineColor( kGreen+2 );
	gauss->SetLineWidth( 2 );
	gauss->SetLineStyle( 1 );
	gauss->Draw( "same" );

	if (pull){
		TF1* gaussRef = new TF1("gaussRef",Gauss,l_min_object,l_max_object,3);
		gaussRef->SetParameter(0,0);
		gaussRef->SetParameter(1,1);
		gaussRef->SetParameter(2,nGenerations*l_step_1D_object);
		gaussRef->SetLineColor( kRed );
		gaussRef->SetLineWidth( 1 );
		gaussRef->SetLineStyle( 3 );
		gaussRef->Draw( "same" );
	}

	if (!pull){
		TLine* injLine = new TLine( injvalue, 0, injvalue, plotMax );
		injLine->SetLineWidth( 1 );
		injLine->SetLineStyle( 2 );
		injLine->SetLineColor( kRed );
		injLine->Draw( "same" );
	}
	char meantext[200], sigmatext[200];
	if (!pull) sprintf(meantext, "#splitline{#mu_{#lambda} = %1.3f #pm %1.3f}{#sigma_{#lambda} = %1.3f #pm %1.3f}",object_histo_mean,object_histo_meanerr,object_histo_sigma,object_histo_sigmaerr);
	if (pull) sprintf(meantext, "#splitline{#mu_{z} = %1.3f #pm %1.3f}{#sigma_{z} = %1.3f #pm %1.3f}",object_histo_mean,object_histo_meanerr,object_histo_sigma,object_histo_sigmaerr);

	TLatex *text = new TLatex(l_min_object+0.05*(l_max_object-l_min_object),plotMax*0.9,meantext);
	text->SetTextSize(0.035);
	text->Draw( "same" );

	if (!pull){
		sprintf(meantext, "#splitline{#lambda^{Inj.} = %1.3f}{#delta_{#lambda} = %1.1f #sigma_{#lambda}}",injvalue,TMath::Abs(object_histo_mean-injvalue)/object_histo_sigma);
		text = new TLatex(l_min_object+0.05*(l_max_object-l_min_object),plotMax*0.75,meantext);
		text->SetTextSize(0.035);
		text->Draw( "same" );
	}
	if (pull){
		sprintf(meantext, "#splitline{#delta_{#mu_{z}} = %1.1f #sigma_{#mu_{z}}}{#delta_{#sigma_{z}} = %1.1f #sigma_{#sigma_{z}}}",object_histo_mean/object_histo_meanerr,TMath::Abs(object_histo_sigma-1)/object_histo_sigmaerr);
		text = new TLatex(l_min_object+0.05*(l_max_object-l_min_object),plotMax*0.75,meantext);
		text->SetTextSize(0.035);
		text->Draw( "same" );
	}
	c1->Print( filename );

	delete gauss ;
	delete c1 ;


}


void PlotRapPt(int ptBin, int cpmBinMin, int cpmBinMax, int rapBin, double yMin, double yMax, char yAxisTitle[200], char filename[200], double lmean[ToyMC::ncpmBins], double lmeanerr_low[ToyMC::ncpmBins], double lmeanerr_high[ToyMC::ncpmBins], bool DrawLine, double lamline, char GraphFolder[200], char GraphName[200], bool SaveGraph, int nState, int nFrame, bool RealData, char* realdatadir, char* TreeID){

	gStyle->SetPalette(1,0);
	gStyle->SetPadBottomMargin(0.12);
	gStyle->SetPadLeftMargin(0.13);
	gStyle->SetPadRightMargin(0.15);

	gStyle->SetTickLength(-0.02, "xyz");
	gStyle->SetLabelOffset(0.02, "x");
	gStyle->SetLabelOffset(0.02, "y");
	gStyle->SetTitleOffset(1.3, "x");
	gStyle->SetTitleOffset(1.4, "y");
	gStyle->SetTitleFillColor(kWhite);

	TCanvas *plotCanvas = new TCanvas("plotCanvas","plotCanvas",1000,800);

	plotCanvas->SetFillColor(kWhite);
	plotCanvas->SetGrid();
	plotCanvas->GetFrame()->SetFillColor(kWhite);
	plotCanvas->GetFrame()->SetBorderSize(0);
	plotCanvas->SetRightMargin(0.05) ;

	TH1F *plotHisto = new TH1F;
	plotHisto = plotCanvas->DrawFrame(onia::cpmRange[cpmBinMin-1],yMin,onia::cpmRange[cpmBinMax],yMax);
	plotHisto->SetXTitle("N_{ch}");
	plotHisto->SetYTitle(yAxisTitle);
	plotHisto->GetYaxis()->SetTitleOffset(1.5);


	double yLine=lamline;

	TLine* RefLine = new TLine( onia::cpmRange[cpmBinMin-1], yLine, onia::cpmRange[cpmBinMax], yLine );
	RefLine->SetLineWidth( 2 );
	RefLine->SetLineStyle( 2 );
	RefLine->SetLineColor( kGreen+2 );
	if(DrawLine) RefLine->Draw( "same" );
	
	int nBinscpm=cpmBinMax-cpmBinMin+1;
//	int nBinspT=ptBinMax-ptBinMin+1;
	double cpmCentre_[nBinscpm];
	double cpmCentreErr_low[nBinscpm];
	double cpmCentreErr_high[nBinscpm];

	
	
	
	int cpm=0;
	  for(int cpmBin = cpmBinMin; cpmBin < cpmBinMax+1; cpmBin++) {
		double mean_cpm=0;
		if(RealData){
			char rootFilename [500];
			sprintf(rootFilename,"%s/data_%s_rap%d_pT%d_cpm%d.root",realdatadir,TreeID,rapBin,ptBin,cpmBin);
			TFile* dataFile = new TFile(rootFilename,"R");
			if(!dataFile) cout<<"dataerror"<<endl;
			TH1D* hist_mean_cpm = (TH1D*)dataFile->Get("mean_cpm");
			if(!hist_mean_cpm) cout<<"othererror"<<endl;
			mean_cpm = hist_mean_cpm->GetBinContent(1);
			//cout<<"rap"<<rapBin<<"_pt"<<ptBin<<" mean_pT: "<<mean_pT<<endl;
			dataFile->Close();
		}
		else
			mean_cpm = ToyMC::cpmCentre[rapBin-1][ptBin-1][cpmBin-1];

		cpmCentre_[cpm]=mean_cpm;
		cpmCentreErr_low[cpm]=TMath::Abs(mean_cpm-onia::cpmRange[cpmBin-1]);
		cpmCentreErr_high[cpm]=TMath::Abs(mean_cpm-onia::cpmRange[cpmBin]);
		cpm++;
		}


	TGraphAsymmErrors *plotGraph = new TGraphAsymmErrors(nBinscpm,cpmCentre_,lmean,cpmCentreErr_low,cpmCentreErr_high,lmeanerr_low,lmeanerr_high);
	plotGraph->SetMarkerColor(ToyMC::MarkerColor[nFrame]);
	plotGraph->SetLineColor(ToyMC::MarkerColor[nFrame]);
	plotGraph->SetMarkerStyle(ToyMC::MarkerStyle[nState][rapBin]);
	plotGraph->SetMarkerSize(ToyMC::MarkerSize[nState][rapBin]);
	plotGraph->SetName(GraphName);
	plotGraph->Draw("P");
	plotGraph->Write();



	char texTex[200];
	if(rapBin==1) sprintf(texTex,"      |y| < 0.6");
	if(rapBin==2) sprintf(texTex,"0.6 < |y| < 1.2");
	if(rapBin==3) sprintf(texTex,"1.2 < |y| < 1.5");
	TLatex *text = new TLatex(onia::cpmRange[cpmBinMax]*0.22,yMin+(yMax-yMin)*0.1,texTex);
	text->SetTextSize(0.035);
	text->Draw( "same" );

	plotCanvas->SetLogx(false);
	plotCanvas->SaveAs(filename);
	plotCanvas->Close();

	delete plotCanvas;

	if(SaveGraph){
		char outfilename[200];
		sprintf(outfilename,"%s/TGraphResults_Psi%dS_temp.root",GraphFolder,nState-3);cout<<outfilename<<endl;
		TFile *outfile = new TFile(outfilename,"UPDATE");

		outfile->cd();
		plotGraph->Draw("P");
		plotGraph->Write();

		outfile->Write();
		outfile->Close();
		delete outfile;
		outfile = NULL;
	}


}

void PlotComparisonRapPt(int ptBin, int cpmBinMin, int cpmBinMax, int rapBin, double yMin, double yMax, char yAxisTitle[200], char filename[200], double lmean1[ToyMC::ncpmBins], double lmean1err_low[ToyMC::ncpmBins], double lmean1err_high[ToyMC::ncpmBins], double lmean2[ToyMC::ncpmBins], double lmean2err_low[ToyMC::ncpmBins], double lmean2err_high[ToyMC::ncpmBins], double lmean3[ToyMC::ncpmBins], double lmean3err_low[ToyMC::ncpmBins], double lmean3err_high[ToyMC::ncpmBins], double lamline, int nState, bool RealData, char* realdatadir, char* TreeID){

	gStyle->SetPalette(1,0);
	gStyle->SetPadBottomMargin(0.12);
	gStyle->SetPadLeftMargin(0.13);
	gStyle->SetPadRightMargin(0.15);

	gStyle->SetTickLength(-0.02, "xyz");
	gStyle->SetLabelOffset(0.02, "x");
	gStyle->SetLabelOffset(0.02, "y");
	gStyle->SetTitleOffset(1.3, "x");
	gStyle->SetTitleOffset(1.4, "y");
	gStyle->SetTitleFillColor(kWhite);

	TCanvas *plotCanvas = new TCanvas("plotCanvas","plotCanvas",1000,800);

	plotCanvas->SetFillColor(kWhite);
	plotCanvas->SetGrid();
	plotCanvas->GetFrame()->SetFillColor(kWhite);
	plotCanvas->GetFrame()->SetBorderSize(0);
	plotCanvas->SetRightMargin(0.05) ;

	TH1F *plotHisto = new TH1F;
	plotHisto = plotCanvas->DrawFrame(onia::cpmRange[cpmBinMin-1],yMin,onia::cpmRange[cpmBinMax],yMax);
	plotHisto->SetXTitle("N_{ch}");
	plotHisto->SetYTitle(yAxisTitle);
	plotHisto->GetYaxis()->SetTitleOffset(1.5);

	TLegend* plotLegend=new TLegend(0.825,0.75,0.95,0.9);
	plotLegend->SetFillColor(kWhite);
	plotLegend->SetTextFont(72);
	plotLegend->SetTextSize(0.035);
	plotLegend->SetBorderSize(1);

	char legendentry[200];

	int nBinscpm=cpmBinMax-cpmBinMin+1;
//	int nBinspT=ptBinMax-ptBinMin+1;
	double cpmCentre_[nBinscpm];
	double cpmCentreErr_low[nBinscpm];
	double cpmCentreErr_high[nBinscpm];

	double yLine=lamline;

	TLine* RefLine = new TLine( onia::cpmRange[cpmBinMin-1], yLine, onia::cpmRange[cpmBinMax], yLine );
	RefLine->SetLineWidth( 2 );
	RefLine->SetLineStyle( 2 );
	RefLine->SetLineColor( kGreen+2 );
	RefLine->Draw( "same" );


	int cpm=0;
	  for(int cpmBin = cpmBinMin; cpmBin < cpmBinMax+1; cpmBin++) {
		double mean_cpm=0;
		if(RealData){
			char rootFilename [500];
			sprintf(rootFilename,"%s/data_%s_rap%d_pT%d_cpm%d.root",realdatadir,TreeID,rapBin,ptBin,cpmBin);
			TFile* dataFile = new TFile(rootFilename,"R");
			if(!dataFile) cout<<"error"<<endl;
			TH1D* hist_mean_cpm = (TH1D*)dataFile->Get("mean_cpm");
			if(!hist_mean_cpm) cout<<"error"<<endl;
			mean_cpm = hist_mean_cpm->GetBinContent(1);
			//cout<<"rap"<<rapBin<<"_pt"<<ptBin<<" mean_pT: "<<mean_pT<<endl;
			dataFile->Close();
		}
		else
			mean_cpm = ToyMC::cpmCentre[rapBin-1][ptBin-1][cpmBin-1];

		cpmCentre_[cpm]= mean_cpm + mean_cpm * 0.02;
		cpmCentreErr_low[cpm]=TMath::Abs(mean_cpm-onia::cpmRange[cpmBin-1]);
		cpmCentreErr_high[cpm]=TMath::Abs(mean_cpm-onia::cpmRange[cpmBin]);
		cpm++;
		}


	TGraphAsymmErrors *plotGraph = new TGraphAsymmErrors(nBinscpm,cpmCentre_,lmean1,cpmCentreErr_low,cpmCentreErr_high,lmean1err_low,lmean1err_high);
	plotGraph->SetMarkerColor(ToyMC::MarkerColor[1]);
	plotGraph->SetLineColor(ToyMC::MarkerColor[1]);
	plotGraph->SetMarkerStyle(ToyMC::MarkerStyle[nState][rapBin]);
	plotGraph->SetMarkerSize(ToyMC::MarkerSize[nState][rapBin]);
	plotGraph->Draw("P");
	sprintf(legendentry,"CS");
	plotLegend->AddEntry(plotGraph,legendentry,"ple");

	cpm=0;
	  for(int cpmBin = cpmBinMin; cpmBin < cpmBinMax+1; cpmBin++) {
		cpmCentre_[cpm]=ToyMC::cpmCentre[rapBin-1][ptBin-1][cpmBin-1];
		cpmCentreErr_low[cpm]=TMath::Abs(ToyMC::cpmCentre[rapBin-1][ptBin-1][cpmBin-1]-onia::cpmRange[cpmBin]);
		cpmCentreErr_high[cpm]=TMath::Abs(ToyMC::cpmCentre[rapBin-1][ptBin-1][cpmBin-1]-onia::cpmRange[cpmBin]);
		cpm++;
		}


	plotGraph = new TGraphAsymmErrors(nBinscpm,cpmCentre_,lmean2,cpmCentreErr_low,cpmCentreErr_high,lmean2err_low,lmean2err_high);
	plotGraph->SetMarkerColor(ToyMC::MarkerColor[2]);
	plotGraph->SetLineColor(ToyMC::MarkerColor[2]);
	plotGraph->SetMarkerStyle(ToyMC::MarkerStyle[nState][rapBin]);
	plotGraph->SetMarkerSize(ToyMC::MarkerSize[nState][rapBin]);
	plotGraph->Draw("P");
	sprintf(legendentry,"HX");
	plotLegend->AddEntry(plotGraph,legendentry,"ple");

	
	cpm=0;
	  for(int cpmBin = cpmBinMin; cpmBin < cpmBinMax+1; cpmBin++) {
		cpmCentre_[cpm]=ToyMC::cpmCentre[rapBin-1][ptBin-1][cpmBin-1]-ToyMC::cpmCentre[rapBin-1][ptBin-1][cpmBin-1]*.02;
		cpmCentreErr_low[cpm]=TMath::Abs(ToyMC::cpmCentre[rapBin-1][ptBin-1][cpmBin-1]-onia::cpmRange[cpmBin-1]);
		cpmCentreErr_high[cpm]=TMath::Abs(ToyMC::cpmCentre[rapBin-1][ptBin-1][cpmBin-1]-onia::cpmRange[cpmBin]);
		cpm++;
		}


	plotGraph = new TGraphAsymmErrors(nBinscpm,cpmCentre_,lmean3,cpmCentreErr_low,cpmCentreErr_high,lmean3err_low,lmean3err_high);
	plotGraph->SetMarkerColor(ToyMC::MarkerColor[3]);
	plotGraph->SetLineColor(ToyMC::MarkerColor[3]);
	plotGraph->SetMarkerStyle(ToyMC::MarkerStyle[nState][rapBin]);
	plotGraph->SetMarkerSize(ToyMC::MarkerSize[nState][rapBin]);
	plotGraph->Draw("P");
	sprintf(legendentry,"PX");
	plotLegend->AddEntry(plotGraph,legendentry,"ple");

	char texTex[200];
	if(rapBin==1) sprintf(texTex,"      |y| < 0.6");
	if(rapBin==2) sprintf(texTex,"0.6 < |y| < 1.2");
	if(rapBin==3) sprintf(texTex,"1.2 < |y| < 1.5");
	TLatex *text = new TLatex(onia::cpmRange[cpmBinMax]*0.22,yMin+(yMax-yMin)*0.1,texTex);
	text->SetTextSize(0.035);
	text->Draw( "same" );

	plotLegend->Draw();

	plotCanvas->SetLogx(false);
	plotCanvas->SaveAs(filename);
	plotCanvas->Close();

	delete plotCanvas;

}


void FindMPV(TH1* PosteriorDist , double& MPV , double& MPVerrorLow, double& MPVerrorHigh, int MPValgo, int nSigma){

	if(MPValgo==1){
		MPV=PosteriorDist->GetMean();
		MPVerrorLow=PosteriorDist->GetRMS();
		MPVerrorHigh=PosteriorDist->GetRMS();
	}

	if(MPValgo==2||MPValgo==3){

		int nBins = PosteriorDist->GetNbinsX();
		int maxbin_PosteriorDist = PosteriorDist->GetMaximumBin();
		double PosteriorDist_initial = PosteriorDist->GetBinCenter(maxbin_PosteriorDist);
		double err_PosteriorDist_initial=PosteriorDist->GetRMS();
		double PosteriorDist_par [3];

		TF1 *gauss;

		int nMaxFits=1;
		if(MPValgo==3) nMaxFits=20;
		for(int iFits=0;iFits<nMaxFits;iFits++){
			gauss = new TF1("f1", "gaus", PosteriorDist_initial-err_PosteriorDist_initial, PosteriorDist_initial+err_PosteriorDist_initial);
			gauss->SetParameters(PosteriorDist_initial,err_PosteriorDist_initial);
			PosteriorDist->Fit(gauss, "R");
			gauss->GetParameters(PosteriorDist_par);
			double ndof = 2*err_PosteriorDist_initial/PosteriorDist->GetBinWidth(1)-3;
			cout<<"chi2/ndf = "<<gauss->GetChisquare()/ndof<<endl;
			PosteriorDist_initial=PosteriorDist_par[1];
			err_PosteriorDist_initial=err_PosteriorDist_initial/2;
			if(gauss->GetChisquare()/ndof<5) break;
		}
		MPV=PosteriorDist_par[1];

		double OneSigmaCL;
		if(nSigma==1) OneSigmaCL=0.682689492137;
		if(nSigma==2) OneSigmaCL=0.954499736104;
		if(nSigma==3) OneSigmaCL=0.997300203937;
		double fullInt=PosteriorDist->Integral(1,nBins);
		//cout<<(1-OneSigmaCL)/2.<<endl;

		for(int i = 1; i < nBins+1; i++){
			//	cout<<i<<" "<<PosteriorDist->Integral(1,i)/fullInt<<endl;
			if(PosteriorDist->Integral(1,i)/fullInt > (1-OneSigmaCL)/2.) {MPVerrorLow=MPV-PosteriorDist->GetBinCenter(i-1); break;}
		}
		for(int i = 1; i < nBins+1; i++){
			//	cout<<i<<" "<<PosteriorDist->Integral(nBins+1-i,nBins)/fullInt<<endl;
			if(PosteriorDist->Integral(nBins+1-i,nBins)/fullInt > (1-OneSigmaCL)/2.) {MPVerrorHigh=PosteriorDist->GetBinCenter(nBins-i)-MPV; break;}
		}

	}


	return;

}

void PlotPosterior(int cpmBin, int ptBin, int rapBin, char xAxisTitle[200], char filename[200], TH1* histo, double MPV, double MPVerrLow, double MPVerrHigh){

	gStyle->SetPalette(1,0);
	gStyle->SetPadBottomMargin(0.12);
	gStyle->SetPadLeftMargin(0.13);
	gStyle->SetPadRightMargin(0.15);

	gStyle->SetTickLength(-0.02, "xyz");
	gStyle->SetLabelOffset(0.02, "x");
	gStyle->SetLabelOffset(0.02, "y");
	gStyle->SetTitleOffset(1.3, "x");
	gStyle->SetTitleOffset(1.4, "y");
	gStyle->SetTitleFillColor(kWhite);

	TLegend* plotLegend=new TLegend(0.775,0.8,0.95,0.9);
	plotLegend->SetFillColor(kWhite);
	plotLegend->SetTextFont(72);
	plotLegend->SetTextSize(0.0175);
	plotLegend->SetBorderSize(1);

	TCanvas *plotCanvas = new TCanvas("plotCanvas","plotCanvas",1000,800);

	plotCanvas->SetFillColor(kWhite);
	//	plotCanvas->SetGrid();
	plotCanvas->GetFrame()->SetFillColor(kWhite);
	plotCanvas->GetFrame()->SetBorderSize(0);
	plotCanvas->SetRightMargin(0.05) ;

	histo->SetStats(kFALSE);
	histo->SetLineColor(kBlack);
	histo->SetYTitle("Posterior Probability");
	histo->SetXTitle(xAxisTitle);
	histo->GetYaxis()->SetTitleOffset(1.5);
	histo->Draw();


	int maxbin_PosteriorDist = histo->GetMaximumBin();
	double maxval = histo->GetBinContent(maxbin_PosteriorDist);

	TLine* MeanLine = new TLine( histo->GetMean(), 0, histo->GetMean(), maxval );
	MeanLine->SetLineWidth( 1 );
	MeanLine->SetLineStyle( 1 );
	MeanLine->SetLineColor( kGreen+2 );
	MeanLine->Draw( "same" );
	TLine* MPVLine = new TLine( MPV, 0, MPV, maxval );
	MPVLine->SetLineWidth( 2.5 );
	MPVLine->SetLineStyle( 1 );
	MPVLine->SetLineColor( kRed );
	MPVLine->Draw( "same" );
	TLine* MPVerrLowLine = new TLine( MPV-MPVerrLow, 0, MPV-MPVerrLow, maxval );
	MPVerrLowLine->SetLineWidth( 2.5 );
	MPVerrLowLine->SetLineStyle( 2 );
	MPVerrLowLine->SetLineColor( kRed );
	MPVerrLowLine->Draw( "same" );
	TLine* MPVerrHighLine = new TLine( MPV+MPVerrHigh, 0, MPV+MPVerrHigh, maxval );
	MPVerrHighLine->SetLineWidth( 2.5 );
	MPVerrHighLine->SetLineStyle( 2 );
	MPVerrHighLine->SetLineColor( kRed );
	MPVerrHighLine->Draw( "same" );
	plotLegend->AddEntry(MeanLine,"Mean","l");
	plotLegend->AddEntry(MPVLine,"MPV","l");
	plotLegend->AddEntry(MPVerrLowLine,"1#sigma high/low","l");

	plotLegend->Draw();
	plotCanvas->SaveAs(filename);
	plotCanvas->Close();

	delete plotCanvas;

}

int main(int argc, char** argv) {


	int nGenerations=1;
	int polScen=1;
	int frame=1;
	Char_t *dirstruct = "ToyDirectory_Default";
	int rapBinMin=1;
	int rapBinMax=1;
	int ptBinMin=1;
	int ptBinMax=1;
	int cpmBinMin=1;
	int cpmBinMax=1;
	int MPValgo=0;
	int nState=0;
	Char_t *TreeID = "TreeID_Default";
	Char_t *realdatadir = "Default"; //Storage Directory
	int nSigma=1;

	bool RealData(false);

	char dirstructPuffer [200];			sprintf(dirstructPuffer,"%s",dirstruct);


	for( int i=0;i < argc; ++i ) {
		if(std::string(argv[i]).find("ptBinMin") != std::string::npos) {char* ptBinMinchar = argv[i]; char* ptBinMinchar2 = strtok (ptBinMinchar, "p"); ptBinMin = atof(ptBinMinchar2); cout<<"ptBinMin = "<<ptBinMin<<endl;}
		if(std::string(argv[i]).find("ptBinMax") != std::string::npos) {char* ptBinMaxchar = argv[i]; char* ptBinMaxchar2 = strtok (ptBinMaxchar, "p"); ptBinMax = atof(ptBinMaxchar2); cout<<"ptBinMax = "<<ptBinMax<<endl;}
		if(std::string(argv[i]).find("cpmBinMin") != std::string::npos) {char* cpmBinMinchar = argv[i]; char* cpmBinMinchar2 = strtok (cpmBinMinchar, "p"); cpmBinMin = atof(cpmBinMinchar2); cout<<"cpmBinMin = "<<cpmBinMin<<endl;}
		if(std::string(argv[i]).find("cpmBinMax") != std::string::npos) {char* cpmBinMaxchar = argv[i]; char* cpmBinMaxchar2 = strtok (cpmBinMaxchar, "p"); cpmBinMax = atof(cpmBinMaxchar2); cout<<"cpmBinMax = "<<cpmBinMax<<endl;}
		if(std::string(argv[i]).find("rapBinMin") != std::string::npos) {char* rapBinMinchar = argv[i]; char* rapBinMinchar2 = strtok (rapBinMinchar, "p"); rapBinMin = atof(rapBinMinchar2); cout<<"rapBinMin = "<<rapBinMin<<endl;}
		if(std::string(argv[i]).find("rapBinMax") != std::string::npos) {char* rapBinMaxchar = argv[i]; char* rapBinMaxchar2 = strtok (rapBinMaxchar, "p"); rapBinMax = atof(rapBinMaxchar2); cout<<"rapBinMax = "<<rapBinMax<<endl;}
		if(std::string(argv[i]).find("frameSig") != std::string::npos) {char* framechar = argv[i]; char* framechar2 = strtok (framechar, "p"); frame = atof(framechar2); cout<<"frame = "<<frame<<endl;}
		if(std::string(argv[i]).find("nGenerations") != std::string::npos) {char* nGenerationschar = argv[i]; char* nGenerationschar2 = strtok (nGenerationschar, "p"); nGenerations = atof(nGenerationschar2); cout<<"nGenerations = "<<nGenerations<<endl;}
		if(std::string(argv[i]).find("polScen") != std::string::npos) {char* polScenchar = argv[i]; char* polScenchar2 = strtok (polScenchar, "p"); polScen = atof(polScenchar2); cout<<"polScen = "<<polScen<<endl;}
		if(std::string(argv[i]).find("MPValgo") != std::string::npos) {char* MPValgochar = argv[i]; char* MPValgochar2 = strtok (MPValgochar, "p"); MPValgo = atof(MPValgochar2); cout<<"MPValgo = "<<MPValgo<<endl;}
		if(std::string(argv[i]).find("nState") != std::string::npos) {char* nStatechar = argv[i]; char* nStatechar2 = strtok (nStatechar, "p"); nState = atof(nStatechar2); cout<<"nState = "<<nState<<endl;}
		if(std::string(argv[i]).find("nSigma") != std::string::npos) {char* nSigmachar = argv[i]; char* nSigmachar2 = strtok (nSigmachar, "p"); nSigma = atof(nSigmachar2); cout<<"nSigma = "<<nSigma<<endl;}

		if(std::string(argv[i]).find("realdata") != std::string::npos) {RealData=true; cout<<"Plotting real data fits"<<endl;}
		if(std::string(argv[i]).find("dirstruct") != std::string::npos) {char* dirstructchar = argv[i]; char* dirstructchar2 = strtok (dirstructchar, "="); dirstruct = dirstructchar2; cout<<"dirstruct = "<<dirstruct<<endl;}
		if(std::string(argv[i]).find("TreeID") != std::string::npos) {char* TreeIDchar = argv[i]; char* TreeIDchar2 = strtok (TreeIDchar, "="); TreeID = TreeIDchar2; cout<<"TreeID = "<<TreeID<<endl;}

		if(std::string(argv[i]).find("realdatadir") != std::string::npos) {char* realdatadirchar = argv[i]; char* realdatadirchar2 = strtok (realdatadirchar, "="); realdatadir = realdatadirchar2; cout<<"realdatadir = "<<realdatadir<<endl;}
	}


	cout<<"";

	gROOT->Reset();
	gROOT->SetStyle("Plain");
	gStyle->SetOptFit(1111);
	gStyle->SetOptStat(0);
	gStyle->SetFrameBorderMode(0);

	char NatFrameName[200];
	sprintf(NatFrameName,"%s",onia::frameLabel[frame-1]);


	char tmpfilename[200];
	sprintf(tmpfilename,"%s/TGraphResults.root",dirstruct);		gSystem->Unlink(tmpfilename);



	bool HX_is_natural;
	bool PX_is_natural;
	if(frame==2) HX_is_natural=true; if(frame==3) PX_is_natural=true; //else CS is the natural frame by default

	char filename [500];

	int NrapBins = rapBinMax-rapBinMin+1;
	int NptBins = ptBinMax-ptBinMin+1;
	int NcpmBins = cpmBinMax-cpmBinMin+1;

	bool emptyBin[NrapBins][NptBins][NcpmBins];

	double pull_lth_CS_mean[NrapBins][NptBins][NcpmBins],	pull_lth_CS_meanerr[NrapBins][NptBins][NcpmBins],	pull_lth_CS_sigma[NrapBins][NptBins][NcpmBins],	pull_lth_CS_sigmaerr[NrapBins][NptBins][NcpmBins];
	double pull_lph_CS_mean[NrapBins][NptBins][NcpmBins],	pull_lph_CS_meanerr[NrapBins][NptBins][NcpmBins],	pull_lph_CS_sigma[NrapBins][NptBins][NcpmBins],	pull_lph_CS_sigmaerr[NrapBins][NptBins][NcpmBins];
	double pull_ltp_CS_mean[NrapBins][NptBins][NcpmBins],	pull_ltp_CS_meanerr[NrapBins][NptBins][NcpmBins],	pull_ltp_CS_sigma[NrapBins][NptBins][NcpmBins],	pull_ltp_CS_sigmaerr[NrapBins][NptBins][NcpmBins];
	double pull_lthstar_CS_mean[NrapBins][NptBins][NcpmBins],	pull_lthstar_CS_meanerr[NrapBins][NptBins][NcpmBins],	pull_lthstar_CS_sigma[NrapBins][NptBins][NcpmBins],	pull_lthstar_CS_sigmaerr[NrapBins][NptBins][NcpmBins];
	double pull_lphstar_CS_mean[NrapBins][NptBins][NcpmBins],	pull_lphstar_CS_meanerr[NrapBins][NptBins][NcpmBins],	pull_lphstar_CS_sigma[NrapBins][NptBins][NcpmBins],	pull_lphstar_CS_sigmaerr[NrapBins][NptBins][NcpmBins];
	double pull_ltilde_CS_mean[NrapBins][NptBins][NcpmBins],	pull_ltilde_CS_meanerr[NrapBins][NptBins][NcpmBins],	pull_ltilde_CS_sigma[NrapBins][NptBins][NcpmBins],	pull_ltilde_CS_sigmaerr[NrapBins][NptBins][NcpmBins];

	double pull_lth_HX_mean[NrapBins][NptBins][NcpmBins],	pull_lth_HX_meanerr[NrapBins][NptBins][NcpmBins],	pull_lth_HX_sigma[NrapBins][NptBins][NcpmBins],	pull_lth_HX_sigmaerr[NrapBins][NptBins][NcpmBins];
	double pull_lph_HX_mean[NrapBins][NptBins][NcpmBins],	pull_lph_HX_meanerr[NrapBins][NptBins][NcpmBins],	pull_lph_HX_sigma[NrapBins][NptBins][NcpmBins],	pull_lph_HX_sigmaerr[NrapBins][NptBins][NcpmBins];
	double pull_ltp_HX_mean[NrapBins][NptBins][NcpmBins],	pull_ltp_HX_meanerr[NrapBins][NptBins][NcpmBins],	pull_ltp_HX_sigma[NrapBins][NptBins][NcpmBins],	pull_ltp_HX_sigmaerr[NrapBins][NptBins][NcpmBins];
	double pull_lthstar_HX_mean[NrapBins][NptBins][NcpmBins],	pull_lthstar_HX_meanerr[NrapBins][NptBins][NcpmBins],	pull_lthstar_HX_sigma[NrapBins][NptBins][NcpmBins],	pull_lthstar_HX_sigmaerr[NrapBins][NptBins][NcpmBins];
	double pull_lphstar_HX_mean[NrapBins][NptBins][NcpmBins],	pull_lphstar_HX_meanerr[NrapBins][NptBins][NcpmBins],	pull_lphstar_HX_sigma[NrapBins][NptBins][NcpmBins],	pull_lphstar_HX_sigmaerr[NrapBins][NptBins][NcpmBins];
	double pull_ltilde_HX_mean[NrapBins][NptBins][NcpmBins],	pull_ltilde_HX_meanerr[NrapBins][NptBins][NcpmBins],	pull_ltilde_HX_sigma[NrapBins][NptBins][NcpmBins],	pull_ltilde_HX_sigmaerr[NrapBins][NptBins][NcpmBins];

	double pull_lth_PX_mean[NrapBins][NptBins][NcpmBins],	pull_lth_PX_meanerr[NrapBins][NptBins][NcpmBins],	pull_lth_PX_sigma[NrapBins][NptBins][NcpmBins],	pull_lth_PX_sigmaerr[NrapBins][NptBins][NcpmBins];
	double pull_lph_PX_mean[NrapBins][NptBins][NcpmBins],	pull_lph_PX_meanerr[NrapBins][NptBins][NcpmBins],	pull_lph_PX_sigma[NrapBins][NptBins][NcpmBins],	pull_lph_PX_sigmaerr[NrapBins][NptBins][NcpmBins];
	double pull_ltp_PX_mean[NrapBins][NptBins][NcpmBins],	pull_ltp_PX_meanerr[NrapBins][NptBins][NcpmBins],	pull_ltp_PX_sigma[NrapBins][NptBins][NcpmBins],	pull_ltp_PX_sigmaerr[NrapBins][NptBins][NcpmBins];
	double pull_lthstar_PX_mean[NrapBins][NptBins][NcpmBins],	pull_lthstar_PX_meanerr[NrapBins][NptBins][NcpmBins],	pull_lthstar_PX_sigma[NrapBins][NptBins][NcpmBins],	pull_lthstar_PX_sigmaerr[NrapBins][NptBins][NcpmBins];
	double pull_lphstar_PX_mean[NrapBins][NptBins][NcpmBins],	pull_lphstar_PX_meanerr[NrapBins][NptBins][NcpmBins],	pull_lphstar_PX_sigma[NrapBins][NptBins][NcpmBins],	pull_lphstar_PX_sigmaerr[NrapBins][NptBins][NcpmBins];
	double pull_ltilde_PX_mean[NrapBins][NptBins][NcpmBins],	pull_ltilde_PX_meanerr[NrapBins][NptBins][NcpmBins],	pull_ltilde_PX_sigma[NrapBins][NptBins][NcpmBins],	pull_ltilde_PX_sigmaerr[NrapBins][NptBins][NcpmBins];

	double param_lth_CS_mean[NrapBins][NptBins][NcpmBins],	param_lth_CS_meanerr[NrapBins][NptBins][NcpmBins],			param_lth_CS_sigma[NrapBins][NptBins][NcpmBins],		param_lth_CS_sigma_high[NrapBins][NptBins][NcpmBins],		param_lth_CS_sigma_low[NrapBins][NptBins][NcpmBins],		param_lth_CS_sigmaerr[NrapBins][NptBins][NcpmBins];
	double param_lph_CS_mean[NrapBins][NptBins][NcpmBins],	param_lph_CS_meanerr[NrapBins][NptBins][NcpmBins],			param_lph_CS_sigma[NrapBins][NptBins][NcpmBins],		param_lph_CS_sigma_high[NrapBins][NptBins][NcpmBins],		param_lph_CS_sigma_low[NrapBins][NptBins][NcpmBins],		param_lph_CS_sigmaerr[NrapBins][NptBins][NcpmBins];
	double param_ltp_CS_mean[NrapBins][NptBins][NcpmBins],	param_ltp_CS_meanerr[NrapBins][NptBins][NcpmBins],			param_ltp_CS_sigma[NrapBins][NptBins][NcpmBins],		param_ltp_CS_sigma_high[NrapBins][NptBins][NcpmBins],		param_ltp_CS_sigma_low[NrapBins][NptBins][NcpmBins],		param_ltp_CS_sigmaerr[NrapBins][NptBins][NcpmBins];
	double param_lthstar_CS_mean[NrapBins][NptBins][NcpmBins],	param_lthstar_CS_meanerr[NrapBins][NptBins][NcpmBins],	param_lthstar_CS_sigma[NrapBins][NptBins][NcpmBins],	param_lthstar_CS_sigma_high[NrapBins][NptBins][NcpmBins],	param_lthstar_CS_sigma_low[NrapBins][NptBins][NcpmBins],	param_lthstar_CS_sigmaerr[NrapBins][NptBins][NcpmBins];
	double param_lphstar_CS_mean[NrapBins][NptBins][NcpmBins],	param_lphstar_CS_meanerr[NrapBins][NptBins][NcpmBins],	param_lphstar_CS_sigma[NrapBins][NptBins][NcpmBins],	param_lphstar_CS_sigma_high[NrapBins][NptBins][NcpmBins],	param_lphstar_CS_sigma_low[NrapBins][NptBins][NcpmBins],	param_lphstar_CS_sigmaerr[NrapBins][NptBins][NcpmBins];
	double param_ltilde_CS_mean[NrapBins][NptBins][NcpmBins],	param_ltilde_CS_meanerr[NrapBins][NptBins][NcpmBins],		param_ltilde_CS_sigma[NrapBins][NptBins][NcpmBins],	param_ltilde_CS_sigma_high[NrapBins][NptBins][NcpmBins],	param_ltilde_CS_sigma_low[NrapBins][NptBins][NcpmBins],	param_ltilde_CS_sigmaerr[NrapBins][NptBins][NcpmBins];

	double param_lth_HX_mean[NrapBins][NptBins][NcpmBins],	param_lth_HX_meanerr[NrapBins][NptBins][NcpmBins],			param_lth_HX_sigma[NrapBins][NptBins][NcpmBins],		param_lth_HX_sigma_high[NrapBins][NptBins][NcpmBins],		param_lth_HX_sigma_low[NrapBins][NptBins][NcpmBins],		param_lth_HX_sigmaerr[NrapBins][NptBins][NcpmBins];
	double param_lph_HX_mean[NrapBins][NptBins][NcpmBins],	param_lph_HX_meanerr[NrapBins][NptBins][NcpmBins],			param_lph_HX_sigma[NrapBins][NptBins][NcpmBins],		param_lph_HX_sigma_high[NrapBins][NptBins][NcpmBins],		param_lph_HX_sigma_low[NrapBins][NptBins][NcpmBins],		param_lph_HX_sigmaerr[NrapBins][NptBins][NcpmBins];
	double param_ltp_HX_mean[NrapBins][NptBins][NcpmBins],	param_ltp_HX_meanerr[NrapBins][NptBins][NcpmBins],			param_ltp_HX_sigma[NrapBins][NptBins][NcpmBins],		param_ltp_HX_sigma_high[NrapBins][NptBins][NcpmBins],		param_ltp_HX_sigma_low[NrapBins][NptBins][NcpmBins],		param_ltp_HX_sigmaerr[NrapBins][NptBins][NcpmBins];
	double param_lthstar_HX_mean[NrapBins][NptBins][NcpmBins],	param_lthstar_HX_meanerr[NrapBins][NptBins][NcpmBins],	param_lthstar_HX_sigma[NrapBins][NptBins][NcpmBins],	param_lthstar_HX_sigma_high[NrapBins][NptBins][NcpmBins],	param_lthstar_HX_sigma_low[NrapBins][NptBins][NcpmBins],	param_lthstar_HX_sigmaerr[NrapBins][NptBins][NcpmBins];
	double param_lphstar_HX_mean[NrapBins][NptBins][NcpmBins],	param_lphstar_HX_meanerr[NrapBins][NptBins][NcpmBins],	param_lphstar_HX_sigma[NrapBins][NptBins][NcpmBins],	param_lphstar_HX_sigma_high[NrapBins][NptBins][NcpmBins],	param_lphstar_HX_sigma_low[NrapBins][NptBins][NcpmBins],	param_lphstar_HX_sigmaerr[NrapBins][NptBins][NcpmBins];
	double param_ltilde_HX_mean[NrapBins][NptBins][NcpmBins],	param_ltilde_HX_meanerr[NrapBins][NptBins][NcpmBins],		param_ltilde_HX_sigma[NrapBins][NptBins][NcpmBins],	param_ltilde_HX_sigma_high[NrapBins][NptBins][NcpmBins],	param_ltilde_HX_sigma_low[NrapBins][NptBins][NcpmBins],	param_ltilde_HX_sigmaerr[NrapBins][NptBins][NcpmBins];

	double param_lth_PX_mean[NrapBins][NptBins][NcpmBins],	param_lth_PX_meanerr[NrapBins][NptBins][NcpmBins],			param_lth_PX_sigma[NrapBins][NptBins][NcpmBins],		param_lth_PX_sigma_high[NrapBins][NptBins][NcpmBins],		param_lth_PX_sigma_low[NrapBins][NptBins][NcpmBins],		param_lth_PX_sigmaerr[NrapBins][NptBins][NcpmBins];
	double param_lph_PX_mean[NrapBins][NptBins][NcpmBins],	param_lph_PX_meanerr[NrapBins][NptBins][NcpmBins],			param_lph_PX_sigma[NrapBins][NptBins][NcpmBins],		param_lph_PX_sigma_high[NrapBins][NptBins][NcpmBins],		param_lph_PX_sigma_low[NrapBins][NptBins][NcpmBins],		param_lph_PX_sigmaerr[NrapBins][NptBins][NcpmBins];
	double param_ltp_PX_mean[NrapBins][NptBins][NcpmBins],	param_ltp_PX_meanerr[NrapBins][NptBins][NcpmBins],			param_ltp_PX_sigma[NrapBins][NptBins][NcpmBins],		param_ltp_PX_sigma_high[NrapBins][NptBins][NcpmBins],		param_ltp_PX_sigma_low[NrapBins][NptBins][NcpmBins],		param_ltp_PX_sigmaerr[NrapBins][NptBins][NcpmBins];
	double param_lthstar_PX_mean[NrapBins][NptBins][NcpmBins],	param_lthstar_PX_meanerr[NrapBins][NptBins][NcpmBins],	param_lthstar_PX_sigma[NrapBins][NptBins][NcpmBins],	param_lthstar_PX_sigma_high[NrapBins][NptBins][NcpmBins],	param_lthstar_PX_sigma_low[NrapBins][NptBins][NcpmBins],	param_lthstar_PX_sigmaerr[NrapBins][NptBins][NcpmBins];
	double param_lphstar_PX_mean[NrapBins][NptBins][NcpmBins],	param_lphstar_PX_meanerr[NrapBins][NptBins][NcpmBins],	param_lphstar_PX_sigma[NrapBins][NptBins][NcpmBins],	param_lphstar_PX_sigma_high[NrapBins][NptBins][NcpmBins],	param_lphstar_PX_sigma_low[NrapBins][NptBins][NcpmBins],	param_lphstar_PX_sigmaerr[NrapBins][NptBins][NcpmBins];
	double param_ltilde_PX_mean[NrapBins][NptBins][NcpmBins],	param_ltilde_PX_meanerr[NrapBins][NptBins][NcpmBins],		param_ltilde_PX_sigma[NrapBins][NptBins][NcpmBins],	param_ltilde_PX_sigma_high[NrapBins][NptBins][NcpmBins],	param_ltilde_PX_sigma_low[NrapBins][NptBins][NcpmBins],	param_ltilde_PX_sigmaerr[NrapBins][NptBins][NcpmBins];


	double lambda_theta_injected_CS[NrapBins][NptBins][NcpmBins];
	double lambda_phi_injected_CS[NrapBins][NptBins][NcpmBins];
	double lambda_thetaphi_injected_CS[NrapBins][NptBins][NcpmBins];
	double lambda_thetastar_injected_CS[NrapBins][NptBins][NcpmBins];
	double lambda_phistar_injected_CS[NrapBins][NptBins][NcpmBins];
	double lambda_tilde_injected_CS[NrapBins][NptBins][NcpmBins];

	double lambda_theta_injected_HX[NrapBins][NptBins][NcpmBins];
	double lambda_phi_injected_HX[NrapBins][NptBins][NcpmBins];
	double lambda_thetaphi_injected_HX[NrapBins][NptBins][NcpmBins];
	double lambda_thetastar_injected_HX[NrapBins][NptBins][NcpmBins];
	double lambda_phistar_injected_HX[NrapBins][NptBins][NcpmBins];
	double lambda_tilde_injected_HX[NrapBins][NptBins][NcpmBins];

	double lambda_theta_injected_PX[NrapBins][NptBins][NcpmBins];
	double lambda_phi_injected_PX[NrapBins][NptBins][NcpmBins];
	double lambda_thetaphi_injected_PX[NrapBins][NptBins][NcpmBins];
	double lambda_thetastar_injected_PX[NrapBins][NptBins][NcpmBins];
	double lambda_phistar_injected_PX[NrapBins][NptBins][NcpmBins];
	double lambda_tilde_injected_PX[NrapBins][NptBins][NcpmBins];

	int nGen_[NrapBins][NptBins][NcpmBins];

	double l_min;
	double l_max;
	double l_step_1D;

	double l_min_pull;
	double l_max_pull;
	double l_step_1D_pull;

	double delta_l;

	int rap=0;
	for(int iRap = rapBinMin; iRap < rapBinMax+1; iRap++){
		rap++;int pt=0;
		for(int iPt = ptBinMin; iPt < ptBinMax+1; iPt++){
			pt++;
		int cpm=0;
		for(int icpm = cpmBinMin; icpm < cpmBinMax+1; icpm++){
			cpm++;	

			///////////////////// Extraction of mean injected parameters from GenResults.root //////////////////////////////////////////////////

			////2013-6-21
			//MPValgo = 3;
			//if(iPt==12) MPValgo = 1;

			l_min = -1.4;
			l_max =  1.4;
			l_step_1D = 0.02;

			l_min_pull = -5;
			l_max_pull =  5;
			l_step_1D_pull = 0.5;

			TH1D* inj_lth_CS = new TH1D( "inj_lth_CS", "", int((l_max-l_min)/l_step_1D), l_min, l_max );
			TH1D* inj_lph_CS = new TH1D( "inj_lph_CS", "", int((l_max-l_min)/l_step_1D), l_min, l_max );
			TH1D* inj_ltp_CS = new TH1D( "inj_ltp_CS", "", int((l_max-l_min)/l_step_1D), l_min, l_max );

			TH1D* inj_lth_HX = new TH1D( "inj_lth_HX", "", int((l_max-l_min)/l_step_1D), l_min, l_max );
			TH1D* inj_lph_HX = new TH1D( "inj_lph_HX", "", int((l_max-l_min)/l_step_1D), l_min, l_max );
			TH1D* inj_ltp_HX = new TH1D( "inj_ltp_HX", "", int((l_max-l_min)/l_step_1D), l_min, l_max );

			TH1D* inj_lth_PX = new TH1D( "inj_lth_PX", "", int((l_max-l_min)/l_step_1D), l_min, l_max );
			TH1D* inj_lph_PX = new TH1D( "inj_lph_PX", "", int((l_max-l_min)/l_step_1D), l_min, l_max );
			TH1D* inj_ltp_PX = new TH1D( "inj_ltp_PX", "", int((l_max-l_min)/l_step_1D), l_min, l_max );


			for(int iGen=1;iGen<nGenerations+1;iGen++){

				sprintf(filename,"%s/rap%d_pT%d_cpm%d/Generation%d/GenResults.root",dirstruct,iRap,iPt,icpm,iGen);
				TFile* GenResultFile = new TFile(filename);

				if(GenResultFile->Get("GenResults") == NULL) {GenResultFile->Close(); continue;}

				TTree* GenResults = (TTree*) GenResultFile->Get("GenResults");


				double lamth_CS; double lamph_CS; double lamtp_CS;
				double lamth_HX; double lamph_HX; double lamtp_HX;
				double lamth_PX; double lamph_PX; double lamtp_PX;

				TBranch        *b_inj_lth_CS;				GenResults->SetBranchAddress("lthCS", &lamth_CS, &b_inj_lth_CS);
				TBranch        *b_inj_lph_CS;				GenResults->SetBranchAddress("lphCS", &lamph_CS, &b_inj_lph_CS);
				TBranch        *b_inj_ltp_CS;				GenResults->SetBranchAddress("ltpCS", &lamtp_CS, &b_inj_ltp_CS);
				TBranch        *b_inj_lth_HX;				GenResults->SetBranchAddress("lthHX", &lamth_HX, &b_inj_lth_HX);
				TBranch        *b_inj_lph_HX;				GenResults->SetBranchAddress("lphHX", &lamph_HX, &b_inj_lph_HX);
				TBranch        *b_inj_ltp_HX;				GenResults->SetBranchAddress("ltpHX", &lamtp_HX, &b_inj_ltp_HX);
				TBranch        *b_inj_lth_PX;				GenResults->SetBranchAddress("lthPX", &lamth_PX, &b_inj_lth_PX);
				TBranch        *b_inj_lph_PX;				GenResults->SetBranchAddress("lphPX", &lamph_PX, &b_inj_lph_PX);
				TBranch        *b_inj_ltp_PX;				GenResults->SetBranchAddress("ltpPX", &lamtp_PX, &b_inj_ltp_PX);

				GenResults->GetEvent( 0 );


				inj_lth_CS->Fill( lamth_CS );
				inj_lph_CS->Fill( lamph_CS );
				inj_ltp_CS->Fill( lamtp_CS );

				inj_lth_HX->Fill( lamth_HX );
				inj_lph_HX->Fill( lamph_HX );
				inj_ltp_HX->Fill( lamtp_HX );

				inj_lth_PX->Fill( lamth_PX );
				inj_lph_PX->Fill( lamph_PX );
				inj_ltp_PX->Fill( lamtp_PX );


				GenResultFile->Close();

			}

			lambda_theta_injected_CS[rap-1][pt-1][cpm-1]=inj_lth_CS->GetMean();
			lambda_phi_injected_CS[rap-1][pt-1][cpm-1]=inj_lph_CS->GetMean();
			lambda_thetaphi_injected_CS[rap-1][pt-1][cpm-1]=inj_ltp_CS->GetMean();
			calcLambdastar(lambda_thetastar_injected_CS[rap-1][pt-1][cpm-1],lambda_phistar_injected_CS[rap-1][pt-1][cpm-1],lambda_theta_injected_CS[rap-1][pt-1][cpm-1],lambda_phi_injected_CS[rap-1][pt-1][cpm-1],lambda_thetaphi_injected_CS[rap-1][pt-1][cpm-1]);
			lambda_tilde_injected_CS[rap-1][pt-1][cpm-1]=(lambda_theta_injected_CS[rap-1][pt-1][cpm-1]+3.*lambda_phi_injected_CS[rap-1][pt-1][cpm-1])/(1-lambda_phi_injected_CS[rap-1][pt-1][cpm-1]);

			lambda_theta_injected_HX[rap-1][pt-1][cpm-1]=inj_lth_HX->GetMean();
			lambda_phi_injected_HX[rap-1][pt-1][cpm-1]=inj_lph_HX->GetMean();
			lambda_thetaphi_injected_HX[rap-1][pt-1][cpm-1]=inj_ltp_HX->GetMean();
			calcLambdastar(lambda_thetastar_injected_HX[rap-1][pt-1][cpm-1],lambda_phistar_injected_HX[rap-1][pt-1][cpm-1],lambda_theta_injected_HX[rap-1][pt-1][cpm-1],lambda_phi_injected_HX[rap-1][pt-1][cpm-1],lambda_thetaphi_injected_HX[rap-1][pt-1][cpm-1]);
			lambda_tilde_injected_HX[rap-1][pt-1][cpm-1]=(lambda_theta_injected_HX[rap-1][pt-1][cpm-1]+3.*lambda_phi_injected_HX[rap-1][pt-1][cpm-1])/(1-lambda_phi_injected_HX[rap-1][pt-1][cpm-1]);

			lambda_theta_injected_PX[rap-1][pt-1][cpm-1]=inj_lth_PX->GetMean();
			lambda_phi_injected_PX[rap-1][pt-1][cpm-1]=inj_lph_PX->GetMean();
			lambda_thetaphi_injected_PX[rap-1][pt-1][cpm-1]=inj_ltp_PX->GetMean();
			calcLambdastar(lambda_thetastar_injected_PX[rap-1][pt-1][cpm-1],lambda_phistar_injected_PX[rap-1][pt-1][cpm-1],lambda_theta_injected_PX[rap-1][pt-1][cpm-1],lambda_phi_injected_PX[rap-1][pt-1][cpm-1],lambda_thetaphi_injected_PX[rap-1][pt-1][cpm-1]);
			lambda_tilde_injected_PX[rap-1][pt-1][cpm-1]=(lambda_theta_injected_PX[rap-1][pt-1][cpm-1]+3.*lambda_phi_injected_PX[rap-1][pt-1][cpm-1])/(1-lambda_phi_injected_PX[rap-1][pt-1][cpm-1]);

			///////////////////// End of Extraction ////////////////////////////////////////////////////////////////////////////////////////////



			l_min = -1.4;
			l_max =  1.4;

			l_min_pull = -10;
			l_max_pull =  10;
			l_step_1D_pull = 2*l_max_pull/50.;

			delta_l = 1.2;
			l_step_1D = 2*delta_l/50.;

			double l_min_lth_CS = lambda_theta_injected_CS[rap-1][pt-1][cpm-1] - delta_l;
			double l_max_lth_CS = lambda_theta_injected_CS[rap-1][pt-1][cpm-1] + delta_l;
			double l_min_lph_CS = lambda_phi_injected_CS[rap-1][pt-1][cpm-1] - delta_l;
			double l_max_lph_CS = lambda_phi_injected_CS[rap-1][pt-1][cpm-1] + delta_l;
			double l_min_ltp_CS = lambda_thetaphi_injected_CS[rap-1][pt-1][cpm-1] - delta_l;
			double l_max_ltp_CS = lambda_thetaphi_injected_CS[rap-1][pt-1][cpm-1] + delta_l;
			double l_min_lthstar_CS = lambda_thetastar_injected_CS[rap-1][pt-1][cpm-1] - delta_l;
			double l_max_lthstar_CS = lambda_thetastar_injected_CS[rap-1][pt-1][cpm-1] + delta_l;
			double l_min_lphstar_CS = lambda_phistar_injected_CS[rap-1][pt-1][cpm-1] - delta_l;
			double l_max_lphstar_CS = lambda_phistar_injected_CS[rap-1][pt-1][cpm-1] + delta_l;
			double l_min_ltilde_CS = lambda_tilde_injected_CS[rap-1][pt-1][cpm-1] - delta_l;
			double l_max_ltilde_CS = lambda_tilde_injected_CS[rap-1][pt-1][cpm-1] + delta_l;

			double l_min_lth_HX = lambda_theta_injected_HX[rap-1][pt-1][cpm-1] - delta_l;
			double l_max_lth_HX = lambda_theta_injected_HX[rap-1][pt-1][cpm-1] + delta_l;
			double l_min_lph_HX = lambda_phi_injected_HX[rap-1][pt-1][cpm-1] - delta_l;
			double l_max_lph_HX = lambda_phi_injected_HX[rap-1][pt-1][cpm-1] + delta_l;
			double l_min_ltp_HX = lambda_thetaphi_injected_HX[rap-1][pt-1][cpm-1] - delta_l;
			double l_max_ltp_HX = lambda_thetaphi_injected_HX[rap-1][pt-1][cpm-1] + delta_l;
			double l_min_lthstar_HX = lambda_thetastar_injected_HX[rap-1][pt-1][cpm-1] - delta_l;
			double l_max_lthstar_HX = lambda_thetastar_injected_HX[rap-1][pt-1][cpm-1] + delta_l;
			double l_min_lphstar_HX = lambda_phistar_injected_HX[rap-1][pt-1][cpm-1] - delta_l;
			double l_max_lphstar_HX = lambda_phistar_injected_HX[rap-1][pt-1][cpm-1] + delta_l;
			double l_min_ltilde_HX = lambda_tilde_injected_HX[rap-1][pt-1][cpm-1] - delta_l;
			double l_max_ltilde_HX = lambda_tilde_injected_HX[rap-1][pt-1][cpm-1] + delta_l;

			double l_min_lth_PX = lambda_theta_injected_PX[rap-1][pt-1][cpm-1] - delta_l;
			double l_max_lth_PX = lambda_theta_injected_PX[rap-1][pt-1][cpm-1] + delta_l;
			double l_min_lph_PX = lambda_phi_injected_PX[rap-1][pt-1][cpm-1] - delta_l;
			double l_max_lph_PX = lambda_phi_injected_PX[rap-1][pt-1][cpm-1] + delta_l;
			double l_min_ltp_PX = lambda_thetaphi_injected_PX[rap-1][pt-1][cpm-1] - delta_l;
			double l_max_ltp_PX = lambda_thetaphi_injected_PX[rap-1][pt-1][cpm-1] + delta_l;
			double l_min_lthstar_PX = lambda_thetastar_injected_PX[rap-1][pt-1][cpm-1] - delta_l;
			double l_max_lthstar_PX = lambda_thetastar_injected_PX[rap-1][pt-1][cpm-1] + delta_l;
			double l_min_lphstar_PX = lambda_phistar_injected_PX[rap-1][pt-1][cpm-1] - delta_l;
			double l_max_lphstar_PX = lambda_phistar_injected_PX[rap-1][pt-1][cpm-1] + delta_l;
			double l_min_ltilde_PX = lambda_tilde_injected_PX[rap-1][pt-1][cpm-1] - delta_l;
			double l_max_ltilde_PX = lambda_tilde_injected_PX[rap-1][pt-1][cpm-1] + delta_l;

			TH1D* param_lth_CS = new TH1D( "param_lth_CS", "", int((l_max_lth_CS-l_min_lth_CS)/l_step_1D), l_min_lth_CS, l_max_lth_CS );
			TH1D* param_lph_CS = new TH1D( "param_lph_CS", "", int((l_max_lph_CS-l_min_lph_CS)/l_step_1D), l_min_lph_CS, l_max_lph_CS );
			TH1D* param_ltp_CS = new TH1D( "param_ltp_CS", "", int((l_max_ltp_CS-l_min_ltp_CS)/l_step_1D), l_min_ltp_CS, l_max_ltp_CS );
			TH1D* param_ltilde_CS = new TH1D( "param_ltilde_CS", "", int((l_max_ltilde_CS-l_min_ltilde_CS)/l_step_1D), l_min_ltilde_CS, l_max_ltilde_CS );
			TH1D* param_lphstar_CS = new TH1D( "param_lphstar_CS", "", int((l_max_lphstar_CS-l_min_lphstar_CS)/l_step_1D), l_min_lphstar_CS, l_max_lphstar_CS );
			TH1D* param_lthstar_CS = new TH1D( "param_lthstar_CS", "", int((l_max_lthstar_CS-l_min_lthstar_CS)/l_step_1D), l_min_lthstar_CS, l_max_lthstar_CS );

			TH1D* param_lth_HX = new TH1D( "param_lth_HX", "", int((l_max_lth_HX-l_min_lth_HX)/l_step_1D), l_min_lth_HX, l_max_lth_HX );
			TH1D* param_lph_HX = new TH1D( "param_lph_HX", "", int((l_max_lph_HX-l_min_lph_HX)/l_step_1D), l_min_lph_HX, l_max_lph_HX );
			TH1D* param_ltp_HX = new TH1D( "param_ltp_HX", "", int((l_max_ltp_HX-l_min_ltp_HX)/l_step_1D), l_min_ltp_HX, l_max_ltp_HX );
			TH1D* param_ltilde_HX = new TH1D( "param_ltilde_HX", "", int((l_max_ltilde_HX-l_min_ltilde_HX)/l_step_1D), l_min_ltilde_HX, l_max_ltilde_HX );
			TH1D* param_lphstar_HX = new TH1D( "param_lphstar_HX", "", int((l_max_lphstar_HX-l_min_lphstar_HX)/l_step_1D), l_min_lphstar_HX, l_max_lphstar_HX );
			TH1D* param_lthstar_HX = new TH1D( "param_lthstar_HX", "", int((l_max_lthstar_HX-l_min_lthstar_HX)/l_step_1D), l_min_lthstar_HX, l_max_lthstar_HX );

			TH1D* param_lth_PX = new TH1D( "param_lth_PX", "", int((l_max_lth_PX-l_min_lth_PX)/l_step_1D), l_min_lth_PX, l_max_lth_PX );
			TH1D* param_lph_PX = new TH1D( "param_lph_PX", "", int((l_max_lph_PX-l_min_lph_PX)/l_step_1D), l_min_lph_PX, l_max_lph_PX );
			TH1D* param_ltp_PX = new TH1D( "param_ltp_PX", "", int((l_max_ltp_PX-l_min_ltp_PX)/l_step_1D), l_min_ltp_PX, l_max_ltp_PX );
			TH1D* param_ltilde_PX = new TH1D( "param_ltilde_PX", "", int((l_max_ltilde_PX-l_min_ltilde_PX)/l_step_1D), l_min_ltilde_PX, l_max_ltilde_PX );
			TH1D* param_lphstar_PX = new TH1D( "param_lphstar_PX", "", int((l_max_lphstar_PX-l_min_lphstar_PX)/l_step_1D), l_min_lphstar_PX, l_max_lphstar_PX );
			TH1D* param_lthstar_PX = new TH1D( "param_lthstar_PX", "", int((l_max_lthstar_PX-l_min_lthstar_PX)/l_step_1D), l_min_lthstar_PX, l_max_lthstar_PX );

			TH1D* pull_lth_CS = new TH1D( "pull_lth_CS", "", int((l_max_pull-l_min_pull)/l_step_1D_pull), l_min_pull, l_max_pull );
			TH1D* pull_lph_CS = new TH1D( "pull_lph_CS", "", int((l_max_pull-l_min_pull)/l_step_1D_pull), l_min_pull, l_max_pull );
			TH1D* pull_ltp_CS = new TH1D( "pull_ltp_CS", "", int((l_max_pull-l_min_pull)/l_step_1D_pull), l_min_pull, l_max_pull );
			TH1D* pull_ltilde_CS = new TH1D( "pull_ltilde_CS", "", int((l_max_pull-l_min_pull)/l_step_1D_pull), l_min_pull, l_max_pull );
			TH1D* pull_lphstar_CS = new TH1D( "pull_lphstar_CS", "", int((l_max_pull-l_min_pull)/l_step_1D_pull), l_min_pull, l_max_pull );
			TH1D* pull_lthstar_CS = new TH1D( "pull_lthstar_CS", "", int((l_max_pull-l_min_pull)/l_step_1D_pull), l_min_pull, l_max_pull );

			TH1D* pull_lth_HX = new TH1D( "pull_lth_HX", "", int((l_max_pull-l_min_pull)/l_step_1D_pull), l_min_pull, l_max_pull );
			TH1D* pull_lph_HX = new TH1D( "pull_lph_HX", "", int((l_max_pull-l_min_pull)/l_step_1D_pull), l_min_pull, l_max_pull );
			TH1D* pull_ltp_HX = new TH1D( "pull_ltp_HX", "", int((l_max_pull-l_min_pull)/l_step_1D_pull), l_min_pull, l_max_pull );
			TH1D* pull_ltilde_HX = new TH1D( "pull_ltilde_HX", "", int((l_max_pull-l_min_pull)/l_step_1D_pull), l_min_pull, l_max_pull );
			TH1D* pull_lphstar_HX = new TH1D( "pull_lphstar_HX", "", int((l_max_pull-l_min_pull)/l_step_1D_pull), l_min_pull, l_max_pull );
			TH1D* pull_lthstar_HX = new TH1D( "pull_lthstar_HX", "", int((l_max_pull-l_min_pull)/l_step_1D_pull), l_min_pull, l_max_pull );

			TH1D* pull_lth_PX = new TH1D( "pull_lth_PX", "", int((l_max_pull-l_min_pull)/l_step_1D_pull), l_min_pull, l_max_pull );
			TH1D* pull_lph_PX = new TH1D( "pull_lph_PX", "", int((l_max_pull-l_min_pull)/l_step_1D_pull), l_min_pull, l_max_pull );
			TH1D* pull_ltp_PX = new TH1D( "pull_ltp_PX", "", int((l_max_pull-l_min_pull)/l_step_1D_pull), l_min_pull, l_max_pull );
			TH1D* pull_ltilde_PX = new TH1D( "pull_ltilde_PX", "", int((l_max_pull-l_min_pull)/l_step_1D_pull), l_min_pull, l_max_pull );
			TH1D* pull_lphstar_PX = new TH1D( "pull_lphstar_PX", "", int((l_max_pull-l_min_pull)/l_step_1D_pull), l_min_pull, l_max_pull );
			TH1D* pull_lthstar_PX = new TH1D( "pull_lthstar_PX", "", int((l_max_pull-l_min_pull)/l_step_1D_pull), l_min_pull, l_max_pull );

			TH1D* err_low_param_lth_CS = new TH1D( "err_low_param_lth_CS", "", 100, 0, 2 );
			TH1D* err_high_param_lth_CS = new TH1D( "err_high_param_lth_CS", "", 100, 0, 2 );
			TH1D* err_param_lth_CS = new TH1D( "err_param_lth_CS", "", 100, 0, 2 );

			TH1D* err_low_param_lph_CS = new TH1D( "err_low_param_lph_CS", "", 100, 0, 2 );
			TH1D* err_high_param_lph_CS = new TH1D( "err_high_param_lph_CS", "", 100, 0, 2 );
			TH1D* err_param_lph_CS = new TH1D( "err_param_lph_CS", "", 100, 0, 2 );

			TH1D* err_low_param_ltp_CS = new TH1D( "err_low_param_ltp_CS", "", 100, 0, 2 );
			TH1D* err_high_param_ltp_CS = new TH1D( "err_high_param_ltp_CS", "", 100, 0, 2 );
			TH1D* err_param_ltp_CS = new TH1D( "err_param_ltp_CS", "", 100, 0, 2 );

			TH1D* err_low_param_ltilde_CS = new TH1D( "err_low_param_ltilde_CS", "", 100, 0, 2 );
			TH1D* err_high_param_ltilde_CS = new TH1D( "err_high_param_ltilde_CS", "", 100, 0, 2 );
			TH1D* err_param_ltilde_CS = new TH1D( "err_param_ltilde_CS", "", 100, 0, 2 );

			TH1D* err_low_param_lphstar_CS = new TH1D( "err_low_param_lphstar_CS", "", 100, 0, 2 );
			TH1D* err_high_param_lphstar_CS = new TH1D( "err_high_param_lphstar_CS", "", 100, 0, 2 );
			TH1D* err_param_lphstar_CS = new TH1D( "err_param_lphstar_CS", "", 100, 0, 2 );

			TH1D* err_low_param_lthstar_CS = new TH1D( "err_low_param_lthstar_CS", "", 100, 0, 2 );
			TH1D* err_high_param_lthstar_CS = new TH1D( "err_high_param_lthstar_CS", "", 100, 0, 2 );
			TH1D* err_param_lthstar_CS = new TH1D( "err_param_lthstar_CS", "", 100, 0, 2 );


			TH1D* err_low_param_lth_HX = new TH1D( "err_low_param_lth_HX", "", 100, 0, 2 );
			TH1D* err_high_param_lth_HX = new TH1D( "err_high_param_lth_HX", "", 100, 0, 2 );
			TH1D* err_param_lth_HX = new TH1D( "err_param_lth_HX", "", 100, 0, 2 );

			TH1D* err_low_param_lph_HX = new TH1D( "err_low_param_lph_HX", "", 100, 0, 2 );
			TH1D* err_high_param_lph_HX = new TH1D( "err_high_param_lph_HX", "", 100, 0, 2 );
			TH1D* err_param_lph_HX = new TH1D( "err_param_lph_HX", "", 100, 0, 2 );

			TH1D* err_low_param_ltp_HX = new TH1D( "err_low_param_ltp_HX", "", 100, 0, 2 );
			TH1D* err_high_param_ltp_HX = new TH1D( "err_high_param_ltp_HX", "", 100, 0, 2 );
			TH1D* err_param_ltp_HX = new TH1D( "err_param_ltp_HX", "", 100, 0, 2 );

			TH1D* err_low_param_ltilde_HX = new TH1D( "err_low_param_ltilde_HX", "", 100, 0, 2 );
			TH1D* err_high_param_ltilde_HX = new TH1D( "err_high_param_ltilde_HX", "", 100, 0, 2 );
			TH1D* err_param_ltilde_HX = new TH1D( "err_param_ltilde_HX", "", 100, 0, 2 );

			TH1D* err_low_param_lphstar_HX = new TH1D( "err_low_param_lphstar_HX", "", 100, 0, 2 );
			TH1D* err_high_param_lphstar_HX = new TH1D( "err_high_param_lphstar_HX", "", 100, 0, 2 );
			TH1D* err_param_lphstar_HX = new TH1D( "err_param_lphstar_HX", "", 100, 0, 2 );

			TH1D* err_low_param_lthstar_HX = new TH1D( "err_low_param_lthstar_HX", "", 100, 0, 2 );
			TH1D* err_high_param_lthstar_HX = new TH1D( "err_high_param_lthstar_HX", "", 100, 0, 2 );
			TH1D* err_param_lthstar_HX = new TH1D( "err_param_lthstar_HX", "", 100, 0, 2 );


			TH1D* err_low_param_lth_PX = new TH1D( "err_low_param_lth_PX", "", 100, 0, 2 );
			TH1D* err_high_param_lth_PX = new TH1D( "err_high_param_lth_PX", "", 100, 0, 2 );
			TH1D* err_param_lth_PX = new TH1D( "err_param_lth_PX", "", 100, 0, 2 );

			TH1D* err_low_param_lph_PX = new TH1D( "err_low_param_lph_PX", "", 100, 0, 2 );
			TH1D* err_high_param_lph_PX = new TH1D( "err_high_param_lph_PX", "", 100, 0, 2 );
			TH1D* err_param_lph_PX = new TH1D( "err_param_lph_PX", "", 100, 0, 2 );

			TH1D* err_low_param_ltp_PX = new TH1D( "err_low_param_ltp_PX", "", 100, 0, 2 );
			TH1D* err_high_param_ltp_PX = new TH1D( "err_high_param_ltp_PX", "", 100, 0, 2 );
			TH1D* err_param_ltp_PX = new TH1D( "err_param_ltp_PX", "", 100, 0, 2 );

			TH1D* err_low_param_ltilde_PX = new TH1D( "err_low_param_ltilde_PX", "", 100, 0, 2 );
			TH1D* err_high_param_ltilde_PX = new TH1D( "err_high_param_ltilde_PX", "", 100, 0, 2 );
			TH1D* err_param_ltilde_PX = new TH1D( "err_param_ltilde_PX", "", 100, 0, 2 );

			TH1D* err_low_param_lphstar_PX = new TH1D( "err_low_param_lphstar_PX", "", 100, 0, 2 );
			TH1D* err_high_param_lphstar_PX = new TH1D( "err_high_param_lphstar_PX", "", 100, 0, 2 );
			TH1D* err_param_lphstar_PX = new TH1D( "err_param_lphstar_PX", "", 100, 0, 2 );

			TH1D* err_low_param_lthstar_PX = new TH1D( "err_low_param_lthstar_PX", "", 100, 0, 2 );
			TH1D* err_high_param_lthstar_PX = new TH1D( "err_high_param_lthstar_PX", "", 100, 0, 2 );
			TH1D* err_param_lthstar_PX = new TH1D( "err_param_lthstar_PX", "", 100, 0, 2 );

			double err_lth_CS; double err_lph_CS; double err_ltp_CS; double err_lthstar_CS; double err_lphstar_CS; double err_ltilde_CS;
			double err_lth_HX; double err_lph_HX; double err_ltp_HX; double err_lthstar_HX; double err_lphstar_HX; double err_ltilde_HX;
			double err_lth_PX; double err_lph_PX; double err_ltp_PX; double err_lthstar_PX; double err_lphstar_PX; double err_ltilde_PX;
			double err_high_lth_CS; double err_high_lph_CS; double err_high_ltp_CS; double err_high_lthstar_CS; double err_high_lphstar_CS; double err_high_ltilde_CS;
			double err_high_lth_HX; double err_high_lph_HX; double err_high_ltp_HX; double err_high_lthstar_HX; double err_high_lphstar_HX; double err_high_ltilde_HX;
			double err_high_lth_PX; double err_high_lph_PX; double err_high_ltp_PX; double err_high_lthstar_PX; double err_high_lphstar_PX; double err_high_ltilde_PX;
			double err_low_lth_CS; double err_low_lph_CS; double err_low_ltp_CS; double err_low_lthstar_CS; double err_low_lphstar_CS; double err_low_ltilde_CS;
			double err_low_lth_HX; double err_low_lph_HX; double err_low_ltp_HX; double err_low_lthstar_HX; double err_low_lphstar_HX; double err_low_ltilde_HX;
			double err_low_lth_PX; double err_low_lph_PX; double err_low_ltp_PX; double err_low_lthstar_PX; double err_low_lphstar_PX; double err_low_ltilde_PX;
			double lth_CS; double lph_CS; double ltp_CS; double lthstar_CS; double lphstar_CS; double ltilde_CS;
			double lth_HX; double lph_HX; double ltp_HX; double lthstar_HX; double lphstar_HX; double ltilde_HX;
			double lth_PX; double lph_PX; double ltp_PX; double lthstar_PX; double lphstar_PX; double ltilde_PX;
			int positivity_CS; int positivity_HX; int positivity_PX;

			nGen_[rap-1][pt-1][cpm-1]=0;

			for(int iGen=1;iGen<nGenerations+1;iGen++){

				sprintf(filename,"%s/rap%d_pT%d_cpm%d/Generation%d/results.root",dirstruct,iRap,iPt,icpm,iGen);
				if(RealData) sprintf(filename,"%s/results_%s_rap%d_pT%d_cpm%d.root",dirstruct,TreeID,iRap,iPt,icpm);

				TFile* results = new TFile(filename);

				if(results->Get("Results") == NULL) {results->Close(); continue;}

				TTree* Results = (TTree*) results->Get("Results");

				nGen_[rap-1][pt-1][cpm-1]++;

				cout<<"";



				TTree* results_CS = (TTree*) results->Get("lambdaCS");
				TBranch        *b_lth_CS;				results_CS->SetBranchAddress("lth", &lth_CS, &b_lth_CS);
				TBranch        *b_lph_CS;				results_CS->SetBranchAddress("lph", &lph_CS, &b_lph_CS);
				TBranch        *b_ltp_CS;				results_CS->SetBranchAddress("ltp", &ltp_CS, &b_ltp_CS);
				TBranch        *b_lthstar_CS;			results_CS->SetBranchAddress("lthstar", &lthstar_CS, &b_lthstar_CS);
				TBranch        *b_lphstar_CS;			results_CS->SetBranchAddress("lphstar", &lphstar_CS, &b_lphstar_CS);
				TBranch        *b_ltilde_CS;			results_CS->SetBranchAddress("ltilde", &ltilde_CS, &b_ltilde_CS);
				TBranch        *b_positivity_CS;		results_CS->SetBranchAddress("positivity", &positivity_CS, &b_positivity_CS);


				TTree* results_HX = (TTree*) results->Get("lambdaHX");
				TBranch        *b_lth_HX;				results_HX->SetBranchAddress("lth", &lth_HX, &b_lth_HX);
				TBranch        *b_lph_HX;				results_HX->SetBranchAddress("lph", &lph_HX, &b_lph_HX);
				TBranch        *b_ltp_HX;				results_HX->SetBranchAddress("ltp", &ltp_HX, &b_ltp_HX);
				TBranch        *b_lthstar_HX;			results_HX->SetBranchAddress("lthstar", &lthstar_HX, &b_lthstar_HX);
				TBranch        *b_lphstar_HX;			results_HX->SetBranchAddress("lphstar", &lphstar_HX, &b_lphstar_HX);
				TBranch        *b_ltilde_HX;			results_HX->SetBranchAddress("ltilde", &ltilde_HX, &b_ltilde_HX);
				TBranch        *b_positivity_HX;		results_HX->SetBranchAddress("positivity", &positivity_HX, &b_positivity_HX);

				TTree* results_PX = (TTree*) results->Get("lambdaPX");
				TBranch        *b_lth_PX;				results_PX->SetBranchAddress("lth", &lth_PX, &b_lth_PX);
				TBranch        *b_lph_PX;				results_PX->SetBranchAddress("lph", &lph_PX, &b_lph_PX);
				TBranch        *b_ltp_PX;				results_PX->SetBranchAddress("ltp", &ltp_PX, &b_ltp_PX);
				TBranch        *b_lthstar_PX;			results_PX->SetBranchAddress("lthstar", &lthstar_PX, &b_lthstar_PX);
				TBranch        *b_lphstar_PX;			results_PX->SetBranchAddress("lphstar", &lphstar_PX, &b_lphstar_PX);
				TBranch        *b_ltilde_PX;			results_PX->SetBranchAddress("ltilde", &ltilde_PX, &b_ltilde_PX);
				TBranch        *b_positivity_PX;		results_PX->SetBranchAddress("positivity", &positivity_PX, &b_positivity_PX);



				//	  cout<<"max "<<results_HX->GetMaximum("lph")<<endl;
				//	  cout<<"min "<<results_HX->GetMinimum("lph")<<endl;
				//	  cout<<"entries "<<results_HX->GetEntries()<<endl;

				///// Extract numerical values of the lambda parameters from the result TTrees /////

				// extremes and binning of lambda plots
				double l_bins = 100;

				if(RealData) l_bins = 200;

				TH1D* h_lth_CS = new TH1D( "h_lth_CS", "", l_bins, results_CS->GetMinimum("lth"), results_CS->GetMaximum("lth") );
				TH1D* h_lph_CS = new TH1D( "h_lph_CS", "", l_bins, results_CS->GetMinimum("lph"), results_CS->GetMaximum("lph") );
				TH1D* h_ltp_CS = new TH1D( "h_ltp_CS", "", l_bins, results_CS->GetMinimum("ltp"), results_CS->GetMaximum("ltp") );
				TH1D* h_ltilde_CS = new TH1D( "h_ltilde_CS", "", l_bins, results_CS->GetMinimum("ltilde"), results_CS->GetMaximum("ltilde") );
				TH1D* h_lphstar_CS = new TH1D( "h_lphstar_CS", "", l_bins, results_CS->GetMinimum("lphstar"), results_CS->GetMaximum("lphstar") );
				TH1D* h_lthstar_CS = new TH1D( "h_lthstar_CS", "", l_bins, results_CS->GetMinimum("lthstar"), results_CS->GetMaximum("lthstar") );

				TH1D* h_lth_HX = new TH1D( "h_lth_HX", "", l_bins, results_HX->GetMinimum("lth"), results_HX->GetMaximum("lth") );
				TH1D* h_lph_HX = new TH1D( "h_lph_HX", "", l_bins, results_HX->GetMinimum("lph"), results_HX->GetMaximum("lph") );
				TH1D* h_ltp_HX = new TH1D( "h_ltp_HX", "", l_bins, results_HX->GetMinimum("ltp"), results_HX->GetMaximum("ltp") );
				TH1D* h_ltilde_HX = new TH1D( "h_ltilde_HX", "", l_bins, results_HX->GetMinimum("ltilde"), results_HX->GetMaximum("ltilde") );
				TH1D* h_lphstar_HX = new TH1D( "h_lphstar_HX", "", l_bins, results_HX->GetMinimum("lphstar"), results_HX->GetMaximum("lphstar") );
				TH1D* h_lthstar_HX = new TH1D( "h_lthstar_HX", "", l_bins, results_HX->GetMinimum("lthstar"), results_HX->GetMaximum("lthstar") );

				TH1D* h_lth_PX = new TH1D( "h_lth_PX", "", l_bins, results_PX->GetMinimum("lth"), results_PX->GetMaximum("lth") );
				TH1D* h_lph_PX = new TH1D( "h_lph_PX", "", l_bins, results_PX->GetMinimum("lph"), results_PX->GetMaximum("lph") );
				TH1D* h_ltp_PX = new TH1D( "h_ltp_PX", "", l_bins, results_PX->GetMinimum("ltp"), results_PX->GetMaximum("ltp") );
				TH1D* h_ltilde_PX = new TH1D( "h_ltilde_PX", "", l_bins, results_PX->GetMinimum("ltilde"), results_PX->GetMaximum("ltilde") );
				TH1D* h_lphstar_PX = new TH1D( "h_lphstar_PX", "", l_bins, results_PX->GetMinimum("lphstar"), results_PX->GetMaximum("lphstar") );
				TH1D* h_lthstar_PX = new TH1D( "h_lthstar_PX", "", l_bins, results_PX->GetMinimum("lthstar"), results_PX->GetMaximum("lthstar") );

				// loop over entries in the ntuples

				results_CS->Draw("lth>>h_lth_CS");
				cout<<"Projected lth"<<endl;
				results_CS->Draw("lph>>h_lph_CS");
				cout<<"Projected lph"<<endl;
				results_CS->Draw("ltp>>h_ltp_CS");
				cout<<"Projected ltp"<<endl;
				results_CS->Draw("ltilde>>h_ltilde_CS");
				cout<<"Projected ltilde"<<endl;
				results_CS->Draw("lthstar>>h_lthstar_CS");
				cout<<"Projected lthstar"<<endl;
				results_CS->Draw("lphstar>>h_lphstar_CS");
				cout<<"Projected lphstar"<<endl;

				results_HX->Draw("lth>>h_lth_HX");
				results_HX->Draw("lph>>h_lph_HX");
				results_HX->Draw("ltp>>h_ltp_HX");
				results_HX->Draw("ltilde>>h_ltilde_HX");
				results_HX->Draw("lthstar>>h_lthstar_HX");
				results_HX->Draw("lphstar>>h_lphstar_HX");

				results_PX->Draw("lth>>h_lth_PX");
				results_PX->Draw("lph>>h_lph_PX");
				results_PX->Draw("ltp>>h_ltp_PX");
				results_PX->Draw("ltilde>>h_ltilde_PX");
				results_PX->Draw("lthstar>>h_lthstar_PX");
				results_PX->Draw("lphstar>>h_lphstar_PX");





				int n_entries = int( results_CS->GetEntries() );
				/*     for ( int i_entry = 0; i_entry < n_entries; i_entry++ ) {
							 results_CS->GetEvent( i_entry );
							 h_lth_CS->Fill( lth_CS );
							 h_lph_CS->Fill( lph_CS );
							 h_ltp_CS->Fill( ltp_CS );
							 h_ltilde_CS->Fill( ltilde_CS );
							 h_lphstar_CS->Fill( lphstar_CS );
							 h_lthstar_CS->Fill( lthstar_CS );
				//        cout<<"ltilde_CS "<<ltilde_CS<<endl;
				}

				n_entries = int( results_HX->GetEntries() );
				for ( int i_entry = 0; i_entry < n_entries; i_entry++ ) {
				results_HX->GetEvent( i_entry );
				h_lth_HX->Fill( lth_HX );
				h_lph_HX->Fill( lph_HX );
				h_ltp_HX->Fill( ltp_HX );
				h_ltilde_HX->Fill( ltilde_HX );
				h_lphstar_HX->Fill( lphstar_HX );
				h_lthstar_HX->Fill( lthstar_HX );
				}

				n_entries = int( results_PX->GetEntries() );
				for ( int i_entry = 0; i_entry < n_entries; i_entry++ ) {
				results_PX->GetEvent( i_entry );
				h_lth_PX->Fill( lth_PX );
				h_lph_PX->Fill( lph_PX );
				h_ltp_PX->Fill( ltp_PX );
				h_ltilde_PX->Fill( ltilde_PX );
				h_lphstar_PX->Fill( lphstar_PX );
				h_lthstar_PX->Fill( lthstar_PX );
				}
				 */
				// end of loop over entries in the ntuples

				//  cout << endl << endl;

				//////////////////////////////////////////////////
				// Implement here the extraction of the result and error from posterior distributions, equally for Toy and Data!
				//////////////////////////////////////////////////



				FindMPV(h_lth_CS,lth_CS,err_low_lth_CS,err_high_lth_CS,MPValgo,nSigma);
				FindMPV(h_lph_CS,lph_CS,err_low_lph_CS,err_high_lph_CS,MPValgo,nSigma);
				FindMPV(h_ltp_CS,ltp_CS,err_low_ltp_CS,err_high_ltp_CS,MPValgo,nSigma);
				FindMPV(h_lthstar_CS,lthstar_CS,err_low_lthstar_CS,err_high_lthstar_CS,MPValgo,nSigma);
				FindMPV(h_lphstar_CS,lphstar_CS,err_low_lphstar_CS,err_high_lphstar_CS,MPValgo,nSigma);
				FindMPV(h_ltilde_CS,ltilde_CS,err_low_ltilde_CS,err_high_ltilde_CS,MPValgo,nSigma);

				FindMPV(h_lth_HX,lth_HX,err_low_lth_HX,err_high_lth_HX,MPValgo,nSigma);
				FindMPV(h_lph_HX,lph_HX,err_low_lph_HX,err_high_lph_HX,MPValgo,nSigma);
				FindMPV(h_ltp_HX,ltp_HX,err_low_ltp_HX,err_high_ltp_HX,MPValgo,nSigma);
				FindMPV(h_lthstar_HX,lthstar_HX,err_low_lthstar_HX,err_high_lthstar_HX,MPValgo,nSigma);
				FindMPV(h_lphstar_HX,lphstar_HX,err_low_lphstar_HX,err_high_lphstar_HX,MPValgo,nSigma);
				FindMPV(h_ltilde_HX,ltilde_HX,err_low_ltilde_HX,err_high_ltilde_HX,MPValgo,nSigma);

				FindMPV(h_lth_PX,lth_PX,err_low_lth_PX,err_high_lth_PX,MPValgo,nSigma);
				FindMPV(h_lph_PX,lph_PX,err_low_lph_PX,err_high_lph_PX,MPValgo,nSigma);
				FindMPV(h_ltp_PX,ltp_PX,err_low_ltp_PX,err_high_ltp_PX,MPValgo,nSigma);
				FindMPV(h_lthstar_PX,lthstar_PX,err_low_lthstar_PX,err_high_lthstar_PX,MPValgo,nSigma);
				FindMPV(h_lphstar_PX,lphstar_PX,err_low_lphstar_PX,err_high_lphstar_PX,MPValgo,nSigma);
				FindMPV(h_ltilde_PX,ltilde_PX,err_low_ltilde_PX,err_high_ltilde_PX,MPValgo,nSigma);

				err_lth_CS=(err_low_lth_CS+err_high_lth_CS)/2.;
				err_lph_CS=(err_low_lph_CS+err_high_lph_CS)/2.;
				err_ltp_CS=(err_low_ltp_CS+err_high_ltp_CS)/2.;
				err_lthstar_CS=(err_low_lthstar_CS+err_high_lthstar_CS)/2.;
				err_lphstar_CS=(err_low_lphstar_CS+err_high_lphstar_CS)/2.;
				err_ltilde_CS=(err_low_ltilde_CS+err_high_ltilde_CS)/2.;

				err_lth_HX=(err_low_lth_HX+err_high_lth_HX)/2.;
				err_lph_HX=(err_low_lph_HX+err_high_lph_HX)/2.;
				err_ltp_HX=(err_low_ltp_HX+err_high_ltp_HX)/2.;
				err_lthstar_HX=(err_low_lthstar_HX+err_high_lthstar_HX)/2.;
				err_lphstar_HX=(err_low_lphstar_HX+err_high_lphstar_HX)/2.;
				err_ltilde_HX=(err_low_ltilde_HX+err_high_ltilde_HX)/2.;

				err_lth_PX=(err_low_lth_PX+err_high_lth_PX)/2.;
				err_lph_PX=(err_low_lph_PX+err_high_lph_PX)/2.;
				err_ltp_PX=(err_low_ltp_PX+err_high_ltp_PX)/2.;
				err_lthstar_PX=(err_low_lthstar_PX+err_high_lthstar_PX)/2.;
				err_lphstar_PX=(err_low_lphstar_PX+err_high_lphstar_PX)/2.;
				err_ltilde_PX=(err_low_ltilde_PX+err_high_ltilde_PX)/2.;



				if(RealData){
					char filename[200];
					sprintf(filename,"%s/Figures/PosteriorDist_lth_CS_%s_rap%d_pt%d_cpm%d.pdf",dirstruct,TreeID,iRap,iPt,icpm,MPValgo);
					PlotPosterior(icpm,iPt, iRap, "#lambda^{CS}_{#theta}", filename, h_lth_CS, lth_CS, err_low_lth_CS, err_high_lth_CS);
					sprintf(filename,"%s/Figures/PosteriorDist_lph_CS_%s_rap%d_pt%d_cpm%d.pdf",dirstruct,TreeID,iRap,iPt,icpm,MPValgo);
					PlotPosterior(icpm,iPt, iRap, "#lambda^{CS}_{#phi}", filename, h_lph_CS, lph_CS, err_low_lph_CS, err_high_lph_CS);
					sprintf(filename,"%s/Figures/PosteriorDist_ltp_CS_%s_rap%d_pt%d_cpm%d.pdf",dirstruct,TreeID,iRap,iPt,icpm,MPValgo);
					PlotPosterior(icpm,iPt, iRap, "#lambda^{CS}_{#theta#phi}", filename, h_ltp_CS, ltp_CS, err_low_ltp_CS, err_high_ltp_CS);
					sprintf(filename,"%s/Figures/PosteriorDist_ltilde_CS_%s_rap%d_pt%d_cpm%d.pdf",dirstruct,TreeID,iRap,iPt,icpm,MPValgo);
					PlotPosterior(icpm,iPt, iRap, "#tilde{#lambda}^{CS}", filename, h_ltilde_CS, ltilde_CS, err_low_ltilde_CS, err_high_ltilde_CS);
					sprintf(filename,"%s/Figures/PosteriorDist_lthstar_CS_%s_rap%d_pt%d_cpm%d.pdf",dirstruct,TreeID,iRap,iPt,icpm,MPValgo);
					PlotPosterior(icpm,iPt, iRap, "#lambda^{*CS}_{#theta}", filename, h_lthstar_CS, lthstar_CS, err_low_lthstar_CS, err_high_lthstar_CS);
					sprintf(filename,"%s/Figures/PosteriorDist_lphstar_CS_%s_rap%d_pt%d_cpm%d.pdf",dirstruct,TreeID,iRap,iPt,icpm,MPValgo);
					PlotPosterior(icpm,iPt, iRap, "#lambda^{*CS}_{#phi}", filename, h_lphstar_CS, lphstar_CS, err_low_lphstar_CS, err_high_lphstar_CS);

					sprintf(filename,"%s/Figures/PosteriorDist_lth_HX_%s_rap%d_pt%d_cpm%d.pdf",dirstruct,TreeID,iRap,iPt,icpm,MPValgo);
					PlotPosterior(icpm,iPt, iRap, "#lambda^{HX}_{#theta}", filename, h_lth_HX, lth_HX, err_low_lth_HX, err_high_lth_HX);
					sprintf(filename,"%s/Figures/PosteriorDist_lph_HX_%s_rap%d_pt%d_cpm%d.pdf",dirstruct,TreeID,iRap,iPt,icpm,MPValgo);
					PlotPosterior(icpm,iPt, iRap, "#lambda^{HX}_{#phi}", filename, h_lph_HX, lph_HX, err_low_lph_HX, err_high_lph_HX);
					sprintf(filename,"%s/Figures/PosteriorDist_ltp_HX_%s_rap%d_pt%d_cpm%d.pdf",dirstruct,TreeID,iRap,iPt,icpm,MPValgo);
					PlotPosterior(icpm,iPt, iRap, "#lambda^{HX}_{#theta#phi}", filename, h_ltp_HX, ltp_HX, err_low_ltp_HX, err_high_ltp_HX);
					sprintf(filename,"%s/Figures/PosteriorDist_ltilde_HX_%s_rap%d_pt%d_cpm%d.pdf",dirstruct,TreeID,iRap,iPt,icpm,MPValgo);
					PlotPosterior(icpm,iPt, iRap, "#tilde{#lambda}^{HX}", filename, h_ltilde_HX, ltilde_HX, err_low_ltilde_HX, err_high_ltilde_HX);
					sprintf(filename,"%s/Figures/PosteriorDist_lthstar_HX_%s_rap%d_pt%d_cpm%d.pdf",dirstruct,TreeID,iRap,iPt,icpm,MPValgo);
					PlotPosterior(icpm,iPt, iRap, "#lambda^{*HX}_{#theta}", filename, h_lthstar_HX, lthstar_HX, err_low_lthstar_HX, err_high_lthstar_HX);
					sprintf(filename,"%s/Figures/PosteriorDist_lphstar_HX_%s_rap%d_pt%d_cpm%d.pdf",dirstruct,TreeID,iRap,iPt,icpm,MPValgo);
					PlotPosterior(icpm,iPt, iRap, "#lambda^{*HX}_{#phi}", filename, h_lphstar_HX, lphstar_HX, err_low_lphstar_HX, err_high_lphstar_HX);

					sprintf(filename,"%s/Figures/PosteriorDist_lth_PX_%s_rap%d_pt%d_cpm%d.pdf",dirstruct,TreeID,iRap,iPt,icpm,MPValgo);
					PlotPosterior(icpm,iPt, iRap, "#lambda^{PX}_{#theta}", filename, h_lth_PX, lth_PX, err_low_lth_PX, err_high_lth_PX);
					sprintf(filename,"%s/Figures/PosteriorDist_lph_PX_%s_rap%d_pt%d_cpm%d.pdf",dirstruct,TreeID,iRap,iPt,icpm,MPValgo);
					PlotPosterior(icpm,iPt, iRap, "#lambda^{PX}_{#phi}", filename, h_lph_PX, lph_PX, err_low_lph_PX, err_high_lph_PX);
					sprintf(filename,"%s/Figures/PosteriorDist_ltp_PX_%s_rap%d_pt%d_cpm%d.pdf",dirstruct,TreeID,iRap,iPt,icpm,MPValgo);
					PlotPosterior(icpm,iPt, iRap, "#lambda^{PX}_{#theta#phi}", filename, h_ltp_PX, ltp_PX, err_low_ltp_PX, err_high_ltp_PX);
					sprintf(filename,"%s/Figures/PosteriorDist_ltilde_PX_%s_rap%d_pt%d_cpm%d.pdf",dirstruct,TreeID,iRap,iPt,icpm,MPValgo);
					PlotPosterior(icpm,iPt, iRap, "#tilde{#lambda}^{PX}", filename, h_ltilde_PX, ltilde_PX, err_low_ltilde_PX, err_high_ltilde_PX);
					sprintf(filename,"%s/Figures/PosteriorDist_lthstar_PX_%s_rap%d_pt%d_cpm%d.pdf",dirstruct,TreeID,iRap,iPt,icpm,MPValgo);
					PlotPosterior(icpm,iPt, iRap, "#lambda^{*PX}_{#theta}", filename, h_lthstar_PX, lthstar_PX, err_low_lthstar_PX, err_high_lthstar_PX);
					sprintf(filename,"%s/Figures/PosteriorDist_lphstar_PX_%s_rap%d_pt%d_cpm%d.pdf",dirstruct,TreeID,iRap,iPt,icpm,MPValgo);
					PlotPosterior(icpm,iPt, iRap, "#lambda^{*PX}_{#phi}", filename, h_lphstar_PX, lphstar_PX, err_low_lphstar_PX, err_high_lphstar_PX);

				}







				param_lth_CS->Fill( lth_CS );
				param_lph_CS->Fill( lph_CS );
				param_ltp_CS->Fill( ltp_CS );
				param_ltilde_CS->Fill( ltilde_CS );
				param_lphstar_CS->Fill( lphstar_CS );
				param_lthstar_CS->Fill( lthstar_CS );

				param_lth_HX->Fill( lth_HX );
				param_lph_HX->Fill( lph_HX );
				param_ltp_HX->Fill( ltp_HX );
				param_ltilde_HX->Fill( ltilde_HX );
				param_lphstar_HX->Fill( lphstar_HX );
				param_lthstar_HX->Fill( lthstar_HX );

				param_lth_PX->Fill( lth_PX );
				param_lph_PX->Fill( lph_PX );
				param_ltp_PX->Fill( ltp_PX );
				param_ltilde_PX->Fill( ltilde_PX );
				param_lphstar_PX->Fill( lphstar_PX );
				param_lthstar_PX->Fill( lthstar_PX );



				err_low_param_lth_CS->Fill( err_low_lth_CS );
				err_high_param_lth_CS->Fill( err_high_lth_CS );
				err_param_lth_CS->Fill( err_lth_CS );

				err_low_param_lph_CS->Fill( err_low_lph_CS );
				err_high_param_lph_CS->Fill( err_high_lph_CS );
				err_param_lph_CS->Fill( err_lph_CS );

				err_low_param_ltp_CS->Fill( err_low_ltp_CS );
				err_high_param_ltp_CS->Fill( err_high_ltp_CS );
				err_param_ltp_CS->Fill( err_ltp_CS );

				err_low_param_ltilde_CS->Fill( err_low_ltilde_CS );
				err_high_param_ltilde_CS->Fill( err_high_ltilde_CS );
				err_param_ltilde_CS->Fill( err_ltilde_CS );

				err_low_param_lphstar_CS->Fill( err_low_lphstar_CS );
				err_high_param_lphstar_CS->Fill( err_high_lphstar_CS );
				err_param_lphstar_CS->Fill( err_lphstar_CS );

				err_low_param_lthstar_CS->Fill( err_low_lthstar_CS );
				err_high_param_lthstar_CS->Fill( err_high_lthstar_CS );
				err_param_lthstar_CS->Fill( err_lthstar_CS );



				err_low_param_lth_HX->Fill( err_low_lth_HX );
				err_high_param_lth_HX->Fill( err_high_lth_HX );
				err_param_lth_HX->Fill( err_lth_HX );

				err_low_param_lph_HX->Fill( err_low_lph_HX );
				err_high_param_lph_HX->Fill( err_high_lph_HX );
				err_param_lph_HX->Fill( err_lph_HX );

				err_low_param_ltp_HX->Fill( err_low_ltp_HX );
				err_high_param_ltp_HX->Fill( err_high_ltp_HX );
				err_param_ltp_HX->Fill( err_ltp_HX );

				err_low_param_ltilde_HX->Fill( err_low_ltilde_HX );
				err_high_param_ltilde_HX->Fill( err_high_ltilde_HX );
				err_param_ltilde_HX->Fill( err_ltilde_HX );

				err_low_param_lphstar_HX->Fill( err_low_lphstar_HX );
				err_high_param_lphstar_HX->Fill( err_high_lphstar_HX );
				err_param_lphstar_HX->Fill( err_lphstar_HX );

				err_low_param_lthstar_HX->Fill( err_low_lthstar_HX );
				err_high_param_lthstar_HX->Fill( err_high_lthstar_HX );
				err_param_lthstar_HX->Fill( err_lthstar_HX );



				err_low_param_lth_PX->Fill( err_low_lth_PX );
				err_high_param_lth_PX->Fill( err_high_lth_PX );
				err_param_lth_PX->Fill( err_lth_PX );

				err_low_param_lph_PX->Fill( err_low_lph_PX );
				err_high_param_lph_PX->Fill( err_high_lph_PX );
				err_param_lph_PX->Fill( err_lph_PX );

				err_low_param_ltp_PX->Fill( err_low_ltp_PX );
				err_high_param_ltp_PX->Fill( err_high_ltp_PX );
				err_param_ltp_PX->Fill( err_ltp_PX );

				err_low_param_ltilde_PX->Fill( err_low_ltilde_PX );
				err_high_param_ltilde_PX->Fill( err_high_ltilde_PX );
				err_param_ltilde_PX->Fill( err_ltilde_PX );

				err_low_param_lphstar_PX->Fill( err_low_lphstar_PX );
				err_high_param_lphstar_PX->Fill( err_high_lphstar_PX );
				err_param_lphstar_PX->Fill( err_lphstar_PX );

				err_low_param_lthstar_PX->Fill( err_low_lthstar_PX );
				err_high_param_lthstar_PX->Fill( err_high_lthstar_PX );
				err_param_lthstar_PX->Fill( err_lthstar_PX );

				pull_lth_CS->Fill( (lth_CS-lambda_theta_injected_CS[rap-1][pt-1][cpm-1])/err_lth_CS );
				pull_lph_CS->Fill( (lph_CS-lambda_phi_injected_CS[rap-1][pt-1][cpm-1])/err_lph_CS );
				pull_ltp_CS->Fill( (ltp_CS-lambda_thetaphi_injected_CS[rap-1][pt-1][cpm-1])/err_ltp_CS );
				pull_lthstar_CS->Fill( (lthstar_CS-lambda_thetastar_injected_CS[rap-1][pt-1][cpm-1])/err_lthstar_CS );
				pull_lphstar_CS->Fill( (lphstar_CS-lambda_phistar_injected_CS[rap-1][pt-1][cpm-1])/err_lphstar_CS );
				pull_ltilde_CS->Fill( (ltilde_CS-lambda_tilde_injected_CS[rap-1][pt-1][cpm-1])/err_ltilde_CS );

				pull_lth_HX->Fill( (lth_HX-lambda_theta_injected_HX[rap-1][pt-1][cpm-1])/err_lth_HX );
				pull_lph_HX->Fill( (lph_HX-lambda_phi_injected_HX[rap-1][pt-1][cpm-1])/err_lph_HX );
				pull_ltp_HX->Fill( (ltp_HX-lambda_thetaphi_injected_HX[rap-1][pt-1][cpm-1])/err_ltp_HX );
				pull_lthstar_HX->Fill( (lthstar_HX-lambda_thetastar_injected_HX[rap-1][pt-1][cpm-1])/err_lthstar_HX );
				pull_lphstar_HX->Fill( (lphstar_HX-lambda_phistar_injected_HX[rap-1][pt-1][cpm-1])/err_lphstar_HX );
				pull_ltilde_HX->Fill( (ltilde_HX-lambda_tilde_injected_HX[rap-1][pt-1][cpm-1])/err_ltilde_HX );

				pull_lth_PX->Fill( (lth_PX-lambda_theta_injected_PX[rap-1][pt-1][cpm-1])/err_lth_PX );
				pull_lph_PX->Fill( (lph_PX-lambda_phi_injected_PX[rap-1][pt-1][cpm-1])/err_lph_PX );
				pull_ltp_PX->Fill( (ltp_PX-lambda_thetaphi_injected_PX[rap-1][pt-1][cpm-1])/err_ltp_PX );
				pull_lthstar_PX->Fill( (lthstar_PX-lambda_thetastar_injected_PX[rap-1][pt-1][cpm-1])/err_lthstar_PX );
				pull_lphstar_PX->Fill( (lphstar_PX-lambda_phistar_injected_PX[rap-1][pt-1][cpm-1])/err_lphstar_PX );
				pull_ltilde_PX->Fill( (ltilde_PX-lambda_tilde_injected_PX[rap-1][pt-1][cpm-1])/err_ltilde_PX );


				results->Close();

			}

			cout<<nGen_[rap-1][pt-1][cpm-1]<<endl;

			emptyBin[rap-1][pt-1][cpm-1]=false;
			if(nGen_[rap-1][pt-1][cpm-1]==0) {emptyBin[rap-1][pt-1][cpm-1]=true; cout<<"empty bin"<<endl;}

			pull_lth_CS_mean[rap-1][pt-1][cpm-1]=pull_lth_CS->GetMean(); pull_lth_CS_meanerr[rap-1][pt-1][cpm-1]=pull_lth_CS->GetMeanError(); pull_lth_CS_sigma[rap-1][pt-1][cpm-1]=pull_lth_CS->GetRMS(); pull_lth_CS_sigmaerr[rap-1][pt-1][cpm-1]=pull_lth_CS->GetRMSError();
			pull_lph_CS_mean[rap-1][pt-1][cpm-1]=pull_lph_CS->GetMean(); pull_lph_CS_meanerr[rap-1][pt-1][cpm-1]=pull_lph_CS->GetMeanError(); pull_lph_CS_sigma[rap-1][pt-1][cpm-1]=pull_lph_CS->GetRMS(); pull_lph_CS_sigmaerr[rap-1][pt-1][cpm-1]=pull_lph_CS->GetRMSError();
			pull_ltp_CS_mean[rap-1][pt-1][cpm-1]=pull_ltp_CS->GetMean(); pull_ltp_CS_meanerr[rap-1][pt-1][cpm-1]=pull_ltp_CS->GetMeanError(); pull_ltp_CS_sigma[rap-1][pt-1][cpm-1]=pull_ltp_CS->GetRMS(); pull_ltp_CS_sigmaerr[rap-1][pt-1][cpm-1]=pull_ltp_CS->GetRMSError();
			pull_lthstar_CS_mean[rap-1][pt-1][cpm-1]=pull_lthstar_CS->GetMean(); pull_lthstar_CS_meanerr[rap-1][pt-1][cpm-1]=pull_lthstar_CS->GetMeanError(); pull_lthstar_CS_sigma[rap-1][pt-1][cpm-1]=pull_lthstar_CS->GetRMS(); pull_lthstar_CS_sigmaerr[rap-1][pt-1][cpm-1]=pull_lthstar_CS->GetRMSError();
			pull_lphstar_CS_mean[rap-1][pt-1][cpm-1]=pull_lphstar_CS->GetMean(); pull_lphstar_CS_meanerr[rap-1][pt-1][cpm-1]=pull_lphstar_CS->GetMeanError(); pull_lphstar_CS_sigma[rap-1][pt-1][cpm-1]=pull_lphstar_CS->GetRMS(); pull_lphstar_CS_sigmaerr[rap-1][pt-1][cpm-1]=pull_lphstar_CS->GetRMSError();
			pull_ltilde_CS_mean[rap-1][pt-1][cpm-1]=pull_ltilde_CS->GetMean(); pull_ltilde_CS_meanerr[rap-1][pt-1][cpm-1]=pull_ltilde_CS->GetMeanError(); pull_ltilde_CS_sigma[rap-1][pt-1][cpm-1]=pull_ltilde_CS->GetRMS(); pull_ltilde_CS_sigmaerr[rap-1][pt-1][cpm-1]=pull_ltilde_CS->GetRMSError();

			pull_lth_HX_mean[rap-1][pt-1][cpm-1]=pull_lth_HX->GetMean(); pull_lth_HX_meanerr[rap-1][pt-1][cpm-1]=pull_lth_HX->GetMeanError(); pull_lth_HX_sigma[rap-1][pt-1][cpm-1]=pull_lth_HX->GetRMS(); pull_lth_HX_sigmaerr[rap-1][pt-1][cpm-1]=pull_lth_HX->GetRMSError();
			pull_lph_HX_mean[rap-1][pt-1][cpm-1]=pull_lph_HX->GetMean(); pull_lph_HX_meanerr[rap-1][pt-1][cpm-1]=pull_lph_HX->GetMeanError(); pull_lph_HX_sigma[rap-1][pt-1][cpm-1]=pull_lph_HX->GetRMS(); pull_lph_HX_sigmaerr[rap-1][pt-1][cpm-1]=pull_lph_HX->GetRMSError();
			pull_ltp_HX_mean[rap-1][pt-1][cpm-1]=pull_ltp_HX->GetMean(); pull_ltp_HX_meanerr[rap-1][pt-1][cpm-1]=pull_ltp_HX->GetMeanError(); pull_ltp_HX_sigma[rap-1][pt-1][cpm-1]=pull_ltp_HX->GetRMS(); pull_ltp_HX_sigmaerr[rap-1][pt-1][cpm-1]=pull_ltp_HX->GetRMSError();
			pull_lthstar_HX_mean[rap-1][pt-1][cpm-1]=pull_lthstar_HX->GetMean(); pull_lthstar_HX_meanerr[rap-1][pt-1][cpm-1]=pull_lthstar_HX->GetMeanError(); pull_lthstar_HX_sigma[rap-1][pt-1][cpm-1]=pull_lthstar_HX->GetRMS(); pull_lthstar_HX_sigmaerr[rap-1][pt-1][cpm-1]=pull_lthstar_HX->GetRMSError();
			pull_lphstar_HX_mean[rap-1][pt-1][cpm-1]=pull_lphstar_HX->GetMean(); pull_lphstar_HX_meanerr[rap-1][pt-1][cpm-1]=pull_lphstar_HX->GetMeanError(); pull_lphstar_HX_sigma[rap-1][pt-1][cpm-1]=pull_lphstar_HX->GetRMS(); pull_lphstar_HX_sigmaerr[rap-1][pt-1][cpm-1]=pull_lphstar_HX->GetRMSError();
			pull_ltilde_HX_mean[rap-1][pt-1][cpm-1]=pull_ltilde_HX->GetMean(); pull_ltilde_HX_meanerr[rap-1][pt-1][cpm-1]=pull_ltilde_HX->GetMeanError(); pull_ltilde_HX_sigma[rap-1][pt-1][cpm-1]=pull_ltilde_HX->GetRMS(); pull_ltilde_HX_sigmaerr[rap-1][pt-1][cpm-1]=pull_ltilde_HX->GetRMSError();

			pull_lth_PX_mean[rap-1][pt-1][cpm-1]=pull_lth_PX->GetMean(); pull_lth_PX_meanerr[rap-1][pt-1][cpm-1]=pull_lth_PX->GetMeanError(); pull_lth_PX_sigma[rap-1][pt-1][cpm-1]=pull_lth_PX->GetRMS(); pull_lth_PX_sigmaerr[rap-1][pt-1][cpm-1]=pull_lth_PX->GetRMSError();
			pull_lph_PX_mean[rap-1][pt-1][cpm-1]=pull_lph_PX->GetMean(); pull_lph_PX_meanerr[rap-1][pt-1][cpm-1]=pull_lph_PX->GetMeanError(); pull_lph_PX_sigma[rap-1][pt-1][cpm-1]=pull_lph_PX->GetRMS(); pull_lph_PX_sigmaerr[rap-1][pt-1][cpm-1]=pull_lph_PX->GetRMSError();
			pull_ltp_PX_mean[rap-1][pt-1][cpm-1]=pull_ltp_PX->GetMean(); pull_ltp_PX_meanerr[rap-1][pt-1][cpm-1]=pull_ltp_PX->GetMeanError(); pull_ltp_PX_sigma[rap-1][pt-1][cpm-1]=pull_ltp_PX->GetRMS(); pull_ltp_PX_sigmaerr[rap-1][pt-1][cpm-1]=pull_ltp_PX->GetRMSError();
			pull_lthstar_PX_mean[rap-1][pt-1][cpm-1]=pull_lthstar_PX->GetMean(); pull_lthstar_PX_meanerr[rap-1][pt-1][cpm-1]=pull_lthstar_PX->GetMeanError(); pull_lthstar_PX_sigma[rap-1][pt-1][cpm-1]=pull_lthstar_PX->GetRMS(); pull_lthstar_PX_sigmaerr[rap-1][pt-1][cpm-1]=pull_lthstar_PX->GetRMSError();
			pull_lphstar_PX_mean[rap-1][pt-1][cpm-1]=pull_lphstar_PX->GetMean(); pull_lphstar_PX_meanerr[rap-1][pt-1][cpm-1]=pull_lphstar_PX->GetMeanError(); pull_lphstar_PX_sigma[rap-1][pt-1][cpm-1]=pull_lphstar_PX->GetRMS(); pull_lphstar_PX_sigmaerr[rap-1][pt-1][cpm-1]=pull_lphstar_PX->GetRMSError();
			pull_ltilde_PX_mean[rap-1][pt-1][cpm-1]=pull_ltilde_PX->GetMean(); pull_ltilde_PX_meanerr[rap-1][pt-1][cpm-1]=pull_ltilde_PX->GetMeanError(); pull_ltilde_PX_sigma[rap-1][pt-1][cpm-1]=pull_ltilde_PX->GetRMS(); pull_ltilde_PX_sigmaerr[rap-1][pt-1][cpm-1]=pull_ltilde_PX->GetRMSError();

			param_lth_CS_mean[rap-1][pt-1][cpm-1]=param_lth_CS->GetMean(); param_lth_CS_meanerr[rap-1][pt-1][cpm-1]=param_lth_CS->GetMeanError();						 param_lth_CS_sigma[rap-1][pt-1][cpm-1]=err_param_lth_CS->GetMean(); 				param_lth_CS_sigma_high[rap-1][pt-1][cpm-1]=err_high_param_lth_CS->GetMean(); 				param_lth_CS_sigma_low[rap-1][pt-1][cpm-1]=err_low_param_lth_CS->GetMean(); 					param_lth_CS_sigmaerr[rap-1][pt-1][cpm-1]=err_param_lth_CS->GetMeanError(); 											if(RealData) {param_lth_CS_mean[rap-1][pt-1][cpm-1]=lth_CS; 			param_lth_CS_sigma[rap-1][pt-1][cpm-1]=err_lth_CS;        		param_lth_CS_sigma_high[rap-1][pt-1][cpm-1]=err_high_lth_CS;        		param_lth_CS_sigma_low[rap-1][pt-1][cpm-1]=err_low_lth_CS;        		}
			param_lph_CS_mean[rap-1][pt-1][cpm-1]=param_lph_CS->GetMean(); param_lph_CS_meanerr[rap-1][pt-1][cpm-1]=param_lph_CS->GetMeanError();						 param_lph_CS_sigma[rap-1][pt-1][cpm-1]=err_param_lph_CS->GetMean(); 				param_lph_CS_sigma_high[rap-1][pt-1][cpm-1]=err_high_param_lph_CS->GetMean(); 				param_lph_CS_sigma_low[rap-1][pt-1][cpm-1]=err_low_param_lph_CS->GetMean(); 					param_lph_CS_sigmaerr[rap-1][pt-1][cpm-1]=err_param_lph_CS->GetMeanError(); 											if(RealData) {param_lph_CS_mean[rap-1][pt-1][cpm-1]=lph_CS; 			param_lph_CS_sigma[rap-1][pt-1][cpm-1]=err_lph_CS;        		param_lph_CS_sigma_high[rap-1][pt-1][cpm-1]=err_high_lph_CS;        		param_lph_CS_sigma_low[rap-1][pt-1][cpm-1]=err_low_lph_CS;        		}
			param_ltp_CS_mean[rap-1][pt-1][cpm-1]=param_ltp_CS->GetMean(); param_ltp_CS_meanerr[rap-1][pt-1][cpm-1]=param_ltp_CS->GetMeanError();						 param_ltp_CS_sigma[rap-1][pt-1][cpm-1]=err_param_ltp_CS->GetMean(); 				param_ltp_CS_sigma_high[rap-1][pt-1][cpm-1]=err_high_param_ltp_CS->GetMean(); 				param_ltp_CS_sigma_low[rap-1][pt-1][cpm-1]=err_low_param_ltp_CS->GetMean(); 					param_ltp_CS_sigmaerr[rap-1][pt-1][cpm-1]=err_param_ltp_CS->GetMeanError(); 											if(RealData) {param_ltp_CS_mean[rap-1][pt-1][cpm-1]=ltp_CS; 			param_ltp_CS_sigma[rap-1][pt-1][cpm-1]=err_ltp_CS;        		param_ltp_CS_sigma_high[rap-1][pt-1][cpm-1]=err_high_ltp_CS;        		param_ltp_CS_sigma_low[rap-1][pt-1][cpm-1]=err_low_ltp_CS;        		}
			param_ltilde_CS_mean[rap-1][pt-1][cpm-1]=param_ltilde_CS->GetMean(); param_ltilde_CS_meanerr[rap-1][pt-1][cpm-1]=param_ltilde_CS->GetMeanError();						 param_ltilde_CS_sigma[rap-1][pt-1][cpm-1]=err_param_ltilde_CS->GetMean(); 				param_ltilde_CS_sigma_high[rap-1][pt-1][cpm-1]=err_high_param_ltilde_CS->GetMean(); 				param_ltilde_CS_sigma_low[rap-1][pt-1][cpm-1]=err_low_param_ltilde_CS->GetMean(); 					param_ltilde_CS_sigmaerr[rap-1][pt-1][cpm-1]=err_param_ltilde_CS->GetMeanError(); 											if(RealData) {param_ltilde_CS_mean[rap-1][pt-1][cpm-1]=ltilde_CS; 			param_ltilde_CS_sigma[rap-1][pt-1][cpm-1]=err_ltilde_CS;        		param_ltilde_CS_sigma_high[rap-1][pt-1][cpm-1]=err_high_ltilde_CS;        		param_ltilde_CS_sigma_low[rap-1][pt-1][cpm-1]=err_low_ltilde_CS;        		}
			param_lphstar_CS_mean[rap-1][pt-1][cpm-1]=param_lphstar_CS->GetMean(); param_lphstar_CS_meanerr[rap-1][pt-1][cpm-1]=param_lphstar_CS->GetMeanError();						 param_lphstar_CS_sigma[rap-1][pt-1][cpm-1]=err_param_lphstar_CS->GetMean(); 				param_lphstar_CS_sigma_high[rap-1][pt-1][cpm-1]=err_high_param_lphstar_CS->GetMean(); 				param_lphstar_CS_sigma_low[rap-1][pt-1][cpm-1]=err_low_param_lphstar_CS->GetMean(); 					param_lphstar_CS_sigmaerr[rap-1][pt-1][cpm-1]=err_param_lphstar_CS->GetMeanError(); 											if(RealData) {param_lphstar_CS_mean[rap-1][pt-1][cpm-1]=lphstar_CS; 			param_lphstar_CS_sigma[rap-1][pt-1][cpm-1]=err_lphstar_CS;        		param_lphstar_CS_sigma_high[rap-1][pt-1][cpm-1]=err_high_lphstar_CS;        		param_lphstar_CS_sigma_low[rap-1][pt-1][cpm-1]=err_low_lphstar_CS;        		}
			param_lthstar_CS_mean[rap-1][pt-1][cpm-1]=param_lthstar_CS->GetMean(); param_lthstar_CS_meanerr[rap-1][pt-1][cpm-1]=param_lthstar_CS->GetMeanError();						 param_lthstar_CS_sigma[rap-1][pt-1][cpm-1]=err_param_lthstar_CS->GetMean(); 				param_lthstar_CS_sigma_high[rap-1][pt-1][cpm-1]=err_high_param_lthstar_CS->GetMean(); 				param_lthstar_CS_sigma_low[rap-1][pt-1][cpm-1]=err_low_param_lthstar_CS->GetMean(); 					param_lthstar_CS_sigmaerr[rap-1][pt-1][cpm-1]=err_param_lthstar_CS->GetMeanError(); 											if(RealData) {param_lthstar_CS_mean[rap-1][pt-1][cpm-1]=lthstar_CS; 			param_lthstar_CS_sigma[rap-1][pt-1][cpm-1]=err_lthstar_CS;        		param_lthstar_CS_sigma_high[rap-1][pt-1][cpm-1]=err_high_lthstar_CS;        		param_lthstar_CS_sigma_low[rap-1][pt-1][cpm-1]=err_low_lthstar_CS;        		}

			param_lth_HX_mean[rap-1][pt-1][cpm-1]=param_lth_HX->GetMean(); param_lth_HX_meanerr[rap-1][pt-1][cpm-1]=param_lth_HX->GetMeanError();						 param_lth_HX_sigma[rap-1][pt-1][cpm-1]=err_param_lth_HX->GetMean(); 				param_lth_HX_sigma_high[rap-1][pt-1][cpm-1]=err_high_param_lth_HX->GetMean(); 				param_lth_HX_sigma_low[rap-1][pt-1][cpm-1]=err_low_param_lth_HX->GetMean(); 					param_lth_HX_sigmaerr[rap-1][pt-1][cpm-1]=err_param_lth_HX->GetMeanError(); 											if(RealData) {param_lth_HX_mean[rap-1][pt-1][cpm-1]=lth_HX; 			param_lth_HX_sigma[rap-1][pt-1][cpm-1]=err_lth_HX;        		param_lth_HX_sigma_high[rap-1][pt-1][cpm-1]=err_high_lth_HX;        		param_lth_HX_sigma_low[rap-1][pt-1][cpm-1]=err_low_lth_HX;        		}
			param_lph_HX_mean[rap-1][pt-1][cpm-1]=param_lph_HX->GetMean(); param_lph_HX_meanerr[rap-1][pt-1][cpm-1]=param_lph_HX->GetMeanError();						 param_lph_HX_sigma[rap-1][pt-1][cpm-1]=err_param_lph_HX->GetMean(); 				param_lph_HX_sigma_high[rap-1][pt-1][cpm-1]=err_high_param_lph_HX->GetMean(); 				param_lph_HX_sigma_low[rap-1][pt-1][cpm-1]=err_low_param_lph_HX->GetMean(); 					param_lph_HX_sigmaerr[rap-1][pt-1][cpm-1]=err_param_lph_HX->GetMeanError(); 											if(RealData) {param_lph_HX_mean[rap-1][pt-1][cpm-1]=lph_HX; 			param_lph_HX_sigma[rap-1][pt-1][cpm-1]=err_lph_HX;        		param_lph_HX_sigma_high[rap-1][pt-1][cpm-1]=err_high_lph_HX;        		param_lph_HX_sigma_low[rap-1][pt-1][cpm-1]=err_low_lph_HX;        		}
			param_ltp_HX_mean[rap-1][pt-1][cpm-1]=param_ltp_HX->GetMean(); param_ltp_HX_meanerr[rap-1][pt-1][cpm-1]=param_ltp_HX->GetMeanError();						 param_ltp_HX_sigma[rap-1][pt-1][cpm-1]=err_param_ltp_HX->GetMean(); 				param_ltp_HX_sigma_high[rap-1][pt-1][cpm-1]=err_high_param_ltp_HX->GetMean(); 				param_ltp_HX_sigma_low[rap-1][pt-1][cpm-1]=err_low_param_ltp_HX->GetMean(); 					param_ltp_HX_sigmaerr[rap-1][pt-1][cpm-1]=err_param_ltp_HX->GetMeanError(); 											if(RealData) {param_ltp_HX_mean[rap-1][pt-1][cpm-1]=ltp_HX; 			param_ltp_HX_sigma[rap-1][pt-1][cpm-1]=err_ltp_HX;        		param_ltp_HX_sigma_high[rap-1][pt-1][cpm-1]=err_high_ltp_HX;        		param_ltp_HX_sigma_low[rap-1][pt-1][cpm-1]=err_low_ltp_HX;        		}
			param_ltilde_HX_mean[rap-1][pt-1][cpm-1]=param_ltilde_HX->GetMean(); param_ltilde_HX_meanerr[rap-1][pt-1][cpm-1]=param_ltilde_HX->GetMeanError();						 param_ltilde_HX_sigma[rap-1][pt-1][cpm-1]=err_param_ltilde_HX->GetMean(); 				param_ltilde_HX_sigma_high[rap-1][pt-1][cpm-1]=err_high_param_ltilde_HX->GetMean(); 				param_ltilde_HX_sigma_low[rap-1][pt-1][cpm-1]=err_low_param_ltilde_HX->GetMean(); 					param_ltilde_HX_sigmaerr[rap-1][pt-1][cpm-1]=err_param_ltilde_HX->GetMeanError(); 											if(RealData) {param_ltilde_HX_mean[rap-1][pt-1][cpm-1]=ltilde_HX; 			param_ltilde_HX_sigma[rap-1][pt-1][cpm-1]=err_ltilde_HX;        		param_ltilde_HX_sigma_high[rap-1][pt-1][cpm-1]=err_high_ltilde_HX;        		param_ltilde_HX_sigma_low[rap-1][pt-1][cpm-1]=err_low_ltilde_HX;        		}
			param_lphstar_HX_mean[rap-1][pt-1][cpm-1]=param_lphstar_HX->GetMean(); param_lphstar_HX_meanerr[rap-1][pt-1][cpm-1]=param_lphstar_HX->GetMeanError();						 param_lphstar_HX_sigma[rap-1][pt-1][cpm-1]=err_param_lphstar_HX->GetMean(); 				param_lphstar_HX_sigma_high[rap-1][pt-1][cpm-1]=err_high_param_lphstar_HX->GetMean(); 				param_lphstar_HX_sigma_low[rap-1][pt-1][cpm-1]=err_low_param_lphstar_HX->GetMean(); 					param_lphstar_HX_sigmaerr[rap-1][pt-1][cpm-1]=err_param_lphstar_HX->GetMeanError(); 											if(RealData) {param_lphstar_HX_mean[rap-1][pt-1][cpm-1]=lphstar_HX; 			param_lphstar_HX_sigma[rap-1][pt-1][cpm-1]=err_lphstar_HX;        		param_lphstar_HX_sigma_high[rap-1][pt-1][cpm-1]=err_high_lphstar_HX;        		param_lphstar_HX_sigma_low[rap-1][pt-1][cpm-1]=err_low_lphstar_HX;        		}
			param_lthstar_HX_mean[rap-1][pt-1][cpm-1]=param_lthstar_HX->GetMean(); param_lthstar_HX_meanerr[rap-1][pt-1][cpm-1]=param_lthstar_HX->GetMeanError();						 param_lthstar_HX_sigma[rap-1][pt-1][cpm-1]=err_param_lthstar_HX->GetMean(); 				param_lthstar_HX_sigma_high[rap-1][pt-1][cpm-1]=err_high_param_lthstar_HX->GetMean(); 				param_lthstar_HX_sigma_low[rap-1][pt-1][cpm-1]=err_low_param_lthstar_HX->GetMean(); 					param_lthstar_HX_sigmaerr[rap-1][pt-1][cpm-1]=err_param_lthstar_HX->GetMeanError(); 											if(RealData) {param_lthstar_HX_mean[rap-1][pt-1][cpm-1]=lthstar_HX; 			param_lthstar_HX_sigma[rap-1][pt-1][cpm-1]=err_lthstar_HX;        		param_lthstar_HX_sigma_high[rap-1][pt-1][cpm-1]=err_high_lthstar_HX;        		param_lthstar_HX_sigma_low[rap-1][pt-1][cpm-1]=err_low_lthstar_HX;        		}

			param_lth_PX_mean[rap-1][pt-1][cpm-1]=param_lth_PX->GetMean(); param_lth_PX_meanerr[rap-1][pt-1][cpm-1]=param_lth_PX->GetMeanError();						 param_lth_PX_sigma[rap-1][pt-1][cpm-1]=err_param_lth_PX->GetMean(); 				param_lth_PX_sigma_high[rap-1][pt-1][cpm-1]=err_high_param_lth_PX->GetMean(); 				param_lth_PX_sigma_low[rap-1][pt-1][cpm-1]=err_low_param_lth_PX->GetMean(); 					param_lth_PX_sigmaerr[rap-1][pt-1][cpm-1]=err_param_lth_PX->GetMeanError(); 											if(RealData) {param_lth_PX_mean[rap-1][pt-1][cpm-1]=lth_PX; 			param_lth_PX_sigma[rap-1][pt-1][cpm-1]=err_lth_PX;        		param_lth_PX_sigma_high[rap-1][pt-1][cpm-1]=err_high_lth_PX;        		param_lth_PX_sigma_low[rap-1][pt-1][cpm-1]=err_low_lth_PX;        		}
			param_lph_PX_mean[rap-1][pt-1][cpm-1]=param_lph_PX->GetMean(); param_lph_PX_meanerr[rap-1][pt-1][cpm-1]=param_lph_PX->GetMeanError();						 param_lph_PX_sigma[rap-1][pt-1][cpm-1]=err_param_lph_PX->GetMean(); 				param_lph_PX_sigma_high[rap-1][pt-1][cpm-1]=err_high_param_lph_PX->GetMean(); 				param_lph_PX_sigma_low[rap-1][pt-1][cpm-1]=err_low_param_lph_PX->GetMean(); 					param_lph_PX_sigmaerr[rap-1][pt-1][cpm-1]=err_param_lph_PX->GetMeanError(); 											if(RealData) {param_lph_PX_mean[rap-1][pt-1][cpm-1]=lph_PX; 			param_lph_PX_sigma[rap-1][pt-1][cpm-1]=err_lph_PX;        		param_lph_PX_sigma_high[rap-1][pt-1][cpm-1]=err_high_lph_PX;        		param_lph_PX_sigma_low[rap-1][pt-1][cpm-1]=err_low_lph_PX;        		}
			param_ltp_PX_mean[rap-1][pt-1][cpm-1]=param_ltp_PX->GetMean(); param_ltp_PX_meanerr[rap-1][pt-1][cpm-1]=param_ltp_PX->GetMeanError();						 param_ltp_PX_sigma[rap-1][pt-1][cpm-1]=err_param_ltp_PX->GetMean(); 				param_ltp_PX_sigma_high[rap-1][pt-1][cpm-1]=err_high_param_ltp_PX->GetMean(); 				param_ltp_PX_sigma_low[rap-1][pt-1][cpm-1]=err_low_param_ltp_PX->GetMean(); 					param_ltp_PX_sigmaerr[rap-1][pt-1][cpm-1]=err_param_ltp_PX->GetMeanError(); 											if(RealData) {param_ltp_PX_mean[rap-1][pt-1][cpm-1]=ltp_PX; 			param_ltp_PX_sigma[rap-1][pt-1][cpm-1]=err_ltp_PX;        		param_ltp_PX_sigma_high[rap-1][pt-1][cpm-1]=err_high_ltp_PX;        		param_ltp_PX_sigma_low[rap-1][pt-1][cpm-1]=err_low_ltp_PX;        		}
			param_ltilde_PX_mean[rap-1][pt-1][cpm-1]=param_ltilde_PX->GetMean(); param_ltilde_PX_meanerr[rap-1][pt-1][cpm-1]=param_ltilde_PX->GetMeanError();						 param_ltilde_PX_sigma[rap-1][pt-1][cpm-1]=err_param_ltilde_PX->GetMean(); 				param_ltilde_PX_sigma_high[rap-1][pt-1][cpm-1]=err_high_param_ltilde_PX->GetMean(); 				param_ltilde_PX_sigma_low[rap-1][pt-1][cpm-1]=err_low_param_ltilde_PX->GetMean(); 					param_ltilde_PX_sigmaerr[rap-1][pt-1][cpm-1]=err_param_ltilde_PX->GetMeanError(); 											if(RealData) {param_ltilde_PX_mean[rap-1][pt-1][cpm-1]=ltilde_PX; 			param_ltilde_PX_sigma[rap-1][pt-1][cpm-1]=err_ltilde_PX;        		param_ltilde_PX_sigma_high[rap-1][pt-1][cpm-1]=err_high_ltilde_PX;        		param_ltilde_PX_sigma_low[rap-1][pt-1][cpm-1]=err_low_ltilde_PX;        		}
			param_lphstar_PX_mean[rap-1][pt-1][cpm-1]=param_lphstar_PX->GetMean(); param_lphstar_PX_meanerr[rap-1][pt-1][cpm-1]=param_lphstar_PX->GetMeanError();						 param_lphstar_PX_sigma[rap-1][pt-1][cpm-1]=err_param_lphstar_PX->GetMean(); 				param_lphstar_PX_sigma_high[rap-1][pt-1][cpm-1]=err_high_param_lphstar_PX->GetMean(); 				param_lphstar_PX_sigma_low[rap-1][pt-1][cpm-1]=err_low_param_lphstar_PX->GetMean(); 					param_lphstar_PX_sigmaerr[rap-1][pt-1][cpm-1]=err_param_lphstar_PX->GetMeanError(); 											if(RealData) {param_lphstar_PX_mean[rap-1][pt-1][cpm-1]=lphstar_PX; 			param_lphstar_PX_sigma[rap-1][pt-1][cpm-1]=err_lphstar_PX;        		param_lphstar_PX_sigma_high[rap-1][pt-1][cpm-1]=err_high_lphstar_PX;        		param_lphstar_PX_sigma_low[rap-1][pt-1][cpm-1]=err_low_lphstar_PX;        		}
			param_lthstar_PX_mean[rap-1][pt-1][cpm-1]=param_lthstar_PX->GetMean(); param_lthstar_PX_meanerr[rap-1][pt-1][cpm-1]=param_lthstar_PX->GetMeanError();						 param_lthstar_PX_sigma[rap-1][pt-1][cpm-1]=err_param_lthstar_PX->GetMean(); 				param_lthstar_PX_sigma_high[rap-1][pt-1][cpm-1]=err_high_param_lthstar_PX->GetMean(); 				param_lthstar_PX_sigma_low[rap-1][pt-1][cpm-1]=err_low_param_lthstar_PX->GetMean(); 					param_lthstar_PX_sigmaerr[rap-1][pt-1][cpm-1]=err_param_lthstar_PX->GetMeanError(); 											if(RealData) {param_lthstar_PX_mean[rap-1][pt-1][cpm-1]=lthstar_PX; 			param_lthstar_PX_sigma[rap-1][pt-1][cpm-1]=err_lthstar_PX;        		param_lthstar_PX_sigma_high[rap-1][pt-1][cpm-1]=err_high_lthstar_PX;        		param_lthstar_PX_sigma_low[rap-1][pt-1][cpm-1]=err_low_lthstar_PX;        		}


			/*	    sigma1Counter_lth_CS[rap-1][pt-1]=0;
							for(int iPull=0;iPull<int((l_max_pull-l_min_pull)/l_step_1D_pull)+1;iPull++){
							if(TMath::Abs(param_lth_CS->GetBinCenter())<1) sigma1Counter_lth_CS[rap-1][pt-1]+=param_lth_CS->GetBinContent();
							}
							sigma1Counter_lth_CS[rap-1][pt-1]=sigma1Counter_lth_CS[rap-1][pt-1]/param_lth_CS->GetSumOfWeights();
							pull_lth_CS_mean[rap-1][pt-1]=pull_lth_CS->GetMean(); pull_lth_CS_meanerr[rap-1][pt-1]=pull_lth_CS->GetMeanError();
							pull_lth_CS_sigma[rap-1][pt-1]=pull_lth_CS->GetRMS();
							pull_lth_CS_sigmaerr[rap-1][pt-1]=pull_lth_CS->GetRMSError();


			 */

			bool UseMedianForToys=true;
			Double_t quantile,prob;
			prob=.5;

			if(!RealData&&UseMedianForToys){

				param_lth_CS->GetQuantiles(1,&quantile,&prob); param_lth_CS_mean[rap-1][pt-1][cpm-1]=quantile;
				param_lph_CS->GetQuantiles(1,&quantile,&prob); param_lph_CS_mean[rap-1][pt-1][cpm-1]=quantile;
				param_ltp_CS->GetQuantiles(1,&quantile,&prob); param_ltp_CS_mean[rap-1][pt-1][cpm-1]=quantile;
				param_ltilde_CS->GetQuantiles(1,&quantile,&prob); param_ltilde_CS_mean[rap-1][pt-1][cpm-1]=quantile;

				param_lth_HX->GetQuantiles(1,&quantile,&prob); param_lth_HX_mean[rap-1][pt-1][cpm-1]=quantile;
				param_lph_HX->GetQuantiles(1,&quantile,&prob); param_lph_HX_mean[rap-1][pt-1][cpm-1]=quantile;
				param_ltp_HX->GetQuantiles(1,&quantile,&prob); param_ltp_HX_mean[rap-1][pt-1][cpm-1]=quantile;
				param_ltilde_HX->GetQuantiles(1,&quantile,&prob); param_ltilde_HX_mean[rap-1][pt-1][cpm-1]=quantile;

				param_lth_PX->GetQuantiles(1,&quantile,&prob); param_lth_PX_mean[rap-1][pt-1][cpm-1]=quantile;
				param_lph_PX->GetQuantiles(1,&quantile,&prob); param_lph_PX_mean[rap-1][pt-1][cpm-1]=quantile;
				param_ltp_PX->GetQuantiles(1,&quantile,&prob); param_ltp_PX_mean[rap-1][pt-1][cpm-1]=quantile;
				param_ltilde_PX->GetQuantiles(1,&quantile,&prob); param_ltilde_PX_mean[rap-1][pt-1][cpm-1]=quantile;

			}



			sprintf(filename,"%s/Figures",dirstruct); gSystem->mkdir(filename);

			if(!RealData){

				sprintf(filename,"%s/Figures/Pull_lth_CS_rap%dpt%dcpm%d.pdf",dirstruct,iRap,iPt,icpm);
				PlotObject(pull_lth_CS, "z[#lambda^{CS}_{#theta}]", l_min_pull, l_max_pull, l_step_1D_pull, pull_lth_CS_mean[rap-1][pt-1][cpm-1], pull_lth_CS_meanerr[rap-1][pt-1][cpm-1], pull_lth_CS_sigma[rap-1][pt-1][cpm-1], pull_lth_CS_sigmaerr[rap-1][pt-1][cpm-1], filename, nGenerations,true,0);

				sprintf(filename,"%s/Figures/Pull_lph_CS_rap%dpt%dcpm%d.pdf",dirstruct,iRap,iPt,icpm);
				PlotObject(pull_lph_CS, "z[#lambda^{CS}_{#phi}]", l_min_pull, l_max_pull, l_step_1D_pull, pull_lph_CS_mean[rap-1][pt-1][cpm-1], pull_lph_CS_meanerr[rap-1][pt-1][cpm-1], pull_lph_CS_sigma[rap-1][pt-1][cpm-1], pull_lph_CS_sigmaerr[rap-1][pt-1][cpm-1], filename, nGenerations,true,0);

				sprintf(filename,"%s/Figures/Pull_ltp_CS_rap%dpt%dcpm%d.pdf",dirstruct,iRap,iPt,icpm);
				PlotObject(pull_ltp_CS, "z[#lambda^{CS}_{#theta#phi}]", l_min_pull, l_max_pull, l_step_1D_pull, pull_ltp_CS_mean[rap-1][pt-1][cpm-1], pull_ltp_CS_meanerr[rap-1][pt-1][cpm-1], pull_ltp_CS_sigma[rap-1][pt-1][cpm-1], pull_ltp_CS_sigmaerr[rap-1][pt-1][cpm-1], filename, nGenerations,true,0);

				sprintf(filename,"%s/Figures/Pull_lthstar_CS_rap%dpt%dcpm%d.pdf",dirstruct,iRap,iPt,icpm);
				PlotObject(pull_lthstar_CS, "z[#lambda^{*CS}_{#theta}]", l_min_pull, l_max_pull, l_step_1D_pull, pull_lthstar_CS_mean[rap-1][pt-1][cpm-1], pull_lthstar_CS_meanerr[rap-1][pt-1][cpm-1], pull_lthstar_CS_sigma[rap-1][pt-1][cpm-1], pull_lthstar_CS_sigmaerr[rap-1][pt-1][cpm-1], filename, nGenerations,true,0);

				sprintf(filename,"%s/Figures/Pull_lphstar_CS_rap%dpt%dcpm%d.pdf",dirstruct,iRap,iPt,icpm);
				PlotObject(pull_lphstar_CS, "z[#lambda^{*CS}_{#phi}]", l_min_pull, l_max_pull, l_step_1D_pull, pull_lphstar_CS_mean[rap-1][pt-1][cpm-1], pull_lphstar_CS_meanerr[rap-1][pt-1][cpm-1], pull_lphstar_CS_sigma[rap-1][pt-1][cpm-1], pull_lphstar_CS_sigmaerr[rap-1][pt-1][cpm-1], filename, nGenerations,true,0);

				sprintf(filename,"%s/Figures/Pull_ltilde_CS_rap%dpt%dcpm%d.pdf",dirstruct,iRap,iPt,icpm);
				PlotObject(pull_ltilde_CS, "z[#tilde{#lambda}^{CS}]", l_min_pull, l_max_pull, l_step_1D_pull, pull_ltilde_CS_mean[rap-1][pt-1][cpm-1], pull_ltilde_CS_meanerr[rap-1][pt-1][cpm-1], pull_ltilde_CS_sigma[rap-1][pt-1][cpm-1], pull_ltilde_CS_sigmaerr[rap-1][pt-1][cpm-1], filename, nGenerations,true,0);


				sprintf(filename,"%s/Figures/Pull_lth_HX_rap%dpt%dcpm%d.pdf",dirstruct,iRap,iPt,icpm);
				PlotObject(pull_lth_HX, "z[#lambda^{HX}_{#theta}]", l_min_pull, l_max_pull, l_step_1D_pull, pull_lth_HX_mean[rap-1][pt-1][cpm-1], pull_lth_HX_meanerr[rap-1][pt-1][cpm-1], pull_lth_HX_sigma[rap-1][pt-1][cpm-1], pull_lth_HX_sigmaerr[rap-1][pt-1][cpm-1], filename, nGenerations,true,0);

				sprintf(filename,"%s/Figures/Pull_lph_HX_rap%dpt%dcpm%d.pdf",dirstruct,iRap,iPt,icpm);
				PlotObject(pull_lph_HX, "z[#lambda^{HX}_{#phi}]", l_min_pull, l_max_pull, l_step_1D_pull, pull_lph_HX_mean[rap-1][pt-1][cpm-1], pull_lph_HX_meanerr[rap-1][pt-1][cpm-1], pull_lph_HX_sigma[rap-1][pt-1][cpm-1], pull_lph_HX_sigmaerr[rap-1][pt-1][cpm-1], filename, nGenerations,true,0);

				sprintf(filename,"%s/Figures/Pull_ltp_HX_rap%dpt%dcpm%d.pdf",dirstruct,iRap,iPt,icpm);
				PlotObject(pull_ltp_HX, "z[#lambda^{HX}_{#theta#phi}]", l_min_pull, l_max_pull, l_step_1D_pull, pull_ltp_HX_mean[rap-1][pt-1][cpm-1], pull_ltp_HX_meanerr[rap-1][pt-1][cpm-1], pull_ltp_HX_sigma[rap-1][pt-1][cpm-1], pull_ltp_HX_sigmaerr[rap-1][pt-1][cpm-1], filename, nGenerations,true,0);

				sprintf(filename,"%s/Figures/Pull_lthstar_HX_rap%dpt%dcpm%d.pdf",dirstruct,iRap,iPt,icpm);
				PlotObject(pull_lthstar_HX, "z[#lambda^{*HX}_{#theta}]", l_min_pull, l_max_pull, l_step_1D_pull, pull_lthstar_HX_mean[rap-1][pt-1][cpm-1], pull_lthstar_HX_meanerr[rap-1][pt-1][cpm-1], pull_lthstar_HX_sigma[rap-1][pt-1][cpm-1], pull_lthstar_HX_sigmaerr[rap-1][pt-1][cpm-1], filename, nGenerations,true,0);

				sprintf(filename,"%s/Figures/Pull_lphstar_HX_rap%dpt%dcpm%d.pdf",dirstruct,iRap,iPt,icpm);
				PlotObject(pull_lphstar_HX, "z[#lambda^{*HX}_{#phi}]", l_min_pull, l_max_pull, l_step_1D_pull, pull_lphstar_HX_mean[rap-1][pt-1][cpm-1], pull_lphstar_HX_meanerr[rap-1][pt-1][cpm-1], pull_lphstar_HX_sigma[rap-1][pt-1][cpm-1], pull_lphstar_HX_sigmaerr[rap-1][pt-1][cpm-1], filename, nGenerations,true,0);

				sprintf(filename,"%s/Figures/Pull_ltilde_HX_rap%dpt%dcpm%d.pdf",dirstruct,iRap,iPt,icpm);
				PlotObject(pull_ltilde_HX, "z[#tilde{#lambda}^{HX}]", l_min_pull, l_max_pull, l_step_1D_pull, pull_ltilde_HX_mean[rap-1][pt-1][cpm-1], pull_ltilde_HX_meanerr[rap-1][pt-1][cpm-1], pull_ltilde_HX_sigma[rap-1][pt-1][cpm-1], pull_ltilde_HX_sigmaerr[rap-1][pt-1][cpm-1], filename, nGenerations,true,0);


				sprintf(filename,"%s/Figures/Pull_lth_PX_rap%dpt%dcpm%d.pdf",dirstruct,iRap,iPt,icpm);
				PlotObject(pull_lth_PX, "z[#lambda^{PX}_{#theta}]", l_min_pull, l_max_pull, l_step_1D_pull, pull_lth_PX_mean[rap-1][pt-1][cpm-1], pull_lth_PX_meanerr[rap-1][pt-1][cpm-1], pull_lth_PX_sigma[rap-1][pt-1][cpm-1], pull_lth_PX_sigmaerr[rap-1][pt-1][cpm-1], filename, nGenerations,true,0);

				sprintf(filename,"%s/Figures/Pull_lph_PX_rap%dpt%dcpm%d.pdf",dirstruct,iRap,iPt,icpm);
				PlotObject(pull_lph_PX, "z[#lambda^{PX}_{#phi}]", l_min_pull, l_max_pull, l_step_1D_pull, pull_lph_PX_mean[rap-1][pt-1][cpm-1], pull_lph_PX_meanerr[rap-1][pt-1][cpm-1], pull_lph_PX_sigma[rap-1][pt-1][cpm-1], pull_lph_PX_sigmaerr[rap-1][pt-1][cpm-1], filename, nGenerations,true,0);

				sprintf(filename,"%s/Figures/Pull_ltp_PX_rap%dpt%dcpm%d.pdf",dirstruct,iRap,iPt,icpm);
				PlotObject(pull_ltp_PX, "z[#lambda^{PX}_{#theta#phi}]", l_min_pull, l_max_pull, l_step_1D_pull, pull_ltp_PX_mean[rap-1][pt-1][cpm-1], pull_ltp_PX_meanerr[rap-1][pt-1][cpm-1], pull_ltp_PX_sigma[rap-1][pt-1][cpm-1], pull_ltp_PX_sigmaerr[rap-1][pt-1][cpm-1], filename, nGenerations,true,0);

				sprintf(filename,"%s/Figures/Pull_lthstar_PX_rap%dpt%dcpm%d.pdf",dirstruct,iRap,iPt,icpm);
				PlotObject(pull_lthstar_PX, "z[#lambda^{*PX}_{#theta}]", l_min_pull, l_max_pull, l_step_1D_pull, pull_lthstar_PX_mean[rap-1][pt-1][cpm-1], pull_lthstar_PX_meanerr[rap-1][pt-1][cpm-1], pull_lthstar_PX_sigma[rap-1][pt-1][cpm-1], pull_lthstar_PX_sigmaerr[rap-1][pt-1][cpm-1], filename, nGenerations,true,0);

				sprintf(filename,"%s/Figures/Pull_lphstar_PX_rap%dpt%dcpm%d.pdf",dirstruct,iRap,iPt,icpm);
				PlotObject(pull_lphstar_PX, "z[#lambda^{*PX}_{#phi}]", l_min_pull, l_max_pull, l_step_1D_pull, pull_lphstar_PX_mean[rap-1][pt-1][cpm-1], pull_lphstar_PX_meanerr[rap-1][pt-1][cpm-1], pull_lphstar_PX_sigma[rap-1][pt-1][cpm-1], pull_lphstar_PX_sigmaerr[rap-1][pt-1][cpm-1], filename, nGenerations,true,0);

				sprintf(filename,"%s/Figures/Pull_ltilde_PX_rap%dpt%dcpm%d.pdf",dirstruct,iRap,iPt,icpm);
				PlotObject(pull_ltilde_PX, "z[#tilde{#lambda}^{PX}]", l_min_pull, l_max_pull, l_step_1D_pull, pull_ltilde_PX_mean[rap-1][pt-1][cpm-1], pull_ltilde_PX_meanerr[rap-1][pt-1][cpm-1], pull_ltilde_PX_sigma[rap-1][pt-1][cpm-1], pull_ltilde_PX_sigmaerr[rap-1][pt-1][cpm-1], filename, nGenerations,true,0);




				sprintf(filename,"%s/Figures/Param_lth_CS_rap%dpt%dcpm%d.pdf",dirstruct,iRap,iPt,icpm);
				PlotObject(param_lth_CS, "#lambda^{CS}_{#theta}", l_min_lth_CS, l_max_lth_CS, l_step_1D, param_lth_CS_mean[rap-1][pt-1][cpm-1], param_lth_CS_meanerr[rap-1][pt-1][cpm-1], param_lth_CS_sigma[rap-1][pt-1][cpm-1], param_lth_CS_sigmaerr[rap-1][pt-1][cpm-1], filename, nGenerations,false,lambda_theta_injected_CS[rap-1][pt-1][cpm-1]);

				sprintf(filename,"%s/Figures/Param_lph_CS_rap%dpt%dcpm%d.pdf",dirstruct,iRap,iPt,icpm);
				PlotObject(param_lph_CS, "#lambda^{CS}_{#phi}", l_min_lph_CS, l_max_lph_CS, l_step_1D, param_lph_CS_mean[rap-1][pt-1][cpm-1], param_lph_CS_meanerr[rap-1][pt-1][cpm-1], param_lph_CS_sigma[rap-1][pt-1][cpm-1], param_lph_CS_sigmaerr[rap-1][pt-1][cpm-1], filename, nGenerations,false,lambda_phi_injected_CS[rap-1][pt-1][cpm-1]);

				sprintf(filename,"%s/Figures/Param_ltp_CS_rap%dpt%dcpm%d.pdf",dirstruct,iRap,iPt,icpm);
				PlotObject(param_ltp_CS, "#lambda^{CS}_{#theta#phi}", l_min_ltp_CS, l_max_ltp_CS, l_step_1D, param_ltp_CS_mean[rap-1][pt-1][cpm-1], param_ltp_CS_meanerr[rap-1][pt-1][cpm-1], param_ltp_CS_sigma[rap-1][pt-1][cpm-1], param_ltp_CS_sigmaerr[rap-1][pt-1][cpm-1], filename, nGenerations,false,lambda_thetaphi_injected_CS[rap-1][pt-1][cpm-1]);

				sprintf(filename,"%s/Figures/Param_lthstar_CS_rap%dpt%dcpm%d.pdf",dirstruct,iRap,iPt,icpm);
				PlotObject(param_lthstar_CS, "#lambda^{*CS}_{#theta}", l_min_lthstar_CS, l_max_lthstar_CS, l_step_1D, param_lthstar_CS_mean[rap-1][pt-1][cpm-1], param_lthstar_CS_meanerr[rap-1][pt-1][cpm-1], param_lthstar_CS_sigma[rap-1][pt-1][cpm-1], param_lthstar_CS_sigmaerr[rap-1][pt-1][cpm-1], filename, nGenerations,false,lambda_thetastar_injected_CS[rap-1][pt-1][cpm-1]);

				sprintf(filename,"%s/Figures/Param_lphstar_CS_rap%dpt%dcpm%d.pdf",dirstruct,iRap,iPt,icpm);
				PlotObject(param_lphstar_CS, "#lambda^{*CS}_{#phi}", l_min_lphstar_CS, l_max_lphstar_CS, l_step_1D, param_lphstar_CS_mean[rap-1][pt-1][cpm-1], param_lphstar_CS_meanerr[rap-1][pt-1][cpm-1], param_lphstar_CS_sigma[rap-1][pt-1][cpm-1], param_lphstar_CS_sigmaerr[rap-1][pt-1][cpm-1], filename, nGenerations,false,lambda_phistar_injected_CS[rap-1][pt-1][cpm-1]);

				sprintf(filename,"%s/Figures/Param_ltilde_CS_rap%dpt%dcpm%d.pdf",dirstruct,iRap,iPt,icpm);
				PlotObject(param_ltilde_CS, "#tilde{#lambda}^{CS}", l_min_ltilde_CS, l_max_ltilde_CS, l_step_1D, param_ltilde_CS_mean[rap-1][pt-1][cpm-1], param_ltilde_CS_meanerr[rap-1][pt-1][cpm-1], param_ltilde_CS_sigma[rap-1][pt-1][cpm-1], param_ltilde_CS_sigmaerr[rap-1][pt-1][cpm-1], filename, nGenerations,false,lambda_tilde_injected_CS[rap-1][pt-1][cpm-1]);


				sprintf(filename,"%s/Figures/Param_lth_HX_rap%dpt%dcpm%d.pdf",dirstruct,iRap,iPt,icpm);
				PlotObject(param_lth_HX, "#lambda^{HX}_{#theta}", l_min_lth_HX, l_max_lth_HX, l_step_1D, param_lth_HX_mean[rap-1][pt-1][cpm-1], param_lth_HX_meanerr[rap-1][pt-1][cpm-1], param_lth_HX_sigma[rap-1][pt-1][cpm-1], param_lth_HX_sigmaerr[rap-1][pt-1][cpm-1], filename, nGenerations,false,lambda_theta_injected_HX[rap-1][pt-1][cpm-1]);

				sprintf(filename,"%s/Figures/Param_lph_HX_rap%dpt%dcpm%d.pdf",dirstruct,iRap,iPt,icpm);
				PlotObject(param_lph_HX, "#lambda^{HX}_{#phi}", l_min_lph_HX, l_max_lph_HX, l_step_1D, param_lph_HX_mean[rap-1][pt-1][cpm-1], param_lph_HX_meanerr[rap-1][pt-1][cpm-1], param_lph_HX_sigma[rap-1][pt-1][cpm-1], param_lph_HX_sigmaerr[rap-1][pt-1][cpm-1], filename, nGenerations,false,lambda_phi_injected_HX[rap-1][pt-1][cpm-1]);

				sprintf(filename,"%s/Figures/Param_ltp_HX_rap%dpt%dcpm%d.pdf",dirstruct,iRap,iPt,icpm);
				PlotObject(param_ltp_HX, "#lambda^{HX}_{#theta#phi}", l_min_ltp_HX, l_max_ltp_HX, l_step_1D, param_ltp_HX_mean[rap-1][pt-1][cpm-1], param_ltp_HX_meanerr[rap-1][pt-1][cpm-1], param_ltp_HX_sigma[rap-1][pt-1][cpm-1], param_ltp_HX_sigmaerr[rap-1][pt-1][cpm-1], filename, nGenerations,false,lambda_thetaphi_injected_HX[rap-1][pt-1][cpm-1]);

				sprintf(filename,"%s/Figures/Param_lthstar_HX_rap%dpt%dcpm%d.pdf",dirstruct,iRap,iPt,icpm);
				PlotObject(param_lthstar_HX, "#lambda^{*HX}_{#theta}", l_min_lthstar_HX, l_max_lthstar_HX, l_step_1D, param_lthstar_HX_mean[rap-1][pt-1][cpm-1], param_lthstar_HX_meanerr[rap-1][pt-1][cpm-1], param_lthstar_HX_sigma[rap-1][pt-1][cpm-1], param_lthstar_HX_sigmaerr[rap-1][pt-1][cpm-1], filename, nGenerations,false,lambda_thetastar_injected_HX[rap-1][pt-1][cpm-1]);

				sprintf(filename,"%s/Figures/Param_lphstar_HX_rap%dpt%dcpm%d.pdf",dirstruct,iRap,iPt,icpm);
				PlotObject(param_lphstar_HX, "#lambda^{*HX}_{#phi}", l_min_lphstar_HX, l_max_lphstar_HX, l_step_1D, param_lphstar_HX_mean[rap-1][pt-1][cpm-1], param_lphstar_HX_meanerr[rap-1][pt-1][cpm-1], param_lphstar_HX_sigma[rap-1][pt-1][cpm-1], param_lphstar_HX_sigmaerr[rap-1][pt-1][cpm-1], filename, nGenerations,false,lambda_phistar_injected_HX[rap-1][pt-1][cpm-1]);

				sprintf(filename,"%s/Figures/Param_ltilde_HX_rap%dpt%dcpm%d.pdf",dirstruct,iRap,iPt,icpm);
				PlotObject(param_ltilde_HX, "#tilde{#lambda}^{HX}", l_min_ltilde_HX, l_max_ltilde_HX, l_step_1D, param_ltilde_HX_mean[rap-1][pt-1][cpm-1], param_ltilde_HX_meanerr[rap-1][pt-1][cpm-1], param_ltilde_HX_sigma[rap-1][pt-1][cpm-1], param_ltilde_HX_sigmaerr[rap-1][pt-1][cpm-1], filename, nGenerations,false,lambda_tilde_injected_HX[rap-1][pt-1][cpm-1]);


				sprintf(filename,"%s/Figures/Param_lth_PX_rap%dpt%dcpm%d.pdf",dirstruct,iRap,iPt,icpm);
				PlotObject(param_lth_PX, "#lambda^{PX}_{#theta}", l_min_lth_PX, l_max_lth_PX, l_step_1D, param_lth_PX_mean[rap-1][pt-1][cpm-1], param_lth_PX_meanerr[rap-1][pt-1][cpm-1], param_lth_PX_sigma[rap-1][pt-1][cpm-1], param_lth_PX_sigmaerr[rap-1][pt-1][cpm-1], filename, nGenerations,false,lambda_theta_injected_PX[rap-1][pt-1][cpm-1]);

				sprintf(filename,"%s/Figures/Param_lph_PX_rap%dpt%dcpm%d.pdf",dirstruct,iRap,iPt,icpm);
				PlotObject(param_lph_PX, "#lambda^{PX}_{#phi}", l_min_lph_PX, l_max_lph_PX, l_step_1D, param_lph_PX_mean[rap-1][pt-1][cpm-1], param_lph_PX_meanerr[rap-1][pt-1][cpm-1], param_lph_PX_sigma[rap-1][pt-1][cpm-1], param_lph_PX_sigmaerr[rap-1][pt-1][cpm-1], filename, nGenerations,false,lambda_phi_injected_PX[rap-1][pt-1][cpm-1]);

				sprintf(filename,"%s/Figures/Param_ltp_PX_rap%dpt%dcpm%d.pdf",dirstruct,iRap,iPt,icpm);
				PlotObject(param_ltp_PX, "#lambda^{PX}_{#theta#phi}", l_min_ltp_PX, l_max_ltp_PX, l_step_1D, param_ltp_PX_mean[rap-1][pt-1][cpm-1], param_ltp_PX_meanerr[rap-1][pt-1][cpm-1], param_ltp_PX_sigma[rap-1][pt-1][cpm-1], param_ltp_PX_sigmaerr[rap-1][pt-1][cpm-1], filename, nGenerations,false,lambda_thetaphi_injected_PX[rap-1][pt-1][cpm-1]);

				sprintf(filename,"%s/Figures/Param_lthstar_PX_rap%dpt%dcpm%d.pdf",dirstruct,iRap,iPt,icpm);
				PlotObject(param_lthstar_PX, "#lambda^{*PX}_{#theta}", l_min_lthstar_PX, l_max_lthstar_PX, l_step_1D, param_lthstar_PX_mean[rap-1][pt-1][cpm-1], param_lthstar_PX_meanerr[rap-1][pt-1][cpm-1], param_lthstar_PX_sigma[rap-1][pt-1][cpm-1], param_lthstar_PX_sigmaerr[rap-1][pt-1][cpm-1], filename, nGenerations,false,lambda_thetastar_injected_PX[rap-1][pt-1][cpm-1]);

				sprintf(filename,"%s/Figures/Param_lphstar_PX_rap%dpt%dcpm%d.pdf",dirstruct,iRap,iPt,icpm);
				PlotObject(param_lphstar_PX, "#lambda^{*PX}_{#phi}", l_min_lphstar_PX, l_max_lphstar_PX, l_step_1D, param_lphstar_PX_mean[rap-1][pt-1][cpm-1], param_lphstar_PX_meanerr[rap-1][pt-1][cpm-1], param_lphstar_PX_sigma[rap-1][pt-1][cpm-1], param_lphstar_PX_sigmaerr[rap-1][pt-1][cpm-1], filename, nGenerations,false,lambda_phistar_injected_PX[rap-1][pt-1][cpm-1]);

				sprintf(filename,"%s/Figures/Param_ltilde_PX_rap%dpt%dcpm%d.pdf",dirstruct,iRap,iPt,icpm);
				PlotObject(param_ltilde_PX, "#tilde{#lambda}^{PX}", l_min_ltilde_PX, l_max_ltilde_PX, l_step_1D, param_ltilde_PX_mean[rap-1][pt-1][cpm-1], param_ltilde_PX_meanerr[rap-1][pt-1][cpm-1], param_ltilde_PX_sigma[rap-1][pt-1][cpm-1], param_ltilde_PX_sigmaerr[rap-1][pt-1][cpm-1], filename, nGenerations,false,lambda_tilde_injected_PX[rap-1][pt-1][cpm-1]);


			}

			delete param_lth_CS ;
			delete param_lph_CS ;
			delete param_ltp_CS ;
			delete param_ltilde_CS ;
			delete param_lphstar_CS ;
			delete param_lthstar_CS  ;
			delete param_lth_HX      ;
			delete param_lph_HX      ;
			delete param_ltp_HX      ;
			delete param_ltilde_HX   ;
			delete param_lphstar_HX  ;
			delete param_lthstar_HX  ;
			delete param_lth_PX      ;
			delete param_lph_PX      ;
			delete param_ltp_PX      ;
			delete param_ltilde_PX   ;
			delete param_lphstar_PX  ;
			delete param_lthstar_PX  ;
			delete pull_lth_CS         ;
			delete pull_lph_CS         ;
			delete pull_ltp_CS         ;
			delete pull_ltilde_CS      ;
			delete pull_lphstar_CS     ;
			delete pull_lthstar_CS     ;
			delete pull_lth_HX         ;
			delete pull_lph_HX         ;
			delete pull_ltp_HX         ;
			delete pull_ltilde_HX      ;
			delete pull_lphstar_HX     ;
			delete pull_lthstar_HX     ;
			delete pull_lth_PX         ;
			delete pull_lph_PX         ;
			delete pull_ltp_PX         ;
			delete pull_ltilde_PX      ;
			delete pull_lphstar_PX     ;
			delete pull_lthstar_PX     ;


		}}}

		//////// Rap/Pt Summary Plots ////////////////


		double lammin;
		double lammax;

		double pmumin=-5.1;
		double pmumax=5.1;
		double psigmin=0;
		double psigmax=2;

		double dsigmin=0;
		double dsigmax=0.25;
		bool DrawLine(false);
		double lamline=999;

		double ddeltamumin=-delta_l;
		double ddeltamumax=delta_l;

		//////////////////////////////////////////////////

		double dmumin=-1.1;
		double dmumax=1.1;
		double thetaphimin=-0.25;
		double thetaphimax=0.25;
		double tildeCompmin=dmumin;
		double tildeCompmax=dmumax;
		/*
			 double dmumin=-2;
			 double dmumax=10;
			 double thetaphimin=-2;
			 double thetaphimax=3;
			 double tildeCompmin=-2;
			 double tildeCompmax=30;
		 */
		//////////////////////////////////////////////////

		double KinDepRap[NptBins][NcpmBins];
		double KinDepRaperr_low[NptBins][NcpmBins];
		double KinDepRaperr_high[NptBins][NcpmBins];
		double KinDepRap2[NptBins][NcpmBins];
		double KinDepRap2err_low[NptBins][NcpmBins];
		double KinDepRap2err_high[NptBins][NcpmBins];
		double KinDepRap3[NptBins][NcpmBins];
		double KinDepRap3err_low[NptBins][NcpmBins];
		double KinDepRap3err_high[NptBins][NcpmBins];

		int nLamPlots=105;
		char axislabel[200];
		char GraphName[200];
		char GraphFolder[200];
		bool SaveGraph(false);
		int nFrame=0;

		sprintf(GraphFolder,"%s",dirstruct);

		if(RealData){
			sprintf(dirstruct,"%s/Figures/%s",dirstruct,TreeID);
		}


		for(int iLam=1; iLam<nLamPlots+1;iLam++){

			DrawLine=false;
			SaveGraph=false;

			if(iLam<13)continue;


			int rap=0;
			for(int rapBin = rapBinMin; rapBin < rapBinMax+1; rapBin++) {
			int pt=0;
				for(int ptBin = ptBinMin; ptBin < ptBinMax+1; ptBin++) {

				if(iLam==13){lammin=dmumin;		lammax=dmumax;	sprintf(filename,"%s/Figures/KinDep_param_CS_mean_lth_rap%d_pt%d.pdf",dirstruct,rapBin,ptBin); 		sprintf(axislabel,"#lambda^{CS}_{#theta}");																sprintf(GraphName,"lth_CS_rap%d_pt%d",rapBin,ptBin); if(RealData) SaveGraph=true;		}
				if(iLam==14){lammin=dsigmin;	lammax=dsigmax;	sprintf(filename,"%s/Figures/KinDep_param_CS_sigma_lth_rap%d_pt%d.pdf",dirstruct,rapBin,ptBin); 	sprintf(axislabel,"#sigma_{#lambda^{CS}_{#theta}}");																}
				if(iLam==15){lammin=thetaphimin;		lammax=thetaphimax;	sprintf(filename,"%s/Figures/KinDep_param_CS_mean_lph_rap%d_pt%d.pdf",dirstruct,rapBin,ptBin); 		sprintf(axislabel,"#lambda^{CS}_{#phi}");																sprintf(GraphName,"lph_CS_rap%d_pt%d",rapBin,ptBin); if(RealData) SaveGraph=true;		}
				if(iLam==16){lammin=dsigmin;	lammax=dsigmax;	sprintf(filename,"%s/Figures/KinDep_param_CS_sigma_lph_rap%d_pt%d.pdf",dirstruct,rapBin,ptBin); 	sprintf(axislabel,"#sigma_{#lambda^{CS}_{#phi}}");																	}
				if(iLam==17){lammin=thetaphimin;		lammax=thetaphimax;	sprintf(filename,"%s/Figures/KinDep_param_CS_mean_ltp_rap%d_pt%d.pdf",dirstruct,rapBin,ptBin); 		sprintf(axislabel,"#lambda^{CS}_{#theta#phi}");															sprintf(GraphName,"ltp_CS_rap%d_pt%d",rapBin,ptBin); if(RealData) SaveGraph=true;		}
				if(iLam==18){lammin=dsigmin;	lammax=dsigmax;	sprintf(filename,"%s/Figures/KinDep_param_CS_sigma_ltp_rap%d_pt%d.pdf",dirstruct,rapBin,ptBin); 	sprintf(axislabel,"#sigma_{#lambda^{CS}_{#theta#phi}}");															}
				if(iLam==19){lammin=dmumin;		lammax=dmumax;	sprintf(filename,"%s/Figures/KinDep_param_CS_mean_lthstar_rap%d_pt%d.pdf",dirstruct,rapBin,ptBin); 		sprintf(axislabel,"#lambda^{*CS}_{#theta}");														sprintf(GraphName,"lthstar_CS_rap%d_pt%d",rapBin,ptBin); if(RealData) SaveGraph=true;		}
				if(iLam==20){lammin=dsigmin;	lammax=dsigmax;	sprintf(filename,"%s/Figures/KinDep_param_CS_sigma_lthstar_rap%d_pt%d.pdf",dirstruct,rapBin,ptBin); 	sprintf(axislabel,"#sigma_{#lambda^{*CS}_{#theta}}");															}
				if(iLam==21){lammin=dmumin;		lammax=dmumax;	sprintf(filename,"%s/Figures/KinDep_param_CS_mean_lphstar_rap%d_pt%d.pdf",dirstruct,rapBin,ptBin); 		sprintf(axislabel,"#lambda^{*CS}_{#phi}");															sprintf(GraphName,"lphstar_CS_rap%d_pt%d",rapBin,ptBin); if(RealData) SaveGraph=true;		}
				if(iLam==22){lammin=dsigmin;	lammax=dsigmax;	sprintf(filename,"%s/Figures/KinDep_param_CS_sigma_lphstar_rap%d_pt%d.pdf",dirstruct,rapBin,ptBin); 	sprintf(axislabel,"#sigma_{#lambda^{*CS}_{#phi}}");																}
				if(iLam==23){lammin=dmumin;		lammax=dmumax;	if(polScen==6){lammin=lammin+2;lammax=lammax+2;}sprintf(filename,"%s/Figures/KinDep_param_CS_mean_ltilde_rap%d_pt%d.pdf",dirstruct,rapBin,ptBin); 		sprintf(axislabel,"#tilde{#lambda}^{CS}");				sprintf(GraphName,"ltilde_CS_rap%d_pt%d",rapBin,ptBin); if(RealData) SaveGraph=true;		}
				if(iLam==24){lammin=dsigmin;	lammax=dsigmax;	if(polScen==6){lammin=lammin+2;lammax=lammax+2;}sprintf(filename,"%s/Figures/KinDep_param_CS_sigma_ltilde_rap%d_pt%d.pdf",dirstruct,rapBin,ptBin); 	sprintf(axislabel,"#sigma_{#tilde{#lambda}^{CS}}");				}

				if(iLam==25){lammin=dmumin;		lammax=dmumax;	sprintf(filename,"%s/Figures/KinDep_param_HX_mean_lth_rap%d_pt%d.pdf",dirstruct,rapBin,ptBin); 		sprintf(axislabel,"#lambda^{HX}_{#theta}");                                                            sprintf(GraphName,"lth_HX_rap%d_pt%d",rapBin,ptBin); if(RealData) SaveGraph=true;		}
				if(iLam==26){lammin=dsigmin;	lammax=dsigmax;	sprintf(filename,"%s/Figures/KinDep_param_HX_sigma_lth_rap%d_pt%d.pdf",dirstruct,rapBin,ptBin); 	sprintf(axislabel,"#sigma_{#lambda^{HX}_{#theta}}");                                                       		}
				if(iLam==27){lammin=thetaphimin;		lammax=thetaphimax;	sprintf(filename,"%s/Figures/KinDep_param_HX_mean_lph_rap%d_pt%d.pdf",dirstruct,rapBin,ptBin); 		sprintf(axislabel,"#lambda^{HX}_{#phi}");                                                              sprintf(GraphName,"lph_HX_rap%d_pt%d",rapBin,ptBin); if(RealData) SaveGraph=true;		}
				if(iLam==28){lammin=dsigmin;	lammax=dsigmax;	sprintf(filename,"%s/Figures/KinDep_param_HX_sigma_lph_rap%d_pt%d.pdf",dirstruct,rapBin,ptBin); 	sprintf(axislabel,"#sigma_{#lambda^{HX}_{#phi}}");                                                         		}
				if(iLam==29){lammin=thetaphimin;		lammax=thetaphimax;	sprintf(filename,"%s/Figures/KinDep_param_HX_mean_ltp_rap%d_pt%d.pdf",dirstruct,rapBin,ptBin); 		sprintf(axislabel,"#lambda^{HX}_{#theta#phi}");                                                        sprintf(GraphName,"ltp_HX_rap%d_pt%d",rapBin,ptBin); if(RealData) SaveGraph=true;		}
				if(iLam==30){lammin=dsigmin;	lammax=dsigmax;	sprintf(filename,"%s/Figures/KinDep_param_HX_sigma_ltp_rap%d_pt%d.pdf",dirstruct,rapBin,ptBin); 	sprintf(axislabel,"#sigma_{#lambda^{HX}_{#theta#phi}}");                                                   		}
				if(iLam==31){lammin=dmumin;		lammax=dmumax;	sprintf(filename,"%s/Figures/KinDep_param_HX_mean_lthstar_rap%d_pt%d.pdf",dirstruct,rapBin,ptBin); 		sprintf(axislabel,"#lambda^{*HX}_{#theta}");                                                       sprintf(GraphName,"lthstar_HX_rap%d_pt%d",rapBin,ptBin); if(RealData) SaveGraph=true;		}
				if(iLam==32){lammin=dsigmin;	lammax=dsigmax;	sprintf(filename,"%s/Figures/KinDep_param_HX_sigma_lthstar_rap%d_pt%d.pdf",dirstruct,rapBin,ptBin); 	sprintf(axislabel,"#sigma_{#lambda^{*HX}_{#theta}}");                                                  		}
				if(iLam==33){lammin=dmumin;		lammax=dmumax;	sprintf(filename,"%s/Figures/KinDep_param_HX_mean_lphstar_rap%d_pt%d.pdf",dirstruct,rapBin,ptBin); 		sprintf(axislabel,"#lambda^{*HX}_{#phi}");                                                         sprintf(GraphName,"lphstar_HX_rap%d_pt%d",rapBin,ptBin); if(RealData) SaveGraph=true;		}
				if(iLam==34){lammin=dsigmin;	lammax=dsigmax;	sprintf(filename,"%s/Figures/KinDep_param_HX_sigma_lphstar_rap%d_pt%d.pdf",dirstruct,rapBin,ptBin); 	sprintf(axislabel,"#sigma_{#lambda^{*HX}_{#phi}}");                                                    		}
				if(iLam==35){lammin=dmumin;		lammax=dmumax;	if(polScen==6){lammin=lammin+2;lammax=lammax+2;}sprintf(filename,"%s/Figures/KinDep_param_HX_mean_ltilde_rap%d_pt%d.pdf",dirstruct,rapBin,ptBin); 		sprintf(axislabel,"#tilde{#lambda}^{HX}");             sprintf(GraphName,"ltilde_HX_rap%d_pt%d",rapBin,ptBin); if(RealData) SaveGraph=true;		}
				if(iLam==36){lammin=dsigmin;	lammax=dsigmax;	if(polScen==6){lammin=lammin+2;lammax=lammax+2;}sprintf(filename,"%s/Figures/KinDep_param_HX_sigma_ltilde_rap%d_pt%d.pdf",dirstruct,rapBin,ptBin); 	sprintf(axislabel,"#sigma_{#tilde{#lambda}^{HX}}");    		}

				if(iLam==37){lammin=dmumin;		lammax=dmumax;	sprintf(filename,"%s/Figures/KinDep_param_PX_mean_lth_rap%d_pt%d.pdf",dirstruct,rapBin,ptBin); 		sprintf(axislabel,"#lambda^{PX}_{#theta}");                                                            sprintf(GraphName,"lth_PX_rap%d_pt%d",rapBin,ptBin); if(RealData) SaveGraph=true;		}
				if(iLam==38){lammin=dsigmin;	lammax=dsigmax;	sprintf(filename,"%s/Figures/KinDep_param_PX_sigma_lth_rap%d_pt%d.pdf",dirstruct,rapBin,ptBin); 	sprintf(axislabel,"#sigma_{#lambda^{PX}_{#theta}}");                                                       		}
				if(iLam==39){lammin=thetaphimin;		lammax=thetaphimax;	sprintf(filename,"%s/Figures/KinDep_param_PX_mean_lph_rap%d_pt%d.pdf",dirstruct,rapBin,ptBin); 		sprintf(axislabel,"#lambda^{PX}_{#phi}");                                                              sprintf(GraphName,"lph_PX_rap%d_pt%d",rapBin,ptBin); if(RealData) SaveGraph=true;		}
				if(iLam==40){lammin=dsigmin;	lammax=dsigmax;	sprintf(filename,"%s/Figures/KinDep_param_PX_sigma_lph_rap%d_pt%d.pdf",dirstruct,rapBin,ptBin); 	sprintf(axislabel,"#sigma_{#lambda^{PX}_{#phi}}");                                                         		}
				if(iLam==41){lammin=thetaphimin;		lammax=thetaphimax;	sprintf(filename,"%s/Figures/KinDep_param_PX_mean_ltp_rap%d_pt%d.pdf",dirstruct,rapBin,ptBin); 		sprintf(axislabel,"#lambda^{PX}_{#theta#phi}");                                                        sprintf(GraphName,"ltp_PX_rap%d_pt%d",rapBin,ptBin); if(RealData) SaveGraph=true;		}
				if(iLam==42){lammin=dsigmin;	lammax=dsigmax;	sprintf(filename,"%s/Figures/KinDep_param_PX_sigma_ltp_rap%d_pt%d.pdf",dirstruct,rapBin,ptBin); 	sprintf(axislabel,"#sigma_{#lambda^{PX}_{#theta#phi}}");                                                   		}
				if(iLam==43){lammin=dmumin;		lammax=dmumax;	sprintf(filename,"%s/Figures/KinDep_param_PX_mean_lthstar_rap%d_pt%d.pdf",dirstruct,rapBin,ptBin); 		sprintf(axislabel,"#lambda^{*PX}_{#theta}");                                                       sprintf(GraphName,"lthstar_PX_rap%d_pt%d",rapBin,ptBin); if(RealData) SaveGraph=true;		}
				if(iLam==44){lammin=dsigmin;	lammax=dsigmax;	sprintf(filename,"%s/Figures/KinDep_param_PX_sigma_lthstar_rap%d_pt%d.pdf",dirstruct,rapBin,ptBin); 	sprintf(axislabel,"#sigma_{#lambda^{*PX}_{#theta}}");                                                  		}
				if(iLam==45){lammin=dmumin;		lammax=dmumax;	sprintf(filename,"%s/Figures/KinDep_param_PX_mean_lphstar_rap%d_pt%d.pdf",dirstruct,rapBin,ptBin); 		sprintf(axislabel,"#lambda^{*PX}_{#phi}");                                                         sprintf(GraphName,"lphstar_PX_rap%d_pt%d",rapBin,ptBin); if(RealData) SaveGraph=true;		}
				if(iLam==46){lammin=dsigmin;	lammax=dsigmax;	sprintf(filename,"%s/Figures/KinDep_param_PX_sigma_lphstar_rap%d_pt%d.pdf",dirstruct,rapBin,ptBin); 	sprintf(axislabel,"#sigma_{#lambda^{*PX}_{#phi}}");                                                    		}
				if(iLam==47){lammin=dmumin;		lammax=dmumax;	if(polScen==6){lammin=lammin+2;lammax=lammax+2;}sprintf(filename,"%s/Figures/KinDep_param_PX_mean_ltilde_rap%d_pt%d.pdf",dirstruct,rapBin,ptBin); 		sprintf(axislabel,"#tilde{#lambda}^{PX}");             sprintf(GraphName,"ltilde_PX_rap%d_pt%d",rapBin,ptBin); if(RealData) SaveGraph=true;		}
				if(iLam==48){lammin=dsigmin;	lammax=dsigmax;	if(polScen==6){lammin=lammin+2;lammax=lammax+2;}sprintf(filename,"%s/Figures/KinDep_param_PX_sigma_ltilde_rap%d_pt%d.pdf",dirstruct,rapBin,ptBin); 	sprintf(axislabel,"#sigma_{#tilde{#lambda}^{PX}}");    		}

				if(iLam==49){lammin=pmumin;		lammax=pmumax;	sprintf(filename,"%s/Figures/KinDep_pull_CS_mean_lth_rap%d_pt%d.pdf",dirstruct,rapBin,ptBin); 			sprintf(axislabel,"#mu_{[z[#lambda^{CS}_{#theta}]}");				DrawLine=true;lamline=0;}
				if(iLam==50){lammin=psigmin;	lammax=psigmax;	sprintf(filename,"%s/Figures/KinDep_pull_CS_sigma_lth_rap%d_pt%d.pdf",dirstruct,rapBin,ptBin); 		sprintf(axislabel,"#sigma_{z[#lambda^{CS}_{#theta}]}");				DrawLine=true;lamline=1;}
				if(iLam==51){lammin=pmumin;		lammax=pmumax;	sprintf(filename,"%s/Figures/KinDep_pull_CS_mean_lph_rap%d_pt%d.pdf",dirstruct,rapBin,ptBin); 			sprintf(axislabel,"#mu_{[z[#lambda^{CS}_{#phi}]}");					DrawLine=true;lamline=0;}
				if(iLam==52){lammin=psigmin;	lammax=psigmax;	sprintf(filename,"%s/Figures/KinDep_pull_CS_sigma_lph_rap%d_pt%d.pdf",dirstruct,rapBin,ptBin); 		sprintf(axislabel,"#sigma_{z[#lambda^{CS}_{#phi}]}");				DrawLine=true;lamline=1;}
				if(iLam==53){lammin=pmumin;		lammax=pmumax;	sprintf(filename,"%s/Figures/KinDep_pull_CS_mean_ltp_rap%d_pt%d.pdf",dirstruct,rapBin,ptBin); 			sprintf(axislabel,"#mu_{[z[#lambda^{CS}_{#theta#phi}]}");			DrawLine=true;lamline=0;}
				if(iLam==54){lammin=psigmin;	lammax=psigmax;	sprintf(filename,"%s/Figures/KinDep_pull_CS_sigma_ltp_rap%d_pt%d.pdf",dirstruct,rapBin,ptBin); 		sprintf(axislabel,"#sigma_{z[#lambda^{CS}_{#theta#phi}]}");			DrawLine=true;lamline=1;}
				if(iLam==55){lammin=pmumin;		lammax=pmumax;	sprintf(filename,"%s/Figures/KinDep_pull_CS_mean_lthstar_rap%d_pt%d.pdf",dirstruct,rapBin,ptBin); 		sprintf(axislabel,"#mu_{[z[#lambda^{*CS}_{#theta}]}");				DrawLine=true;lamline=0;}
				if(iLam==56){lammin=psigmin;	lammax=psigmax;	sprintf(filename,"%s/Figures/KinDep_pull_CS_sigma_lthstar_rap%d_pt%d.pdf",dirstruct,rapBin,ptBin); 	sprintf(axislabel,"#sigma_{z[#lambda^{*CS}_{#theta}]}");			DrawLine=true;lamline=1;}
				if(iLam==57){lammin=pmumin;		lammax=pmumax;	sprintf(filename,"%s/Figures/KinDep_pull_CS_mean_lphstar_rap%d_pt%d.pdf",dirstruct,rapBin,ptBin); 		sprintf(axislabel,"#mu_{[z[#lambda^{*CS}_{#phi}]}");				DrawLine=true;lamline=0;}
				if(iLam==58){lammin=psigmin;	lammax=psigmax;	sprintf(filename,"%s/Figures/KinDep_pull_CS_sigma_lphstar_rap%d_pt%d.pdf",dirstruct,rapBin,ptBin); 	sprintf(axislabel,"#sigma_{z[#lambda^{*CS}_{#phi}]}");				DrawLine=true;lamline=1;}
				if(iLam==59){lammin=pmumin;		lammax=pmumax;	sprintf(filename,"%s/Figures/KinDep_pull_CS_mean_ltilde_rap%d_pt%d.pdf",dirstruct,rapBin,ptBin); 		sprintf(axislabel,"#mu_{[z[#tilde{#lambda}^{CS}]}");				DrawLine=true;lamline=0;}
				if(iLam==60){lammin=psigmin;	lammax=psigmax;	sprintf(filename,"%s/Figures/KinDep_pull_CS_sigma_ltilde_rap%d_pt%d.pdf",dirstruct,rapBin,ptBin); 		sprintf(axislabel,"#sigma_{z[#tilde{#lambda}^{CS}]}");				DrawLine=true;lamline=1;}

				if(iLam==61){lammin=pmumin;		lammax=pmumax;	sprintf(filename,"%s/Figures/KinDep_pull_HX_mean_lth_rap%d_pt%d.pdf",dirstruct,rapBin,ptBin); 			sprintf(axislabel,"#mu_{[z[#lambda^{HX}_{#theta}]}");				DrawLine=true;lamline=0;}
				if(iLam==62){lammin=psigmin;	lammax=psigmax;	sprintf(filename,"%s/Figures/KinDep_pull_HX_sigma_lth_rap%d_pt%d.pdf",dirstruct,rapBin,ptBin); 		sprintf(axislabel,"#sigma_{z[#lambda^{HX}_{#theta}]}");				DrawLine=true;lamline=1;}
				if(iLam==63){lammin=pmumin;		lammax=pmumax;	sprintf(filename,"%s/Figures/KinDep_pull_HX_mean_lph_rap%d_pt%d.pdf",dirstruct,rapBin,ptBin); 			sprintf(axislabel,"#mu_{[z[#lambda^{HX}_{#phi}]}");					DrawLine=true;lamline=0;}
				if(iLam==64){lammin=psigmin;	lammax=psigmax;	sprintf(filename,"%s/Figures/KinDep_pull_HX_sigma_lph_rap%d_pt%d.pdf",dirstruct,rapBin,ptBin); 		sprintf(axislabel,"#sigma_{z[#lambda^{HX}_{#phi}]}");				DrawLine=true;lamline=1;}
				if(iLam==65){lammin=pmumin;		lammax=pmumax;	sprintf(filename,"%s/Figures/KinDep_pull_HX_mean_ltp_rap%d_pt%d.pdf",dirstruct,rapBin,ptBin); 			sprintf(axislabel,"#mu_{[z[#lambda^{HX}_{#theta#phi}]}");			DrawLine=true;lamline=0;}
				if(iLam==66){lammin=psigmin;	lammax=psigmax;	sprintf(filename,"%s/Figures/KinDep_pull_HX_sigma_ltp_rap%d_pt%d.pdf",dirstruct,rapBin,ptBin); 		sprintf(axislabel,"#sigma_{z[#lambda^{HX}_{#theta#phi}]}");			DrawLine=true;lamline=1;}
				if(iLam==67){lammin=pmumin;		lammax=pmumax;	sprintf(filename,"%s/Figures/KinDep_pull_HX_mean_lthstar_rap%d_pt%d.pdf",dirstruct,rapBin,ptBin); 		sprintf(axislabel,"#mu_{[z[#lambda^{*HX}_{#theta}]}");				DrawLine=true;lamline=0;}
				if(iLam==68){lammin=psigmin;	lammax=psigmax;	sprintf(filename,"%s/Figures/KinDep_pull_HX_sigma_lthstar_rap%d_pt%d.pdf",dirstruct,rapBin,ptBin); 	sprintf(axislabel,"#sigma_{z[#lambda^{*HX}_{#theta}]}");			DrawLine=true;lamline=1;}
				if(iLam==69){lammin=pmumin;		lammax=pmumax;	sprintf(filename,"%s/Figures/KinDep_pull_HX_mean_lphstar_rap%d_pt%d.pdf",dirstruct,rapBin,ptBin); 		sprintf(axislabel,"#mu_{[z[#lambda^{*HX}_{#phi}]}");				DrawLine=true;lamline=0;}
				if(iLam==70){lammin=psigmin;	lammax=psigmax;	sprintf(filename,"%s/Figures/KinDep_pull_HX_sigma_lphstar_rap%d_pt%d.pdf",dirstruct,rapBin,ptBin); 	sprintf(axislabel,"#sigma_{z[#lambda^{*HX}_{#phi}]}");				DrawLine=true;lamline=1;}
				if(iLam==71){lammin=pmumin;		lammax=pmumax;	sprintf(filename,"%s/Figures/KinDep_pull_HX_mean_ltilde_rap%d_pt%d.pdf",dirstruct,rapBin,ptBin); 		sprintf(axislabel,"#mu_{[z[#tilde{#lambda}^{HX}]}");				DrawLine=true;lamline=0;}
				if(iLam==72){lammin=psigmin;	lammax=psigmax;	sprintf(filename,"%s/Figures/KinDep_pull_HX_sigma_ltilde_rap%d_pt%d.pdf",dirstruct,rapBin,ptBin); 		sprintf(axislabel,"#sigma_{z[#tilde{#lambda}^{HX}]}");				DrawLine=true;lamline=1;}

				if(iLam==73){lammin=pmumin;		lammax=pmumax;	sprintf(filename,"%s/Figures/KinDep_pull_PX_mean_lth_rap%d_pt%d.pdf",dirstruct,rapBin,ptBin); 			sprintf(axislabel,"#mu_{[z[#lambda^{PX}_{#theta}]}");				DrawLine=true;lamline=0;}
				if(iLam==74){lammin=psigmin;	lammax=psigmax;	sprintf(filename,"%s/Figures/KinDep_pull_PX_sigma_lth_rap%d_pt%d.pdf",dirstruct,rapBin,ptBin); 		sprintf(axislabel,"#sigma_{z[#lambda^{PX}_{#theta}]}");				DrawLine=true;lamline=1;}
				if(iLam==75){lammin=pmumin;		lammax=pmumax;	sprintf(filename,"%s/Figures/KinDep_pull_PX_mean_lph_rap%d_pt%d.pdf",dirstruct,rapBin,ptBin); 			sprintf(axislabel,"#mu_{[z[#lambda^{PX}_{#phi}]}");					DrawLine=true;lamline=0;}
				if(iLam==76){lammin=psigmin;	lammax=psigmax;	sprintf(filename,"%s/Figures/KinDep_pull_PX_sigma_lph_rap%d_pt%d.pdf",dirstruct,rapBin,ptBin); 		sprintf(axislabel,"#sigma_{z[#lambda^{PX}_{#phi}]}");				DrawLine=true;lamline=1;}
				if(iLam==77){lammin=pmumin;		lammax=pmumax;	sprintf(filename,"%s/Figures/KinDep_pull_PX_mean_ltp_rap%d_pt%d.pdf",dirstruct,rapBin,ptBin); 			sprintf(axislabel,"#mu_{[z[#lambda^{PX}_{#theta#phi}]}");			DrawLine=true;lamline=0;}
				if(iLam==78){lammin=psigmin;	lammax=psigmax;	sprintf(filename,"%s/Figures/KinDep_pull_PX_sigma_ltp_rap%d_pt%d.pdf",dirstruct,rapBin,ptBin); 		sprintf(axislabel,"#sigma_{z[#lambda^{PX}_{#theta#phi}]}");			DrawLine=true;lamline=1;}
				if(iLam==79){lammin=pmumin;		lammax=pmumax;	sprintf(filename,"%s/Figures/KinDep_pull_PX_mean_lthstar_rap%d_pt%d.pdf",dirstruct,rapBin,ptBin); 		sprintf(axislabel,"#mu_{[z[#lambda^{*PX}_{#theta}]}");				DrawLine=true;lamline=0;}
				if(iLam==80){lammin=psigmin;	lammax=psigmax;	sprintf(filename,"%s/Figures/KinDep_pull_PX_sigma_lthstar_rap%d_pt%d.pdf",dirstruct,rapBin,ptBin); 	sprintf(axislabel,"#sigma_{z[#lambda^{*PX}_{#theta}]}");			DrawLine=true;lamline=1;}
				if(iLam==81){lammin=pmumin;		lammax=pmumax;	sprintf(filename,"%s/Figures/KinDep_pull_PX_mean_lphstar_rap%d_pt%d.pdf",dirstruct,rapBin,ptBin); 		sprintf(axislabel,"#mu_{[z[#lambda^{*PX}_{#phi}]}");				DrawLine=true;lamline=0;}
				if(iLam==82){lammin=psigmin;	lammax=psigmax;	sprintf(filename,"%s/Figures/KinDep_pull_PX_sigma_lphstar_rap%d_pt%d.pdf",dirstruct,rapBin,ptBin); 	sprintf(axislabel,"#sigma_{z[#lambda^{*PX}_{#phi}]}");				DrawLine=true;lamline=1;}
				if(iLam==83){lammin=pmumin;		lammax=pmumax;	sprintf(filename,"%s/Figures/KinDep_pull_PX_mean_ltilde_rap%d_pt%d.pdf",dirstruct,rapBin,ptBin); 		sprintf(axislabel,"#mu_{[z[#tilde{#lambda}^{PX}]}");				DrawLine=true;lamline=0;}
				if(iLam==84){lammin=psigmin;	lammax=psigmax;	sprintf(filename,"%s/Figures/KinDep_pull_PX_sigma_ltilde_rap%d_pt%d.pdf",dirstruct,rapBin,ptBin); 		sprintf(axislabel,"#sigma_{z[#tilde{#lambda}^{PX}]}");				DrawLine=true;lamline=1;}


				if(iLam==85) {lammin=dmumin;	lammax=dmumax;	sprintf(filename,"%s/Figures/KinDep_param_CS_deltamean_lth_rap%d_pt%d.pdf",dirstruct,rapBin,ptBin); 			sprintf(axislabel,"#delta_{#lambda^{CS}_{#theta}}");		DrawLine=true;lamline=0;	sprintf(GraphName,"lth_CS_rap%d_pt%d",rapBin,ptBin); if(!RealData) SaveGraph=true;		}
				if(iLam==86) {lammin=thetaphimin;		lammax=thetaphimax;	sprintf(filename,"%s/Figures/KinDep_param_CS_deltamean_lph_rap%d_pt%d.pdf",dirstruct,rapBin,ptBin); 			sprintf(axislabel,"#delta_{#lambda^{CS}_{#phi}}");			DrawLine=true;lamline=0;	sprintf(GraphName,"lph_CS_rap%d_pt%d",rapBin,ptBin); if(!RealData) SaveGraph=true;		}
				if(iLam==87) {lammin=thetaphimin;		lammax=thetaphimax;	sprintf(filename,"%s/Figures/KinDep_param_CS_deltamean_ltp_rap%d_pt%d.pdf",dirstruct,rapBin,ptBin); 			sprintf(axislabel,"#delta_{#lambda^{CS}_{#theta#phi}}");	DrawLine=true;lamline=0;	sprintf(GraphName,"ltp_CS_rap%d_pt%d",rapBin,ptBin); if(!RealData) SaveGraph=true;		}
				if(iLam==88) {lammin=dmumin;	lammax=dmumax;	sprintf(filename,"%s/Figures/KinDep_param_CS_deltamean_lthstar_rap%d_pt%d.pdf",dirstruct,rapBin,ptBin); 		sprintf(axislabel,"#delta_{#lambda^{*CS}_{#theta}}");		DrawLine=true;lamline=0;	sprintf(GraphName,"lthstar_CS_rap%d_pt%d",rapBin,ptBin); if(!RealData) SaveGraph=true;	}
				if(iLam==89) {lammin=dmumin;	lammax=dmumax;	sprintf(filename,"%s/Figures/KinDep_param_CS_deltamean_lphstar_rap%d_pt%d.pdf",dirstruct,rapBin,ptBin); 		sprintf(axislabel,"#delta_{#lambda^{*CS}_{#phi}}");			DrawLine=true;lamline=0;	sprintf(GraphName,"lphstar_CS_rap%d_pt%d",rapBin,ptBin); if(!RealData) SaveGraph=true;	}
				if(iLam==90) {lammin=dmumin;	lammax=dmumax;	sprintf(filename,"%s/Figures/KinDep_param_CS_deltamean_ltilde_rap%d_pt%d.pdf",dirstruct,rapBin,ptBin); 		sprintf(axislabel,"#delta_{#tilde{#lambda}^{CS}}");			DrawLine=true;lamline=0;	sprintf(GraphName,"ltilde_CS_rap%d_pt%d",rapBin,ptBin); if(!RealData) SaveGraph=true;	}

				if(iLam==91) {lammin=dmumin;	lammax=dmumax;	sprintf(filename,"%s/Figures/KinDep_param_HX_deltamean_lth_rap%d_pt%d.pdf",dirstruct,rapBin,ptBin); 			sprintf(axislabel,"#delta_{#lambda^{HX}_{#theta}}");		DrawLine=true;lamline=0;	sprintf(GraphName,"lth_HX_rap%d_pt%d",rapBin,ptBin); if(!RealData) SaveGraph=true;		}
				if(iLam==92) {lammin=thetaphimin;		lammax=thetaphimax;	sprintf(filename,"%s/Figures/KinDep_param_HX_deltamean_lph_rap%d_pt%d.pdf",dirstruct,rapBin,ptBin); 			sprintf(axislabel,"#delta_{#lambda^{HX}_{#phi}}");			DrawLine=true;lamline=0;	sprintf(GraphName,"lph_HX_rap%d_pt%d",rapBin,ptBin); if(!RealData) SaveGraph=true;		}
				if(iLam==93) {lammin=thetaphimin;		lammax=thetaphimax;	sprintf(filename,"%s/Figures/KinDep_param_HX_deltamean_ltp_rap%d_pt%d.pdf",dirstruct,rapBin,ptBin); 			sprintf(axislabel,"#delta_{#lambda^{HX}_{#theta#phi}}");	DrawLine=true;lamline=0;	sprintf(GraphName,"ltp_HX_rap%d_pt%d",rapBin,ptBin); if(!RealData) SaveGraph=true;		}
				if(iLam==94) {lammin=dmumin;	lammax=dmumax;	sprintf(filename,"%s/Figures/KinDep_param_HX_deltamean_lthstar_rap%d_pt%d.pdf",dirstruct,rapBin,ptBin); 		sprintf(axislabel,"#delta_{#lambda^{*HX}_{#theta}}");		DrawLine=true;lamline=0;	sprintf(GraphName,"lthstar_HX_rap%d_pt%d",rapBin,ptBin); if(!RealData) SaveGraph=true;	}
				if(iLam==95) {lammin=dmumin;	lammax=dmumax;	sprintf(filename,"%s/Figures/KinDep_param_HX_deltamean_lphstar_rap%d_pt%d.pdf",dirstruct,rapBin,ptBin); 		sprintf(axislabel,"#delta_{#lambda^{*HX}_{#phi}}");			DrawLine=true;lamline=0;	sprintf(GraphName,"lphstar_HX_rap%d_pt%d",rapBin,ptBin); if(!RealData) SaveGraph=true;	}
				if(iLam==96) {lammin=dmumin;	lammax=dmumax;	sprintf(filename,"%s/Figures/KinDep_param_HX_deltamean_ltilde_rap%d_pt%d.pdf",dirstruct,rapBin,ptBin); 		sprintf(axislabel,"#delta_{#tilde{#lambda}^{HX}}");			DrawLine=true;lamline=0;	sprintf(GraphName,"ltilde_HX_rap%d_pt%d",rapBin,ptBin); if(!RealData) SaveGraph=true;	}

				if(iLam==97) {lammin=dmumin;	lammax=dmumax;	sprintf(filename,"%s/Figures/KinDep_param_PX_deltamean_lth_rap%d_pt%d.pdf",dirstruct,rapBin,ptBin); 			sprintf(axislabel,"#delta_{#lambda^{PX}_{#theta}}");		DrawLine=true;lamline=0;	sprintf(GraphName,"lth_PX_rap%d_pt%d",rapBin,ptBin); if(!RealData) SaveGraph=true;		}
				if(iLam==98) {lammin=thetaphimin;		lammax=thetaphimax;	sprintf(filename,"%s/Figures/KinDep_param_PX_deltamean_lph_rap%d_pt%d.pdf",dirstruct,rapBin,ptBin); 			sprintf(axislabel,"#delta_{#lambda^{PX}_{#phi}}");			DrawLine=true;lamline=0;	sprintf(GraphName,"lph_PX_rap%d_pt%d",rapBin,ptBin); if(!RealData) SaveGraph=true;		}
				if(iLam==99) {lammin=thetaphimin;		lammax=thetaphimax;	sprintf(filename,"%s/Figures/KinDep_param_PX_deltamean_ltp_rap%d_pt%d.pdf",dirstruct,rapBin,ptBin); 			sprintf(axislabel,"#delta_{#lambda^{PX}_{#theta#phi}}");	DrawLine=true;lamline=0;	sprintf(GraphName,"ltp_PX_rap%d_pt%d",rapBin,ptBin); if(!RealData) SaveGraph=true;		}
				if(iLam==100){lammin=dmumin;	lammax=dmumax;	sprintf(filename,"%s/Figures/KinDep_param_PX_deltamean_lthstar_rap%d_pt%d.pdf",dirstruct,rapBin,ptBin); 	sprintf(axislabel,"#delta_{#lambda^{*PX}_{#theta}}");			DrawLine=true;lamline=0;	sprintf(GraphName,"lthstar_PX_rap%d_pt%d",rapBin,ptBin); if(!RealData) SaveGraph=true;	}
				if(iLam==101){lammin=dmumin;	lammax=dmumax;	sprintf(filename,"%s/Figures/KinDep_param_PX_deltamean_lphstar_rap%d_pt%d.pdf",dirstruct,rapBin,ptBin); 		sprintf(axislabel,"#delta_{#lambda^{*PX}_{#phi}}");			DrawLine=true;lamline=0;	sprintf(GraphName,"lphstar_PX_rap%d_pt%d",rapBin,ptBin); if(!RealData) SaveGraph=true;	}
				if(iLam==102){lammin=dmumin;	lammax=dmumax;	sprintf(filename,"%s/Figures/KinDep_param_PX_deltamean_ltilde_rap%d_pt%d.pdf",dirstruct,rapBin,ptBin); 		sprintf(axislabel,"#delta_{#tilde{#lambda}^{PX}}");			DrawLine=true;lamline=0;	sprintf(GraphName,"ltilde_PX_rap%d_pt%d",rapBin,ptBin); if(!RealData) SaveGraph=true;	}

				if(iLam==103){lammin=tildeCompmin;	lammax=tildeCompmax;	if(polScen==6){lammin=lammin+2;lammax=lammax+2;}sprintf(filename,"%s/Figures/KinDep_param_mean_Comparison_ltilde_rap%d_pt%d.pdf",dirstruct,rapBin,ptBin); 		sprintf(axislabel,"#tilde{#lambda}"); lamline=lambda_tilde_injected_CS[rapBin-1][ptBin-1][0];}
				if(iLam==104){lammin=dmumin;	lammax=dmumax;	sprintf(filename,"%s/Figures/KinDep_param_mean_Comparison_lthstar_rap%d_pt%d.pdf",dirstruct,rapBin,ptBin); 	sprintf(axislabel,"#lambda*_{#theta}"); lamline=lambda_thetastar_injected_CS[rapBin-1][ptBin-1][0];}
				if(iLam==105){lammin=dmumin;	lammax=dmumax;	sprintf(filename,"%s/Figures/KinDep_param_mean_Comparison_lphstar_rap%d_pt%d.pdf",dirstruct,rapBin,ptBin); 	sprintf(axislabel,"#lambda*_{#phi}"); lamline=lambda_phistar_injected_CS[rapBin-1][ptBin-1][0];}

				
				int cpm=0;
				for(int cpmBin = cpmBinMin; cpmBin < cpmBinMax+1; cpmBin++) {


					if(iLam==13){ KinDepRap[pt][cpm]= param_lth_CS_mean[rap][pt][cpm]; 		KinDepRaperr_low[pt][cpm]=param_lth_CS_meanerr[rap][pt][cpm];		 	if(RealData) {KinDepRaperr_low[pt][cpm]= param_lth_CS_sigma_low[rap][pt][cpm];		KinDepRaperr_high[pt][cpm]= param_lth_CS_sigma_high[rap][pt][cpm];}    }
					if(iLam==14){ KinDepRap[pt][cpm]= param_lth_CS_sigma[rap][pt][cpm]; 		KinDepRaperr_low[pt][cpm]=param_lth_CS_sigmaerr[rap][pt][cpm];		}
					if(iLam==15){ KinDepRap[pt][cpm]= param_lph_CS_mean[rap][pt][cpm]; 		KinDepRaperr_low[pt][cpm]=param_lph_CS_meanerr[rap][pt][cpm]; 		if(RealData) {KinDepRaperr_low[pt][cpm]= param_lph_CS_sigma_low[rap][pt][cpm];		KinDepRaperr_high[pt][cpm]= param_lph_CS_sigma_high[rap][pt][cpm];}    }
					if(iLam==16){ KinDepRap[pt][cpm]= param_lph_CS_sigma[rap][pt][cpm]; 		KinDepRaperr_low[pt][cpm]=param_lph_CS_sigmaerr[rap][pt][cpm];		}
					if(iLam==17){ KinDepRap[pt][cpm]= param_ltp_CS_mean[rap][pt][cpm]; 		KinDepRaperr_low[pt][cpm]=param_ltp_CS_meanerr[rap][pt][cpm]; 		if(RealData) {KinDepRaperr_low[pt][cpm]= param_ltp_CS_sigma_low[rap][pt][cpm];		KinDepRaperr_high[pt][cpm]= param_ltp_CS_sigma_high[rap][pt][cpm];}    }
					if(iLam==18){ KinDepRap[pt][cpm]= param_ltp_CS_sigma[rap][pt][cpm]; 		KinDepRaperr_low[pt][cpm]=param_ltp_CS_sigmaerr[rap][pt][cpm];		}
					if(iLam==19){ KinDepRap[pt][cpm]= param_lthstar_CS_mean[rap][pt][cpm]; 	KinDepRaperr_low[pt][cpm]=param_lthstar_CS_meanerr[rap][pt][cpm]; 	if(RealData) {KinDepRaperr_low[pt][cpm]= param_lthstar_CS_sigma_low[rap][pt][cpm];	KinDepRaperr_high[pt][cpm]= param_lthstar_CS_sigma_high[rap][pt][cpm];}}
					if(iLam==20){ KinDepRap[pt][cpm]= param_lthstar_CS_sigma[rap][pt][cpm]; 	KinDepRaperr_low[pt][cpm]=param_lthstar_CS_sigmaerr[rap][pt][cpm];	}
					if(iLam==21){ KinDepRap[pt][cpm]= param_lphstar_CS_mean[rap][pt][cpm]; 	KinDepRaperr_low[pt][cpm]=param_lphstar_CS_meanerr[rap][pt][cpm]; 	if(RealData) {KinDepRaperr_low[pt][cpm]= param_lphstar_CS_sigma_low[rap][pt][cpm];	KinDepRaperr_high[pt][cpm]= param_lphstar_CS_sigma_high[rap][pt][cpm];}}
					if(iLam==22){ KinDepRap[pt][cpm]= param_lphstar_CS_sigma[rap][pt][cpm]; 	KinDepRaperr_low[pt][cpm]=param_lphstar_CS_sigmaerr[rap][pt][cpm];	}
					if(iLam==23){ KinDepRap[pt][cpm]= param_ltilde_CS_mean[rap][pt][cpm]; 	KinDepRaperr_low[pt][cpm]=param_ltilde_CS_meanerr[rap][pt][cpm]; 		if(RealData) {KinDepRaperr_low[pt][cpm]= param_ltilde_CS_sigma_low[rap][pt][cpm];		KinDepRaperr_high[pt][cpm]= param_ltilde_CS_sigma_high[rap][pt][cpm];} }
					if(iLam==24){ KinDepRap[pt][cpm]= param_ltilde_CS_sigma[rap][pt][cpm]; 	KinDepRaperr_low[pt][cpm]=param_ltilde_CS_sigmaerr[rap][pt][cpm];		}

					if(iLam==25){ KinDepRap[pt][cpm]= param_lth_HX_mean[rap][pt][cpm]; 		KinDepRaperr_low[pt][cpm]=param_lth_HX_meanerr[rap][pt][cpm]; 		if(RealData) {KinDepRaperr_low[pt][cpm]= param_lth_HX_sigma_low[rap][pt][cpm];		KinDepRaperr_high[pt][cpm]= param_lth_HX_sigma_high[rap][pt][cpm];}    }
					if(iLam==26){ KinDepRap[pt][cpm]= param_lth_HX_sigma[rap][pt][cpm]; 		KinDepRaperr_low[pt][cpm]=param_lth_HX_sigmaerr[rap][pt][cpm];		}
					if(iLam==27){ KinDepRap[pt][cpm]= param_lph_HX_mean[rap][pt][cpm]; 		KinDepRaperr_low[pt][cpm]=param_lph_HX_meanerr[rap][pt][cpm]; 		if(RealData) {KinDepRaperr_low[pt][cpm]= param_lph_HX_sigma_low[rap][pt][cpm];		KinDepRaperr_high[pt][cpm]= param_lph_HX_sigma_high[rap][pt][cpm];}    }
					if(iLam==28){ KinDepRap[pt][cpm]= param_lph_HX_sigma[rap][pt][cpm]; 		KinDepRaperr_low[pt][cpm]=param_lph_HX_sigmaerr[rap][pt][cpm];		}
					if(iLam==29){ KinDepRap[pt][cpm]= param_ltp_HX_mean[rap][pt][cpm]; 		KinDepRaperr_low[pt][cpm]=param_ltp_HX_meanerr[rap][pt][cpm]; 		if(RealData) {KinDepRaperr_low[pt][cpm]= param_ltp_HX_sigma_low[rap][pt][cpm];		KinDepRaperr_high[pt][cpm]= param_ltp_HX_sigma_high[rap][pt][cpm];}    }
					if(iLam==30){ KinDepRap[pt][cpm]= param_ltp_HX_sigma[rap][pt][cpm]; 		KinDepRaperr_low[pt][cpm]=param_ltp_HX_sigmaerr[rap][pt][cpm];		}
					if(iLam==31){ KinDepRap[pt][cpm]= param_lthstar_HX_mean[rap][pt][cpm]; 	KinDepRaperr_low[pt][cpm]=param_lthstar_HX_meanerr[rap][pt][cpm]; 	if(RealData) {KinDepRaperr_low[pt][cpm]= param_lthstar_HX_sigma_low[rap][pt][cpm];	KinDepRaperr_high[pt][cpm]= param_lthstar_HX_sigma_high[rap][pt][cpm];}}
					if(iLam==32){ KinDepRap[pt][cpm]= param_lthstar_HX_sigma[rap][pt][cpm]; 	KinDepRaperr_low[pt][cpm]=param_lthstar_HX_sigmaerr[rap][pt][cpm];	}
					if(iLam==33){ KinDepRap[pt][cpm]= param_lphstar_HX_mean[rap][pt][cpm]; 	KinDepRaperr_low[pt][cpm]=param_lphstar_HX_meanerr[rap][pt][cpm]; 	if(RealData) {KinDepRaperr_low[pt][cpm]= param_lphstar_HX_sigma_low[rap][pt][cpm];	KinDepRaperr_high[pt][cpm]= param_lphstar_HX_sigma_high[rap][pt][cpm];}}
					if(iLam==34){ KinDepRap[pt][cpm]= param_lphstar_HX_sigma[rap][pt][cpm]; 	KinDepRaperr_low[pt][cpm]=param_lphstar_HX_sigmaerr[rap][pt][cpm];	}
					if(iLam==35){ KinDepRap[pt][cpm]= param_ltilde_HX_mean[rap][pt][cpm]; 	KinDepRaperr_low[pt][cpm]=param_ltilde_HX_meanerr[rap][pt][cpm]; 		if(RealData) {KinDepRaperr_low[pt][cpm]= param_ltilde_HX_sigma_low[rap][pt][cpm];		KinDepRaperr_high[pt][cpm]= param_ltilde_HX_sigma_high[rap][pt][cpm];} }
					if(iLam==36){ KinDepRap[pt][cpm]= param_ltilde_HX_sigma[rap][pt][cpm]; 	KinDepRaperr_low[pt][cpm]=param_ltilde_HX_sigmaerr[rap][pt][cpm];		}

					if(iLam==37){ KinDepRap[pt][cpm]= param_lth_PX_mean[rap][pt][cpm]; 		KinDepRaperr_low[pt][cpm]=param_lth_PX_meanerr[rap][pt][cpm]; 		if(RealData) {KinDepRaperr_low[pt][cpm]= param_lth_PX_sigma_low[rap][pt][cpm];		KinDepRaperr_high[pt][cpm]= param_lth_PX_sigma_high[rap][pt][cpm];}    }
					if(iLam==38){ KinDepRap[pt][cpm]= param_lth_PX_sigma[rap][pt][cpm]; 		KinDepRaperr_low[pt][cpm]=param_lth_PX_sigmaerr[rap][pt][cpm];		}
					if(iLam==39){ KinDepRap[pt][cpm]= param_lph_PX_mean[rap][pt][cpm]; 		KinDepRaperr_low[pt][cpm]=param_lph_PX_meanerr[rap][pt][cpm]; 		if(RealData) {KinDepRaperr_low[pt][cpm]= param_lph_PX_sigma_low[rap][pt][cpm];		KinDepRaperr_high[pt][cpm]= param_lph_PX_sigma_high[rap][pt][cpm];}    }
					if(iLam==40){ KinDepRap[pt][cpm]= param_lph_PX_sigma[rap][pt][cpm]; 		KinDepRaperr_low[pt][cpm]=param_lph_PX_sigmaerr[rap][pt][cpm];		}
					if(iLam==41){ KinDepRap[pt][cpm]= param_ltp_PX_mean[rap][pt][cpm]; 		KinDepRaperr_low[pt][cpm]=param_ltp_PX_meanerr[rap][pt][cpm]; 		if(RealData) {KinDepRaperr_low[pt][cpm]= param_ltp_PX_sigma_low[rap][pt][cpm];		KinDepRaperr_high[pt][cpm]= param_ltp_PX_sigma_high[rap][pt][cpm];}    }
					if(iLam==42){ KinDepRap[pt][cpm]= param_ltp_PX_sigma[rap][pt][cpm]; 		KinDepRaperr_low[pt][cpm]=param_ltp_PX_sigmaerr[rap][pt][cpm];		}
					if(iLam==43){ KinDepRap[pt][cpm]= param_lthstar_PX_mean[rap][pt][cpm]; 	KinDepRaperr_low[pt][cpm]=param_lthstar_PX_meanerr[rap][pt][cpm];	 	if(RealData) {KinDepRaperr_low[pt][cpm]= param_lthstar_PX_sigma_low[rap][pt][cpm];	KinDepRaperr_high[pt][cpm]= param_lthstar_PX_sigma_high[rap][pt][cpm];}}
					if(iLam==44){ KinDepRap[pt][cpm]= param_lthstar_PX_sigma[rap][pt][cpm]; 	KinDepRaperr_low[pt][cpm]=param_lthstar_PX_sigmaerr[rap][pt][cpm];	}
					if(iLam==45){ KinDepRap[pt][cpm]= param_lphstar_PX_mean[rap][pt][cpm]; 	KinDepRaperr_low[pt][cpm]=param_lphstar_PX_meanerr[rap][pt][cpm]; 	if(RealData) {KinDepRaperr_low[pt][cpm]= param_lphstar_PX_sigma_low[rap][pt][cpm];	KinDepRaperr_high[pt][cpm]= param_lphstar_PX_sigma_high[rap][pt][cpm];}}
					if(iLam==46){ KinDepRap[pt][cpm]= param_lphstar_PX_sigma[rap][pt][cpm]; 	KinDepRaperr_low[pt][cpm]=param_lphstar_PX_sigmaerr[rap][pt][cpm];	}
					if(iLam==47){ KinDepRap[pt][cpm]= param_ltilde_PX_mean[rap][pt][cpm]; 	KinDepRaperr_low[pt][cpm]=param_ltilde_PX_meanerr[rap][pt][cpm]; 		if(RealData) {KinDepRaperr_low[pt][cpm]= param_ltilde_PX_sigma_low[rap][pt][cpm];		KinDepRaperr_high[pt][cpm]= param_ltilde_PX_sigma_high[rap][pt][cpm];} }
					if(iLam==48){ KinDepRap[pt][cpm]= param_ltilde_PX_sigma[rap][pt][cpm]; 	KinDepRaperr_low[pt][cpm]=param_ltilde_PX_sigmaerr[rap][pt][cpm];		}

					if(iLam==49){ KinDepRap[pt][cpm]= pull_lth_CS_mean[rap][pt][cpm]; 		KinDepRaperr_low[pt][cpm]=pull_lth_CS_meanerr[rap][pt][cpm];			}
					if(iLam==50){ KinDepRap[pt][cpm]= pull_lth_CS_sigma[rap][pt][cpm]; 		KinDepRaperr_low[pt][cpm]=pull_lth_CS_sigmaerr[rap][pt][cpm];			}
					if(iLam==51){ KinDepRap[pt][cpm]= pull_lph_CS_mean[rap][pt][cpm]; 		KinDepRaperr_low[pt][cpm]=pull_lph_CS_meanerr[rap][pt][cpm];			}
					if(iLam==52){ KinDepRap[pt][cpm]= pull_lph_CS_sigma[rap][pt][cpm]; 		KinDepRaperr_low[pt][cpm]=pull_lph_CS_sigmaerr[rap][pt][cpm];			}
					if(iLam==53){ KinDepRap[pt][cpm]= pull_ltp_CS_mean[rap][pt][cpm]; 		KinDepRaperr_low[pt][cpm]=pull_ltp_CS_meanerr[rap][pt][cpm];			}
					if(iLam==54){ KinDepRap[pt][cpm]= pull_ltp_CS_sigma[rap][pt][cpm]; 		KinDepRaperr_low[pt][cpm]=pull_ltp_CS_sigmaerr[rap][pt][cpm];			}
					if(iLam==55){ KinDepRap[pt][cpm]= pull_lthstar_CS_mean[rap][pt][cpm]; 	KinDepRaperr_low[pt][cpm]=pull_lthstar_CS_meanerr[rap][pt][cpm];		}
					if(iLam==56){ KinDepRap[pt][cpm]= pull_lthstar_CS_sigma[rap][pt][cpm]; 	KinDepRaperr_low[pt][cpm]=pull_lthstar_CS_sigmaerr[rap][pt][cpm];		}
					if(iLam==57){ KinDepRap[pt][cpm]= pull_lphstar_CS_mean[rap][pt][cpm]; 	KinDepRaperr_low[pt][cpm]=pull_lphstar_CS_meanerr[rap][pt][cpm];		}
					if(iLam==58){ KinDepRap[pt][cpm]= pull_lphstar_CS_sigma[rap][pt][cpm]; 	KinDepRaperr_low[pt][cpm]=pull_lphstar_CS_sigmaerr[rap][pt][cpm];		}
					if(iLam==59){ KinDepRap[pt][cpm]= pull_ltilde_CS_mean[rap][pt][cpm]; 		KinDepRaperr_low[pt][cpm]=pull_ltilde_CS_meanerr[rap][pt][cpm];		}
					if(iLam==60){ KinDepRap[pt][cpm]= pull_ltilde_CS_sigma[rap][pt][cpm]; 	KinDepRaperr_low[pt][cpm]=pull_ltilde_CS_sigmaerr[rap][pt][cpm];		}

					if(iLam==61){ KinDepRap[pt][cpm]= pull_lth_HX_mean[rap][pt][cpm]; 		KinDepRaperr_low[pt][cpm]=pull_lth_HX_meanerr[rap][pt][cpm];			}
					if(iLam==62){ KinDepRap[pt][cpm]= pull_lth_HX_sigma[rap][pt][cpm]; 		KinDepRaperr_low[pt][cpm]=pull_lth_HX_sigmaerr[rap][pt][cpm];			}
					if(iLam==63){ KinDepRap[pt][cpm]= pull_lph_HX_mean[rap][pt][cpm]; 		KinDepRaperr_low[pt][cpm]=pull_lph_HX_meanerr[rap][pt][cpm];			}
					if(iLam==64){ KinDepRap[pt][cpm]= pull_lph_HX_sigma[rap][pt][cpm]; 		KinDepRaperr_low[pt][cpm]=pull_lph_HX_sigmaerr[rap][pt][cpm];			}
					if(iLam==65){ KinDepRap[pt][cpm]= pull_ltp_HX_mean[rap][pt][cpm]; 		KinDepRaperr_low[pt][cpm]=pull_ltp_HX_meanerr[rap][pt][cpm];			}
					if(iLam==66){ KinDepRap[pt][cpm]= pull_ltp_HX_sigma[rap][pt][cpm]; 		KinDepRaperr_low[pt][cpm]=pull_ltp_HX_sigmaerr[rap][pt][cpm];			}
					if(iLam==67){ KinDepRap[pt][cpm]= pull_lthstar_HX_mean[rap][pt][cpm]; 	KinDepRaperr_low[pt][cpm]=pull_lthstar_HX_meanerr[rap][pt][cpm];		}
					if(iLam==68){ KinDepRap[pt][cpm]= pull_lthstar_HX_sigma[rap][pt][cpm]; 	KinDepRaperr_low[pt][cpm]=pull_lthstar_HX_sigmaerr[rap][pt][cpm];		}
					if(iLam==69){ KinDepRap[pt][cpm]= pull_lphstar_HX_mean[rap][pt][cpm]; 	KinDepRaperr_low[pt][cpm]=pull_lphstar_HX_meanerr[rap][pt][cpm];		}
					if(iLam==70){ KinDepRap[pt][cpm]= pull_lphstar_HX_sigma[rap][pt][cpm]; 	KinDepRaperr_low[pt][cpm]=pull_lphstar_HX_sigmaerr[rap][pt][cpm];		}
					if(iLam==71){ KinDepRap[pt][cpm]= pull_ltilde_HX_mean[rap][pt][cpm]; 		KinDepRaperr_low[pt][cpm]=pull_ltilde_HX_meanerr[rap][pt][cpm];		}
					if(iLam==72){ KinDepRap[pt][cpm]= pull_ltilde_HX_sigma[rap][pt][cpm]; 	KinDepRaperr_low[pt][cpm]=pull_ltilde_HX_sigmaerr[rap][pt][cpm];		}

					if(iLam==73){ KinDepRap[pt][cpm]= pull_lth_PX_mean[rap][pt][cpm]; 		KinDepRaperr_low[pt][cpm]=pull_lth_PX_meanerr[rap][pt][cpm];			}
					if(iLam==74){ KinDepRap[pt][cpm]= pull_lth_PX_sigma[rap][pt][cpm]; 		KinDepRaperr_low[pt][cpm]=pull_lth_PX_sigmaerr[rap][pt][cpm];			}
					if(iLam==75){ KinDepRap[pt][cpm]= pull_lph_PX_mean[rap][pt][cpm]; 		KinDepRaperr_low[pt][cpm]=pull_lph_PX_meanerr[rap][pt][cpm];			}
					if(iLam==76){ KinDepRap[pt][cpm]= pull_lph_PX_sigma[rap][pt][cpm]; 		KinDepRaperr_low[pt][cpm]=pull_lph_PX_sigmaerr[rap][pt][cpm];			}
					if(iLam==77){ KinDepRap[pt][cpm]= pull_ltp_PX_mean[rap][pt][cpm]; 		KinDepRaperr_low[pt][cpm]=pull_ltp_PX_meanerr[rap][pt][cpm];			}
					if(iLam==78){ KinDepRap[pt][cpm]= pull_ltp_PX_sigma[rap][pt][cpm]; 		KinDepRaperr_low[pt][cpm]=pull_ltp_PX_sigmaerr[rap][pt][cpm];			}
					if(iLam==79){ KinDepRap[pt][cpm]= pull_lthstar_PX_mean[rap][pt][cpm]; 	KinDepRaperr_low[pt][cpm]=pull_lthstar_PX_meanerr[rap][pt][cpm];		}
					if(iLam==80){ KinDepRap[pt][cpm]= pull_lthstar_PX_sigma[rap][pt][cpm]; 	KinDepRaperr_low[pt][cpm]=pull_lthstar_PX_sigmaerr[rap][pt][cpm];		}
					if(iLam==81){ KinDepRap[pt][cpm]= pull_lphstar_PX_mean[rap][pt][cpm]; 	KinDepRaperr_low[pt][cpm]=pull_lphstar_PX_meanerr[rap][pt][cpm];		}
					if(iLam==82){ KinDepRap[pt][cpm]= pull_lphstar_PX_sigma[rap][pt][cpm]; 	KinDepRaperr_low[pt][cpm]=pull_lphstar_PX_sigmaerr[rap][pt][cpm];		}
					if(iLam==83){ KinDepRap[pt][cpm]= pull_ltilde_PX_mean[rap][pt][cpm]; cout<<"blah blah kin dep is "<<KinDepRap[pt][cpm]<<endl;		KinDepRaperr_low[pt][cpm]=pull_ltilde_PX_meanerr[rap][pt][cpm];		cout<<"the mean is "<<pull_ltilde_PX_mean[rap][pt][cpm]<<endl;}
					if(iLam==84){ KinDepRap[pt][cpm]= pull_ltilde_PX_sigma[rap][pt][cpm]; 	KinDepRaperr_low[pt][cpm]=pull_ltilde_PX_sigmaerr[rap][pt][cpm];		 }


					if(iLam==85){ KinDepRap[pt][cpm]= param_lth_CS_mean[rap][pt][cpm]-lambda_theta_injected_CS[rap][pt][cpm]; 			KinDepRaperr_low[pt][cpm]=param_lth_CS_sigma[rap][pt][cpm];}
					if(iLam==86){ KinDepRap[pt][cpm]= param_lph_CS_mean[rap][pt][cpm]-lambda_phi_injected_CS[rap][pt][cpm]; 			KinDepRaperr_low[pt][cpm]=param_lph_CS_sigma[rap][pt][cpm];}
					if(iLam==87){ KinDepRap[pt][cpm]= param_ltp_CS_mean[rap][pt][cpm]-lambda_thetaphi_injected_CS[rap][pt][cpm]; 		KinDepRaperr_low[pt][cpm]=param_ltp_CS_sigma[rap][pt][cpm];}
					if(iLam==88){ KinDepRap[pt][cpm]= param_lthstar_CS_mean[rap][pt][cpm]-lambda_thetastar_injected_CS[rap][pt][cpm]; 	KinDepRaperr_low[pt][cpm]=param_lthstar_CS_sigma[rap][pt][cpm];}
					if(iLam==89){ KinDepRap[pt][cpm]= param_lphstar_CS_mean[rap][pt][cpm]-lambda_phistar_injected_CS[rap][pt][cpm]; 	KinDepRaperr_low[pt][cpm]=param_lphstar_CS_sigma[rap][pt][cpm];}
					if(iLam==90){ KinDepRap[pt][cpm]= param_ltilde_CS_mean[rap][pt][cpm]-lambda_tilde_injected_CS[rap][pt][cpm]; 		KinDepRaperr_low[pt][cpm]=param_ltilde_CS_sigma[rap][pt][cpm];}

					if(iLam==91){ KinDepRap[pt][cpm]= param_lth_HX_mean[rap][pt][cpm]-lambda_theta_injected_HX[rap][pt][cpm]; 			KinDepRaperr_low[pt][cpm]=param_lth_HX_sigma[rap][pt][cpm];}
					if(iLam==92){ KinDepRap[pt][cpm]= param_lph_HX_mean[rap][pt][cpm]-lambda_phi_injected_HX[rap][pt][cpm]; 			KinDepRaperr_low[pt][cpm]=param_lph_HX_sigma[rap][pt][cpm];}
					if(iLam==93){ KinDepRap[pt][cpm]= param_ltp_HX_mean[rap][pt][cpm]-lambda_thetaphi_injected_HX[rap][pt][cpm]; 		KinDepRaperr_low[pt][cpm]=param_ltp_HX_sigma[rap][pt][cpm];}
					if(iLam==94){ KinDepRap[pt][cpm]= param_lthstar_HX_mean[rap][pt][cpm]-lambda_thetastar_injected_HX[rap][pt][cpm]; 	KinDepRaperr_low[pt][cpm]=param_lthstar_HX_sigma[rap][pt][cpm];}
					if(iLam==95){ KinDepRap[pt][cpm]= param_lphstar_HX_mean[rap][pt][cpm]-lambda_phistar_injected_HX[rap][pt][cpm]; 	KinDepRaperr_low[pt][cpm]=param_lphstar_HX_sigma[rap][pt][cpm];}
					if(iLam==96){ KinDepRap[pt][cpm]= param_ltilde_HX_mean[rap][pt][cpm]-lambda_tilde_injected_HX[rap][pt][cpm]; 		KinDepRaperr_low[pt][cpm]=param_ltilde_HX_sigma[rap][pt][cpm];}

					if(iLam==97){ KinDepRap[pt][cpm]= param_lth_PX_mean[rap][pt][cpm]-lambda_theta_injected_PX[rap][pt][cpm]; 			KinDepRaperr_low[pt][cpm]=param_lth_PX_sigma[rap][pt][cpm];}
					if(iLam==98){ KinDepRap[pt][cpm]= param_lph_PX_mean[rap][pt][cpm]-lambda_phi_injected_PX[rap][pt][cpm]; 			KinDepRaperr_low[pt][cpm]=param_lph_PX_sigma[rap][pt][cpm];}
					if(iLam==99){ KinDepRap[pt][cpm]= param_ltp_PX_mean[rap][pt][cpm]-lambda_thetaphi_injected_PX[rap][pt][cpm]; 		KinDepRaperr_low[pt][cpm]=param_ltp_PX_sigma[rap][pt][cpm];}
					if(iLam==100){ KinDepRap[pt][cpm]= param_lthstar_PX_mean[rap][pt][cpm]-lambda_thetastar_injected_PX[rap][pt][cpm]; KinDepRaperr_low[pt][cpm]=param_lthstar_PX_sigma[rap][pt][cpm];}
					if(iLam==101){ KinDepRap[pt][cpm]= param_lphstar_PX_mean[rap][pt][cpm]-lambda_phistar_injected_PX[rap][pt][cpm]; 	KinDepRaperr_low[pt][cpm]=param_lphstar_PX_sigma[rap][pt][cpm];}
					if(iLam==102){ KinDepRap[pt][cpm]= param_ltilde_PX_mean[rap][pt][cpm]-lambda_tilde_injected_PX[rap][pt][cpm]; 		KinDepRaperr_low[pt][cpm]=param_ltilde_PX_sigma[rap][pt][cpm];}

					if(iLam==103){  KinDepRap[pt][cpm]= param_ltilde_CS_mean[rap][pt][cpm]; 		 KinDepRaperr_low[pt][cpm]=param_ltilde_CS_sigma_low[rap][pt][cpm];	 KinDepRaperr_high[pt][cpm]=param_ltilde_CS_sigma_high[rap][pt][cpm];	 KinDepRap2[pt][cpm]= param_ltilde_HX_mean[rap][pt][cpm]; 		 KinDepRap2err_low[pt][cpm]=param_ltilde_HX_sigma_low[rap][pt][cpm];	 KinDepRap2err_high[pt][cpm]=param_ltilde_HX_sigma_high[rap][pt][cpm];	KinDepRap3[pt][cpm]= param_ltilde_PX_mean[rap][pt][cpm]; 		 	 KinDepRap3err_low[pt][cpm]=param_ltilde_PX_sigma_low[rap][pt][cpm];	KinDepRap3err_high[pt][cpm]=param_ltilde_PX_sigma_high[rap][pt][cpm];}
					if(iLam==104){ KinDepRap[pt][cpm]= param_lthstar_CS_mean[rap][pt][cpm]; 		KinDepRaperr_low[pt][cpm]=param_lthstar_CS_sigma_low[rap][pt][cpm];	KinDepRaperr_high[pt][cpm]=param_lthstar_CS_sigma_high[rap][pt][cpm];	KinDepRap2[pt][cpm]= param_lthstar_HX_mean[rap][pt][cpm]; 		KinDepRap2err_low[pt][cpm]=param_lthstar_HX_sigma_low[rap][pt][cpm];	KinDepRap2err_high[pt][cpm]=param_lthstar_HX_sigma_high[rap][pt][cpm];	KinDepRap3[pt][cpm]= param_lthstar_PX_mean[rap][pt][cpm]; 		KinDepRap3err_low[pt][cpm]=param_lthstar_PX_sigma_low[rap][pt][cpm];	KinDepRap3err_high[pt][cpm]=param_lthstar_PX_sigma_high[rap][pt][cpm];}
					if(iLam==105){ KinDepRap[pt][cpm]= param_lphstar_CS_mean[rap][pt][cpm]; 		KinDepRaperr_low[pt][cpm]=param_lphstar_CS_sigma_low[rap][pt][cpm];	KinDepRaperr_high[pt][cpm]=param_lphstar_CS_sigma_high[rap][pt][cpm];	KinDepRap2[pt][cpm]= param_lphstar_HX_mean[rap][pt][cpm]; 		KinDepRap2err_low[pt][cpm]=param_lphstar_HX_sigma_low[rap][pt][cpm];	KinDepRap2err_high[pt][cpm]=param_lphstar_HX_sigma_high[rap][pt][cpm];	KinDepRap3[pt][cpm]= param_lphstar_PX_mean[rap][pt][cpm]; 		KinDepRap3err_low[pt][cpm]=param_lphstar_PX_sigma_low[rap][pt][cpm];	KinDepRap3err_high[pt][cpm]=param_lphstar_PX_sigma_high[rap][pt][cpm];}

//					if(emptyBin[rap][pt]){ KinDepRap[pt][cpm]=999; KinDepRap2[pt][cpm]=999; KinDepRap3[pt][cpm]=999;}

					if(iLam<103 && !RealData) {KinDepRaperr_high[pt][cpm]=KinDepRaperr_low[pt][cpm];}
					
					
				//	cout<<"blah blah mean is "<<pull_ltilde_PX_mean[rap][pt][cpm]<<endl;
					
					cpm++;
					
					}
					

				if(iLam>12&&iLam<25||iLam>48&&iLam<61||iLam>84&&iLam<91) nFrame=1;
				if(iLam>24&&iLam<37||iLam>60&&iLam<73||iLam>90&&iLam<97) nFrame=2;
				if(iLam>36&&iLam<49||iLam>72&&iLam<85||iLam>96&&iLam<103) nFrame=3;
				if(iLam<49 && iLam % 2 == 0 && RealData) continue;
				if(iLam>48 && iLam < 103 && RealData) continue;
				if(iLam>102 && RealData) lamline=1000;
				if(iLam<103) PlotRapPt(ptBin, cpmBinMin, cpmBinMax, rapBin, lammin, lammax, axislabel, filename, KinDepRap[pt], KinDepRaperr_low[pt], KinDepRaperr_high[pt], DrawLine, lamline, GraphFolder, GraphName, SaveGraph, nState, nFrame, RealData, realdatadir, TreeID);
//				void PlotRapPt(int ptBinMin, int ptBinMax, int cpmBinMin, int cpmBinMax, int rapBin, double yMin, double yMax, char yAxisTitle[200], char filename[200], double lmean[ToyMC::nPtBins][ToyMC::ncpmBins], double lmeanerr_low[ToyMC::nPtBins][ToyMC::ncpmBins], double lmeanerr_high[ToyMC::nPtBins][ToyMC::ncpmBins], bool DrawLine, double lamline, char GraphFolder[200], char GraphName[200], bool SaveGraph, int nState, int nFrame, bool RealData, char* realdatadir, char* TreeID){
//error: cannot convert 'double*' to 'double (*)[2]' for argument '10' to 'void PlotRapPt(int, int, int, int, int, double, double, char*, char*, double (*)[2], double (*)[2], double (*)[2], bool, double, char*, char*, bool, int, int, bool, char*, char*)'
				
				
				else PlotComparisonRapPt(ptBin, cpmBinMin, cpmBinMax, rapBin, lammin, lammax, axislabel, filename, KinDepRap[pt], KinDepRaperr_low[pt], KinDepRaperr_high[pt], KinDepRap2[pt], KinDepRap2err_low[pt], KinDepRap2err_high[pt], KinDepRap3[pt], KinDepRap3err_low[pt], KinDepRap3err_high[pt], lamline, nState, RealData, realdatadir, TreeID);
				
				pt++;
				}
				rap++;
			}



		}



		///////////////// TABLE PRODUCTION /////////////////////////////////////////////////////////////////////////


		char framerap[200];
		int nTables=3;

		char NumFileName[200];
		sprintf(NumFileName,"ToyNumericalResults_%s.tex",dirstruct);
		if(RealData) sprintf(NumFileName,"ToyNumericalResults.tex");
		FILE *NumFile = fopen(NumFileName,"w");

		fprintf(NumFile, "\n");
		fprintf(NumFile,"\\documentclass{article}\n\\usepackage[applemac]{inputenc}\n\\usepackage{amsmath}\n\\usepackage{textcomp}\n\\pagestyle{plain}\n\\usepackage{graphicx}\n\\usepackage{multicol}\n\\usepackage{geometry}\n\\geometry{a3paper,left=20mm,right=20mm, top=1.5cm, bottom=1.5cm}\n\\usepackage{subfigure}\n\\usepackage{booktabs}\n\\usepackage{setspace}\n\n\n\n\\begin{document}\n");

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
		double lthstar_tab;
		double lthstarerr_tab;
		double lthstarerr_high_tab;
		double lphstar_tab;
		double lphstarerr_tab;
		double lphstarerr_high_tab;

		if(!RealData){


			double meandeviation[3][3][6]={NULL};//index 1...table,2...frame,3...parameter
			int n_=0;
			rap=0;
			for(int rapBin = rapBinMin; rapBin < rapBinMax+1; rapBin++) {
				int pt=0;
				for(int ptBin = ptBinMin; ptBin < ptBinMax+1; ptBin++) {
				int cpm=0;
				  for(int cpmBin = cpmBinMin; cpmBin < cpmBinMax+1; cpmBin++) {
/////////did not edit the if statement!!!!
	/*				if (!emptyBin[rap][pt][cpm]) {
						meandeviation[0][0][0]=meandeviation[0][0][0]+TMath::Abs((param_lth_CS_mean[rap][pt]-lambda_theta_injected_CS[rap][pt]))/param_lth_CS_sigma[rap][pt];
						meandeviation[0][0][1]=meandeviation[0][0][1]+TMath::Abs((param_lph_CS_mean[rap][pt]-lambda_phi_injected_CS[rap][pt]))/param_lph_CS_sigma[rap][pt];
						meandeviation[0][0][2]=meandeviation[0][0][2]+TMath::Abs((param_ltp_CS_mean[rap][pt]-lambda_thetaphi_injected_CS[rap][pt]))/param_ltp_CS_sigma[rap][pt];
						meandeviation[0][0][3]=meandeviation[0][0][3]+TMath::Abs((param_ltilde_CS_mean[rap][pt]-lambda_tilde_injected_CS[rap][pt]))/param_ltilde_CS_sigma[rap][pt];
						meandeviation[0][0][4]=meandeviation[0][0][4]+TMath::Abs((param_lthstar_CS_mean[rap][pt]-lambda_thetastar_injected_CS[rap][pt]))/param_lthstar_CS_sigma[rap][pt];
						meandeviation[0][0][5]=meandeviation[0][0][5]+TMath::Abs((param_lphstar_CS_mean[rap][pt]-lambda_phistar_injected_CS[rap][pt]))/param_lphstar_CS_sigma[rap][pt];

						meandeviation[0][1][0]=meandeviation[0][1][0]+TMath::Abs((param_lth_HX_mean[rap][pt]-lambda_theta_injected_HX[rap][pt]))/param_lth_HX_sigma[rap][pt];
						meandeviation[0][1][1]=meandeviation[0][1][1]+TMath::Abs((param_lph_HX_mean[rap][pt]-lambda_phi_injected_HX[rap][pt]))/param_lph_HX_sigma[rap][pt];
						meandeviation[0][1][2]=meandeviation[0][1][2]+TMath::Abs((param_ltp_HX_mean[rap][pt]-lambda_thetaphi_injected_HX[rap][pt]))/param_ltp_HX_sigma[rap][pt];
						meandeviation[0][1][3]=meandeviation[0][1][3]+TMath::Abs((param_ltilde_HX_mean[rap][pt]-lambda_tilde_injected_HX[rap][pt]))/param_ltilde_HX_sigma[rap][pt];
						meandeviation[0][1][4]=meandeviation[0][1][4]+TMath::Abs((param_lthstar_HX_mean[rap][pt]-lambda_thetastar_injected_HX[rap][pt]))/param_lthstar_HX_sigma[rap][pt];
						meandeviation[0][1][5]=meandeviation[0][1][5]+TMath::Abs((param_lphstar_HX_mean[rap][pt]-lambda_phistar_injected_HX[rap][pt]))/param_lphstar_HX_sigma[rap][pt];

						meandeviation[0][2][0]=meandeviation[0][2][0]+TMath::Abs((param_lth_PX_mean[rap][pt]-lambda_theta_injected_PX[rap][pt]))/param_lth_PX_sigma[rap][pt];
						meandeviation[0][2][1]=meandeviation[0][2][1]+TMath::Abs((param_lph_PX_mean[rap][pt]-lambda_phi_injected_PX[rap][pt]))/param_lph_PX_sigma[rap][pt];
						meandeviation[0][2][2]=meandeviation[0][2][2]+TMath::Abs((param_ltp_PX_mean[rap][pt]-lambda_thetaphi_injected_PX[rap][pt]))/param_ltp_PX_sigma[rap][pt];
						meandeviation[0][2][3]=meandeviation[0][2][3]+TMath::Abs((param_ltilde_PX_mean[rap][pt]-lambda_tilde_injected_PX[rap][pt]))/param_ltilde_PX_sigma[rap][pt];
						meandeviation[0][2][4]=meandeviation[0][2][4]+TMath::Abs((param_lthstar_PX_mean[rap][pt]-lambda_thetastar_injected_PX[rap][pt]))/param_lthstar_PX_sigma[rap][pt];
						meandeviation[0][2][5]=meandeviation[0][2][5]+TMath::Abs((param_lphstar_PX_mean[rap][pt]-lambda_phistar_injected_PX[rap][pt]))/param_lphstar_PX_sigma[rap][pt];


						meandeviation[1][0][0]=meandeviation[1][0][0]+TMath::Abs(pull_lth_CS_mean[rap][pt])/pull_lth_CS_meanerr[rap][pt];
						meandeviation[1][0][1]=meandeviation[1][0][1]+TMath::Abs(pull_lph_CS_mean[rap][pt])/pull_lph_CS_meanerr[rap][pt];
						meandeviation[1][0][2]=meandeviation[1][0][2]+TMath::Abs(pull_ltp_CS_mean[rap][pt])/pull_ltp_CS_meanerr[rap][pt];
						meandeviation[1][0][3]=meandeviation[1][0][3]+TMath::Abs(pull_ltilde_CS_mean[rap][pt])/pull_ltilde_CS_meanerr[rap][pt];
						meandeviation[1][0][4]=meandeviation[1][0][4]+TMath::Abs(pull_lthstar_CS_mean[rap][pt])/pull_lthstar_CS_meanerr[rap][pt];
						meandeviation[1][0][5]=meandeviation[1][0][5]+TMath::Abs(pull_lphstar_CS_mean[rap][pt])/pull_lphstar_CS_meanerr[rap][pt];

						meandeviation[1][1][0]=meandeviation[1][1][0]+TMath::Abs(pull_lth_HX_mean[rap][pt])/pull_lth_HX_meanerr[rap][pt];
						meandeviation[1][1][1]=meandeviation[1][1][1]+TMath::Abs(pull_lph_HX_mean[rap][pt])/pull_lph_HX_meanerr[rap][pt];
						meandeviation[1][1][2]=meandeviation[1][1][2]+TMath::Abs(pull_ltp_HX_mean[rap][pt])/pull_ltp_HX_meanerr[rap][pt];
						meandeviation[1][1][3]=meandeviation[1][1][3]+TMath::Abs(pull_ltilde_HX_mean[rap][pt])/pull_ltilde_HX_meanerr[rap][pt];
						meandeviation[1][1][4]=meandeviation[1][1][4]+TMath::Abs(pull_lthstar_HX_mean[rap][pt])/pull_lthstar_HX_meanerr[rap][pt];
						meandeviation[1][1][5]=meandeviation[1][1][5]+TMath::Abs(pull_lphstar_HX_mean[rap][pt])/pull_lphstar_HX_meanerr[rap][pt];

						meandeviation[1][2][0]=meandeviation[1][2][0]+TMath::Abs(pull_lth_PX_mean[rap][pt])/pull_lth_PX_meanerr[rap][pt];
						meandeviation[1][2][1]=meandeviation[1][2][1]+TMath::Abs(pull_lph_PX_mean[rap][pt])/pull_lph_PX_meanerr[rap][pt];
						meandeviation[1][2][2]=meandeviation[1][2][2]+TMath::Abs(pull_ltp_PX_mean[rap][pt])/pull_ltp_PX_meanerr[rap][pt];
						meandeviation[1][2][3]=meandeviation[1][2][3]+TMath::Abs(pull_ltilde_PX_mean[rap][pt])/pull_ltilde_PX_meanerr[rap][pt];
						meandeviation[1][2][4]=meandeviation[1][2][4]+TMath::Abs(pull_lthstar_PX_mean[rap][pt])/pull_lthstar_PX_meanerr[rap][pt];
						meandeviation[1][2][5]=meandeviation[1][2][5]+TMath::Abs(pull_lphstar_PX_mean[rap][pt])/pull_lphstar_PX_meanerr[rap][pt];


						meandeviation[2][0][0]=meandeviation[2][0][0]+TMath::Abs(pull_lth_CS_sigma[rap][pt]-1)/pull_lth_CS_sigmaerr[rap][pt];
						meandeviation[2][0][1]=meandeviation[2][0][1]+TMath::Abs(pull_lph_CS_sigma[rap][pt]-1)/pull_lph_CS_sigmaerr[rap][pt];
						meandeviation[2][0][2]=meandeviation[2][0][2]+TMath::Abs(pull_ltp_CS_sigma[rap][pt]-1)/pull_ltp_CS_sigmaerr[rap][pt];
						meandeviation[2][0][3]=meandeviation[2][0][3]+TMath::Abs(pull_ltilde_CS_sigma[rap][pt]-1)/pull_ltilde_CS_sigmaerr[rap][pt];
						meandeviation[2][0][4]=meandeviation[2][0][4]+TMath::Abs(pull_lthstar_CS_sigma[rap][pt]-1)/pull_lthstar_CS_sigmaerr[rap][pt];
						meandeviation[2][0][5]=meandeviation[2][0][5]+TMath::Abs(pull_lphstar_CS_sigma[rap][pt]-1)/pull_lphstar_CS_sigmaerr[rap][pt];

						meandeviation[2][1][0]=meandeviation[2][1][0]+TMath::Abs(pull_lth_HX_sigma[rap][pt]-1)/pull_lth_HX_sigmaerr[rap][pt];
						meandeviation[2][1][1]=meandeviation[2][1][1]+TMath::Abs(pull_lph_HX_sigma[rap][pt]-1)/pull_lph_HX_sigmaerr[rap][pt];
						meandeviation[2][1][2]=meandeviation[2][1][2]+TMath::Abs(pull_ltp_HX_sigma[rap][pt]-1)/pull_ltp_HX_sigmaerr[rap][pt];
						meandeviation[2][1][3]=meandeviation[2][1][3]+TMath::Abs(pull_ltilde_HX_sigma[rap][pt]-1)/pull_ltilde_HX_sigmaerr[rap][pt];
						meandeviation[2][1][4]=meandeviation[2][1][4]+TMath::Abs(pull_lthstar_HX_sigma[rap][pt]-1)/pull_lthstar_HX_sigmaerr[rap][pt];
						meandeviation[2][1][5]=meandeviation[2][1][5]+TMath::Abs(pull_lphstar_HX_sigma[rap][pt]-1)/pull_lphstar_HX_sigmaerr[rap][pt];

						meandeviation[2][2][0]=meandeviation[2][2][0]+TMath::Abs(pull_lth_PX_sigma[rap][pt]-1)/pull_lth_PX_sigmaerr[rap][pt];
						meandeviation[2][2][1]=meandeviation[2][2][1]+TMath::Abs(pull_lph_PX_sigma[rap][pt]-1)/pull_lph_PX_sigmaerr[rap][pt];
						meandeviation[2][2][2]=meandeviation[2][2][2]+TMath::Abs(pull_ltp_PX_sigma[rap][pt]-1)/pull_ltp_PX_sigmaerr[rap][pt];
						meandeviation[2][2][3]=meandeviation[2][2][3]+TMath::Abs(pull_ltilde_PX_sigma[rap][pt]-1)/pull_ltilde_PX_sigmaerr[rap][pt];
						meandeviation[2][2][4]=meandeviation[2][2][4]+TMath::Abs(pull_lthstar_PX_sigma[rap][pt]-1)/pull_lthstar_PX_sigmaerr[rap][pt];
						meandeviation[2][2][5]=meandeviation[2][2][5]+TMath::Abs(pull_lphstar_PX_sigma[rap][pt]-1)/pull_lphstar_PX_sigmaerr[rap][pt];

						//							cout<<meandeviation[0][0][5]<<endl;

						n_++;
					}
					else cout<<"empty bin"<<endl;
	*/
					cpm++;
					}
					pt++;
				}
				rap++;
			}
			cout<<n_<<endl;
			for(int i=0;i<3;i++){
				for(int j=0;j<3;j++){
					for(int k=0;k<6;k++){
						meandeviation[i][j][k]=meandeviation[i][j][k]/n_;
					}
				}
			}

			for(int iTab=1; iTab<nTables+1;iTab++){

				if(iTab==1){
					fprintf(NumFile, "\\begin{table}[!h]\n\\centering \\caption{Mean Deviation of the mean of the distribution of the parameter estimates $\\lambda_{i}$ from the injected parameter: $n\\sigma(\\delta_{\\lambda_{i}})=\\sum_{j=1}^{N}\\frac{1}{N}\\frac{|\\mu^j_{\\lambda_{i}}-\\lambda^{Truth}_{i}|}{\\sigma^j_{\\lambda_{i}}}$. Averaged over all N kinematic bins j, in units of $\\sigma_{\\lambda_{i}}$.}\n\\begin{tabular}{|cccccc|}\n\\hline\n");
					fprintf(NumFile, "$n\\sigma(\\delta_{\\lambda_{\\vartheta}})$ & $n\\sigma(\\delta_{\\lambda_{\\varphi}})$ &  $n\\sigma(\\delta_{\\lambda_{\\vartheta \\varphi}})$ & $n\\sigma(\\delta_{\\tilde{\\lambda}})$ & $n\\sigma(\\delta_{\\lambda^*_{\\vartheta}})$ & $n\\sigma(\\delta_{\\lambda^*_{\\varphi}})$ \\\\\n");
				}

				if(iTab==2){
					fprintf(NumFile, "\\begin{table}[!h]\n\\centering \\caption{Mean Deviation of the mean of the distribution of the standard score $z(\\lambda_{i})$ from 0: $n\\sigma(\\mu_{z(\\lambda_{i})})=\\sum_{j=1}^{N}\\frac{1}{N}\\frac{|\\mu_{z(\\lambda_{i})}|}{\\sigma^j_{\\mu_{z(\\lambda_{i})}}}$. Averaged over all N kinematic bins j, in units of $\\sigma_{\\mu_{z(\\lambda_{i})}}$.}\n\\begin{tabular}{|cccccc|}\n\\hline\n");
					fprintf(NumFile, "$n\\sigma(\\mu_{z(\\lambda_{\\vartheta})})$ & $n\\sigma(\\mu_{z(\\lambda_{\\varphi})})$ &  $n\\sigma(\\mu_{z(\\lambda_{\\vartheta \\varphi})})$ & $n\\sigma(\\mu_{z(\\tilde{\\lambda})})$ & $n\\sigma(\\mu_{z(\\lambda^*_{\\vartheta})})$ & $n\\sigma(\\mu_{z(\\lambda^*_{\\varphi})})$ \\\\\n");
				}

				if(iTab==3){
					fprintf(NumFile, "\\begin{table}[!h]\n\\centering \\caption{Mean Deviation of the r.m.s. of the distribution of the standard score $z(\\lambda_{i})$ from 1: $n\\sigma(\\sigma_{z(\\lambda_{i})})=\\sum_{j=1}^{N}\\frac{1}{N}\\frac{|\\sigma_{z(\\lambda_{i})}-1|}{\\sigma^j_{\\sigma_{z(\\lambda_{i})}}}$. Averaged over all N kinematic bins j, in units of $\\sigma_{\\sigma_{z(\\lambda_{i})}}$.}\n\\begin{tabular}{|cccccc|}\n\\hline\n");
					fprintf(NumFile, "$n\\sigma(\\sigma_{z(\\lambda_{\\vartheta})})$ & $n\\sigma(\\sigma_{z(\\lambda_{\\varphi})})$ &  $n\\sigma(\\sigma_{z(\\lambda_{\\vartheta \\varphi})})$ & $n\\sigma(\\sigma_{z(\\tilde{\\lambda})})$ & $n\\sigma(\\sigma_{z(\\lambda^*_{\\vartheta})})$ & $n\\sigma(\\sigma_{z(\\lambda^*_{\\varphi})})$ \\\\\n");
				}

				for(int iFrame=1; iFrame<4; iFrame++){


					if(iFrame==1) sprintf(framerap,"\\hline \\multicolumn{6}{|c|}{CS frame}\\\\ \\hline \\rule{0pt}{4mm}\n");
					if(iFrame==2) sprintf(framerap,"\\hline \\multicolumn{6}{|c|}{HX frame}\\\\ \\hline \\rule{0pt}{4mm}\n");
					if(iFrame==3) sprintf(framerap,"\\hline \\multicolumn{6}{|c|}{PX frame}\\\\ \\hline \\rule{0pt}{4mm}\n");

					fprintf(NumFile,framerap);
					fprintf(NumFile, "$%1.1f$ & $%1.1f$ & $%1.1f$ & $%1.1f$ & $%1.1f$ & $%1.1f$ \\\\\n",meandeviation[iTab-1][iFrame-1][0], meandeviation[iTab-1][iFrame-1][1], meandeviation[iTab-1][iFrame-1][2], meandeviation[iTab-1][iFrame-1][3], meandeviation[iTab-1][iFrame-1][4], meandeviation[iTab-1][iFrame-1][5]);



				}

				fprintf(NumFile, "\\hline\n");
				fprintf(NumFile, "\\end{tabular}\n");
				fprintf(NumFile, "\\label{tab:syst_acceptance}\n");
				fprintf(NumFile, "\\end{table}\n");
				fprintf(NumFile, "\n");

			}




			for(int iTab=1; iTab<nTables+1;iTab++){

				fprintf(NumFile, "\n\n\n\n");

				if(iTab==1){
					fprintf(NumFile, "\\begin{table}[!h]\n\\centering \\caption{Mean Deviation of the distribution of the parameter estimates $\\lambda_{i}$: $\\delta_{\\lambda_{i}}=(\\mu_{\\lambda_{i}}-\\lambda^{Truth}_{i}) \\pm  \\sigma_{\\lambda_{i}}$}\n\\begin{tabular}{|c|cccccc|}\n\\hline\n");
					fprintf(NumFile, "$p_{T}$ [GeV] & $\\delta_{\\lambda_{\\vartheta}}$ & $\\delta_{\\lambda_{\\varphi}}$ &  $\\delta_{\\lambda_{\\vartheta \\varphi}}$ & $\\delta_{\\tilde{\\lambda}}$ & $\\delta_{\\lambda^*_{\\vartheta}}$ & $\\delta_{\\lambda^*_{\\varphi}}$ \\\\\n");
				}
				if(iTab==2){
					fprintf(NumFile, "\\begin{table}[!h]\n\\centering \\caption{Mean of the distribution of the standard score $z(\\lambda_{i})$: $\\mu_{z(\\lambda_{i})} \\pm \\sigma_{\\mu_{z(\\lambda_{i})}}$}\n\\begin{tabular}{|c|cccccc|}\n\\hline\n");
					fprintf(NumFile, "$p_{T}$ [GeV] & $\\mu_{z(\\lambda_{\\vartheta})}$ & $\\mu_{z(\\lambda_{\\varphi})}$ &  $\\mu_{z(\\lambda_{\\vartheta \\varphi})}$ & $\\mu_{z(\\tilde{\\lambda})}$ & $\\mu_{z(\\lambda^*_{\\vartheta})}$ & $\\mu_{z(\\lambda^*_{\\varphi})}$ \\\\\n");
				}
				if(iTab==3){
					fprintf(NumFile, "\\begin{table}[!h]\n\\centering \\caption{R.m.s. of the distribution of the standard score $z(\\lambda_{i})$: $\\sigma_{z(\\lambda_{i})} \\pm \\sigma_{\\sigma_{z(\\lambda_{i})}}$}\n\\begin{tabular}{|c|cccccc|}\n\\hline\n");
					fprintf(NumFile, "$p_{T}$ [GeV] & $\\sigma_{z(\\lambda_{\\vartheta})}$ & $\\sigma_{z(\\lambda_{\\varphi})}$ &  $\\sigma_{z(\\lambda_{\\vartheta \\varphi})}$ & $\\sigma_{z(\\tilde{\\lambda})}$ & $\\sigma_{z(\\lambda^*_{\\vartheta})}$ & $\\sigma_{z(\\lambda^*_{\\varphi})}$ \\\\\n");
				}

				for(int iFrame=1; iFrame<4; iFrame++){


					int rap=0;
					for(int rapBin = rapBinMin; rapBin < rapBinMax+1; rapBin++) {

						if(iFrame==1) {sprintf(framerap,"\\hline \\multicolumn{7}{|c|}{CS frame, $%1.1f < |y| < %1.1f$}\\\\ \\hline \\rule{0pt}{4mm}\n",onia::rapForPTRange[rapBin-1],onia::rapForPTRange[rapBin]); fprintf(NumFile,framerap);}
						if(iFrame==2) {sprintf(framerap,"\\hline \\multicolumn{7}{|c|}{HX frame, $%1.1f < |y| < %1.1f$}\\\\ \\hline \\rule{0pt}{4mm}\n",onia::rapForPTRange[rapBin-1],onia::rapForPTRange[rapBin]); fprintf(NumFile,framerap);}
						if(iFrame==3) {sprintf(framerap,"\\hline \\multicolumn{7}{|c|}{PX frame, $%1.1f < |y| < %1.1f$}\\\\ \\hline \\rule{0pt}{4mm}\n",onia::rapForPTRange[rapBin-1],onia::rapForPTRange[rapBin]); fprintf(NumFile,framerap);}

						int pt=0;
						for(int ptBin = ptBinMin; ptBin < ptBinMax+1; ptBin++) {
							int cpm=0;
							for(int cpmBin = cpmBinMin; cpmBin < cpmBinMax+1; cpmBin++) {


							if(iTab==1){
								if(iFrame==1){
									lth_tab=param_lth_CS_mean[rap][pt][cpm]			-lambda_theta_injected_CS[rap][pt][cpm];  	ltherr_tab=param_lth_CS_sigma[rap][pt][cpm];
									lph_tab=param_lph_CS_mean[rap][pt][cpm]			-lambda_phi_injected_CS[rap][pt][cpm]; 		lpherr_tab=param_lph_CS_sigma[rap][pt][cpm];
									ltp_tab=param_ltp_CS_mean[rap][pt][cpm]			-lambda_thetaphi_injected_CS[rap][pt][cpm]; 	ltperr_tab=param_ltp_CS_sigma[rap][pt][cpm];
									ltilde_tab=param_ltilde_CS_mean[rap][pt][cpm]	-lambda_tilde_injected_CS[rap][pt][cpm]; 	ltildeerr_tab=param_ltilde_CS_sigma[rap][pt][cpm];
									lthstar_tab=param_lthstar_CS_mean[rap][pt][cpm]	-lambda_thetastar_injected_CS[rap][pt][cpm];	lthstarerr_tab=param_lthstar_CS_sigma[rap][pt][cpm];
									lphstar_tab=param_lphstar_CS_mean[rap][pt][cpm]	-lambda_phistar_injected_CS[rap][pt][cpm];	lphstarerr_tab=param_lphstar_CS_sigma[rap][pt][cpm];
								}
								if(iFrame==2){
									lth_tab=param_lth_HX_mean[rap][pt][cpm]			-lambda_theta_injected_HX[rap][pt][cpm];  	ltherr_tab=param_lth_HX_sigma[rap][pt][cpm];
									lph_tab=param_lph_HX_mean[rap][pt][cpm]			-lambda_phi_injected_HX[rap][pt][cpm]; 		lpherr_tab=param_lph_HX_sigma[rap][pt][cpm];
									ltp_tab=param_ltp_HX_mean[rap][pt][cpm]			-lambda_thetaphi_injected_HX[rap][pt][cpm]; 	ltperr_tab=param_ltp_HX_sigma[rap][pt][cpm];
									ltilde_tab=param_ltilde_HX_mean[rap][pt][cpm]	-lambda_tilde_injected_HX[rap][pt][cpm]; 	ltildeerr_tab=param_ltilde_HX_sigma[rap][pt][cpm];
									lthstar_tab=param_lthstar_HX_mean[rap][pt][cpm]	-lambda_thetastar_injected_HX[rap][pt][cpm];	lthstarerr_tab=param_lthstar_HX_sigma[rap][pt][cpm];
									lphstar_tab=param_lphstar_HX_mean[rap][pt][cpm]	-lambda_phistar_injected_HX[rap][pt][cpm];	lphstarerr_tab=param_lphstar_HX_sigma[rap][pt][cpm];
								}
								if(iFrame==3){
									lth_tab=param_lth_PX_mean[rap][pt][cpm]			-lambda_theta_injected_PX[rap][pt][cpm];  	ltherr_tab=param_lth_PX_sigma[rap][pt][cpm];
									lph_tab=param_lph_PX_mean[rap][pt][cpm]			-lambda_phi_injected_PX[rap][pt][cpm]; 		lpherr_tab=param_lph_PX_sigma[rap][pt][cpm];
									ltp_tab=param_ltp_PX_mean[rap][pt][cpm]			-lambda_thetaphi_injected_PX[rap][pt][cpm]; 	ltperr_tab=param_ltp_PX_sigma[rap][pt][cpm];
									ltilde_tab=param_ltilde_PX_mean[rap][pt][cpm]	-lambda_tilde_injected_PX[rap][pt][cpm]; 	ltildeerr_tab=param_ltilde_PX_sigma[rap][pt][cpm];
									lthstar_tab=param_lthstar_PX_mean[rap][pt][cpm]	-lambda_thetastar_injected_PX[rap][pt][cpm];	lthstarerr_tab=param_lthstar_PX_sigma[rap][pt][cpm];
									lphstar_tab=param_lphstar_PX_mean[rap][pt][cpm]	-lambda_phistar_injected_PX[rap][pt][cpm];	lphstarerr_tab=param_lphstar_PX_sigma[rap][pt][cpm];
								}
							}

							if(iTab==2){
								if(iFrame==1){
									lth_tab=pull_lth_CS_mean[rap][pt][cpm]; ltherr_tab=pull_lth_CS_meanerr[rap][pt][cpm];
									lph_tab=pull_lph_CS_mean[rap][pt][cpm]; lpherr_tab=pull_lph_CS_meanerr[rap][pt][cpm];
									ltp_tab=pull_ltp_CS_mean[rap][pt][cpm]; ltperr_tab=pull_ltp_CS_meanerr[rap][pt][cpm];
									ltilde_tab=pull_ltilde_CS_mean[rap][pt][cpm]; ltildeerr_tab=pull_ltilde_CS_meanerr[rap][pt][cpm];
									lthstar_tab=pull_lthstar_CS_mean[rap][pt][cpm]; lthstarerr_tab=pull_lthstar_CS_meanerr[rap][pt][cpm];
									lphstar_tab=pull_lphstar_CS_mean[rap][pt][cpm]; lphstarerr_tab=pull_lphstar_CS_meanerr[rap][pt][cpm];
								}
								if(iFrame==2){
									lth_tab=pull_lth_HX_mean[rap][pt][cpm]; ltherr_tab=pull_lth_HX_meanerr[rap][pt][cpm];
									lph_tab=pull_lph_HX_mean[rap][pt][cpm]; lpherr_tab=pull_lph_HX_meanerr[rap][pt][cpm];
									ltp_tab=pull_ltp_HX_mean[rap][pt][cpm]; ltperr_tab=pull_ltp_HX_meanerr[rap][pt][cpm];
									ltilde_tab=pull_ltilde_HX_mean[rap][pt][cpm]; ltildeerr_tab=pull_ltilde_HX_meanerr[rap][pt][cpm];
									lthstar_tab=pull_lthstar_HX_mean[rap][pt][cpm]; lthstarerr_tab=pull_lthstar_HX_meanerr[rap][pt][cpm];
									lphstar_tab=pull_lphstar_HX_mean[rap][pt][cpm]; lphstarerr_tab=pull_lphstar_HX_meanerr[rap][pt][cpm];
								}
								if(iFrame==3){
									lth_tab=pull_lth_PX_mean[rap][pt][cpm]; ltherr_tab=pull_lth_PX_meanerr[rap][pt][cpm];
									lph_tab=pull_lph_PX_mean[rap][pt][cpm]; lpherr_tab=pull_lph_PX_meanerr[rap][pt][cpm];
									ltp_tab=pull_ltp_PX_mean[rap][pt][cpm]; ltperr_tab=pull_ltp_PX_meanerr[rap][pt][cpm];
									ltilde_tab=pull_ltilde_PX_mean[rap][pt][cpm]; ltildeerr_tab=pull_ltilde_PX_meanerr[rap][pt][cpm];
									lthstar_tab=pull_lthstar_PX_mean[rap][pt][cpm]; lthstarerr_tab=pull_lthstar_PX_meanerr[rap][pt][cpm];
									lphstar_tab=pull_lphstar_PX_mean[rap][pt][cpm]; lphstarerr_tab=pull_lphstar_PX_meanerr[rap][pt][cpm];
								}
							}

							if(iTab==3){
								if(iFrame==1){
									lth_tab=pull_lth_CS_sigma[rap][pt][cpm]; ltherr_tab=pull_lth_CS_sigmaerr[rap][pt][cpm];
									lph_tab=pull_lph_CS_sigma[rap][pt][cpm]; lpherr_tab=pull_lph_CS_sigmaerr[rap][pt][cpm];
									ltp_tab=pull_ltp_CS_sigma[rap][pt][cpm]; ltperr_tab=pull_ltp_CS_sigmaerr[rap][pt][cpm];
									ltilde_tab=pull_ltilde_CS_sigma[rap][pt][cpm]; ltildeerr_tab=pull_ltilde_CS_sigmaerr[rap][pt][cpm];
									lthstar_tab=pull_lthstar_CS_sigma[rap][pt][cpm]; lthstarerr_tab=pull_lthstar_CS_sigmaerr[rap][pt][cpm];
									lphstar_tab=pull_lphstar_CS_sigma[rap][pt][cpm]; lphstarerr_tab=pull_lphstar_CS_sigmaerr[rap][pt][cpm];
								}
								if(iFrame==2){
									lth_tab=pull_lth_HX_sigma[rap][pt][cpm]; ltherr_tab=pull_lth_HX_sigmaerr[rap][pt][cpm];
									lph_tab=pull_lph_HX_sigma[rap][pt][cpm]; lpherr_tab=pull_lph_HX_sigmaerr[rap][pt][cpm];
									ltp_tab=pull_ltp_HX_sigma[rap][pt][cpm]; ltperr_tab=pull_ltp_HX_sigmaerr[rap][pt][cpm];
									ltilde_tab=pull_ltilde_HX_sigma[rap][pt][cpm]; ltildeerr_tab=pull_ltilde_HX_sigmaerr[rap][pt][cpm];
									lthstar_tab=pull_lthstar_HX_sigma[rap][pt][cpm]; lthstarerr_tab=pull_lthstar_HX_sigmaerr[rap][pt][cpm];
									lphstar_tab=pull_lphstar_HX_sigma[rap][pt][cpm]; lphstarerr_tab=pull_lphstar_HX_sigmaerr[rap][pt][cpm];
								}
								if(iFrame==3){
									lth_tab=pull_lth_PX_sigma[rap][pt][cpm]; ltherr_tab=pull_lth_PX_sigmaerr[rap][pt][cpm];
									lph_tab=pull_lph_PX_sigma[rap][pt][cpm]; lpherr_tab=pull_lph_PX_sigmaerr[rap][pt][cpm];
									ltp_tab=pull_ltp_PX_sigma[rap][pt][cpm]; ltperr_tab=pull_ltp_PX_sigmaerr[rap][pt][cpm];
									ltilde_tab=pull_ltilde_PX_sigma[rap][pt][cpm]; ltildeerr_tab=pull_ltilde_PX_sigmaerr[rap][pt][cpm];
									lthstar_tab=pull_lthstar_PX_sigma[rap][pt][cpm]; lthstarerr_tab=pull_lthstar_PX_sigmaerr[rap][pt][cpm];
									lphstar_tab=pull_lphstar_PX_sigma[rap][pt][cpm]; lphstarerr_tab=pull_lphstar_PX_sigmaerr[rap][pt][cpm];
								}
							}



							if(emptyBin[rap][pt][cpm])fprintf(NumFile, "%1.0f--%1.0f   &  $nan$  & $nan$  &  $nan$ &  $nan$  & $nan$  &  $nan$ \\\\\n",onia::pTRange[rapBin][ptBin-1],onia::pTRange[rapBin][ptBin]);

							else fprintf(NumFile, "%1.0f--%1.0f   &  $%1.3f \\pm %1.3f$  & $%1.3f \\pm %1.3f$  &  $%1.3f \\pm %1.3f$ &  $%1.3f \\pm %1.3f$  & $%1.3f \\pm %1.3f$  &  $%1.3f \\pm %1.3f$ \\\\\n",onia::pTRange[rapBin][ptBin-1],onia::pTRange[rapBin][ptBin],lth_tab,ltherr_tab   ,lph_tab,lpherr_tab   ,ltp_tab,ltperr_tab   ,ltilde_tab,ltildeerr_tab   ,lthstar_tab,lthstarerr_tab   ,lphstar_tab,lphstarerr_tab);

							cpm++;
							}
							pt++;
						}



						rap++;

					}//end rapBin

				}//end iFrame



				fprintf(NumFile, "\\hline\n");
				fprintf(NumFile, "\\end{tabular}\n");
				fprintf(NumFile, "\\label{tab:syst_acceptance}\n");
				fprintf(NumFile, "\\end{table}\n");
				fprintf(NumFile, "\n");

			}//end iTab
		}
		else{


			nTables=1;
			for(int iTab=1; iTab<nTables+1;iTab++){

				fprintf(NumFile, "\n\n\n\n");

				if(iTab==1){
					fprintf(NumFile, "\\begin{table}[!h]\n\\centering \\caption{Parameter Results}\n\\begin{tabular}{|c|cccccc|}\n\\hline\n");
					fprintf(NumFile, "$p_{T}$ [GeV] & $\\lambda_{\\vartheta}$ & $\\lambda_{\\varphi}$ &  $\\lambda_{\\vartheta \\varphi}$ & $\\tilde{\\lambda}$ & $\\lambda^*_{\\vartheta}$ & $\\lambda^*_{\\varphi}$ \\\\\n");
				}

				for(int iFrame=1; iFrame<4; iFrame++){


					int rap=0;
					for(int rapBin = rapBinMin; rapBin < rapBinMax+1; rapBin++) {

						if(iFrame==1) {sprintf(framerap,"\\hline \\multicolumn{7}{|c|}{CS frame, $%1.1f < |y| < %1.1f$}\\\\ \\hline \\rule{0pt}{4mm}\n",onia::rapForPTRange[rapBin-1],onia::rapForPTRange[rapBin]); fprintf(NumFile,framerap);}
						if(iFrame==2) {sprintf(framerap,"\\hline \\multicolumn{7}{|c|}{HX frame, $%1.1f < |y| < %1.1f$}\\\\ \\hline \\rule{0pt}{4mm}\n",onia::rapForPTRange[rapBin-1],onia::rapForPTRange[rapBin]); fprintf(NumFile,framerap);}
						if(iFrame==3) {sprintf(framerap,"\\hline \\multicolumn{7}{|c|}{PX frame, $%1.1f < |y| < %1.1f$}\\\\ \\hline \\rule{0pt}{4mm}\n",onia::rapForPTRange[rapBin-1],onia::rapForPTRange[rapBin]); fprintf(NumFile,framerap);}

						int pt=0;
						for(int ptBin = ptBinMin; ptBin < ptBinMax+1; ptBin++) {
						  int cpm=0;
						  for(int cpmBin = cpmBinMin; cpmBin < cpmBinMax+1; cpmBin++) {


							if(iTab==1){
								if(iFrame==1){
									lth_tab=param_lth_CS_mean[rap][pt][cpm]			;   ltherr_tab=param_lth_CS_sigma_low[rap][pt][cpm];            ltherr_high_tab=param_lth_CS_sigma_high[rap][pt][cpm];
									lph_tab=param_lph_CS_mean[rap][pt][cpm]			;   lpherr_tab=param_lph_CS_sigma_low[rap][pt][cpm];            lpherr_high_tab=param_lph_CS_sigma_high[rap][pt][cpm];
									ltp_tab=param_ltp_CS_mean[rap][pt][cpm]			;   ltperr_tab=param_ltp_CS_sigma_low[rap][pt][cpm];            ltperr_high_tab=param_ltp_CS_sigma_high[rap][pt][cpm];
									ltilde_tab=param_ltilde_CS_mean[rap][pt][cpm]	;   ltildeerr_tab=param_ltilde_CS_sigma_low[rap][pt][cpm];      ltildeerr_high_tab=param_ltilde_CS_sigma_high[rap][pt][cpm];
									lthstar_tab=param_lthstar_CS_mean[rap][pt][cpm]	;   lthstarerr_tab=param_lthstar_CS_sigma_low[rap][pt][cpm];    lthstarerr_high_tab=param_lthstar_CS_sigma_high[rap][pt][cpm];
									lphstar_tab=param_lphstar_CS_mean[rap][pt][cpm]	;   lphstarerr_tab=param_lphstar_CS_sigma_low[rap][pt][cpm];    lphstarerr_high_tab=param_lphstar_CS_sigma_high[rap][pt][cpm];
								}
								if(iFrame==2){
									lth_tab=param_lth_HX_mean[rap][pt][cpm]			;   ltherr_tab=param_lth_HX_sigma_low[rap][pt][cpm];            ltherr_high_tab=param_lth_HX_sigma_high[rap][pt][cpm];
									lph_tab=param_lph_HX_mean[rap][pt][cpm]			;   lpherr_tab=param_lph_HX_sigma_low[rap][pt][cpm];            lpherr_high_tab=param_lph_HX_sigma_high[rap][pt][cpm];
									ltp_tab=param_ltp_HX_mean[rap][pt][cpm]			;   ltperr_tab=param_ltp_HX_sigma_low[rap][pt][cpm];            ltperr_high_tab=param_ltp_HX_sigma_high[rap][pt][cpm];
									ltilde_tab=param_ltilde_HX_mean[rap][pt][cpm]	;   ltildeerr_tab=param_ltilde_HX_sigma_low[rap][pt][cpm];      ltildeerr_high_tab=param_ltilde_HX_sigma_high[rap][pt][cpm];
									lthstar_tab=param_lthstar_HX_mean[rap][pt][cpm]	;   lthstarerr_tab=param_lthstar_HX_sigma_low[rap][pt][cpm];    lthstarerr_high_tab=param_lthstar_HX_sigma_high[rap][pt][cpm];
									lphstar_tab=param_lphstar_HX_mean[rap][pt][cpm]	;   lphstarerr_tab=param_lphstar_HX_sigma_low[rap][pt][cpm];    lphstarerr_high_tab=param_lphstar_HX_sigma_high[rap][pt][cpm];
								}
								if(iFrame==3){
									lth_tab=param_lth_PX_mean[rap][pt][cpm]			;   ltherr_tab=param_lth_PX_sigma_low[rap][pt][cpm];            ltherr_high_tab=param_lth_PX_sigma_high[rap][pt][cpm];
									lph_tab=param_lph_PX_mean[rap][pt][cpm]			;   lpherr_tab=param_lph_PX_sigma_low[rap][pt][cpm];            lpherr_high_tab=param_lph_PX_sigma_high[rap][pt][cpm];
									ltp_tab=param_ltp_PX_mean[rap][pt][cpm]			;   ltperr_tab=param_ltp_PX_sigma_low[rap][pt][cpm];            ltperr_high_tab=param_ltp_PX_sigma_high[rap][pt][cpm];
									ltilde_tab=param_ltilde_PX_mean[rap][pt][cpm]	;   ltildeerr_tab=param_ltilde_PX_sigma_low[rap][pt][cpm];      ltildeerr_high_tab=param_ltilde_PX_sigma_high[rap][pt][cpm];
									lthstar_tab=param_lthstar_PX_mean[rap][pt][cpm]	;   lthstarerr_tab=param_lthstar_PX_sigma_low[rap][pt][cpm];    lthstarerr_high_tab=param_lthstar_PX_sigma_high[rap][pt][cpm];
									lphstar_tab=param_lphstar_PX_mean[rap][pt][cpm]	;   lphstarerr_tab=param_lphstar_PX_sigma_low[rap][pt][cpm];    lphstarerr_high_tab=param_lphstar_PX_sigma_high[rap][pt][cpm];
								}
							}

							if(emptyBin[rap][pt][cpm])fprintf(NumFile, "%1.0f--%1.0f   &  $nan$  & $nan$  &  $nan$ &  $nan$  & $nan$  &  $nan$ \\\\\n",onia::pTRange[rapBin][ptBin-1],onia::pTRange[rapBin][ptBin]);

							else fprintf(NumFile, "%1.0f--%1.0f   &  $%1.3f _{-%1.3f}^{+%1.3f} $  & $%1.3f _{-%1.3f}^{+%1.3f}$  &  $%1.3f _{-%1.3f}^{+%1.3f}$ &  $%1.3f _{-%1.3f}^{+%1.3f}$  & $%1.3f _{-%1.3f}^{+%1.3f}$  &  $%1.3f _{-%1.3f}^{+%1.3f}$ \\\\\n",onia::pTRange[rapBin][ptBin-1],onia::pTRange[rapBin][ptBin],lth_tab,ltherr_tab,ltherr_high_tab   ,lph_tab,lpherr_tab,lpherr_high_tab   ,ltp_tab,ltperr_tab,ltperr_high_tab   ,ltilde_tab,ltildeerr_tab,ltildeerr_high_tab   ,lthstar_tab,lthstarerr_tab,lthstarerr_high_tab   ,lphstar_tab,lphstarerr_tab,lphstarerr_high_tab);

							cpm++;
							}
							pt++;
						}



						rap++;

					}//end rapBin

				}//end iFrame



				fprintf(NumFile, "\\hline\n");
				fprintf(NumFile, "\\end{tabular}\n");
				fprintf(NumFile, "\\label{tab:syst_acceptance}\n");
				fprintf(NumFile, "\\end{table}\n");
				fprintf(NumFile, "\n");

			}//end iTab



		}

		fprintf(NumFile, "\\end{document}");

		fclose(NumFile);

		if(RealData){
			sprintf(dirstruct,"%s",dirstructPuffer);
		}
		cout<<"End of Plots"<<endl;
		return 0;
}
