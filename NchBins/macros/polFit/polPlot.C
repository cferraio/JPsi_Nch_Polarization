#include "Riostream.h"
#include "TSystem.h"
#include "TString.h"
#include "TROOT.h"
#include "TTree.h"
#include "TMath.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "TF1.h"
#include "TF2.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TProfile.h"
#include "TCanvas.h"
#include "TColor.h"
#include "TLine.h"
#include "TPolyLine.h"
#include "TEllipse.h"
#include "TStyle.h"
#include "TFile.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TPaveStats.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TFrame.h"

// binning of costh - phi plots
const int nbin_cth = 80; //40.;//80
const int nbin_ph  = 72; //36.;//72

// extremes and binning of lambda plots
double lth_min = -1.1;
double lth_max =  1.1;
double lth_step_1D = 0.02;
double lth_step_2D = 0.02;

double lph_min = -1.1;
double lph_max =  1.1;
double lph_step_1D = 0.02;
double lph_step_2D = 0.02;

double ltp_min = -1.1;
double ltp_max =  1.1;
double ltp_step_1D = 0.02;
double ltp_step_2D = 0.02;

double ContourHeightForIntegral;
double FractionInsideContour;

// function to calculate height of the contour of a 2D distribution for a certain confidence level
inline double contourHeight2D ( TH2D *h, double confidenceLevel ) {
	int Nx = h->GetXaxis()->GetNbins();
	int Ny = h->GetYaxis()->GetNbins();

	double totSum = h->GetSum();
	double targetSum = confidenceLevel * totSum;
	double maxHeight = h->GetMaximum();
	double step = 0.001*maxHeight;

	double tempHeight = 0.;
	double tempSum = totSum;

	while ( tempSum > targetSum && tempHeight < maxHeight ) {
		tempHeight += step;
		tempSum = 0.;
		for ( int ix = 0; ix < Nx; ix++ ) {
			for ( int iy = 0; iy < Ny; iy++ ) {
				double binContent = h->GetBinContent(ix,iy);
				if ( binContent > tempHeight ) tempSum += binContent;
			}
		}
	}
	//  cout<<"tempHeight = "<<tempHeight<<endl;
	return tempHeight;
}

// function to set the 99% and 68% C.L. contours of a 2D histogram
inline void setContourHistogram ( TH2D *h ) {
	double cont0 = contourHeight2D( h, 0.997 );
	double cont1 = contourHeight2D( h, 0.955 );
	double cont2 = contourHeight2D( h, 0.683 );
	/*  h->SetContour(3);
			h->SetContourLevel(0,cont0);
			h->SetContourLevel(1,cont1);
			h->SetContourLevel(2,cont2);
	 */
	h->SetContour(2);
	h->SetContourLevel(0,cont0);
	h->SetContourLevel(1,cont2);

	/*  h->SetContour(1);
	//  h->SetContourLevel(0,cont1);
	h->SetContourLevel(0,cont2);
	ContourHeightForIntegral=cont2;

	// Calc area in contour:
	int iCalcNbins=0;
	int iCalcNbinsInContour=0;
	for(int iCalcIntX=0;iCalcIntX<h->GetNbinsX();iCalcIntX++){
	for(int iCalcIntY=1;iCalcIntY<h->GetNbinsY();iCalcIntY++){
	iCalcNbins++;
	if(h->GetBinContent(iCalcIntX,iCalcIntY)>ContourHeightForIntegral) iCalcNbinsInContour++;
	}
	}
	FractionInsideContour=double(iCalcNbinsInContour)/double(iCalcNbins);
	 */

}

inline void setContourHistogram1Sig ( TH2D *h ) {
	double cont0 = contourHeight2D( h, 0.683 );
	h->SetContour(1);
	h->SetContourLevel(0,cont0);
}
inline void setContourHistogram3Sig ( TH2D *h ) {
	double cont0 = contourHeight2D( h, 0.997 );
	h->SetContour(1);
	h->SetContourLevel(0,cont0);
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

void PlotPosterior(int ptBin, int rapBin, char xAxisTitle[200], char filename[200], TH1* histo, double MPV, double MPVerrLow, double MPVerrHigh){

	double plotborder = 1.1; // how much space to leave beyond the maximum of the plots
	double plotMax = plotborder * histo->GetMaximum();

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

	TLegend* plotLegend=new TLegend(0.725,0.7,0.95,0.9);
	plotLegend->SetFillColor(kWhite);
	//		plotLegend->SetTextFont(72);
	plotLegend->SetTextSize(0.03);
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
	histo->SetMaximum(plotMax);
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

void polPlot(Char_t *dirstruct = "OutputDirectory_Default",
		Char_t *TreeBinID = "TreeBinID_Default",
		bool RealData=false,
		int MPValgo=0,
		bool scalePlots=false,
		int nTotalFits=1,
		int nState=999,
		double ptlow=999.,
		double pthigh=999.,
		double raplow=999.,
		double raphigh=999.){

	cout<<"/////////////////////////////////"<<endl;
	cout<<"running polPlot.C ........///////"<<endl;
	cout<<"/////////////////////////////////"<<endl;

	gROOT->Reset();
	gROOT->SetStyle("Plain");
	gStyle->SetOptFit(1111);
	gStyle->SetOptStat(0);

	char filename [500];

	cout<<"nState: "<<nState<<endl;

	// input file:

	sprintf(filename,"%s/results.root",dirstruct);
	if(RealData) sprintf(filename,"%s/results_%s.root",dirstruct,TreeBinID);
	cout<<filename<<endl;
	TFile* results = new TFile(filename);

	sprintf(dirstruct,"%s/Figures",dirstruct); gSystem->mkdir(dirstruct);

	int n_entries, n_step;

	// input ntuples with PDFs of the anisotropy parameters

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


	// create CS histograms

	int l_bins=100;
	if(RealData) l_bins=200;

	TH1D* h_lth_CS = new TH1D( "h_lth_CS", "", l_bins, lambdaCS->GetMinimum("lth"), lambdaCS->GetMaximum("lth") );
	TH1D* h_lph_CS = new TH1D( "h_lph_CS", "", l_bins, lambdaCS->GetMinimum("lph"), lambdaCS->GetMaximum("lph") );
	TH1D* h_ltp_CS = new TH1D( "h_ltp_CS", "", l_bins, lambdaCS->GetMinimum("ltp"), lambdaCS->GetMaximum("ltp") );
	TH1D* h_ltilde_CS = new TH1D( "h_ltilde_CS", "", l_bins, lambdaCS->GetMinimum("ltilde"), lambdaCS->GetMaximum("ltilde") );

	TH2D* h_lph_vs_lth_CS = new TH2D( "h_lph_vs_lth_CS", "", int((lth_max-lth_min)/lth_step_2D), lth_min, lth_max, int((lph_max-lph_min)/lph_step_2D), lph_min, lph_max );
	TH2D* h_ltp_vs_lth_CS = new TH2D( "h_ltp_vs_lth_CS", "", int((lth_max-lth_min)/lth_step_2D), lth_min, lth_max, int((ltp_max-ltp_min)/ltp_step_2D), ltp_min, ltp_max );
	TH2D* h_ltp_vs_lph_CS = new TH2D( "h_ltp_vs_lph_CS", "", int((lph_max-lph_min)/lph_step_2D), lph_min, lph_max, int((ltp_max-ltp_min)/ltp_step_2D), ltp_min, ltp_max );

	TH2D* h_lphstar_vs_lthstar_CS = new TH2D( "h_lphstar_vs_lthstar_CS", "", int((lth_max-lth_min)/lth_step_2D), lth_min, lth_max, int((lph_max-lph_min)/lph_step_2D), lph_min, lph_max );


	// loop over entries in the CS ntuple

	n_entries = int( lambdaCS->GetEntries() );

	cout << endl;
	cout << "Reading distribution of CS parameters (" << n_entries << ") entries"<< endl;
	cout << "-------------------------------------------------------------" << endl;
	cout << "Progress: ";

	n_step = n_entries/50;

	for ( int i_entry = 0; i_entry < n_entries; i_entry++ ) {

		if ( i_entry%n_step == 0 ) cout << "X";

		lambdaCS->GetEvent( i_entry );

		h_lth_CS->Fill( lth_CS );
		h_lph_CS->Fill( lph_CS );
		h_ltp_CS->Fill( ltp_CS );

		h_lph_vs_lth_CS->Fill( lth_CS, lph_CS );
		h_ltp_vs_lth_CS->Fill( lth_CS, ltp_CS );
		h_ltp_vs_lph_CS->Fill( lph_CS, ltp_CS );

		h_ltilde_CS->Fill( ltilde_CS );
		h_lphstar_vs_lthstar_CS->Fill( lthstar_CS, lphstar_CS );

	}

	// end of loop over entries in the CS ntuple

	cout << endl << endl;



	// create HX histograms

	TH1D* h_lth_HX = new TH1D( "h_lth_HX", "", l_bins, lambdaHX->GetMinimum("lth"), lambdaHX->GetMaximum("lth") );
	TH1D* h_lph_HX = new TH1D( "h_lph_HX", "", l_bins, lambdaHX->GetMinimum("lph"), lambdaHX->GetMaximum("lph") );
	TH1D* h_ltp_HX = new TH1D( "h_ltp_HX", "", l_bins, lambdaHX->GetMinimum("ltp"), lambdaHX->GetMaximum("ltp") );
	TH1D* h_ltilde_HX = new TH1D( "h_ltilde_HX", "", l_bins, lambdaHX->GetMinimum("ltilde"), lambdaHX->GetMaximum("ltilde") );

	TH2D* h_lph_vs_lth_HX = new TH2D( "h_lph_vs_lth_HX", "", int((lth_max-lth_min)/lth_step_2D), lth_min, lth_max, int((lph_max-lph_min)/lph_step_2D), lph_min, lph_max );
	TH2D* h_ltp_vs_lth_HX = new TH2D( "h_ltp_vs_lth_HX", "", int((lth_max-lth_min)/lth_step_2D), lth_min, lth_max, int((ltp_max-ltp_min)/ltp_step_2D), ltp_min, ltp_max );
	TH2D* h_ltp_vs_lph_HX = new TH2D( "h_ltp_vs_lph_HX", "", int((lph_max-lph_min)/lph_step_2D), lph_min, lph_max, int((ltp_max-ltp_min)/ltp_step_2D), ltp_min, ltp_max );

	TH2D* h_lphstar_vs_lthstar_HX = new TH2D( "h_lphstar_vs_lthstar_HX", "", int((lth_max-lth_min)/lth_step_2D), lth_min, lth_max, int((lph_max-lph_min)/lph_step_2D), lph_min, lph_max );


	// loop over entries in the HX ntuple

	n_entries = int( lambdaHX->GetEntries() );

	cout << endl;
	cout << "Reading distribution of HX parameters (" << n_entries << ") entries"<< endl;
	cout << "-------------------------------------------------------------" << endl;
	cout << "Progress: ";

	n_step = n_entries/50;

	for ( int i_entry = 0; i_entry < n_entries; i_entry++ ) {

		if ( i_entry%n_step == 0 ) cout << "X";

		lambdaHX->GetEvent( i_entry );

		h_lth_HX->Fill( lth_HX );
		h_lph_HX->Fill( lph_HX );
		h_ltp_HX->Fill( ltp_HX );

		h_lph_vs_lth_HX->Fill( lth_HX, lph_HX );
		h_ltp_vs_lth_HX->Fill( lth_HX, ltp_HX );
		h_ltp_vs_lph_HX->Fill( lph_HX, ltp_HX );

		h_ltilde_HX->Fill( ltilde_HX );
		h_lphstar_vs_lthstar_HX->Fill( lthstar_HX, lphstar_HX );

	}
	// end of loop over entries in the HX ntuple

	cout << endl << endl;


	// create PX histograms

	TH1D* h_lth_PX = new TH1D( "h_lth_PX", "", l_bins, lambdaPX->GetMinimum("lth"), lambdaPX->GetMaximum("lth") );
	TH1D* h_lph_PX = new TH1D( "h_lph_PX", "", l_bins, lambdaPX->GetMinimum("lph"), lambdaPX->GetMaximum("lph") );
	TH1D* h_ltp_PX = new TH1D( "h_ltp_PX", "", l_bins, lambdaPX->GetMinimum("ltp"), lambdaPX->GetMaximum("ltp") );
	TH1D* h_ltilde_PX = new TH1D( "h_ltilde_PX", "", l_bins, lambdaPX->GetMinimum("ltilde"), lambdaPX->GetMaximum("ltilde") );

	TH2D* h_lph_vs_lth_PX = new TH2D( "h_lph_vs_lth_PX", "", int((lth_max-lth_min)/lth_step_2D), lth_min, lth_max, int((lph_max-lph_min)/lph_step_2D), lph_min, lph_max );
	TH2D* h_ltp_vs_lth_PX = new TH2D( "h_ltp_vs_lth_PX", "", int((lth_max-lth_min)/lth_step_2D), lth_min, lth_max, int((ltp_max-ltp_min)/ltp_step_2D), ltp_min, ltp_max );
	TH2D* h_ltp_vs_lph_PX = new TH2D( "h_ltp_vs_lph_PX", "", int((lph_max-lph_min)/lph_step_2D), lph_min, lph_max, int((ltp_max-ltp_min)/ltp_step_2D), ltp_min, ltp_max );

	TH2D* h_lphstar_vs_lthstar_PX = new TH2D( "h_lphstar_vs_lthstar_PX", "", int((lth_max-lth_min)/lth_step_2D), lth_min, lth_max, int((lph_max-lph_min)/lph_step_2D), lph_min, lph_max );

	// loop over entries in the PX ntuple

	n_entries = int( lambdaPX->GetEntries() );

	cout << endl;
	cout << "Reading distribution of PX parameters (" << n_entries << ") entries"<< endl;
	cout << "-------------------------------------------------------------" << endl;
	cout << "Progress: ";

	n_step = n_entries/50;

	for ( int i_entry = 0; i_entry < n_entries; i_entry++ ) {

		if ( i_entry%n_step == 0 ) cout << "X";

		lambdaPX->GetEvent( i_entry );

		h_lth_PX->Fill( lth_PX );
		h_lph_PX->Fill( lph_PX );
		h_ltp_PX->Fill( ltp_PX );

		h_lph_vs_lth_PX->Fill( lth_PX, lph_PX );
		h_ltp_vs_lth_PX->Fill( lth_PX, ltp_PX );
		h_ltp_vs_lph_PX->Fill( lph_PX, ltp_PX );

		h_ltilde_PX->Fill( ltilde_PX );
		h_lphstar_vs_lthstar_PX->Fill( lthstar_PX, lphstar_PX );

	}
	// end of loop over entries in the PX ntuple



	cout << endl << endl;

	bool CMSprelim=false;
	bool TreeBinIDTextAdd=false;

	TLatex *CentralsText1_2D;
	TLatex *CentralsText2_2D;
	TLatex *CentralsText3_2D;
	TLatex *CentralsTextNumbering_2D;
	TLatex *TreeBinIDText1_2D;
	TLatex *TreeBinIDText2_2D;
	TLatex *TreeBinIDText3_2D;
	double xCentrals2D;
	char text[200];
	double CentralsFontSize;



	TGraph* graph_lph_vs_lth_CS_1Sig       = NULL;
	TGraph* graph_lph_vs_lth_CS_3Sig       = NULL;
	TGraph* graph_ltp_vs_lph_CS_1Sig       = NULL;
	TGraph* graph_ltp_vs_lph_CS_3Sig       = NULL;
	TGraph* graph_ltp_vs_lth_CS_1Sig       = NULL;
	TGraph* graph_ltp_vs_lth_CS_3Sig       = NULL;

	TCanvas *cBuffer = new TCanvas("cBuffer", "cBuffer", 10, 28, 580,571);
	cBuffer->Range(-237.541,-66.47556,187.377,434.8609);
	cBuffer->SetFillColor(0);
	cBuffer->SetBorderMode(0);
	cBuffer->SetBorderSize(0);
	cBuffer->SetLeftMargin(0.1354167);
	cBuffer->SetRightMargin(0.01736111);
	cBuffer->SetTopMargin(0.01841621);
	cBuffer->SetBottomMargin(0.1325967);
	cBuffer->SetFrameBorderMode(0);

	Int_t i, j;
	Double_t x0, y0, z0;
	Double_t contours[6];
	TObjArray *conts;
	TList* contLevel = NULL;
	TGraph* curv     = NULL;
	TGraph* gc       = NULL;
	Int_t TotalConts = 0;

	setContourHistogram1Sig ( h_lph_vs_lth_CS );
	h_lph_vs_lth_CS->Draw("CONT Z LIST");
	cBuffer->Update();
	conts = (TObjArray*)gROOT->GetListOfSpecials()->FindObject("contours");
	TotalConts = conts->GetSize();
	for(i = 0; i < TotalConts; i++) contLevel = (TList*)conts->At(i);
	for(i = 0; i < TotalConts; i++){
		contLevel = (TList*)conts->At(i);
		curv = (TGraph*)contLevel->First();
		for(j = 0; j < contLevel->GetSize(); j++){
			gc = (TGraph*)curv->Clone();
			curv = (TGraph*)contLevel->After(curv);
		}
	}
	cBuffer->Update();
	graph_lph_vs_lth_CS_1Sig= (TGraph*)gc->Clone();

	setContourHistogram3Sig ( h_lph_vs_lth_CS );
	h_lph_vs_lth_CS->Draw("CONT Z LIST");
	cBuffer->Update();
	conts = (TObjArray*)gROOT->GetListOfSpecials()->FindObject("contours");
	TotalConts = conts->GetSize();
	for(i = 0; i < TotalConts; i++) contLevel = (TList*)conts->At(i);
	for(i = 0; i < TotalConts; i++){
		contLevel = (TList*)conts->At(i);
		curv = (TGraph*)contLevel->First();
		for(j = 0; j < contLevel->GetSize(); j++){
			gc = (TGraph*)curv->Clone();
			curv = (TGraph*)contLevel->After(curv);
		}
	}
	cBuffer->Update();
	graph_lph_vs_lth_CS_3Sig= (TGraph*)gc->Clone();







	setContourHistogram1Sig ( h_ltp_vs_lph_CS );
	h_ltp_vs_lph_CS->Draw("CONT Z LIST");
	cBuffer->Update();
	conts = (TObjArray*)gROOT->GetListOfSpecials()->FindObject("contours");
	TotalConts = conts->GetSize();
	for(i = 0; i < TotalConts; i++) contLevel = (TList*)conts->At(i);
	for(i = 0; i < TotalConts; i++){
		contLevel = (TList*)conts->At(i);
		curv = (TGraph*)contLevel->First();
		for(j = 0; j < contLevel->GetSize(); j++){
			gc = (TGraph*)curv->Clone();
			curv = (TGraph*)contLevel->After(curv);
		}
	}
	cBuffer->Update();
	graph_ltp_vs_lph_CS_1Sig= (TGraph*)gc->Clone();

	setContourHistogram3Sig ( h_ltp_vs_lph_CS );
	h_ltp_vs_lph_CS->Draw("CONT Z LIST");
	cBuffer->Update();
	conts = (TObjArray*)gROOT->GetListOfSpecials()->FindObject("contours");
	TotalConts = conts->GetSize();
	for(i = 0; i < TotalConts; i++) contLevel = (TList*)conts->At(i);
	for(i = 0; i < TotalConts; i++){
		contLevel = (TList*)conts->At(i);
		curv = (TGraph*)contLevel->First();
		for(j = 0; j < contLevel->GetSize(); j++){
			gc = (TGraph*)curv->Clone();
			curv = (TGraph*)contLevel->After(curv);
		}
	}
	cBuffer->Update();
	graph_ltp_vs_lph_CS_3Sig= (TGraph*)gc->Clone();





	setContourHistogram1Sig ( h_ltp_vs_lth_CS );
	h_ltp_vs_lth_CS->Draw("CONT Z LIST");
	cBuffer->Update();
	conts = (TObjArray*)gROOT->GetListOfSpecials()->FindObject("contours");
	TotalConts = conts->GetSize();
	for(i = 0; i < TotalConts; i++) contLevel = (TList*)conts->At(i);
	for(i = 0; i < TotalConts; i++){
		contLevel = (TList*)conts->At(i);
		curv = (TGraph*)contLevel->First();
		for(j = 0; j < contLevel->GetSize(); j++){
			gc = (TGraph*)curv->Clone();
			curv = (TGraph*)contLevel->After(curv);
		}
	}
	cBuffer->Update();
	graph_ltp_vs_lth_CS_1Sig= (TGraph*)gc->Clone();

	setContourHistogram3Sig ( h_ltp_vs_lth_CS );
	h_ltp_vs_lth_CS->Draw("CONT Z LIST");
	cBuffer->Update();
	conts = (TObjArray*)gROOT->GetListOfSpecials()->FindObject("contours");
	TotalConts = conts->GetSize();
	for(i = 0; i < TotalConts; i++) contLevel = (TList*)conts->At(i);
	for(i = 0; i < TotalConts; i++){
		contLevel = (TList*)conts->At(i);
		curv = (TGraph*)contLevel->First();
		for(j = 0; j < contLevel->GetSize(); j++){
			gc = (TGraph*)curv->Clone();
			curv = (TGraph*)contLevel->After(curv);
		}
	}
	cBuffer->Update();
	graph_ltp_vs_lth_CS_3Sig= (TGraph*)gc->Clone();






	// canvas

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


	double NewXAxisMin=lth_min;
	double NewXAxisMax=lth_max;
	double NewYAxisMin=lph_min;
	double NewYAxisMax=lph_max;

	double CentralsTextyMin=0;
	double CentralsTextxMin=0;
	double CentralsTextyDelta=0.05;
	double CentralsTextxDelta=0.05;

	////////////////////////////////////////////////////////////////
	// 2D plot lph vs lth

	bool TomChange=true;
	bool TomZoom=false;
	bool TomZoomButNeverthelessSpendHoursToWorkOnThatBoundaries=false;

	char DrawContourStyle[200];
	sprintf(DrawContourStyle,"cont2,same");
	int LineWidth1Sig=4;
	int LineWidth3Sig=2;
	int LineStyleCS=2;
	int LineStyleHX=3; //1
	int LineStylePX=1; //3
	int LineColorCS=418;
	int LineColorHX=632;
	int LineColorPX=600;


	//TomChange
	if(TomChange&&TomZoom){
		lth_min=-0.7;
		lth_max=0.7;
		lph_min=-0.5;
		lph_max=0.5;
	}

	TH2D* h_lph_vs_lth = new TH2D( "h_lph_vs_lth", "", int((lth_max-lth_min)/lth_step_2D), lth_min, lth_max, int((lph_max-lph_min)/lph_step_2D), lph_min, lph_max );
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
	if(!TomZoom)background_ph_vs_th->Draw("f same");

	double x_triangle_ph_vs_th[3] = {-1.,  1., 1.};
	double y_triangle_ph_vs_th[3] = { 0., -1., 1.};
	TPolyLine *triangle_ph_vs_th = new TPolyLine( 3, x_triangle_ph_vs_th, y_triangle_ph_vs_th );
	triangle_ph_vs_th->SetFillColor(kWhite);
	triangle_ph_vs_th->SetLineStyle(0);
	if(!TomZoom)triangle_ph_vs_th->Draw("f same");

	TLine* lh_p1 = new TLine( lth_min, 1., lth_max, 1. );
	lh_p1->SetLineWidth( 1 );
	lh_p1->SetLineStyle( 2 );
	lh_p1->SetLineColor( kGray+2 );
	if(!TomZoom)lh_p1->Draw( "same" );

	TLine* lh_m1 = new TLine( lth_min, -1., lth_max, -1. );
	lh_m1->SetLineWidth( 1 );
	lh_m1->SetLineStyle( 2 );
	lh_m1->SetLineColor( kGray+2 );
	if(!TomZoom)lh_m1->Draw( "same" );

	TLine* lh_0 = new TLine( lth_min, 0., lth_max, 0. );
	lh_0->SetLineWidth( 1 );
	lh_0->SetLineStyle( 2 );
	lh_0->SetLineColor( kGray+2 );
	lh_0->Draw( "same" );

	TLine* lv_p1 = new TLine( 1., lph_min, 1., lph_max );
	lv_p1->SetLineWidth( 1 );
	lv_p1->SetLineStyle( 2 );
	lv_p1->SetLineColor( kGray+2 );
	if(!TomZoom)lv_p1->Draw( "same" );

	TLine* lv_m1 = new TLine( -1., lph_min, -1., lph_max );
	lv_m1->SetLineWidth( 1 );
	lv_m1->SetLineStyle( 2 );
	lv_m1->SetLineColor( kGray+2 );
	if(!TomZoom)lv_m1->Draw( "same" );

	TLine* lv_0 = new TLine( 0., lph_min, 0., lph_max );
	lv_0->SetLineWidth( 1 );
	lv_0->SetLineStyle( 2 );
	lv_0->SetLineColor( kGray+2 );
	lv_0->Draw( "same" );


	// CS frame

	/*  setContourHistogram ( h_lph_vs_lth_CS );
			h_lph_vs_lth_CS->SetLineColor( kGreen+2 );
			h_lph_vs_lth_CS->SetLineWidth( LineWidth3Sig );
			h_lph_vs_lth_CS->SetLineStyle( LineStyleCS  );
			h_lph_vs_lth_CS->Draw( DrawContourStyle );
	 */


	graph_lph_vs_lth_CS_1Sig->SetLineWidth(LineWidth1Sig);
	graph_lph_vs_lth_CS_1Sig->SetLineColor(LineColorCS);
	graph_lph_vs_lth_CS_1Sig->SetLineStyle(LineStyleCS);
	graph_lph_vs_lth_CS_1Sig->Draw("C");
	graph_lph_vs_lth_CS_3Sig->SetLineWidth(LineWidth3Sig);
	graph_lph_vs_lth_CS_3Sig->SetLineColor(LineColorCS);
	graph_lph_vs_lth_CS_3Sig->SetLineStyle(LineStyleCS);
	graph_lph_vs_lth_CS_3Sig->Draw("C");


	h_lph_vs_lth_CS->SetLineColor( kGreen+2 );
	h_lph_vs_lth_CS->SetLineWidth( LineWidth1Sig );
	h_lph_vs_lth_CS->SetLineStyle( LineStyleCS  );
	//  h_lph_vs_lth_CS->Draw( DrawContourStyle );
	TH2D *h_lph_vs_lth_CS_clone=(TH2D*)h_lph_vs_lth_CS->Clone();
	setContourHistogram3Sig ( h_lph_vs_lth_CS_clone );
	h_lph_vs_lth_CS_clone->SetLineColor( kGreen+2 );
	h_lph_vs_lth_CS_clone->SetLineWidth( LineWidth3Sig );
	h_lph_vs_lth_CS_clone->SetLineStyle( LineStyleCS  );
	//  h_lph_vs_lth_CS_clone->Draw( DrawContourStyle );

	// HX frame

	/*  setContourHistogram ( h_lph_vs_lth_HX );
			h_lph_vs_lth_HX->SetLineColor( kBlue );
			h_lph_vs_lth_HX->SetLineWidth( LineWidth3Sig );
			h_lph_vs_lth_HX->SetLineStyle( LineStyleHX );
			h_lph_vs_lth_HX->Draw( DrawContourStyle );
	 */
	setContourHistogram1Sig ( h_lph_vs_lth_HX );
	h_lph_vs_lth_HX->SetLineColor( kBlue );
	h_lph_vs_lth_HX->SetLineWidth( LineWidth1Sig );
	h_lph_vs_lth_HX->SetLineStyle( LineStyleHX  );
	if(!TomChange) h_lph_vs_lth_HX->Draw( DrawContourStyle );
	TH2D *h_lph_vs_lth_HX_clone=(TH2D*)h_lph_vs_lth_HX->Clone();
	setContourHistogram3Sig ( h_lph_vs_lth_HX_clone );
	h_lph_vs_lth_HX_clone->SetLineColor( kBlue );
	h_lph_vs_lth_HX_clone->SetLineWidth( LineWidth3Sig );
	h_lph_vs_lth_HX_clone->SetLineStyle( LineStyleHX  );
	if(!TomChange) h_lph_vs_lth_HX_clone->Draw( DrawContourStyle );

	// PX frame

	/*  setContourHistogram ( h_lph_vs_lth_PX );
			h_lph_vs_lth_PX->SetLineColor( kRed );
			h_lph_vs_lth_PX->SetLineWidth( LineWidth3Sig );
			h_lph_vs_lth_PX->SetLineStyle( LineStylePX );
			if(!TomChange) h_lph_vs_lth_PX->Draw( DrawContourStyle );
	 */

	setContourHistogram1Sig ( h_lph_vs_lth_PX );
	h_lph_vs_lth_PX->SetLineColor( kRed ); 
	h_lph_vs_lth_PX->SetLineWidth( LineWidth1Sig );
	h_lph_vs_lth_PX->SetLineStyle( LineStylePX  );
	h_lph_vs_lth_PX->Draw( DrawContourStyle );
	TH2D *h_lph_vs_lth_PX_clone=(TH2D*)h_lph_vs_lth_PX->Clone();
	setContourHistogram3Sig ( h_lph_vs_lth_PX_clone );
	h_lph_vs_lth_PX_clone->SetLineColor( kRed );
	h_lph_vs_lth_PX_clone->SetLineWidth( LineWidth3Sig );
	h_lph_vs_lth_PX_clone->SetLineStyle( LineStylePX  );
	h_lph_vs_lth_PX_clone->Draw( DrawContourStyle );

	char legendentry[200];

	TLegend* plotLegend;
	if(!TomChange){
		plotLegend=new TLegend(0.15,0.8,0.25,0.95);
		plotLegend->SetFillColor(kWhite);
		plotLegend->SetTextSize(0.035);
		plotLegend->SetBorderSize(1);
		sprintf(legendentry,"CS");
		plotLegend->AddEntry(h_lph_vs_lth_CS,legendentry,"l");
		sprintf(legendentry,"HX");
		plotLegend->AddEntry(h_lph_vs_lth_HX,legendentry,"l");
		sprintf(legendentry,"PX");
		plotLegend->AddEntry(h_lph_vs_lth_PX,legendentry,"l");
		plotLegend->Draw(); plotLegend->Draw();
	}

	if(TomChange){
		plotLegend=new TLegend(0.195,0.735,0.475,0.93);
		if(TomZoom) 	 plotLegend=new TLegend(0.65,0.735,0.93,0.93);
		plotLegend->SetFillColor(kWhite);
		plotLegend->SetTextSize(0.035);
		plotLegend->SetBorderSize(1);

		TH2D* h_lph_vs_lth_CS_legend1 = (TH2D*)h_lph_vs_lth_CS->Clone();
		TH2D* h_lph_vs_lth_CS_legend2 = (TH2D*)h_lph_vs_lth_CS_clone->Clone();
		TH2D* h_lph_vs_lth_PX_legend1 = (TH2D*)h_lph_vs_lth_PX->Clone();
		TH2D* h_lph_vs_lth_PX_legend2 = (TH2D*)h_lph_vs_lth_PX_clone->Clone();

		sprintf(legendentry,"CS, 68.3%% CL");
		plotLegend->AddEntry(h_lph_vs_lth_CS_legend1,legendentry,"l");
		sprintf(legendentry,"CS, 99.7%% CL");
		plotLegend->AddEntry(h_lph_vs_lth_CS_legend2,legendentry,"l");
		sprintf(legendentry,"PX, 68.3%% CL");
		plotLegend->AddEntry(h_lph_vs_lth_PX_legend1,legendentry,"l");
		sprintf(legendentry,"PX, 99.7%% CL");
		plotLegend->AddEntry(h_lph_vs_lth_PX_legend2,legendentry,"l");
		plotLegend->Draw();
	}

	if(!TomChange){
		cout<<"DRAW CMS preliminary Latex"<<endl;
		xCentrals2D=-0.75;
		CentralsFontSize=0.035;
		sprintf(text,"CMS Preliminary");
		if(!CMSprelim) sprintf(text,"CMS");
		CentralsText1_2D = new TLatex(xCentrals2D,0.91,text);
		CentralsText1_2D->SetTextSize(CentralsFontSize);
		sprintf(text,"L = 4.9 fb^{-1}");
		CentralsText2_2D = new TLatex(xCentrals2D,0.65,text);
		//	    if(!CMSprelim) CentralsText2_2D = new TLatex(xCentrals2D,0.72,text);
		CentralsText2_2D->SetTextSize(CentralsFontSize);
		sprintf(text,"pp    #sqrt{s} = 7 TeV");
		CentralsText3_2D = new TLatex(xCentrals2D,0.78,text);
		//	    if(!CMSprelim) CentralsText3_2D = new TLatex(xCentrals2D,0.85,text);
		CentralsText3_2D->SetTextSize(CentralsFontSize);
	}

	if(TomChange){
		cout<<"DRAW CMS preliminary Latex"<<endl;
		xCentrals2D=lth_min+(lth_max-lth_min)*0.10;
		CentralsFontSize=0.035;
		sprintf(text,"CMS Preliminary");
		if(!CMSprelim) sprintf(text,"CMS");
		CentralsText1_2D = new TLatex(xCentrals2D,lph_min+(lph_max-lph_min)*0.21,text);
		CentralsText1_2D->SetTextSize(CentralsFontSize);
		sprintf(text,"L = 4.9 fb^{-1}");
		CentralsText2_2D = new TLatex(xCentrals2D,lph_min+(lph_max-lph_min)*0.11,text);
		CentralsText2_2D->SetTextSize(CentralsFontSize);
		sprintf(text,"pp    #sqrt{s} = 7 TeV");
		CentralsText3_2D = new TLatex(xCentrals2D,lph_min+(lph_max-lph_min)*0.16,text);
		CentralsText3_2D->SetTextSize(CentralsFontSize);
		if(TomZoom){
			cout<<"DRAW CMS preliminary Latex"<<endl;
			xCentrals2D=lth_min+(lth_max-lth_min)*0.65;
			CentralsFontSize=0.035;
			sprintf(text,"CMS Preliminary");
			if(!CMSprelim) sprintf(text,"CMS");
			CentralsText1_2D = new TLatex(xCentrals2D,lph_min+(lph_max-lph_min)*0.21,text);
			CentralsText1_2D->SetTextSize(CentralsFontSize);
			sprintf(text,"L = 4.9 fb^{-1}");
			CentralsText2_2D = new TLatex(xCentrals2D,lph_min+(lph_max-lph_min)*0.11,text);
			CentralsText2_2D->SetTextSize(CentralsFontSize);
			sprintf(text,"pp    #sqrt{s} = 7 TeV");
			CentralsText3_2D = new TLatex(xCentrals2D,lph_min+(lph_max-lph_min)*0.16,text);
			CentralsText3_2D->SetTextSize(CentralsFontSize);
		}

		sprintf(text,"(a)");
		CentralsTextNumbering_2D = new TLatex(xCentrals2D,lph_min+(lph_max-lph_min)*0.31,text);
		CentralsTextNumbering_2D->SetTextSize(CentralsFontSize);
		CentralsTextNumbering_2D->Draw( "same" );

	}
	CentralsText1_2D->Draw( "same" );
	CentralsText2_2D->Draw( "same" );
	CentralsText3_2D->Draw( "same" );

	sprintf(text,"#psi(%dS)",nState);
	TreeBinIDText1_2D = new TLatex(xCentrals2D,-0.73,text);
	TreeBinIDText1_2D->SetTextSize(CentralsFontSize);
	if(raplow<0.3) sprintf(text,"|y| < 0.6");
	if(raplow>0.3 && raplow<0.9) sprintf(text,"0.6 < |y| < 1.2");
	if(raplow>0.9) sprintf(text,"1.2 < |y| < 1.8");

	TreeBinIDText2_2D = new TLatex(xCentrals2D,-0.84,text);
	TreeBinIDText2_2D->SetTextSize(CentralsFontSize);
	sprintf(text,"%1.0f < p_{T} < %1.0f GeV",ptlow,pthigh);
	TreeBinIDText3_2D = new TLatex(xCentrals2D,-0.95,text);
	TreeBinIDText3_2D->SetTextSize(CentralsFontSize);

	if(TreeBinIDTextAdd){
		TreeBinIDText1_2D->Draw( "same" );
		TreeBinIDText2_2D->Draw( "same" );
		TreeBinIDText3_2D->Draw( "same" );
	}

	sprintf(filename,"%s/lph_vs_lth_%s.pdf",dirstruct,TreeBinID);
	c1->Print( filename );
	sprintf(filename,"%s/lph_vs_lth_%s.C",dirstruct,TreeBinID);
	c1->Print( filename );
	sprintf(filename,"%s/lph_vs_lth_%s.jpg",dirstruct,TreeBinID);
	c1->Print( filename );




	////////////////////////////////////////////////////////////////
	// 2D plot lphstar vs lthstar

	TH2D* h_lphstar_vs_lthstar = (TH2D*)h_lphstar_vs_lthstar_CS->Clone();
	h_lphstar_vs_lthstar->SetName("h_lphstar_vs_lthstar");
	h_lphstar_vs_lthstar->Reset();

	h_lphstar_vs_lthstar->GetXaxis()->SetTitle("#lambda_{#vartheta}^{*}");
	h_lphstar_vs_lthstar->GetXaxis()->SetLabelOffset(0.028);
	h_lphstar_vs_lthstar->GetXaxis()->SetTitleSize(0.05);
	h_lphstar_vs_lthstar->GetXaxis()->SetTickLength(-0.03);
	h_lphstar_vs_lthstar->GetXaxis()->SetTitleOffset(1.15);
	h_lphstar_vs_lthstar->GetYaxis()->SetTitle("#lambda_{#varphi}^{*}");
	h_lphstar_vs_lthstar->GetYaxis()->SetLabelOffset(0.032);
	h_lphstar_vs_lthstar->GetYaxis()->SetTitleSize(0.05);
	h_lphstar_vs_lthstar->GetYaxis()->SetTickLength(-0.03);
	h_lphstar_vs_lthstar->GetYaxis()->SetTitleOffset(1.3);
	h_lphstar_vs_lthstar->Draw("");

	background_ph_vs_th->Draw("f same");
	triangle_ph_vs_th->Draw("f same");

	lh_p1->Draw( "same" );
	lh_m1->Draw( "same" );
	lh_0->Draw( "same" );
	lv_p1->Draw( "same" );
	lv_m1->Draw( "same" );
	lv_0->Draw( "same" );


	// CS frame

	setContourHistogram ( h_lphstar_vs_lthstar_CS );
	h_lphstar_vs_lthstar_CS->SetLineColor( kGreen+2 );
	h_lphstar_vs_lthstar_CS->SetLineWidth( LineWidth3Sig );
	h_lphstar_vs_lthstar_CS->SetLineStyle( LineStyleCS   );
	h_lphstar_vs_lthstar_CS->Draw( "cont2, same" );

	// HX frame

	setContourHistogram ( h_lphstar_vs_lthstar_HX );
	h_lphstar_vs_lthstar_HX->SetLineColor( kBlue );
	h_lphstar_vs_lthstar_HX->SetLineWidth( LineWidth3Sig );
	h_lphstar_vs_lthstar_HX->SetLineStyle( LineStyleHX   );
	h_lphstar_vs_lthstar_HX->Draw( "cont2, same" );

	// PX frame

	setContourHistogram ( h_lphstar_vs_lthstar_PX );
	h_lphstar_vs_lthstar_PX->SetLineColor( kRed );
	h_lphstar_vs_lthstar_PX->SetLineWidth( LineWidth3Sig );
	h_lphstar_vs_lthstar_PX->SetLineStyle( LineStylePX   );
	h_lphstar_vs_lthstar_PX->Draw( "cont2, same" );

	plotLegend->Draw();


	CentralsText1_2D->Draw( "same" );
	CentralsText2_2D->Draw( "same" );
	CentralsText3_2D->Draw( "same" );

	if(TreeBinIDTextAdd){
		TreeBinIDText1_2D->Draw( "same" );
		TreeBinIDText2_2D->Draw( "same" );
		TreeBinIDText3_2D->Draw( "same" );
	}

	sprintf(filename,"%s/lphstar_vs_lthstar_%s.pdf",dirstruct,TreeBinID);
	c1->Print( filename );
	sprintf(filename,"%s/lphstar_vs_lthstar_%s.C",dirstruct,TreeBinID);
	c1->Print( filename );
	sprintf(filename,"%s/lphstar_vs_lthstar_%s.jpg",dirstruct,TreeBinID);
	c1->Print( filename );




	////////////////////////////////////////////////////////////////
	// 2D plot ltp vs lth

	if(TomChange&&TomZoom){
		lth_min=-0.7;
		lth_max=0.7;
		ltp_min=-0.5;
		ltp_max=0.5;
	}

	TH2D* h_ltp_vs_lth = new TH2D( "h_ltp_vs_lth", "", int((lth_max-lth_min)/lth_step_2D), lth_min, lth_max, int((ltp_max-ltp_min)/ltp_step_2D), ltp_min, ltp_max );
	//  TH2D* h_ltp_vs_lth = (TH2D*)h_ltp_vs_lth_CS->Clone();
	h_ltp_vs_lth->SetName("h_ltp_vs_lth");
	h_ltp_vs_lth->Reset();

	h_ltp_vs_lth->GetXaxis()->SetTitle("#lambda_{#vartheta}");
	h_ltp_vs_lth->GetXaxis()->SetLabelOffset(0.028);
	h_ltp_vs_lth->GetXaxis()->SetTitleSize(0.05);
	h_ltp_vs_lth->GetXaxis()->SetTickLength(-0.03);
	h_ltp_vs_lth->GetXaxis()->SetTitleOffset(1.15);
	h_ltp_vs_lth->GetYaxis()->SetTitle("#lambda_{#vartheta#varphi}");
	h_ltp_vs_lth->GetYaxis()->SetLabelOffset(0.032);
	h_ltp_vs_lth->GetYaxis()->SetTitleSize(0.05);
	h_ltp_vs_lth->GetYaxis()->SetTickLength(-0.03);
	h_ltp_vs_lth->GetYaxis()->SetTitleOffset(1.3);
	h_ltp_vs_lth->Draw("");

	double x_background_tp_vs_th[4] = { lth_min+0.01, lth_min+0.01, lth_max-0.01, lth_max-0.01 };
	double y_background_tp_vs_th[4] = { ltp_min+0.01, ltp_max-0.01, ltp_max-0.01, ltp_min+0.01 };
	TPolyLine *background_tp_vs_th = new TPolyLine( 4, x_background_tp_vs_th, y_background_tp_vs_th );
	background_tp_vs_th->SetFillColor(kGray);
	background_tp_vs_th->SetLineStyle(0);
	if(!TomZoom) background_tp_vs_th->Draw("f same");

	TEllipse *ellipse_tp_vs_th = new TEllipse( 0., 0., 1., sqrt(2.)/2. );
	ellipse_tp_vs_th->SetFillColor(kWhite);
	ellipse_tp_vs_th->SetLineStyle(0);
	if(!TomZoom) ellipse_tp_vs_th->Draw("f same");


	lh_p1 = new TLine( lth_min, sqrt(2.)/2., lth_max, sqrt(2.)/2. );
	lh_p1->SetLineWidth( 1 );
	lh_p1->SetLineStyle( 2 );
	lh_p1->SetLineColor( kGray+2 );
	if(!TomZoom) lh_p1->Draw( "same" );

	lh_m1 = new TLine( lth_min, -sqrt(2.)/2., lth_max, -sqrt(2.)/2. );
	lh_m1->SetLineWidth( 1 );
	lh_m1->SetLineStyle( 2 );
	lh_m1->SetLineColor( kGray+2 );
	if(!TomZoom) lh_m1->Draw( "same" );

	lh_0 = new TLine( lth_min, 0., lth_max, 0. );
	lh_0->SetLineWidth( 1 );
	lh_0->SetLineStyle( 2 );
	lh_0->SetLineColor( kGray+2 );
	lh_0->Draw( "same" );

	lv_p1 = new TLine( 1., ltp_min, 1., ltp_max );
	lv_p1->SetLineWidth( 1 );
	lv_p1->SetLineStyle( 2 );
	lv_p1->SetLineColor( kGray+2 );
	if(!TomZoom) lv_p1->Draw( "same" );

	lv_m1 = new TLine( -1., ltp_min, -1., ltp_max );
	lv_m1->SetLineWidth( 1 );
	lv_m1->SetLineStyle( 2 );
	lv_m1->SetLineColor( kGray+2 );
	if(!TomZoom) lv_m1->Draw( "same" );

	lv_0 = new TLine( 0., ltp_min, 0., ltp_max );
	lv_0->SetLineWidth( 1 );
	lv_0->SetLineStyle( 2 );
	lv_0->SetLineColor( kGray+2 );
	lv_0->Draw( "same" );


	// CS frame

	/*  setContourHistogram ( h_ltp_vs_lth_CS );
			h_ltp_vs_lth_CS->SetLineColor( kGreen+2 );
			h_ltp_vs_lth_CS->SetLineWidth( LineWidth3Sig );
			h_ltp_vs_lth_CS->SetLineStyle( LineStyleCS  );
			h_ltp_vs_lth_CS->Draw( DrawContourStyle );
	 */

	graph_ltp_vs_lth_CS_1Sig->SetLineWidth(LineWidth1Sig);
	graph_ltp_vs_lth_CS_1Sig->SetLineColor(LineColorCS);
	graph_ltp_vs_lth_CS_1Sig->SetLineStyle(LineStyleCS);
	graph_ltp_vs_lth_CS_1Sig->Draw("C");
	graph_ltp_vs_lth_CS_3Sig->SetLineWidth(LineWidth3Sig);
	graph_ltp_vs_lth_CS_3Sig->SetLineColor(LineColorCS);
	graph_ltp_vs_lth_CS_3Sig->SetLineStyle(LineStyleCS);
	graph_ltp_vs_lth_CS_3Sig->Draw("C");

	setContourHistogram1Sig ( h_ltp_vs_lth_CS );
	h_ltp_vs_lth_CS->SetLineColor( kGreen+2 );
	h_ltp_vs_lth_CS->SetLineWidth( LineWidth1Sig );
	h_ltp_vs_lth_CS->SetLineStyle( LineStyleCS  );
	//  h_ltp_vs_lth_CS->Draw( DrawContourStyle );
	TH2D *h_ltp_vs_lth_CS_clone=(TH2D*)h_ltp_vs_lth_CS->Clone();
	setContourHistogram3Sig ( h_ltp_vs_lth_CS_clone );
	h_ltp_vs_lth_CS_clone->SetLineColor( kGreen+2 );
	h_ltp_vs_lth_CS_clone->SetLineWidth( LineWidth3Sig );
	h_ltp_vs_lth_CS_clone->SetLineStyle( LineStyleCS  );
	//  h_ltp_vs_lth_CS_clone->Draw( DrawContourStyle );

	// HX frame

	/*  setContourHistogram ( h_ltp_vs_lth_HX );
			h_ltp_vs_lth_HX->SetLineColor( kBlue );
			h_ltp_vs_lth_HX->SetLineWidth( LineWidth3Sig );
			h_ltp_vs_lth_HX->SetLineStyle( LineStyleHX );
			h_ltp_vs_lth_HX->Draw( DrawContourStyle );
	 */
	setContourHistogram1Sig ( h_ltp_vs_lth_HX );
	h_ltp_vs_lth_HX->SetLineColor( kBlue );
	h_ltp_vs_lth_HX->SetLineWidth( LineWidth1Sig );
	h_ltp_vs_lth_HX->SetLineStyle( LineStyleHX  );
	if(!TomChange) h_ltp_vs_lth_HX->Draw( DrawContourStyle );
	TH2D *h_ltp_vs_lth_HX_clone=(TH2D*)h_ltp_vs_lth_HX->Clone();
	setContourHistogram3Sig ( h_ltp_vs_lth_HX_clone );
	h_ltp_vs_lth_HX_clone->SetLineColor( kBlue );
	h_ltp_vs_lth_HX_clone->SetLineWidth( LineWidth3Sig );
	h_ltp_vs_lth_HX_clone->SetLineStyle( LineStyleHX  );
	if(!TomChange) h_ltp_vs_lth_HX_clone->Draw( DrawContourStyle );

	// PX frame

	/*  setContourHistogram ( h_ltp_vs_lth_PX );
			h_ltp_vs_lth_PX->SetLineColor( kRed );
			h_ltp_vs_lth_PX->SetLineWidth( LineWidth3Sig );
			h_ltp_vs_lth_PX->SetLineStyle( LineStylePX );
			if(!TomChange) h_ltp_vs_lth_PX->Draw( DrawContourStyle );
	 */

	setContourHistogram1Sig ( h_ltp_vs_lth_PX );
	h_ltp_vs_lth_PX->SetLineColor( kRed );
	h_ltp_vs_lth_PX->SetLineWidth( LineWidth1Sig );
	h_ltp_vs_lth_PX->SetLineStyle( LineStylePX  );
	h_ltp_vs_lth_PX->Draw( DrawContourStyle );
	TH2D *h_ltp_vs_lth_PX_clone=(TH2D*)h_ltp_vs_lth_PX->Clone();
	setContourHistogram3Sig ( h_ltp_vs_lth_PX_clone );
	h_ltp_vs_lth_PX_clone->SetLineColor( kRed );
	h_ltp_vs_lth_PX_clone->SetLineWidth( LineWidth3Sig );
	h_ltp_vs_lth_PX_clone->SetLineStyle( LineStylePX  );
	h_ltp_vs_lth_PX_clone->Draw( DrawContourStyle );

	if(!TomChange){
		plotLegend=new TLegend(0.15,0.825,0.25,0.975);
		plotLegend->SetFillColor(kWhite);
		plotLegend->SetTextSize(0.035);
		plotLegend->SetBorderSize(1);
		sprintf(legendentry,"CS");
		plotLegend->AddEntry(h_lph_vs_lth_CS,legendentry,"l");
		sprintf(legendentry,"HX");
		plotLegend->AddEntry(h_lph_vs_lth_HX,legendentry,"l");
		sprintf(legendentry,"PX");
		plotLegend->AddEntry(h_lph_vs_lth_PX,legendentry,"l");
		plotLegend->Draw(); plotLegend->Draw();
		plotLegend->Draw();
	}

	if(TomChange){
		plotLegend=new TLegend(0.195,0.735,0.475,0.93);
		if(TomZoom) 	 plotLegend=new TLegend(0.65,0.735,0.93,0.93);
		plotLegend->SetFillColor(kWhite);
		plotLegend->SetTextSize(0.035);
		plotLegend->SetBorderSize(1);

		TH2D* h_lph_vs_lth_CS_legend1 = (TH2D*)h_lph_vs_lth_CS->Clone();
		TH2D* h_lph_vs_lth_CS_legend2 = (TH2D*)h_lph_vs_lth_CS_clone->Clone();
		TH2D* h_lph_vs_lth_PX_legend1 = (TH2D*)h_lph_vs_lth_PX->Clone();
		TH2D* h_lph_vs_lth_PX_legend2 = (TH2D*)h_lph_vs_lth_PX_clone->Clone();

		sprintf(legendentry,"CS, 68.3%% CL");
		plotLegend->AddEntry(h_lph_vs_lth_CS_legend1,legendentry,"l");
		sprintf(legendentry,"CS, 99.7%% CL");
		plotLegend->AddEntry(h_lph_vs_lth_CS_legend2,legendentry,"l");
		sprintf(legendentry,"PX, 68.3%% CL");
		plotLegend->AddEntry(h_lph_vs_lth_PX_legend1,legendentry,"l");
		sprintf(legendentry,"PX, 99.7%% CL");
		plotLegend->AddEntry(h_lph_vs_lth_PX_legend2,legendentry,"l");
		plotLegend->Draw();
	}

	if(!TomChange){
		cout<<"DRAW CMS preliminary Latex"<<endl;
		xCentrals2D=-0.75;
		CentralsFontSize=0.035;
		sprintf(text,"CMS Preliminary");
		if(!CMSprelim) sprintf(text,"CMS");
		CentralsText1_2D = new TLatex(xCentrals2D,0.98,text);
		CentralsText1_2D->SetTextSize(CentralsFontSize);
		sprintf(text,"L = 4.9 fb^{-1}");
		CentralsText2_2D = new TLatex(xCentrals2D,0.74,text);
		//  if(!CMSprelim) CentralsText2_2D = new TLatex(xCentrals2D,0.81,text);
		CentralsText2_2D->SetTextSize(CentralsFontSize);
		sprintf(text,"pp    #sqrt{s} = 7 TeV");
		CentralsText3_2D = new TLatex(xCentrals2D,0.86,text);
		//  if(!CMSprelim) CentralsText3_2D = new TLatex(xCentrals2D,0.93,text);
		CentralsText3_2D->SetTextSize(CentralsFontSize);
	}

	if(TomChange){
		cout<<"DRAW CMS preliminary Latex"<<endl;
		xCentrals2D=lth_min+(lth_max-lth_min)*0.1;
		CentralsFontSize=0.035;
		sprintf(text,"CMS Preliminary");
		if(!CMSprelim) sprintf(text,"CMS");
		CentralsText1_2D = new TLatex(xCentrals2D,ltp_min+(ltp_max-ltp_min)*0.13,text);
		CentralsText1_2D->SetTextSize(CentralsFontSize);
		sprintf(text,"L = 4.9 fb^{-1}");
		CentralsText2_2D = new TLatex(xCentrals2D,ltp_min+(ltp_max-ltp_min)*0.03,text);
		CentralsText2_2D->SetTextSize(CentralsFontSize);
		sprintf(text,"pp    #sqrt{s} = 7 TeV");
		CentralsText3_2D = new TLatex(xCentrals2D,ltp_min+(ltp_max-ltp_min)*0.08,text);
		CentralsText3_2D->SetTextSize(CentralsFontSize);
		if(TomZoom){
			cout<<"DRAW CMS preliminary Latex"<<endl;
			xCentrals2D=lth_min+(lth_max-lth_min)*0.65;
			CentralsFontSize=0.035;
			sprintf(text,"CMS Preliminary");
			if(!CMSprelim) sprintf(text,"CMS");
			CentralsText1_2D = new TLatex(xCentrals2D,ltp_min+(ltp_max-ltp_min)*0.21,text);
			CentralsText1_2D->SetTextSize(CentralsFontSize);
			sprintf(text,"L = 4.9 fb^{-1}");
			CentralsText2_2D = new TLatex(xCentrals2D,ltp_min+(ltp_max-ltp_min)*0.11,text);
			CentralsText2_2D->SetTextSize(CentralsFontSize);
			sprintf(text,"pp    #sqrt{s} = 7 TeV");
			CentralsText3_2D = new TLatex(xCentrals2D,ltp_min+(ltp_max-ltp_min)*0.16,text);
			CentralsText3_2D->SetTextSize(CentralsFontSize);
		}
	}

	CentralsText1_2D->Draw( "same" );
	CentralsText2_2D->Draw( "same" );
	CentralsText3_2D->Draw( "same" );

	sprintf(text,"#psi(%dS)",nState);
	TreeBinIDText1_2D = new TLatex(xCentrals2D,-0.80,text);
	TreeBinIDText1_2D->SetTextSize(CentralsFontSize);
	if(raplow<0.3) sprintf(text,"|y| < 0.6");
	if(raplow>0.3 && raplow<0.9) sprintf(text,"0.6 < |y| < 1.2");
	if(raplow>0.9) sprintf(text,"1.2 < |y| < 1.8");
	TreeBinIDText2_2D = new TLatex(xCentrals2D,-0.91,text);
	TreeBinIDText2_2D->SetTextSize(CentralsFontSize);
	sprintf(text,"%1.0f < p_{T} < %1.0f GeV",ptlow,pthigh);
	TreeBinIDText3_2D = new TLatex(xCentrals2D,-1.02,text);
	TreeBinIDText3_2D->SetTextSize(CentralsFontSize);

	if(TreeBinIDTextAdd){
		TreeBinIDText1_2D->Draw( "same" );
		TreeBinIDText2_2D->Draw( "same" );
		TreeBinIDText3_2D->Draw( "same" );
	}

	sprintf(filename,"%s/ltp_vs_lth_%s.pdf",dirstruct,TreeBinID);
	c1->Print( filename );
	sprintf(filename,"%s/ltp_vs_lth_%s.C",dirstruct,TreeBinID);
	c1->Print( filename );
	sprintf(filename,"%s/ltp_vs_lth_%s.jpg",dirstruct,TreeBinID);
	c1->Print( filename );


	////////////////////////////////////////////////////////////////
	// 2D plot ltp vs lph

	if(TomChange&&TomZoom){
		ltp_min=-0.5;
		ltp_max=0.5;
		lph_min=-0.5;
		lph_max=0.5;
	}

	TH2D* h_ltp_vs_lph = new TH2D( "h_ltp_vs_lph", "", int((lph_max-lph_min)/lph_step_2D), lph_min, lph_max, int((ltp_max-ltp_min)/ltp_step_2D), ltp_min, ltp_max );
	// TH2D* h_ltp_vs_lph = (TH2D*)h_ltp_vs_lph_CS->Clone();
	h_ltp_vs_lph->SetName("h_ltp_vs_lph");
	h_ltp_vs_lph->Reset();

	h_ltp_vs_lph->GetXaxis()->SetTitle("#lambda_{#varphi}");
	h_ltp_vs_lph->GetXaxis()->SetLabelOffset(0.028);
	h_ltp_vs_lph->GetXaxis()->SetTitleSize(0.05);
	h_ltp_vs_lph->GetXaxis()->SetTickLength(-0.03);
	h_ltp_vs_lph->GetXaxis()->SetTitleOffset(1.15);
	h_ltp_vs_lph->GetYaxis()->SetTitle("#lambda_{#vartheta#varphi}");
	h_ltp_vs_lph->GetYaxis()->SetLabelOffset(0.032);
	h_ltp_vs_lph->GetYaxis()->SetTitleSize(0.05);
	h_ltp_vs_lph->GetYaxis()->SetTickLength(-0.03);
	h_ltp_vs_lph->GetYaxis()->SetTitleOffset(1.3);
	h_ltp_vs_lph->Draw("");

	double x_background_tp_vs_ph[4] = { lph_min+0.01, lph_min+0.01, lph_max-0.01, lph_max-0.01 };
	double y_background_tp_vs_ph[4] = { ltp_min+0.01, ltp_max-0.01, ltp_max-0.01, ltp_min+0.01 };
	TPolyLine *background_tp_vs_ph = new TPolyLine( 4, x_background_tp_vs_ph, y_background_tp_vs_ph );
	background_tp_vs_ph->SetFillColor(kGray);
	background_tp_vs_ph->SetLineStyle(0);
	if(!TomZoom) background_tp_vs_ph->Draw("f same");

	TEllipse *ellipse_tp_vs_ph = new TEllipse( -0.5, 0., 0.5, sqrt(2.)/2. );
	ellipse_tp_vs_ph->SetFillColor(kWhite);
	ellipse_tp_vs_ph->SetLineStyle(0);
	if(!TomZoom) ellipse_tp_vs_ph->Draw("f same");

	double x_triangle_tp_vs_th[3] = {-1./3., -1./3., 1.};
	double y_triangle_tp_vs_th[3] = { 2./3., -2./3., 0.};
	TPolyLine *triangle_tp_vs_th = new TPolyLine( 3, x_triangle_tp_vs_th, y_triangle_tp_vs_th );
	triangle_tp_vs_th->SetFillColor(kWhite);
	triangle_tp_vs_th->SetLineStyle(0);
	if(!TomZoom) triangle_tp_vs_th->Draw("f same");

	lh_p1 = new TLine( lph_min, sqrt(2.)/2., lph_max, sqrt(2.)/2. );
	lh_p1->SetLineWidth( 1 );
	lh_p1->SetLineStyle( 2 );
	lh_p1->SetLineColor( kGray+2 );
	if(!TomZoom) lh_p1->Draw( "same" );

	lh_m1 = new TLine( lph_min, -sqrt(2.)/2., lph_max, -sqrt(2.)/2. );
	lh_m1->SetLineWidth( 1 );
	lh_m1->SetLineStyle( 2 );
	lh_m1->SetLineColor( kGray+2 );
	if(!TomZoom) lh_m1->Draw( "same" );

	lh_0 = new TLine( lph_min, 0., lph_max, 0. );
	lh_0->SetLineWidth( 1 );
	lh_0->SetLineStyle( 2 );
	lh_0->SetLineColor( kGray+2 );
	lh_0->Draw( "same" );

	lv_p1 = new TLine( 1., ltp_min, 1., ltp_max );
	lv_p1->SetLineWidth( 1 );
	lv_p1->SetLineStyle( 2 );
	lv_p1->SetLineColor( kGray+2 );
	if(!TomZoom) lv_p1->Draw( "same" );

	lv_m1 = new TLine( -1., ltp_min, -1., ltp_max );
	lv_m1->SetLineWidth( 1 );
	lv_m1->SetLineStyle( 2 );
	lv_m1->SetLineColor( kGray+2 );
	if(!TomZoom) lv_m1->Draw( "same" );

	lv_0 = new TLine( 0., ltp_min, 0., ltp_max );
	lv_0->SetLineWidth( 1 );
	lv_0->SetLineStyle( 2 );
	lv_0->SetLineColor( kGray+2 );
	lv_0->Draw( "same" );


	// CS frame

	/*  setContourHistogram ( h_ltp_vs_lph_CS );
			h_ltp_vs_lph_CS->SetLineColor( kGreen+2 );
			h_ltp_vs_lph_CS->SetLineWidth( LineWidth3Sig );
			h_ltp_vs_lph_CS->SetLineStyle( LineStyleCS  );
			h_ltp_vs_lph_CS->Draw( DrawContourStyle );
	 */

	graph_ltp_vs_lph_CS_1Sig->SetLineWidth(LineWidth1Sig);
	graph_ltp_vs_lph_CS_1Sig->SetLineColor(LineColorCS);
	graph_ltp_vs_lph_CS_1Sig->SetLineStyle(LineStyleCS);
	graph_ltp_vs_lph_CS_1Sig->Draw("C");
	graph_ltp_vs_lph_CS_3Sig->SetLineWidth(LineWidth3Sig);
	graph_ltp_vs_lph_CS_3Sig->SetLineColor(LineColorCS);
	graph_ltp_vs_lph_CS_3Sig->SetLineStyle(LineStyleCS);
	graph_ltp_vs_lph_CS_3Sig->Draw("C");

	setContourHistogram1Sig ( h_ltp_vs_lph_CS );
	h_ltp_vs_lph_CS->SetLineColor( kGreen+2 );
	h_ltp_vs_lph_CS->SetLineWidth( LineWidth1Sig );
	h_ltp_vs_lph_CS->SetLineStyle( LineStyleCS  );
	//  h_ltp_vs_lph_CS->Draw( DrawContourStyle );
	TH2D *h_ltp_vs_lph_CS_clone=(TH2D*)h_ltp_vs_lph_CS->Clone();
	setContourHistogram3Sig ( h_ltp_vs_lph_CS_clone );
	h_ltp_vs_lph_CS_clone->SetLineColor( kGreen+2 );
	h_ltp_vs_lph_CS_clone->SetLineWidth( LineWidth3Sig );
	h_ltp_vs_lph_CS_clone->SetLineStyle( LineStyleCS  );
	//  h_ltp_vs_lph_CS_clone->Draw( DrawContourStyle );

	// HX frame

	/*  setContourHistogram ( h_ltp_vs_lph_HX );
			h_ltp_vs_lph_HX->SetLineColor( kBlue);
			h_ltp_vs_lph_HX->SetLineWidth( LineWidth3Sig );
			h_ltp_vs_lph_HX->SetLineStyle( LineStyleHX );
			h_ltp_vs_lph_HX->Draw( DrawContourStyle );
	 */
	setContourHistogram1Sig ( h_ltp_vs_lph_HX );
	h_ltp_vs_lph_HX->SetLineColor( kBlue );
	h_ltp_vs_lph_HX->SetLineWidth( LineWidth1Sig );
	h_ltp_vs_lph_HX->SetLineStyle( LineStyleHX  );
	if(!TomChange) h_ltp_vs_lph_HX->Draw( DrawContourStyle );
	TH2D *h_ltp_vs_lph_HX_clone=(TH2D*)h_ltp_vs_lph_HX->Clone();
	setContourHistogram3Sig ( h_ltp_vs_lph_HX_clone );
	h_ltp_vs_lph_HX_clone->SetLineColor( kBlue );
	h_ltp_vs_lph_HX_clone->SetLineWidth( LineWidth3Sig );
	h_ltp_vs_lph_HX_clone->SetLineStyle( LineStyleHX  );
	if(!TomChange) h_ltp_vs_lph_HX_clone->Draw( DrawContourStyle );

	// PX frame

	/*  setContourHistogram ( h_ltp_vs_lph_PX );
			h_ltp_vs_lph_PX->SetLineColor( kRed );
			h_ltp_vs_lph_PX->SetLineWidth( LineWidth3Sig );
			h_ltp_vs_lph_PX->SetLineStyle( LineStylePX );
			if(!TomChange) h_ltp_vs_lph_PX->Draw( DrawContourStyle );
	 */

	setContourHistogram1Sig ( h_ltp_vs_lph_PX );
	h_ltp_vs_lph_PX->SetLineColor( kRed );
	h_ltp_vs_lph_PX->SetLineWidth( LineWidth1Sig );
	h_ltp_vs_lph_PX->SetLineStyle( LineStylePX  );
	h_ltp_vs_lph_PX->Draw( DrawContourStyle );
	TH2D *h_ltp_vs_lph_PX_clone=(TH2D*)h_ltp_vs_lph_PX->Clone();
	setContourHistogram3Sig ( h_ltp_vs_lph_PX_clone );
	h_ltp_vs_lph_PX_clone->SetLineColor( kRed );
	h_ltp_vs_lph_PX_clone->SetLineWidth( LineWidth3Sig );
	h_ltp_vs_lph_PX_clone->SetLineStyle( LineStylePX  );
	h_ltp_vs_lph_PX_clone->Draw( DrawContourStyle );


	if(!TomChange){
		plotLegend->Draw();
	}

	if(TomChange){
		plotLegend=new TLegend(0.65,0.735,0.93,0.93);
		plotLegend->SetFillColor(kWhite);
		plotLegend->SetTextSize(0.035);
		plotLegend->SetBorderSize(1);

		TH2D* h_lph_vs_lth_CS_legend1 = (TH2D*)h_lph_vs_lth_CS->Clone();
		TH2D* h_lph_vs_lth_CS_legend2 = (TH2D*)h_lph_vs_lth_CS_clone->Clone();
		TH2D* h_lph_vs_lth_PX_legend1 = (TH2D*)h_lph_vs_lth_PX->Clone();
		TH2D* h_lph_vs_lth_PX_legend2 = (TH2D*)h_lph_vs_lth_PX_clone->Clone();

		sprintf(legendentry,"CS, 68.3%% CL");
		plotLegend->AddEntry(h_lph_vs_lth_CS_legend1,legendentry,"l");
		sprintf(legendentry,"CS, 99.7%% CL");
		plotLegend->AddEntry(h_lph_vs_lth_CS_legend2,legendentry,"l");
		sprintf(legendentry,"PX, 68.3%% CL");
		plotLegend->AddEntry(h_lph_vs_lth_PX_legend1,legendentry,"l");
		sprintf(legendentry,"PX, 99.7%% CL");
		plotLegend->AddEntry(h_lph_vs_lth_PX_legend2,legendentry,"l");
		plotLegend->Draw();
	}

	if(TomChange){
		cout<<"DRAW CMS preliminary Latex"<<endl;
		xCentrals2D=lph_min+(lph_max-lph_min)*0.65;
		CentralsFontSize=0.035;
		sprintf(text,"CMS Preliminary");
		if(!CMSprelim) sprintf(text,"CMS");
		CentralsText1_2D = new TLatex(xCentrals2D,ltp_min+(ltp_max-ltp_min)*0.13,text);
		CentralsText1_2D->SetTextSize(CentralsFontSize);
		sprintf(text,"L = 4.9 fb^{-1}");
		CentralsText2_2D = new TLatex(xCentrals2D,ltp_min+(ltp_max-ltp_min)*0.03,text);
		CentralsText2_2D->SetTextSize(CentralsFontSize);
		sprintf(text,"pp    #sqrt{s} = 7 TeV");
		CentralsText3_2D = new TLatex(xCentrals2D,ltp_min+(ltp_max-ltp_min)*0.08,text);
		CentralsText3_2D->SetTextSize(CentralsFontSize);
		if(TomZoom){
			cout<<"DRAW CMS preliminary Latex"<<endl;
			xCentrals2D=lph_min+(lph_max-lph_min)*0.65;
			CentralsFontSize=0.035;
			sprintf(text,"CMS Preliminary");
			if(!CMSprelim) sprintf(text,"CMS");
			CentralsText1_2D = new TLatex(xCentrals2D,ltp_min+(ltp_max-ltp_min)*0.21,text);
			CentralsText1_2D->SetTextSize(CentralsFontSize);
			sprintf(text,"L = 4.9 fb^{-1}");
			CentralsText2_2D = new TLatex(xCentrals2D,ltp_min+(ltp_max-ltp_min)*0.11,text);
			CentralsText2_2D->SetTextSize(CentralsFontSize);
			sprintf(text,"pp    #sqrt{s} = 7 TeV");
			CentralsText3_2D = new TLatex(xCentrals2D,ltp_min+(ltp_max-ltp_min)*0.16,text);
			CentralsText3_2D->SetTextSize(CentralsFontSize);
		}
		sprintf(text,"(b)");
		CentralsTextNumbering_2D = new TLatex(lph_min+(lph_max-lph_min)*0.1,0.875,text);
		CentralsTextNumbering_2D->SetTextSize(CentralsFontSize);
		CentralsTextNumbering_2D->Draw( "same" );


	}

	plotLegend->Draw();


	CentralsText1_2D->Draw( "same" );
	CentralsText2_2D->Draw( "same" );
	CentralsText3_2D->Draw( "same" );

	if(TreeBinIDTextAdd){
		TreeBinIDText1_2D->Draw( "same" );
		TreeBinIDText2_2D->Draw( "same" );
		TreeBinIDText3_2D->Draw( "same" );
	}

	sprintf(filename,"%s/ltp_vs_lph_%s.pdf",dirstruct,TreeBinID);
	c1->Print( filename );
	sprintf(filename,"%s/ltp_vs_lph_%s.C",dirstruct,TreeBinID);
	c1->Print( filename );
	sprintf(filename,"%s/ltp_vs_lph_%s.jpg",dirstruct,TreeBinID);
	c1->Print( filename );




	bool PlotJochenStuff=false;

	if(PlotJochenStuff){

		c1 = new TCanvas("c1", "c1", 10, 28, 650,571);
		c1->Range(-237.541,-66.47556,187.377,434.8609);
		c1->SetFillColor(0);
		c1->SetBorderMode(0);
		c1->SetBorderSize(0);
		c1->SetTopMargin(0.01841621);
		c1->SetBottomMargin(0.1325967);
		c1->SetFrameBorderMode(0);
		c1->SetLeftMargin(0.2);
		c1->SetRightMargin(0.075);




		TFile* TestFile1;

		sprintf(filename,"/afs/hephy.at/scratch/k/knuenz/CMSSW_4_2_4_patch2/src/UpsilonPol/macros/polFit/h_lth_lph/Data_TheGreatRun_10B_May10_TESTS_DEFAULT.root");
		TestFile1 = new TFile(filename);
		TH2D* h_lph_vs_lth_Test1 = (TH2D*)TestFile1->Get("h_lth_lph");

		sprintf(filename,"/afs/hephy.at/scratch/k/knuenz/CMSSW_4_2_4_patch2/src/UpsilonPol/macros/polFit/h_lth_lph/Data_TheGreatRun_10B_May10_TESTS_LamthRand.root");
		TestFile1 = new TFile(filename);
		TH2D* h_lph_vs_lth_Test2 = (TH2D*)TestFile1->Get("h_lth_lph");

		sprintf(filename,"/afs/hephy.at/scratch/k/knuenz/CMSSW_4_2_4_patch2/src/UpsilonPol/macros/polFit/h_lth_lph/Data_TheGreatRun_10B_May10_TESTS_LamthP05.root");
		TestFile1 = new TFile(filename);
		TH2D* h_lph_vs_lth_Test3 = (TH2D*)TestFile1->Get("h_lth_lph");

		sprintf(filename,"/afs/hephy.at/scratch/k/knuenz/CMSSW_4_2_4_patch2/src/UpsilonPol/macros/polFit/h_lth_lph/Data_TheGreatRun_10B_May10_TESTS_LamthM05.root");
		TestFile1 = new TFile(filename);
		TH2D* h_lph_vs_lth_Test4 = (TH2D*)TestFile1->Get("h_lth_lph");

		sprintf(filename,"/afs/hephy.at/scratch/k/knuenz/CMSSW_4_2_4_patch2/src/UpsilonPol/macros/polFit/h_lth_lph/Data_TheGreatRun_10B_May10_TESTS_SigmaBI0p05.root");
		TestFile1 = new TFile(filename);
		TH2D* h_lph_vs_lth_Test5 = (TH2D*)TestFile1->Get("h_lth_lph");

		sprintf(filename,"/afs/hephy.at/scratch/k/knuenz/CMSSW_4_2_4_patch2/src/UpsilonPol/macros/polFit/h_lth_lph/Data_TheGreatRun_10B_May10_TESTS_SigmaBI0p2.root");
		TestFile1 = new TFile(filename);
		TH2D* h_lph_vs_lth_Test6 = (TH2D*)TestFile1->Get("h_lth_lph");

		sprintf(filename,"/afs/hephy.at/scratch/k/knuenz/CMSSW_4_2_4_patch2/src/UpsilonPol/macros/polFit/h_lth_lph/h_lth_lph_Toy_May9_mal4AfterBurnIn.root");
		TestFile1 = new TFile(filename);
		TH2D* h_lph_vs_lth_Test7 = (TH2D*)TestFile1->Get("h_lth_lph");


		//	  TH2D* h_lph_vs_lth_Test   = new TH2D( "h_lph_vs_lth_Test", "", 50, -0.15,0.2,50,-.02,.075);
		TH2D* h_lph_vs_lth_Test   = new TH2D( "h_lph_vs_lth_Test", "", 50, -0.3,0.5,50,-.125,.15);

		h_lph_vs_lth_Test->Reset();
		h_lph_vs_lth_Test->SetName("h_lph_vs_lth");
		h_lph_vs_lth_Test->SetStats(0);
		h_lph_vs_lth_Test->SetTitle("");
		h_lph_vs_lth_Test->GetXaxis()->SetTitle("#lambda_{#vartheta}^{PX}");
		h_lph_vs_lth_Test->GetXaxis()->SetLabelOffset(0.028);
		h_lph_vs_lth_Test->GetXaxis()->SetTitleSize(0.05);
		h_lph_vs_lth_Test->GetXaxis()->SetTickLength(-0.03);
		h_lph_vs_lth_Test->GetXaxis()->SetTitleOffset(1.25);
		h_lph_vs_lth_Test->GetYaxis()->SetTitle("#lambda_{#varphi}^{PX}");
		h_lph_vs_lth_Test->GetYaxis()->SetLabelOffset(0.032);
		h_lph_vs_lth_Test->GetYaxis()->SetTitleSize(0.05);
		h_lph_vs_lth_Test->GetYaxis()->SetTickLength(-0.03);
		h_lph_vs_lth_Test->GetYaxis()->SetTitleOffset(1.75);
		h_lph_vs_lth_Test->Draw("");




		double IntegralInsideContour[100];
		// PX frame

		setContourHistogram ( h_lph_vs_lth_Test1 );
		h_lph_vs_lth_Test1->SetStats(0);
		h_lph_vs_lth_Test1->SetTitle("");
		h_lph_vs_lth_Test1->SetLineColor( kGreen+2 );
		h_lph_vs_lth_Test1->SetLineWidth( 2 );
		h_lph_vs_lth_Test1->SetLineStyle( 1 );
		h_lph_vs_lth_Test1->Draw( "cont2, same" );
		IntegralInsideContour[1]=FractionInsideContour;

		setContourHistogram ( h_lph_vs_lth_Test2 );
		h_lph_vs_lth_Test2->SetStats(0);
		h_lph_vs_lth_Test2->SetTitle("");
		h_lph_vs_lth_Test2->SetLineColor( kRed );
		h_lph_vs_lth_Test2->SetLineWidth( 2 );
		h_lph_vs_lth_Test2->SetLineStyle( 1 );
		h_lph_vs_lth_Test2->Draw( "cont2, same" );
		IntegralInsideContour[2]=FractionInsideContour;

		setContourHistogram ( h_lph_vs_lth_Test3 );
		h_lph_vs_lth_Test3->SetStats(0);
		h_lph_vs_lth_Test3->SetTitle("");
		h_lph_vs_lth_Test3->SetLineColor( kBlue );
		h_lph_vs_lth_Test3->SetLineWidth( 2 );
		h_lph_vs_lth_Test3->SetLineStyle( 1 );
		h_lph_vs_lth_Test3->Draw( "cont2, same" );
		IntegralInsideContour[3]=FractionInsideContour;

		setContourHistogram ( h_lph_vs_lth_Test4 );
		h_lph_vs_lth_Test4->SetStats(0);
		h_lph_vs_lth_Test4->SetTitle("");
		h_lph_vs_lth_Test4->SetLineColor( kMagenta );
		h_lph_vs_lth_Test4->SetLineWidth( 2 );
		h_lph_vs_lth_Test4->SetLineStyle( 1 );
		h_lph_vs_lth_Test4->Draw( "cont2, same" );
		IntegralInsideContour[4]=FractionInsideContour;

		setContourHistogram ( h_lph_vs_lth_Test5 );
		h_lph_vs_lth_Test5->SetStats(0);
		h_lph_vs_lth_Test5->SetTitle("");
		h_lph_vs_lth_Test5->SetLineColor( kOrange );
		h_lph_vs_lth_Test5->SetLineWidth( 2 );
		h_lph_vs_lth_Test5->SetLineStyle( 1 );
		h_lph_vs_lth_Test5->Draw( "cont2, same" );
		IntegralInsideContour[5]=FractionInsideContour;

		setContourHistogram ( h_lph_vs_lth_Test6 );
		h_lph_vs_lth_Test6->SetStats(0);
		h_lph_vs_lth_Test6->SetTitle("");
		h_lph_vs_lth_Test6->SetLineColor( kGreen-7 );
		h_lph_vs_lth_Test6->SetLineWidth( 2 );
		h_lph_vs_lth_Test6->SetLineStyle( 1 );
		h_lph_vs_lth_Test6->Draw( "cont2, same" );
		IntegralInsideContour[6]=FractionInsideContour;

		setContourHistogram ( h_lph_vs_lth_Test7 );
		h_lph_vs_lth_Test7->SetStats(0);
		h_lph_vs_lth_Test7->SetTitle("");
		h_lph_vs_lth_Test7->SetLineColor( kCyan );
		h_lph_vs_lth_Test7->SetLineWidth( 2 );
		h_lph_vs_lth_Test7->SetLineStyle( 1 );
		//	  h_lph_vs_lth_Test7->Draw( "cont2, same" );
		IntegralInsideContour[7]=FractionInsideContour;




		plotLegend=new TLegend(0.225,0.625,0.9,0.95);
		plotLegend->SetFillColor(kWhite);
		//		 plotLegend->SetTextFont(72);
		plotLegend->SetTextSize(0.03);
		plotLegend->SetBorderSize(1);
		sprintf(legendentry,"#lambda_{#theta}^{PX,start} = 0.0, #sigma^{BurnIn}=0.1,Int_{cont.}=%1.4f",IntegralInsideContour[1]);
		plotLegend->AddEntry(h_lph_vs_lth_Test1,legendentry,"l");
		sprintf(legendentry,"#lambda_{#theta}^{PX,start} = rand., #sigma^{BurnIn}=0.1,Int_{cont.}=%1.4f",IntegralInsideContour[2]);
		plotLegend->AddEntry(h_lph_vs_lth_Test2,legendentry,"l");
		sprintf(legendentry,"#lambda_{#theta}^{PX,start} = 0.5, #sigma^{BurnIn}=0.1,Int_{cont.}=%1.4f",IntegralInsideContour[3]);
		plotLegend->AddEntry(h_lph_vs_lth_Test3,legendentry,"l");
		sprintf(legendentry,"#lambda_{#theta}^{PX,start} = -0.5, #sigma^{BurnIn}=0.1,Int_{cont.}=%1.4f",IntegralInsideContour[4]);
		plotLegend->AddEntry(h_lph_vs_lth_Test4,legendentry,"l");
		sprintf(legendentry,"#lambda_{#theta}^{PX,start} = 0.0, #sigma^{BurnIn}=0.05,Int_{cont.}=%1.4f",IntegralInsideContour[5]);
		plotLegend->AddEntry(h_lph_vs_lth_Test5,legendentry,"l");
		sprintf(legendentry,"#lambda_{#theta}^{PX,start} = 0.0, #sigma^{BurnIn}=0.2,Int_{cont.}=%1.4f",IntegralInsideContour[6]);
		plotLegend->AddEntry(h_lph_vs_lth_Test6,legendentry,"l");
		/*		 sprintf(legendentry,"#lambda_{#theta}^{PX,start} = 0.5, #sigma^{BurnIn}=0.1, #sigma^{AfterBurnIn}=1.2*#sigma_{test}");
					 plotLegend->AddEntry(h_lph_vs_lth_Test7,legendentry,"l");
		 */
		plotLegend->Draw(); plotLegend->Draw();


		sprintf(filename,"%s/lph_vs_lth_%s_BurnInTest.pdf",dirstruct,TreeBinID);
		c1->Print( filename );










		sprintf(filename,"/afs/hephy.at/scratch/k/knuenz/CMSSW_4_2_4_patch2/src/UpsilonPol/macros/polFit/h_lth_lph/Data_TheGreatRun_10B_May10_TESTS_DEFAULT.root");
		TestFile1 = new TFile(filename);
		TH2D* h_lph_vs_lth_TestABI1 = (TH2D*)TestFile1->Get("h_lth_lph");

		sprintf(filename,"/afs/hephy.at/scratch/k/knuenz/CMSSW_4_2_4_patch2/src/UpsilonPol/macros/polFit/h_lth_lph/Data_TheGreatRun_10B_May10_TESTS_SigmaABI0p5.root");
		TestFile1 = new TFile(filename);
		TH2D* h_lph_vs_lth_TestABI2 = (TH2D*)TestFile1->Get("h_lth_lph");

		sprintf(filename,"/afs/hephy.at/scratch/k/knuenz/CMSSW_4_2_4_patch2/src/UpsilonPol/macros/polFit/h_lth_lph/Data_TheGreatRun_10B_May10_TESTS_SigmaABI1p0.root");
		TestFile1 = new TFile(filename);
		TH2D* h_lph_vs_lth_TestABI3 = (TH2D*)TestFile1->Get("h_lth_lph");

		sprintf(filename,"/afs/hephy.at/scratch/k/knuenz/CMSSW_4_2_4_patch2/src/UpsilonPol/macros/polFit/h_lth_lph/Data_TheGreatRun_10B_May10_TESTS_SigmaABI1p5.root");
		TestFile1 = new TFile(filename);
		TH2D* h_lph_vs_lth_TestABI4 = (TH2D*)TestFile1->Get("h_lth_lph");

		sprintf(filename,"/afs/hephy.at/scratch/k/knuenz/CMSSW_4_2_4_patch2/src/UpsilonPol/macros/polFit/h_lth_lph/Data_TheGreatRun_10B_May10_TESTS_SigmaABI2p0.root");
		TestFile1 = new TFile(filename);
		TH2D* h_lph_vs_lth_TestABI5 = (TH2D*)TestFile1->Get("h_lth_lph");

		sprintf(filename,"/afs/hephy.at/scratch/k/knuenz/CMSSW_4_2_4_patch2/src/UpsilonPol/macros/polFit/h_lth_lph/Data_TheGreatRun_10B_May10_TESTS_SigmaABI3p0.root");
		TestFile1 = new TFile(filename);
		TH2D* h_lph_vs_lth_TestABI6 = (TH2D*)TestFile1->Get("h_lth_lph");

		sprintf(filename,"/afs/hephy.at/scratch/k/knuenz/CMSSW_4_2_4_patch2/src/UpsilonPol/macros/polFit/h_lth_lph/h_lth_lph_Toy_May9_PXp05_new.root");
		TestFile1 = new TFile(filename);
		TH2D* h_lph_vs_lth_TestABI7 = (TH2D*)TestFile1->Get("h_lth_lph");


		TH2D* h_lph_vs_lth_TestABI   = new TH2D( "h_lph_vs_lth_TestABI", "", 50, -0.3,0.5,50,-.125,.15);

		h_lph_vs_lth_TestABI->Reset();
		h_lph_vs_lth_TestABI->SetName("h_lph_vs_lth");
		h_lph_vs_lth_TestABI->SetStats(0);
		h_lph_vs_lth_TestABI->SetTitle("");
		h_lph_vs_lth_TestABI->GetXaxis()->SetTitle("#lambda_{#vartheta}^{PX}");
		h_lph_vs_lth_TestABI->GetXaxis()->SetLabelOffset(0.028);
		h_lph_vs_lth_TestABI->GetXaxis()->SetTitleSize(0.05);
		h_lph_vs_lth_TestABI->GetXaxis()->SetTickLength(-0.03);
		h_lph_vs_lth_TestABI->GetXaxis()->SetTitleOffset(1.25);
		h_lph_vs_lth_TestABI->GetYaxis()->SetTitle("#lambda_{#varphi}^{PX}");
		h_lph_vs_lth_TestABI->GetYaxis()->SetLabelOffset(0.032);
		h_lph_vs_lth_TestABI->GetYaxis()->SetTitleSize(0.05);
		h_lph_vs_lth_TestABI->GetYaxis()->SetTickLength(-0.03);
		h_lph_vs_lth_TestABI->GetYaxis()->SetTitleOffset(1.75);
		h_lph_vs_lth_TestABI->Draw("");




		// PX frame

		setContourHistogram ( h_lph_vs_lth_TestABI1 );
		h_lph_vs_lth_TestABI1->SetStats(0);
		h_lph_vs_lth_TestABI1->SetTitle("");
		h_lph_vs_lth_TestABI1->SetLineColor( kGreen+2 );
		h_lph_vs_lth_TestABI1->SetLineWidth( 2 );
		h_lph_vs_lth_TestABI1->SetLineStyle( 1 );
		h_lph_vs_lth_TestABI1->Draw( "cont2, same" );
		IntegralInsideContour[1]=FractionInsideContour;

		setContourHistogram ( h_lph_vs_lth_TestABI2 );
		h_lph_vs_lth_TestABI2->SetStats(0);
		h_lph_vs_lth_TestABI2->SetTitle("");
		h_lph_vs_lth_TestABI2->SetLineColor( kRed );
		h_lph_vs_lth_TestABI2->SetLineWidth( 2 );
		h_lph_vs_lth_TestABI2->SetLineStyle( 1 );
		h_lph_vs_lth_TestABI2->Draw( "cont2, same" );
		IntegralInsideContour[2]=FractionInsideContour;

		setContourHistogram ( h_lph_vs_lth_TestABI3 );
		h_lph_vs_lth_TestABI3->SetStats(0);
		h_lph_vs_lth_TestABI3->SetTitle("");
		h_lph_vs_lth_TestABI3->SetLineColor( kBlue );
		h_lph_vs_lth_TestABI3->SetLineWidth( 2 );
		h_lph_vs_lth_TestABI3->SetLineStyle( 1 );
		h_lph_vs_lth_TestABI3->Draw( "cont2, same" );
		IntegralInsideContour[3]=FractionInsideContour;

		setContourHistogram ( h_lph_vs_lth_TestABI4 );
		h_lph_vs_lth_TestABI4->SetStats(0);
		h_lph_vs_lth_TestABI4->SetTitle("");
		h_lph_vs_lth_TestABI4->SetLineColor( kMagenta );
		h_lph_vs_lth_TestABI4->SetLineWidth( 2 );
		h_lph_vs_lth_TestABI4->SetLineStyle( 1 );
		h_lph_vs_lth_TestABI4->Draw( "cont2, same" );
		IntegralInsideContour[4]=FractionInsideContour;

		setContourHistogram ( h_lph_vs_lth_TestABI5 );
		h_lph_vs_lth_TestABI5->SetStats(0);
		h_lph_vs_lth_TestABI5->SetTitle("");
		h_lph_vs_lth_TestABI5->SetLineColor( kOrange );
		h_lph_vs_lth_TestABI5->SetLineWidth( 2 );
		h_lph_vs_lth_TestABI5->SetLineStyle( 1 );
		h_lph_vs_lth_TestABI5->Draw( "cont2, same" );
		IntegralInsideContour[5]=FractionInsideContour;

		setContourHistogram ( h_lph_vs_lth_TestABI6 );
		h_lph_vs_lth_TestABI6->SetStats(0);
		h_lph_vs_lth_TestABI6->SetTitle("");
		h_lph_vs_lth_TestABI6->SetLineColor( kGreen-7 );
		h_lph_vs_lth_TestABI6->SetLineWidth( 2 );
		h_lph_vs_lth_TestABI6->SetLineStyle( 1 );
		h_lph_vs_lth_TestABI6->Draw( "cont2, same" );
		IntegralInsideContour[6]=FractionInsideContour;

		setContourHistogram ( h_lph_vs_lth_TestABI7 );
		h_lph_vs_lth_TestABI7->SetStats(0);
		h_lph_vs_lth_TestABI7->SetTitle("");
		h_lph_vs_lth_TestABI7->SetLineColor( kCyan );
		h_lph_vs_lth_TestABI7->SetLineWidth( 2 );
		h_lph_vs_lth_TestABI7->SetLineStyle( 1 );
		//	  h_lph_vs_lth_TestABI7->Draw( "cont2, same" );
		IntegralInsideContour[7]=FractionInsideContour;




		plotLegend=new TLegend(0.225,0.625,0.9,0.95);
		plotLegend->SetFillColor(kWhite);
		//		 plotLegend->SetTextFont(72);
		plotLegend->SetTextSize(0.03);
		plotLegend->SetBorderSize(1);
		sprintf(legendentry,"#sigma^{AfterBurnIn}=0.3*#sigma_{test},Int_{cont.}=%1.4f",IntegralInsideContour[1]);
		plotLegend->AddEntry(h_lph_vs_lth_TestABI1,legendentry,"l");
		sprintf(legendentry,"#sigma^{AfterBurnIn}=0.5*#sigma_{test},Int_{cont.}=%1.4f",IntegralInsideContour[2]);
		plotLegend->AddEntry(h_lph_vs_lth_TestABI2,legendentry,"l");
		sprintf(legendentry,"#sigma^{AfterBurnIn}=1.0*#sigma_{test},Int_{cont.}=%1.4f",IntegralInsideContour[3]);
		plotLegend->AddEntry(h_lph_vs_lth_TestABI3,legendentry,"l");
		sprintf(legendentry,"#sigma^{AfterBurnIn}=1.5*#sigma_{test},Int_{cont.}=%1.4f",IntegralInsideContour[4]);
		plotLegend->AddEntry(h_lph_vs_lth_TestABI4,legendentry,"l");
		sprintf(legendentry,"#sigma^{AfterBurnIn}=2.0*#sigma_{test},Int_{cont.}=%1.4f",IntegralInsideContour[5]);
		plotLegend->AddEntry(h_lph_vs_lth_TestABI5,legendentry,"l");
		sprintf(legendentry,"#sigma^{AfterBurnIn}=3.0*#sigma_{test},Int_{cont.}=%1.4f",IntegralInsideContour[6]);
		plotLegend->AddEntry(h_lph_vs_lth_TestABI6,legendentry,"l");
		//		 sprintf(legendentry,"#lambda_{#theta}^{PX,start} = 0.5, #sigma^{BurnIn}=0.1, #sigma^{AfterBurnIn}=0.3*#sigma_{test}, check");
		//		 plotLegend->AddEntry(h_lph_vs_lth_TestABI7,legendentry,"l");

		plotLegend->Draw(); plotLegend->Draw();


		sprintf(filename,"%s/lph_vs_lth_%s_AfterBurnInTest.pdf",dirstruct,TreeBinID);
		c1->Print( filename );


















		sprintf(filename,"/afs/hephy.at/scratch/k/knuenz/CMSSW_4_2_4_patch2/src/UpsilonPol/macros/polFit/h_lth_lph/Toy_May9_BurnInConditionTests_PXp05_BI5000_NoBGtry.root");
		TestFile1 = new TFile(filename);
		TH2D* h_lph_vs_lth_TestSanity1 = (TH2D*)TestFile1->Get("h_lth_lph");

		sprintf(filename,"/afs/hephy.at/scratch/k/knuenz/CMSSW_4_2_4_patch2/src/UpsilonPol/macros/polFit/h_lth_lph/Toy_May9_BurnInConditionTests_PXp05_BI5000_NoBGtry2.root");
		TestFile1 = new TFile(filename);
		TH2D* h_lph_vs_lth_TestSanity2 = (TH2D*)TestFile1->Get("h_lth_lph");

		sprintf(filename,"/afs/hephy.at/scratch/k/knuenz/CMSSW_4_2_4_patch2/src/UpsilonPol/macros/polFit/h_lth_lph/Toy_May9_BurnInConditionTests_PXp05_BI5000_NoBGtry3.root");
		TestFile1 = new TFile(filename);
		TH2D* h_lph_vs_lth_TestSanity3 = (TH2D*)TestFile1->Get("h_lth_lph");

		sprintf(filename,"/afs/hephy.at/scratch/k/knuenz/CMSSW_4_2_4_patch2/src/UpsilonPol/macros/polFit/h_lth_lph/Toy_May9_BurnInConditionTests_PXp05_BI5000_NoBGtry4.root");
		TestFile1 = new TFile(filename);
		TH2D* h_lph_vs_lth_TestSanity4 = (TH2D*)TestFile1->Get("h_lth_lph");

		sprintf(filename,"/afs/hephy.at/scratch/k/knuenz/CMSSW_4_2_4_patch2/src/UpsilonPol/macros/polFit/h_lth_lph/Toy_May9_BurnInConditionTests_PXp05_BI5000_NoBGtry5.root");
		TestFile1 = new TFile(filename);
		TH2D* h_lph_vs_lth_TestSanity5 = (TH2D*)TestFile1->Get("h_lth_lph");

		sprintf(filename,"/afs/hephy.at/scratch/k/knuenz/CMSSW_4_2_4_patch2/src/UpsilonPol/macros/polFit/h_lth_lph/h_lth_lph_Toy_May9_2AfterBurnIn.root");
		TestFile1 = new TFile(filename);
		TH2D* h_lph_vs_lth_TestSanity6 = (TH2D*)TestFile1->Get("h_lth_lph");

		sprintf(filename,"/afs/hephy.at/scratch/k/knuenz/CMSSW_4_2_4_patch2/src/UpsilonPol/macros/polFit/h_lth_lph/h_lth_lph_Toy_May9_PXp05_new.root");
		TestFile1 = new TFile(filename);
		TH2D* h_lph_vs_lth_TestSanity7 = (TH2D*)TestFile1->Get("h_lth_lph");


		TH2D* h_lph_vs_lth_TestSanity   = new TH2D( "h_lph_vs_lth_TestSanity", "", 50, -0.15,0.2,50,-.02,.075);

		h_lph_vs_lth_TestSanity->Reset();
		h_lph_vs_lth_TestSanity->SetName("h_lph_vs_lth");
		h_lph_vs_lth_TestSanity->SetStats(0);
		h_lph_vs_lth_TestSanity->SetTitle("");
		h_lph_vs_lth_TestSanity->GetXaxis()->SetTitle("#lambda_{#vartheta}^{PX}");
		h_lph_vs_lth_TestSanity->GetXaxis()->SetLabelOffset(0.028);
		h_lph_vs_lth_TestSanity->GetXaxis()->SetTitleSize(0.05);
		h_lph_vs_lth_TestSanity->GetXaxis()->SetTickLength(-0.03);
		h_lph_vs_lth_TestSanity->GetXaxis()->SetTitleOffset(1.25);
		h_lph_vs_lth_TestSanity->GetYaxis()->SetTitle("#lambda_{#varphi}^{PX}");
		h_lph_vs_lth_TestSanity->GetYaxis()->SetLabelOffset(0.032);
		h_lph_vs_lth_TestSanity->GetYaxis()->SetTitleSize(0.05);
		h_lph_vs_lth_TestSanity->GetYaxis()->SetTickLength(-0.03);
		h_lph_vs_lth_TestSanity->GetYaxis()->SetTitleOffset(1.75);
		h_lph_vs_lth_TestSanity->Draw("");




		// PX frame

		setContourHistogram ( h_lph_vs_lth_TestSanity1 );
		h_lph_vs_lth_TestSanity1->SetStats(0);
		h_lph_vs_lth_TestSanity1->SetTitle("");
		h_lph_vs_lth_TestSanity1->SetLineColor( kGreen+2 );
		h_lph_vs_lth_TestSanity1->SetLineWidth( 2 );
		h_lph_vs_lth_TestSanity1->SetLineStyle( 1 );
		h_lph_vs_lth_TestSanity1->Draw( "cont2, same" );

		setContourHistogram ( h_lph_vs_lth_TestSanity2 );
		h_lph_vs_lth_TestSanity2->SetStats(0);
		h_lph_vs_lth_TestSanity2->SetTitle("");
		h_lph_vs_lth_TestSanity2->SetLineColor( kRed );
		h_lph_vs_lth_TestSanity2->SetLineWidth( 2 );
		h_lph_vs_lth_TestSanity2->SetLineStyle( 1 );
		h_lph_vs_lth_TestSanity2->Draw( "cont2, same" );

		setContourHistogram ( h_lph_vs_lth_TestSanity3 );
		h_lph_vs_lth_TestSanity3->SetStats(0);
		h_lph_vs_lth_TestSanity3->SetTitle("");
		h_lph_vs_lth_TestSanity3->SetLineColor( kBlue );
		h_lph_vs_lth_TestSanity3->SetLineWidth( 2 );
		h_lph_vs_lth_TestSanity3->SetLineStyle( 1 );
		h_lph_vs_lth_TestSanity3->Draw( "cont2, same" );

		setContourHistogram ( h_lph_vs_lth_TestSanity4 );
		h_lph_vs_lth_TestSanity4->SetStats(0);
		h_lph_vs_lth_TestSanity4->SetTitle("");
		h_lph_vs_lth_TestSanity4->SetLineColor( kMagenta );
		h_lph_vs_lth_TestSanity4->SetLineWidth( 2 );
		h_lph_vs_lth_TestSanity4->SetLineStyle( 1 );
		h_lph_vs_lth_TestSanity4->Draw( "cont2, same" );

		setContourHistogram ( h_lph_vs_lth_TestSanity5 );
		h_lph_vs_lth_TestSanity5->SetStats(0);
		h_lph_vs_lth_TestSanity5->SetTitle("");
		h_lph_vs_lth_TestSanity5->SetLineColor( kOrange );
		h_lph_vs_lth_TestSanity5->SetLineWidth( 2 );
		h_lph_vs_lth_TestSanity5->SetLineStyle( 1 );
		h_lph_vs_lth_TestSanity5->Draw( "cont2, same" );

		setContourHistogram ( h_lph_vs_lth_TestSanity6 );
		h_lph_vs_lth_TestSanity6->SetStats(0);
		h_lph_vs_lth_TestSanity6->SetTitle("");
		h_lph_vs_lth_TestSanity6->SetLineColor( kGreen-7 );
		h_lph_vs_lth_TestSanity6->SetLineWidth( 2 );
		h_lph_vs_lth_TestSanity6->SetLineStyle( 1 );
		//	  h_lph_vs_lth_TestSanity6->Draw( "cont2, same" );

		setContourHistogram ( h_lph_vs_lth_TestSanity7 );
		h_lph_vs_lth_TestSanity7->SetStats(0);
		h_lph_vs_lth_TestSanity7->SetTitle("");
		h_lph_vs_lth_TestSanity7->SetLineColor( kCyan );
		h_lph_vs_lth_TestSanity7->SetLineWidth( 2 );
		h_lph_vs_lth_TestSanity7->SetLineStyle( 1 );
		//	  h_lph_vs_lth_TestSanity7->Draw( "cont2, same" );




		plotLegend=new TLegend(0.225,0.625,0.9,0.95);
		plotLegend->SetFillColor(kWhite);
		//		 plotLegend->SetTextFont(72);
		plotLegend->SetTextSize(0.03);
		plotLegend->SetBorderSize(1);
		sprintf(legendentry,"#lambda_{#theta}^{PX,start} = 0.5, #sigma^{BurnIn}=0.1, #sigma^{AfterBurnIn}=0.3*#sigma_{test},1");
		plotLegend->AddEntry(h_lph_vs_lth_TestSanity1,legendentry,"l");
		sprintf(legendentry,"#lambda_{#theta}^{PX,start} = 0.5, #sigma^{BurnIn}=0.1, #sigma^{AfterBurnIn}=0.3*#sigma_{test},2");
		plotLegend->AddEntry(h_lph_vs_lth_TestSanity2,legendentry,"l");
		sprintf(legendentry,"#lambda_{#theta}^{PX,start} = 0.5, #sigma^{BurnIn}=0.1, #sigma^{AfterBurnIn}=0.3*#sigma_{test},3");
		plotLegend->AddEntry(h_lph_vs_lth_TestSanity3,legendentry,"l");
		sprintf(legendentry,"#lambda_{#theta}^{PX,start} = 0.5, #sigma^{BurnIn}=0.1, #sigma^{AfterBurnIn}=0.3*#sigma_{test},4");
		plotLegend->AddEntry(h_lph_vs_lth_TestSanity4,legendentry,"l");
		sprintf(legendentry,"#lambda_{#theta}^{PX,start} = 0.5, #sigma^{BurnIn}=0.1, #sigma^{AfterBurnIn}=0.3*#sigma_{test},5");
		plotLegend->AddEntry(h_lph_vs_lth_TestSanity5,legendentry,"l");
		//		 sprintf(legendentry,"#lambda_{#theta}^{PX,start} = 0.5, #sigma^{BurnIn}=0.1, #sigma^{AfterBurnIn}=0.3*#sigma_{test}");
		//		 plotLegend->AddEntry(h_lph_vs_lth_TestSanity6,legendentry,"l");
		//		 sprintf(legendentry,"#lambda_{#theta}^{PX,start} = 0.5, #sigma^{BurnIn}=0.1, #sigma^{AfterBurnIn}=0.3*#sigma_{test}, check");
		//		 plotLegend->AddEntry(h_lph_vs_lth_TestSanity7,legendentry,"l");

		plotLegend->Draw(); plotLegend->Draw();


		sprintf(filename,"%s/lph_vs_lth_%s_SanityCheck.pdf",dirstruct,TreeBinID);
		c1->Print( filename );


















		sprintf(filename,"/afs/hephy.at/scratch/k/knuenz/CMSSW_4_2_4_patch2/src/UpsilonPol/macros/polFit/h_lth_lph/h_lth_lph_Data_3S_rap1_pT9_May9_MergedUpToFit50.root");
		TestFile1 = new TFile(filename);
		TH2D* h_lph_vs_lth_Test1_nFit = (TH2D*)TestFile1->Get("h_lth_lph");

		sprintf(filename,"/afs/hephy.at/scratch/k/knuenz/CMSSW_4_2_4_patch2/src/UpsilonPol/macros/polFit/h_lth_lph/h_lth_lph_Data_3S_rap1_pT9_May9_MergedUpToFit40.root");
		TestFile1 = new TFile(filename);
		TH2D* h_lph_vs_lth_Test2_nFit = (TH2D*)TestFile1->Get("h_lth_lph");

		sprintf(filename,"/afs/hephy.at/scratch/k/knuenz/CMSSW_4_2_4_patch2/src/UpsilonPol/macros/polFit/h_lth_lph/h_lth_lph_Data_3S_rap1_pT9_May9_MergedUpToFit30.root");
		TestFile1 = new TFile(filename);
		TH2D* h_lph_vs_lth_Test3_nFit = (TH2D*)TestFile1->Get("h_lth_lph");

		sprintf(filename,"/afs/hephy.at/scratch/k/knuenz/CMSSW_4_2_4_patch2/src/UpsilonPol/macros/polFit/h_lth_lph/h_lth_lph_Data_3S_rap1_pT9_May9_MergedUpToFit20.root");
		TestFile1 = new TFile(filename);
		TH2D* h_lph_vs_lth_Test4_nFit = (TH2D*)TestFile1->Get("h_lth_lph");

		sprintf(filename,"/afs/hephy.at/scratch/k/knuenz/CMSSW_4_2_4_patch2/src/UpsilonPol/macros/polFit/h_lth_lph/h_lth_lph_Data_3S_rap1_pT9_May9_MergedUpToFit10.root");
		TestFile1 = new TFile(filename);
		TH2D* h_lph_vs_lth_Test5_nFit = (TH2D*)TestFile1->Get("h_lth_lph");

		sprintf(filename,"/afs/hephy.at/scratch/k/knuenz/CMSSW_4_2_4_patch2/src/UpsilonPol/macros/polFit/h_lth_lph/h_lth_lph_Data_3S_rap1_pT9_May9_MergedUpToFit5.root");
		TestFile1 = new TFile(filename);
		TH2D* h_lph_vs_lth_Test6_nFit = (TH2D*)TestFile1->Get("h_lth_lph");

		sprintf(filename,"/afs/hephy.at/scratch/k/knuenz/CMSSW_4_2_4_patch2/src/UpsilonPol/macros/polFit/h_lth_lph/h_lth_lph_Data_3S_rap1_pT9_May9_MergedUpToFit1.root");
		TestFile1 = new TFile(filename);
		TH2D* h_lph_vs_lth_Test7_nFit = (TH2D*)TestFile1->Get("h_lth_lph");


		TH2D* h_lph_vs_lth_Test_nFit   = new TH2D( "h_lph_vs_lth_Test_nFit", "", 50, -0.3,0.5,50,-.125,.15);

		h_lph_vs_lth_Test_nFit->Reset();
		h_lph_vs_lth_Test_nFit->SetName("h_lph_vs_lth");
		h_lph_vs_lth_Test_nFit->SetStats(0);
		h_lph_vs_lth_Test_nFit->SetTitle("");
		h_lph_vs_lth_Test_nFit->GetXaxis()->SetTitle("#lambda_{#vartheta}^{PX}");
		h_lph_vs_lth_Test_nFit->GetXaxis()->SetLabelOffset(0.028);
		h_lph_vs_lth_Test_nFit->GetXaxis()->SetTitleSize(0.05);
		h_lph_vs_lth_Test_nFit->GetXaxis()->SetTickLength(-0.03);
		h_lph_vs_lth_Test_nFit->GetXaxis()->SetTitleOffset(1.25);
		h_lph_vs_lth_Test_nFit->GetYaxis()->SetTitle("#lambda_{#varphi}^{PX}");
		h_lph_vs_lth_Test_nFit->GetYaxis()->SetLabelOffset(0.032);
		h_lph_vs_lth_Test_nFit->GetYaxis()->SetTitleSize(0.05);
		h_lph_vs_lth_Test_nFit->GetYaxis()->SetTickLength(-0.03);
		h_lph_vs_lth_Test_nFit->GetYaxis()->SetTitleOffset(1.75);
		h_lph_vs_lth_Test_nFit->Draw("");




		// PX frame

		setContourHistogram ( h_lph_vs_lth_Test1_nFit );
		h_lph_vs_lth_Test1_nFit->SetStats(0);
		h_lph_vs_lth_Test1_nFit->SetTitle("");
		h_lph_vs_lth_Test1_nFit->SetLineColor( kGreen+2 );
		h_lph_vs_lth_Test1_nFit->SetLineWidth( 2 );
		h_lph_vs_lth_Test1_nFit->SetLineStyle( 1 );
		h_lph_vs_lth_Test1_nFit->Draw( "cont2, same" );

		setContourHistogram ( h_lph_vs_lth_Test2_nFit );
		h_lph_vs_lth_Test2_nFit->SetStats(0);
		h_lph_vs_lth_Test2_nFit->SetTitle("");
		h_lph_vs_lth_Test2_nFit->SetLineColor( kRed );
		h_lph_vs_lth_Test2_nFit->SetLineWidth( 2 );
		h_lph_vs_lth_Test2_nFit->SetLineStyle( 1 );
		h_lph_vs_lth_Test2_nFit->Draw( "cont2, same" );

		setContourHistogram ( h_lph_vs_lth_Test3_nFit );
		h_lph_vs_lth_Test3_nFit->SetStats(0);
		h_lph_vs_lth_Test3_nFit->SetTitle("");
		h_lph_vs_lth_Test3_nFit->SetLineColor( kBlue );
		h_lph_vs_lth_Test3_nFit->SetLineWidth( 2 );
		h_lph_vs_lth_Test3_nFit->SetLineStyle( 1 );
		h_lph_vs_lth_Test3_nFit->Draw( "cont2, same" );

		setContourHistogram ( h_lph_vs_lth_Test4_nFit );
		h_lph_vs_lth_Test4_nFit->SetStats(0);
		h_lph_vs_lth_Test4_nFit->SetTitle("");
		h_lph_vs_lth_Test4_nFit->SetLineColor( kMagenta );
		h_lph_vs_lth_Test4_nFit->SetLineWidth( 2 );
		h_lph_vs_lth_Test4_nFit->SetLineStyle( 1 );
		h_lph_vs_lth_Test4_nFit->Draw( "cont2, same" );

		setContourHistogram ( h_lph_vs_lth_Test5_nFit );
		h_lph_vs_lth_Test5_nFit->SetStats(0);
		h_lph_vs_lth_Test5_nFit->SetTitle("");
		h_lph_vs_lth_Test5_nFit->SetLineColor( kOrange );
		h_lph_vs_lth_Test5_nFit->SetLineWidth( 2 );
		h_lph_vs_lth_Test5_nFit->SetLineStyle( 1 );
		h_lph_vs_lth_Test5_nFit->Draw( "cont2, same" );

		setContourHistogram ( h_lph_vs_lth_Test6_nFit );
		h_lph_vs_lth_Test6_nFit->SetStats(0);
		h_lph_vs_lth_Test6_nFit->SetTitle("");
		h_lph_vs_lth_Test6_nFit->SetLineColor( kGreen-7 );
		h_lph_vs_lth_Test6_nFit->SetLineWidth( 2 );
		h_lph_vs_lth_Test6_nFit->SetLineStyle( 1 );
		h_lph_vs_lth_Test6_nFit->Draw( "cont2, same" );

		setContourHistogram ( h_lph_vs_lth_Test7_nFit );
		h_lph_vs_lth_Test7_nFit->SetStats(0);
		h_lph_vs_lth_Test7_nFit->SetTitle("");
		h_lph_vs_lth_Test7_nFit->SetLineColor( kCyan );
		h_lph_vs_lth_Test7_nFit->SetLineWidth( 2 );
		h_lph_vs_lth_Test7_nFit->SetLineStyle( 1 );
		h_lph_vs_lth_Test7_nFit->Draw( "cont2, same" );

		plotLegend=new TLegend(0.225,0.7,0.45,0.95);
		plotLegend->SetFillColor(kWhite);
		//		 plotLegend->SetTextFont(72);
		plotLegend->SetTextSize(0.03);
		plotLegend->SetBorderSize(1);
		sprintf(legendentry,"n_{fit} = 50");
		plotLegend->AddEntry(h_lph_vs_lth_Test1_nFit,legendentry,"l");
		sprintf(legendentry,"n_{fit} = 40");
		plotLegend->AddEntry(h_lph_vs_lth_Test2_nFit,legendentry,"l");
		sprintf(legendentry,"n_{fit} = 30");
		plotLegend->AddEntry(h_lph_vs_lth_Test3_nFit,legendentry,"l");
		sprintf(legendentry,"n_{fit} = 20");
		plotLegend->AddEntry(h_lph_vs_lth_Test4_nFit,legendentry,"l");
		sprintf(legendentry,"n_{fit} = 10");
		plotLegend->AddEntry(h_lph_vs_lth_Test5_nFit,legendentry,"l");
		sprintf(legendentry,"n_{fit} = 5");
		plotLegend->AddEntry(h_lph_vs_lth_Test6_nFit,legendentry,"l");
		sprintf(legendentry,"n_{fit} = 1");
		plotLegend->AddEntry(h_lph_vs_lth_Test7_nFit,legendentry,"l");

		plotLegend->Draw(); plotLegend->Draw();


		sprintf(filename,"%s/lph_vs_lth_%s_nFitTest.pdf",dirstruct,TreeBinID);
		c1->Print( filename );


	}













	///////////////////////////////////////////////////////////////////////////////
	// plot PDFs for lambdatilde in three frames

	c1 = new TCanvas("c1", "c1", 10, 28, 588,563);
	c1->Range(-1.370833,-114.012,1.029167,745.8285);
	c1->SetFillColor(0);
	c1->SetBorderMode(0);
	c1->SetBorderSize(0);
	c1->SetLeftMargin(0.1545139);
	c1->SetRightMargin(0.01215278);
	c1->SetTopMargin(0.01841621);
	c1->SetBottomMargin(0.1325967);
	c1->SetFrameBorderMode(0);

	TH1D* h_ltilde = (TH1D*)h_ltilde_CS->Clone();
	h_ltilde->SetName("h_ltilde");
	h_ltilde->Reset();

	h_ltilde_CS->Scale( 1./h_ltilde_CS->Integral() );
	h_ltilde_HX->Scale( 1./h_ltilde_HX->Integral() );
	h_ltilde_PX->Scale( 1./h_ltilde_PX->Integral() );

	double plotborder = 1.3; // how much space to leave beyond the maximum of the plots

	double plotMax = plotborder * h_ltilde_CS->GetMaximum();
	if ( plotborder * h_ltilde_HX->GetMaximum() > plotMax ) plotMax = plotborder * h_ltilde_HX->GetMaximum();
	if ( plotborder * h_ltilde_PX->GetMaximum() > plotMax ) plotMax = plotborder * h_ltilde_PX->GetMaximum();

	h_ltilde->GetXaxis()->SetTitle("#tilde{#lambda}");
	h_ltilde->GetXaxis()->SetLabelOffset(0.028);
	h_ltilde->GetXaxis()->SetTitleSize(0.05);
	h_ltilde->GetXaxis()->SetTickLength(-0.03);
	h_ltilde->GetXaxis()->SetTitleOffset(1.20);
	h_ltilde->GetYaxis()->SetTitle("parameter PDF [a.u.]");
	h_ltilde->GetYaxis()->SetLabelOffset(0.032);
	h_ltilde->GetYaxis()->SetTitleSize(0.05);
	h_ltilde->GetYaxis()->SetTickLength(-0.03);
	h_ltilde->GetYaxis()->SetTitleOffset(1.55);
	h_ltilde->SetMinimum(0.);
	h_ltilde->SetMaximum(plotMax);
	h_ltilde->Draw("");

	// CS frame

	h_ltilde_CS->SetLineColor( kBlue );
	h_ltilde_CS->SetLineWidth( 2 );
	h_ltilde_CS->SetLineStyle( 1 );
	h_ltilde_CS->Draw( "L same" );

	// HX frame

	h_ltilde_HX->SetLineColor( kRed );
	h_ltilde_HX->SetLineWidth( 2 );
	h_ltilde_HX->SetLineStyle( 1 );
	h_ltilde_HX->Draw( "L same" );

	// PX frame

	h_ltilde_PX->SetLineColor( kGreen+2 );
	h_ltilde_PX->SetLineWidth( 2 );
	h_ltilde_PX->SetLineStyle( 1 );
	h_ltilde_PX->Draw( "L same" );

	plotLegend=new TLegend(0.175,0.81,0.275,0.96);
	plotLegend->SetFillColor(kWhite);
	//  plotLegend->SetTextFont(72);
	plotLegend->SetTextSize(0.035);
	plotLegend->SetBorderSize(1);
	sprintf(legendentry,"CS");
	plotLegend->AddEntry(h_lph_vs_lth_CS,legendentry,"l");
	sprintf(legendentry,"HX");
	plotLegend->AddEntry(h_lph_vs_lth_HX,legendentry,"l");
	sprintf(legendentry,"PX");
	plotLegend->AddEntry(h_lph_vs_lth_PX,legendentry,"l");
	plotLegend->Draw(); plotLegend->Draw();
	plotLegend->Draw();

	plotLegend->Draw();

	sprintf(filename,"%s/ltilde_%s.pdf",dirstruct,TreeBinID);
	c1->Print( filename );


	///////////////////////////////////////////////////////////////////////////
	// extract best-fit values from the 1D histograms (taking mean)

	/*  double lth_CS_best = h_lth_CS->GetMean();
			double lph_CS_best = h_lph_CS->GetMean();
			double ltp_CS_best = h_ltp_CS->GetMean();

			double lth_HX_best = h_lth_HX->GetMean();
			double lph_HX_best = h_lph_HX->GetMean();
			double ltp_HX_best = h_ltp_HX->GetMean();

			double lth_PX_best = h_lth_PX->GetMean();
			double lph_PX_best = h_lph_PX->GetMean();
			double ltp_PX_best = h_ltp_PX->GetMean();

			double ltilde_CS_best = h_ltilde_CS->GetMean();
			double ltilde_HX_best = h_ltilde_HX->GetMean();
			double ltilde_PX_best = h_ltilde_PX->GetMean();
	 */

	double MPVerrorLow, MPVerrorHigh;
	double lth_CS_best, lph_CS_best, ltp_CS_best, ltilde_CS_best, lth_HX_best, lph_HX_best, ltp_HX_best, ltilde_HX_best, lth_PX_best, lph_PX_best, ltp_PX_best, ltilde_PX_best;
	bool PlotPosteriorDist=true;
	int nSigma=1;

	cout<<"fitting lth_CS_best"<<endl;
	FindMPV(h_lth_CS , lth_CS_best , MPVerrorLow, MPVerrorHigh, MPValgo, nSigma);
	sprintf(filename,"%s/PosteriorDist_lth_CS_%s.pdf",dirstruct,TreeBinID);
	if(PlotPosteriorDist)PlotPosterior(999, 999, "#lambda^{CS}_{#theta}", filename, h_lth_CS, lth_CS_best, MPVerrorLow, MPVerrorHigh);
	cout<<"fitting lph_CS_best"<<endl;
	FindMPV(h_lph_CS , lph_CS_best , MPVerrorLow, MPVerrorHigh, MPValgo, nSigma);
	sprintf(filename,"%s/PosteriorDist_lph_CS_%s.pdf",dirstruct,TreeBinID);
	if(PlotPosteriorDist)PlotPosterior(999, 999, "#lambda^{CS}_{#phi}", filename, h_lph_CS, lph_CS_best, MPVerrorLow, MPVerrorHigh);
	cout<<"fitting ltp_CS_best"<<endl;
	FindMPV(h_ltp_CS , ltp_CS_best , MPVerrorLow, MPVerrorHigh, MPValgo, nSigma);
	sprintf(filename,"%s/PosteriorDist_ltp_CS_%s.pdf",dirstruct,TreeBinID);
	if(PlotPosteriorDist)PlotPosterior(999, 999, "#lambda^{CS}_{#theta #phi}", filename, h_ltp_CS, ltp_CS_best, MPVerrorLow, MPVerrorHigh);
	cout<<"fitting ltilde_CS_best"<<endl;
	FindMPV(h_ltilde_CS , ltilde_CS_best , MPVerrorLow, MPVerrorHigh, MPValgo, nSigma);
	sprintf(filename,"%s/PosteriorDist_ltilde_CS_%s.pdf",dirstruct,TreeBinID);
	if(PlotPosteriorDist)PlotPosterior(999, 999, "#tilde{#lambda}^{CS}", filename, h_ltilde_CS, ltilde_CS_best, MPVerrorLow, MPVerrorHigh);

	cout<<"fitting lth_HX_best"<<endl;
	FindMPV(h_lth_HX , lth_HX_best , MPVerrorLow, MPVerrorHigh, MPValgo, nSigma);
	sprintf(filename,"%s/PosteriorDist_lth_HX_%s.pdf",dirstruct,TreeBinID);
	if(PlotPosteriorDist)PlotPosterior(999, 999, "#lambda^{HX}_{#theta}", filename, h_lth_HX, lth_HX_best, MPVerrorLow, MPVerrorHigh);
	cout<<"fitting lph_HX_best"<<endl;
	FindMPV(h_lph_HX , lph_HX_best , MPVerrorLow, MPVerrorHigh, MPValgo, nSigma);
	sprintf(filename,"%s/PosteriorDist_lph_HX_%s.pdf",dirstruct,TreeBinID);
	if(PlotPosteriorDist)PlotPosterior(999, 999, "#lambda^{HX}_{#phi}", filename, h_lph_HX, lph_HX_best, MPVerrorLow, MPVerrorHigh);
	cout<<"fitting ltp_HX_best"<<endl;
	FindMPV(h_ltp_HX , ltp_HX_best , MPVerrorLow, MPVerrorHigh, MPValgo, nSigma);
	sprintf(filename,"%s/PosteriorDist_ltp_HX_%s.pdf",dirstruct,TreeBinID);
	if(PlotPosteriorDist)PlotPosterior(999, 999, "#lambda^{HX}_{#theta #phi}", filename, h_ltp_HX, ltp_HX_best, MPVerrorLow, MPVerrorHigh);
	cout<<"fitting ltilde_HX_best"<<endl;
	FindMPV(h_ltilde_HX , ltilde_HX_best , MPVerrorLow, MPVerrorHigh, MPValgo, nSigma);
	sprintf(filename,"%s/PosteriorDist_ltilde_HX_%s.pdf",dirstruct,TreeBinID);
	if(PlotPosteriorDist)PlotPosterior(999, 999, "#tilde{#lambda}^{HX}", filename, h_ltilde_HX, ltilde_HX_best, MPVerrorLow, MPVerrorHigh);

	cout<<"fitting lth_PX_best"<<endl;
	FindMPV(h_lth_PX , lth_PX_best , MPVerrorLow, MPVerrorHigh, MPValgo, nSigma);
	sprintf(filename,"%s/PosteriorDist_lth_PX_%s.pdf",dirstruct,TreeBinID);
	if(PlotPosteriorDist)PlotPosterior(999, 999, "#lambda^{PX}_{#theta}", filename, h_lth_PX, lth_PX_best, MPVerrorLow, MPVerrorHigh);
	cout<<"fitting lph_PX_best"<<endl;
	FindMPV(h_lph_PX , lph_PX_best , MPVerrorLow, MPVerrorHigh, MPValgo, nSigma);
	sprintf(filename,"%s/PosteriorDist_lph_PX_%s.pdf",dirstruct,TreeBinID);
	if(PlotPosteriorDist)PlotPosterior(999, 999, "#lambda^{PX}_{#phi}", filename, h_lph_PX, lph_PX_best, MPVerrorLow, MPVerrorHigh);
	cout<<"fitting ltp_PX_best"<<endl;
	FindMPV(h_ltp_PX , ltp_PX_best , MPVerrorLow, MPVerrorHigh, MPValgo, nSigma);
	sprintf(filename,"%s/PosteriorDist_ltp_PX_%s.pdf",dirstruct,TreeBinID);
	if(PlotPosteriorDist)PlotPosterior(999, 999, "#lambda^{PX}_{#theta #phi}", filename, h_ltp_PX, ltp_PX_best, MPVerrorLow, MPVerrorHigh);
	cout<<"fitting ltilde_PX_best"<<endl;
	FindMPV(h_ltilde_PX , ltilde_PX_best , MPVerrorLow, MPVerrorHigh, MPValgo, nSigma);
	sprintf(filename,"%s/PosteriorDist_ltilde_PX_%s.pdf",dirstruct,TreeBinID);
	if(PlotPosteriorDist)PlotPosterior(999, 999, "#tilde{#lambda}^{PX}", filename, h_ltilde_PX, ltilde_PX_best, MPVerrorLow, MPVerrorHigh);

	cout<<"lth_CS_best    "<<lth_CS_best<<endl;
	cout<<"lph_CS_best    "<<lph_CS_best<<endl;
	cout<<"ltp_CS_best    "<<ltp_CS_best<<endl;
	cout<<"ltilde_CS_best "<<ltilde_CS_best<<endl;

	cout<<"lth_HX_best    "<<lth_HX_best<<endl;
	cout<<"lph_HX_best    "<<lph_HX_best<<endl;
	cout<<"ltp_HX_best    "<<ltp_HX_best<<endl;
	cout<<"ltilde_HX_best "<<ltilde_HX_best<<endl;

	cout<<"lth_PX_best    "<<lth_PX_best<<endl;
	cout<<"lph_PX_best    "<<lph_PX_best<<endl;
	cout<<"ltp_PX_best    "<<ltp_PX_best<<endl;
	cout<<"ltilde_PX_best "<<ltilde_PX_best<<endl;


	double ltilde_best = (ltilde_CS_best + ltilde_HX_best + ltilde_PX_best)/3.;

	// input ntuple with anguar distribution of background-subtracted dilepton events

	TTree* angles = (TTree*)results->Get("angles");

	double costh_CS;  angles->SetBranchAddress("costh_CS",     &costh_CS );
	double phi_CS;    angles->SetBranchAddress("phi_CS",       &phi_CS   );
	double phith_CS;  angles->SetBranchAddress("phithCS",      &phith_CS );

	double costh_HX;  angles->SetBranchAddress("costh_HX",     &costh_HX );
	double phi_HX;    angles->SetBranchAddress("phi_HX",       &phi_HX   );
	double phith_HX;  angles->SetBranchAddress("phith_HX",     &phith_HX );

	double costh_PX;  angles->SetBranchAddress("costh_PX",     &costh_PX );
	double phi_PX;    angles->SetBranchAddress("phi_PX",       &phi_PX   );
	double phith_PX;  angles->SetBranchAddress("phith_PX",     &phith_PX );

	double cosalpha;  angles->SetBranchAddress("cosalpha",     &cosalpha );
	double epsilon;   angles->SetBranchAddress("epsilon",      &epsilon  );


	// input histograms with angular dependencies of the parameter-independent PDF terms


	TH1D*  PDF_1_vs_cth_CS  = (TH1D*)results->Get("PDF_1_vs_cth_CS");
	TH1D*  PDF_1_vs_ph_CS   = (TH1D*)results->Get("PDF_1_vs_ph_CS");
	TH1D*  PDF_1_vs_phth_CS = (TH1D*)results->Get("PDF_1_vs_phth_CS");

	TH1D*  PDF_cth2_vs_cth_CS  = (TH1D*)results->Get("PDF_cth2_vs_cth_CS");
	TH1D*  PDF_cth2_vs_ph_CS   = (TH1D*)results->Get("PDF_cth2_vs_ph_CS");
	TH1D*  PDF_cth2_vs_phth_CS = (TH1D*)results->Get("PDF_cth2_vs_phth_CS");

	TH1D*  PDF_sth2c2ph_vs_cth_CS  = (TH1D*)results->Get("PDF_sth2c2ph_vs_cth_CS");
	TH1D*  PDF_sth2c2ph_vs_ph_CS   = (TH1D*)results->Get("PDF_sth2c2ph_vs_ph_CS");
	TH1D*  PDF_sth2c2ph_vs_phth_CS = (TH1D*)results->Get("PDF_sth2c2ph_vs_phth_CS");

	TH1D*  PDF_s2thcph_vs_cth_CS  = (TH1D*)results->Get("PDF_s2thcph_vs_cth_CS");
	TH1D*  PDF_s2thcph_vs_ph_CS   = (TH1D*)results->Get("PDF_s2thcph_vs_ph_CS");
	TH1D*  PDF_s2thcph_vs_phth_CS = (TH1D*)results->Get("PDF_s2thcph_vs_phth_CS");


	TH1D*  PDF_1_vs_cth_HX   = (TH1D*)results->Get("PDF_1_vs_cth_HX");
	TH1D*  PDF_1_vs_ph_HX    = (TH1D*)results->Get("PDF_1_vs_ph_HX");
	TH1D*  PDF_1_vs_phth_HX  = (TH1D*)results->Get("PDF_1_vs_phth_HX");

	TH1D*  PDF_cth2_vs_cth_HX  = (TH1D*)results->Get("PDF_cth2_vs_cth_HX");
	TH1D*  PDF_cth2_vs_ph_HX   = (TH1D*)results->Get("PDF_cth2_vs_ph_HX");
	TH1D*  PDF_cth2_vs_phth_HX = (TH1D*)results->Get("PDF_cth2_vs_phth_HX");

	TH1D*  PDF_sth2c2ph_vs_cth_HX  = (TH1D*)results->Get("PDF_sth2c2ph_vs_cth_HX");
	TH1D*  PDF_sth2c2ph_vs_ph_HX   = (TH1D*)results->Get("PDF_sth2c2ph_vs_ph_HX");
	TH1D*  PDF_sth2c2ph_vs_phth_HX = (TH1D*)results->Get("PDF_sth2c2ph_vs_phth_HX");

	TH1D*  PDF_s2thcph_vs_cth_HX  = (TH1D*)results->Get("PDF_s2thcph_vs_cth_HX");
	TH1D*  PDF_s2thcph_vs_ph_HX   = (TH1D*)results->Get("PDF_s2thcph_vs_ph_HX");
	TH1D*  PDF_s2thcph_vs_phth_HX = (TH1D*)results->Get("PDF_s2thcph_vs_phth_HX");


	TH1D*  PDF_1_vs_cth_PX   = (TH1D*)results->Get("PDF_1_vs_cth_PX");
	TH1D*  PDF_1_vs_ph_PX    = (TH1D*)results->Get("PDF_1_vs_ph_PX");
	TH1D*  PDF_1_vs_phth_PX  = (TH1D*)results->Get("PDF_1_vs_phth_PX");

	TH1D*  PDF_cth2_vs_cth_PX  = (TH1D*)results->Get("PDF_cth2_vs_cth_PX");
	TH1D*  PDF_cth2_vs_ph_PX   = (TH1D*)results->Get("PDF_cth2_vs_ph_PX");
	TH1D*  PDF_cth2_vs_phth_PX = (TH1D*)results->Get("PDF_cth2_vs_phth_PX");

	TH1D*  PDF_sth2c2ph_vs_cth_PX  = (TH1D*)results->Get("PDF_sth2c2ph_vs_cth_PX");
	TH1D*  PDF_sth2c2ph_vs_ph_PX   = (TH1D*)results->Get("PDF_sth2c2ph_vs_ph_PX");
	TH1D*  PDF_sth2c2ph_vs_phth_PX = (TH1D*)results->Get("PDF_sth2c2ph_vs_phth_PX");

	TH1D*  PDF_s2thcph_vs_cth_PX  = (TH1D*)results->Get("PDF_s2thcph_vs_cth_PX");
	TH1D*  PDF_s2thcph_vs_ph_PX   = (TH1D*)results->Get("PDF_s2thcph_vs_ph_PX");
	TH1D*  PDF_s2thcph_vs_phth_PX = (TH1D*)results->Get("PDF_s2thcph_vs_phth_PX");


	TH1D*  PDF_1_vs_calpha        = (TH1D*)results->Get("PDF_1_vs_calpha");
	TH1D*  PDF_cth2_vs_calpha     = (TH1D*)results->Get("PDF_cth2_vs_calpha");
	TH1D*  PDF_sth2c2ph_vs_calpha = (TH1D*)results->Get("PDF_sth2c2ph_vs_calpha");
	TH1D*  PDF_s2thcph_vs_calpha  = (TH1D*)results->Get("PDF_s2thcph_vs_calpha");


	double scaleFits=1./double(nTotalFits);
	if(scalePlots){
		cout<<"scaling PDF's by "<<scaleFits<<endl;

		PDF_1_vs_cth_CS->Scale(scaleFits);
		PDF_1_vs_ph_CS->Scale(scaleFits);
		PDF_1_vs_phth_CS->Scale(scaleFits);
		PDF_cth2_vs_cth_CS->Scale(scaleFits);
		PDF_cth2_vs_ph_CS->Scale(scaleFits);
		PDF_cth2_vs_phth_CS->Scale(scaleFits);
		PDF_sth2c2ph_vs_cth_CS->Scale(scaleFits);
		PDF_sth2c2ph_vs_ph_CS->Scale(scaleFits);
		PDF_sth2c2ph_vs_phth_CS->Scale(scaleFits);
		PDF_s2thcph_vs_cth_CS->Scale(scaleFits);
		PDF_s2thcph_vs_ph_CS->Scale(scaleFits);
		PDF_s2thcph_vs_phth_CS->Scale(scaleFits);

		PDF_1_vs_cth_HX->Scale(scaleFits);
		PDF_1_vs_ph_HX->Scale(scaleFits);
		PDF_1_vs_phth_HX->Scale(scaleFits);
		PDF_cth2_vs_cth_HX->Scale(scaleFits);
		PDF_cth2_vs_ph_HX->Scale(scaleFits);
		PDF_cth2_vs_phth_HX->Scale(scaleFits);
		PDF_sth2c2ph_vs_cth_HX->Scale(scaleFits);
		PDF_sth2c2ph_vs_ph_HX->Scale(scaleFits);
		PDF_sth2c2ph_vs_phth_HX->Scale(scaleFits);
		PDF_s2thcph_vs_cth_HX->Scale(scaleFits);
		PDF_s2thcph_vs_ph_HX->Scale(scaleFits);
		PDF_s2thcph_vs_phth_HX->Scale(scaleFits);

		PDF_1_vs_cth_PX->Scale(scaleFits);
		PDF_1_vs_ph_PX->Scale(scaleFits);
		PDF_1_vs_phth_PX->Scale(scaleFits);
		PDF_cth2_vs_cth_PX->Scale(scaleFits);
		PDF_cth2_vs_ph_PX->Scale(scaleFits);
		PDF_cth2_vs_phth_PX->Scale(scaleFits);
		PDF_sth2c2ph_vs_cth_PX->Scale(scaleFits);
		PDF_sth2c2ph_vs_ph_PX->Scale(scaleFits);
		PDF_sth2c2ph_vs_phth_PX->Scale(scaleFits);
		PDF_s2thcph_vs_cth_PX->Scale(scaleFits);
		PDF_s2thcph_vs_ph_PX->Scale(scaleFits);
		PDF_s2thcph_vs_phth_PX->Scale(scaleFits);

		PDF_1_vs_calpha->Scale(scaleFits);
		PDF_cth2_vs_calpha->Scale(scaleFits);
		PDF_sth2c2ph_vs_calpha->Scale(scaleFits);
		PDF_s2thcph_vs_calpha->Scale(scaleFits);

	}

	// create data histograms

	TH1D* DATA_cth_CS   = new TH1D( "DATA_cth_CS",  "", nbin_cth, -1., 1. );
	TH1D* DATA_ph_CS    = new TH1D( "DATA_ph_CS",   "", nbin_ph, -180., 180. );
	TH1D* DATA_phth_CS  = new TH1D( "DATA_phth_CS", "", nbin_ph, -180., 180. );

	TH1D* DATA_cth_HX   = new TH1D( "DATA_cth_HX",  "", nbin_cth, -1., 1. );
	TH1D* DATA_ph_HX    = new TH1D( "DATA_ph_HX",   "", nbin_ph, -180., 180. );
	TH1D* DATA_phth_HX  = new TH1D( "DATA_phth_HX", "", nbin_ph, -180., 180. );

	TH1D* DATA_cth_PX   = new TH1D( "DATA_cth_PX",  "", nbin_cth, -1., 1. );
	TH1D* DATA_ph_PX    = new TH1D( "DATA_ph_PX",   "", nbin_ph, -180., 180. );
	TH1D* DATA_phth_PX  = new TH1D( "DATA_phth_PX", "", nbin_ph, -180., 180. );

	TH1D* DATA_calpha   = new TH1D( "DATA_calpha",  "", nbin_cth, -1., 1. );


	// loop over entries in the ntuple of angles

	n_entries = int( angles->GetEntries() );

	cout << endl;
	cout << "Reading angular distributions (" << n_entries << ") entries"<< endl;
	cout << "-------------------------------------------------------------" << endl;
	cout << "Progress: ";

	n_step = n_entries/50;

	for ( int i_entry = 0; i_entry < n_entries; i_entry++ ) {

		if ( i_entry%n_step == 0 ) cout << "X";

		angles->GetEvent( i_entry );

		DATA_cth_CS->Fill( costh_CS );
		DATA_ph_CS->Fill( phi_CS );
		DATA_phth_CS->Fill( phith_CS );

		DATA_cth_HX->Fill( costh_HX );
		DATA_ph_HX->Fill( phi_HX );
		DATA_phth_HX->Fill( phith_HX );

		DATA_cth_PX->Fill( costh_PX );
		DATA_ph_PX->Fill( phi_PX );
		DATA_phth_PX->Fill( phith_PX );

		DATA_calpha->Fill( cosalpha );

	}
	// end of loop over entries in the ntuple of angles

	cout << endl << endl;

	if(scalePlots){
		cout<<"scaling Data by "<<scaleFits<<endl;
		DATA_cth_CS->Scale(scaleFits);
		DATA_ph_CS->Scale(scaleFits);
		DATA_phth_CS->Scale(scaleFits);

		DATA_cth_HX->Scale(scaleFits);
		DATA_ph_HX->Scale(scaleFits);
		DATA_phth_HX->Scale(scaleFits);

		DATA_cth_PX->Scale(scaleFits);
		DATA_ph_PX->Scale(scaleFits);
		DATA_phth_PX->Scale(scaleFits);

		DATA_calpha->Scale(scaleFits);
	}


	double lth_phiplots=0.0; //0.5;
	double lth_phithplots=0.5;
	//double lthplot_zeroLine=0.5;
	double lthplot_zeroLine=0.; // for Fiducial-test

	//bool plotZeros=false;
	bool plotZeros=true; // for Fiducial-test

	//plotborder = 1.5;
	plotborder = 1.7; // for Fiducial-test

	// plot angular PDFs (calculated combining parameter-independent terms with best-fit parameters)
	// normalized to data distribution

	TCanvas* c3 = new TCanvas("c3", "c3", 10, 28, 588,563);
	c3->Range(-1.370833,-114.012,1.029167,745.8285);
	c3->SetFillColor(0);
	c3->SetBorderMode(0);
	c3->SetBorderSize(0);
	c3->SetLeftMargin(0.1545139);
	c3->SetRightMargin(0.01215278);
	c3->SetTopMargin(0.01841621);
	c3->SetBottomMargin(0.1325967);
	c3->SetFrameBorderMode(0);


	// costh_CS

	TH1D* PDF_cth_CS = (TH1D*)PDF_1_vs_cth_CS->Clone();
	PDF_cth_CS->SetName("PDF_cth_CS");
	PDF_cth_CS->Add( PDF_1_vs_cth_CS, PDF_cth2_vs_cth_CS, 1., lth_CS_best );
	PDF_cth_CS->Add( PDF_sth2c2ph_vs_cth_CS, lph_CS_best );
	PDF_cth_CS->Add( PDF_s2thcph_vs_cth_CS, ltp_CS_best );
	PDF_cth_CS->Scale( DATA_cth_CS->Integral() * PDF_cth_CS->GetNbinsX() / ( PDF_cth_CS->Integral()*DATA_cth_CS->GetNbinsX() ) );

	TH1D* PDF_cth_CS_p1 = (TH1D*)PDF_1_vs_cth_CS->Clone();
	PDF_cth_CS_p1->SetName("PDF_cth_CS_p1");
	PDF_cth_CS_p1->Add( PDF_1_vs_cth_CS, PDF_cth2_vs_cth_CS, 1., 1. );
	PDF_cth_CS_p1->Scale( DATA_cth_CS->Integral() * PDF_cth_CS_p1->GetNbinsX() / ( PDF_cth_CS_p1->Integral()*DATA_cth_CS->GetNbinsX() ) );

	TH1D* PDF_cth_CS_m1 = (TH1D*)PDF_1_vs_cth_CS->Clone();
	PDF_cth_CS_m1->SetName("PDF_cth_CS_m1");
	PDF_cth_CS_m1->Add( PDF_1_vs_cth_CS, PDF_cth2_vs_cth_CS, 1., -1. );
	PDF_cth_CS_m1->Scale( DATA_cth_CS->Integral() * PDF_cth_CS_m1->GetNbinsX() / ( PDF_cth_CS_m1->Integral()*DATA_cth_CS->GetNbinsX() ) );

	TH1D* PDF_cth_CS_0 = (TH1D*)PDF_1_vs_cth_CS->Clone();
	PDF_cth_CS_0->SetName("PDF_cth_CS_0");
	PDF_cth_CS_0->Add( PDF_1_vs_cth_CS, PDF_cth2_vs_cth_CS, 1., lthplot_zeroLine );
	PDF_cth_CS_0->Scale( DATA_cth_CS->Integral() * PDF_cth_CS_0->GetNbinsX() / ( PDF_cth_CS_0->Integral()*DATA_cth_CS->GetNbinsX() ) );



	plotMax = plotborder * PDF_cth_CS->GetMaximum();
	if ( plotborder * PDF_cth_CS_p1->GetMaximum() > plotMax ) plotMax = plotborder * PDF_cth_CS_p1->GetMaximum();
	if ( plotborder * PDF_cth_CS_m1->GetMaximum() > plotMax ) plotMax = plotborder * PDF_cth_CS_m1->GetMaximum();

	PDF_cth_CS->GetXaxis()->SetTitle("cos#vartheta_{CS}");
	PDF_cth_CS->GetXaxis()->SetLabelOffset(0.028);
	PDF_cth_CS->GetXaxis()->SetTitleSize(0.05);
	PDF_cth_CS->GetXaxis()->SetTickLength(-0.03);
	PDF_cth_CS->GetXaxis()->SetTitleOffset(1.20);
	PDF_cth_CS->GetYaxis()->SetTitle("event PDF [a.u.]");
	PDF_cth_CS->GetYaxis()->SetLabelOffset(0.032);
	PDF_cth_CS->GetYaxis()->SetTitleSize(0.05);
	PDF_cth_CS->GetYaxis()->SetTickLength(-0.03);
	PDF_cth_CS->GetYaxis()->SetTitleOffset(1.55);
	PDF_cth_CS->SetMinimum(0.);
	PDF_cth_CS->SetLineColor(kRed);
	PDF_cth_CS->SetLineStyle(1);
	PDF_cth_CS->SetLineWidth(2);
	PDF_cth_CS->SetMaximum(plotMax);
	PDF_cth_CS->Draw("L");

	PDF_cth_CS_p1->SetLineColor(kRed);
	PDF_cth_CS_p1->SetLineStyle(2);
	PDF_cth_CS_p1->SetLineWidth(2);
	PDF_cth_CS_p1->Draw("L same");

	PDF_cth_CS_m1->SetLineColor(kRed);
	PDF_cth_CS_m1->SetLineStyle(3);
	PDF_cth_CS_m1->SetLineWidth(2);
	PDF_cth_CS_m1->Draw("L same");

	PDF_cth_CS_0->SetLineColor(kBlue); //(kRed);
	PDF_cth_CS_0->SetLineStyle(10);
	PDF_cth_CS_0->SetLineWidth(2);
	if(plotZeros) PDF_cth_CS_0->Draw("L same");

	DATA_cth_CS->SetLineColor(kBlack);
	DATA_cth_CS->Draw("E same");

	TLegend* plotLegend2;

	plotLegend2=new TLegend(0.2,0.8,0.75,0.95);
	plotLegend2->SetFillColor(kWhite);
	//plotLegend2->SetTextFont(72);
	plotLegend2->SetTextSize(0.035);
	plotLegend2->SetBorderSize(1);

	sprintf(legendentry,"Best fit curve");
	plotLegend2->AddEntry(PDF_cth_CS,legendentry,"l");
	sprintf(legendentry,"#lambda_{#theta} = +1,  #lambda_{#phi} = 0,  #lambda_{#theta#phi} = 0");
	plotLegend2->AddEntry(PDF_cth_CS_p1,legendentry,"l");
	sprintf(legendentry,"#lambda_{#theta} = -1,  #lambda_{#phi} = 0,  #lambda_{#theta#phi} = 0");
	plotLegend2->AddEntry(PDF_cth_CS_m1,legendentry,"l");
	if(plotZeros){
		sprintf(legendentry,"#lambda_{#theta} = 0,  #lambda_{#phi} = 0,  #lambda_{#theta#phi} = 0");
		plotLegend2->AddEntry(PDF_cth_CS_0,legendentry,"l");
	}
	plotLegend2->Draw(); plotLegend2->Draw();

	sprintf(filename,"%s/fit_CS_costh_%s.pdf",dirstruct,TreeBinID);
	c3->Print( filename );


	// phi_CS

	TH1D* PDF_ph_CS = (TH1D*)PDF_1_vs_ph_CS->Clone();
	PDF_ph_CS->SetName("PDF_ph_CS");
	PDF_ph_CS->Add( PDF_1_vs_ph_CS, PDF_cth2_vs_ph_CS, 1., lth_CS_best );
	PDF_ph_CS->Add( PDF_sth2c2ph_vs_ph_CS, lph_CS_best );
	PDF_ph_CS->Add( PDF_s2thcph_vs_ph_CS, ltp_CS_best );
	PDF_ph_CS->Scale( DATA_ph_CS->Integral() * PDF_ph_CS->GetNbinsX() / ( PDF_ph_CS->Integral()*DATA_ph_CS->GetNbinsX() ) );

	TH1D* PDF_ph_CS_p1 = (TH1D*)PDF_1_vs_ph_CS->Clone();
	PDF_ph_CS_p1->SetName("PDF_ph_CS_p1");
	PDF_ph_CS_p1->Add( PDF_1_vs_ph_CS, PDF_cth2_vs_ph_CS, 1., lth_phiplots );
	PDF_ph_CS_p1->Add( PDF_sth2c2ph_vs_ph_CS, 1. );
	PDF_ph_CS_p1->Scale( DATA_ph_CS->Integral() * PDF_ph_CS_p1->GetNbinsX() / ( PDF_ph_CS_p1->Integral()*DATA_ph_CS->GetNbinsX() ) );

	TH1D* PDF_ph_CS_m1 = (TH1D*)PDF_1_vs_ph_CS->Clone();
	PDF_ph_CS_m1->SetName("PDF_ph_CS_m1");
	PDF_ph_CS_m1->Add( PDF_1_vs_ph_CS, PDF_cth2_vs_ph_CS, 1., lth_phiplots );
	PDF_ph_CS_m1->Add( PDF_sth2c2ph_vs_ph_CS, -1. );
	PDF_ph_CS_m1->Scale( DATA_ph_CS->Integral() * PDF_ph_CS_m1->GetNbinsX() / ( PDF_ph_CS_m1->Integral()*DATA_ph_CS->GetNbinsX() ) );

	TH1D* PDF_ph_CS_0 = (TH1D*)PDF_1_vs_ph_CS->Clone();
	PDF_ph_CS_0->SetName("PDF_ph_CS_0");
	PDF_ph_CS_0->Add( PDF_1_vs_ph_CS, PDF_cth2_vs_ph_CS, 1., lth_phiplots );
	PDF_ph_CS_0->Add( PDF_sth2c2ph_vs_ph_CS, 0. );
	PDF_ph_CS_0->Scale( DATA_ph_CS->Integral() * PDF_ph_CS_0->GetNbinsX() / ( PDF_ph_CS_0->Integral()*DATA_ph_CS->GetNbinsX() ) );



	plotMax = plotborder * PDF_ph_CS->GetMaximum();
	if ( plotborder * PDF_ph_CS_p1->GetMaximum() > plotMax ) plotMax = plotborder * PDF_ph_CS_p1->GetMaximum();
	if ( plotborder * PDF_ph_CS_m1->GetMaximum() > plotMax ) plotMax = plotborder * PDF_ph_CS_m1->GetMaximum();

	PDF_ph_CS->GetXaxis()->SetTitle("#varphi_{CS}");
	PDF_ph_CS->GetXaxis()->SetLabelOffset(0.028);
	PDF_ph_CS->GetXaxis()->SetTitleSize(0.05);
	PDF_ph_CS->GetXaxis()->SetTickLength(-0.03);
	PDF_ph_CS->GetXaxis()->SetTitleOffset(1.20);
	PDF_ph_CS->GetYaxis()->SetTitle("event PDF [a.u.]");
	PDF_ph_CS->GetYaxis()->SetLabelOffset(0.032);
	PDF_ph_CS->GetYaxis()->SetTitleSize(0.05);
	PDF_ph_CS->GetYaxis()->SetTickLength(-0.03);
	PDF_ph_CS->GetYaxis()->SetTitleOffset(1.55);
	PDF_ph_CS->SetMinimum(0.);
	PDF_ph_CS->SetLineColor(kRed);
	PDF_ph_CS->SetLineStyle(1);
	PDF_ph_CS->SetLineWidth(2);
	PDF_ph_CS->SetMaximum(plotMax);
	PDF_ph_CS->Draw("L");

	PDF_ph_CS_p1->SetLineColor(kRed);
	PDF_ph_CS_p1->SetLineStyle(2);
	PDF_ph_CS_p1->SetLineWidth(2);
	PDF_ph_CS_p1->Draw("L same");

	PDF_ph_CS_m1->SetLineColor(kRed);
	PDF_ph_CS_m1->SetLineStyle(3);
	PDF_ph_CS_m1->SetLineWidth(2);
	PDF_ph_CS_m1->Draw("L same");

	PDF_ph_CS_0->SetLineColor(kBlue); //(kRed);
	PDF_ph_CS_0->SetLineStyle(10);
	PDF_ph_CS_0->SetLineWidth(2);
	if(plotZeros) PDF_ph_CS_0->Draw("L same");

	DATA_ph_CS->SetLineColor(kBlack);
	DATA_ph_CS->Draw("E same");

	plotLegend2=new TLegend(0.2,0.8,0.75,0.95);
	plotLegend2->SetFillColor(kWhite);
	//plotLegend2->SetTextFont(72);
	plotLegend2->SetTextSize(0.035);
	plotLegend2->SetBorderSize(1);

	sprintf(legendentry,"Best fit curve");
	plotLegend2->AddEntry(PDF_ph_CS,legendentry,"l");
	//sprintf(legendentry,"#lambda_{#theta} = +0.5,  #lambda_{#phi} = +0.5,  #lambda_{#theta#phi} = 0");
	sprintf(legendentry,"#lambda_{#theta} = 0,  #lambda_{#phi} = +0.5,  #lambda_{#theta#phi} = 0");
	plotLegend2->AddEntry(PDF_ph_CS_p1,legendentry,"l");
	//sprintf(legendentry,"#lambda_{#theta} = +0.5,  #lambda_{#phi} = -0.5,  #lambda_{#theta#phi} = 0");
	sprintf(legendentry,"#lambda_{#theta} = 0,  #lambda_{#phi} = -0.5,  #lambda_{#theta#phi} = 0");
	plotLegend2->AddEntry(PDF_ph_CS_m1,legendentry,"l");
	if(plotZeros){
		sprintf(legendentry,"#lambda_{#theta} = 0,  #lambda_{#phi} = 0,  #lambda_{#theta#phi} = 0");
		plotLegend2->AddEntry(PDF_ph_CS_0,legendentry,"l");
	}
	plotLegend2->Draw(); plotLegend2->Draw();

	sprintf(filename,"%s/fit_CS_phi_%s.pdf",dirstruct,TreeBinID);
	c3->Print( filename );


	// phith_CS

	TH1D* PDF_phth_CS = (TH1D*)PDF_1_vs_phth_CS->Clone();
	PDF_phth_CS->SetName("PDF_phth_CS");
	PDF_phth_CS->Add( PDF_1_vs_phth_CS, PDF_cth2_vs_phth_CS, 1., lth_CS_best );
	PDF_phth_CS->Add( PDF_sth2c2ph_vs_phth_CS, lph_CS_best );
	PDF_phth_CS->Add( PDF_s2thcph_vs_phth_CS, ltp_CS_best );
	PDF_phth_CS->Scale( DATA_phth_CS->Integral() * PDF_phth_CS->GetNbinsX() / ( PDF_phth_CS->Integral()*DATA_phth_CS->GetNbinsX() ) );

	TH1D* PDF_phth_CS_p1 = (TH1D*)PDF_1_vs_phth_CS->Clone();
	PDF_phth_CS_p1->SetName("PDF_phth_CS_p1");
	PDF_phth_CS_p1->Add( PDF_1_vs_phth_CS, PDF_sth2c2ph_vs_phth_CS, 1., lth_phithplots );
	PDF_phth_CS_p1->Add( PDF_s2thcph_vs_phth_CS, sqrt(2.)/2. );
	PDF_phth_CS_p1->Scale( DATA_phth_CS->Integral() * PDF_phth_CS_p1->GetNbinsX() / ( PDF_phth_CS_p1->Integral()*DATA_phth_CS->GetNbinsX() ) );

	TH1D* PDF_phth_CS_m1 = (TH1D*)PDF_1_vs_phth_CS->Clone();
	PDF_phth_CS_m1->SetName("PDF_phth_CS_m1");
	PDF_phth_CS_m1->Add( PDF_1_vs_phth_CS, PDF_sth2c2ph_vs_phth_CS, 1., lth_phithplots );
	PDF_phth_CS_m1->Add( PDF_s2thcph_vs_phth_CS, -sqrt(2.)/2. );
	PDF_phth_CS_m1->Scale( DATA_phth_CS->Integral() * PDF_phth_CS_m1->GetNbinsX() / ( PDF_phth_CS_m1->Integral()*DATA_phth_CS->GetNbinsX() ) );

	TH1D* PDF_phth_CS_0 = (TH1D*)PDF_1_vs_phth_CS->Clone();
	PDF_phth_CS_0->SetName("PDF_phth_CS_0");
	PDF_phth_CS_0->Add( PDF_1_vs_phth_CS, PDF_sth2c2ph_vs_phth_CS, 1., lth_phithplots );
	PDF_phth_CS_0->Add( PDF_s2thcph_vs_phth_CS, 0. );
	PDF_phth_CS_0->Scale( DATA_phth_CS->Integral() * PDF_phth_CS_0->GetNbinsX() / ( PDF_phth_CS_0->Integral()*DATA_phth_CS->GetNbinsX() ) );



	plotMax = plotborder * PDF_phth_CS->GetMaximum();
	if ( plotborder * PDF_phth_CS_p1->GetMaximum() > plotMax ) plotMax = plotborder * PDF_phth_CS_p1->GetMaximum();
	if ( plotborder * PDF_phth_CS_m1->GetMaximum() > plotMax ) plotMax = plotborder * PDF_phth_CS_m1->GetMaximum();

	PDF_phth_CS->GetXaxis()->SetTitle("#tilde{#varphi}_{CS}");
	PDF_phth_CS->GetXaxis()->SetLabelOffset(0.028);
	PDF_phth_CS->GetXaxis()->SetTitleSize(0.05);
	PDF_phth_CS->GetXaxis()->SetTickLength(-0.03);
	PDF_phth_CS->GetXaxis()->SetTitleOffset(1.20);
	PDF_phth_CS->GetYaxis()->SetTitle("event PDF [a.u.]");
	PDF_phth_CS->GetYaxis()->SetLabelOffset(0.032);
	PDF_phth_CS->GetYaxis()->SetTitleSize(0.05);
	PDF_phth_CS->GetYaxis()->SetTickLength(-0.03);
	PDF_phth_CS->GetYaxis()->SetTitleOffset(1.55);
	PDF_phth_CS->SetMinimum(0.);
	PDF_phth_CS->SetLineColor(kRed);
	PDF_phth_CS->SetLineStyle(1);
	PDF_phth_CS->SetLineWidth(2);
	PDF_phth_CS->SetMaximum(plotMax);
	PDF_phth_CS->Draw("L");

	PDF_phth_CS_p1->SetLineColor(kRed);
	PDF_phth_CS_p1->SetLineStyle(2);
	PDF_phth_CS_p1->SetLineWidth(2);
	PDF_phth_CS_p1->Draw("L same");

	PDF_phth_CS_m1->SetLineColor(kRed);
	PDF_phth_CS_m1->SetLineStyle(3);
	PDF_phth_CS_p1->SetLineWidth(2);
	PDF_phth_CS_m1->Draw("L same");

	PDF_phth_CS_0->SetLineColor(kBlue); //(kRed);
	PDF_phth_CS_0->SetLineStyle(10);
	PDF_phth_CS_0->SetLineWidth(2);
	if(plotZeros) PDF_phth_CS_0->Draw("L same");

	DATA_phth_CS->SetLineColor(kBlack);
	DATA_phth_CS->Draw("E same");

	plotLegend2=new TLegend(0.2,0.8,0.75,0.95);
	plotLegend2->SetFillColor(kWhite);
	//plotLegend2->SetTextFont(72);
	plotLegend2->SetTextSize(0.035);
	plotLegend2->SetBorderSize(1);

	sprintf(legendentry,"Best fit curve");
	plotLegend2->AddEntry(PDF_phth_CS,legendentry,"l");
	sprintf(legendentry,"#lambda_{#theta} = +0.5,  #lambda_{#phi} = 0,  #lambda_{#theta#phi} = +#sqrt{2}/2");
	plotLegend2->AddEntry(PDF_phth_CS_p1,legendentry,"l");
	sprintf(legendentry,"#lambda_{#theta} = +0.5,  #lambda_{#phi} = 0,  #lambda_{#theta#phi} = -#sqrt{2}/2");
	plotLegend2->AddEntry(PDF_phth_CS_m1,legendentry,"l");
	if(plotZeros){
		sprintf(legendentry,"#lambda_{#theta} = 0,  #lambda_{#phi} = 0,  #lambda_{#theta#phi} = 0");
		plotLegend2->AddEntry(PDF_phth_CS_0,legendentry,"l");
	}
	plotLegend2->Draw(); plotLegend2->Draw();

	sprintf(filename,"%s/fit_CS_phith_%s.pdf",dirstruct,TreeBinID);
	c3->Print( filename );


	// costh_HX

	TH1D* PDF_cth_HX = (TH1D*)PDF_1_vs_cth_HX->Clone();
	PDF_cth_HX->SetName("PDF_cth_HX");
	PDF_cth_HX->Add( PDF_1_vs_cth_HX, PDF_cth2_vs_cth_HX, 1., lth_HX_best );
	PDF_cth_HX->Add( PDF_sth2c2ph_vs_cth_HX, lph_HX_best );
	PDF_cth_HX->Add( PDF_s2thcph_vs_cth_HX, ltp_HX_best );
	PDF_cth_HX->Scale( DATA_cth_HX->Integral() * PDF_cth_HX->GetNbinsX() / ( PDF_cth_HX->Integral()*DATA_cth_HX->GetNbinsX() ) );

	TH1D* PDF_cth_HX_p1 = (TH1D*)PDF_1_vs_cth_HX->Clone();
	PDF_cth_HX_p1->SetName("PDF_cth_HX_p1");
	PDF_cth_HX_p1->Add( PDF_1_vs_cth_HX, PDF_cth2_vs_cth_HX, 1., 1. );
	PDF_cth_HX_p1->Scale( DATA_cth_HX->Integral() * PDF_cth_HX_p1->GetNbinsX() / ( PDF_cth_HX_p1->Integral()*DATA_cth_HX->GetNbinsX() ) );

	TH1D* PDF_cth_HX_m1 = (TH1D*)PDF_1_vs_cth_HX->Clone();
	PDF_cth_HX_m1->SetName("PDF_cth_HX_m1");
	PDF_cth_HX_m1->Add( PDF_1_vs_cth_HX, PDF_cth2_vs_cth_HX, 1., -1. );
	PDF_cth_HX_m1->Scale( DATA_cth_HX->Integral() * PDF_cth_HX_m1->GetNbinsX() / ( PDF_cth_HX_m1->Integral()*DATA_cth_HX->GetNbinsX() ) );

	TH1D* PDF_cth_HX_0 = (TH1D*)PDF_1_vs_cth_HX->Clone();
	PDF_cth_HX_0->SetName("PDF_cth_HX_0");
	PDF_cth_HX_0->Add( PDF_1_vs_cth_HX, PDF_cth2_vs_cth_HX, 1., lthplot_zeroLine );
	PDF_cth_HX_0->Scale( DATA_cth_HX->Integral() * PDF_cth_HX_0->GetNbinsX() / ( PDF_cth_HX_0->Integral()*DATA_cth_HX->GetNbinsX() ) );

	TH1D* PDF_cth_HX_0p1 = (TH1D*)PDF_1_vs_cth_HX->Clone();
	PDF_cth_HX_0p1->SetName("PDF_cth_HX_0p1");
	PDF_cth_HX_0p1->Add( PDF_1_vs_cth_HX, PDF_cth2_vs_cth_HX, 1., 0.1 );
	PDF_cth_HX_0p1->Scale( DATA_cth_HX->Integral() * PDF_cth_HX_0p1->GetNbinsX() / ( PDF_cth_HX_0p1->Integral()*DATA_cth_HX->GetNbinsX() ) );


	plotMax = plotborder * PDF_cth_HX->GetMaximum();
	if ( plotborder * PDF_cth_HX_p1->GetMaximum() > plotMax ) plotMax = plotborder * PDF_cth_HX_p1->GetMaximum();
	if ( plotborder * PDF_cth_HX_m1->GetMaximum() > plotMax ) plotMax = plotborder * PDF_cth_HX_m1->GetMaximum();

	PDF_cth_HX->GetXaxis()->SetTitle("cos#vartheta_{HX}");
	PDF_cth_HX->GetXaxis()->SetLabelOffset(0.028);
	PDF_cth_HX->GetXaxis()->SetTitleSize(0.05);
	PDF_cth_HX->GetXaxis()->SetTickLength(-0.03);
	PDF_cth_HX->GetXaxis()->SetTitleOffset(1.20);
	PDF_cth_HX->GetYaxis()->SetTitle("event PDF [a.u.]");
	PDF_cth_HX->GetYaxis()->SetLabelOffset(0.032);
	PDF_cth_HX->GetYaxis()->SetTitleSize(0.05);
	PDF_cth_HX->GetYaxis()->SetTickLength(-0.03);
	PDF_cth_HX->GetYaxis()->SetTitleOffset(1.55);
	PDF_cth_HX->SetMinimum(0.);
	PDF_cth_HX->SetLineColor(kRed);
	PDF_cth_HX->SetLineStyle(1);
	PDF_cth_HX->SetLineWidth(2);
	PDF_cth_HX->SetMaximum(plotMax);
	PDF_cth_HX->Draw("L");

	PDF_cth_HX_p1->SetLineColor(kRed);
	PDF_cth_HX_p1->SetLineStyle(2);
	PDF_cth_HX_p1->SetLineWidth(2);
	PDF_cth_HX_p1->Draw("L same");

	PDF_cth_HX_m1->SetLineColor(kRed);
	PDF_cth_HX_m1->SetLineStyle(3);
	PDF_cth_HX_m1->SetLineWidth(2);
	PDF_cth_HX_m1->Draw("L same");

	PDF_cth_HX_0->SetLineColor(kBlue); //(kRed);
	PDF_cth_HX_0->SetLineStyle(10);
	PDF_cth_HX_0->SetLineWidth(2);
	if(plotZeros) PDF_cth_HX_0->Draw("L same");

	DATA_cth_HX->SetLineColor(kBlack);
	DATA_cth_HX->Draw("E same");

	plotLegend2=new TLegend(0.2,0.8,0.75,0.95);
	plotLegend2->SetFillColor(kWhite);
	//plotLegend2->SetTextFont(72);
	plotLegend2->SetTextSize(0.035);
	plotLegend2->SetBorderSize(1);

	sprintf(legendentry,"Best fit curve");
	plotLegend2->AddEntry(PDF_cth_HX,legendentry,"l");
	sprintf(legendentry,"#lambda_{#theta} = +1,  #lambda_{#phi} = 0,  #lambda_{#theta#phi} = 0");
	plotLegend2->AddEntry(PDF_cth_HX_p1,legendentry,"l");
	sprintf(legendentry,"#lambda_{#theta} = -1,  #lambda_{#phi} = 0,  #lambda_{#theta#phi} = 0");
	plotLegend2->AddEntry(PDF_cth_HX_m1,legendentry,"l");
	if(plotZeros){
		sprintf(legendentry,"#lambda_{#theta} = 0,  #lambda_{#phi} = 0,  #lambda_{#theta#phi} = 0");
		plotLegend2->AddEntry(PDF_cth_HX_0,legendentry,"l");
	}
	plotLegend2->Draw(); plotLegend2->Draw();

	sprintf(filename,"%s/fit_HX_costh_%s.pdf",dirstruct,TreeBinID);
	c3->Print( filename );

	//2013-6-14
	//PDF_cth_HX->Sumw2(); PDF_cth_HX_0->Sumw2(); DATA_cth_HX->Sumw2(); 

	TH1D* ratio_best_To_zero_cth_HX = (TH1D*)PDF_cth_HX->Clone();
	ratio_best_To_zero_cth_HX->SetName("ratio_best_To_zero_cth_HX");
	ratio_best_To_zero_cth_HX->Divide(PDF_cth_HX_0);

	TH1D* ratio_data_To_zero_cth_HX = (TH1D*)DATA_cth_HX->Clone();
	ratio_data_To_zero_cth_HX->SetName("ratio_data_To_zero_cth_HX");
	ratio_data_To_zero_cth_HX->Divide(PDF_cth_HX_0);

	TH1D* ratio_data_To_0p1_cth_HX = (TH1D*)DATA_cth_HX->Clone();
	ratio_data_To_0p1_cth_HX->SetName("ratio_data_To_0p1_cth_HX");
	ratio_data_To_0p1_cth_HX->Divide(PDF_cth_HX_0p1);



	ratio_best_To_zero_cth_HX->GetXaxis()->SetTitle("cos#vartheta_{HX}");
	ratio_best_To_zero_cth_HX->GetXaxis()->SetLabelOffset(0.028);
	ratio_best_To_zero_cth_HX->GetXaxis()->SetTitleSize(0.05);
	ratio_best_To_zero_cth_HX->GetXaxis()->SetTickLength(-0.03);
	ratio_best_To_zero_cth_HX->GetXaxis()->SetTitleOffset(1.20);
	ratio_best_To_zero_cth_HX->GetYaxis()->SetTitle(" best / zero");
	ratio_best_To_zero_cth_HX->GetYaxis()->SetLabelOffset(0.032);
	ratio_best_To_zero_cth_HX->GetYaxis()->SetTitleSize(0.05);
	ratio_best_To_zero_cth_HX->GetYaxis()->SetTickLength(-0.03);
	ratio_best_To_zero_cth_HX->GetYaxis()->SetTitleOffset(1.55);
	ratio_best_To_zero_cth_HX->SetMinimum(0.);
	ratio_best_To_zero_cth_HX->SetMarkerStyle(20);
	ratio_best_To_zero_cth_HX->SetMarkerSize(0.5);
	ratio_best_To_zero_cth_HX->SetMarkerColor(kRed);
	ratio_best_To_zero_cth_HX->SetLineColor(kRed);
	ratio_best_To_zero_cth_HX->SetLineStyle(1);
	ratio_best_To_zero_cth_HX->SetLineWidth(2);
	ratio_best_To_zero_cth_HX->SetMaximum(2.);
	ratio_best_To_zero_cth_HX->Draw("p"); //pe2

	TLine *line0 = new TLine(-1, 1., 1., 1);
	line0->SetLineWidth( 2 );
	line0->SetLineStyle( 1 );
	line0->SetLineColor( kBlack );
	line0->Draw( "same" );

	sprintf(filename,"%s/fit_HX_costh_ratio_bestTOzero_%s.pdf",dirstruct,TreeBinID);
	c3->Print( filename );

	ratio_data_To_zero_cth_HX->GetXaxis()->SetTitle("cos#vartheta_{HX}");
	ratio_data_To_zero_cth_HX->GetXaxis()->SetLabelOffset(0.028);
	ratio_data_To_zero_cth_HX->GetXaxis()->SetTitleSize(0.05);
	ratio_data_To_zero_cth_HX->GetXaxis()->SetTickLength(-0.03);
	ratio_data_To_zero_cth_HX->GetXaxis()->SetTitleOffset(1.20);
	ratio_data_To_zero_cth_HX->GetYaxis()->SetTitle(" data / zero");
	ratio_data_To_zero_cth_HX->GetYaxis()->SetLabelOffset(0.032);
	ratio_data_To_zero_cth_HX->GetYaxis()->SetTitleSize(0.05);
	ratio_data_To_zero_cth_HX->GetYaxis()->SetTickLength(-0.03);
	ratio_data_To_zero_cth_HX->GetYaxis()->SetTitleOffset(1.55);
	ratio_data_To_zero_cth_HX->SetMinimum(0.);
	ratio_data_To_zero_cth_HX->SetMarkerStyle(20);
	ratio_data_To_zero_cth_HX->SetMarkerSize(0.5);
	ratio_data_To_zero_cth_HX->SetMarkerColor(kRed);
	ratio_data_To_zero_cth_HX->SetLineColor(kRed);
	ratio_data_To_zero_cth_HX->SetLineStyle(1);
	ratio_data_To_zero_cth_HX->SetLineWidth(2);
	ratio_data_To_zero_cth_HX->SetMaximum(2.);
	ratio_data_To_zero_cth_HX->Draw("p"); //pe2
	line0->Draw( "same" );

	sprintf(filename,"%s/fit_HX_costh_ratio_dataTOzero_%s.pdf",dirstruct,TreeBinID);
	c3->Print( filename );

	ratio_data_To_0p1_cth_HX->GetXaxis()->SetTitle("cos#vartheta_{HX}");
	ratio_data_To_0p1_cth_HX->GetXaxis()->SetLabelOffset(0.028);
	ratio_data_To_0p1_cth_HX->GetXaxis()->SetTitleSize(0.05);
	ratio_data_To_0p1_cth_HX->GetXaxis()->SetTickLength(-0.03);
	ratio_data_To_0p1_cth_HX->GetXaxis()->SetTitleOffset(1.20);
	ratio_data_To_0p1_cth_HX->GetYaxis()->SetTitle(" data / 0.1 curve");
	ratio_data_To_0p1_cth_HX->GetYaxis()->SetLabelOffset(0.032);
	ratio_data_To_0p1_cth_HX->GetYaxis()->SetTitleSize(0.05);
	ratio_data_To_0p1_cth_HX->GetYaxis()->SetTickLength(-0.03);
	ratio_data_To_0p1_cth_HX->GetYaxis()->SetTitleOffset(1.55);
	ratio_data_To_0p1_cth_HX->SetMinimum(0.);
	ratio_data_To_0p1_cth_HX->SetMarkerStyle(20);
	ratio_data_To_0p1_cth_HX->SetMarkerSize(0.5);
	ratio_data_To_0p1_cth_HX->SetMarkerColor(kRed);
	ratio_data_To_0p1_cth_HX->SetLineColor(kRed);
	ratio_data_To_0p1_cth_HX->SetLineStyle(1);
	ratio_data_To_0p1_cth_HX->SetLineWidth(2);
	ratio_data_To_0p1_cth_HX->SetMaximum(2.);
	ratio_data_To_0p1_cth_HX->Draw("p"); //pe2
	line0->Draw( "same" );

	sprintf(filename,"%s/fit_HX_costh_ratio_data_TO_0p1_%s.pdf",dirstruct,TreeBinID);
	c3->Print( filename );


	// phi_HX


	TH1D* PDF_ph_HX = (TH1D*)PDF_1_vs_ph_HX->Clone();
	PDF_ph_HX->SetName("PDF_ph_HX");
	PDF_ph_HX->Add( PDF_1_vs_ph_HX, PDF_cth2_vs_ph_HX, 1., lth_HX_best );
	PDF_ph_HX->Add( PDF_sth2c2ph_vs_ph_HX, lph_HX_best );
	PDF_ph_HX->Add( PDF_s2thcph_vs_ph_HX, ltp_HX_best );
	PDF_ph_HX->Scale( DATA_ph_HX->Integral() * PDF_ph_HX->GetNbinsX() / ( PDF_ph_HX->Integral()*DATA_ph_HX->GetNbinsX() ) );

	TH1D* PDF_ph_HX_p1 = (TH1D*)PDF_1_vs_ph_HX->Clone();
	PDF_ph_HX_p1->SetName("PDF_ph_HX_p1");
	PDF_ph_HX_p1->Add( PDF_1_vs_ph_HX, PDF_cth2_vs_ph_HX, 1., lth_phiplots );
	PDF_ph_HX_p1->Add( PDF_sth2c2ph_vs_ph_HX, 1. );
	PDF_ph_HX_p1->Scale( DATA_ph_HX->Integral() * PDF_ph_HX_p1->GetNbinsX() / ( PDF_ph_HX_p1->Integral()*DATA_ph_HX->GetNbinsX() ) );

	TH1D* PDF_ph_HX_m1 = (TH1D*)PDF_1_vs_ph_HX->Clone();
	PDF_ph_HX_m1->SetName("PDF_ph_HX_m1");
	PDF_ph_HX_m1->Add( PDF_1_vs_ph_HX, PDF_cth2_vs_ph_HX, 1., lth_phiplots );
	PDF_ph_HX_m1->Add( PDF_sth2c2ph_vs_ph_HX, -1. );
	PDF_ph_HX_m1->Scale( DATA_ph_HX->Integral() * PDF_ph_HX_m1->GetNbinsX() / ( PDF_ph_HX_m1->Integral()*DATA_ph_HX->GetNbinsX() ) );

	TH1D* PDF_ph_HX_0 = (TH1D*)PDF_1_vs_ph_HX->Clone();
	PDF_ph_HX_0->SetName("PDF_ph_HX_0");
	PDF_ph_HX_0->Add( PDF_1_vs_ph_HX, PDF_cth2_vs_ph_HX, 1., lth_phiplots );
	PDF_ph_HX_0->Add( PDF_sth2c2ph_vs_ph_HX, 0. );
	PDF_ph_HX_0->Scale( DATA_ph_HX->Integral() * PDF_ph_HX_0->GetNbinsX() / ( PDF_ph_HX_0->Integral()*DATA_ph_HX->GetNbinsX() ) );



	plotMax = plotborder * PDF_ph_HX->GetMaximum();
	if ( plotborder * PDF_ph_HX_p1->GetMaximum() > plotMax ) plotMax = plotborder * PDF_ph_HX_p1->GetMaximum();
	if ( plotborder * PDF_ph_HX_m1->GetMaximum() > plotMax ) plotMax = plotborder * PDF_ph_HX_m1->GetMaximum();

	PDF_ph_HX->GetXaxis()->SetTitle("#varphi_{HX}");
	PDF_ph_HX->GetXaxis()->SetLabelOffset(0.028);
	PDF_ph_HX->GetXaxis()->SetTitleSize(0.05);
	PDF_ph_HX->GetXaxis()->SetTickLength(-0.03);
	PDF_ph_HX->GetXaxis()->SetTitleOffset(1.20);
	PDF_ph_HX->GetYaxis()->SetTitle("event PDF [a.u.]");
	PDF_ph_HX->GetYaxis()->SetLabelOffset(0.032);
	PDF_ph_HX->GetYaxis()->SetTitleSize(0.05);
	PDF_ph_HX->GetYaxis()->SetTickLength(-0.03);
	PDF_ph_HX->GetYaxis()->SetTitleOffset(1.55);
	PDF_ph_HX->SetMinimum(0.);
	PDF_ph_HX->SetLineColor(kRed);
	PDF_ph_HX->SetLineStyle(1);
	PDF_ph_HX->SetLineWidth(2);
	PDF_ph_HX->SetMaximum(plotMax);
	PDF_ph_HX->Draw("L");

	PDF_ph_HX_p1->SetLineColor(kRed);
	PDF_ph_HX_p1->SetLineStyle(2);
	PDF_ph_HX_p1->SetLineWidth(2);
	PDF_ph_HX_p1->Draw("L same");

	PDF_ph_HX_m1->SetLineColor(kRed);
	PDF_ph_HX_m1->SetLineStyle(3);
	PDF_ph_HX_m1->SetLineWidth(2);
	PDF_ph_HX_m1->Draw("L same");

	PDF_ph_HX_0->SetLineColor(kBlue); //(kRed);
	PDF_ph_HX_0->SetLineStyle(10);
	PDF_ph_HX_0->SetLineWidth(2);
	if(plotZeros) PDF_ph_HX_0->Draw("L same");

	DATA_ph_HX->SetLineColor(kBlack);
	DATA_ph_HX->Draw("E same");

	plotLegend2=new TLegend(0.2,0.8,0.75,0.95);
	plotLegend2->SetFillColor(kWhite);
	//plotLegend2->SetTextFont(72);
	plotLegend2->SetTextSize(0.035);
	plotLegend2->SetBorderSize(1);

	sprintf(legendentry,"Best fit curve");
	plotLegend2->AddEntry(PDF_ph_HX,legendentry,"l");
	//sprintf(legendentry,"#lambda_{#theta} = +0.5,  #lambda_{#phi} = +0.5,  #lambda_{#theta#phi} = 0");
	sprintf(legendentry,"#lambda_{#theta} = 0,  #lambda_{#phi} = +0.5,  #lambda_{#theta#phi} = 0");
	plotLegend2->AddEntry(PDF_ph_HX_p1,legendentry,"l");
	//sprintf(legendentry,"#lambda_{#theta} = +0.5,  #lambda_{#phi} = -0.5,  #lambda_{#theta#phi} = 0");
	sprintf(legendentry,"#lambda_{#theta} = 0,  #lambda_{#phi} = -0.5,  #lambda_{#theta#phi} = 0");
	plotLegend2->AddEntry(PDF_ph_HX_m1,legendentry,"l");
	if(plotZeros){
		sprintf(legendentry,"#lambda_{#theta} = 0,  #lambda_{#phi} = 0,  #lambda_{#theta#phi} = 0");
		plotLegend2->AddEntry(PDF_ph_HX_0,legendentry,"l");
	}
	plotLegend2->Draw(); plotLegend2->Draw();

	sprintf(filename,"%s/fit_HX_phi_%s.pdf",dirstruct,TreeBinID);
	c3->Print( filename );

	// 2013-6-17
	char costh_phi_fileName[500];
	sprintf(costh_phi_fileName,"%s/costh_phi_HX_%s.root",dirstruct,TreeBinID);
	TFile *costh_phi_file = new TFile(costh_phi_fileName,"RECREATE");
	costh_phi_file->cd();
	PDF_cth_HX_0->Write();
	PDF_ph_HX_0->Write();
	ratio_best_To_zero_cth_HX->Write();
	ratio_data_To_zero_cth_HX->Write();
	ratio_data_To_0p1_cth_HX->Write();
	costh_phi_file->Close();

	// phith_HX

	TH1D* PDF_phth_HX = (TH1D*)PDF_1_vs_phth_HX->Clone();
	PDF_phth_HX->SetName("PDF_phth_HX");
	PDF_phth_HX->Add( PDF_1_vs_phth_HX, PDF_cth2_vs_phth_HX, 1., lth_HX_best );
	PDF_phth_HX->Add( PDF_sth2c2ph_vs_phth_HX, lph_HX_best );
	PDF_phth_HX->Add( PDF_s2thcph_vs_phth_HX, ltp_HX_best );
	PDF_phth_HX->Scale( DATA_phth_HX->Integral() * PDF_phth_HX->GetNbinsX() / ( PDF_phth_HX->Integral()*DATA_phth_HX->GetNbinsX() ) );

	TH1D* PDF_phth_HX_p1 = (TH1D*)PDF_1_vs_phth_HX->Clone();
	PDF_phth_HX_p1->SetName("PDF_phth_HX_p1");
	PDF_phth_HX_p1->Add( PDF_1_vs_phth_HX, PDF_sth2c2ph_vs_phth_HX, 1., lth_phithplots );
	PDF_phth_HX_p1->Add( PDF_s2thcph_vs_phth_HX, sqrt(2.)/2. );
	PDF_phth_HX_p1->Scale( DATA_phth_HX->Integral() * PDF_phth_HX_p1->GetNbinsX() / ( PDF_phth_HX_p1->Integral()*DATA_phth_HX->GetNbinsX() ) );

	TH1D* PDF_phth_HX_m1 = (TH1D*)PDF_1_vs_phth_HX->Clone();
	PDF_phth_HX_m1->SetName("PDF_phth_HX_m1");
	PDF_phth_HX_m1->Add( PDF_1_vs_phth_HX, PDF_sth2c2ph_vs_phth_HX, 1., lth_phithplots );
	PDF_phth_HX_m1->Add( PDF_s2thcph_vs_phth_HX, -sqrt(2.)/2. );
	PDF_phth_HX_m1->Scale( DATA_phth_HX->Integral() * PDF_phth_HX_m1->GetNbinsX() / ( PDF_phth_HX_m1->Integral()*DATA_phth_HX->GetNbinsX() ) );

	TH1D* PDF_phth_HX_0 = (TH1D*)PDF_1_vs_phth_HX->Clone();
	PDF_phth_HX_0->SetName("PDF_phth_HX_0");
	PDF_phth_HX_0->Add( PDF_1_vs_phth_HX, PDF_sth2c2ph_vs_phth_HX, 1., lth_phithplots );
	PDF_phth_HX_0->Add( PDF_s2thcph_vs_phth_HX, 0. );
	PDF_phth_HX_0->Scale( DATA_phth_HX->Integral() * PDF_phth_HX_0->GetNbinsX() / ( PDF_phth_HX_0->Integral()*DATA_phth_HX->GetNbinsX() ) );




	plotMax = plotborder * PDF_phth_HX->GetMaximum();
	if ( plotborder * PDF_phth_HX_p1->GetMaximum() > plotMax ) plotMax = plotborder * PDF_phth_HX_p1->GetMaximum();
	if ( plotborder * PDF_phth_HX_m1->GetMaximum() > plotMax ) plotMax = plotborder * PDF_phth_HX_m1->GetMaximum();

	PDF_phth_HX->GetXaxis()->SetTitle("#tilde{#varphi}_{HX}");
	PDF_phth_HX->GetXaxis()->SetLabelOffset(0.028);
	PDF_phth_HX->GetXaxis()->SetTitleSize(0.05);
	PDF_phth_HX->GetXaxis()->SetTickLength(-0.03);
	PDF_phth_HX->GetXaxis()->SetTitleOffset(1.20);
	PDF_phth_HX->GetYaxis()->SetTitle("event PDF [a.u.]");
	PDF_phth_HX->GetYaxis()->SetLabelOffset(0.032);
	PDF_phth_HX->GetYaxis()->SetTitleSize(0.05);
	PDF_phth_HX->GetYaxis()->SetTickLength(-0.03);
	PDF_phth_HX->GetYaxis()->SetTitleOffset(1.55);
	PDF_phth_HX->SetMinimum(0.);
	PDF_phth_HX->SetLineColor(kRed);
	PDF_phth_HX->SetLineStyle(1);
	PDF_phth_HX->SetLineWidth(2);
	PDF_phth_HX->SetMaximum(plotMax);
	PDF_phth_HX->Draw("L");

	PDF_phth_HX_p1->SetLineColor(kRed);
	PDF_phth_HX_p1->SetLineStyle(2);
	PDF_phth_HX_p1->SetLineWidth(2);
	PDF_phth_HX_p1->Draw("L same");

	PDF_phth_HX_m1->SetLineColor(kRed);
	PDF_phth_HX_m1->SetLineStyle(3);
	PDF_phth_HX_m1->SetLineWidth(2);
	PDF_phth_HX_m1->Draw("L same");

	PDF_phth_HX_0->SetLineColor(kBlue); //kRed);
	PDF_phth_HX_0->SetLineStyle(10);
	PDF_phth_HX_0->SetLineWidth(2);
	if(plotZeros) PDF_phth_HX_0->Draw("L same");

	DATA_phth_HX->SetLineColor(kBlack);
	DATA_phth_HX->Draw("E same");

	plotLegend2=new TLegend(0.2,0.8,0.75,0.95);
	plotLegend2->SetFillColor(kWhite);
	//plotLegend2->SetTextFont(72);
	plotLegend2->SetTextSize(0.035);
	plotLegend2->SetBorderSize(1);

	sprintf(legendentry,"Best fit curve");
	plotLegend2->AddEntry(PDF_phth_HX,legendentry,"l");
	sprintf(legendentry,"#lambda_{#theta} = +0.5,  #lambda_{#phi} = 0,  #lambda_{#theta#phi} = +#sqrt{2}/2");
	plotLegend2->AddEntry(PDF_phth_HX_p1,legendentry,"l");
	sprintf(legendentry,"#lambda_{#theta} = +0.5,  #lambda_{#phi} = 0,  #lambda_{#theta#phi} = -#sqrt{2}/2");
	plotLegend2->AddEntry(PDF_phth_HX_m1,legendentry,"l");
	if(plotZeros){
		sprintf(legendentry,"#lambda_{#theta} = 0,  #lambda_{#phi} = 0,  #lambda_{#theta#phi} = 0");
		plotLegend2->AddEntry(PDF_phth_HX_0,legendentry,"l");
	}
	plotLegend2->Draw(); plotLegend2->Draw();

	sprintf(filename,"%s/fit_HX_phith_%s.pdf",dirstruct,TreeBinID);
	c3->Print( filename );



	// costh_PX

	TH1D* PDF_cth_PX = (TH1D*)PDF_1_vs_cth_PX->Clone();
	PDF_cth_PX->SetName("PDF_cth_PX");
	PDF_cth_PX->Add( PDF_1_vs_cth_PX, PDF_cth2_vs_cth_PX, 1., lth_PX_best );
	PDF_cth_PX->Add( PDF_sth2c2ph_vs_cth_PX, lph_PX_best );
	PDF_cth_PX->Add( PDF_s2thcph_vs_cth_PX, ltp_PX_best );
	PDF_cth_PX->Scale( DATA_cth_PX->Integral() * PDF_cth_PX->GetNbinsX() / ( PDF_cth_PX->Integral()*DATA_cth_PX->GetNbinsX() ) );

	TH1D* PDF_cth_PX_p1 = (TH1D*)PDF_1_vs_cth_PX->Clone();
	PDF_cth_PX_p1->SetName("PDF_cth_PX_p1");
	PDF_cth_PX_p1->Add( PDF_1_vs_cth_PX, PDF_cth2_vs_cth_PX, 1., 1. );
	PDF_cth_PX_p1->Scale( DATA_cth_PX->Integral() * PDF_cth_PX_p1->GetNbinsX() / ( PDF_cth_PX_p1->Integral()*DATA_cth_PX->GetNbinsX() ) );

	TH1D* PDF_cth_PX_m1 = (TH1D*)PDF_1_vs_cth_PX->Clone();
	PDF_cth_PX_m1->SetName("PDF_cth_PX_m1");
	PDF_cth_PX_m1->Add( PDF_1_vs_cth_PX, PDF_cth2_vs_cth_PX, 1., -1. );
	PDF_cth_PX_m1->Scale( DATA_cth_PX->Integral() * PDF_cth_PX_m1->GetNbinsX() / ( PDF_cth_PX_m1->Integral()*DATA_cth_PX->GetNbinsX() ) );

	TH1D* PDF_cth_PX_0 = (TH1D*)PDF_1_vs_cth_PX->Clone();
	PDF_cth_PX_0->SetName("PDF_cth_PX_0");
	PDF_cth_PX_0->Add( PDF_1_vs_cth_PX, PDF_cth2_vs_cth_PX, 1., lthplot_zeroLine );
	PDF_cth_PX_0->Scale( DATA_cth_PX->Integral() * PDF_cth_PX_0->GetNbinsX() / ( PDF_cth_PX_0->Integral()*DATA_cth_PX->GetNbinsX() ) );



	plotMax = plotborder * PDF_cth_PX->GetMaximum();
	if ( plotborder * PDF_cth_PX_p1->GetMaximum() > plotMax ) plotMax = plotborder * PDF_cth_PX_p1->GetMaximum();
	if ( plotborder * PDF_cth_PX_m1->GetMaximum() > plotMax ) plotMax = plotborder * PDF_cth_PX_m1->GetMaximum();

	PDF_cth_PX->GetXaxis()->SetTitle("cos#vartheta_{PX}");
	PDF_cth_PX->GetXaxis()->SetLabelOffset(0.028);
	PDF_cth_PX->GetXaxis()->SetTitleSize(0.05);
	PDF_cth_PX->GetXaxis()->SetTickLength(-0.03);
	PDF_cth_PX->GetXaxis()->SetTitleOffset(1.20);
	PDF_cth_PX->GetYaxis()->SetTitle("event PDF [a.u.]");
	PDF_cth_PX->GetYaxis()->SetLabelOffset(0.032);
	PDF_cth_PX->GetYaxis()->SetTitleSize(0.05);
	PDF_cth_PX->GetYaxis()->SetTickLength(-0.03);
	PDF_cth_PX->GetYaxis()->SetTitleOffset(1.55);
	PDF_cth_PX->SetMinimum(0.);
	PDF_cth_PX->SetLineColor(kRed);
	PDF_cth_PX->SetLineStyle(1);
	PDF_cth_PX->SetLineWidth(2);
	PDF_cth_PX->SetMaximum(plotMax);
	PDF_cth_PX->Draw("L");

	PDF_cth_PX_p1->SetLineColor(kRed);
	PDF_cth_PX_p1->SetLineStyle(2);
	PDF_cth_PX_p1->SetLineWidth(2);
	PDF_cth_PX_p1->Draw("L same");

	PDF_cth_PX_m1->SetLineColor(kRed);
	PDF_cth_PX_m1->SetLineStyle(3);
	PDF_cth_PX_m1->SetLineWidth(2);
	PDF_cth_PX_m1->Draw("L same");

	PDF_cth_PX_0->SetLineColor(kBlue); //(kRed);
	PDF_cth_PX_0->SetLineStyle(10);
	PDF_cth_PX_0->SetLineWidth(2);
	if(plotZeros) PDF_cth_PX_0->Draw("L same");

	DATA_cth_PX->SetLineColor(kBlack);
	DATA_cth_PX->Draw("E same");

	plotLegend2=new TLegend(0.2,0.8,0.75,0.95);
	plotLegend2->SetFillColor(kWhite);
	//plotLegend2->SetTextFont(72);
	plotLegend2->SetTextSize(0.035);
	plotLegend2->SetBorderSize(1);

	sprintf(legendentry,"Best fit curve");
	plotLegend2->AddEntry(PDF_cth_PX,legendentry,"l");
	sprintf(legendentry,"#lambda_{#theta} = +1,  #lambda_{#phi} = 0,  #lambda_{#theta#phi} = 0");
	plotLegend2->AddEntry(PDF_cth_PX_p1,legendentry,"l");
	sprintf(legendentry,"#lambda_{#theta} = -1,  #lambda_{#phi} = 0,  #lambda_{#theta#phi} = 0");
	plotLegend2->AddEntry(PDF_cth_PX_m1,legendentry,"l");
	if(plotZeros){
		sprintf(legendentry,"#lambda_{#theta} = 0,  #lambda_{#phi} = 0,  #lambda_{#theta#phi} = 0");
		plotLegend2->AddEntry(PDF_cth_PX_0,legendentry,"l");
	}
	plotLegend2->Draw(); plotLegend2->Draw();

	sprintf(filename,"%s/fit_PX_costh_%s.pdf",dirstruct,TreeBinID);
	c3->Print( filename );


	// phi_PX

	TH1D* PDF_ph_PX = (TH1D*)PDF_1_vs_ph_PX->Clone();
	PDF_ph_PX->SetName("PDF_ph_PX");
	PDF_ph_PX->Add( PDF_1_vs_ph_PX, PDF_cth2_vs_ph_PX, 1., lth_PX_best );
	PDF_ph_PX->Add( PDF_sth2c2ph_vs_ph_PX, lph_PX_best );
	PDF_ph_PX->Add( PDF_s2thcph_vs_ph_PX, ltp_PX_best );
	PDF_ph_PX->Scale( DATA_ph_PX->Integral() * PDF_ph_PX->GetNbinsX() / ( PDF_ph_PX->Integral()*DATA_ph_PX->GetNbinsX() ) );

	TH1D* PDF_ph_PX_p1 = (TH1D*)PDF_1_vs_ph_PX->Clone();
	PDF_ph_PX_p1->SetName("PDF_ph_PX_p1");
	PDF_ph_PX_p1->Add( PDF_1_vs_ph_PX, PDF_cth2_vs_ph_PX, 1., lth_phiplots );
	PDF_ph_PX_p1->Add( PDF_sth2c2ph_vs_ph_PX, 1. );
	PDF_ph_PX_p1->Scale( DATA_ph_PX->Integral() * PDF_ph_PX_p1->GetNbinsX() / ( PDF_ph_PX_p1->Integral()*DATA_ph_PX->GetNbinsX() ) );

	TH1D* PDF_ph_PX_m1 = (TH1D*)PDF_1_vs_ph_PX->Clone();
	PDF_ph_PX_m1->SetName("PDF_ph_PX_m1");
	PDF_ph_PX_m1->Add( PDF_1_vs_ph_PX, PDF_cth2_vs_ph_PX, 1., lth_phiplots );
	PDF_ph_PX_m1->Add( PDF_sth2c2ph_vs_ph_PX, -1. );
	PDF_ph_PX_m1->Scale( DATA_ph_PX->Integral() * PDF_ph_PX_m1->GetNbinsX() / ( PDF_ph_PX_m1->Integral()*DATA_ph_PX->GetNbinsX() ) );

	TH1D* PDF_ph_PX_0 = (TH1D*)PDF_1_vs_ph_PX->Clone();
	PDF_ph_PX_0->SetName("PDF_ph_PX_0");
	PDF_ph_PX_0->Add( PDF_1_vs_ph_PX, PDF_cth2_vs_ph_PX, 1., lth_phiplots );
	PDF_ph_PX_0->Add( PDF_sth2c2ph_vs_ph_PX, 0. );
	PDF_ph_PX_0->Scale( DATA_ph_PX->Integral() * PDF_ph_PX_0->GetNbinsX() / ( PDF_ph_PX_0->Integral()*DATA_ph_PX->GetNbinsX() ) );



	plotMax = plotborder * PDF_ph_PX->GetMaximum();
	if ( plotborder * PDF_ph_PX_p1->GetMaximum() > plotMax ) plotMax = plotborder * PDF_ph_PX_p1->GetMaximum();
	if ( plotborder * PDF_ph_PX_m1->GetMaximum() > plotMax ) plotMax = plotborder * PDF_ph_PX_m1->GetMaximum();

	PDF_ph_PX->GetXaxis()->SetTitle("#varphi_{PX}");
	PDF_ph_PX->GetXaxis()->SetLabelOffset(0.028);
	PDF_ph_PX->GetXaxis()->SetTitleSize(0.05);
	PDF_ph_PX->GetXaxis()->SetTickLength(-0.03);
	PDF_ph_PX->GetXaxis()->SetTitleOffset(1.20);
	PDF_ph_PX->GetYaxis()->SetTitle("event PDF [a.u.]");
	PDF_ph_PX->GetYaxis()->SetLabelOffset(0.032);
	PDF_ph_PX->GetYaxis()->SetTitleSize(0.05);
	PDF_ph_PX->GetYaxis()->SetTickLength(-0.03);
	PDF_ph_PX->GetYaxis()->SetTitleOffset(1.55);
	PDF_ph_PX->SetMinimum(0.);
	PDF_ph_PX->SetLineColor(kRed);
	PDF_ph_PX->SetLineStyle(1);
	PDF_ph_PX->SetLineWidth(2);
	PDF_ph_PX->SetMaximum(plotMax);
	PDF_ph_PX->Draw("L");

	PDF_ph_PX_p1->SetLineColor(kRed);
	PDF_ph_PX_p1->SetLineStyle(2);
	PDF_ph_PX_p1->SetLineWidth(2);
	PDF_ph_PX_p1->Draw("L same");

	PDF_ph_PX_m1->SetLineColor(kRed);
	PDF_ph_PX_m1->SetLineStyle(3);
	PDF_ph_PX_m1->SetLineWidth(2);
	PDF_ph_PX_m1->Draw("L same");

	PDF_ph_PX_0->SetLineColor(kBlue); //(kRed);
	PDF_ph_PX_0->SetLineStyle(10);
	PDF_ph_PX_0->SetLineWidth(2);
	if(plotZeros) PDF_ph_PX_0->Draw("L same");

	DATA_ph_PX->SetLineColor(kBlack);
	DATA_ph_PX->Draw("E same");

	plotLegend2=new TLegend(0.2,0.8,0.75,0.95);
	plotLegend2->SetFillColor(kWhite);
	//plotLegend2->SetTextFont(72);
	plotLegend2->SetTextSize(0.035);
	plotLegend2->SetBorderSize(1);

	sprintf(legendentry,"Best fit curve");
	plotLegend2->AddEntry(PDF_ph_PX,legendentry,"l");
	//sprintf(legendentry,"#lambda_{#theta} = +0.5,  #lambda_{#phi} = +0.5,  #lambda_{#theta#phi} = 0");
	sprintf(legendentry,"#lambda_{#theta} = 0,  #lambda_{#phi} = +0.5,  #lambda_{#theta#phi} = 0");
	plotLegend2->AddEntry(PDF_ph_PX_p1,legendentry,"l");
	//sprintf(legendentry,"#lambda_{#theta} = +0.5,  #lambda_{#phi} = -0.5,  #lambda_{#theta#phi} = 0");
	sprintf(legendentry,"#lambda_{#theta} = 0,  #lambda_{#phi} = -0.5,  #lambda_{#theta#phi} = 0");
	plotLegend2->AddEntry(PDF_ph_PX_m1,legendentry,"l");
	if(plotZeros){
		sprintf(legendentry,"#lambda_{#theta} = 0,  #lambda_{#phi} = 0,  #lambda_{#theta#phi} = 0");
		plotLegend2->AddEntry(PDF_ph_PX_0,legendentry,"l");
	}
	plotLegend2->Draw(); plotLegend2->Draw();

	sprintf(filename,"%s/fit_PX_phi_%s.pdf",dirstruct,TreeBinID);
	c3->Print( filename );


	// phith_PX

	TH1D* PDF_phth_PX = (TH1D*)PDF_1_vs_phth_PX->Clone();
	PDF_phth_PX->SetName("PDF_phth_PX");
	PDF_phth_PX->Add( PDF_1_vs_phth_PX, PDF_cth2_vs_phth_PX, 1., lth_PX_best );
	PDF_phth_PX->Add( PDF_sth2c2ph_vs_phth_PX, lph_PX_best );
	PDF_phth_PX->Add( PDF_s2thcph_vs_phth_PX, ltp_PX_best );
	PDF_phth_PX->Scale( DATA_phth_PX->Integral() * PDF_phth_PX->GetNbinsX() / ( PDF_phth_PX->Integral()*DATA_phth_PX->GetNbinsX() ) );

	TH1D* PDF_phth_PX_p1 = (TH1D*)PDF_1_vs_phth_PX->Clone();
	PDF_phth_PX_p1->SetName("PDF_phth_PX_p1");
	PDF_phth_PX_p1->Add( PDF_1_vs_phth_PX, PDF_sth2c2ph_vs_phth_PX, 1., lth_phithplots );
	PDF_phth_PX_p1->Add( PDF_s2thcph_vs_phth_PX, sqrt(2.)/2. );
	PDF_phth_PX_p1->Scale( DATA_phth_PX->Integral() * PDF_phth_PX_p1->GetNbinsX() / ( PDF_phth_PX_p1->Integral()*DATA_phth_PX->GetNbinsX() ) );

	TH1D* PDF_phth_PX_m1 = (TH1D*)PDF_1_vs_phth_PX->Clone();
	PDF_phth_PX_m1->SetName("PDF_phth_PX_m1");
	PDF_phth_PX_m1->Add( PDF_1_vs_phth_PX, PDF_sth2c2ph_vs_phth_PX, 1., lth_phithplots );
	PDF_phth_PX_m1->Add( PDF_s2thcph_vs_phth_PX, -sqrt(2.)/2. );
	PDF_phth_PX_m1->Scale( DATA_phth_PX->Integral() * PDF_phth_PX_m1->GetNbinsX() / ( PDF_phth_PX_m1->Integral()*DATA_phth_PX->GetNbinsX() ) );

	TH1D* PDF_phth_PX_0 = (TH1D*)PDF_1_vs_phth_PX->Clone();
	PDF_phth_PX_0->SetName("PDF_phth_PX_0");
	PDF_phth_PX_0->Add( PDF_1_vs_phth_PX, PDF_sth2c2ph_vs_phth_PX, 1., lth_phithplots );
	PDF_phth_PX_0->Add( PDF_s2thcph_vs_phth_PX, 0. );
	PDF_phth_PX_0->Scale( DATA_phth_PX->Integral() * PDF_phth_PX_0->GetNbinsX() / ( PDF_phth_PX_0->Integral()*DATA_phth_PX->GetNbinsX() ) );



	plotMax = plotborder * PDF_phth_PX->GetMaximum();
	if ( plotborder * PDF_phth_PX_p1->GetMaximum() > plotMax ) plotMax = plotborder * PDF_phth_PX_p1->GetMaximum();
	if ( plotborder * PDF_phth_PX_m1->GetMaximum() > plotMax ) plotMax = plotborder * PDF_phth_PX_m1->GetMaximum();

	PDF_phth_PX->GetXaxis()->SetTitle("#tilde{#varphi}_{PX}");
	PDF_phth_PX->GetXaxis()->SetLabelOffset(0.028);
	PDF_phth_PX->GetXaxis()->SetTitleSize(0.05);
	PDF_phth_PX->GetXaxis()->SetTickLength(-0.03);
	PDF_phth_PX->GetXaxis()->SetTitleOffset(1.20);
	PDF_phth_PX->GetYaxis()->SetTitle("event PDF [a.u.]");
	PDF_phth_PX->GetYaxis()->SetLabelOffset(0.032);
	PDF_phth_PX->GetYaxis()->SetTitleSize(0.05);
	PDF_phth_PX->GetYaxis()->SetTickLength(-0.03);
	PDF_phth_PX->GetYaxis()->SetTitleOffset(1.55);
	PDF_phth_PX->SetMinimum(0.);
	PDF_phth_PX->SetLineColor(kRed);
	PDF_phth_PX->SetLineStyle(1);
	PDF_phth_PX->SetLineWidth(2);
	PDF_phth_PX->SetMaximum(plotMax);
	PDF_phth_PX->Draw("L");

	PDF_phth_PX_p1->SetLineColor(kRed);
	PDF_phth_PX_p1->SetLineStyle(2);
	PDF_phth_PX_p1->SetLineWidth(2);
	PDF_phth_PX_p1->Draw("L same");

	PDF_phth_PX_m1->SetLineColor(kRed);
	PDF_phth_PX_m1->SetLineStyle(3);
	PDF_phth_PX_m1->SetLineWidth(2);
	PDF_phth_PX_m1->Draw("L same");

	PDF_phth_PX_0->SetLineColor(kBlue); //(kRed);
	PDF_phth_PX_0->SetLineStyle(10);
	PDF_phth_PX_0->SetLineWidth(2);
	if(plotZeros) PDF_phth_PX_0->Draw("L same");

	DATA_phth_PX->SetLineColor(kBlack);
	DATA_phth_PX->Draw("E same");

	plotLegend2=new TLegend(0.2,0.8,0.75,0.95);
	plotLegend2->SetFillColor(kWhite);
	//plotLegend2->SetTextFont(72);
	plotLegend2->SetTextSize(0.035);
	plotLegend2->SetBorderSize(1);

	sprintf(legendentry,"Best fit curve");
	plotLegend2->AddEntry(PDF_phth_PX,legendentry,"l");
	sprintf(legendentry,"#lambda_{#theta} = +0.5,  #lambda_{#phi} = 0,  #lambda_{#theta#phi} = +#sqrt{2}/2");
	plotLegend2->AddEntry(PDF_phth_PX_p1,legendentry,"l");
	sprintf(legendentry,"#lambda_{#theta} = +0.5,  #lambda_{#phi} = 0,  #lambda_{#theta#phi} = -#sqrt{2}/2");
	plotLegend2->AddEntry(PDF_phth_PX_m1,legendentry,"l");
	if(plotZeros){
		sprintf(legendentry,"#lambda_{#theta} = 0,  #lambda_{#phi} = 0,  #lambda_{#theta#phi} = 0");
		plotLegend2->AddEntry(PDF_phth_PX_0,legendentry,"l");
	}
	plotLegend2->Draw(); plotLegend2->Draw();

	sprintf(filename,"%s/fit_PX_phith_%s.pdf",dirstruct,TreeBinID);
	c3->Print( filename );




	// cosalpha

	TH1D* PDF_calpha = (TH1D*)PDF_1_vs_calpha->Clone();
	PDF_calpha->SetName("PDF_calpha");
	PDF_calpha->Add( PDF_1_vs_calpha, PDF_cth2_vs_calpha, 1., lth_PX_best );
	PDF_calpha->Add( PDF_sth2c2ph_vs_calpha, lph_PX_best );
	PDF_calpha->Add( PDF_s2thcph_vs_calpha, ltp_PX_best );
	PDF_calpha->Scale( DATA_calpha->Integral() * PDF_calpha->GetNbinsX() / ( PDF_calpha->Integral()*DATA_calpha->GetNbinsX() ) );

	TH1D* PDF_calpha_p1 = (TH1D*)PDF_1_vs_calpha->Clone();
	PDF_calpha_p1->SetName("PDF_calpha_p1");
	PDF_calpha_p1->Add( PDF_1_vs_calpha, PDF_cth2_vs_calpha, 1., 1. );
	PDF_calpha_p1->Add( PDF_sth2c2ph_vs_calpha, 1. );
	PDF_calpha_p1->Scale( DATA_calpha->Integral() * PDF_calpha_p1->GetNbinsX() / ( PDF_calpha_p1->Integral()*DATA_calpha->GetNbinsX() ) );

	TH1D* PDF_calpha_m1 = (TH1D*)PDF_1_vs_calpha->Clone();
	PDF_calpha_m1->SetName("PDF_calpha_m1");
	PDF_calpha_m1->Add( PDF_1_vs_calpha, PDF_cth2_vs_calpha, 1., -1. );
	PDF_calpha_m1->Scale( DATA_calpha->Integral() * PDF_calpha_m1->GetNbinsX() / ( PDF_calpha_m1->Integral()*DATA_calpha->GetNbinsX() ) );



	plotMax = plotborder * PDF_calpha->GetMaximum();
	if ( plotborder * PDF_calpha_p1->GetMaximum() > plotMax ) plotMax = plotborder * PDF_calpha_p1->GetMaximum();
	if ( plotborder * PDF_calpha_m1->GetMaximum() > plotMax ) plotMax = plotborder * PDF_calpha_m1->GetMaximum();

	PDF_calpha->GetXaxis()->SetTitle("cos#alpha");
	PDF_calpha->GetXaxis()->SetLabelOffset(0.028);
	PDF_calpha->GetXaxis()->SetTitleSize(0.05);
	PDF_calpha->GetXaxis()->SetTickLength(-0.03);
	PDF_calpha->GetXaxis()->SetTitleOffset(1.20);
	PDF_calpha->GetYaxis()->SetTitle("event PDF [a.u.]");
	PDF_calpha->GetYaxis()->SetLabelOffset(0.032);
	PDF_calpha->GetYaxis()->SetTitleSize(0.05);
	PDF_calpha->GetYaxis()->SetTickLength(-0.03);
	PDF_calpha->GetYaxis()->SetTitleOffset(1.55);
	PDF_calpha->SetMinimum(0.);
	PDF_calpha->SetLineColor(kRed);
	PDF_calpha->SetLineStyle(1);
	PDF_calpha->SetLineWidth(2);
	PDF_calpha->SetMaximum(plotMax);
	PDF_calpha->Draw("L");

	PDF_calpha_p1->SetLineColor(kRed);
	PDF_calpha_p1->SetLineStyle(2);
	PDF_calpha_p1->Draw("L same");

	PDF_calpha_m1->SetLineColor(kRed);
	PDF_calpha_m1->SetLineStyle(3);
	PDF_calpha_m1->Draw("L same");

	DATA_calpha->SetLineColor(kBlack);
	DATA_calpha->Draw("E same");

	plotLegend2=new TLegend(0.2,0.8,0.75,0.95);
	plotLegend2->SetFillColor(kWhite);
	//plotLegend2->SetTextFont(72);
	plotLegend2->SetTextSize(0.035);
	plotLegend2->SetBorderSize(1);

	sprintf(legendentry,"Best fit curve");
	plotLegend2->AddEntry(PDF_calpha,legendentry,"l");
	sprintf(legendentry,"#lambda_{#theta} = +1,  #lambda_{#phi} = 0,  #lambda_{#theta#phi} = 0");
	plotLegend2->AddEntry(PDF_calpha_p1,legendentry,"l");
	sprintf(legendentry,"#lambda_{#theta} = -1,  #lambda_{#phi} = 0,  #lambda_{#theta#phi} = 0");
	plotLegend2->AddEntry(PDF_calpha_m1,legendentry,"l");
	plotLegend2->Draw(); plotLegend2->Draw();

	sprintf(filename,"%s/fit_cosalpha_%s.pdf",dirstruct,TreeBinID);
	c3->Print( filename );




	TLine  *bkgLine;


	double minBkgPlot=0;

	c3->SetRightMargin(0.03);

	TH1* background_rap_test = (TH1*)results->Get("background_rap_test");
	background_rap_test->Scale(scaleFits);
	background_rap_test->GetXaxis()->SetTitle("|y|");
	background_rap_test->GetYaxis()->SetTitle("Sampled / Original");
	background_rap_test->GetYaxis()->SetTitleOffset(1.55);
	background_rap_test->Draw();
	background_rap_test->SetMinimum(minBkgPlot);

	bkgLine = new TLine( background_rap_test->GetBinCenter(1)-background_rap_test->GetBinWidth(1)/2., 1., background_rap_test->GetBinCenter(background_rap_test->GetNbinsX())+background_rap_test->GetBinWidth(background_rap_test->GetNbinsX())/2.,1.);
	bkgLine->SetLineWidth( 1 );
	bkgLine->SetLineStyle( 2 );
	bkgLine->SetLineColor( kBlack);
	bkgLine->Draw( "same" );

	sprintf(filename,"%s/fit_background_rap_test_%s.pdf",dirstruct,TreeBinID);
	c3->Print( filename );

	TH1* background_mass_test = (TH1*)results->Get("background_mass_test");
	background_mass_test->Scale(scaleFits);
	background_mass_test->GetXaxis()->SetTitle("mass");
	background_mass_test->GetYaxis()->SetTitle("Sampled / Original");
	background_mass_test->GetYaxis()->SetTitleOffset(1.55);
	background_mass_test->Draw();
	background_mass_test->SetMinimum(minBkgPlot);
	background_mass_test->SetNdivisions(505, "X");

	bkgLine = new TLine( background_mass_test->GetBinCenter(1)-background_mass_test->GetBinWidth(1)/2., 1., background_mass_test->GetBinCenter(background_mass_test->GetNbinsX())+background_mass_test->GetBinWidth(background_mass_test->GetNbinsX())/2.,1.);
	bkgLine->SetLineWidth( 1 );
	bkgLine->SetLineStyle( 2 );
	bkgLine->SetLineColor( kBlack);
	bkgLine->Draw( "same" );

	sprintf(filename,"%s/fit_background_mass_test_%s.pdf",dirstruct,TreeBinID);
	c3->Print( filename );

	TH1* background_pT_test = (TH1*)results->Get("background_pT_test");
	background_pT_test->Scale(scaleFits);
	background_pT_test->GetXaxis()->SetTitle("p_{T}");
	background_pT_test->GetYaxis()->SetTitle("Sampled / Original");
	background_pT_test->GetYaxis()->SetTitleOffset(1.55);
	background_pT_test->Draw();
	background_pT_test->SetMinimum(minBkgPlot);

	bkgLine = new TLine( background_pT_test->GetBinCenter(1)-background_pT_test->GetBinWidth(1)/2., 1., background_pT_test->GetBinCenter(background_pT_test->GetNbinsX())+background_pT_test->GetBinWidth(background_pT_test->GetNbinsX())/2.,1.);
	bkgLine->SetLineWidth( 1 );
	bkgLine->SetLineStyle( 2 );
	bkgLine->SetLineColor( kBlack);
	bkgLine->Draw( "same" );

	sprintf(filename,"%s/fit_background_pT_test_%s.pdf",dirstruct,TreeBinID);
	c3->Print( filename );

	TH1* background_phiPX_test = (TH1*)results->Get("background_phiPX_test");
	background_phiPX_test->Scale(scaleFits);
	background_phiPX_test->GetXaxis()->SetTitle("#phi_{PX}");
	background_phiPX_test->GetYaxis()->SetTitle("Sampled / Original");
	background_phiPX_test->GetYaxis()->SetTitleOffset(1.55);
	background_phiPX_test->Draw();
	background_phiPX_test->SetMinimum(minBkgPlot);

	bkgLine = new TLine( background_phiPX_test->GetBinCenter(1)-background_phiPX_test->GetBinWidth(1)/2., 1., background_phiPX_test->GetBinCenter(background_phiPX_test->GetNbinsX())+background_phiPX_test->GetBinWidth(background_phiPX_test->GetNbinsX())/2.,1.);
	bkgLine->SetLineWidth( 1 );
	bkgLine->SetLineStyle( 2 );
	bkgLine->SetLineColor( kBlack);
	bkgLine->Draw( "same" );

	sprintf(filename,"%s/fit_background_phiPX_test_%s.pdf",dirstruct,TreeBinID);
	c3->Print( filename );

	TH1* background_costhPX_test = (TH1*)results->Get("background_costhPX_test");
	background_costhPX_test->Scale(scaleFits);
	background_costhPX_test->GetXaxis()->SetTitle("cos #theta_{PX}");
	background_costhPX_test->GetYaxis()->SetTitle("Sampled / Original");
	background_costhPX_test->GetYaxis()->SetTitleOffset(1.55);
	background_costhPX_test->Draw();
	background_costhPX_test->SetMinimum(minBkgPlot);

	bkgLine = new TLine( background_costhPX_test->GetBinCenter(1)-background_costhPX_test->GetBinWidth(1)/2., 1., background_costhPX_test->GetBinCenter(background_costhPX_test->GetNbinsX())+background_costhPX_test->GetBinWidth(background_costhPX_test->GetNbinsX())/2.,1.);
	bkgLine->SetLineWidth( 1 );
	bkgLine->SetLineStyle( 2 );
	bkgLine->SetLineColor( kBlack);
	bkgLine->Draw( "same" );

	sprintf(filename,"%s/fit_background_costhPX_test_%s.pdf",dirstruct,TreeBinID);
	c3->Print( filename );




	c3->SetRightMargin(0.03);

	double titleOffset=1.8;


	TH1* total_rap = (TH1*)results->Get("total_rap");
	total_rap->Scale(scaleFits);
	total_rap->GetXaxis()->SetTitle("|y|");
	total_rap->GetYaxis()->SetTitle("Entries [a.u.]");
	total_rap->GetYaxis()->SetTitleOffset(titleOffset);
	total_rap->Draw();
	total_rap->SetMinimum(minBkgPlot);
	sprintf(filename,"%s/fit_total_rap_%s.pdf",dirstruct,TreeBinID);
	c3->Print( filename );

	TH1* total_mass = (TH1*)results->Get("total_mass");
	total_mass->Scale(scaleFits);
	total_mass->GetXaxis()->SetTitle("mass");
	total_mass->GetYaxis()->SetTitle("Entries [a.u.]");
	total_mass->GetYaxis()->SetTitleOffset(titleOffset);
	total_mass->Draw();
	total_mass->SetMinimum(minBkgPlot);
	sprintf(filename,"%s/fit_total_mass_%s.pdf",dirstruct,TreeBinID);
	c3->Print( filename );

	TH1* total_pT = (TH1*)results->Get("total_pT");
	total_pT->Scale(scaleFits);
	total_pT->GetXaxis()->SetTitle("p_{T}");
	total_pT->GetYaxis()->SetTitle("Entries [a.u.]");
	total_pT->GetYaxis()->SetTitleOffset(titleOffset);
	total_pT->Draw();
	total_pT->SetMinimum(minBkgPlot);
	sprintf(filename,"%s/fit_total_pT_%s.pdf",dirstruct,TreeBinID);
	c3->Print( filename );

	TH1* total_phiPX = (TH1*)results->Get("total_phiPX");
	total_phiPX->Scale(scaleFits);
	total_phiPX->GetXaxis()->SetTitle("#phi_{PX}");
	total_phiPX->GetYaxis()->SetTitle("Entries [a.u.]");
	total_phiPX->GetYaxis()->SetTitleOffset(titleOffset);
	total_phiPX->Draw();
	total_phiPX->SetMinimum(minBkgPlot);
	sprintf(filename,"%s/fit_total_phiPX_%s.pdf",dirstruct,TreeBinID);
	c3->Print( filename );

	TH1* total_costhPX = (TH1*)results->Get("total_costhPX");
	total_costhPX->Scale(scaleFits);
	total_costhPX->GetXaxis()->SetTitle("cos #theta_{PX}");
	total_costhPX->GetYaxis()->SetTitle("Entries [a.u.]");
	total_costhPX->GetYaxis()->SetTitleOffset(titleOffset);
	total_costhPX->Draw();
	total_costhPX->SetMinimum(minBkgPlot);
	sprintf(filename,"%s/fit_total_costhPX_%s.pdf",dirstruct,TreeBinID);
	c3->Print( filename );

	bool DrawAccEff=false;
	if(DrawAccEff){
		TH1D* AccEff = (TH1D*)results->Get("AccEff");
		AccEff->GetXaxis()->SetTitle("Acceptance*#epsilon_{#mu#mu}");
		AccEff->GetYaxis()->SetTitleOffset(titleOffset);
		AccEff->Draw();
		sprintf(filename,"%s/AccEff_%s.pdf",dirstruct,TreeBinID);
		c3->Print( filename );

		TH1D* Acc = (TH1D*)results->Get("Acc");
		Acc->GetXaxis()->SetTitle("Acceptance");
		Acc->GetYaxis()->SetTitleOffset(titleOffset);
		Acc->Draw();
		sprintf(filename,"%s/Acc_%s.pdf",dirstruct,TreeBinID);
		c3->Print( filename );

		TH1D* Eff = (TH1D*)results->Get("Eff");
		Eff->GetXaxis()->SetTitle("#epsilon_{#mu#mu}");
		Eff->GetYaxis()->SetTitleOffset(titleOffset);
		Eff->Draw();
		sprintf(filename,"%s/DilepEfficiency_%s.pdf",dirstruct,TreeBinID);
		c3->Print( filename );
	}

	cout << endl << endl;

}
