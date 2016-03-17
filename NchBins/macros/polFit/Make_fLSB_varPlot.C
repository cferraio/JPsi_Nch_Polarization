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
#include "TH3.h"
#include "TGraphAsymmErrors.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLatex.h"
#include "TFrame.h"
#include "TLegend.h"
#include "TLine.h"

void Make_fLSB_varPlot(){

	  gROOT->Reset();
	  gROOT->SetBatch();

	  char savename[200];

	  char PlotDir[200];
	  sprintf(PlotDir,"May25_28p_46p_30p");
	  char PlotDirFull[200];
	  sprintf(PlotDirFull,"FigBuffer/fLSB_varPlots",PlotDir);
	  gSystem->mkdir(PlotDirFull);
	  sprintf(PlotDirFull,"FigBuffer/fLSB_varPlots/%s",PlotDir);
	  gSystem->mkdir(PlotDirFull);


	  double massMin=8.85;
	  double massMax=11.;
	  double fLSBMin=0.;
	  double fLSBMax=1.;

	  double kLin=1./(massMin-massMax);
	  double dLin=1.-kLin*massMin;

	  TF1 *fLSBlin = new TF1( "fLSBlin", "[0]*x[0]+[1]", massMin, massMax);
	  fLSBlin->SetParameter(0., kLin);
	  fLSBlin->SetParameter(1., dLin);

	  fLSBlin->SetLineColor(kRed);
	  fLSBlin->SetLineWidth(2.);
	  fLSBlin->SetLineStyle(1);

	  double mass1S;
	  double mass2S;
	  double mass3S;

	  double fLSBcentral1S=0.72;
	  double fLSBcentral2S=0.46;
	  double fLSBcentral3S=0.30;

	  double fLSBvar1S=0.28;
	  double fLSBvar2S=0.46;
	  double fLSBvar3S=0.30;

	  int nScan=1000;
	  for(int iScan=0;iScan<nScan+1;iScan++){
		  double massScan=massMin+(massMax-massMin)*iScan/nScan;
		  if(fLSBlin->Eval(massScan)<fLSBcentral1S) {mass1S=massScan; break;}
	  }
	  for(int iScan=0;iScan<nScan+1;iScan++){
		  double massScan=massMin+(massMax-massMin)*iScan/nScan;
		  if(fLSBlin->Eval(massScan)<fLSBcentral2S) {mass2S=massScan; break;}
	  }
	  for(int iScan=0;iScan<nScan+1;iScan++){
		  double massScan=massMin+(massMax-massMin)*iScan/nScan;
		  if(fLSBlin->Eval(massScan)<fLSBcentral3S) {mass3S=massScan; break;}
	  }

	  mass1S=(fLSBcentral1S-dLin)/kLin;
	  mass2S=(fLSBcentral2S-dLin)/kLin;
	  mass3S=(fLSBcentral3S-dLin)/kLin;

	  cout<<"mass1S "<<mass1S<<endl;
	  cout<<"mass2S "<<mass2S<<endl;
	  cout<<"mass3S "<<mass3S<<endl;

	  double mass1SvarMin=(fLSBcentral1S+fLSBvar1S-dLin)/kLin;
	  double mass1SvarMax=(fLSBcentral1S-fLSBvar1S-dLin)/kLin;
	  double mass2SvarMin=(fLSBcentral2S+fLSBvar2S-dLin)/kLin;
	  double mass2SvarMax=(fLSBcentral2S-fLSBvar2S-dLin)/kLin;
	  double mass3SvarMin=(fLSBcentral3S+fLSBvar3S-dLin)/kLin;
	  double mass3SvarMax=(fLSBcentral3S-fLSBvar3S-dLin)/kLin;


	  int NsColors[4]={0,418,600,616};

		TCanvas *fLSBCanvas = new TCanvas("fLSBCanvas","fLSBCanvas",1200,800);
	  	  gStyle->SetPalette(1);
	  	  gPad->SetFillColor(kWhite);

		fLSBCanvas->SetFillColor(kWhite);
//		fLSBCanvas->SetGrid();
		fLSBCanvas->GetFrame()->SetFillColor(kWhite);
		fLSBCanvas->GetFrame()->SetBorderSize(0);
		fLSBCanvas->SetRightMargin(0.05) ;
	  	fLSBCanvas->SetLeftMargin(0.1);
	    fLSBCanvas->SetTopMargin(0.05);
	    fLSBCanvas->SetBottomMargin(0.125);


		TLegend* plotLegend=new TLegend(0.65,0.775,0.95,0.95);
		plotLegend->SetFillColor(kWhite);
		plotLegend->SetTextFont(72);
		plotLegend->SetTextSize(0.04);
		plotLegend->SetBorderSize(1);
		char legendentry[200];

		double lineWidth=3;


		TH1F *fLSBHisto = new TH1F;
		fLSBHisto = fLSBCanvas->DrawFrame(massMin,fLSBMin,massMax,fLSBMax);
		fLSBHisto->SetXTitle("dimuon mass [GeV]");
		fLSBHisto->SetYTitle("f_{LSB}");
		fLSBHisto->GetYaxis()->SetTitleOffset(1.25);

		fLSBlin->Draw("same");


		TLine* Line1S = new TLine( mass1S, fLSBcentral1S-fLSBvar1S, mass1S ,fLSBcentral1S+fLSBvar1S);
		Line1S->SetLineWidth( 2 );
		Line1S->SetLineStyle( 1 );
		Line1S->SetLineColor( NsColors[1] );
		Line1S->Draw( "same" );

		TLine* Line2S = new TLine( mass2S, fLSBcentral2S-fLSBvar2S, mass2S ,fLSBcentral2S+fLSBvar2S);
		Line2S->SetLineWidth( 2 );
		Line2S->SetLineStyle( 1 );
		Line2S->SetLineColor( NsColors[2] );
		Line2S->Draw( "same" );

		TLine* Line3S = new TLine( mass3S, fLSBcentral3S-fLSBvar3S, mass3S ,fLSBcentral3S+fLSBvar3S);
		Line3S->SetLineWidth( 2 );
		Line3S->SetLineStyle( 1 );
		Line3S->SetLineColor( NsColors[3] );
		Line3S->Draw( "same" );

		TLine* LineRectangle1S = new TLine( mass1SvarMin, fLSBcentral1S+fLSBvar1S, mass1SvarMax ,fLSBcentral1S+fLSBvar1S);
		LineRectangle1S->SetLineWidth( 2 );
		LineRectangle1S->SetLineStyle( 2 );
		LineRectangle1S->SetLineColor( NsColors[1] );
		LineRectangle1S->Draw( "same" );

		LineRectangle1S = new TLine( mass1SvarMin, fLSBcentral1S-fLSBvar1S, mass1SvarMax ,fLSBcentral1S-fLSBvar1S);
		LineRectangle1S->SetLineWidth( 2 );
		LineRectangle1S->SetLineStyle( 2 );
		LineRectangle1S->SetLineColor( NsColors[1] );
		LineRectangle1S->Draw( "same" );

		LineRectangle1S = new TLine( mass1SvarMin, fLSBcentral1S-fLSBvar1S, mass1SvarMin ,fLSBcentral1S+fLSBvar1S);
		LineRectangle1S->SetLineWidth( 2 );
		LineRectangle1S->SetLineStyle( 2 );
		LineRectangle1S->SetLineColor( NsColors[1] );
		LineRectangle1S->Draw( "same" );

		LineRectangle1S = new TLine( mass1SvarMax, fLSBcentral1S-fLSBvar1S, mass1SvarMax ,fLSBcentral1S+fLSBvar1S);
		LineRectangle1S->SetLineWidth( 2 );
		LineRectangle1S->SetLineStyle( 2 );
		LineRectangle1S->SetLineColor( NsColors[1] );
		LineRectangle1S->Draw( "same" );



		TLine* LineRectangle2S = new TLine( mass2SvarMin, fLSBcentral2S+fLSBvar2S, mass2SvarMax ,fLSBcentral2S+fLSBvar2S);
		LineRectangle2S->SetLineWidth( 2 );
		LineRectangle2S->SetLineStyle( 2 );
		LineRectangle2S->SetLineColor( NsColors[2] );
		LineRectangle2S->Draw( "same" );

		LineRectangle2S = new TLine( mass2SvarMin, fLSBcentral2S-fLSBvar2S, mass2SvarMax ,fLSBcentral2S-fLSBvar2S);
		LineRectangle2S->SetLineWidth( 2 );
		LineRectangle2S->SetLineStyle( 2 );
		LineRectangle2S->SetLineColor( NsColors[2] );
		LineRectangle2S->Draw( "same" );

		LineRectangle2S = new TLine( mass2SvarMin, fLSBcentral2S-fLSBvar2S, mass2SvarMin ,fLSBcentral2S+fLSBvar2S);
		LineRectangle2S->SetLineWidth( 2 );
		LineRectangle2S->SetLineStyle( 2 );
		LineRectangle2S->SetLineColor( NsColors[2] );
		LineRectangle2S->Draw( "same" );

		LineRectangle2S = new TLine( mass2SvarMax, fLSBcentral2S-fLSBvar2S, mass2SvarMax ,fLSBcentral2S+fLSBvar2S);
		LineRectangle2S->SetLineWidth( 2 );
		LineRectangle2S->SetLineStyle( 2 );
		LineRectangle2S->SetLineColor( NsColors[2] );
		LineRectangle2S->Draw( "same" );



		TLine* LineRectangle3S = new TLine( mass3SvarMin, fLSBcentral3S+fLSBvar3S, mass3SvarMax ,fLSBcentral3S+fLSBvar3S);
		LineRectangle3S->SetLineWidth( 2 );
		LineRectangle3S->SetLineStyle( 2 );
		LineRectangle3S->SetLineColor( NsColors[3] );
		LineRectangle3S->Draw( "same" );

		LineRectangle3S = new TLine( mass3SvarMin, fLSBcentral3S-fLSBvar3S, mass3SvarMax ,fLSBcentral3S-fLSBvar3S);
		LineRectangle3S->SetLineWidth( 2 );
		LineRectangle3S->SetLineStyle( 2 );
		LineRectangle3S->SetLineColor( NsColors[3] );
		LineRectangle3S->Draw( "same" );

		LineRectangle3S = new TLine( mass3SvarMin, fLSBcentral3S-fLSBvar3S, mass3SvarMin ,fLSBcentral3S+fLSBvar3S);
		LineRectangle3S->SetLineWidth( 2 );
		LineRectangle3S->SetLineStyle( 2 );
		LineRectangle3S->SetLineColor( NsColors[3] );
		LineRectangle3S->Draw( "same" );

		LineRectangle3S = new TLine( mass3SvarMax, fLSBcentral3S-fLSBvar3S, mass3SvarMax ,fLSBcentral3S+fLSBvar3S);
		LineRectangle3S->SetLineWidth( 2 );
		LineRectangle3S->SetLineStyle( 2 );
		LineRectangle3S->SetLineColor( NsColors[3] );
		LineRectangle3S->Draw( "same" );

		sprintf(legendentry,"Linear Assumption");
		plotLegend->AddEntry(fLSBlin,legendentry,"l");
		sprintf(legendentry,"1S variation");
		plotLegend->AddEntry(Line1S,legendentry,"l");
		sprintf(legendentry,"2S variation");
		plotLegend->AddEntry(Line2S,legendentry,"l");
		sprintf(legendentry,"3S variation");
		plotLegend->AddEntry(Line3S,legendentry,"l");

		plotLegend->Draw();

  	  sprintf(savename,"%s/fLSB_vs_mass.pdf",PlotDirFull);
  	  fLSBCanvas->SaveAs(savename);




}
