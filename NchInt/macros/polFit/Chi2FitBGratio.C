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
#include "TLegend.h"
#include "TFrame.h"
#include "TLine.h"

void Chi2FitBGratio(){

	  gROOT->Reset();
	  gROOT->SetBatch();


////// Settings /////////

	  char DataID[200];
	  sprintf(DataID,"SetOfCuts11_HighCtau3_BGratioTest");
	  int nSigma=1;

	  char PlotDir[200];
	  sprintf(PlotDir,"Chi2_May21_2D_AbsCosthMax1_Rap0Pt0_1sigma");

	  int FitDimension=1;
	  double AbsCosthMax=1.;
	  double RelativeErrorMax=100.;

	  int rapMin=0;
	  int rapMax=0;
	  int ptMin=0;
	  int ptMax=0;
	  int frameMin=1;
	  int frameMax=3;
	  int stateMin=1;
	  int stateMax=3;

	  const int nFracL=201;
	  double fracLmin=0.;
	  double fracLmax=1.;

	  bool plotStuff=false;
	  bool StatesInd=true;
	  bool SetModelErrorZero=false;

/// Plotting settings:

	  double yMax=2.5;
	  double yMin=0.5;
	  double minLineWidth=1.;
	  int minLineStyle=2;

/////////////////////////

	  char PlotDirFull[200];
	  sprintf(PlotDirFull,"FigBuffer/BGratioFits/%s",PlotDir);
	  gSystem->mkdir(PlotDirFull);
	  char SystDirFull[200];
	  sprintf(SystDirFull,"Systematics/TheGreatRun_BKGratio/%s",PlotDir);
	  gSystem->mkdir(SystDirFull);

	  char frameIdentity[50];
	    char text[200];

	  const int nStatesact=stateMax-stateMin;
	  const int npTact = ptMax-ptMin;//11
	  const int nrapact = rapMax-rapMin;//3
	  const int nframeact = frameMax-frameMin;//4

	  const int nStates=3;
	  const int npT = 11;//11
	  const int nrap = 3;//3
	  const int nframe = 4;//4

	  double fracLvar[nFracL][npT][nrap][nframe][nStates];
	  double FitConst[nFracL][npT][nrap][nframe][nStates];
	  double err_FitConst[nFracL][npT][nrap][nframe][nStates];
	  double FitConst_Chi2[nFracL][npT][nrap][nframe][nStates];
	  double FitConst_NDF[nFracL][npT][nrap][nframe][nStates];


	  for( int nState  = stateMin; nState  < stateMax+1; nState++ ){



		  TFile* inFile;
		  TFile* inFile2;



  TF1 *fcosthphi1D = new TF1( "fcosthphi1D", "[0]", -AbsCosthMax, AbsCosthMax);
  fcosthphi1D->SetParameter(0,1.);

  TF2 *fcosthphi2D = new TF2( "fcosthphi2D", "[0]+(0*x[0]+0*x[1])", -AbsCosthMax, AbsCosthMax, -180., 180. );
  fcosthphi2D->SetParameter(0,1.);


  double  lth[npT][nrap][nframe],   dlth[npT][nrap][nframe];
  double  lph[npT][nrap][nframe],   dlph[npT][nrap][nframe];
  double  lthph[npT][nrap][nframe], dlthph[npT][nrap][nframe];

  char histoName[50];

  TH2D* histoBGmodelL;
  TH2D* histoBGmodelR;
  TH2D* histoBGmodel;
  TH2D* histoBGdata;
  TH2D* histoBGratio;

  TH1D* histoBGmodel1D;
  TH1D* histoBGdata1D;
  TH1D* histoBGratio1D;


  TCanvas *c3;


  for( int iframe  = frameMin; iframe  < frameMax+1; iframe++ ){

	  if(iframe==1) sprintf(frameIdentity,"CS");
	  if(iframe==2) sprintf(frameIdentity,"HX");
	  if(iframe==3) sprintf(frameIdentity,"PX");

  for( int irap  = rapMin; irap  < rapMax+1; irap++ ){
  for( int ipT   = ptMin; ipT   < ptMax+1;  ipT++  ){

	  char inFileName[500];
	  sprintf(inFileName,"/afs/hephy.at/scratch/k/knuenz/CMSSW_4_2_4_patch2/src/UpsilonPol/macros/DataFiles/%s/AllStates_%d.00Sigma_FracLSB50Percent/data_%dSUps_rap%d_pT%d.root",DataID,nSigma,nState,irap,ipT);
//	  sprintf(inFileName,"/afs/hephy.at/scratch/k/knuenz/CMSSW_4_2_4_patch2/src/UpsilonPol/macros/DataFiles/SetOfCuts11_May20Centrals/AllStates_%d.00Sigma_FracLSB50Percent/data_%dSUps_rap%d_pT%d.root",nSigma,nState,irap,ipT);
	    inFile = new TFile(inFileName,"READ");
	  sprintf(inFileName,"/afs/hephy.at/scratch/k/knuenz/CMSSW_4_2_4_patch2/src/UpsilonPol/macros/DataFiles/%s/AllStates_%d.00Sigma_FracLSB50Percent/Figures/Ups%dS/AngDistHist.root",DataID,nSigma,nState);
	  inFile2 = new TFile(inFileName,"READ");



  if(iframe==1)    sprintf ( histoName, "CosthPhi_CS_rap%d_pT%d", irap, ipT );
  if(iframe==2)    sprintf ( histoName, "CosthPhi_HX_rap%d_pT%d", irap, ipT );
  if(iframe==3)    sprintf ( histoName, "CosthPhi_PHX_rap%d_pT%d", irap, ipT );

    histoBGdata  = (TH2D*)inFile2->Get( histoName );
    histoBGdata->Sumw2();


    if(iframe==1)    sprintf ( histoName, "hBG_cosThetaPhi_CS_L");
    if(iframe==2)    sprintf ( histoName, "hBG_cosThetaPhi_HX_L");
    if(iframe==3)    sprintf ( histoName, "hBG_cosThetaPhi_PHX_L");

  cout << histoName << endl;

  histoBGmodelL  = (TH2D*)inFile->Get( histoName );
  histoBGmodelL->Sumw2();

    if(iframe==1)    sprintf ( histoName, "hBG_cosThetaPhi_CS_R");
    if(iframe==2)    sprintf ( histoName, "hBG_cosThetaPhi_HX_R");
    if(iframe==3)    sprintf ( histoName, "hBG_cosThetaPhi_PHX_R");

  cout << histoName << endl;

  histoBGmodelR  = (TH2D*)inFile->Get( histoName );
  histoBGmodelR->Sumw2();


  for(int iFracL=0;iFracL<nFracL;iFracL++){

	  double fracL=fracLmin+iFracL*fracLmax/double(nFracL-1);

  TH2D* histoBGmodelRclone  = (TH2D*)histoBGmodelR->Clone( "histoBGmodelRclone" );
  TH2D* histoBGmodelLclone  = (TH2D*)histoBGmodelL->Clone( "histoBGmodelLclone" );

  histoBGmodelRclone->Sumw2();
  histoBGmodelLclone->Sumw2();

  if(iFracL>0) histoBGmodel->Reset();
  histoBGmodelLclone->Scale(fracL/histoBGmodelLclone->Integral());
  histoBGmodelRclone->Scale((1.-fracL)/histoBGmodelRclone->Integral());
  histoBGmodel = (TH2D *) histoBGmodelLclone->Clone("histoBGmodel");
  histoBGmodel->Add(histoBGmodelRclone);

  histoBGmodel->Sumw2();

////////// Cut on costh //////////////////////


    histoBGmodel1D=histoBGmodel->ProjectionX("histoBGmodel1D", 0, -1,"e");
    histoBGmodel1D->Sumw2();

    if(SetModelErrorZero){
	    for(int iX=0;iX<histoBGmodel1D->GetNbinsX()+1;iX++){
	        	int globalBin=histoBGmodel1D->GetBin(iX);
	        	histoBGmodel1D->SetBinError(iX,0.);
	        }
    }

	histoBGdata1D=histoBGdata->ProjectionX("histoBGdata1D", 0, -1,"e");

    histoBGratio1D=(TH1D*)histoBGdata1D->Clone("histoBGratio1D");
    histoBGratio1D->Divide(histoBGmodel1D);
    histoBGratio1D->Sumw2();

    int numBinsWithContent=0;
	    for(int iX=0;iX<histoBGratio1D->GetNbinsX()+1;iX++){
	        	int globalBin=histoBGratio1D->GetBin(iX);
	        	if(TMath::Abs(histoBGratio1D->GetBinCenter(globalBin))>AbsCosthMax||histoBGratio1D->GetBinError(globalBin)/histoBGratio1D->GetBinContent(globalBin)>RelativeErrorMax) {
	        		histoBGratio1D->SetBinContent(iX,0.);
	        		histoBGratio1D->SetBinError(iX,0.);
	        	}
	        	else numBinsWithContent++;
	        }


	    histoBGratio=(TH2D*)histoBGdata->Clone("histoBGratio");
	    histoBGratio->Divide(histoBGmodel);

	int numBinsWithContent2D=0;
	for(int iX=0;iX<histoBGratio->GetNbinsX()+1;iX++){
	    for(int iY=0;iY<histoBGratio->GetNbinsY()+1;iY++){
	    	int globalBin=histoBGratio->GetBin(iX);
	    	if(TMath::Abs(histoBGratio->GetBinCenter(globalBin))>AbsCosthMax||histoBGratio->GetBinError(globalBin)/histoBGratio->GetBinContent(globalBin)>RelativeErrorMax) {
	    		histoBGratio->SetBinContent(iX,iY,0.);
	    		histoBGratio->SetBinError(iX,iY,0.);
	    	}
        	else numBinsWithContent2D++;
	    }
	}


//////////////////////////////////////////////////////////////////////////


	    if(FitDimension==1){
		    histoBGratio1D->Scale(double(numBinsWithContent)/histoBGratio1D->Integral());
		    histoBGratio1D->Sumw2();
		    histoBGratio1D->Print();
		    histoBGratio1D->Fit("fcosthphi1D","0");

		      fracLvar[iFracL][ipT][irap][iframe][nState] = fracL;
		      FitConst[iFracL][ipT][irap][iframe][nState] = fcosthphi1D->GetParameter(0);
		      err_FitConst[iFracL][ipT][irap][iframe][nState] = fcosthphi1D->GetParError(0);
		      FitConst_Chi2[iFracL][ipT][irap][iframe][nState] = fcosthphi1D->GetChisquare();
		      FitConst_NDF[iFracL][ipT][irap][iframe][nState] = fcosthphi1D->GetNDF();
	    }
	    if(FitDimension==2){
			histoBGratio->Sumw2();
			histoBGratio->Scale(double(numBinsWithContent2D)/histoBGratio->Integral());
		    histoBGratio->Print();
		    histoBGratio->Fit("fcosthphi2D","0");

		      fracLvar[iFracL][ipT][irap][iframe][nState] = fracL;
		      FitConst[iFracL][ipT][irap][iframe][nState] = fcosthphi2D->GetParameter(0);
		      err_FitConst[iFracL][ipT][irap][iframe][nState] = fcosthphi2D->GetParError(0);
		      FitConst_Chi2[iFracL][ipT][irap][iframe][nState] = fcosthphi2D->GetChisquare();
		      FitConst_NDF[iFracL][ipT][irap][iframe][nState] = fcosthphi2D->GetNDF();
	    }


    	      cout<<"FitConst "<<FitConst[iFracL][ipT][irap][iframe][nState]<<endl;
    	      cout<<"err_FitConst "<<err_FitConst[iFracL][ipT][irap][iframe][nState]<<endl;
    	      cout<<"FitConst_Chi2 "<<FitConst_Chi2[iFracL][ipT][irap][iframe][nState]<<endl;
    	      cout<<"FitConst_NDF "<<FitConst_NDF[iFracL][ipT][irap][iframe][nState]<<endl;


    		gStyle->SetPalette(1);

    		if(plotStuff){
    		TCanvas *c2 = new TCanvas("c2","c2",1200,800);
    		gPad->SetFillColor(kWhite);

  		    histoBGratio1D->SetMarkerStyle(20);
  		    histoBGratio1D->SetTitle(0);
  		    histoBGratio1D->SetStats(0);
  		    histoBGratio1D->Draw();
 		    fcosthphi1D->SetLineColor(kGreen+2);
  		    fcosthphi1D->SetLineWidth(2.);
  		    fcosthphi1D->SetLineStyle(1);
  		    fcosthphi1D->Draw("same");

  		    char text[200];
  		    sprintf(text,"Const = %1.3f #pm %1.3f",FitConst[iFracL][ipT][irap][iframe][nState],err_FitConst[iFracL][ipT][irap][iframe][nState]);
  		    TLatex *fitBGratio1Dtex = new TLatex(-0.8,histoBGratio1D->GetMaximum()*0.2,text);
  		    fitBGratio1Dtex->SetTextSize(0.05)                                                                                                                                                                                                                                             ;
  		    fitBGratio1Dtex->Draw( "same" )                                                                                                                                                                                                                                                 ;

    		char savenamePlot[200];
  		    sprintf(savenamePlot,"FigBuffer/BGratioFits/%s/histoBGratioWithFit_%s_%dSUps_rap%d_pT%d_nFracL%d.pdf",PlotDir,frameIdentity,nState,irap,ipT,iFracL);
  		    c2->SaveAs(savenamePlot);
    		}


  }//FracLloop

    inFile->Close();
    inFile2->Close();

  }//ptLoop
  }//rapLoop
  }//FrameLoop
  }//nState LOOP


		double FracLCenter[nFracL];
		double err_FracLCenter[nFracL];
		double Chi2Func1S_CS[nFracL];
		double err_Chi2Func1S_CS[nFracL];
		double Chi2Func1S_HX[nFracL];
		double err_Chi2Func1S_HX[nFracL];
		double Chi2Func1S_PX[nFracL];
		double err_Chi2Func1S_PX[nFracL];

		double Chi2Func2S_CS[nFracL];
		double err_Chi2Func2S_CS[nFracL];
		double Chi2Func2S_HX[nFracL];
		double err_Chi2Func2S_HX[nFracL];
		double Chi2Func2S_PX[nFracL];
		double err_Chi2Func2S_PX[nFracL];

		double Chi2Func3S_CS[nFracL];
		double err_Chi2Func3S_CS[nFracL];
		double Chi2Func3S_HX[nFracL];
		double err_Chi2Func3S_HX[nFracL];
		double Chi2Func3S_PX[nFracL];
		double err_Chi2Func3S_PX[nFracL];

		double Chi2Min1S_CS[nFracL];

		double Chi2Func_nSnFrame[nFracL];
		int NDF_nSnFrame[nFracL];

		double minChi2FracL[nframe][nStates+1]={{1000,1000,1000,1000},{1000,1000,1000,1000},{1000,1000,1000,1000},{1000,1000,1000,1000}};
		double Chi2Buffer[nframe][nStates+1]={{1000,1000,1000,1000},{1000,1000,1000,1000},{1000,1000,1000,1000},{1000,1000,1000,1000}};
		double Chi2Min[nframe][nStates+1];
		int NDFMin[nframe][nStates+1];

	for(int iFracL=0;iFracL<nFracL;iFracL++){
		  for( int nState  = stateMin; nState  < stateMax+1; nState++ ){
		  for( int iframe  = frameMin; iframe  < frameMax+1; iframe++ ){


			  FracLCenter[iFracL]=fracLvar[iFracL][ptMin][rapMin][frameMin][stateMin];
			  err_FracLCenter[iFracL]=0;

			  Chi2Func_nSnFrame[iFracL]=0;
			  NDF_nSnFrame[iFracL]=0;

			  for( int irap  = rapMin; irap  < rapMax+1; irap++ ){
				  for( int ipT   = ptMin; ipT   < ptMax+1;  ipT++  ){

						  err_Chi2Func1S_CS[iFracL]=0;
						  err_Chi2Func1S_HX[iFracL]=0;
						  err_Chi2Func1S_PX[iFracL]=0;
						  err_Chi2Func2S_CS[iFracL]=0;
						  err_Chi2Func2S_HX[iFracL]=0;
						  err_Chi2Func2S_PX[iFracL]=0;
						  err_Chi2Func3S_CS[iFracL]=0;
						  err_Chi2Func3S_HX[iFracL]=0;
						  err_Chi2Func3S_PX[iFracL]=0;

						  Chi2Func_nSnFrame[iFracL]+=FitConst_Chi2[iFracL][ipT][irap][iframe][nState];
						  NDF_nSnFrame[iFracL]+=FitConst_NDF[iFracL][ipT][irap][iframe][nState];

					  }//ptLoop
					  }//rapLoop

			  if(iframe==1&&nState==1) Chi2Func1S_CS[iFracL]=Chi2Func_nSnFrame[iFracL]/double(NDF_nSnFrame[iFracL]);
			  if(iframe==2&&nState==1) Chi2Func1S_HX[iFracL]=Chi2Func_nSnFrame[iFracL]/double(NDF_nSnFrame[iFracL]);
			  if(iframe==3&&nState==1) Chi2Func1S_PX[iFracL]=Chi2Func_nSnFrame[iFracL]/double(NDF_nSnFrame[iFracL]);
			  if(iframe==1&&nState==2) Chi2Func2S_CS[iFracL]=Chi2Func_nSnFrame[iFracL]/double(NDF_nSnFrame[iFracL]);
			  if(iframe==2&&nState==2) Chi2Func2S_HX[iFracL]=Chi2Func_nSnFrame[iFracL]/double(NDF_nSnFrame[iFracL]);
			  if(iframe==3&&nState==2) Chi2Func2S_PX[iFracL]=Chi2Func_nSnFrame[iFracL]/double(NDF_nSnFrame[iFracL]);
			  if(iframe==1&&nState==3) Chi2Func3S_CS[iFracL]=Chi2Func_nSnFrame[iFracL]/double(NDF_nSnFrame[iFracL]);
			  if(iframe==2&&nState==3) Chi2Func3S_HX[iFracL]=Chi2Func_nSnFrame[iFracL]/double(NDF_nSnFrame[iFracL]);
			  if(iframe==3&&nState==3) Chi2Func3S_PX[iFracL]=Chi2Func_nSnFrame[iFracL]/double(NDF_nSnFrame[iFracL]);

			  if(Chi2Func_nSnFrame[iFracL]/double(NDF_nSnFrame[iFracL])<Chi2Buffer[iframe][nState]) {Chi2Buffer[iframe][nState]=Chi2Func_nSnFrame[iFracL]/double(NDF_nSnFrame[iFracL]); minChi2FracL[iframe][nState]=FracLCenter[iFracL]; Chi2Min[iframe][nState]=Chi2Func_nSnFrame[iFracL]; NDFMin[iframe][nState]=NDF_nSnFrame[iFracL];}

					  }//FrameLoop
					  }//nState LOOP
					  }//FracLloop

		TGraphAsymmErrors* Chi2Func1Sgraph_CS;
		TGraphAsymmErrors* Chi2Func2Sgraph_CS;
		TGraphAsymmErrors* Chi2Func3Sgraph_CS;

		Chi2Func1Sgraph_CS = new TGraphAsymmErrors(nFracL,FracLCenter,Chi2Func1S_CS,err_FracLCenter,err_FracLCenter,err_Chi2Func1S_CS,err_Chi2Func1S_CS);
		Chi2Func2Sgraph_CS = new TGraphAsymmErrors(nFracL,FracLCenter,Chi2Func2S_CS,err_FracLCenter,err_FracLCenter,err_Chi2Func2S_CS,err_Chi2Func2S_CS);
		Chi2Func3Sgraph_CS = new TGraphAsymmErrors(nFracL,FracLCenter,Chi2Func3S_CS,err_FracLCenter,err_FracLCenter,err_Chi2Func3S_CS,err_Chi2Func3S_CS);

		TGraphAsymmErrors* Chi2Func1Sgraph_HX;
		TGraphAsymmErrors* Chi2Func2Sgraph_HX;
		TGraphAsymmErrors* Chi2Func3Sgraph_HX;

		Chi2Func1Sgraph_HX = new TGraphAsymmErrors(nFracL,FracLCenter,Chi2Func1S_HX,err_FracLCenter,err_FracLCenter,err_Chi2Func1S_HX,err_Chi2Func1S_HX);
		Chi2Func2Sgraph_HX = new TGraphAsymmErrors(nFracL,FracLCenter,Chi2Func2S_HX,err_FracLCenter,err_FracLCenter,err_Chi2Func2S_HX,err_Chi2Func2S_HX);
		Chi2Func3Sgraph_HX = new TGraphAsymmErrors(nFracL,FracLCenter,Chi2Func3S_HX,err_FracLCenter,err_FracLCenter,err_Chi2Func3S_HX,err_Chi2Func3S_HX);

		TGraphAsymmErrors* Chi2Func1Sgraph_PX;
		TGraphAsymmErrors* Chi2Func2Sgraph_PX;
		TGraphAsymmErrors* Chi2Func3Sgraph_PX;

		Chi2Func1Sgraph_PX = new TGraphAsymmErrors(nFracL,FracLCenter,Chi2Func1S_PX,err_FracLCenter,err_FracLCenter,err_Chi2Func1S_PX,err_Chi2Func1S_PX);
		Chi2Func2Sgraph_PX = new TGraphAsymmErrors(nFracL,FracLCenter,Chi2Func2S_PX,err_FracLCenter,err_FracLCenter,err_Chi2Func2S_PX,err_Chi2Func2S_PX);
		Chi2Func3Sgraph_PX = new TGraphAsymmErrors(nFracL,FracLCenter,Chi2Func3S_PX,err_FracLCenter,err_FracLCenter,err_Chi2Func3S_PX,err_Chi2Func3S_PX);



		   gStyle->SetPalette(1,0);
		   gStyle->SetPadBottomMargin(0.12);
		   gStyle->SetPadLeftMargin(0.13);
		   gStyle->SetPadRightMargin(0.15);
		   gStyle->SetPadTopMargin(0.05);

		   gStyle->SetTickLength(-0.02, "xyz");
		   gStyle->SetLabelOffset(0.02, "x");
		   gStyle->SetLabelOffset(0.02, "y");
		   gStyle->SetTitleOffset(1.3, "x");
		   gStyle->SetTitleOffset(1.4, "y");
		   gStyle->SetTitleFillColor(kWhite);

		TCanvas *plotCanvas = new TCanvas("plotCanvas","plotCanvas",1000,800);

		plotCanvas->SetFillColor(kWhite);
		plotCanvas->GetFrame()->SetFillColor(kWhite);
		plotCanvas->GetFrame()->SetBorderSize(0);
		plotCanvas->SetRightMargin(0.05) ;


		TLegend* plotLegend;
		TH1F *plotHisto;

		int fLSBLineColor[nframe][nStates+1]={{0,632,600,418},{0,632,600,418},{0,632,600,418},{0,632,600,418}};
		int fLSBMarkerColor[nframe][nStates+1]={{0,632,600,418},{0,632,600,418},{0,632,600,418},{0,632,600,418}};
		double fLSBMarkerSize=1.;
		int fLSBMarkerStyle=1;

		Chi2Func1Sgraph_CS->SetMarkerColor(fLSBMarkerColor[1][1]);
		Chi2Func1Sgraph_CS->SetMarkerSize(fLSBMarkerSize);
		Chi2Func1Sgraph_CS->SetMarkerStyle(fLSBMarkerStyle);
		Chi2Func1Sgraph_CS->SetLineColor(fLSBLineColor[1][1]);

		Chi2Func1Sgraph_HX->SetMarkerColor(fLSBMarkerColor[2][1]);
		Chi2Func1Sgraph_HX->SetMarkerSize(fLSBMarkerSize);
		Chi2Func1Sgraph_HX->SetMarkerStyle(fLSBMarkerStyle);
		Chi2Func1Sgraph_HX->SetLineColor(fLSBLineColor[2][1]);

		Chi2Func1Sgraph_PX->SetMarkerColor(fLSBMarkerColor[3][1]);
		Chi2Func1Sgraph_PX->SetMarkerSize(fLSBMarkerSize);
		Chi2Func1Sgraph_PX->SetMarkerStyle(fLSBMarkerStyle);
		Chi2Func1Sgraph_PX->SetLineColor(fLSBLineColor[3][1]);

		Chi2Func2Sgraph_CS->SetMarkerColor(fLSBMarkerColor[1][2]);
		Chi2Func2Sgraph_CS->SetMarkerSize(fLSBMarkerSize);
		Chi2Func2Sgraph_CS->SetMarkerStyle(fLSBMarkerStyle);
		Chi2Func2Sgraph_CS->SetLineColor(fLSBLineColor[1][2]);

		Chi2Func2Sgraph_HX->SetMarkerColor(fLSBMarkerColor[2][2]);
		Chi2Func2Sgraph_HX->SetMarkerSize(fLSBMarkerSize);
		Chi2Func2Sgraph_HX->SetMarkerStyle(fLSBMarkerStyle);
		Chi2Func2Sgraph_HX->SetLineColor(fLSBLineColor[2][2]);

		Chi2Func2Sgraph_PX->SetMarkerColor(fLSBMarkerColor[3][2]);
		Chi2Func2Sgraph_PX->SetMarkerSize(fLSBMarkerSize);
		Chi2Func2Sgraph_PX->SetMarkerStyle(fLSBMarkerStyle);
		Chi2Func2Sgraph_PX->SetLineColor(fLSBLineColor[3][2]);

		Chi2Func3Sgraph_CS->SetMarkerColor(fLSBMarkerColor[1][3]);
		Chi2Func3Sgraph_CS->SetMarkerSize(fLSBMarkerSize);
		Chi2Func3Sgraph_CS->SetMarkerStyle(fLSBMarkerStyle);
		Chi2Func3Sgraph_CS->SetLineColor(fLSBLineColor[1][3]);

		Chi2Func3Sgraph_HX->SetMarkerColor(fLSBMarkerColor[2][3]);
		Chi2Func3Sgraph_HX->SetMarkerSize(fLSBMarkerSize);
		Chi2Func3Sgraph_HX->SetMarkerStyle(fLSBMarkerStyle);
		Chi2Func3Sgraph_HX->SetLineColor(fLSBLineColor[2][3]);

		Chi2Func3Sgraph_PX->SetMarkerColor(fLSBMarkerColor[3][3]);
		Chi2Func3Sgraph_PX->SetMarkerSize(fLSBMarkerSize);
		Chi2Func3Sgraph_PX->SetMarkerStyle(fLSBMarkerStyle);
		Chi2Func3Sgraph_PX->SetLineColor(fLSBLineColor[3][3]);

		char drawGraphStyle[200];
		sprintf(drawGraphStyle,"l");

		char savenamePlot[200];
		double FontSizeChi2=0.03;
		int nStateBuff;

		  for( int iframe  = frameMin; iframe  < frameMax+1; iframe++ ){

			  if(iframe==1) sprintf(frameIdentity,"CS");
			  if(iframe==2) sprintf(frameIdentity,"HX");
			  if(iframe==3) sprintf(frameIdentity,"PX");

				plotLegend=new TLegend(0.55,0.75,0.95,0.95);
				plotLegend->SetFillColor(0);
				plotLegend->SetTextFont(72);
				plotLegend->SetTextSize(0.04);
				plotLegend->SetBorderSize(1);
				char legendentry[200];


				plotHisto = new TH1F;
				plotHisto = plotCanvas->DrawFrame(fracLmin,yMin,fracLmax,yMax);
				plotHisto->SetXTitle("f_{LSB}");
				plotHisto->SetYTitle("#chi^{2} / ndf");
				plotHisto->GetYaxis()->SetTitleOffset(1.5);

				if(iframe==1){
			    Chi2Func1Sgraph_CS->Draw(drawGraphStyle);
				if(StatesInd){
				Chi2Func2Sgraph_CS->Draw(drawGraphStyle);
			    Chi2Func3Sgraph_CS->Draw(drawGraphStyle);
				}
				sprintf(legendentry,"#Upsilon(1S) mass region");
				plotLegend->AddEntry(Chi2Func1Sgraph_CS,legendentry,"l");
				sprintf(legendentry,"#Upsilon(2S) mass region");
				plotLegend->AddEntry(Chi2Func2Sgraph_CS,legendentry,"l");
				sprintf(legendentry,"#Upsilon(3S) mass region");
				plotLegend->AddEntry(Chi2Func3Sgraph_CS,legendentry,"l");


				}
				if(iframe==2){
			    Chi2Func1Sgraph_HX->Draw(drawGraphStyle);
				if(StatesInd){
			    Chi2Func2Sgraph_HX->Draw(drawGraphStyle);
			    Chi2Func3Sgraph_HX->Draw(drawGraphStyle);
				}
				sprintf(legendentry,"#Upsilon(1S) mass region");
				plotLegend->AddEntry(Chi2Func1Sgraph_HX,legendentry,"l");
				sprintf(legendentry,"#Upsilon(2S) mass region");
				plotLegend->AddEntry(Chi2Func2Sgraph_HX,legendentry,"l");
				sprintf(legendentry,"#Upsilon(3S) mass region");
				plotLegend->AddEntry(Chi2Func3Sgraph_HX,legendentry,"l");
				}
				if(iframe==3){
			    Chi2Func1Sgraph_PX->Draw(drawGraphStyle);
				if(StatesInd){
			    Chi2Func2Sgraph_PX->Draw(drawGraphStyle);
			    Chi2Func3Sgraph_PX->Draw(drawGraphStyle);
				}



				sprintf(legendentry,"#Upsilon(1S) mass region");
				plotLegend->AddEntry(Chi2Func1Sgraph_PX,legendentry,"l");
				sprintf(legendentry,"#Upsilon(2S) mass region");
				plotLegend->AddEntry(Chi2Func2Sgraph_PX,legendentry,"l");
				sprintf(legendentry,"#Upsilon(3S) mass region");
				plotLegend->AddEntry(Chi2Func3Sgraph_PX,legendentry,"l");
				}


				if(StatesInd){
				nStateBuff=1;
			    TLine* Chi2Min1S = new TLine( minChi2FracL[iframe][nStateBuff], yMin , minChi2FracL[iframe][nStateBuff] ,yMax);
			    Chi2Min1S->SetLineWidth( minLineWidth );
			    Chi2Min1S->SetLineStyle( minLineStyle );
			    Chi2Min1S->SetLineColor( fLSBLineColor[iframe][nStateBuff] );
			    Chi2Min1S->Draw( "same" );
			    sprintf(text,"#color[2]{#Upsilon(1S) minimum:} f_{LSB} = %1.2f, #chi^{2} / ndf = %1.3f / %d = %1.3f,  #chi^{2}-prob. = %1.3f", minChi2FracL[iframe][nStateBuff], Chi2Min[iframe][nStateBuff], NDFMin[iframe][nStateBuff], Chi2Min[iframe][nStateBuff]/double(NDFMin[iframe][nStateBuff]), TMath::Prob(Chi2Min[iframe][nStateBuff], NDFMin[iframe][nStateBuff]));
			    TLatex *textChi2_1S = new TLatex(0.02,yMin+(yMax-yMin)*0.125,text);
			    textChi2_1S->SetTextSize(FontSizeChi2)                                                                                                                                                                                                                                             ;
			    textChi2_1S->Draw( "same" )                                                                                                                                                                                                                                                 ;

			    nStateBuff=2;
			    TLine* Chi2Min2S = new TLine( minChi2FracL[iframe][nStateBuff], yMin, minChi2FracL[iframe][nStateBuff] ,yMax);
			    Chi2Min2S->SetLineWidth( minLineWidth );
			    Chi2Min2S->SetLineStyle( minLineStyle );
			    Chi2Min2S->SetLineColor( fLSBLineColor[iframe][nStateBuff] );
			    Chi2Min2S->Draw( "same" );
			    sprintf(text,"#color[4]{#Upsilon(2S) minimum:} f_{LSB} = %1.2f, #chi^{2} / ndf = %1.3f / %d = %1.3f,  #chi^{2}-prob. = %1.3f", minChi2FracL[iframe][nStateBuff], Chi2Min[iframe][nStateBuff], NDFMin[iframe][nStateBuff], Chi2Min[iframe][nStateBuff]/double(NDFMin[iframe][nStateBuff]), TMath::Prob(Chi2Min[iframe][nStateBuff], NDFMin[iframe][nStateBuff]));
			    TLatex *textChi2_2S = new TLatex(0.02,yMin+(yMax-yMin)*0.075,text);
			    textChi2_2S->SetTextSize(FontSizeChi2)                                                                                                                                                                                                                                             ;
			    textChi2_2S->Draw( "same" )                                                                                                                                                                                                                                                 ;

			    nStateBuff=3;
			    TLine* Chi2Min3S = new TLine( minChi2FracL[iframe][nStateBuff], yMin, minChi2FracL[iframe][nStateBuff] ,yMax);
			    Chi2Min3S->SetLineWidth( minLineWidth );
			    Chi2Min3S->SetLineStyle( minLineStyle );
			    Chi2Min3S->SetLineColor( fLSBLineColor[iframe][nStateBuff] );
			    Chi2Min3S->Draw( "same" );
			    sprintf(text,"#color[3]{#Upsilon(3S) minimum:} f_{LSB} = %1.2f, #chi^{2} / ndf = %1.3f / %d = %1.3f,  #chi^{2}-prob. = %1.3f", minChi2FracL[iframe][nStateBuff], Chi2Min[iframe][nStateBuff], NDFMin[iframe][nStateBuff], Chi2Min[iframe][nStateBuff]/double(NDFMin[iframe][nStateBuff]), TMath::Prob(Chi2Min[iframe][nStateBuff], NDFMin[iframe][nStateBuff]));
			    TLatex *textChi2_3S = new TLatex(0.02,yMin+(yMax-yMin)*0.025,text);
			    textChi2_3S->SetTextSize(FontSizeChi2)                                                                                                                                                                                                                                             ;
			    textChi2_3S->Draw( "same" )                                                                                                                                                                                                                                                 ;
				plotLegend->Draw();
				}
				if(!StatesInd){
					nStateBuff=1;
				    TLine* Chi2Min1S = new TLine( minChi2FracL[iframe][nStateBuff], yMin, minChi2FracL[iframe][nStateBuff] ,yMax);
				    Chi2Min1S->SetLineWidth( minLineWidth );
				    Chi2Min1S->SetLineStyle( minLineStyle );
				    Chi2Min1S->SetLineColor( fLSBLineColor[iframe][nStateBuff] );
				    Chi2Min1S->Draw( "same" );
				    sprintf(text,"#color[2]{Minimum:} f_{LSB} = %1.2f, #chi^{2} / ndf = %1.3f / %d = %1.3f,  #chi^{2}-prob. = %1.3f", minChi2FracL[iframe][nStateBuff], Chi2Min[iframe][nStateBuff], NDFMin[iframe][nStateBuff], Chi2Min[iframe][nStateBuff]/double(NDFMin[iframe][nStateBuff]), TMath::Prob(Chi2Min[iframe][nStateBuff], NDFMin[iframe][nStateBuff]));
				    TLatex *textChi2_1S = new TLatex(0.02,yMin+(yMax-yMin)*0.05,text);
				    textChi2_1S->SetTextSize(FontSizeChi2)                                                                                                                                                                                                                                             ;
				    textChi2_1S->Draw( "same" )                                                                                                                                                                                                                                                 ;
				}
	  		    sprintf(savenamePlot,"FigBuffer/BGratioFits/%s/Chi2Func_%s.pdf",PlotDir,frameIdentity);
	  		    plotCanvas->SaveAs(savenamePlot);

		  }



}
