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

void Chi2FitBGSBratio(){

	  gROOT->Reset();
	  gROOT->SetBatch();


	  char DataID[200];
	  sprintf(DataID,"SetOfCuts11_May20Centrals");
	  int nSigma=1;

	  char PlotDir[200];
	  sprintf(PlotDir,"Chi2_SBratio_May24_2D_AbsCosthMax1_with1DPlots");

	  int FitDimension=2;
	  double AbsCosthMax=1;
	  double RelativeErrorMax=100.;

	  char PlotDirFull[200];
	  sprintf(PlotDirFull,"FigBuffer/BGratioFits/%s",PlotDir);
	  gSystem->mkdir(PlotDirFull);
	  char SystDirFull[200];
	  sprintf(SystDirFull,"Systematics/TheGreatRun_BKGratio/%s",PlotDir);
	  gSystem->mkdir(SystDirFull);


	  char frameIdentity[50];
	  bool SetModelErrorZero=false;

	  int rapMin=1;
	  int rapMax=2;
	  int ptMin=6;
	  int ptMax=10;
	  int frameMin=3;
	  int frameMax=3;
	  int stateMin=1;
	  int stateMax=1;

	  const int nStates=3;
	  const int npT = 11;//11
	  const int nrap = 3;//3
	  const int nframe = 4;//4

	  int nFracL=1;
	  int iFracL=0;

	  double fracLvar[nFracL][npT][nrap][nframe][nStates];
	  double FitConst[nFracL][npT][nrap][nframe][nStates];
	  double err_FitConst[nFracL][npT][nrap][nframe][nStates];
	  double FitConst_Chi2[nFracL][npT][nrap][nframe][nStates];
	  double FitConst_NDF[nFracL][npT][nrap][nframe][nStates];

	  for( int nState  = stateMin; nState  < stateMax+1; nState++ ){



		  TFile* inFile;


		  TF1 *fcosthphi1D = new TF1( "fcosthphi1D", "[0]", -AbsCosthMax, AbsCosthMax);
		  fcosthphi1D->SetParameter(0,1.);

		  TF2 *fcosthphi2D = new TF2( "fcosthphi2D", "[0]+(0*x[0]+0*x[1])", -AbsCosthMax, AbsCosthMax, -180., 180. );
		  fcosthphi2D->SetParameter(0,1.);


  char histoName[50];

  TH2D* histoBGmodel;
  TH2D* histoBGdata;
  TH2D* histoBGratio;
  TH2D* histoBGmodelR;
  TH2D* histoBGmodelL;
  TH1D* histoBGmodelR1D;
  TH1D* histoBGmodelL1D;

  TH1D* histoBGmodel1D;
  TH1D* histoBGdata1D;
  TH1D* histoBGratio1D;

  TH1D* histoBGratio1D_rap1_pT0;
  TH1D* histoBGratio1D_rap2_pT0;

  TF1* Clone_rap1_pT0;
  TF1* Clone_rap2_pT0;

  TLatex *fitBGratio1DtexSumm1;
  TLatex *fitBGratio1DtexSumm2;

  TCanvas *c3;

  for( int iframe  = frameMin; iframe  < frameMax+1; iframe++ ){

	  if(iframe==1) sprintf(frameIdentity,"CS");
	  if(iframe==2) sprintf(frameIdentity,"HX");
	  if(iframe==3) sprintf(frameIdentity,"PX");

  for( int irap  = rapMin; irap  < rapMax+1; irap++ ){
  for( int ipT   = ptMin; ipT   < ptMax+1;  ipT++  ){

	  char inFileName[500];
	  sprintf(inFileName,"/afs/hephy.at/scratch/k/knuenz/CMSSW_4_2_4_patch2/src/UpsilonPol/macros/DataFiles/%s/AllStates_%d.00Sigma_FracLSB50Percent/data_%dSUps_rap%d_pT%d.root",DataID,nSigma,nState,irap,ipT);
	    inFile = new TFile(inFileName,"READ");

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

	  TH2D* histoBGmodelRclone  = (TH2D*)histoBGmodelR->Clone( "histoBGmodelRclone" );
	  TH2D* histoBGmodelLclone  = (TH2D*)histoBGmodelL->Clone( "histoBGmodelLclone" );

	  histoBGmodelRclone->Sumw2();
	  histoBGmodelLclone->Sumw2();


	////////// Cut on costh //////////////////////


	  histoBGmodelR1D=histoBGmodelRclone->ProjectionX("histoBGmodelR1D", 0, -1,"e");
	  histoBGmodelR1D->Sumw2();

	    if(SetModelErrorZero){
		    for(int iX=0;iX<histoBGmodelR1D->GetNbinsX()+1;iX++){
		        	int globalBin=histoBGmodelR1D->GetBin(iX);
		        	histoBGmodelR1D->SetBinError(iX,0.);
		        }
	    }

	    histoBGmodelL1D=histoBGmodelLclone->ProjectionX("histoBGmodelL1D", 0, -1,"e");

	    histoBGratio1D=(TH1D*)histoBGmodelL1D->Clone("histoBGratio1D");
	    histoBGratio1D->Divide(histoBGmodelR1D);
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


		    histoBGratio=(TH2D*)histoBGmodelL->Clone("histoBGratio");
		    histoBGratio->Divide(histoBGmodelR);

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

			      fracLvar[iFracL][ipT][irap][iframe][nState] = 0.;
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

			      fracLvar[iFracL][ipT][irap][iframe][nState] = 0.;
			      FitConst[iFracL][ipT][irap][iframe][nState] = fcosthphi2D->GetParameter(0);
			      err_FitConst[iFracL][ipT][irap][iframe][nState] = fcosthphi2D->GetParError(0);
			      FitConst_Chi2[iFracL][ipT][irap][iframe][nState] = fcosthphi2D->GetChisquare();
			      FitConst_NDF[iFracL][ipT][irap][iframe][nState] = fcosthphi2D->GetNDF();
		    }


	    	      cout<<"FitConst "<<FitConst[iFracL][ipT][irap][iframe][nState]<<endl;
	    	      cout<<"err_FitConst "<<err_FitConst[iFracL][ipT][irap][iframe][nState]<<endl;
	    	      cout<<"FitConst_Chi2 "<<FitConst_Chi2[iFracL][ipT][irap][iframe][nState]<<endl;
	    	      cout<<"FitConst_NDF "<<FitConst_NDF[iFracL][ipT][irap][iframe][nState]<<endl;


	    char savenamePlot[200];
	    gStyle->SetPalette(1);


	    TCanvas *c2 = new TCanvas("c2","c2",1200,800);
	    gPad->SetFillColor(kWhite);

	    c2->cd();


	    if(FitDimension==2){
	    histoBGratio->SetTitle(0);
	    histoBGratio->SetStats(0);
	    histoBGratio->Draw("colz");
		char text[200];
		sprintf(text,"#chi^{2} / NDF = %1.3f / %1.0f = %1.3f",FitConst_Chi2[iFracL][ipT][irap][iframe][nState],FitConst_NDF[iFracL][ipT][irap][iframe][nState],FitConst_Chi2[iFracL][ipT][irap][iframe][nState]/FitConst_NDF[iFracL][ipT][irap][iframe][nState]);
		TLatex *fitBGratio1Dtex = new TLatex(-0.8,190.,text);
		fitBGratio1Dtex->SetTextSize(0.03)                                                                                                                                                                                                                                             ;
		fitBGratio1Dtex->Draw( "same" );
	    sprintf(savenamePlot,"FigBuffer/BGratioFits/%s/histoBGratio_LSBoverRSB2D_%s_rap%d_pT%d.pdf",PlotDir,frameIdentity,irap,ipT);
	    c2->SaveAs(savenamePlot);

	    histoBGratio1D->SetTitle(0);
	    histoBGratio1D->SetStats(0);
	    histoBGratio1D->Draw("colz");
	    fcosthphi1D->SetLineColor(kGreen+2);
	    fcosthphi1D->SetLineWidth(2.);
	    fcosthphi1D->SetLineStyle(1);
		fcosthphi1D->SetParameter(0,fcosthphi1D->GetParameter(0));
	    fcosthphi1D->Draw("same");
		sprintf(text,"#chi^{2} / NDF = %1.3f / %1.0f = %1.3f",FitConst_Chi2[iFracL][ipT][irap][iframe][nState],FitConst_NDF[iFracL][ipT][irap][iframe][nState],FitConst_Chi2[iFracL][ipT][irap][iframe][nState]/FitConst_NDF[iFracL][ipT][irap][iframe][nState]);
		TLatex *fitBGratio1Dtex2 = new TLatex(-0.8,histoBGratio1D->GetMaximum()*0.2,text);
		fitBGratio1Dtex2->SetTextSize(0.03)                                                                                                                                                                                                                                             ;
		fitBGratio1Dtex2->Draw( "same" )                                                                                                                                                                                                                                                 ;
		sprintf(savenamePlot,"FigBuffer/BGratioFits/%s/histoBGratio_LSBoverRSB1D_%s_rap%d_pT%d.pdf",PlotDir,frameIdentity,irap,ipT);
		c2->SaveAs(savenamePlot);

	    }
	    if(FitDimension==1){
	    histoBGratio1D->SetTitle(0);
	    histoBGratio1D->SetStats(0);
	    histoBGratio1D->Draw("colz");
	    fcosthphi1D->SetLineColor(kGreen+2);
	    fcosthphi1D->SetLineWidth(2.);
	    fcosthphi1D->SetLineStyle(1);
	    fcosthphi1D->Draw("same");
		char text[200];
		sprintf(text,"#chi^{2} / NDF = %1.3f / %1.0f = %1.3f",FitConst_Chi2[iFracL][ipT][irap][iframe][nState],FitConst_NDF[iFracL][ipT][irap][iframe][nState],FitConst_Chi2[iFracL][ipT][irap][iframe][nState]/FitConst_NDF[iFracL][ipT][irap][iframe][nState]);
		TLatex *fitBGratio1Dtex = new TLatex(-0.8,histoBGratio1D->GetMaximum()*0.2,text);
		fitBGratio1Dtex->SetTextSize(0.03)                                                                                                                                                                                                                                             ;
		fitBGratio1Dtex->Draw( "same" )                                                                                                                                                                                                                                                 ;
	    sprintf(savenamePlot,"FigBuffer/BGratioFits/%s/histoBGratio_LSBoverRSB1D_%s_rap%d_pT%d.pdf",PlotDir,frameIdentity,irap,ipT);
	    c2->SaveAs(savenamePlot);
	    }


    inFile->Close();

  }
  }



  }






	char filename[200];
	sprintf(filename,"/afs/hephy.at/scratch/k/knuenz/CMSSW_4_2_4_patch2/src/UpsilonPol/macros/polFit/Systematics/TheGreatRun_BKGratio/%s/TGraphResults_%dSUps.root",PlotDir,nState);
	gSystem->Unlink(filename);
	cout<<filename<<endl;
	char GraphName[200];


	for(int iLam = 1; iLam<19; iLam++){

	for(int rapBin = 1; rapBin < nrap; rapBin++){

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

	int nframe;

	if(iLam>0&&iLam<7) nframe=1;
	if(iLam>6&&iLam<13) nframe=2;
	if(iLam>12&&iLam<19) nframe=3;


	const int nBinspT=10;
	double ptCentre_[nBinspT]={5.747, 6.50359, 7.60069, 8.47588, 9.59466, 10.9576, 13.7561, 17.7305, 23.5741, 35.862};
	double pTrange[nBinspT+1]={5., 6., 7., 8., 9., 10., 12., 16., 20., 30., 50.};
	double ptCentreErr_low[nBinspT];
	double ptCentreErr_high[nBinspT];
	double lmean[nBinspT];
	double errlmean[nBinspT];

	int pt=0;
	for(int ptBin = 1; ptBin < npT; ptBin++) {


		ptCentreErr_low[pt]=ptCentre_[pt]-pTrange[pt];
		ptCentreErr_high[pt]=pTrange[pt+1]-ptCentre_[pt];

		if(FitConst_NDF[0][ptBin][rapBin][nframe][1]<0.5)FitConst_NDF[0][ptBin][rapBin][nframe][1]=100000;
		lmean[pt]=FitConst_Chi2[0][ptBin][rapBin][nframe][1]/FitConst_NDF[0][ptBin][rapBin][nframe][1]; errlmean[pt]=0;


	pt++;
	}


	TGraphAsymmErrors *graphSyst = new TGraphAsymmErrors(nBinspT,ptCentre_,lmean,ptCentreErr_low,ptCentreErr_high,errlmean,errlmean);
	graphSyst->SetMarkerColor(kRed);
	graphSyst->SetLineColor(kRed);
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







}
