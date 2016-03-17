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

void FitRho(){

	  gROOT->Reset();
	  gROOT->SetBatch();

	  TFile* inFile = new TFile("/afs/hephy.at/scratch/k/knuenz/CMSSW_4_2_4_patch2/src/UpsilonPol/macros/polFit/EffFiles/rhofactor.root");
	  TFile* inFile = new TFile("/afs/hephy.at/scratch/k/knuenz/CMSSW_4_2_4_patch2/src/UpsilonPol/macros/polFit/EffFiles/rhoFactor_Ups3S_MCTruthEff_29March2012_finalPTBins_TPVCuts_EachEvUsed40x_10x12cosThetaPhiBins_maxVariation.root");


  TF2 *fcosthphi = new TF2( "fcosthphi", "[0]*(1.+[1]*x[0]*x[0]+[2]*(1.-x[0]*x[0])*cos(2.*x[1]*0.0174532925)+[3]*2.*x[0]*sqrt(1.-x[0]*x[0])*cos(x[1]*0.0174532925))", -1., 1., -180., 180. );
  fcosthphi->SetParameters(1., 0.0, 0.0, 0.0);

  const int npT = 11;
  const int nrap = 3;
  const int nframe = 4;

  double  lth[npT][nrap][nframe],   dlth[npT][nrap][nframe];
  double  lph[npT][nrap][nframe],   dlph[npT][nrap][nframe];
  double  lthph[npT][nrap][nframe], dlthph[npT][nrap][nframe];

  char histoName[50];

  for( int iframe  = 1; iframe  < nframe; iframe++ ){
  for( int irap  = 1; irap  < nrap; irap++ ){
  for( int ipT   = 1; ipT   < npT;  ipT++  ){


	  if(iframe==1)    sprintf ( histoName, "hRhoWErr_tot_CS_rap%u_pT%u", irap, ipT );
	  if(iframe==2)    sprintf ( histoName, "hRhoWErr_tot_HX_rap%u_pT%u", irap, ipT );
	  if(iframe==3)    sprintf ( histoName, "hRhoWErr_tot_PHX_rap%u_pT%u", irap, ipT );

    cout << histoName << endl;

    TH2D* histo  = (TH2D*)inFile->Get( histoName );

    histo->Fit("fcosthphi");

     lth[ipT][irap][iframe]   = fcosthphi->GetParameter(1);
    dlth[ipT][irap][iframe]   = fcosthphi->GetParError(1);
     lph[ipT][irap][iframe]   = fcosthphi->GetParameter(2);
    dlph[ipT][irap][iframe]   = fcosthphi->GetParError(2);
     lthph[ipT][irap][iframe] = fcosthphi->GetParameter(3);
    dlthph[ipT][irap][iframe] = fcosthphi->GetParError(3);

  }
  }
  }

  cout << endl;

  for( int iframe  = 1; iframe  < nframe; iframe++ ){
  for( int irap  = 1; irap  < nrap; irap++ ){
  for( int ipT   = 1; ipT   < npT;  ipT++  ){

      cout << "pT" << ipT << ", rap" << irap  << " in frame "<<iframe << ":  "
           <<          lth[ipT][irap][iframe]   << " +- " << dlth[ipT][irap][iframe]
           << "   " << lph[ipT][irap][iframe]   << " +- " << dlph[ipT][irap][iframe]
           << "   " << lthph[ipT][irap][iframe] << " +- " << dlthph[ipT][irap][iframe]
           << endl;
  }
  }
  }






	char filename[200];
	sprintf(filename,"/afs/hephy.at/scratch/k/knuenz/CMSSW_4_2_4_patch2/src/UpsilonPol/macros/polFit/Systematics/TheGreatRun_Rho/RhoFactorFit/TGraphResults_3SUps.root");
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

		if(iLam==1) {lmean[pt]=lth[ptBin][rapBin][1]; errlmean[pt]=dlth[ptBin][rapBin][1];}
		if(iLam==2) {lmean[pt]=lph[ptBin][rapBin][1];errlmean[pt]=dlph[ptBin][rapBin][1];}
		if(iLam==3) {lmean[pt]=lthph[ptBin][rapBin][1];errlmean[pt]=dlthph[ptBin][rapBin][1];}
		if(iLam==4 || iLam==5) {lmean[pt]=0;errlmean[pt]=0;}
		if(iLam==6) {lmean[pt]=(lth[ptBin][rapBin][1]+3*lph[ptBin][rapBin][1])/(1-lph[ptBin][rapBin][1]);}

		if(iLam==7) {lmean[pt]=lth[ptBin][rapBin][2]; errlmean[pt]=dlth[ptBin][rapBin][2];}
		if(iLam==8) {lmean[pt]=lph[ptBin][rapBin][2];errlmean[pt]=dlph[ptBin][rapBin][2];}
		if(iLam==9) {lmean[pt]=lthph[ptBin][rapBin][2];errlmean[pt]=dlthph[ptBin][rapBin][2];}
		if(iLam==10 || iLam==11) {lmean[pt]=0;errlmean[pt]=0;}
		if(iLam==12) {lmean[pt]=(lth[ptBin][rapBin][2]+3*lph[ptBin][rapBin][2])/(1-lph[ptBin][rapBin][2]);}

		if(iLam==13) {lmean[pt]=lth[ptBin][rapBin][3]; errlmean[pt]=dlth[ptBin][rapBin][3];}
		if(iLam==14) {lmean[pt]=lph[ptBin][rapBin][3];errlmean[pt]=dlph[ptBin][rapBin][3];}
		if(iLam==15) {lmean[pt]=lthph[ptBin][rapBin][3];errlmean[pt]=dlthph[ptBin][rapBin][3];}
		if(iLam==16 || iLam==17) {lmean[pt]=0;errlmean[pt]=0;}
		if(iLam==18) {lmean[pt]=(lth[ptBin][rapBin][3]+3*lph[ptBin][rapBin][3])/(1-lph[ptBin][rapBin][3]);}

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
