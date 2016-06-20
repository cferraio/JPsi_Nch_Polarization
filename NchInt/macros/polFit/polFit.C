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
#include "TGraph2D.h"
#include "TGraph.h"
#include "TEfficiency.h"

//#include "effsAndCuts.h"

// number of random extractions for the sampling of the parameter space
// (the number of entries in the results ntuple will be smaller than this,
// due to efficiency < 1 of the sampling method and the rejected extractions of the burn-in period):
//const int n_sampledPoints = 12000;  // this INCLUDES the burn-in below.
                                    // 12000 is already perfectly good for obtaining central values and errors
                                    // more itetarions necessary for smooth 2D distributions

// number of discarded intial random extractions (burn-in period):
int n_burnIn = 10000; // do not change this, 2000 is default
int Toy_n_burnIn = 2000; // do not change this, 2000 is default 
//int Toy_n_burnIn = 10000; // do not change this, 2000 is default // for Gen-Fiducial test
bool RandomSeed=true;
bool RandomLambdaTheta=false;

// Proposal width to be multiplied by the RMS of the test distribution from the burn in, default: 0.3
double proposalWidthTimesRMS = 1.0;//0.3 default
double proposalWidthBurnIn = 0.1;//0.1 default

// number of MC events for acceptance calculation
const double n_MCevents_ACC  = 1000000;
bool subtractBG=true;

// number of MC events for the calculation of the likelihood normalization
// (preliminary integration of the parameter-independent parts of the likelihood)
const double n_MCevents_NORM = 1000000;

// Numerical inputs to calculate decay angles:
const double pbeam_ = 7000.; // exact number irrelevant as long as pbeam >> Mprot
const double Mprot_ = 0.9382720;
const double Mlepton_ = 0.10566;  // (muon)
const double gPI_ = TMath::Pi();
const double Ebeam_ = sqrt( pbeam_*pbeam_ + Mprot_*Mprot_ );
TLorentzVector beam1_LAB_( 0., 0., pbeam_, Ebeam_ );
TLorentzVector beam2_LAB_( 0., 0., -pbeam_, Ebeam_ );

// polarization assumed for the calculation of the acceptance as a function of pT, y, mass
const double lambda_theta_ACC = 0.0;
const double lambda_phi_ACC =  0.0;
const double lambda_thetaphi_ACC = 0.0;

// reference frame
const bool HX_is_natural_ACC = true;  // put both to false to generate in the CS frame
const bool PX_is_natural_ACC = false;


// minimum value of the dilepton efficiency
const double min_dileptonEff = 0.01;



bool isMuonInAcceptance(int iCut, double pT, double eta);
double singleLeptonEfficiency( double& pT, double& eta, int nEff, TFile *fInEff, TH2D* hEvalEff, bool MCeff, TEfficiency* TEff);
void EvaluateEffFileName(int nEff, char EffFileName [200], bool singleLeptonEff);
double DiLeptonEfficiency( double& Dilepton_pT, double& Dilepton_rap, int nDileptonEff, TFile *fInDileptonEff, bool MCeff);
double EvaluateRhoFactor( double& costh, double& phi, int nEff, TFile* fInRhoFactor, double rap, double pT, bool StatVarRho);
void FindMPV(TH1* PosteriorDist , double& MPV , double& MPVerrorLow, double& MPVerrorHigh, int MPValgo, int nSigma);
double DenominatorAmapEfficiency( double& pT, double& eta, int nDenominatorAmap, TFile *fInEff_nDenominatorAmap, TH2D* hEvalEff_nDenominatorAmap, bool MCeff, TEfficiency* TEff_nDenominatorAmap);
double EvaluateAmap( double& costh_Amap, double& phi_Amap, int nAmap, TFile* fInAmap, double rap, double pT);

////////////////////////////////////////////////////////////////////////////////////////////////////
// Function randomly sampling the lambda values of the "next iteration"
// given the current values, according to a Gaussian "proposal pdf"
// including or not the positivity constraints of the angular distribution:

void extractFromProposalPDF( double& lth_candidate,     double& lph_candidate,     double& ltp_candidate,
                             double& lth,               double& lph,               double& ltp,
                             double& proposalWidth_lth, double& proposalWidth_lph, double& proposalWidth_ltp ) {

  if ( proposalWidth_lth < 0. || proposalWidth_lph < 0. || proposalWidth_ltp < 0. ) {
       do {
           lth_candidate = gRandom->Gaus( lth, proposalWidthBurnIn ); // Gaussian proposal pdf with "large" sigma and positivity constraints
           lph_candidate = gRandom->Gaus( lph, proposalWidthBurnIn ); // 0.10 all default
           ltp_candidate = gRandom->Gaus( ltp, proposalWidthBurnIn );
       }

/*    do {
        lth_candidate = gRandom->Uniform(-1,1); // Gaussian proposal pdf with "large" sigma and positivity constraints
        lph_candidate = gRandom->Uniform(-1,1);
        ltp_candidate = gRandom->Uniform(-0.707106781186548,0.707106781186548);
    }
*/	while ( TMath::Abs( lph_candidate ) > 0.5*( 1 + lth_candidate ) || lth_candidate*lth_candidate + 2.*ltp_candidate*ltp_candidate > 1
            || TMath::Abs( ltp_candidate ) > 0.5*( 1 - lph_candidate )
            || (  (1.+2.*lph_candidate)*(1.+2.*lph_candidate) + 2.*ltp_candidate*ltp_candidate > 1 && lph_candidate < -1./3. ) );
  }
  else {
      lth_candidate = gRandom->Gaus ( lth, proposalWidth_lth ); // Gaussian proposal pdf with possibly smaller sigmas
      lph_candidate = gRandom->Gaus ( lph, proposalWidth_lph ); // and no positivity constraints
      ltp_candidate = gRandom->Gaus ( ltp, proposalWidth_ltp );
  }

}

// starting values of the widths of the proposal functions to scan the parameter space.
// Negative = the burn-in period starts with "large" sigma and positivity constraints

const double proposalWidth_lth_CS_start = -1.;
const double proposalWidth_lph_CS_start = -1.;
const double proposalWidth_ltp_CS_start = -1.;

const double proposalWidth_lth_HX_start = -1.;
const double proposalWidth_lph_HX_start = -1.;
const double proposalWidth_ltp_HX_start = -1.;

const double proposalWidth_lth_PX_start = -1.;
const double proposalWidth_lph_PX_start = -1.;
const double proposalWidth_ltp_PX_start = -1.;


// calculation of the frame-invariants lambdatheta^star and lambdaphi^star:

void calcLambdastar( double& lthstar, double& lphstar,
                     double& lth,     double& lph,     double& ltp ) {

  double LamPlus = 0.25 * ( lth - lph + sqrt( pow( lth - lph, 2. ) + 4. * ltp*ltp ) );
  double LamMnus = 0.25 * ( lth - lph - sqrt( pow( lth - lph, 2. ) + 4. * ltp*ltp ) );

  lthstar = ( lth - 3.*LamPlus ) / ( 1. + LamPlus );
  lphstar = ( lph + LamPlus )    / ( 1. + LamPlus );

  double lphstarMnus = ( lph + LamMnus ) / ( 1. + LamMnus );

  if ( TMath::Abs(lphstarMnus) < TMath::Abs(lphstar) ) {
     lthstar = ( lth - 3.*LamMnus ) / ( 1. + LamMnus );
     lphstar = lphstarMnus;
  }
}




///////////// FindMPV //////////////////////////////////////////////////////////////////////////////

/*
void FindMPV(TH1* PosteriorDist , double& MPV , double& MPVerrorLow, double& MPVerrorHigh, int MPValgo){
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
double OneSigmaCL = 0.682689492137;
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
*/

////////////////////////////////////////////////////////////////////////////////////////////////////
///////////// main program /////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////

void polFit(int n_sampledPoints=1,
		int FidCuts=0,
		int nEff=1,
		int nDileptonEff=1,
		int nRhoFactor=1,
		Char_t *dirstruct = "OutputDirectory_Default",
		Char_t *realdatadir = "RealDataDirectory_Default",
		Char_t *TreeBinID = "TreeBinID_Default",
		Char_t *TreeBinID_dataFile = "TreeBinID_Default",
		bool RealData=false,
		Char_t *effDir = "effDir_Default",
		bool MCeff=false,
		bool MCDileptoneff=false,
		int rapBin=999,
		int pTbin=999,
		bool NewAccCalc=true,
		int MPValgo=999,
		bool useAmapApproach=false,
		int nAmap=999,
		int nDenominatorAmap=999,
		bool StatVarTotBGfraction=false,
		bool StatVarTotBGmodel=false,
		bool StatVarRho=false,
		bool cutDeltaREllDpt=false){

	cout<<"/////////////////////////////////"<<endl;
	cout<<"running polFit.C ........////////"<<endl;
	cout<<"/////////////////////////////////"<<endl;

	bool PhiHX_test = true;

  gROOT->Reset();

	double DeltaREllDptValue = 0.18;

  delete gRandom;

  if(RandomSeed) gRandom = new TRandom3(0);  // better random generator
  else gRandom = new TRandom3(23101987);

  bool ForceEpsSmallerOne=false;

  if(!RealData) n_burnIn=Toy_n_burnIn;

//Get single Lepton Efficiency File name
  char EffFile[200];
  char EffFileName[200];
  EvaluateEffFileName(nEff,EffFileName,true);
  sprintf(EffFile,"%s/%s",effDir,EffFileName);

	char EffType[200];
	if(MCeff) sprintf(EffType,"hEff_MC_central");
	else sprintf(EffType,"hEff_DATA_central");

  TFile *fInEff = new TFile(EffFile);
  TH1* hEff=(TH1*) fInEff->Get(EffType);

//Get DiLepton Efficiency File name

  EvaluateEffFileName(nDileptonEff,EffFileName,false);
  sprintf(EffFile,"%s/%s",effDir,EffFileName);

  TFile *fInDileptonEff = new TFile(EffFile);
  TH1* hDileptonEff=(TH1*) fInDileptonEff->Get(EffType);

  EvaluateEffFileName(nRhoFactor,EffFileName,false);
  sprintf(EffFile,"%s/%s",effDir,EffFileName);
  TFile *fInRhoFactor = new TFile(EffFile);
	//cout<<"EffFile: "<<EffFile<<endl;


  cout<<"Filling Efficiency Evaluation Histogram"<<endl;
  TH2D*  hEvalEff1D;
  TH2D*  hEvalEff2D;
  char hEvalEffName[200];

  double AlmostZero=1e-8;

  if(fInEff->Get(EffType) != NULL){
 int pTBinsTotal = hEff -> GetNbinsX();
  int etaBinsTotal = hEff -> GetNbinsY();
  int pTBinsNew = 2000;
  int etaBinsNew = 200;
  hEvalEff1D   = new TH2D( "hEvalEff1D", "hEvalEff1D", etaBinsNew, hEff->GetXaxis()->GetXmin(),1.6, pTBinsNew,  hEff->GetYaxis()->GetXmin(), 100);
  hEvalEff2D   = new TH2D( "hEvalEff2D", "hEvalEff2D", etaBinsNew, hEff->GetXaxis()->GetXmin(),1.6, pTBinsNew,  hEff->GetYaxis()->GetXmin(), 100);
  double eff;
  double effBuffer;
  char graphName[200], graphName2D[200];

//  if(fInEff->Get("gEff2D_MC") != NULL){

  if(MCeff) sprintf(graphName2D,"gEff2D_MC");
  else sprintf(graphName2D,"gEff2D_DATA");
//  TGraph2D *graph2D = new TGraph2D(*((TGraph2D *) fInEff->Get(graphName2D)));

  for(int etaBin=0;etaBin<etaBinsNew;etaBin++){
	  int currentEtaBin = hEff->GetXaxis()->FindBin(hEvalEff1D->GetXaxis()->GetBinCenter(etaBin+1));
	  if(MCeff) sprintf(graphName,"gEff_MC_PT_AETA%d",currentEtaBin-1);
	  else sprintf(graphName,"gEff_DATA_PT_AETA%d",currentEtaBin-1);
	  TGraphAsymmErrors *graph = new TGraphAsymmErrors(*((TGraphAsymmErrors *) fInEff->Get(graphName)));
	  double pTMaxCenter=-999;
	  double EtaConst=-999;
	  double pTcenterGraphM1=-999;
	  for(int pTBinGraph=0;pTBinGraph<100;pTBinGraph++){
		double pTcenterGraph;
		double effGraph;
		graph->GetPoint(pTBinGraph,pTcenterGraph,effGraph);
		if(TMath::Abs(pTcenterGraph-pTcenterGraphM1)<AlmostZero) {pTMaxCenter=pTcenterGraph; EtaConst=effGraph; break;}
		pTcenterGraphM1=pTcenterGraph;
	  }

		  for(int pTBin=0;pTBin<pTBinsNew;pTBin++){
			  eff = graph->Eval(hEvalEff1D->GetYaxis()->GetBinCenter(pTBin+1));
		  	  if(eff<0) eff=0;
		  	  if(hEvalEff1D->GetYaxis()->GetBinCenter(pTBin+1)>pTMaxCenter) eff=EtaConst;
		  	  hEvalEff1D->SetBinContent(etaBin+1,pTBin+1,eff); effBuffer=eff;
//			  eff = graph2D->Interpolate(hEvalEff1D->GetYaxis()->GetBinCenter(pTBin+1),hEvalEff1D->GetXaxis()->GetBinCenter(etaBin+1));
		  	  if(eff<AlmostZero && currentEtaBin==1) eff=effBuffer;
		  	  if(eff<AlmostZero && hEvalEff1D->GetYaxis()->GetBinCenter(pTBin+1)>pTMaxCenter*0.5) eff=effBuffer;
			  if(eff<0) eff=0;
		  	  hEvalEff2D->SetBinContent(etaBin+1,pTBin+1,eff);
  }
  }
  sprintf(hEvalEffName,"%s/EvalHisto1D.root",dirstruct);
//  hEvalEff1D->SaveAs(hEvalEffName);
  sprintf(hEvalEffName,"%s/EvalHisto2D.root",dirstruct);
	//  hEvalEff2D->SaveAs(hEvalEffName);
	}




	if( nEff>1019 && nEff<1099 ){
		const int etaBinsTotal = 8;
		double etaBinningParametrized[etaBinsTotal+1]={0.,0.2,0.3,0.6,0.8,1.0,1.2,1.4,1.6};
		int pTBinsNew = 2000;
		int etaBinsNew = 200;
		hEvalEff1D   = new TH2D( "hEvalEff1D", "hEvalEff1D", etaBinsNew, 0,1.6, pTBinsNew,  0, 100);
		double eff;
		double effBuffer;
		char graphName[200];


		int currentEtaBin;
		for(int etaBin=0;etaBin<etaBinsNew;etaBin++){
			for(int etaSearch=0;etaSearch<etaBinsTotal;etaSearch++){
				if(hEvalEff1D->GetXaxis()->GetBinCenter(etaBin+1)>etaBinningParametrized[etaSearch] && hEvalEff1D->GetXaxis()->GetBinCenter(etaBin+1)<etaBinningParametrized[etaSearch+1])
					currentEtaBin=etaSearch+1;
			}
			//cout<<"debug--currentEtaBin: "<<currentEtaBin<<endl;
			sprintf(graphName,"gEff_DATA_PT_AETA%d",currentEtaBin-1);
			TGraphAsymmErrors *graph = new TGraphAsymmErrors(*((TGraphAsymmErrors *) fInEff->Get(graphName)));

			for(int pTBin=0;pTBin<pTBinsNew;pTBin++){
				eff = graph->Eval(hEvalEff1D->GetYaxis()->GetBinCenter(pTBin+1));
				if(eff<0) eff=0;
				hEvalEff1D->SetBinContent(etaBin+1,pTBin+1,eff); effBuffer=eff;
			}

		}
		sprintf(hEvalEffName,"%s/EvalHisto1D.root",dirstruct);
		//  hEvalEff1D->SaveAs(hEvalEffName);
	}


	if( nEff==1101||nEff==1102 ){
		const int etaBinsTotal = 16;
		double etaBinningParametrized[etaBinsTotal+1]={0.,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6};
		int pTBinsNew = 2000;
		int etaBinsNew = 200;
		hEvalEff1D   = new TH2D( "hEvalEff1D", "hEvalEff1D", etaBinsNew, 0,1.6, pTBinsNew,  0, 100);
		double eff;
		double effBuffer;
		char graphName[200];


		int currentEtaBin;
		for(int etaBin=0;etaBin<etaBinsNew;etaBin++){
			for(int etaSearch=0;etaSearch<etaBinsTotal;etaSearch++){
				if(hEvalEff1D->GetXaxis()->GetBinCenter(etaBin+1)>etaBinningParametrized[etaSearch] && hEvalEff1D->GetXaxis()->GetBinCenter(etaBin+1)<etaBinningParametrized[etaSearch+1])
					currentEtaBin=etaSearch+1;
			}
			sprintf(graphName,"gEff_DATA_PT_AETA%d",currentEtaBin-1);
			TGraphAsymmErrors *graph = new TGraphAsymmErrors(*((TGraphAsymmErrors *) fInEff->Get(graphName)));

			for(int pTBin=0;pTBin<pTBinsNew;pTBin++){
				eff = graph->Eval(hEvalEff1D->GetYaxis()->GetBinCenter(pTBin+1));
				if(eff<0) eff=0;
				hEvalEff1D->SetBinContent(etaBin+1,pTBin+1,eff); effBuffer=eff;
			}

		}
		sprintf(hEvalEffName,"%s/EvalHisto1D.root",dirstruct);
		//  hEvalEff1D->SaveAs(hEvalEffName);
	}





	TEfficiency* TEff;
	TH2D* hEvalEff;
	if(nEff==105 || nEff==106){
		sprintf(EffType,"totEff_MCTRUTH_pT_eta");
		TEfficiency* TEff2=(TEfficiency*) fInEff->Get(EffType);
		TH1* hEffTOT=(TH1*)TEff2->GetTotalHistogram();
		TH1* hEffPASS=(TH1*)TEff2->GetPassedHistogram();
		hEffPASS->Divide(hEffTOT);
		hEvalEff=(TH2D*)hEffPASS->Clone("hEffPASS");
		TEff=(TEfficiency*) fInEff->Get(EffType);
	}


	if(nEff!=105 && nEff!=106 && nEff!=1)
		hEvalEff = (TH2D*)hEvalEff1D->Clone("hEvalEff");
	if(nEff > 10000) hEvalEff = (TH2D*)hEvalEff2D->Clone("hEvalEff");

	sprintf(hEvalEffName,"%s/EvalHisto.root",dirstruct);
	//hEvalEff->SaveAs(hEvalEffName);





	cout<<"Get Amap"<<endl;

	////// Get Amap

	EvaluateEffFileName(nAmap,EffFileName,false);
	sprintf(EffFile,"%s/%s",effDir,EffFileName);

	TFile *fInAmap = new TFile(EffFile);





	cout<<"Get Amap Denominator"<<endl;


	////// Get Amap Denominator

	TEfficiency* TEff_nDenominatorAmap;
	TH2D* hEvalEff_nDenominatorAmap;

	TH2D*  hEvalEff1D_nDenominatorAmap;

	EvaluateEffFileName(nDenominatorAmap,EffFileName,true);
	sprintf(EffFile,"%s/%s",effDir,EffFileName);

	TFile *fInEff_nDenominatorAmap = new TFile(EffFile);

	if( nDenominatorAmap==1101 ){
		const int etaBinsTotal = 16;
		double etaBinningParametrized[etaBinsTotal+1]={0.,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6};
		int pTBinsNew = 2000;
		int etaBinsNew = 200;
		hEvalEff1D_nDenominatorAmap   = new TH2D( "hEvalEff1D_nDenominatorAmap", "hEvalEff1D_nDenominatorAmap", etaBinsNew, 0,1.6, pTBinsNew,  0, 100);
		double eff;
		double effBuffer;
		char graphName[200];


		int currentEtaBin;
		for(int etaBin=0;etaBin<etaBinsNew;etaBin++){
			for(int etaSearch=0;etaSearch<etaBinsTotal;etaSearch++){
				if(hEvalEff1D_nDenominatorAmap->GetXaxis()->GetBinCenter(etaBin+1)>etaBinningParametrized[etaSearch] && hEvalEff1D_nDenominatorAmap->GetXaxis()->GetBinCenter(etaBin+1)<etaBinningParametrized[etaSearch+1])
					currentEtaBin=etaSearch+1;
			}
			sprintf(graphName,"gEff_DATA_PT_AETA%d",currentEtaBin-1);
			TGraphAsymmErrors *graph = new TGraphAsymmErrors(*((TGraphAsymmErrors *) fInEff_nDenominatorAmap->Get(graphName)));

			for(int pTBin=0;pTBin<pTBinsNew;pTBin++){
				eff = graph->Eval(hEvalEff1D_nDenominatorAmap->GetYaxis()->GetBinCenter(pTBin+1));
				if(eff<0) eff=0;
				hEvalEff1D_nDenominatorAmap->SetBinContent(etaBin+1,pTBin+1,eff); effBuffer=eff;
			}

		}
	}

	if(nDenominatorAmap!=105 && nDenominatorAmap!=106 && nDenominatorAmap!=1)
		hEvalEff_nDenominatorAmap = (TH2D*)hEvalEff1D_nDenominatorAmap->Clone("hEvalEff_nDenominatorAmap");


	if(nDenominatorAmap==105 || nDenominatorAmap==106){
		sprintf(EffType,"totEff_MCTRUTH_pT_eta");
		TEfficiency* TEff2=(TEfficiency*) fInEff_nDenominatorAmap->Get(EffType);
		TH1* hEffTOT=(TH1*)TEff2->GetTotalHistogram();
		TH1* hEffPASS=(TH1*)TEff2->GetPassedHistogram();
		hEffPASS->Divide(hEffTOT);
		hEvalEff_nDenominatorAmap=(TH2D*)hEffPASS->Clone("hEffPASS");
		TEff_nDenominatorAmap=(TEfficiency*) fInEff_nDenominatorAmap->Get(EffType);
	}








  // input file with the background histograms
  // and the ntuple of dilepton events

  char filename [500];
  sprintf(filename,"%s/data.root",dirstruct);
  if(RealData) sprintf(filename,"%s/data_%s.root",realdatadir,TreeBinID_dataFile);
  TFile* dataFile = new TFile(filename);

	cout<<TreeBinID_dataFile<<endl;
	cout<<TreeBinID<<endl;

  TH2D* background_costhphiPX  = (TH2D*)dataFile->Get("background_costhphiPHX");
  TH3D* background_pTrapMass   = (TH3D*)dataFile->Get("background_pTrapMass");

  TH1D* background_costhPX = background_costhphiPX->ProjectionX( "background_costhPX" );
  TH1D* background_phiPX   = background_costhphiPX->ProjectionY( "background_phiPX" );

  TH1D* background_pT      = background_pTrapMass->ProjectionX( "background_pT" );
  TH1D* background_rap     = background_pTrapMass->ProjectionY( "background_rap" );
  TH1D* background_mass    = background_pTrapMass->ProjectionZ( "background_mass" );

  TH1D* background_fraction = (TH1D*)dataFile->Get("background_fraction");

  background_costhphiPX->Sumw2();
  background_costhPX->Sumw2();
  background_phiPX->Sumw2();
  background_pTrapMass->Sumw2();
  background_pT->Sumw2();
  background_rap->Sumw2();
  background_mass->Sumw2();

  TTree* data               = (TTree*)dataFile->Get("selectedData");

  TLorentzVector* lepP = 0;  data->SetBranchAddress( "lepP",  &lepP );  // lepton 4-vectors
  TLorentzVector* lepN = 0;  data->SetBranchAddress( "lepN",  &lepN );


  double  f_background = background_fraction->GetBinContent( 1 );
	cout << "f_background: " << f_background << endl;

	// apply statistical fluctuations on f_background with guassian distribution
	if(StatVarTotBGfraction && f_background>0.){
		double  f_backgroundPre    = background_fraction->GetBinContent( 1 );
		double  f_backgroundPreErr = background_fraction->GetBinError(1);
		cout << "f_backgroundErr: " << f_backgroundPreErr << endl;

		do {
			f_background = gRandom->Gaus( f_backgroundPre, f_backgroundPreErr );
		} while( f_background <=0. || f_background >= 2. * f_backgroundPre );

	  cout << "after apply statistical fluctuation: \n"
		 	<< "f_background: " << f_background << endl;
	}



  //////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////    CALCULATIONS    /////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////////
  //
  //  Likelihood(lambdas) = Product_{i=1...n_events} PDF(lepton_momenta[i], lambdas)
  //  where PDF(lepton_momenta[i], lambdas) = probability distribution of the dilepton events
  //
  //  PDF[i] =     W_signal_theor(costheta[i], phi[i], lambdas) * epsilon(lepton_momenta[i])
  //                                / normalization
  //
  //  where epsilon[i] = product of two single-lepton efficiencies
  //
  //  W_signal_theor[i]  is proportional to   1 + lambdatheta    * costheta[i]^2
  //                                            + lambdaphi      * sintheta[i]^2 cos2phi[i]
  //                                            + lambdathetaphi * sin2theta[i] cosphi[i]
  //
  //  This is the PDF for the signal events. The background is subtracted preliminary.
  //  The PDF must be normalized so that
  //         Sum_{i=1...n_events} PDF[i] = constant independent of parameters.
  //  In general the normalization *depends* on the lambda parameters, so that it
  //  should be re-calculated (by re-scanning all events) many times during the "fit",
  //  each time that the values of the parameters change. Fortunately, in our case
  //  PDF[i] is a linear function of the lambdas, so that its normalization can
  //  be calculated as an analytical function of the lambdas, therefore *before*
  //  starting the scan of the parameter space.
  //  We preliminary calculate the parts of the PDF that do not depend on the lambdas
  //  event by event, filling arrays with one element per event:
  //
  //  w_1[i]        = epsilon[i] * 1,
  //  w_cth2[i]     = epsilon[i] * costheta[i]^2,
  //  w_sth2c2ph[i] = epsilon[i] * sintheta[i]^2 cos2phi[i],
  //  w_s2thcph[i]  = epsilon[i] * sin2theta[i] cosphi[i],       i = single *real* envent
  //
  //  To calculate the normalization of the PDF as a function of the parameters, we integrate
  //  numerically these PDF terms.
  //  To perform the integration we generate events according to the efficiency-and-acceptance-corrected
  //  pT, eta and mass distributions (in the considered cell). Histograms
  //  for the dilepton variables corrected by efficiency are prepared preliminarly while
  //  we scan over all events in the ntuple. Acceptance histograms (in pT, eta and mass) are calculated
  //  by scanning uniformly the pT, y plane and the mass variable, generating muons and applying
  //  the fiducial acceptance cuts.
  //  The integration must be uniform over costh and phi in the acceptance region. This condition
  //  is realized by generating the decay distributions of the dileptons as uniform.
  //  We cannot use the data themselves (corrected by efficiency and acceptance) to calculate the integrals,
  //  just because the dilepton in the data are not unpolarized.
  //
  //  sum_1        = Sum_{k} w_1[k]
  //  sum_cth2     = Sum_{k} w_cth2[k]
  //  sum_sth2c2ph = Sum_{k} w_sth2c2ph[k]
  //  sum_s2thcph  = Sum_{k} w_s2thcph[k]        k = single *generated* (unpolarized) envent
  //
  /////////////////////////////////////////////////////////////////////////////////////////////////


  ///////////////////////////////////////////////////////////
  // Preliminary normalization calculations
  ///////////////////////////////////////////////////////////

  // boundaries of the kinematic cell taken from the background histograms:

  double pTdilepton_min = background_pT->GetXaxis()->GetXmin();
  double pTdilepton_max = background_pT->GetXaxis()->GetXmax();
  double rapdilepton_min = background_rap->GetXaxis()->GetXmin();
  double rapdilepton_max = background_rap->GetXaxis()->GetXmax();
  double mass_min = background_mass->GetXaxis()->GetXmin();
  double mass_max = background_mass->GetXaxis()->GetXmax();

	cout << "rapdilepton_min " << rapdilepton_min << endl;
	cout << "rapdilepton_max " << rapdilepton_max << endl;

  sprintf(filename,"%s/results.root",dirstruct);
  if(RealData) sprintf(filename,"%s/results_%s.root",dirstruct,TreeBinID);
  TFile* resultsFile = new TFile(filename, "RECREATE", "results");


  ////////////////////////////////////////////////////////////////////////////
  // calculate acceptance (wrt fiducial cuts) as a function of pT, y, mass

  const int nbinpT   = 10;
  const int nbinrap  = 10;
  const int nbinmass = 10;

  TH3D* acceptance = new TH3D( "acceptance", "", nbinpT,   pTdilepton_min,   pTdilepton_max,
                                                 nbinrap,  rapdilepton_min,  rapdilepton_max,
                                                 nbinmass, mass_min, mass_max  );

  // The pT-rap-mass dependence of the acceptance may be sensitive (at least at low pT)
  // to the polarization hypothesis, to be specified in "effsAndCuts.h"

  // Monte Carlo calculation of the acceptance

  int n_step = n_MCevents_ACC/5;  // to visualize progress of the event scan (50 steps)
  int n_step_=1;

  cout << endl;
  cout << "Calculation of acceptance vs pT, y and mass" << endl;
  cout << "(" << n_MCevents_ACC << " dilepton events)"<< endl;
  cout << "------------------------------------------------------------" << endl;
  cout << "Progress: "<<endl;


  ///////////////// cycle of MC events ////////////////////////
  for(int i_MCevent = 1; i_MCevent <= n_MCevents_ACC; i_MCevent++){

    if (i_MCevent%n_step == 0) {cout << n_step_*20 <<" % "<<endl; n_step_++;}

    // generation of dilepton in the pp CM, flat in pT, rapidity and mass

    // pT, rapidity, mass:

    double pT   = gRandom->Uniform(pTdilepton_min,  pTdilepton_max);
    double rap  = gRandom->Uniform(rapdilepton_min, rapdilepton_max);
    double mass = gRandom->Uniform(mass_min,        mass_max);

    // pL:

    double rap_sign = gRandom->Uniform(-1., 1.); rap_sign /= TMath::Abs(rap_sign);
    rap *= rap_sign;
    double mT = sqrt( mass*mass + pT*pT );
    double pL1 = 0.5 *mT * exp(rap);
    double pL2 = - 0.5 *mT * exp(-rap);
    double pL = pL1 + pL2;

    // Phi:

    double Phi = 2. * gPI_ * gRandom->Uniform(1.);

    // 4-vector:

    TLorentzVector dilepton;
    dilepton.SetXYZM( pT * cos(Phi) , pT * sin(Phi), pL, mass );


    // generation of polarization (generic reference frame)

    double lambda_theta    = lambda_theta_ACC;
    double lambda_phi      = lambda_phi_ACC;
    double lambda_thetaphi = lambda_thetaphi_ACC;
    bool HX_is_natural = HX_is_natural_ACC;
    bool PX_is_natural = PX_is_natural_ACC;

    double costhphidistr_max = 1. + TMath::Abs(lambda_phi) + TMath::Abs(lambda_thetaphi);
    double costhphidistr_rnd;
    double costhphidistr;
    double costh_gen;
    double sinth_gen;
    double phi_gen;

    if ( lambda_theta > 0. ) costhphidistr_max += lambda_theta;

    do { costh_gen = -1. + 2. * gRandom->Uniform(1.);
         phi_gen   = 2. * gPI_ * gRandom->Uniform(1.);
         sinth_gen = sqrt( 1. - costh_gen*costh_gen );
         costhphidistr_rnd = costhphidistr_max * gRandom->Uniform(1.);
         costhphidistr = 1. + lambda_theta    * costh_gen*costh_gen
                            + lambda_phi      * sinth_gen*sinth_gen * cos(2.*phi_gen)
                            + lambda_thetaphi * 2.* sinth_gen*costh_gen * cos(phi_gen);
       } while ( costhphidistr_rnd > costhphidistr );

    // lepton momentum in the dilepton rest frame:

    double p_lepton_DILEP = sqrt( 0.25*mass*mass - Mlepton_*Mlepton_ );

    TLorentzVector lepton_DILEP;

    lepton_DILEP.SetXYZM( p_lepton_DILEP * sinth_gen * cos(phi_gen),
                          p_lepton_DILEP * sinth_gen * sin(phi_gen),
                          p_lepton_DILEP * costh_gen,
                          Mlepton_ );

    // reference directions for leptons in the dilepton rest frame:

    TVector3 lab_to_dilep = -dilepton.BoostVector();

    TLorentzVector beam1_DILEP = beam1_LAB_;
    beam1_DILEP.Boost(lab_to_dilep);         // beam1 in the dilepton rest frame
    TLorentzVector beam2_DILEP = beam2_LAB_;
    beam2_DILEP.Boost(lab_to_dilep);         // beam2 in the dilepton rest frame

    TVector3 beam1_direction     = beam1_DILEP.Vect().Unit();
    TVector3 beam2_direction     = beam2_DILEP.Vect().Unit();
    TVector3 dilep_direction     = dilepton.Vect().Unit();
    TVector3 beam1_beam2_bisect  = ( beam1_direction - beam2_direction ).Unit();

    // Y axis = the normal to the plane formed by
    // the directions of the colliding hadrons

    TVector3 Yaxis = ( beam1_direction.Cross( beam2_direction ) ).Unit();

    // flip of y axis with rapidity

    if ( rap < 0 ) Yaxis = - Yaxis;

    TVector3 perpendicular_to_beam = ( beam1_beam2_bisect.Cross( Yaxis ) ).Unit();

    // transform (rotation) lepton momentum components from generation frame
    // to the frame with x,y,z axes as in the laboratory

    TVector3 oldZaxis = beam1_beam2_bisect;
    if ( HX_is_natural ) oldZaxis = dilep_direction;
    if ( PX_is_natural ) oldZaxis = perpendicular_to_beam;

    TVector3 oldYaxis = Yaxis;
    TVector3 oldXaxis = oldYaxis.Cross(oldZaxis);

    TRotation rotation;
    rotation.RotateAxes(oldXaxis, oldYaxis, oldZaxis);
                     // transforms coordinates from the "old" frame to the "xyz" frame

    TLorentzVector lepton_DILEP_xyz = lepton_DILEP;

    lepton_DILEP_xyz.Transform(rotation);
                     // lepton_DILEP_xyz is the lepton in the dilepton rest frame
                     // wrt to the lab axes

    // lepton 4-vectors in the LAB frame:

    TVector3 dilep_to_lab = dilepton.BoostVector();

    TLorentzVector lepP = lepton_DILEP_xyz;
    lepP.Boost(dilep_to_lab);
    TLorentzVector lepN = lepP;
    lepN.SetPxPyPzE(-lepton_DILEP_xyz.Px(),-lepton_DILEP_xyz.Py(),-lepton_DILEP_xyz.Pz(),lepton_DILEP_xyz.E());
    lepN.Boost(dilep_to_lab);

    double lepP_pT  = lepP.Pt();
    double lepN_pT  = lepN.Pt();

    double lepP_eta = lepP.PseudoRapidity();
    double lepN_eta = lepN.PseudoRapidity();


    // CS frame angles:

    TVector3 newZaxis = beam1_beam2_bisect;
    TVector3 newYaxis = Yaxis;
    TVector3 newXaxis = newYaxis.Cross( newZaxis );

    TRotation rotation2;
    rotation2.SetToIdentity();
    rotation2.RotateAxes( newXaxis, newYaxis, newZaxis );
    rotation2.Invert();  // transforms coordinates from the "xyz" frame to the new frame
    TVector3 lepton_DILEP_rotated = lepton_DILEP.Vect();
    lepton_DILEP_rotated.Transform(rotation2);

    double costh_CS = lepton_DILEP_rotated.CosTheta();
    double phi_CS   = lepton_DILEP_rotated.Phi() * 180. / gPI_;
    double phith_CS;
    if ( costh_CS < 0. ) phith_CS = phi_CS - 135.;
    if ( costh_CS > 0. ) phith_CS = phi_CS - 45.;
    if ( phith_CS < -180. ) phith_CS = 360. + phith_CS;


    // HELICITY frame angles:

    newZaxis = dilep_direction;
    newYaxis = Yaxis;
    newXaxis = newYaxis.Cross( newZaxis );

    rotation2.SetToIdentity();
    rotation2.RotateAxes( newXaxis, newYaxis, newZaxis );
    rotation2.Invert();
    lepton_DILEP_rotated = lepton_DILEP.Vect();
    lepton_DILEP_rotated.Transform(rotation2);

    double costh_HX = lepton_DILEP_rotated.CosTheta();
    double phi_HX   = lepton_DILEP_rotated.Phi() * 180. / gPI_;
    double phith_HX;
    if ( costh_HX < 0. ) phith_HX = phi_HX - 135.;
    if ( costh_HX > 0. ) phith_HX = phi_HX - 45.;
    if ( phith_HX < -180. ) phith_HX = 360. + phith_HX;

		//PhiHX test
		//if(PhiHX_test) { if(phi_HX>80. || phi_HX<91.) continue; }

    // PERPENDICULAR HELICITY frame angles:

    newZaxis = perpendicular_to_beam;
    newYaxis = Yaxis;
    newXaxis = newYaxis.Cross( newZaxis );

    rotation2.SetToIdentity();
    rotation2.RotateAxes( newXaxis, newYaxis, newZaxis );
    rotation2.Invert();
    lepton_DILEP_rotated = lepton_DILEP.Vect();
    lepton_DILEP_rotated.Transform(rotation2);

    double costh_PX = lepton_DILEP_rotated.CosTheta();
    double phi_PX   = lepton_DILEP_rotated.Phi() * 180. / gPI_;
    double phith_PX;
    if ( costh_PX < 0. ) phith_PX = phi_PX - 135.;
    if ( costh_PX > 0. ) phith_PX = phi_PX - 45.;
    if ( phith_PX < -180. ) phith_PX = 360. + phith_PX;

//    cout<<"costh_CS = "<<costh_CS<<", phi_CS = "<<phi_CS<<endl;
//    cout<<"costh_HX = "<<costh_HX<<", phi_HX = "<<phi_HX<<endl;
//    cout<<"costh_PX = "<<costh_PX<<", phi_PX = "<<phi_PX<<endl;

    //bool isEventAccepted = isMuonInAcceptance( FidCuts-1, lepP_pT, lepP_eta )
    //                    * isMuonInAcceptance( FidCuts-1, lepN_pT, lepN_eta ); 

		double deltaPhi = lepN.Phi() - lepP.Phi();
		if(deltaPhi > TMath::Pi()) deltaPhi -= 2.*TMath::Pi();
		else if(deltaPhi < -TMath::Pi()) deltaPhi += 2.*TMath::Pi();
		double deltaEta = lepN_eta - lepP_eta; 
		double deltaPT  = lepN_pT  - lepP_pT ;
		double deltaREll = TMath::Sqrt(deltaEta*deltaEta + TMath::Power(1.2*deltaPhi,2));
		double deltaREllDpt = TMath::Sqrt(deltaREll*deltaREll + TMath::Power(0.00157*deltaPT,2));

		bool accept = isMuonInAcceptance( FidCuts-1, lepP_pT, lepP_eta )
		 	* isMuonInAcceptance( FidCuts-1, lepN_pT, lepN_eta );
		bool isEventAccepted = false;
		if(accept && (!cutDeltaREllDpt || (cutDeltaREllDpt && deltaREllDpt > DeltaREllDptValue) )) isEventAccepted = true;

		//if(deltaREllDpt<DeltaREllDptValue)
		// 	cout<<"1: deltaREllDpt: "<<deltaREllDpt<<endl;
		//if(accept && !isEventAccepted){
		//	cout<<"accept: "<<accept<<endl;
		//	cout<<"isEventAccepted: "<<isEventAccepted<<endl;
		//}


    double epsilon = singleLeptonEfficiency( lepP_pT, lepP_eta, nEff, fInEff, hEvalEff, MCeff, TEff)
                   * singleLeptonEfficiency( lepN_pT, lepN_eta, nEff, fInEff, hEvalEff, MCeff, TEff);

    double costh_DILEff;
  //const double beta = 3.7;  // Psi(2S)
  //const double pTsq = 20.;  // Psi(2S)
    double phi_DILEff;
    if(nDileptonEff>200 && nDileptonEff<211) {costh_DILEff=costh_CS; phi_DILEff=phi_CS;}
    if(nDileptonEff>210 && nDileptonEff<221) {costh_DILEff=costh_HX; phi_DILEff=phi_HX;}
    if(nDileptonEff>220 && nDileptonEff<231) {costh_DILEff=costh_PX; phi_DILEff=phi_PX;}
    double DileptonEff = DiLeptonEfficiency( costh_DILEff, phi_DILEff, nDileptonEff, fInDileptonEff, MCDileptoneff);

    double costh_RhoFactor=0;
    double phi_RhoFactor=0;
    if(nRhoFactor>300 && nRhoFactor<311) {costh_RhoFactor=costh_CS; phi_RhoFactor=phi_CS;}
    if(nRhoFactor>310 && nRhoFactor<321) {costh_RhoFactor=costh_HX; phi_RhoFactor=phi_HX;}
    if(nRhoFactor>320 && nRhoFactor<331) {costh_RhoFactor=costh_PX; phi_RhoFactor=phi_PX;}
	double RhoFactor = EvaluateRhoFactor( costh_RhoFactor, phi_RhoFactor, nRhoFactor, fInRhoFactor, rap, pT, StatVarRho);

	//cout<<"epsilon: "<<epsilon<<endl;
	//cout<<"RhoFactor: "<<RhoFactor<<endl;


	if(!useAmapApproach) epsilon = epsilon*DileptonEff*RhoFactor;

	if(useAmapApproach){
	    double costh_Amap=0;
	    double phi_Amap=0;
	    if(nAmap>10000 && nAmap<20000) {costh_Amap=costh_CS; phi_Amap=phi_CS;}
	    if(nAmap>20000 && nAmap<30000) {costh_Amap=costh_HX; phi_Amap=phi_HX;}
	    if(nAmap>30000 && nAmap<40000) {costh_Amap=costh_PX; phi_Amap=phi_PX;}
		double AmapValue=EvaluateAmap( costh_Amap, phi_Amap, nAmap, fInAmap, rap, pT);
		double AmapDenominator= DenominatorAmapEfficiency( lepP_pT, lepP_eta, nDenominatorAmap, fInEff_nDenominatorAmap, hEvalEff_nDenominatorAmap, true, TEff_nDenominatorAmap) * DenominatorAmapEfficiency( lepN_pT, lepN_eta, nDenominatorAmap, fInEff_nDenominatorAmap, hEvalEff_nDenominatorAmap, true, TEff_nDenominatorAmap);
		epsilon = epsilon*AmapValue/AmapDenominator;

/*		cout<<"muonP pT = "<<lepP_pT<<" muonP eta = "<<lepP_eta<<endl;
		cout<<"muonN pT = "<<lepN_pT<<" muonN eta = "<<lepN_eta<<endl;
		cout<<"costh    = "<<costh_Amap<<" phi = "<<phi_Amap<<endl;
		cout<<"Ups pT = "<<pT<<" Ups rap = "<<rap<<endl;
		cout<<"muonP eff = "<<singleLeptonEfficiency( lepP_pT, lepP_eta, nEff, fInEff, hEvalEff, MCeff, TEff)<<endl;
		cout<<"muonN eff = "<<singleLeptonEfficiency( lepN_pT, lepN_eta, nEff, fInEff, hEvalEff, MCeff, TEff)<<endl;
		cout<<"muonP Denominator eff = "<<DenominatorAmapEfficiency( lepP_pT, lepP_eta, nDenominatorAmap, fInEff_nDenominatorAmap, hEvalEff_nDenominatorAmap, true, TEff_nDenominatorAmap)<<endl;
		cout<<"muonN Denominator eff = "<<DenominatorAmapEfficiency( lepN_pT, lepN_eta, nDenominatorAmap, fInEff_nDenominatorAmap, hEvalEff_nDenominatorAmap, true, TEff_nDenominatorAmap)<<endl;
		cout<<"AmapValue = "<<EvaluateAmap( costh_Amap, phi_Amap, nAmap, fInAmap, rap, pT)<<endl;
		cout<<"dimuon efficiency = "<<epsilon<<endl;
		cout<<"dimuon efficiency Check 0 = "<<epsilon-singleLeptonEfficiency( lepP_pT, lepP_eta, nEff, fInEff, hEvalEff, MCeff, TEff)*singleLeptonEfficiency( lepN_pT, lepN_eta, nEff, fInEff, hEvalEff, MCeff, TEff)*EvaluateAmap( costh_Amap, phi_Amap, nAmap, fInAmap, rap, pT)/(DenominatorAmapEfficiency( lepP_pT, lepP_eta, nDenominatorAmap, fInEff_nDenominatorAmap, hEvalEff_nDenominatorAmap, true, TEff_nDenominatorAmap)*DenominatorAmapEfficiency( lepN_pT, lepN_eta, nDenominatorAmap, fInEff_nDenominatorAmap, hEvalEff_nDenominatorAmap, true, TEff_nDenominatorAmap))<<endl;
*/

	}


	if(epsilon>1&&ForceEpsSmallerOne) {epsilon=1;}

    if ( isEventAccepted ) {
        if(!NewAccCalc) acceptance->Fill( pT, TMath::Abs(rap), mass );
        if(NewAccCalc) acceptance->Fill( pT, TMath::Abs(rap), mass ,epsilon );
    }

  } // end loop of MC acceptance calculation

  cout << endl;




  //////////////////////////////////////////////////////////////////////////////////
  // fill histograms with the kinematics of all events (signal+background)
  // for the subsequent background subtraction

  TH2D* total_costhphiPX  = (TH2D*)background_costhphiPX->Clone();
  total_costhphiPX->SetName("total_costhphiPX");
  total_costhphiPX->Reset();

  TH1D* total_costhPX  = (TH1D*)background_costhPX->Clone();
  total_costhPX->SetName("total_costhPX");
  total_costhPX->Reset();

  TH1D* total_phiPX    = (TH1D*)background_phiPX->Clone();
  total_phiPX->SetName("total_phiPX");
  total_phiPX->Reset();

  TH3D* total_pTrapMass = (TH3D*)background_pTrapMass->Clone();
  total_pTrapMass->SetName("total_pTrapMass");
  total_pTrapMass->Reset();

  TH1D* total_pT       = (TH1D*)background_pT->Clone();
  total_pT->SetName("total_pT");
  total_pT->Reset();

  TH1D* total_rap      = (TH1D*)background_rap->Clone();
  total_rap->SetName("total_rap");
  total_rap->Reset();

  TH1D* total_mass     = (TH1D*)background_mass->Clone();
  total_mass->SetName("total_mass");
  total_mass->Reset();


  int n_events = int( data->GetEntries() );

  cout << endl;
  cout << "Reading input ntuple (" << n_events << " dilepton events)" << endl;
  cout << "to fill preliminary histograms"<< endl;
  cout << "------------------------------------------------------------" << endl;
  cout << "Progress: "<<endl;

  n_step = n_events/5;  // to visualize progress of the event scan (50 steps)
  n_step_=1;


  // loop over dilepton events in the ntuple

  for ( int i_event = 1; i_event <= n_events; i_event++ ) {

    if (i_event%n_step == 0) {cout << n_step_*20 <<" % "<<endl; n_step_++;}

    data->GetEvent( i_event-1 );

    double lepP_pT  = lepP->Pt();
    double lepN_pT  = lepN->Pt();

    double lepP_eta = lepP->PseudoRapidity();
    double lepN_eta = lepN->PseudoRapidity();

    // dilepton 4-vector:

    TLorentzVector dilepton = *lepP + *lepN;
    double pT   = dilepton.Pt();
    double rap  = dilepton.Rapidity();
    double mass = dilepton.M();

    // calculation of decay angles in three polarization frames

    // reference directions to calculate angles:

    TVector3 lab_to_dilep = -dilepton.BoostVector();

    TLorentzVector beam1_DILEP = beam1_LAB_;
    beam1_DILEP.Boost(lab_to_dilep);         // beam1 in the dilepton rest frame
    TLorentzVector beam2_DILEP = beam2_LAB_;
    beam2_DILEP.Boost(lab_to_dilep);         // beam2 in the dilepton rest frame

    TVector3 beam1_direction     = beam1_DILEP.Vect().Unit();
    TVector3 beam2_direction     = beam2_DILEP.Vect().Unit();
    TVector3 dilep_direction     = dilepton.Vect().Unit();
    TVector3 beam1_beam2_bisect  = ( beam1_direction - beam2_direction ).Unit();


    // all polarization frames have the same Y axis = the normal to the plane formed by
    // the directions of the colliding hadrons:

    TVector3 Yaxis = ( beam1_direction.Cross( beam2_direction ) ).Unit();

    // flip of y axis with rapidity:

    if ( rap < 0. ) Yaxis = - Yaxis;

    TVector3 perpendicular_to_beam = ( beam1_beam2_bisect.Cross( Yaxis ) ).Unit();


    // positive lepton in the dilepton rest frame:

    TLorentzVector lepton_DILEP = *lepP;
    lepton_DILEP.Boost(lab_to_dilep);

    // CS frame angles:

    TVector3 newZaxis = beam1_beam2_bisect;
    TVector3 newYaxis = Yaxis;
    TVector3 newXaxis = newYaxis.Cross( newZaxis );

    TRotation rotation;
    rotation.SetToIdentity();
    rotation.RotateAxes( newXaxis, newYaxis, newZaxis );
    rotation.Invert();  // transforms coordinates from the "xyz" frame to the new frame
    TVector3 lepton_DILEP_rotated = lepton_DILEP.Vect();
    lepton_DILEP_rotated.Transform(rotation);

    double costh_CS = lepton_DILEP_rotated.CosTheta();
    double phi_CS   = lepton_DILEP_rotated.Phi() * 180. / gPI_;
    double phith_CS;
    if ( costh_CS < 0. ) phith_CS = phi_CS - 135.;
    if ( costh_CS > 0. ) phith_CS = phi_CS - 45.;
    if ( phith_CS < -180. ) phith_CS = 360. + phith_CS;


    // HELICITY frame angles:

    newZaxis = dilep_direction;
    newYaxis = Yaxis;
    newXaxis = newYaxis.Cross( newZaxis );

    rotation.SetToIdentity();
    rotation.RotateAxes( newXaxis, newYaxis, newZaxis );
    rotation.Invert();
    lepton_DILEP_rotated = lepton_DILEP.Vect();
    lepton_DILEP_rotated.Transform(rotation);

    double costh_HX = lepton_DILEP_rotated.CosTheta();
    double phi_HX   = lepton_DILEP_rotated.Phi() * 180. / gPI_;
    double phith_HX;
    if ( costh_HX < 0. ) phith_HX = phi_HX - 135.;
    if ( costh_HX > 0. ) phith_HX = phi_HX - 45.;
    if ( phith_HX < -180. ) phith_HX = 360. + phith_HX;

		//PhiHX test
		//if(PhiHX_test) { if(phi_HX>80. || phi_HX<91.) continue; }

    // PERPENDICULAR HELICITY frame angles:

    newZaxis = perpendicular_to_beam;
    newYaxis = Yaxis;
    newXaxis = newYaxis.Cross( newZaxis );

    rotation.SetToIdentity();
    rotation.RotateAxes( newXaxis, newYaxis, newZaxis );
    rotation.Invert();
    lepton_DILEP_rotated = lepton_DILEP.Vect();
    lepton_DILEP_rotated.Transform(rotation);

    double costh_PX = lepton_DILEP_rotated.CosTheta();
    double phi_PX   = lepton_DILEP_rotated.Phi() * 180. / gPI_;
    double phith_PX;
    if ( costh_PX < 0. ) phith_PX = phi_PX - 135.;
    if ( costh_PX > 0. ) phith_PX = phi_PX - 45.;
    if ( phith_PX < -180. ) phith_PX = 360. + phith_PX;


    // efficiency:

    double epsilon = singleLeptonEfficiency( lepP_pT, lepP_eta, nEff, fInEff, hEvalEff, MCeff, TEff)
                   * singleLeptonEfficiency( lepN_pT, lepN_eta, nEff, fInEff, hEvalEff, MCeff, TEff);

    double costh_DILEff;
    double phi_DILEff;
    if(nDileptonEff>200 && nDileptonEff<211) {costh_DILEff=costh_CS; phi_DILEff=phi_CS;}
    if(nDileptonEff>210 && nDileptonEff<221) {costh_DILEff=costh_HX; phi_DILEff=phi_HX;}
    if(nDileptonEff>220 && nDileptonEff<231) {costh_DILEff=costh_PX; phi_DILEff=phi_PX;}
    double DileptonEff = DiLeptonEfficiency( costh_DILEff, phi_DILEff, nDileptonEff, fInDileptonEff, MCDileptoneff);

    double costh_RhoFactor=0;
    double phi_RhoFactor=0;
    if(nRhoFactor>300 && nRhoFactor<311) {costh_RhoFactor=costh_CS; phi_RhoFactor=phi_CS;}
    if(nRhoFactor>310 && nRhoFactor<321) {costh_RhoFactor=costh_HX; phi_RhoFactor=phi_HX;}
    if(nRhoFactor>320 && nRhoFactor<331) {costh_RhoFactor=costh_PX; phi_RhoFactor=phi_PX;}
	double RhoFactor = EvaluateRhoFactor( costh_RhoFactor, phi_RhoFactor, nRhoFactor, fInRhoFactor, rap, pT, StatVarRho);


	if(!useAmapApproach) epsilon = epsilon*DileptonEff*RhoFactor;

	if(useAmapApproach){
	    double costh_Amap=0;
	    double phi_Amap=0;
	    if(nAmap>10000 && nAmap<20000) {costh_Amap=costh_CS; phi_Amap=phi_CS;}
	    if(nAmap>20000 && nAmap<30000) {costh_Amap=costh_HX; phi_Amap=phi_HX;}
	    if(nAmap>30000 && nAmap<40000) {costh_Amap=costh_PX; phi_Amap=phi_PX;}
		double AmapValue=EvaluateAmap( costh_Amap, phi_Amap, nAmap, fInAmap, rap, pT);
		double AmapDenominator= DenominatorAmapEfficiency( lepP_pT, lepP_eta, nDenominatorAmap, fInEff_nDenominatorAmap, hEvalEff_nDenominatorAmap, true, TEff_nDenominatorAmap) * DenominatorAmapEfficiency( lepN_pT, lepN_eta, nDenominatorAmap, fInEff_nDenominatorAmap, hEvalEff_nDenominatorAmap, true, TEff_nDenominatorAmap);
		epsilon = epsilon*AmapValue/AmapDenominator;
	}


	if(epsilon>1&&ForceEpsSmallerOne) {epsilon=1;}

    // acceptance:

    //bool isEventAccepted = isMuonInAcceptance( FidCuts-1, lepP_pT, lepP_eta )
    //                     * isMuonInAcceptance( FidCuts-1, lepN_pT, lepN_eta );
		
		double deltaPhi = lepN->Phi() - lepP->Phi();
		if(deltaPhi > TMath::Pi()) deltaPhi -= 2.*TMath::Pi();
		else if(deltaPhi < -TMath::Pi()) deltaPhi += 2.*TMath::Pi();
		double deltaEta = lepN_eta - lepP_eta; 
		double deltaPT  = lepN_pT  - lepP_pT ;
		double deltaREll = TMath::Sqrt(deltaEta*deltaEta + TMath::Power(1.2*deltaPhi,2));
		double deltaREllDpt = TMath::Sqrt(deltaREll*deltaREll + TMath::Power(0.00157*deltaPT,2));

		bool accept = isMuonInAcceptance( FidCuts-1, lepP_pT, lepP_eta )
		 	* isMuonInAcceptance( FidCuts-1, lepN_pT, lepN_eta );
		bool isEventAccepted = false;
		if(accept && (!cutDeltaREllDpt || (cutDeltaREllDpt && deltaREllDpt > DeltaREllDptValue) )) isEventAccepted = true;

		//if(deltaREllDpt<DeltaREllDptValue)
		// 	cout<<"2: deltaREllDpt: "<<deltaREllDpt<<endl;
		//if(accept && !isEventAccepted){
		//	cout<<"accept: "<<accept<<endl;
		//	cout<<"isEventAccepted: "<<isEventAccepted<<endl;
		//}

    if ( epsilon > min_dileptonEff && isEventAccepted ) {

      // fill histograms
      total_costhphiPX->Fill( costh_PX, phi_PX, 1. );
      total_costhPX->Fill( costh_PX, 1. );
      total_phiPX->Fill( phi_PX, 1. );
      total_pTrapMass->Fill( pT, TMath::Abs( rap ), mass, 1. );
      total_pT->Fill( pT, 1. );
      total_rap->Fill( TMath::Abs( rap ), 1. );
      total_mass->Fill( mass, 1. );
    }

  } // end of loop over dilepton events in the ntuple

  cout << endl;

  // Normalize histograms

  background_costhphiPX->Scale( f_background*n_events / background_costhphiPX->Integral() );
  background_costhPX->Scale( f_background*n_events / background_costhPX->Integral() );
  background_phiPX->Scale( f_background*n_events / background_phiPX->Integral() );

  background_pTrapMass->Scale( f_background*n_events / background_pTrapMass->Integral() );
  background_pT->Scale( f_background*n_events / background_pT->Integral() );
  background_rap->Scale( f_background*n_events / background_rap->Integral() );
  background_mass->Scale( f_background*n_events / background_mass->Integral() );


  total_costhphiPX->Scale( n_events / total_costhphiPX->Integral() );
  total_costhPX->Scale( n_events / total_costhPX->Integral() );
  total_phiPX->Scale( n_events / total_phiPX->Integral() );

  total_pTrapMass->Scale( n_events / total_pTrapMass->Integral() );
  total_pT->Scale( n_events / total_pT->Integral() );
  total_rap->Scale( n_events / total_rap->Integral() );
  total_mass->Scale( n_events / total_mass->Integral() );



  ////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // define histograms with |y|-pT-mass background-subtracted and efficiency*acceptance-corrected
  // distributions of the dileptons, to generate events in the subsequent MC integration of the likelihood

  TH3D* pTrapMass_gen = new TH3D( "pTrapMass_gen", "", nbinpT,   pTdilepton_min,  pTdilepton_max,
                                                       nbinrap,  rapdilepton_min, rapdilepton_max,
                                                       nbinmass, mass_min,        mass_max        );
  pTrapMass_gen->Sumw2();

  TH1D* AccEff = new TH1D( "AccEff", "AccEff", 100,   0,  1000);
  TH1D* Acc = new TH1D( "Acc", "Acc", 100,   0,  1000);
  TH1D* Eff = new TH1D( "Eff", "Eff", 100,   0,  1);


  //////////////////////////////////////////////////////////////////////////////
  // create output file with ntuples of results
  // (distributions of angles and parameters)


  // output ntuple with angular distribution of background-subtracted events

  TTree* angles = new TTree("angles","angles");

  double costh_CS;  angles->Branch("costh_CS",     &costh_CS,     "costh_CS/D");
  double phi_CS;    angles->Branch("phi_CS",       &phi_CS,       "phi_CS/D"  );
  double phith_CS;  angles->Branch("phithCS",      &phith_CS,     "phith_CS/D");

  double costh_HX;  angles->Branch("costh_HX",     &costh_HX,     "costh_HX/D");
  double phi_HX;    angles->Branch("phi_HX",       &phi_HX,       "phi_HX/D"  );
  double phith_HX;  angles->Branch("phith_HX",     &phith_HX,     "phith_HX/D");

  double costh_PX;  angles->Branch("costh_PX",     &costh_PX,     "costh_PX/D");
  double phi_PX;    angles->Branch("phi_PX",       &phi_PX,       "phi_PX/D"  );
  double phith_PX;  angles->Branch("phith_PX",     &phith_PX,     "phith_PX/D");

  double cosalpha;  angles->Branch("cosalpha",     &cosalpha,     "cosalpha/D");

  double epsilon;   angles->Branch("epsilon",      &epsilon,      "epsilon/D");


  //////////////////////////////////////////////////////////////////////////////////////////
  // Read all events of the data ntuple and
  // -- subtract background
  // -- calculate the parameter-independent pieces of the event PDF and store them
  //    in large arrays, so that after this loop the ntuple will not be accessed anymore
  // -- create output ntuple with decay angles in all frames for the background-subtracted
  //    events and the efficiency weight
  // -- fill pT, y and mass histograms with efficiency*acceptance corrected events

  cout << endl;
  cout << "Reading input ntuple (" << n_events << " dilepton events)," << endl;
  cout << "subtracting background"<< endl;
  cout << "and calculating parameter-independent PDF terms"<< endl;
  cout << "------------------------------------------------------------" << endl;
  cout << "Progress: "<< endl;


  int i_signalEvent = 0;
  int i_EffRejectedSignalEvent=0;

  n_step = n_events/5;  // to visualize progress of the event scan (50 steps)
  n_step_=1;


  // arrays to be filled, one element per event

  double* w_1           = new double[n_events];

  double* w_cth2_CS     = new double[n_events];
  double* w_sth2c2ph_CS = new double[n_events];
  double* w_s2thcph_CS  = new double[n_events];

  double* w_cth2_HX     = new double[n_events];
  double* w_sth2c2ph_HX = new double[n_events];
  double* w_s2thcph_HX  = new double[n_events];

  double* w_cth2_PX     = new double[n_events];
  double* w_sth2c2ph_PX = new double[n_events];
  double* w_s2thcph_PX  = new double[n_events];


  // backgrond TEST distributions (go into output file)
  // should be flat (possible border fluctuations - they are divisions)
  // indicating that the background distributions are well reproduced
  // (and subtracted) by the likelihood-ratio method

  TH1D* background_costhPX_test  = (TH1D*)background_costhPX->Clone();
  background_costhPX_test->SetName("background_costhPX_test");
  background_costhPX_test->Reset();

  TH1D* background_phiPX_test    = (TH1D*)background_phiPX->Clone();
  background_phiPX_test->SetName("background_phiPX_test");
  background_phiPX_test->Reset();

  TH1D* background_pT_test  = (TH1D*)background_pT->Clone();
  background_pT_test->SetName("background_pT_test");
  background_pT_test->Reset();

  TH1D* background_rap_test  = (TH1D*)background_rap->Clone();
  background_rap_test->SetName("background_rap_test");
  background_rap_test->Reset();

  TH1D* background_mass_test  = (TH1D*)background_mass->Clone();
  background_mass_test->SetName("background_mass_test");
  background_mass_test->Reset();

  TH1D*  SubtractedBG_test   = new TH1D( "SubtractedBG_test",       "", 100, 0.7, 1.3 );

  int n_background_subtracted = 0;

  ///////////////////////////////////////////////////////////////////////////////////////////////
  // loop over dilepton events in the ntuple

  for ( int i_event = 1; i_event <= n_events; i_event++ ) {


	    if (i_event%n_step == 0) {cout << n_step_*20 <<" % "<<endl; n_step_++;}

    data->GetEvent( i_event-1 );

    double lepP_pT  = lepP->Pt();
    double lepN_pT  = lepN->Pt();

    double lepP_eta = lepP->PseudoRapidity();
    double lepN_eta = lepN->PseudoRapidity();

    // dilepton 4-vector:

    TLorentzVector dilepton = *lepP + *lepN;
    double pT   = dilepton.Pt();
    double rap  = dilepton.Rapidity();
    double mass = dilepton.M();

    // calculation of decay angles in three polarization frames

    // reference directions to calculate angles:

    TVector3 lab_to_dilep = -dilepton.BoostVector();

    TLorentzVector beam1_DILEP = beam1_LAB_;
    beam1_DILEP.Boost(lab_to_dilep);         // beam1 in the dilepton rest frame
    TLorentzVector beam2_DILEP = beam2_LAB_;
    beam2_DILEP.Boost(lab_to_dilep);         // beam2 in the dilepton rest frame

    TVector3 beam1_direction     = beam1_DILEP.Vect().Unit();
    TVector3 beam2_direction     = beam2_DILEP.Vect().Unit();
    TVector3 dilep_direction     = dilepton.Vect().Unit();
    TVector3 beam1_beam2_bisect  = ( beam1_direction - beam2_direction ).Unit();

    // all polarization frames have the same Y axis = the normal to the plane formed by
    // the directions of the colliding hadrons:

    TVector3 Yaxis = ( beam1_direction.Cross( beam2_direction ) ).Unit();

    // flip of y axis with rapidity:

    if ( rap < 0. ) Yaxis = - Yaxis;

    TVector3 perpendicular_to_beam = ( beam1_beam2_bisect.Cross( Yaxis ) ).Unit();

    // positive lepton in the dilepton rest frame:

    TLorentzVector lepton_DILEP = *lepP;
    lepton_DILEP.Boost(lab_to_dilep);


    // CS frame angles:

    TVector3 newZaxis = beam1_beam2_bisect;
    TVector3 newYaxis = Yaxis;
    TVector3 newXaxis = newYaxis.Cross( newZaxis );

    TRotation rotation;
    rotation.SetToIdentity();
    rotation.RotateAxes( newXaxis, newYaxis, newZaxis );
    rotation.Invert();  // transforms coordinates from the "xyz" frame to the new frame
    TVector3 lepton_DILEP_rotated = lepton_DILEP.Vect();
    lepton_DILEP_rotated.Transform(rotation);

    costh_CS = lepton_DILEP_rotated.CosTheta();
    phi_CS   = lepton_DILEP_rotated.Phi() * 180. / gPI_;
    if ( costh_CS < 0. ) phith_CS = phi_CS - 135.;
    if ( costh_CS > 0. ) phith_CS = phi_CS - 45.;
    if ( phith_CS < -180. ) phith_CS = 360. + phith_CS;


    // HELICITY frame angles:

    newZaxis = dilep_direction;
    newYaxis = Yaxis;
    newXaxis = newYaxis.Cross( newZaxis );

    rotation.SetToIdentity();
    rotation.RotateAxes( newXaxis, newYaxis, newZaxis );
    rotation.Invert();
    lepton_DILEP_rotated = lepton_DILEP.Vect();
    lepton_DILEP_rotated.Transform(rotation);

    costh_HX = lepton_DILEP_rotated.CosTheta();
    phi_HX   = lepton_DILEP_rotated.Phi() * 180. / gPI_;
    if ( costh_HX < 0. ) phith_HX = phi_HX - 135.;
    if ( costh_HX > 0. ) phith_HX = phi_HX - 45.;
    if ( phith_HX < -180. ) phith_HX = 360. + phith_HX;

		//PhiHX test
		//if(PhiHX_test) { if(phi_HX>80. || phi_HX<91.) continue; }

    // PERPENDICULAR HELICITY frame angles:

    newZaxis = perpendicular_to_beam;
    newYaxis = Yaxis;
    newXaxis = newYaxis.Cross( newZaxis );

    rotation.SetToIdentity();
    rotation.RotateAxes( newXaxis, newYaxis, newZaxis );
    rotation.Invert();
    lepton_DILEP_rotated = lepton_DILEP.Vect();
    lepton_DILEP_rotated.Transform(rotation);

    costh_PX = lepton_DILEP_rotated.CosTheta();
    phi_PX   = lepton_DILEP_rotated.Phi() * 180. / gPI_;
    if ( costh_PX < 0. ) phith_PX = phi_PX - 135.;
    if ( costh_PX > 0. ) phith_PX = phi_PX - 45.;
    if ( phith_PX < -180. ) phith_PX = 360. + phith_PX;


    // invariant polarization angle

    cosalpha = sqrt( 1. - pow(costh_PX, 2.) ) * sin( lepton_DILEP_rotated.Phi() );



    // background subtraction as a function of costh and phi: reject events having distributions
    // compatible with the corresponding background distributions.

    int ibin_costhPX = background_costhPX->GetXaxis()->FindBin( costh_PX );
    int ibin_phiPX   = background_phiPX->GetXaxis()->FindBin( phi_PX );
    int ibin_pT      = background_pT->GetXaxis()->FindBin( pT );
    int ibin_rap     = background_rap->GetXaxis()->FindBin( TMath::Abs( rap ) );
    int ibin_mass    = background_mass->GetXaxis()->FindBin( mass );

    // background likelihood

    double likelihood_BG  = background_costhphiPX->GetBinContent( ibin_costhPX, ibin_phiPX )
                          * background_pTrapMass->GetBinContent( ibin_pT, ibin_rap, ibin_mass );
		
	  // apply statistical fluctuations on likelihood_BG with poisson distribution
		if(StatVarTotBGmodel){
			double binContent_costhphiPX= background_costhphiPX->GetBinContent( ibin_costhPX, ibin_phiPX );
			double binError_costhphiPX  = background_costhphiPX->GetBinError  ( ibin_costhPX, ibin_phiPX );
			double binContent_pTrapMass = background_pTrapMass->GetBinContent( ibin_pT, ibin_rap, ibin_mass );
			double binError_pTrapMass   = background_pTrapMass->GetBinError   ( ibin_pT, ibin_rap, ibin_mass );
			//cout << "binContent_costhphiPX:" << binContent_costhphiPX << " binError_costhphiPX:" << binError_costhphiPX << endl;
			//cout << "binContent_pTrapMass:" << binContent_pTrapMass << " binError_pTrapMass:" << binError_pTrapMass << endl;


			double nPoisson_costhphiPX  = TMath::Power( binContent_costhphiPX / binError_costhphiPX , 2 );
			double nPoisson_pTrapMass   = TMath::Power( binContent_pTrapMass  / binError_pTrapMass  , 2 );
			//cout << "nPoisson_costhphiPX: " << nPoisson_costhphiPX << " nPoisson_pTrapMass:" << nPoisson_pTrapMass << endl;

			double rnd_costhphiPX = gRandom -> Poisson( nPoisson_costhphiPX );
			double rnd_pTrapMass  = gRandom -> Poisson( nPoisson_pTrapMass  );
			//cout << "rnd_costhphiPX: " << rnd_costhphiPX << " rnd_pTrapMass:" << rnd_pTrapMass << endl;

			likelihood_BG = likelihood_BG * rnd_costhphiPX * rnd_pTrapMass / ( nPoisson_costhphiPX * nPoisson_pTrapMass ) ;

		}

    // total likelihood

    double likelihood_SpB = total_costhphiPX->GetBinContent( ibin_costhPX, ibin_phiPX )
                          * total_pTrapMass->GetBinContent( ibin_pT, ibin_rap, ibin_mass );

    if (    likelihood_BG > likelihood_SpB * f_background
         || likelihood_BG > likelihood_SpB * gRandom->Uniform( f_background )  ) {

       ++n_background_subtracted;
       background_costhPX_test->Fill(costh_PX);
       background_phiPX_test->Fill(phi_PX);
       background_pT_test->Fill( pT );
       background_rap_test->Fill( TMath::Abs(rap) );
       background_mass_test->Fill( mass );
       //if(subtractBG)
       continue;
    }


    // efficiency:

    epsilon = singleLeptonEfficiency( lepP_pT, lepP_eta, nEff, fInEff, hEvalEff, MCeff, TEff)
            * singleLeptonEfficiency( lepN_pT, lepN_eta, nEff, fInEff, hEvalEff, MCeff, TEff);

    double costh_DILEff;
    double phi_DILEff;
    if(nDileptonEff>200 && nDileptonEff<211) {costh_DILEff=costh_CS; phi_DILEff=phi_CS;}
    if(nDileptonEff>210 && nDileptonEff<221) {costh_DILEff=costh_HX; phi_DILEff=phi_HX;}
    if(nDileptonEff>220 && nDileptonEff<231) {costh_DILEff=costh_PX; phi_DILEff=phi_PX;}
    double DileptonEff = DiLeptonEfficiency( costh_DILEff, phi_DILEff, nDileptonEff, fInDileptonEff , MCDileptoneff);

    double costh_RhoFactor=0;
    double phi_RhoFactor=0;
    if(nRhoFactor>300 && nRhoFactor<311) {costh_RhoFactor=costh_CS; phi_RhoFactor=phi_CS;}
    if(nRhoFactor>310 && nRhoFactor<321) {costh_RhoFactor=costh_HX; phi_RhoFactor=phi_HX;}
    if(nRhoFactor>320 && nRhoFactor<331) {costh_RhoFactor=costh_PX; phi_RhoFactor=phi_PX;}
	double RhoFactor = EvaluateRhoFactor( costh_RhoFactor, phi_RhoFactor, nRhoFactor, fInRhoFactor, rap, pT, StatVarRho);


	if(!useAmapApproach) epsilon = epsilon*DileptonEff*RhoFactor;

	if(useAmapApproach){
//		cout<<"useAmapApproach"<<endl;
	    double costh_Amap=0;
	    double phi_Amap=0;
	    if(nAmap>10000 && nAmap<20000) {costh_Amap=costh_CS; phi_Amap=phi_CS;}
	    if(nAmap>20000 && nAmap<30000) {costh_Amap=costh_HX; phi_Amap=phi_HX;}
	    if(nAmap>30000 && nAmap<40000) {costh_Amap=costh_PX; phi_Amap=phi_PX;}
//		cout<<"before EvalAmap"<<endl;
		double AmapValue=EvaluateAmap( costh_Amap, phi_Amap, nAmap, fInAmap, rap, pT);
//		cout<<"AmapValue "<<AmapValue<<endl;
//		cout<<"after EvalAmap"<<endl;
		double AmapDenominator= DenominatorAmapEfficiency( lepP_pT, lepP_eta, nDenominatorAmap, fInEff_nDenominatorAmap, hEvalEff_nDenominatorAmap, true, TEff_nDenominatorAmap) * DenominatorAmapEfficiency( lepN_pT, lepN_eta, nDenominatorAmap, fInEff_nDenominatorAmap, hEvalEff_nDenominatorAmap, true, TEff_nDenominatorAmap);
		epsilon = epsilon*AmapValue/AmapDenominator;
	}


	if(epsilon>1&&ForceEpsSmallerOne) {epsilon=1;}

    // acceptance:

    //bool isEventAccepted = isMuonInAcceptance( FidCuts-1, lepP_pT, lepP_eta )
    //                     * isMuonInAcceptance( FidCuts-1, lepN_pT, lepN_eta );

		double deltaPhi = lepN->Phi() - lepP->Phi();
		if(deltaPhi > TMath::Pi()) deltaPhi -= 2.*TMath::Pi();
		else if(deltaPhi < -TMath::Pi()) deltaPhi += 2.*TMath::Pi();
		double deltaEta = lepN_eta - lepP_eta; 
		double deltaPT  = lepN_pT  - lepP_pT ;
		double deltaREll = TMath::Sqrt(deltaEta*deltaEta + TMath::Power(1.2*deltaPhi,2));
		double deltaREllDpt = TMath::Sqrt(deltaREll*deltaREll + TMath::Power(0.00157*deltaPT,2));

		bool accept = isMuonInAcceptance( FidCuts-1, lepP_pT, lepP_eta )
		 	* isMuonInAcceptance( FidCuts-1, lepN_pT, lepN_eta );
		bool isEventAccepted = false;
		if(accept && (!cutDeltaREllDpt || (cutDeltaREllDpt && deltaREllDpt > DeltaREllDptValue) )) isEventAccepted = true;

		//if(deltaREllDpt<DeltaREllDptValue)
		// 	cout<<"3: deltaREllDpt: "<<deltaREllDpt<<endl;
		//if(accept && !isEventAccepted){
		//	cout<<"accept: "<<accept<<endl;
		//	cout<<"isEventAccepted: "<<isEventAccepted<<endl;
		//}

//    double additionalPtCut=0;
//    if( lepP_pT < additionalPtCut || lepN_pT < additionalPtCut ) continue;

    if ( epsilon > min_dileptonEff && isEventAccepted ) {

      ++i_signalEvent;

      // store calculated values in the arrays

      w_1[i_signalEvent]           = epsilon;

      w_cth2_CS[i_signalEvent]     = epsilon * costh_CS*costh_CS;
      w_sth2c2ph_CS[i_signalEvent] = epsilon * (1. - costh_CS*costh_CS) * cos(gPI_/90.*phi_CS);
      w_s2thcph_CS[i_signalEvent]  = epsilon * 2.*costh_CS*sqrt(1.-costh_CS*costh_CS) * cos(gPI_/180.*phi_CS);

      w_cth2_HX[i_signalEvent]     = epsilon * costh_HX*costh_HX;
      w_sth2c2ph_HX[i_signalEvent] = epsilon * (1. - costh_HX*costh_HX) * cos(gPI_/90.*phi_HX);
      w_s2thcph_HX[i_signalEvent]  = epsilon * 2.*costh_HX*sqrt(1.-costh_HX*costh_HX) * cos(gPI_/180.*phi_HX);

      w_cth2_PX[i_signalEvent]     = epsilon * costh_PX*costh_PX;
      w_sth2c2ph_PX[i_signalEvent] = epsilon * (1. - costh_PX*costh_PX) * cos(gPI_/90.*phi_PX);
      w_s2thcph_PX[i_signalEvent]  = epsilon * 2.*costh_PX*sqrt(1.-costh_PX*costh_PX) * cos(gPI_/180.*phi_PX);


      // fill ntuple with decay angles

      angles->Fill();

      // fill pT, y and mass efficiency*acceptance-corrected histograms

      int ibin_pT   = acceptance->GetXaxis()->FindBin( pT );
      int ibin_rap  = acceptance->GetYaxis()->FindBin( TMath::Abs(rap) );
      int ibin_mass = acceptance->GetZaxis()->FindBin( mass );

      double acc  = acceptance->GetBinContent( ibin_pT, ibin_rap, ibin_mass );

      double effTimesAcc;
      effTimesAcc = epsilon * acc;
      if(NewAccCalc) effTimesAcc = acc;

      if ( effTimesAcc > 0.01 ) {
    	  pTrapMass_gen->Fill( pT, TMath::Abs( rap ), mass, 1./effTimesAcc );
      }
      AccEff->Fill(effTimesAcc);
      Acc->Fill(acc);
      Eff->Fill(epsilon);


    }
    if(epsilon < min_dileptonEff){i_EffRejectedSignalEvent++;cout<<"Event rejected because of dimuonEff<1%"<<endl;}
  } // end of loop over dilepton events in the ntuple

  cout << endl;


  double n_signalEvents = (double)i_signalEvent;

  cout<<"Fraction of rejected signal events due to eff<1%: "<<
		(double)i_EffRejectedSignalEvent/(n_signalEvents+(double)i_EffRejectedSignalEvent)<<endl;

  double f_background_actual = n_background_subtracted / double( n_events );
	cout << "n_events: " << n_events << " n_background_subtracted: " << n_background_subtracted<<endl;
  cout << "Background fraction: " << f_background << ", subtracted: " << f_background_actual << endl;
	cout << "Ratio: " << f_background_actual/f_background << endl;


  SubtractedBG_test->Fill(f_background_actual/f_background);

  background_costhPX_test->Divide(background_costhPX);
  background_phiPX_test->Divide(background_phiPX);

  background_mass_test->Divide(background_mass);
  background_pT_test->Divide(background_pT);
  background_rap_test->Divide(background_rap);


  ////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////
  ////// Monte Carlo integration of the parameter-independent pieces of the PDF
  ////// *within the acceptance*, using the efficiency-corrected pT-y distribution

  n_step = n_MCevents_NORM/5;
  n_step_=1;

  cout << endl;
  cout << "MC integration of the parameter-inpependent pieces" << endl;
  cout << "(" << n_MCevents_NORM << " dilepton events)"<< endl;
  cout << "------------------------------------------------------------" << endl;
  cout << "Progress: "<<endl;


  double sum_1           = 0.;

  double sum_cth2_CS     = 0.;
  double sum_sth2c2ph_CS = 0.;
  double sum_s2thcph_CS  = 0.;

  double sum_cth2_HX     = 0.;
  double sum_sth2c2ph_HX = 0.;
  double sum_s2thcph_HX  = 0.;

  double sum_cth2_PX     = 0.;
  double sum_sth2c2ph_PX = 0.;
  double sum_s2thcph_PX  = 0.;

  int n_acceptedEvents = 0;


  // histograms with the angular dependencies of the parameter-independent PDF terms,
  // to calculate the event PDFs to be compared with the data distributions in the final plots

  int nbin_cth = 80.;
  int nbin_ph  = 72.;

  TH1D*  PDF_1_vs_cth_CS   = new TH1D( "PDF_1_vs_cth_CS",       "", nbin_cth, -1., 1. );
  TH1D*  PDF_1_vs_ph_CS    = new TH1D( "PDF_1_vs_ph_CS",        "", nbin_ph, -180., 180. );
  TH1D*  PDF_1_vs_phth_CS  = new TH1D( "PDF_1_vs_phth_CS",      "", nbin_ph, -180., 180. );

  TH1D*  PDF_cth2_vs_cth_CS  = new TH1D( "PDF_cth2_vs_cth_CS",       "", nbin_cth, -1., 1. );
  TH1D*  PDF_cth2_vs_ph_CS   = new TH1D( "PDF_cth2_vs_ph_CS",        "", nbin_ph, -180., 180. );
  TH1D*  PDF_cth2_vs_phth_CS = new TH1D( "PDF_cth2_vs_phth_CS",      "", nbin_ph, -180., 180. );

  TH1D*  PDF_sth2c2ph_vs_cth_CS  = new TH1D( "PDF_sth2c2ph_vs_cth_CS",       "", nbin_cth, -1., 1. );
  TH1D*  PDF_sth2c2ph_vs_ph_CS   = new TH1D( "PDF_sth2c2ph_vs_ph_CS",        "", nbin_ph, -180., 180. );
  TH1D*  PDF_sth2c2ph_vs_phth_CS = new TH1D( "PDF_sth2c2ph_vs_phth_CS",      "", nbin_ph, -180., 180. );

  TH1D*  PDF_s2thcph_vs_cth_CS  = new TH1D( "PDF_s2thcph_vs_cth_CS",       "", nbin_cth, -1., 1. );
  TH1D*  PDF_s2thcph_vs_ph_CS   = new TH1D( "PDF_s2thcph_vs_ph_CS",        "", nbin_ph, -180., 180. );
  TH1D*  PDF_s2thcph_vs_phth_CS = new TH1D( "PDF_s2thcph_vs_phth_CS",      "", nbin_ph, -180., 180. );


  TH1D*  PDF_1_vs_cth_HX   = new TH1D( "PDF_1_vs_cth_HX",       "", nbin_cth, -1., 1. );
  TH1D*  PDF_1_vs_ph_HX    = new TH1D( "PDF_1_vs_ph_HX",        "", nbin_ph, -180., 180. );
  TH1D*  PDF_1_vs_phth_HX  = new TH1D( "PDF_1_vs_phth_HX",      "", nbin_ph, -180., 180. );

  TH1D*  PDF_cth2_vs_cth_HX  = new TH1D( "PDF_cth2_vs_cth_HX",       "", nbin_cth, -1., 1. );
  TH1D*  PDF_cth2_vs_ph_HX   = new TH1D( "PDF_cth2_vs_ph_HX",        "", nbin_ph, -180., 180. );
  TH1D*  PDF_cth2_vs_phth_HX = new TH1D( "PDF_cth2_vs_phth_HX",      "", nbin_ph, -180., 180. );

  TH1D*  PDF_sth2c2ph_vs_cth_HX  = new TH1D( "PDF_sth2c2ph_vs_cth_HX",       "", nbin_cth, -1., 1. );
  TH1D*  PDF_sth2c2ph_vs_ph_HX   = new TH1D( "PDF_sth2c2ph_vs_ph_HX",        "", nbin_ph, -180., 180. );
  TH1D*  PDF_sth2c2ph_vs_phth_HX = new TH1D( "PDF_sth2c2ph_vs_phth_HX",      "", nbin_ph, -180., 180. );

  TH1D*  PDF_s2thcph_vs_cth_HX  = new TH1D( "PDF_s2thcph_vs_cth_HX",       "", nbin_cth, -1., 1. );
  TH1D*  PDF_s2thcph_vs_ph_HX   = new TH1D( "PDF_s2thcph_vs_ph_HX",        "", nbin_ph, -180., 180. );
  TH1D*  PDF_s2thcph_vs_phth_HX = new TH1D( "PDF_s2thcph_vs_phth_HX",      "", nbin_ph, -180., 180. );


  TH1D*  PDF_1_vs_cth_PX   = new TH1D( "PDF_1_vs_cth_PX",       "", nbin_cth, -1., 1. );
  TH1D*  PDF_1_vs_ph_PX    = new TH1D( "PDF_1_vs_ph_PX",        "", nbin_ph, -180., 180. );
  TH1D*  PDF_1_vs_phth_PX  = new TH1D( "PDF_1_vs_phth_PX",      "", nbin_ph, -180., 180. );

  TH1D*  PDF_cth2_vs_cth_PX  = new TH1D( "PDF_cth2_vs_cth_PX",       "", nbin_cth, -1., 1. );
  TH1D*  PDF_cth2_vs_ph_PX   = new TH1D( "PDF_cth2_vs_ph_PX",        "", nbin_ph, -180., 180. );
  TH1D*  PDF_cth2_vs_phth_PX = new TH1D( "PDF_cth2_vs_phth_PX",      "", nbin_ph, -180., 180. );

  TH1D*  PDF_sth2c2ph_vs_cth_PX  = new TH1D( "PDF_sth2c2ph_vs_cth_PX",       "", nbin_cth, -1., 1. );
  TH1D*  PDF_sth2c2ph_vs_ph_PX   = new TH1D( "PDF_sth2c2ph_vs_ph_PX",        "", nbin_ph, -180., 180. );
  TH1D*  PDF_sth2c2ph_vs_phth_PX = new TH1D( "PDF_sth2c2ph_vs_phth_PX",      "", nbin_ph, -180., 180. );

  TH1D*  PDF_s2thcph_vs_cth_PX  = new TH1D( "PDF_s2thcph_vs_cth_PX",       "", nbin_cth, -1., 1. );
  TH1D*  PDF_s2thcph_vs_ph_PX   = new TH1D( "PDF_s2thcph_vs_ph_PX",        "", nbin_ph, -180., 180. );
  TH1D*  PDF_s2thcph_vs_phth_PX = new TH1D( "PDF_s2thcph_vs_phth_PX",      "", nbin_ph, -180., 180. );


  TH1D*  PDF_1_vs_calpha        = new TH1D( "PDF_1_vs_calpha",             "", nbin_cth, -1., 1. );
  TH1D*  PDF_cth2_vs_calpha     = new TH1D( "PDF_cth2_vs_calpha",          "", nbin_cth, -1., 1. );
  TH1D*  PDF_sth2c2ph_vs_calpha = new TH1D( "PDF_sth2c2ph_vs_calpha",      "", nbin_cth, -1., 1. );
  TH1D*  PDF_s2thcph_vs_calpha  = new TH1D( "PDF_s2thcph_vs_calpha",       "", nbin_cth, -1., 1. );


  ///////////////// cycle of MC events ////////////////////////
  for(int i_MCevent = 1; i_MCevent <= n_MCevents_NORM; i_MCevent++){

	    if (i_MCevent%n_step == 0) {cout << n_step_*20 <<" % "<<endl; n_step_++;}

    // generation of dilepton in the pp CM

    // pT, rapidity, mass:

    double pT;
    double rap;
    double mass;

    pTrapMass_gen->GetRandom3( pT, rap, mass );

    // pL:

    double rap_sign = gRandom->Uniform(-1., 1.); rap_sign /= TMath::Abs(rap_sign);
    rap *= rap_sign;
    double mT = sqrt( mass*mass + pT*pT );
    double pL1 = 0.5 *mT * exp(rap);
    double pL2 = - 0.5 *mT * exp(-rap);
    double pL = pL1 + pL2;

    // Phi:

    double Phi = 2. * gPI_ * gRandom->Uniform(1.);

    // 4-vector:

    TLorentzVector dilepton;
    dilepton.SetXYZM( pT * cos(Phi) , pT * sin(Phi), pL, mass );


    // random extraction of decay angles (generic reference frame),
    // with isotropic distribution!

    double costh_gen = -1. + 2. * gRandom->Uniform(1.);
    double sinth_gen = sqrt( 1. - costh_gen*costh_gen );
    double phi_gen = 2. * gPI_ * gRandom->Uniform(1.);

    // lepton momentum in the dilepton rest frame:

    double p_lepton_DILEP = sqrt( 0.25*mass*mass - Mlepton_*Mlepton_ );

    TLorentzVector lepton_DILEP;

    lepton_DILEP.SetXYZM( p_lepton_DILEP * sinth_gen * cos(phi_gen),
                          p_lepton_DILEP * sinth_gen * sin(phi_gen),
                          p_lepton_DILEP * costh_gen,
                          Mlepton_ );


    // reference directions to calculate angles:

    TVector3 lab_to_dilep = -dilepton.BoostVector();

    TLorentzVector beam1_DILEP = beam1_LAB_;
    beam1_DILEP.Boost(lab_to_dilep);         // beam1 in the dilepton rest frame
    TLorentzVector beam2_DILEP = beam2_LAB_;
    beam2_DILEP.Boost(lab_to_dilep);         // beam2 in the dilepton rest frame

    TVector3 beam1_direction     = beam1_DILEP.Vect().Unit();
    TVector3 beam2_direction     = beam2_DILEP.Vect().Unit();
    TVector3 dilep_direction     = dilepton.Vect().Unit();
    TVector3 beam1_beam2_bisect  = ( beam1_direction - beam2_direction ).Unit();


    // all polarization frames have the same Y axis = the normal to the plane formed by
    // the directions of the colliding hadrons

    TVector3 Yaxis = ( beam1_direction.Cross( beam2_direction ) ).Unit();

    // flip of y axis with rapidity

    if ( rap < 0 ) Yaxis = - Yaxis;

    TVector3 perpendicular_to_beam = ( beam1_beam2_bisect.Cross( Yaxis ) ).Unit();


    // step 1: transform (rotation) lepton momentum components from generation frame
    // to the frame with x,y,z axes as in the laboratory

    TVector3 oldZaxis = beam1_beam2_bisect;

    TVector3 oldYaxis = Yaxis;
    TVector3 oldXaxis = oldYaxis.Cross(oldZaxis);

    TRotation rotation;
    rotation.RotateAxes(oldXaxis, oldYaxis, oldZaxis);
                     // transforms coordinates from the "old" frame to the "xyz" frame

    TLorentzVector lepton_DILEP_xyz = lepton_DILEP;

    lepton_DILEP_xyz.Transform(rotation);
                     // lepton_DILEP_xyz is the lepton in the dilepton rest frame
                     // wrt to the lab axes

    // CS frame

    TVector3 newZaxis = beam1_beam2_bisect;
    TVector3 newYaxis = Yaxis;
    TVector3 newXaxis = newYaxis.Cross( newZaxis );

    rotation.SetToIdentity();
    rotation.RotateAxes( newXaxis, newYaxis, newZaxis );
    rotation.Invert();  // transforms coordinates from the "xyz" frame to the new frame

    TVector3 lepton_DILEP_rotated = lepton_DILEP_xyz.Vect();

    lepton_DILEP_rotated.Transform(rotation);

    double costh_CS = lepton_DILEP_rotated.CosTheta();
    double phi_CS = lepton_DILEP_rotated.Phi() * 180. / gPI_;
    double phith_CS;
    if ( costh_CS < 0. ) phith_CS = phi_CS - 135.;
    if ( costh_CS > 0. ) phith_CS = phi_CS - 45.;
    if ( phith_CS < -180. ) phith_CS = 360. + phith_CS;


    // HELICITY frame

    newZaxis = dilep_direction;
    newYaxis = Yaxis;
    newXaxis = newYaxis.Cross( newZaxis );

    rotation.SetToIdentity();
    rotation.RotateAxes( newXaxis, newYaxis, newZaxis );
    rotation.Invert();

    lepton_DILEP_rotated = lepton_DILEP_xyz.Vect();

    lepton_DILEP_rotated.Transform(rotation);

    double costh_HX = lepton_DILEP_rotated.CosTheta();
    double phi_HX = lepton_DILEP_rotated.Phi() * 180. / gPI_;
    double phith_HX;
    if ( costh_HX < 0. ) phith_HX = phi_HX - 135.;
    if ( costh_HX > 0. ) phith_HX = phi_HX - 45.;
    if ( phith_HX < -180. ) phith_HX = 360. + phith_HX;

		//PhiHX test
		//if(PhiHX_test) { if(phi_HX>80. || phi_HX<91.) continue; }

    // PERPENDICULAR HELICITY frame

    newZaxis = perpendicular_to_beam;
    newYaxis = Yaxis;
    newXaxis = newYaxis.Cross( newZaxis );

    rotation.SetToIdentity();
    rotation.RotateAxes( newXaxis, newYaxis, newZaxis );
    rotation.Invert();

    lepton_DILEP_rotated = lepton_DILEP_xyz.Vect();

    lepton_DILEP_rotated.Transform(rotation);

    double costh_PX = lepton_DILEP_rotated.CosTheta();
    double phi_PX = lepton_DILEP_rotated.Phi() * 180. / gPI_;
    double phith_PX;
    if ( costh_PX < 0. ) phith_PX = phi_PX - 135.;
    if ( costh_PX > 0. ) phith_PX = phi_PX - 45.;
    if ( phith_PX < -180. ) phith_PX = 360. + phith_PX;


    // invariant polarization angle

    double cosalpha = sqrt( 1. - pow(costh_PX, 2.) ) * sin( lepton_DILEP_rotated.Phi() );



    // lepton 4-vectors in the LAB frame:

    TVector3 dilep_to_lab = dilepton.BoostVector();

    *lepP = lepton_DILEP_xyz;
    lepP->Boost(dilep_to_lab);
    lepN->SetPxPyPzE(-lepton_DILEP_xyz.Px(),-lepton_DILEP_xyz.Py(),-lepton_DILEP_xyz.Pz(),lepton_DILEP_xyz.E());
    lepN->Boost(dilep_to_lab);

    double lepP_pT  = lepP->Pt();
    double lepN_pT  = lepN->Pt();

    double lepP_eta = lepP->PseudoRapidity();
    double lepN_eta = lepN->PseudoRapidity();

    //bool isEventAccepted = isMuonInAcceptance( FidCuts-1, lepP_pT, lepP_eta )
    //                     * isMuonInAcceptance( FidCuts-1, lepN_pT, lepN_eta );

		double deltaPhi = lepN->Phi() - lepP->Phi();
		if(deltaPhi > TMath::Pi()) deltaPhi -= 2.*TMath::Pi();
		else if(deltaPhi < -TMath::Pi()) deltaPhi += 2.*TMath::Pi();
		double deltaEta = lepN_eta - lepP_eta; 
		double deltaPT  = lepN_pT  - lepP_pT ;
		double deltaREll = TMath::Sqrt(deltaEta*deltaEta + TMath::Power(1.2*deltaPhi,2));
		double deltaREllDpt = TMath::Sqrt(deltaREll*deltaREll + TMath::Power(0.00157*deltaPT,2));

		bool accept = isMuonInAcceptance( FidCuts-1, lepP_pT, lepP_eta )
		 	* isMuonInAcceptance( FidCuts-1, lepN_pT, lepN_eta );
		bool isEventAccepted = false;
		if(accept && (!cutDeltaREllDpt || (cutDeltaREllDpt && deltaREllDpt > DeltaREllDptValue) )) isEventAccepted = true;

		//if(deltaREllDpt<DeltaREllDptValue)
		// 	cout<<"4: deltaREllDpt: "<<deltaREllDpt<<endl;
		//if(accept && !isEventAccepted){
		//	cout<<"accept: "<<accept<<endl;
		//	cout<<"isEventAccepted: "<<isEventAccepted<<endl;
		//}

    double epsilon = singleLeptonEfficiency( lepP_pT, lepP_eta, nEff, fInEff, hEvalEff, MCeff, TEff)
                   * singleLeptonEfficiency( lepN_pT, lepN_eta, nEff, fInEff, hEvalEff, MCeff, TEff);

    double costh_DILEff;
    double phi_DILEff;
    if(nDileptonEff>200 && nDileptonEff<211) {costh_DILEff=costh_CS; phi_DILEff=phi_CS;}
    if(nDileptonEff>210 && nDileptonEff<221) {costh_DILEff=costh_HX; phi_DILEff=phi_HX;}
    if(nDileptonEff>220 && nDileptonEff<231) {costh_DILEff=costh_PX; phi_DILEff=phi_PX;}
    double DileptonEff = DiLeptonEfficiency( costh_DILEff, phi_DILEff, nDileptonEff, fInDileptonEff , MCDileptoneff);


    double costh_RhoFactor;
    double phi_RhoFactor;
    if(nRhoFactor>300 && nRhoFactor<311) {costh_RhoFactor=costh_CS; phi_RhoFactor=phi_CS;}
    if(nRhoFactor>310 && nRhoFactor<321) {costh_RhoFactor=costh_HX; phi_RhoFactor=phi_HX;}
    if(nRhoFactor>320 && nRhoFactor<331) {costh_RhoFactor=costh_PX; phi_RhoFactor=phi_PX;}
	double RhoFactor = EvaluateRhoFactor( costh_RhoFactor, phi_RhoFactor, nRhoFactor, fInRhoFactor, rap, pT, StatVarRho);


	if(!useAmapApproach) epsilon = epsilon*DileptonEff*RhoFactor;

	if(useAmapApproach){
	    double costh_Amap=0;
	    double phi_Amap=0;
	    if(nAmap>10000 && nAmap<20000) {costh_Amap=costh_CS; phi_Amap=phi_CS;}
	    if(nAmap>20000 && nAmap<30000) {costh_Amap=costh_HX; phi_Amap=phi_HX;}
	    if(nAmap>30000 && nAmap<40000) {costh_Amap=costh_PX; phi_Amap=phi_PX;}
		double AmapValue=EvaluateAmap( costh_Amap, phi_Amap, nAmap, fInAmap, rap, pT);
		double AmapDenominator= DenominatorAmapEfficiency( lepP_pT, lepP_eta, nDenominatorAmap, fInEff_nDenominatorAmap, hEvalEff_nDenominatorAmap, true, TEff_nDenominatorAmap) * DenominatorAmapEfficiency( lepN_pT, lepN_eta, nDenominatorAmap, fInEff_nDenominatorAmap, hEvalEff_nDenominatorAmap, true, TEff_nDenominatorAmap);
		epsilon = epsilon*AmapValue/AmapDenominator;
	}


	if(epsilon>1&&ForceEpsSmallerOne) {epsilon=1;}

    if ( epsilon > min_dileptonEff && isEventAccepted ) {

        double temp;

        //  update the sums and fill histograms for the PDF calculations

        ++n_acceptedEvents;

        temp = epsilon;
        sum_1           += temp;

        PDF_1_vs_cth_CS->Fill( costh_CS, temp );
        PDF_1_vs_ph_CS->Fill( phi_CS, temp );
        PDF_1_vs_phth_CS->Fill( phith_CS, temp );

        PDF_1_vs_cth_HX->Fill( costh_HX, temp );
        PDF_1_vs_ph_HX->Fill( phi_HX, temp );
        PDF_1_vs_phth_HX->Fill( phith_HX, temp );

        PDF_1_vs_cth_PX->Fill( costh_PX, temp );
        PDF_1_vs_ph_PX->Fill( phi_PX, temp );
        PDF_1_vs_phth_PX->Fill( phith_PX, temp );

        PDF_1_vs_calpha->Fill( cosalpha, temp );

        temp = epsilon * costh_CS*costh_CS;
        sum_cth2_CS     += temp;
        PDF_cth2_vs_cth_CS->Fill( costh_CS, temp );
        PDF_cth2_vs_ph_CS->Fill( phi_CS, temp );
        PDF_cth2_vs_phth_CS->Fill( phith_CS, temp );

        temp = epsilon * (1.- costh_CS*costh_CS) * cos(gPI_/90.*phi_CS);
        sum_sth2c2ph_CS += temp;
        PDF_sth2c2ph_vs_cth_CS->Fill( costh_CS, temp );
        PDF_sth2c2ph_vs_ph_CS->Fill( phi_CS, temp );
        PDF_sth2c2ph_vs_phth_CS->Fill( phith_CS, temp );

        temp = epsilon * 2.*costh_CS*sqrt(1.-costh_CS*costh_CS) * cos(gPI_/180.*phi_CS);
        sum_s2thcph_CS  += temp;
        PDF_s2thcph_vs_cth_CS->Fill( costh_CS, temp );
        PDF_s2thcph_vs_ph_CS->Fill( phi_CS, temp );
        PDF_s2thcph_vs_phth_CS->Fill( phith_CS, temp );


        temp = epsilon * costh_HX*costh_HX;
        sum_cth2_HX     += temp;
        PDF_cth2_vs_cth_HX->Fill( costh_HX, temp );
        PDF_cth2_vs_ph_HX->Fill( phi_HX, temp );
        PDF_cth2_vs_phth_HX->Fill( phith_HX, temp );

        temp = epsilon * (1.- costh_HX*costh_HX) * cos(gPI_/90.*phi_HX);
        sum_sth2c2ph_HX += temp;
        PDF_sth2c2ph_vs_cth_HX->Fill( costh_HX, temp );
        PDF_sth2c2ph_vs_ph_HX->Fill( phi_HX, temp );
        PDF_sth2c2ph_vs_phth_HX->Fill( phith_HX, temp );

        temp = epsilon * 2.*costh_HX*sqrt(1.-costh_HX*costh_HX) * cos(gPI_/180.*phi_HX);
        sum_s2thcph_HX  += temp;
        PDF_s2thcph_vs_cth_HX->Fill( costh_HX, temp );
        PDF_s2thcph_vs_ph_HX->Fill( phi_HX, temp );
        PDF_s2thcph_vs_phth_HX->Fill( phith_HX, temp );


        temp = epsilon * costh_PX*costh_PX;
        sum_cth2_PX     += temp;
        PDF_cth2_vs_cth_PX->Fill( costh_PX, temp );
        PDF_cth2_vs_ph_PX->Fill( phi_PX, temp );
        PDF_cth2_vs_phth_PX->Fill( phith_PX, temp );
        PDF_cth2_vs_calpha->Fill( cosalpha, temp );

        temp = epsilon * (1.- costh_PX*costh_PX) * cos(gPI_/90.*phi_PX);
        sum_sth2c2ph_PX += temp;
        PDF_sth2c2ph_vs_cth_PX->Fill( costh_PX, temp );
        PDF_sth2c2ph_vs_ph_PX->Fill( phi_PX, temp );
        PDF_sth2c2ph_vs_phth_PX->Fill( phith_PX, temp );
        PDF_sth2c2ph_vs_calpha->Fill( cosalpha, temp );

        temp = epsilon * 2.*costh_PX*sqrt(1.-costh_PX*costh_PX) * cos(gPI_/180.*phi_PX);
        sum_s2thcph_PX  += temp;
        PDF_s2thcph_vs_cth_PX->Fill( costh_PX, temp );
        PDF_s2thcph_vs_ph_PX->Fill( phi_PX, temp );
        PDF_s2thcph_vs_phth_PX->Fill( phith_PX, temp );
        PDF_s2thcph_vs_calpha->Fill( cosalpha, temp );

    }


  } // end loop of MC integration

  cout << endl;



  sum_1           *= n_signalEvents / n_acceptedEvents;

  sum_cth2_CS     *= n_signalEvents / n_acceptedEvents;
  sum_sth2c2ph_CS *= n_signalEvents / n_acceptedEvents;
  sum_s2thcph_CS  *= n_signalEvents / n_acceptedEvents;

  sum_cth2_HX     *= n_signalEvents / n_acceptedEvents;
  sum_sth2c2ph_HX *= n_signalEvents / n_acceptedEvents;
  sum_s2thcph_HX  *= n_signalEvents / n_acceptedEvents;

  sum_cth2_PX     *= n_signalEvents / n_acceptedEvents;
  sum_sth2c2ph_PX *= n_signalEvents / n_acceptedEvents;
  sum_s2thcph_PX  *= n_signalEvents / n_acceptedEvents;


  // input data file is no more used and is closed

  dataFile->Close();


  // output ntuples for the CS, HX, PX (perpendicular helicity) frames
  // lth = lambdatheta, lph = lambdaphi, ltp = lambdathetaphi

  TTree* results_CS = new TTree("lambdaCS","lambdaCS");
  double lth_CS;        results_CS->Branch("lth",         &lth_CS,         "lth/D");
  double lph_CS;        results_CS->Branch("lph",         &lph_CS,         "lph/D");
  double ltp_CS;        results_CS->Branch("ltp",         &ltp_CS,         "ltp/D");
  double lthstar_CS;    results_CS->Branch("lthstar",     &lthstar_CS,     "lthstar/D");
  double lphstar_CS;    results_CS->Branch("lphstar",     &lphstar_CS,     "lphstar/D");
  double ltilde_CS;     results_CS->Branch("ltilde",      &ltilde_CS,      "ltilde/D");
  int    positivity_CS; results_CS->Branch("positivity",  &positivity_CS,  "positivity/I");

  TTree* results_HX = new TTree("lambdaHX","lambdaHX");
  double lth_HX;        results_HX->Branch("lth",         &lth_HX,         "lth/D");
  double lph_HX;        results_HX->Branch("lph",         &lph_HX,         "lph/D");
  double ltp_HX;        results_HX->Branch("ltp",         &ltp_HX,         "ltp/D");
  double lthstar_HX;    results_HX->Branch("lthstar",     &lthstar_HX,     "lthstar/D");
  double lphstar_HX;    results_HX->Branch("lphstar",     &lphstar_HX,     "lphstar/D");
  double ltilde_HX;     results_HX->Branch("ltilde",      &ltilde_HX,      "ltilde/D");
  int    positivity_HX; results_HX->Branch("positivity",  &positivity_HX,  "positivity/I");

  TTree* results_PX = new TTree("lambdaPX","lambdaPX");
  double lth_PX;        results_PX->Branch("lth",         &lth_PX,         "lth/D");
  double lph_PX;        results_PX->Branch("lph",         &lph_PX,         "lph/D");
  double ltp_PX;        results_PX->Branch("ltp",         &ltp_PX,         "ltp/D");
  double lthstar_PX;    results_PX->Branch("lthstar",     &lthstar_PX,     "lthstar/D");
  double lphstar_PX;    results_PX->Branch("lphstar",     &lphstar_PX,     "lphstar/D");
  double ltilde_PX;     results_PX->Branch("ltilde",      &ltilde_PX,      "ltilde/D");
  int    positivity_PX; results_PX->Branch("positivity",  &positivity_PX,  "positivity/I");


  ////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////
  ////// sampling of the parameter space ("fit") with the Markov chain

  // starting point for the "chain" of sampled values

/*  lth_CS = -0.5;
  lph_CS = 0.5;
  ltp_CS = 0.0;
  lth_HX = 0.5;
  lph_HX = 0.0;
  ltp_HX = 0.0;
  lth_PX = 0.5;
  lph_PX = 0.0;
  ltp_PX = 0.0;
*/
  lth_CS = 0.0;
  lph_CS = 0.0;
  ltp_CS = 0.0;

  lth_HX = 0.0;
  lph_HX = 0.0;
  ltp_HX = 0.0;

  lth_PX = 0.0;
  lph_PX = 0.0;
  ltp_PX = 0.0;

  if(RandomLambdaTheta){
	  lth_CS = gRandom->Uniform(1.)-0.5;
	  lph_CS = 0.0;
	  ltp_CS = 0.0;

	  lth_HX = gRandom->Uniform(1.)-0.5;
	  lph_HX = 0.0;
	  ltp_HX = 0.0;

	  lth_PX = gRandom->Uniform(1.)-0.5;
	  lph_PX = 0.0;
	  ltp_PX = 0.0;
  }

  cout<<"starting values:"<<endl;
  cout<<"lth_CS = "<<lth_CS<<endl;
  cout<<"lph_CS = "<<lph_CS<<endl;
  cout<<"ltp_CS = "<<ltp_CS<<endl;
  cout<<"lth_HX = "<<lth_HX<<endl;
  cout<<"lph_HX = "<<lph_HX<<endl;
  cout<<"ltp_HX = "<<ltp_HX<<endl;
  cout<<"lth_PX = "<<lth_PX<<endl;
  cout<<"lph_PX = "<<lph_PX<<endl;
  cout<<"ltp_PX = "<<ltp_PX<<endl;



  double loglikelihood_CS = -1.e30;  // intial (arbitrary) values
  double loglikelihood_HX = -1.e30;
  double loglikelihood_PX = -1.e30;

  // test histograms of the distributions after the burn-in period,
  // for a subsequent adjustment of the algorithm

  TH1D* test_lth_CS = new TH1D( "test_lth_CS", "", 3000, -1.5, 1.5 );
  TH1D* test_lph_CS = new TH1D( "test_lph_CS", "", 3000, -1.5, 1.5 );
  TH1D* test_ltp_CS = new TH1D( "test_ltp_CS", "", 3000, -1.0, 1.0 );

  TH1D* test_lth_HX = new TH1D( "test_lth_HX", "", 3000, -1.5, 1.5 );
  TH1D* test_lph_HX = new TH1D( "test_lph_HX", "", 3000, -1.5, 1.5 );
  TH1D* test_ltp_HX = new TH1D( "test_ltp_HX", "", 3000, -1.0, 1.0 );

  TH1D* test_lth_PX = new TH1D( "test_lth_PX", "", 3000, -1.5, 1.5 );
  TH1D* test_lph_PX = new TH1D( "test_lph_PX", "", 3000, -1.5, 1.5 );
  TH1D* test_ltp_PX = new TH1D( "test_ltp_PX", "", 3000, -1.0, 1.0 );


  cout << endl;
  cout << "Sampling of the parameter space (" << n_burnIn << "+"
                         << n_sampledPoints - n_burnIn << " iterations)" << endl;
  cout << "------------------------------------------------------------" << endl;
  cout << "Progress: "<< endl;

  // Initial values of the sigmas of the proposal functions.
  // Speed and efficiency of the algorithm are extremely sensistive to the
  // sigma of the proposal function. The intial sigma should be large enough, because
  // a large sigma helps finding qickly where the distribution is centred. However, if the
  // sigma is too large, the efficiency of the method is very small, i.e. many events are
  // rejected and it is necessary to increase a lot the number of extractions.
  // For these reasons, during the burn-in period a large sigma is used for each parameter
  // and after this period the sigma is adjusted looking at the distribution produced so far.
  // Optimal: sigma(proposal) ~0.3 * sigma(target)


  double proposalWidth_lth_CS = proposalWidth_lth_CS_start;
  double proposalWidth_lph_CS = proposalWidth_lph_CS_start;
  double proposalWidth_ltp_CS = proposalWidth_ltp_CS_start;

  double proposalWidth_lth_HX = proposalWidth_lth_HX_start;
  double proposalWidth_lph_HX = proposalWidth_lph_HX_start;
  double proposalWidth_ltp_HX = proposalWidth_ltp_HX_start;

  double proposalWidth_lth_PX = proposalWidth_lth_PX_start;
  double proposalWidth_lph_PX = proposalWidth_lph_PX_start;
  double proposalWidth_ltp_PX = proposalWidth_ltp_PX_start;

  int rejectPointCS=0;
  int rejectPointHX=0;
  int rejectPointPX=0;

  n_step = n_sampledPoints/5;  // visualize progress of the parameter sampling
  n_step_=1;

  int actuallySampledPoints=0;
 int BuffStep=0;

  for(int i_sampledPoint = 1; i_sampledPoint <= n_sampledPoints; i_sampledPoint++){

	    if (i_sampledPoint%n_step == 0&&i_sampledPoint>BuffStep) {cout << n_step_*20 <<" % "<<endl; n_step_++; BuffStep=i_sampledPoint; }

	    // initialize positivity condition to positive

     positivity_CS = 1;
     positivity_HX = 1;
     positivity_PX = 1;

     // random extraction of parameter values according to the Metropolis-Hastings
     // implementation of the Markov chain

     double lth_CS_candidate, lph_CS_candidate, ltp_CS_candidate;
     double lth_HX_candidate, lph_HX_candidate, ltp_HX_candidate;
     double lth_PX_candidate, lph_PX_candidate, ltp_PX_candidate;

     // sampling of the parameter space using Gaussian "proposal PDF"
     // centred in the "previous" parameter values and generating
     // "candidate" values for the next iteration

     extractFromProposalPDF( lth_CS_candidate,     lph_CS_candidate,     ltp_CS_candidate,
                             lth_CS,               lph_CS,               ltp_CS,
                             proposalWidth_lth_CS, proposalWidth_lph_CS, proposalWidth_ltp_CS );
     extractFromProposalPDF( lth_HX_candidate,     lph_HX_candidate,     ltp_HX_candidate,
                             lth_HX,               lph_HX,               ltp_HX,
                             proposalWidth_lth_HX, proposalWidth_lph_HX, proposalWidth_ltp_HX );
     extractFromProposalPDF( lth_PX_candidate,     lph_PX_candidate,     ltp_PX_candidate,
                             lth_PX,               lph_PX,               ltp_PX,
                             proposalWidth_lth_PX, proposalWidth_lph_PX, proposalWidth_ltp_PX );


     // calculation of the likelihoods to decide to keep or reject the candidates

     double loglikelihood_CS_candidate = 0; // calculate log(likelihood) to avoid that
     double loglikelihood_HX_candidate = 0; // the likelihood goes to zero in the product over all events
     double loglikelihood_PX_candidate = 0; // because of limited numerical precision

     const double min_likelihood = 1.e-100;


     // loop over arrays of selected dilepton events to calculate the likelihoods
     /////////////////////////////////////////////////////////////////////////////////////
     for ( int i = 0; i < n_signalEvents; i++ ) {

       double singleEventLikelihood_CS = (  w_1[i] + lth_CS_candidate * w_cth2_CS[i]
                                                   + lph_CS_candidate * w_sth2c2ph_CS[i]
                                                   + ltp_CS_candidate * w_s2thcph_CS[i]  )
                                      /  (  sum_1  + lth_CS_candidate * sum_cth2_CS
                                                   + lph_CS_candidate * sum_sth2c2ph_CS
                                                   + ltp_CS_candidate * sum_s2thcph_CS   ) ;
       if ( singleEventLikelihood_CS < min_likelihood ) singleEventLikelihood_CS = min_likelihood;
       loglikelihood_CS_candidate += log( singleEventLikelihood_CS );

       double singleEventLikelihood_HX = (  w_1[i] + lth_HX_candidate * w_cth2_HX[i]
                                                   + lph_HX_candidate * w_sth2c2ph_HX[i]
                                                   + ltp_HX_candidate * w_s2thcph_HX[i]  )
                                      /  (  sum_1  + lth_HX_candidate * sum_cth2_HX
                                                   + lph_HX_candidate * sum_sth2c2ph_HX
                                                   + ltp_HX_candidate * sum_s2thcph_HX   ) ;
       if ( singleEventLikelihood_HX < min_likelihood ) singleEventLikelihood_HX = min_likelihood;
       loglikelihood_HX_candidate += log( singleEventLikelihood_HX );

       double singleEventLikelihood_PX = (  w_1[i] + lth_PX_candidate * w_cth2_PX[i]
                                                   + lph_PX_candidate * w_sth2c2ph_PX[i]
                                                   + ltp_PX_candidate * w_s2thcph_PX[i]  )
                                      /  (  sum_1  + lth_PX_candidate * sum_cth2_PX
                                                   + lph_PX_candidate * sum_sth2c2ph_PX
                                                   + ltp_PX_candidate * sum_s2thcph_PX   ) ;
       if ( singleEventLikelihood_PX < min_likelihood ) singleEventLikelihood_PX = min_likelihood;
       loglikelihood_PX_candidate += log( singleEventLikelihood_PX );

     } // end of loop over arrays of selected dilepton events


     // apply Metropolis-Hastings algorithm: the candidate parameter values are
     // kept or rejected depending on the likelihood ratio wrt the "previous" values.
     // The likelihood ratio, if smaller than 1, represents the "probability" with which
     // the new values must be taken as next "good" values of the "chain" and as starting
     // point for the next iteration.
     // If the ratio is > 1 the new values are always taken as good.
     // The condition "likelihood ratio > or < 1" is translated into the condition
     // "log(likelihood) difference > or < 0"

     double loglikelihood_CS_difference = loglikelihood_CS_candidate - loglikelihood_CS;
     if(  loglikelihood_CS_difference > 0.  ||  log( gRandom->Uniform(1.) ) < loglikelihood_CS_difference  ) {
         lth_CS = lth_CS_candidate; lph_CS = lph_CS_candidate; ltp_CS = ltp_CS_candidate;
         loglikelihood_CS = loglikelihood_CS_candidate;
         if ( i_sampledPoint > n_burnIn ) {  // reject first n_burnIn extractions in the chain

            if ( TMath::Abs( lph_CS ) > 0.5*( 1 + lth_CS ) || lth_CS*lth_CS + 2.*ltp_CS*ltp_CS > 1
            || TMath::Abs( ltp_CS ) > 0.5*( 1 - lph_CS )
            || (  (1.+2.*lph_CS)*(1.+2.*lph_CS) + 2.*ltp_CS*ltp_CS > 1 && lph_CS < -1./3. ) )
            positivity_CS = 0; // apply positivity constraint in a flag-variable of the ntuple

            calcLambdastar( lthstar_CS, lphstar_CS, lth_CS, lph_CS, ltp_CS );
            ltilde_CS = (lth_CS + 3.*lph_CS)/(1.-lph_CS);
            results_CS->Fill();
         }
         else if ( i_sampledPoint > n_burnIn/2 ) { test_lth_CS->Fill( lth_CS );
                                                   test_lph_CS->Fill( lph_CS );
                                                   test_ltp_CS->Fill( ltp_CS ); }
     }
     else {if(i_sampledPoint > n_burnIn) rejectPointCS++; }
     if ( i_sampledPoint == n_burnIn +1 ) { proposalWidth_lth_CS = TMath::Max( proposalWidthTimesRMS * test_lth_CS->GetRMS(), 0.001 );
                                            proposalWidth_lph_CS = TMath::Max( proposalWidthTimesRMS * test_lph_CS->GetRMS(), 0.001 );
                                            proposalWidth_ltp_CS = TMath::Max( proposalWidthTimesRMS * test_ltp_CS->GetRMS(), 0.001 ); }
     // fill test histograms in the second half of the burn-in period (when the algorithm
     // should already have found where the bulk of the distribution is) in order to estimate
     // the sigmas of the output distributions and to adjust the sigmas of the proposal functions
     // accordingly (optimal: between 1/4 and 1/3 of the target distributions)


     double loglikelihood_HX_difference = loglikelihood_HX_candidate - loglikelihood_HX;
     if(  loglikelihood_HX_difference > 0.  ||  log( gRandom->Uniform(1.) ) < loglikelihood_HX_difference  ) {
         lth_HX = lth_HX_candidate; lph_HX = lph_HX_candidate; ltp_HX = ltp_HX_candidate;
         loglikelihood_HX = loglikelihood_HX_candidate;
         if ( i_sampledPoint > n_burnIn ) {

            if ( TMath::Abs( lph_HX ) > 0.5*( 1 + lth_HX ) || lth_HX*lth_HX + 2.*ltp_HX*ltp_HX > 1
            || TMath::Abs( ltp_HX ) > 0.5*( 1 - lph_HX )
            || (  (1.+2.*lph_HX)*(1.+2.*lph_HX) + 2.*ltp_HX*ltp_HX > 1 && lph_HX < -1./3. ) )
            positivity_HX = 0;

            calcLambdastar( lthstar_HX, lphstar_HX, lth_HX, lph_HX, ltp_HX );
            ltilde_HX = (lth_HX + 3.*lph_HX)/(1.-lph_HX);
            results_HX->Fill();
         }
         else if ( i_sampledPoint > n_burnIn/2 ) { test_lth_HX->Fill( lth_HX );
                                                   test_lph_HX->Fill( lph_HX );
                                                   test_ltp_HX->Fill( ltp_HX ); }
     }
     else {if(i_sampledPoint > n_burnIn) rejectPointHX++; }
     if ( i_sampledPoint == n_burnIn +1 ) { proposalWidth_lth_HX = TMath::Max( proposalWidthTimesRMS * test_lth_HX->GetRMS(), 0.001 );
                                            proposalWidth_lph_HX = TMath::Max( proposalWidthTimesRMS * test_lph_HX->GetRMS(), 0.001 );
                                            proposalWidth_ltp_HX = TMath::Max( proposalWidthTimesRMS * test_ltp_HX->GetRMS(), 0.001 ); }


     double loglikelihood_PX_difference = loglikelihood_PX_candidate - loglikelihood_PX;
     if(  loglikelihood_PX_difference > 0.  ||  log( gRandom->Uniform(1.) ) < loglikelihood_PX_difference  ) {
         lth_PX = lth_PX_candidate; lph_PX = lph_PX_candidate; ltp_PX = ltp_PX_candidate;
         loglikelihood_PX = loglikelihood_PX_candidate;
         if ( i_sampledPoint > n_burnIn ) {

            if ( TMath::Abs( lph_PX ) > 0.5*( 1 + lth_PX ) || lth_PX*lth_PX + 2.*ltp_PX*ltp_PX > 1
            || TMath::Abs( ltp_PX ) > 0.5*( 1 - lph_PX )
            || (  (1.+2.*lph_PX)*(1.+2.*lph_PX) + 2.*ltp_PX*ltp_PX > 1 && lph_PX < -1./3. ) )
            positivity_PX = 0;

            calcLambdastar( lthstar_PX, lphstar_PX, lth_PX, lph_PX, ltp_PX );
            ltilde_PX = (lth_PX + 3.*lph_PX)/(1.-lph_PX);
            results_PX->Fill();
         }
         else if ( i_sampledPoint > n_burnIn/2 ) { test_lth_PX->Fill( lth_PX );
                                                   test_lph_PX->Fill( lph_PX );
                                                   test_ltp_PX->Fill( ltp_PX ); }
     }
     else {if(i_sampledPoint > n_burnIn) {rejectPointPX++; i_sampledPoint=i_sampledPoint-1;}}
     if ( i_sampledPoint == n_burnIn +1 ) { proposalWidth_lth_PX = TMath::Max( proposalWidthTimesRMS * test_lth_PX->GetRMS(), 0.001 );
                                            proposalWidth_lph_PX = TMath::Max( proposalWidthTimesRMS * test_lph_PX->GetRMS(), 0.001 );
                                            proposalWidth_ltp_PX = TMath::Max( proposalWidthTimesRMS * test_ltp_PX->GetRMS(), 0.001 ); }



     actuallySampledPoints++;

  } // end of parameter sampling loop

  cout << endl << endl;

  cout<<"Number of steps of Marcov Chain (without burn-in): "<<actuallySampledPoints-n_burnIn<<endl;

  cout<<"Fraction of rejected Points CS: "<<double(rejectPointCS)/double(actuallySampledPoints-n_burnIn)<<endl;
  cout<<"Fraction of rejected Points HX: "<<double(rejectPointHX)/double(actuallySampledPoints-n_burnIn)<<endl;
  cout<<"Fraction of rejected Points PX: "<<double(rejectPointPX)/double(actuallySampledPoints-n_burnIn)<<endl;

  TH1D*  MetropolisHastingsAcceptanceCS   = new TH1D( "MetropolisHastingsAcceptanceCS",       "", 100, 0, 1 );
  TH1D*  MetropolisHastingsAcceptanceHX   = new TH1D( "MetropolisHastingsAcceptanceHX",       "", 100, 0, 1 );
  TH1D*  MetropolisHastingsAcceptancePX   = new TH1D( "MetropolisHastingsAcceptancePX",       "", 100, 0, 1 );

  MetropolisHastingsAcceptanceCS->Fill(1.-double(rejectPointCS)/double(actuallySampledPoints-n_burnIn));
  MetropolisHastingsAcceptanceHX->Fill(1.-double(rejectPointHX)/double(actuallySampledPoints-n_burnIn));
  MetropolisHastingsAcceptancePX->Fill(1.-double(rejectPointPX)/double(actuallySampledPoints-n_burnIn));

  ///// Extract numerical values of the lambda parameters from the result TTrees /////

    // extremes and binning of lambda plots
    const double l_min = -3;
    const double l_max =  3;
    const double l_step_1D = 0.02;

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

    int n_entries = int( results_CS->GetEntries() );
    for ( int i_entry = 0; i_entry < n_entries; i_entry++ ) {
      results_CS->GetEvent( i_entry );
      h_lth_CS->Fill( lth_CS );
      h_lph_CS->Fill( lph_CS );
      h_ltp_CS->Fill( ltp_CS );
      h_ltilde_CS->Fill( ltilde_CS );
      h_lphstar_CS->Fill( lphstar_CS );
      h_lthstar_CS->Fill( lthstar_CS );
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

    // end of loop over entries in the ntuples

  //  cout << endl << endl;


    TTree* Results = new TTree("Results","Results");
    Results->Branch("lthCS",         &lth_CS,         "lthCS/D");
    Results->Branch("lphCS",         &lph_CS,         "lphCS/D");
    Results->Branch("ltpCS",         &ltp_CS,         "ltpCS/D");
    Results->Branch("lthstarCS",     &lthstar_CS,     "lthstarCS/D");
    Results->Branch("lphstarCS",     &lphstar_CS,     "lphstarCS/D");
    Results->Branch("ltildeCS",      &ltilde_CS,      "ltildeCS/D");
    Results->Branch("positivityCS",  &positivity_CS,  "positivityCS/I");
    Results->Branch("lthHX",         &lth_HX,         "lthHX/D");
    Results->Branch("lphHX",         &lph_HX,         "lphHX/D");
    Results->Branch("ltpHX",         &ltp_HX,         "ltpHX/D");
    Results->Branch("lthstarHX",     &lthstar_HX,     "lthstarHX/D");
    Results->Branch("lphstarHX",     &lphstar_HX,     "lphstarHX/D");
    Results->Branch("ltildeHX",      &ltilde_HX,      "ltildeHX/D");
    Results->Branch("positivityHX",  &positivity_HX,  "positivityHX/I");
    Results->Branch("lthPX",         &lth_PX,         "lthPX/D");
    Results->Branch("lphPX",         &lph_PX,         "lphPX/D");
    Results->Branch("ltpPX",         &ltp_PX,         "ltpPX/D");
    Results->Branch("lthstarPX",     &lthstar_PX,     "lthstarPX/D");
    Results->Branch("lphstarPX",     &lphstar_PX,     "lphstarPX/D");
    Results->Branch("ltildePX",      &ltilde_PX,      "ltildePX/D");
    Results->Branch("positivityPX",  &positivity_PX,  "positivityPX/I");

    double err_lth_CS;
    double err_lph_CS;
    double err_ltp_CS;
    double err_lthstar_CS;
    double err_lphstar_CS;
    double err_ltilde_CS;
    double err_lth_HX;
    double err_lph_HX;
    double err_ltp_HX;
    double err_lthstar_HX;
    double err_lphstar_HX;
    double err_ltilde_HX;
    double err_lth_PX;
    double err_lph_PX;
    double err_ltp_PX;
    double err_lthstar_PX;
    double err_lphstar_PX;
    double err_ltilde_PX;

    Results->Branch("err_lthCS",         &err_lth_CS,         "err_lthCS/D");
    Results->Branch("err_lphCS",         &err_lph_CS,         "err_lphCS/D");
    Results->Branch("err_ltpCS",         &err_ltp_CS,         "err_ltpCS/D");
    Results->Branch("err_lthstarCS",     &err_lthstar_CS,     "err_lthstarCS/D");
    Results->Branch("err_lphstarCS",     &err_lphstar_CS,     "err_lphstarCS/D");
    Results->Branch("err_ltildeCS",      &err_ltilde_CS,      "err_ltildeCS/D");
    Results->Branch("err_lthHX",         &err_lth_HX,         "err_lthHX/D");
    Results->Branch("err_lphHX",         &err_lph_HX,         "err_lphHX/D");
    Results->Branch("err_ltpHX",         &err_ltp_HX,         "err_ltpHX/D");
    Results->Branch("err_lthstarHX",     &err_lthstar_HX,     "err_lthstarHX/D");
    Results->Branch("err_lphstarHX",     &err_lphstar_HX,     "err_lphstarHX/D");
    Results->Branch("err_ltildeHX",      &err_ltilde_HX,      "err_ltildeHX/D");
    Results->Branch("err_lthPX",         &err_lth_PX,         "err_lthPX/D");
    Results->Branch("err_lphPX",         &err_lph_PX,         "err_lphPX/D");
    Results->Branch("err_ltpPX",         &err_ltp_PX,         "err_ltpPX/D");
    Results->Branch("err_lthstarPX",     &err_lthstar_PX,     "err_lthstarPX/D");
    Results->Branch("err_lphstarPX",     &err_lphstar_PX,     "err_lphstarPX/D");
    Results->Branch("err_ltildePX",      &err_ltilde_PX,      "err_ltildePX/D");

/*    lth_CS=h_lth_CS->GetMean();
    err_lth_CS=h_lth_CS->GetRMS();
    lph_CS=h_lph_CS->GetMean();
    err_lph_CS=h_lph_CS->GetRMS();
    ltp_CS=h_ltp_CS->GetMean();
    err_ltp_CS=h_ltp_CS->GetRMS();
    lthstar_CS=h_lthstar_CS->GetMean();
    err_lthstar_CS=h_lthstar_CS->GetRMS();
    lphstar_CS=h_lphstar_CS->GetMean();
    err_lphstar_CS=h_lphstar_CS->GetRMS();
    ltilde_CS=h_ltilde_CS->GetMean();
    err_ltilde_CS=h_ltilde_CS->GetRMS();
    lth_HX=h_lth_HX->GetMean();
    err_lth_HX=h_lth_HX->GetRMS();
    lph_HX=h_lph_HX->GetMean();
    err_lph_HX=h_lph_HX->GetRMS();
    ltp_HX=h_ltp_HX->GetMean();
    err_ltp_HX=h_ltp_HX->GetRMS();
    lthstar_HX=h_lthstar_HX->GetMean();
    err_lthstar_HX=h_lthstar_HX->GetRMS();
    lphstar_HX=h_lphstar_HX->GetMean();
    err_lphstar_HX=h_lphstar_HX->GetRMS();
    ltilde_HX=h_ltilde_HX->GetMean();
    err_ltilde_HX=h_ltilde_HX->GetRMS();
    lth_PX=h_lth_PX->GetMean();
    err_lth_PX=h_lth_PX->GetRMS();
    lph_PX=h_lph_PX->GetMean();
    err_lph_PX=h_lph_PX->GetRMS();
    ltp_PX=h_ltp_PX->GetMean();
    err_ltp_PX=h_ltp_PX->GetRMS();
    lthstar_PX=h_lthstar_PX->GetMean();
    err_lthstar_PX=h_lthstar_PX->GetRMS();
    lphstar_PX=h_lphstar_PX->GetMean();
    err_lphstar_PX=h_lphstar_PX->GetRMS();
    ltilde_PX=h_ltilde_PX->GetMean();
    err_ltilde_PX=h_ltilde_PX->GetRMS();
*/
    double MPVerrorLow,MPVerrorHigh;
    int nSigma=1;

    FindMPV(h_lth_CS , lth_CS , MPVerrorLow, MPVerrorHigh, MPValgo, nSigma); err_lth_CS=(MPVerrorLow+MPVerrorHigh)/2.;
    FindMPV(h_lph_CS , lph_CS , MPVerrorLow, MPVerrorHigh, MPValgo, nSigma); err_lph_CS=(MPVerrorLow+MPVerrorHigh)/2.;
    FindMPV(h_ltp_CS , ltp_CS , MPVerrorLow, MPVerrorHigh, MPValgo, nSigma); err_ltp_CS=(MPVerrorLow+MPVerrorHigh)/2.;
    FindMPV(h_ltilde_CS , ltilde_CS , MPVerrorLow, MPVerrorHigh, MPValgo, nSigma);  err_ltilde_CS=(MPVerrorLow+MPVerrorHigh)/2.;
    FindMPV(h_lthstar_CS , lthstar_CS , MPVerrorLow, MPVerrorHigh, MPValgo, nSigma); err_lthstar_CS=(MPVerrorLow+MPVerrorHigh)/2.;
    FindMPV(h_lphstar_CS , lphstar_CS , MPVerrorLow, MPVerrorHigh, MPValgo, nSigma); err_lphstar_CS=(MPVerrorLow+MPVerrorHigh)/2.;

    FindMPV(h_lth_HX , lth_HX , MPVerrorLow, MPVerrorHigh, MPValgo, nSigma); err_lth_HX=(MPVerrorLow+MPVerrorHigh)/2.;
    FindMPV(h_lph_HX , lph_HX , MPVerrorLow, MPVerrorHigh, MPValgo, nSigma); err_lph_HX=(MPVerrorLow+MPVerrorHigh)/2.;
    FindMPV(h_ltp_HX , ltp_HX , MPVerrorLow, MPVerrorHigh, MPValgo, nSigma); err_ltp_HX=(MPVerrorLow+MPVerrorHigh)/2.;
    FindMPV(h_ltilde_HX , ltilde_HX , MPVerrorLow, MPVerrorHigh, MPValgo, nSigma);  err_ltilde_HX=(MPVerrorLow+MPVerrorHigh)/2.;
    FindMPV(h_lthstar_HX , lthstar_HX , MPVerrorLow, MPVerrorHigh, MPValgo, nSigma); err_lthstar_HX=(MPVerrorLow+MPVerrorHigh)/2.;
    FindMPV(h_lphstar_HX , lphstar_HX , MPVerrorLow, MPVerrorHigh, MPValgo, nSigma); err_lphstar_HX=(MPVerrorLow+MPVerrorHigh)/2.;

    FindMPV(h_lth_PX , lth_PX , MPVerrorLow, MPVerrorHigh, MPValgo, nSigma); err_lth_PX=(MPVerrorLow+MPVerrorHigh)/2.;
    FindMPV(h_lph_PX , lph_PX , MPVerrorLow, MPVerrorHigh, MPValgo, nSigma); err_lph_PX=(MPVerrorLow+MPVerrorHigh)/2.;
    FindMPV(h_ltp_PX , ltp_PX , MPVerrorLow, MPVerrorHigh, MPValgo, nSigma); err_ltp_PX=(MPVerrorLow+MPVerrorHigh)/2.;
    FindMPV(h_ltilde_PX , ltilde_PX , MPVerrorLow, MPVerrorHigh, MPValgo, nSigma);  err_ltilde_PX=(MPVerrorLow+MPVerrorHigh)/2.;
    FindMPV(h_lthstar_PX , lthstar_PX , MPVerrorLow, MPVerrorHigh, MPValgo, nSigma); err_lthstar_PX=(MPVerrorLow+MPVerrorHigh)/2.;
    FindMPV(h_lphstar_PX , lphstar_PX , MPVerrorLow, MPVerrorHigh, MPValgo, nSigma); err_lphstar_PX=(MPVerrorLow+MPVerrorHigh)/2.;


    cout<<"lth_CS     = "<<lth_CS<<" +/- "<<err_lth_CS<<endl;
    cout<<"lph_CS     = "<<lph_CS<<" +/- "<<err_lph_CS<<endl;
    cout<<"ltp_CS     = "<<ltp_CS<<" +/- "<<err_ltp_CS<<endl;
    cout<<"ltilde_CS  = "<<ltilde_CS<<" +/- "<<err_ltilde_CS<<endl;
    cout<<"lthstar_CS = "<<lthstar_CS<<" +/- "<<err_lthstar_CS<<endl;
    cout<<"lphstar_CS = "<<lphstar_CS<<" +/- "<<err_lphstar_CS<<endl;

    cout<<"lth_HX     = "<<lth_HX<<" +/- "<<err_lth_HX<<endl;
    cout<<"lph_HX     = "<<lph_HX<<" +/- "<<err_lph_HX<<endl;
    cout<<"ltp_HX     = "<<ltp_HX<<" +/- "<<err_ltp_HX<<endl;
    cout<<"ltilde_HX  = "<<ltilde_HX<<" +/- "<<err_ltilde_HX<<endl;
    cout<<"lthstar_HX = "<<lthstar_HX<<" +/- "<<err_lthstar_HX<<endl;
    cout<<"lphstar_HX = "<<lphstar_HX<<" +/- "<<err_lphstar_HX<<endl;

    cout<<"lth_PX     = "<<lth_PX<<" +/- "<<err_lth_PX<<endl;
    cout<<"lph_PX     = "<<lph_PX<<" +/- "<<err_lph_PX<<endl;
    cout<<"ltp_PX     = "<<ltp_PX<<" +/- "<<err_ltp_PX<<endl;
    cout<<"ltilde_PX  = "<<ltilde_PX<<" +/- "<<err_ltilde_PX<<endl;
    cout<<"lthstar_PX = "<<lthstar_PX<<" +/- "<<err_lthstar_PX<<endl;
    cout<<"lphstar_PX = "<<lphstar_PX<<" +/- "<<err_lphstar_PX<<endl;


    Results->Fill();


  test_lth_CS->Delete();
  test_lph_CS->Delete();
  test_ltp_CS->Delete();

  test_lth_HX->Delete();
  test_lph_HX->Delete();
  test_ltp_HX->Delete();

  test_lth_PX->Delete();
  test_lph_PX->Delete();
  test_ltp_PX->Delete();

  resultsFile->Write();

} // end of main