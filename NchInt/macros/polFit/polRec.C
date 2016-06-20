#include "Riostream.h"
#include "TRandom3.h"
#include "TSystem.h"
#include "TString.h"
#include "TROOT.h"
#include "TTree.h"
#include "TMath.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "TF1.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TFile.h"
#include "TGraphAsymmErrors.h"
#include "TGraph2D.h"
#include "TGraph.h"
#include "TEfficiency.h"

//#include "genDefs.h"
//#include "effsAndCuts.h"

bool isMuonInAcceptance(int iCut, double pT, double eta);
double singleLeptonEfficiency( double& pT, double& eta, int nEff, TFile *fInEff, TH2D* hEvalEff, bool MCeff, TEfficiency* TEff);
void EvaluateEffFileName(int nEff, char EffFileName [200], bool singleLeptonEff);
double DiLeptonEfficiency( double& Dilepton_pT, double& Dilepton_rap, int nDileptonEff, TFile *fInDileptonEff, bool MCeff);
double EvaluateRhoFactor( double& costh, double& phi, int nEff, TFile* fInRhoFactor, double rap, double pT, bool StatVarRho);
double DenominatorAmapEfficiency( double& pT, double& eta, int nDenominatorAmap, TFile *fInEff_nDenominatorAmap, TH2D* hEvalEff_nDenominatorAmap, bool MCeff, TEfficiency* TEff_nDenominatorAmap);
double EvaluateAmap( double& costh_Amap, double& phi_Amap, int nAmap, TFile* fInAmap, double rap, double pT);

void polRec(double rapdilepton_min = 1,
		double rapdilepton_max = 1,
		double pTdilepton_min = 1,
		double pTdilepton_max = 1,
		double mass_signal_peak  =  1,
		double mass_signal_sigma =  1,
		double n_sigmas_signal = 1,
		int nEff=1,
		int nDileptonEff=1,
		int nRhoFactor=1,
		int FidCuts=0,
		Char_t *dirstruct = "ToyDirectory_Default",
		bool applyFidCuts=false,
		Char_t *effDir = "effDir_Default",
		bool MCeff=false,
		bool MCDileptoneff=false,
		int rapBin=999,
		int pTbin=999,
		bool useAmapApproach=false,
		int nAmap=999,
		int nDenominatorAmap=999,
		bool StatVarTotBGfraction=false,
		bool StatVarRho=false){

	cout<<"/////////////////////////////////"<<endl;
	cout<<"running polRec.C ........////////"<<endl;
	cout<<"/////////////////////////////////"<<endl;

	gROOT->Reset();

	bool ForceEpsSmallerOne=false;

	//Get single Lepton Efficiency File name
	char EffFile[200];
	char EffFileName[200];
	EvaluateEffFileName(nEff,EffFileName,true);
	sprintf(EffFile,"%s/%s",effDir,EffFileName);

	char EffType[200];
	if(MCeff) sprintf(EffType,"hEff_MC_central");
	else sprintf(EffType,"hEff_DATA_central");
	cout<<"EffFile: "<<EffFile<<endl;

	TFile *fInEff = new TFile(EffFile);
	TH1* hEff=(TH1*) fInEff->Get(EffType);

	//Get DiLepton Efficiency File name

	EvaluateEffFileName(nDileptonEff,EffFileName,false);
	sprintf(EffFile,"%s/%s",effDir,EffFileName);

	TFile *fInDileptonEff = new TFile(EffFile);
	TH1* hDileptonEff=(TH1*) fInDileptonEff->Get(EffType);

	//Get rho factor File name
	EvaluateEffFileName(nRhoFactor,EffFileName,false);
	sprintf(EffFile,"%s/%s",effDir,EffFileName);
	TFile *fInRhoFactor = new TFile(EffFile);


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
			} //pTBin
		} //etaBin
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
		TEff=(TEfficiency*) fInEff->Get(EffType);
		TH1* hEffTOT=(TH1*)TEff->GetTotalHistogram();
		TH1* hEffPASS=(TH1*)TEff->GetPassedHistogram();
		hEffPASS->Divide(hEffTOT);
		hEvalEff=(TH2D*)hEffPASS->Clone("hEffPASS");
	}

	if(nEff!=105 && nEff!=106 && nEff!=1 && nEff!=1)
		hEvalEff = (TH2D*)hEvalEff1D->Clone("hEvalEff");
	if(nEff > 10000) hEvalEff = (TH2D*)hEvalEff2D->Clone("hEvalEff");

	//sprintf(hEvalEffName,"%s/EvalHisto.root",dirstruct);
	//hEvalEff->SaveAs(hEvalEffName);





	////// Get Amap

	EvaluateEffFileName(nAmap,EffFileName,false);
	sprintf(EffFile,"%s/%s",effDir,EffFileName);

	TFile *fInAmap = new TFile(EffFile);







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




















	double mass_min = mass_signal_peak - n_sigmas_signal*mass_signal_sigma;
	double mass_max = mass_signal_peak + n_sigmas_signal*mass_signal_sigma;

	char filename [500];

	// get ntuple of generated data from file

	sprintf(filename,"%s/genData.root",dirstruct);
	TFile* genFile = new TFile(filename,"READ");
	TTree* genData = (TTree*)genFile->Get("genData");


	// create output data file

	sprintf(filename,"%s/data.root",dirstruct);
	TFile* dataFile = new TFile(filename, "RECREATE", "dataFile");


	// output ntuple

	TTree* data = new TTree("selectedData","selectedData");


	// histograms to be output to the opened file

	// background distributions

	//const int nbinth   = 15;
	//const int nbinph   = 15;
	//// binning of costh and phi same as data Bg histogram
	const int nbinth   = ToyMC::binCosth[rapBin-1][pTbin-1];
	const int nbinph   = ToyMC::binPhi[rapBin-1][pTbin-1];
	const int nbinpT   =  7;
	const int nbinrap  =  2;
	const int nbinmass =  7;
	std::cout << "nbinth: " << nbinth << std::endl;
	std::cout << "nbinph: " << nbinph << std::endl;

	TH2D* background_costhphiPX = new TH2D( "background_costhphiPHX", "", nbinth,    -1.,    1.,
			nbinph,   -180.,  180.  );
	TH3D* background_pTrapMass  = new TH3D( "background_pTrapMass",   "", nbinpT,   pTdilepton_min,  pTdilepton_max,
			nbinrap,  rapdilepton_min, rapdilepton_max,
			nbinmass, mass_min,        mass_max         );
	TH1D* background_fraction = new TH1D( "background_fraction", "", 1, 0., 1. );

	// this is a temporary histogram to calculate BG fraction after acceptence and efficiency cuts
	TH1D* isBGdistribution    = new TH1D( "isBGdistribution", "", 2, 0., 2. );


	// structure of input ntuple

	TLorentzVector* lepP_gen = 0;    genData->SetBranchAddress( "lepP",    &lepP_gen );
	TLorentzVector* lepN_gen = 0;    genData->SetBranchAddress( "lepN",    &lepN_gen );

	double costh_CS;  genData->SetBranchAddress( "costh_CS",     &costh_CS );
	double phi_CS;    genData->SetBranchAddress( "phi_CS",       &phi_CS   );
	double phith_CS;  genData->SetBranchAddress( "phith_CS",     &phith_CS );

	double costh_HX;  genData->SetBranchAddress( "costh_HX",     &costh_HX );
	double phi_HX;    genData->SetBranchAddress( "phi_HX",       &phi_HX   );
	double phith_HX;  genData->SetBranchAddress( "phith_HX",     &phith_HX );

	double costh_PX;  genData->SetBranchAddress( "costh_PX",     &costh_PX );
	double phi_PX;    genData->SetBranchAddress( "phi_PX",       &phi_PX   );
	double phith_PX;  genData->SetBranchAddress( "phith_PX",     &phith_PX );

	double mass;      genData->SetBranchAddress( "mass",         &mass     );
	double pT;        genData->SetBranchAddress( "pT",           &pT       );
	double rap;       genData->SetBranchAddress( "rap",          &rap      );

	int isBG;         genData->SetBranchAddress( "isBG",         &isBG     );


	// structure of output ntuple

	TLorentzVector* lepP = new TLorentzVector(0.,0.,0.,0.);  data->Branch( "lepP", "TLorentzVector", &lepP );
	TLorentzVector* lepN = new TLorentzVector(0.,0.,0.,0.);  data->Branch( "lepN", "TLorentzVector", &lepN );


	// loop over events in the input ntuple

	int numEvts = int( genData->GetEntries() );


	cout << endl;
	cout << "Reading " << numEvts << " dilepton events"<< endl;
	cout << "------------------------------------------------------------" << endl;
	cout << "Progress: "<<endl;

	int n_step = numEvts/5;
	int n_step_=1;
	int rejected=0;

	for ( int evtNumber = 0; evtNumber < numEvts; evtNumber++ ) {

		if ((evtNumber+1)%n_step == 0) {cout << n_step_*20 <<" % "<<endl; n_step_++;}

		genData->GetEvent( evtNumber );




		// select data in acceptance and apply efficiency

		double lepP_pT  = lepP_gen->Pt();
		double lepN_pT  = lepN_gen->Pt();
		double lepP_eta  = lepP_gen->PseudoRapidity();
		double lepN_eta  = lepN_gen->PseudoRapidity();

		bool isEventAccepted = isMuonInAcceptance( FidCuts-1, lepP_pT, lepP_eta )
			* isMuonInAcceptance( FidCuts-1, lepN_pT, lepN_eta );

		//    if (lepP_gen->Phi()>0 && lepP_gen->Phi()<TMath::Pi()/8.) isEventAccepted=false;
		//    if (lepN_gen->Phi()>0 && lepN_gen->Phi()<TMath::Pi()/8.) isEventAccepted=false;

		if ( !isEventAccepted ) {rejected++; continue;}

		double effP = singleLeptonEfficiency( lepP_pT, lepP_eta, nEff, fInEff , hEvalEff , MCeff , TEff);
		double effN = singleLeptonEfficiency( lepN_pT, lepN_eta, nEff, fInEff , hEvalEff , MCeff , TEff);

		double costh_DILEff;
		double phi_DILEff;
		if(nDileptonEff>200 && nDileptonEff<211) {costh_DILEff=costh_CS; phi_DILEff=phi_CS;}
		if(nDileptonEff>210 && nDileptonEff<221) {costh_DILEff=costh_HX; phi_DILEff=phi_HX;}
		if(nDileptonEff>220 && nDileptonEff<231) {costh_DILEff=costh_PX; phi_DILEff=phi_PX;}
		double DileptonEff = DiLeptonEfficiency( costh_DILEff, phi_DILEff, nDileptonEff, fInDileptonEff, MCDileptoneff);

		double costh_RhoFactor;
		double phi_RhoFactor;
		if(nRhoFactor>300 && nRhoFactor<311) {costh_RhoFactor=costh_CS; phi_RhoFactor=phi_CS;}
		if(nRhoFactor>310 && nRhoFactor<321) {costh_RhoFactor=costh_HX; phi_RhoFactor=phi_HX;}
		if(nRhoFactor>320 && nRhoFactor<331) {costh_RhoFactor=costh_PX; phi_RhoFactor=phi_PX;}
		double RhoFactor = EvaluateRhoFactor( costh_RhoFactor, phi_RhoFactor, nRhoFactor, fInRhoFactor, rap, pT, StatVarRho);

		double epsilon = effP*effN;


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

		/*    double rndmeffP = gRandom->Uniform(1.);
					double rndmeffN = gRandom->Uniform(1.);
					double rndmDileptoneff = gRandom->Uniform(1.);
					double rndmRhoFactor = gRandom->Uniform(1.);
		 */
		double rndmepsilon = gRandom->Uniform(1.);

		if ( rndmepsilon > epsilon ) {rejected++; continue;}

		// fill background histograms and output ntuple

		isBGdistribution->Fill( isBG );
		if ( isBG ) {
			background_costhphiPX->Fill( costh_PX, phi_PX );
			background_pTrapMass->Fill( pT, TMath::Abs(rap), mass );
		}

		*lepP = *lepP_gen;
		*lepN = *lepN_gen;

		data->Fill();

	} // end of RD ntuple loop

	cout << endl << endl;



	// background fraction

	double f_BG = isBGdistribution->GetMean();


	background_fraction->SetBinContent( 1, f_BG );
	if(StatVarTotBGfraction)
		background_fraction->SetBinError( 1,  ToyMC::fracBackgrounderr[rapBin-1][pTbin-1]);
	else
		background_fraction->SetBinError( 1,  0.);

	double effFrac=(double(numEvts-rejected))/double(numEvts);


	cout<<"Effective Background Fraction:           "<<f_BG<<endl;
	cout<<"Fraction of Events Surviving Efficiency: "<<effFrac<<endl;
	cout<<"Surviving Signal Events:                 "<<isBGdistribution->GetEntries()*(1-f_BG)<<endl;

	// end

	genFile->Close();
	dataFile->Write();
	dataFile->Close();

	}