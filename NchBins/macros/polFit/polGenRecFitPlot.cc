#include "rootIncludes.inc"
#include "commonVar.h"
#include "ToyMC.h"
#include "effsAndCuts.h"

#include "TSystem.h"
#include "TROOT.h"
#include "polGen.C"
#include "polRec.C"
#include "polFit.C"
#include "polPlot.C"
#include "TGraphAsymmErrors.h"
#include "TFile.h"

#include <time.h>


//====================================

int main(int argc, char** argv) {


	int nGenerations=999;
	int polScenSig=9999;
	int frameSig=999;
	int polScenBkg=999;
	int frameBkg=999;
	int rapBinMin=999;
	int rapBinMax=999;
	int ptBinMin=999;
	int ptBinMax=999;
	int cpmBinMin=999;
	int cpmBinMax=999;	
	int nEff=999;
	int nRhoFactor=999;
	int nDileptonEff=999;
	int nRecEff=999;
	int nRecDileptonEff=999;
	int nRecRhoFactor=999;
	int FidCuts=999;
	int nSample=999;
	int ConstEvents=999;
	int nSkipGen=999;
	int ThisGen=999;
	int MPValgo=0;
	int nState=0;
	int nAmap=999;
	int nDenominatorAmap=999;

	bool ConstEvents_(false);
	bool gen(false);
	bool rec(false);
	bool fit(false);
	bool plot(false);
	bool RealData(false);
	bool UseDifferingEff(false);
	bool MCeff(false);
	bool MCReceff(false);
	bool MCDileptoneff(false);
	bool MCDileptonReceff(false);
	bool scalePlots(false);
	bool NewAccCalc=true;
	bool deletePseudoData=false;
	bool useAmapApproach=false;
	bool useBatch=false;
	bool StatVarTotBGfraction=false;
	bool StatVarTotBGmodel=false;
	bool StatVarRho=false;
	bool cutDeltaREllDpt=false;

	double n_sigmas_signal = 3.;

	Char_t *storagedir = "Default"; //Storage Directory
	Char_t *basedir = "Default"; //Code Directory
	Char_t *JobID = "Default";
	Char_t *realdatadir = "Default"; //Storage Directory
	Char_t *TreeID = "ToyMC"; //Storage Directory
	Char_t *TreeBinID_dataFile = "ToyMC"; //Storage Directory

	for( int i=0;i < argc; ++i ) {
		if(std::string(argv[i]).find("ptBinMin") != std::string::npos) {char* ptBinMinchar = argv[i]; char* ptBinMinchar2 = strtok (ptBinMinchar, "p"); ptBinMin = atof(ptBinMinchar2); cout<<"ptBinMin = "<<ptBinMin<<endl;}
		if(std::string(argv[i]).find("ptBinMax") != std::string::npos) {char* ptBinMaxchar = argv[i]; char* ptBinMaxchar2 = strtok (ptBinMaxchar, "p"); ptBinMax = atof(ptBinMaxchar2); cout<<"ptBinMax = "<<ptBinMax<<endl;}
		if(std::string(argv[i]).find("cpmBinMin") != std::string::npos) {char* cpmBinMinchar = argv[i]; char* cpmBinMinchar2 = strtok (cpmBinMinchar, "p"); cpmBinMin = atof(cpmBinMinchar2); cout<<"cpmBinMin = "<<cpmBinMin<<endl;}
		if(std::string(argv[i]).find("cpmBinMax") != std::string::npos) {char* cpmBinMaxchar = argv[i]; char* cpmBinMaxchar2 = strtok (cpmBinMaxchar, "p"); cpmBinMax = atof(cpmBinMaxchar2); cout<<"cpmBinMax = "<<cpmBinMax<<endl;}
		if(std::string(argv[i]).find("rapBinMin") != std::string::npos) {char* rapBinMinchar = argv[i]; char* rapBinMinchar2 = strtok (rapBinMinchar, "p"); rapBinMin = atof(rapBinMinchar2); cout<<"rapBinMin = "<<rapBinMin<<endl;}
		if(std::string(argv[i]).find("rapBinMax") != std::string::npos) {char* rapBinMaxchar = argv[i]; char* rapBinMaxchar2 = strtok (rapBinMaxchar, "p"); rapBinMax = atof(rapBinMaxchar2); cout<<"rapBinMax = "<<rapBinMax<<endl;}
		if(std::string(argv[i]).find("nGenerations") != std::string::npos) {char* nGenerationschar = argv[i]; char* nGenerationschar2 = strtok (nGenerationschar, "p"); nGenerations = atof(nGenerationschar2); cout<<"nGenerations = "<<nGenerations<<endl;}
		if(std::string(argv[i]).find("frameSig") != std::string::npos) {char* framecharSig = argv[i]; char* framecharSig2 = strtok (framecharSig, "p"); frameSig = atof(framecharSig2); cout<<"frameSig = "<<frameSig<<endl;}
		if(std::string(argv[i]).find("polScenSig") != std::string::npos) {char* polScencharSig = argv[i]; char* polScencharSig2 = strtok (polScencharSig, "p"); polScenSig = atof(polScencharSig2); cout<<"polScenSig = "<<polScenSig<<endl;}
		if(std::string(argv[i]).find("frameBkg") != std::string::npos) {char* framecharBkg = argv[i]; char* framecharBkg2 = strtok (framecharBkg, "p"); frameBkg = atof(framecharBkg2); cout<<"frameBkg = "<<frameBkg<<endl;}
		if(std::string(argv[i]).find("polScenBkg") != std::string::npos) {char* polScencharBkg = argv[i]; char* polScencharBkg2 = strtok (polScencharBkg, "p"); polScenBkg = atof(polScencharBkg2); cout<<"polScenBkg = "<<polScenBkg<<endl;}
		if(std::string(argv[i]).find("nEff") != std::string::npos) {char* nEffchar = argv[i]; char* nEffchar2 = strtok (nEffchar, "p"); nEff = atof(nEffchar2); cout<<"nEff = "<<nEff<<endl;}
		if(std::string(argv[i]).find("nDiEff") != std::string::npos) {char* nDileptonEffchar = argv[i]; char* nDileptonEffchar2 = strtok (nDileptonEffchar, "p"); nDileptonEff = atof(nDileptonEffchar2); cout<<"nDileptonEff = "<<nDileptonEff<<endl;}
		if(std::string(argv[i]).find("nRhoFactor") != std::string::npos) {char* nRhoFactorchar = argv[i]; char* nRhoFactorchar2 = strtok (nRhoFactorchar, "p"); nRhoFactor = atof(nRhoFactorchar2); cout<<"nRhoFactor = "<<nRhoFactor<<endl;}
		if(std::string(argv[i]).find("FidCuts") != std::string::npos) {char* FidCutschar = argv[i]; char* FidCutschar2 = strtok (FidCutschar, "p"); FidCuts = atof(FidCutschar2); cout<<"FidCuts = "<<FidCuts<<endl;}
		if(std::string(argv[i]).find("nSample") != std::string::npos) {char* nSamplechar = argv[i]; char* nSamplechar2 = strtok (nSamplechar, "p"); nSample = atof(nSamplechar2); cout<<"nSample = "<<nSample<<endl;}
		if(std::string(argv[i]).find("ConstEvents") != std::string::npos) {char* ConstEventschar = argv[i]; char* ConstEventschar2 = strtok (ConstEventschar, "p"); ConstEvents = atof(ConstEventschar2); cout<<"ConstEvents = "<<ConstEvents<<endl;}
		if(std::string(argv[i]).find("nSkipGen") != std::string::npos) {char* nSkipGenchar = argv[i]; char* nSkipGenchar2 = strtok (nSkipGenchar, "p"); nSkipGen = atof(nSkipGenchar2); cout<<"nSkipGen = "<<nSkipGen<<endl;}
		if(std::string(argv[i]).find("ThisGen") != std::string::npos) {char* ThisGenchar = argv[i]; char* ThisGenchar2 = strtok (ThisGenchar, "p"); ThisGen = atof(ThisGenchar2); cout<<"ThisGen = "<<ThisGen<<endl;}
		if(std::string(argv[i]).find("nAmap") != std::string::npos) {char* nAmapchar = argv[i]; char* nAmapchar2 = strtok (nAmapchar, "p"); nAmap = atof(nAmapchar2); cout<<"nAmap = "<<nAmap<<endl;}
		if(std::string(argv[i]).find("nDenominatorAmap") != std::string::npos) {char* nDenominatorAmapchar = argv[i]; char* nDenominatorAmapchar2 = strtok (nDenominatorAmapchar, "p"); nDenominatorAmap = atof(nDenominatorAmapchar2); cout<<"nDenominatorAmap = "<<nDenominatorAmap<<endl;}

		if(std::string(argv[i]).find("nSigma") != std::string::npos) {char* nSigmachar = argv[i]; char* nSigmachar2 = strtok (nSigmachar, "p"); n_sigmas_signal = atof(nSigmachar2); cout<<"nSigma = "<<n_sigmas_signal<<endl;}
		if(std::string(argv[i]).find("nState") != std::string::npos) {char* nStatechar = argv[i]; char* nStatechar2 = strtok (nStatechar, "p"); nState = atof(nStatechar2); cout<<"nState = "<<nState<<endl;}

		if(std::string(argv[i]).find("NewAccCalc=false") != std::string::npos) {NewAccCalc=false; cout<<"use old acceptance correction method"<<endl;}
		if(std::string(argv[i]).find("UseConstEv=true") != std::string::npos) {ConstEvents_=true; cout<<"use constant number of reconstructed events"<<endl;}
		if(std::string(argv[i]).find("deletePseudoData=true") != std::string::npos) {deletePseudoData=true; cout<<"deletePseudoData"<<endl;}
		if(std::string(argv[i]).find("gen=true") != std::string::npos) {gen=true; cout<<"run polGen.C"<<endl;}
		if(std::string(argv[i]).find("rec=true") != std::string::npos) {rec=true; cout<<"run polRec.C"<<endl;}
		if(std::string(argv[i]).find("fit=true") != std::string::npos) {fit=true; cout<<"run polFit.C"<<endl;}
		if(std::string(argv[i]).find("plot=true") != std::string::npos) {plot=true; cout<<"run polPlot.C"<<endl;}
		if(std::string(argv[i]).find("UseDifferingEff=true") != std::string::npos) {UseDifferingEff=true; cout<<"Using differing efficiency definitions for generation and fitting"<<endl;}
		if(std::string(argv[i]).find("nRecEff") != std::string::npos) {char* nRecEffchar = argv[i]; char* nRecEffchar2 = strtok (nRecEffchar, "p"); nRecEff = atof(nRecEffchar2); cout<<"nRecEff = "<<nRecEff<<endl;}
		if(std::string(argv[i]).find("nRecDiEff") != std::string::npos) {char* nRecDileptonEffchar = argv[i]; char* nRecDileptonEffchar2 = strtok (nRecDileptonEffchar, "p"); nRecDileptonEff = atof(nRecDileptonEffchar2); cout<<"nRecDileptonEff = "<<nRecDileptonEff<<endl;}
		if(std::string(argv[i]).find("nRecRhoFactor") != std::string::npos) {char* nRecRhoFactorchar = argv[i]; char* nRecRhoFactorchar2 = strtok (nRecRhoFactorchar, "p"); nRecRhoFactor = atof(nRecRhoFactorchar2); cout<<"nRecRhoFactor = "<<nRecRhoFactor<<endl;}
		if(std::string(argv[i]).find("scalePlots=true") != std::string::npos) {scalePlots=true; cout<<"run polGen.C"<<endl;}
		if(std::string(argv[i]).find("useAmapApproach=true") != std::string::npos) {useAmapApproach=true; cout<<"use new A-map approach to calculate dimuon efficiencies"<<endl;}
		if(std::string(argv[i]).find("useBatch=1") != std::string::npos) {useBatch=true; cout<<"use batch submission system"<<endl;}
		if(std::string(argv[i]).find("StatVarTotBGfraction=1") != std::string::npos) {StatVarTotBGfraction=true; cout<<"apply statistical fluctuations on f_background"<<endl;}
		if(std::string(argv[i]).find("StatVarTotBGmodel=1") != std::string::npos) {StatVarTotBGmodel=true; cout<<"apply statistical fluctuations on Bg model"<<endl;}
		if(std::string(argv[i]).find("StatVarRho=1") != std::string::npos) {StatVarRho=true; cout<<"apply statistical fluctuations on rho factor"<<endl;}
		if(std::string(argv[i]).find("cutDeltaREllDpt=1") != std::string::npos) {cutDeltaREllDpt=true; cout<<"cut DeltaREllDpt"<<endl;}

		if(std::string(argv[i]).find("JobID") != std::string::npos) {char* JobIDchar = argv[i]; char* JobIDchar2 = strtok (JobIDchar, "="); JobID = JobIDchar2; cout<<"JobID = "<<JobID<<endl;}
		if(std::string(argv[i]).find("basedir") != std::string::npos) {char* basedirchar = argv[i]; char* basedirchar2 = strtok (basedirchar, "="); basedir = basedirchar2; cout<<"basedir = "<<basedir<<endl;}
		if(std::string(argv[i]).find("storagedir") != std::string::npos) {char* storagedirchar = argv[i]; char* storagedirchar2 = strtok (storagedirchar, "="); storagedir = storagedirchar2; cout<<"storagedir = "<<storagedir<<endl;}

		if(std::string(argv[i]).find("realdatadir") != std::string::npos) {RealData=true; cout<<"Fit to real data"<<endl; char* realdatadirchar = argv[i]; char* realdatadirchar2 = strtok (realdatadirchar, "="); realdatadir = realdatadirchar2; cout<<"realdatadir = "<<realdatadir<<endl;}
		if(std::string(argv[i]).find("TreeID") != std::string::npos) {char* TreeIDchar = argv[i]; char* TreeIDchar2 = strtok (TreeIDchar, "="); TreeID = TreeIDchar2; cout<<"TreeID = "<<TreeID<<endl;}
		if(std::string(argv[i]).find("UseMCeff=true") != std::string::npos) {MCeff=true; cout<<"use MC efficiency"<<endl;}
		if(std::string(argv[i]).find("UseMCReceff=true") != std::string::npos) {MCReceff=true; cout<<"use MC efficiency for generation"<<endl;}
		if(std::string(argv[i]).find("UseMCDileptoneff=true") != std::string::npos) {MCDileptoneff=true; cout<<"use MC Dilepton efficiency"<<endl;}
		if(std::string(argv[i]).find("UseMCDileptonReceff=true") != std::string::npos) {MCDileptonReceff=true; cout<<"use MC Dilepton efficiency for generation"<<endl;}
		if(std::string(argv[i]).find("MPValgo") != std::string::npos) {char* MPValgochar = argv[i]; char* MPValgochar2 = strtok (MPValgochar, "p"); MPValgo = atof(MPValgochar2); cout<<"MPValgo = "<<MPValgo<<endl;}

	}


	double mass_signal_peak;
	double mass_signal_sigma;

	if(nState < 4 && rapBinMin==1) mass_signal_sigma=0.075;
	if(nState < 4 && rapBinMin==2) mass_signal_sigma=0.1;

	if(nState == 4 && rapBinMin==1) mass_signal_sigma=0.025;
	if(nState == 4 && rapBinMin==2) mass_signal_sigma=0.035;

	if(nState == 5 && rapBinMin==1) mass_signal_sigma=0.026;
	if(nState == 5 && rapBinMin==2) mass_signal_sigma=0.038;
	if(nState == 5 && rapBinMin==3) mass_signal_sigma=0.048;


	double lambda_theta_sig_;
	double lambda_phi_sig_;
	double lambda_thetaphi_sig_;

	double lambda_theta_bkg_ = ToyMC::ScenarioBkg[0][polScenBkg-1];
	double lambda_phi_bkg_ = ToyMC::ScenarioBkg[1][polScenBkg-1];
	double lambda_thetaphi_bkg_ = ToyMC::ScenarioBkg[2][polScenBkg-1];


	bool injectRealDataPolarization(false);
	if(polScenSig==999) injectRealDataPolarization=true;

	if(!injectRealDataPolarization){
		lambda_theta_sig_ = ToyMC::ScenarioSig[0][polScenSig-1];
		lambda_phi_sig_ = ToyMC::ScenarioSig[1][polScenSig-1];
		lambda_thetaphi_sig_ = ToyMC::ScenarioSig[2][polScenSig-1];
	}

	if(injectRealDataPolarization){
		cout<<"injectRealDataPolarization"<<endl;

		char filename[1000];
		sprintf(filename,"/afs/hephy.at/scratch/k/knuenz/CMSSW_4_2_4_patch2/src/UpsilonPol/macros/polFit/Systematics/TotalSyst/May20Centrals_CentralsFromAlteredPPDMay20_1SigmaStatError_FracHighCorrected/TGraphResults_Psi%dS.root",nState);
		cout<<filename<<endl;
		TFile *infile1 = new TFile(filename,"READ");

		char GraphName[1000];


		int nRapBins = 2; if(nState == 5) nRapBins=3;

		for(int rapBin = 1; rapBin < nRapBins+1; rapBin++){



			TGraphAsymmErrors* graph_lth;
			TGraphAsymmErrors* graph_lph;
			TGraphAsymmErrors* graph_ltp;
			TGraphAsymmErrors* graph_lthstar;
			TGraphAsymmErrors* graph_lphstar;
			TGraphAsymmErrors* graph_ltilde;

			if(frameSig==1)  {
				cout<<GraphName<<endl;
				sprintf(GraphName,"lth_CS_rap%d",rapBin);
				graph_lth = (TGraphAsymmErrors*) infile1->Get(GraphName);
				sprintf(GraphName,"lph_CS_rap%d",rapBin);
				graph_lph = (TGraphAsymmErrors*) infile1->Get(GraphName);
				sprintf(GraphName,"ltp_CS_rap%d",rapBin);
				graph_ltp = (TGraphAsymmErrors*) infile1->Get(GraphName);
				sprintf(GraphName,"lthstar_CS_rap%d",rapBin);
				graph_lthstar = (TGraphAsymmErrors*) infile1->Get(GraphName);
				sprintf(GraphName,"lphstar_CS_rap%d",rapBin);
				graph_lphstar = (TGraphAsymmErrors*) infile1->Get(GraphName);
				sprintf(GraphName,"ltilde_CS_rap%d",rapBin);
				graph_ltilde = (TGraphAsymmErrors*) infile1->Get(GraphName);
			}

			if(frameSig==2)  {
				sprintf(GraphName,"lth_HX_rap%d",rapBin);
				graph_lth = (TGraphAsymmErrors*) infile1->Get(GraphName);
				sprintf(GraphName,"lph_HX_rap%d",rapBin);
				graph_lph = (TGraphAsymmErrors*) infile1->Get(GraphName);
				sprintf(GraphName,"ltp_HX_rap%d",rapBin);
				graph_ltp = (TGraphAsymmErrors*) infile1->Get(GraphName);
				sprintf(GraphName,"lthstar_HX_rap%d",rapBin);
				graph_lthstar = (TGraphAsymmErrors*) infile1->Get(GraphName);
				sprintf(GraphName,"lphstar_HX_rap%d",rapBin);
				graph_lphstar = (TGraphAsymmErrors*) infile1->Get(GraphName);
				sprintf(GraphName,"ltilde_HX_rap%d",rapBin);
				graph_ltilde = (TGraphAsymmErrors*) infile1->Get(GraphName);
			}

			if(frameSig==3)  {
				sprintf(GraphName,"lth_PX_rap%d",rapBin);
				graph_lth = (TGraphAsymmErrors*) infile1->Get(GraphName);
				sprintf(GraphName,"lph_PX_rap%d",rapBin);
				graph_lph = (TGraphAsymmErrors*) infile1->Get(GraphName);
				sprintf(GraphName,"ltp_PX_rap%d",rapBin);
				graph_ltp = (TGraphAsymmErrors*) infile1->Get(GraphName);
				sprintf(GraphName,"lthstar_PX_rap%d",rapBin);
				graph_lthstar = (TGraphAsymmErrors*) infile1->Get(GraphName);
				sprintf(GraphName,"lphstar_PX_rap%d",rapBin);
				graph_lphstar = (TGraphAsymmErrors*) infile1->Get(GraphName);
				sprintf(GraphName,"ltilde_PX_rap%d",rapBin);
				graph_ltilde = (TGraphAsymmErrors*) infile1->Get(GraphName);
			}

			cout<<"TGraphs of all parameters loaded for frame "<<frameSig<<endl;



			double cpmCentre;
			double lth_lmean;
			double lph_lmean;
			double ltp_lmean;



			graph_lth->GetPoint(cpmBinMin-1,cpmCentre,lth_lmean);

			graph_lph->GetPoint(cpmBinMin-1,cpmCentre,lph_lmean);

			graph_ltp->GetPoint(cpmBinMin-1,cpmCentre,ltp_lmean);


			cout<<"Values of all parameters loaded for frame "<<frameSig<<endl;



			if(rapBin==rapBinMin){
				lambda_theta_sig_ = lth_lmean;
				lambda_phi_sig_ = lph_lmean;
				lambda_thetaphi_sig_ = ltp_lmean;
			}

		}

		cout<<"Using Real Data Results as Input for toyMC-samples, rap"<<rapBinMin<<"_pT"<<ptBinMin<<"_cpm"<<cpmBinMin<<", injected in frame "<<frameSig<<endl;
		cout<<"lth = "<<lambda_theta_sig_<<endl;
		cout<<"lph = "<<lambda_phi_sig_<<endl;
		cout<<"ltp = "<<lambda_thetaphi_sig_<<endl;

	}



	if(nState==1) mass_signal_peak=9.5;
	if(nState==2) mass_signal_peak=10.;
	if(nState==3) mass_signal_peak=10.4;
	if(nState==4) mass_signal_peak=3.097;
	if(nState==5) mass_signal_peak=3.686;

	if(!UseDifferingEff) {nRecEff=nEff; nRecDileptonEff=nDileptonEff; MCReceff=MCeff; MCDileptonReceff=MCDileptoneff; nRecRhoFactor=nRhoFactor;}

	Char_t *OutputDirectory;
	Char_t *TreeBinID;
	double f_BG;
	int n_events;
	char basestruct[1000],substruct[1000], dirstruct[1000], rapptstruct[1000], filenameFrom[1000], filenameTo[1000] , tmpfilename[1000], TreeBinID_[1000], effDir[1000];

	sprintf(basestruct,"%s/%s",storagedir,JobID);gSystem->mkdir(basestruct);
	sprintf(substruct,"%s/Sig_frame%dscen%d_Bkg_frame%dscen%d",basestruct,frameSig,polScenSig,frameBkg,polScenBkg); if(!RealData) gSystem->mkdir(substruct);
	sprintf(effDir,"%s/macros/polFit/EffFiles",basedir);

	cout<<"storagedir: "<<storagedir<<endl;
	cout<<"JobID: "<<JobID<<endl;
	cout<<"basestruct: "<<basestruct<<endl;
	cout<<"substruct: "<<substruct<<endl;

	time_t seconds; seconds = time (NULL); double time_0=seconds; double time_1;


	int iRap = rapBinMin;
	int iPt = ptBinMin;
	int icpm = cpmBinMin;
	sprintf(TreeBinID_,"%s_rap%d_pT%d_cpm%d",TreeID,iRap,iPt,icpm);
	char TreeBinID_dataFileChar[200];
	sprintf(TreeBinID_dataFileChar,"%s",TreeBinID_);
	TreeBinID_dataFile=TreeBinID_dataFileChar;
	if(useBatch) sprintf(TreeBinID_,"Fit%d_%s_rap%d_pT%d_cpm%d",ThisGen,TreeID,iRap,iPt,icpm);
	TreeBinID=TreeBinID_;

	cout<<TreeBinID_dataFile<<endl;
	cout<<TreeBinID<<endl;


	double cpmlow=onia::cpmRange[icpm-1];
	double cpmhigh=onia::cpmRange[icpm];
	double ptlow=onia::pTRange[iRap][iPt-1];
	double pthigh=onia::pTRange[iRap][iPt];
	double raplow=onia::rapForPTRange[iRap-1];
	double raphigh=onia::rapForPTRange[iRap];

	sprintf(rapptstruct,"%s/rap%d_pT%d_cpm%d",substruct,iRap,iPt,icpm); if(!RealData) gSystem->mkdir(rapptstruct);

	if(gen) {sprintf(filenameFrom,"%s/polGen.C",substruct);					sprintf(filenameTo,"%s/polGen.C",rapptstruct); 							gSystem->CopyFile(filenameFrom,filenameTo,kTRUE);}
	if(rec) {sprintf(filenameFrom,"%s/polRec.C",substruct);					sprintf(filenameTo,"%s/polRec.C",rapptstruct); 							gSystem->CopyFile(filenameFrom,filenameTo,kTRUE);}
	if(fit) {sprintf(filenameFrom,"%s/polFit.C",substruct);					sprintf(filenameTo,"%s/polFit.C",rapptstruct); 							gSystem->CopyFile(filenameFrom,filenameTo,kTRUE);}
	if(plot) {sprintf(filenameFrom,"%s/polPlot.C",substruct);				sprintf(filenameTo,"%s/polPlot.C",rapptstruct); 						gSystem->CopyFile(filenameFrom,filenameTo,kTRUE);}

	sprintf(filenameFrom,"%s/effsAndCuts.h",substruct);					sprintf(filenameTo,"%s/effsAndCuts.h",rapptstruct); 							gSystem->CopyFile(filenameFrom,filenameTo,kTRUE);

	gSystem->cd(rapptstruct);

	//if(gen) gROOT->ProcessLine(".L polGen.C+");
	//if(rec) gROOT->ProcessLine(".L polRec.C+");
	//if(fit) gROOT->ProcessLine(".L polFit.C+");
	//if(plot) gROOT->ProcessLine(".L polPlot.C+");


	/// Extract number of signal and background events to be generated, as well as f_BG to be generated to result in desired effective f_BG:

	int nTargetEvents;
	nTargetEvents = ToyMC::numEvents[iRap-1][iPt-1][icpm-1];

	if(gen){

		int numEvCheck = 500000;
		f_BG = ToyMC::fracBackground[iRap-1][iPt-1][icpm-1];
		sprintf(tmpfilename,"%s/data.root",rapptstruct);
		TFile* dataFile = new TFile(tmpfilename, "READ");
		cout<<"f_BG: "<<f_BG<<endl;

		if(dataFile->Get("isBGdistribution")==NULL){
			OutputDirectory=rapptstruct;
			polGen(raplow,raphigh,ptlow,pthigh,mass_signal_peak,mass_signal_sigma,n_sigmas_signal,numEvCheck,f_BG,lambda_theta_sig_,lambda_phi_sig_,lambda_thetaphi_sig_,lambda_theta_bkg_,lambda_phi_bkg_,lambda_thetaphi_bkg_,frameSig,frameBkg,-999,nState,OutputDirectory);
			if(rec)polRec(raplow,raphigh,ptlow,pthigh,mass_signal_peak,mass_signal_sigma,n_sigmas_signal,nRecEff,nRecDileptonEff,nRecRhoFactor,FidCuts,OutputDirectory, true, effDir, MCReceff, MCDileptonReceff, iRap, iPt, icpm, useAmapApproach, nAmap, nDenominatorAmap);
			sprintf(tmpfilename,"%s/genData.root",rapptstruct);			gSystem->Unlink(tmpfilename);
			sprintf(tmpfilename,"%s/GenResults.root",rapptstruct);		gSystem->Unlink(tmpfilename);
		}

		sprintf(tmpfilename,"%s/data.root",rapptstruct);
		dataFile = new TFile(tmpfilename, "READ");
		TH1D* isBG_distribution = (TH1D*)dataFile->Get("isBGdistribution");

		double sigFact = isBG_distribution->GetBinContent(1)/(numEvCheck*(1-f_BG));
		double bkgFact = isBG_distribution->GetBinContent(2)/(numEvCheck*f_BG);

		dataFile->Close();


		if(ConstEvents_) nTargetEvents = ConstEvents;

		n_events = nTargetEvents/sigFact+(nTargetEvents/(1-f_BG)-nTargetEvents)/bkgFact;

		f_BG = (n_events-nTargetEvents/sigFact)/n_events;
		cout<<"n_events: "<<n_events<<endl;

	}
	/// Start actual Generation and Fits:
	cout<<"Start actual Generation and Fits....."<<endl;

	int iGen = ThisGen;
	int nTotalFits = nSkipGen+nGenerations;

	seconds = time (NULL); time_1=seconds;
	//sprintf(dirstruct,"%s/Generation%d",rapptstruct,iGen+nSkipGen); if(!RealData) gSystem->mkdir(dirstruct);
	sprintf(dirstruct,"%s/Generation%d",rapptstruct,iGen); if(!RealData) gSystem->mkdir(dirstruct);
	OutputDirectory=dirstruct;
	if(RealData) OutputDirectory=basestruct;

	cout<<"nState: "<<nState<<endl;
	cout<<"OutputDirectory: "<<OutputDirectory<<endl;
	cout<<"basestruct: "<<basestruct<<endl;

	if(gen)polGen(raplow,raphigh,ptlow,pthigh,mass_signal_peak,mass_signal_sigma,n_sigmas_signal,n_events,f_BG,lambda_theta_sig_,lambda_phi_sig_,lambda_thetaphi_sig_,lambda_theta_bkg_,lambda_phi_bkg_,lambda_thetaphi_bkg_,frameSig,frameBkg,iGen,nState,OutputDirectory);
	if(rec)polRec(raplow,raphigh,ptlow,pthigh,mass_signal_peak,mass_signal_sigma,n_sigmas_signal,nRecEff,nRecDileptonEff,nRecRhoFactor,FidCuts,OutputDirectory, false, effDir, MCReceff, MCDileptonReceff, iRap, iPt, icpm, useAmapApproach, nAmap, nDenominatorAmap, StatVarTotBGfraction, StatVarRho);
	if(fit)polFit(nSample,FidCuts, nEff, nDileptonEff, nRhoFactor, OutputDirectory, realdatadir, TreeBinID, TreeBinID_dataFile, RealData, effDir, MCeff, MCDileptoneff, iRap, iPt, NewAccCalc, MPValgo, useAmapApproach, nAmap, nDenominatorAmap, StatVarTotBGfraction, StatVarTotBGmodel, StatVarRho, cutDeltaREllDpt);
	if(plot)polPlot(OutputDirectory, TreeBinID, RealData, MPValgo, scalePlots, nTotalFits, nState, ptlow, pthigh, raplow, raphigh);

	//sprintf(dirstruct,"%s/Generation%d",rapptstruct,iGen+nSkipGen);
	sprintf(dirstruct,"%s/Generation%d",rapptstruct,iGen);
	if(deletePseudoData){
		sprintf(tmpfilename,"%s/genData.root",dirstruct);			gSystem->Unlink(tmpfilename);
		sprintf(tmpfilename,"%s/data.root",dirstruct);				gSystem->Unlink(tmpfilename);
		sprintf(tmpfilename,"%s/efficiency.root",dirstruct);		gSystem->Unlink(tmpfilename);
	}

	seconds = time (NULL);

	if(fit) cout<<"Proccessing time for this generation: "<<seconds-time_1<<" s"<<" (Corresponding to "<<(seconds-time_1)/60<<" min)"<<endl;
	if(fit) cout<<"Per signal event in final sample:     "<<(seconds-time_1)/nTargetEvents*1000<<" ms"<<endl;

	if(gen)  {sprintf(tmpfilename,"%s/polGen.C",rapptstruct);			gSystem->Unlink(tmpfilename);		sprintf(tmpfilename,"%s/polGen_C.d",rapptstruct);			gSystem->Unlink(tmpfilename);		sprintf(tmpfilename,"%s/polGen_C.so",rapptstruct);			gSystem->Unlink(tmpfilename);}
	if(rec)  {sprintf(tmpfilename,"%s/polRec.C",rapptstruct);			gSystem->Unlink(tmpfilename);		sprintf(tmpfilename,"%s/polRec_C.d",rapptstruct);			gSystem->Unlink(tmpfilename);		sprintf(tmpfilename,"%s/polRec_C.so",rapptstruct);			gSystem->Unlink(tmpfilename);}
	if(fit)  {sprintf(tmpfilename,"%s/polFit.C",rapptstruct);			gSystem->Unlink(tmpfilename);		sprintf(tmpfilename,"%s/polFit_C.d",rapptstruct);			gSystem->Unlink(tmpfilename);		sprintf(tmpfilename,"%s/polFit_C.so",rapptstruct);			gSystem->Unlink(tmpfilename);}
	if(plot)  {sprintf(tmpfilename,"%s/polPlot.C",rapptstruct);			gSystem->Unlink(tmpfilename);		sprintf(tmpfilename,"%s/polPlot_C.d",rapptstruct);			gSystem->Unlink(tmpfilename);		sprintf(tmpfilename,"%s/polPlot_C.so",rapptstruct);			gSystem->Unlink(tmpfilename);}
	if(fit && RealData)  {sprintf(tmpfilename,"%s/polFit_C.d",basestruct);			gSystem->Unlink(tmpfilename);		sprintf(tmpfilename,"%s/polFit_C.so",basestruct);			gSystem->Unlink(tmpfilename);}
	if(plot && RealData)  {sprintf(tmpfilename,"%s/polPlot_C.d",basestruct);			gSystem->Unlink(tmpfilename);		sprintf(tmpfilename,"%s/polPlot_C.so",basestruct);			gSystem->Unlink(tmpfilename);}

	sprintf(tmpfilename,"%s/effsAndCuts.h",rapptstruct);			gSystem->Unlink(tmpfilename);

	return 0;
}


