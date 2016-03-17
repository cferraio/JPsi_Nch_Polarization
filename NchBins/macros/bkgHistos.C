#include "calculatePar.cc"
#include "calcPol.C"

#include <string>
#include <iostream>
#include <sstream>

using namespace RooFit;
using namespace std;

int order(int n);
int findEvenNum(double number);
vector<double> calculateInte(RooWorkspace *ws, RooDataSet *dataJpsictErr, double ctCutMin, double ctCutMax);
vector<double> calculateInteSB(RooWorkspace *ws, RooDataSet *dataJpsictErr, double ctCutMin, double ctCutMax);
TH2D *subtract2D(TH2D* hist1, TH2D* hist2);
TH3D *subtract3D(TH3D* hist1, TH3D* hist2);
TH2D* ReSetBin(TH2D* hist, int nBinX, int nBinY, const std::stringstream& name, const std::stringstream& title);

double calcuFracL(RooWorkspace *ws, double mean, double sigma);

//---------------------------------------------------------------------------------------------------------------
void bkgHistos(const std::string infilename, int rapBin, int ptBin, int cpmBin, int nState, bool folding, bool MC, bool doCtauUncer, bool PolLSB, bool PolRSB, bool PolNP, int ctauScen, int FracLSB, bool forceBinning, bool normApproach, bool scaleFracBg, char *polDataPath){

	const std::string
		datafilename = "tmpFiles/selEvents_data.root",
								 treename = "selectedData",
								 wsname = "ws_masslifetime";

	// input
	TFile *datafile = TFile::Open(datafilename.c_str());
	if(!datafile){
		std::cout << "Inputfile missing" << std::endl;
		return;
	}
	TTree *intree = (TTree *)datafile->Get(treename.c_str());
	TLorentzVector *lepP = 0, *lepN = 0, *jpsi = 0;
	double cpmval = 0;

	TFile *fitfile = TFile::Open(infilename.c_str());
	if(!fitfile){
		std::cout << "fitfile is missing" << std::endl;
		return;
	}
	RooWorkspace *ws = (RooWorkspace*)fitfile->Get(wsname.c_str());
	if(!ws){
		std::cout << "workspace not found" << std::endl;
		return;
	}

	gStyle->SetPadRightMargin(0.2);
	gROOT->SetStyle("Plain");
	gStyle->SetTitleBorderSize(0);

	// create output
	std::stringstream outfilename;
	outfilename << "tmpFiles/data_Psi" << nState-3 << "S_rap" << rapBin << "_pT" << ptBin << "_cpm" << cpmBin << ".root";
	TFile *output = TFile::Open(outfilename.str().c_str(), "RECREATE");

	TTree *outtree = intree->CloneTree(0);

	// background histos
	TH2D *hBG_cosThetaPhiL[onia::kNbFrames];
	TH2D *hBG_cosThetaPhiR[onia::kNbFrames];
	TH2D *hBG_cosThetaPhi[onia::kNbFrames];
	TH2D *hNPBG_cosThetaPhi[onia::kNbFrames];
	TH2D *hNPS_cosThetaPhi[onia::kNbFrames];
	TH2D *hBGinNP_cosThetaPhi[onia::kNbFrames];
	TH2D *hBGinNP_cosThetaPhiL[onia::kNbFrames];
	TH2D *hBGinNP_cosThetaPhiR[onia::kNbFrames];
	TH2D *hTBG_cosThetaPhi[onia::kNbFrames];
	TH2D *hSR_cosThetaPhiL[onia::kNbFrames];
	TH2D *hSR_cosThetaPhiR[onia::kNbFrames];
	TH2D *hSR_cosThetaPhi[onia::kNbFrames];

	for(int iFrame = 0; iFrame < onia::kNbFrames; iFrame++){
		//book the 2D (cosTheta, phi) histos for the L and R mass sideband
		std::stringstream nameL, nameR, nameNP, nameBGinNPL, nameBGinNPR, nameSRL, nameSRR, nameSR, title;
		nameL << "hBG_cosThetaPhi_" << onia::frameLabel[iFrame] << "_L";
		nameR << "hBG_cosThetaPhi_" << onia::frameLabel[iFrame] << "_R";
		nameNP << "hNPBG_cosThetaPhi_" << onia::frameLabel[iFrame];
		nameBGinNPL << "hBGinNP_cosThetaPhi_" << onia::frameLabel[iFrame] << "_L";
		nameBGinNPR << "hBGinNP_cosThetaPhi_" << onia::frameLabel[iFrame] << "_R";
		nameSRL << "hSR_cosThetaPhi_" << onia::frameLabel[iFrame] << "_L";
		nameSRR << "hSR_cosThetaPhi_" << onia::frameLabel[iFrame] << "_R";
		nameSR << "hSR_cosThetaPhi_" << onia::frameLabel[iFrame];
		title << ";cos#theta_{"<< onia::frameLabel[iFrame] << "};#phi_{" << onia::frameLabel[iFrame] << "} [deg]";
		hBG_cosThetaPhiL[iFrame] = new TH2D(nameL.str().c_str(), title.str().c_str(),
				onia::kNbBinsCosT, onia::cosTMin, onia::cosTMax, onia::kNbBinsPhiPol, onia::phiPolMin, onia::phiPolMax);
		hBG_cosThetaPhiL[iFrame]->Sumw2();
		hBG_cosThetaPhiR[iFrame] = new TH2D(nameR.str().c_str(), title.str().c_str(),
				onia::kNbBinsCosT, onia::cosTMin, onia::cosTMax, onia::kNbBinsPhiPol, onia::phiPolMin, onia::phiPolMax);
		hBG_cosThetaPhiR[iFrame]->Sumw2();

		hNPBG_cosThetaPhi[iFrame] = new TH2D(nameNP.str().c_str(), title.str().c_str(),
				onia::kNbBinsCosT, onia::cosTMin, onia::cosTMax, onia::kNbBinsPhiPol, onia::phiPolMin, onia::phiPolMax);
		hNPBG_cosThetaPhi[iFrame]->Sumw2();

		hBGinNP_cosThetaPhiL[iFrame] = new TH2D(nameBGinNPL.str().c_str(), title.str().c_str(),
				onia::kNbBinsCosT, onia::cosTMin, onia::cosTMax, onia::kNbBinsPhiPol, onia::phiPolMin, onia::phiPolMax);
		hBGinNP_cosThetaPhiL[iFrame]->Sumw2();
		hBGinNP_cosThetaPhiR[iFrame] = new TH2D(nameBGinNPR.str().c_str(), title.str().c_str(),
				onia::kNbBinsCosT, onia::cosTMin, onia::cosTMax, onia::kNbBinsPhiPol, onia::phiPolMin, onia::phiPolMax);
		hBGinNP_cosThetaPhiR[iFrame]->Sumw2();
		hSR_cosThetaPhi[iFrame] = new TH2D(nameSR.str().c_str(), title.str().c_str(),
				onia::kNbBinsCosT, onia::cosTMin, onia::cosTMax, onia::kNbBinsPhiPol, onia::phiPolMin, onia::phiPolMax);
		hSR_cosThetaPhi[iFrame]->Sumw2();

	} // iFrame

	// mean pT and y histos (background-subtracted)
	int nBins = 100;
	TH1D* cpm_L   = new TH1D( "cpmLSB", "cpmLSB", nBins, onia::cpmRange[cpmBin-1], onia::cpmRange[cpmBin]);
	TH1D* cpm_R   = new TH1D( "cpmRSB", "cpmRSB", nBins, onia::cpmRange[cpmBin-1], onia::cpmRange[cpmBin]);
	TH1D* cpm_highct_L   = new TH1D( "cpm_highct_LSB", "cpm_highct_LSB", nBins, onia::cpmRange[cpmBin-1], onia::cpmRange[cpmBin]);
	TH1D* cpm_highct_R   = new TH1D( "cpm_highcta_RSB", "cpm_highct_RSB", nBins, onia::cpmRange[cpmBin-1], onia::cpmRange[cpmBin]);
	TH1D* cpm_NP   = new TH1D( "cpmNP", "cpmNP", nBins, onia::cpmRange[cpmBin-1], onia::cpmRange[cpmBin]);
	TH1D* cpm_PSR   = new TH1D( "cpmPSR", "cpmPSR", nBins, onia::cpmRange[cpmBin-1], onia::cpmRange[cpmBin]);
	
	TH1D* pT_L   = new TH1D( "pTLSB", "pTLSB", nBins, onia::pTRange[rapBin-1][ptBin-1], onia::pTRange[rapBin-1][ptBin]);
	TH1D* pT_R   = new TH1D( "pTRSB", "pTRSB", nBins, onia::pTRange[rapBin-1][ptBin-1], onia::pTRange[rapBin-1][ptBin]);
	TH1D* pT_highct_L   = new TH1D( "pT_highct_LSB", "pT_highct_LSB", nBins, onia::pTRange[rapBin-1][ptBin-1], onia::pTRange[rapBin-1][ptBin]);
	TH1D* pT_highct_R   = new TH1D( "pT_highcta_RSB", "pT_highct_RSB", nBins, onia::pTRange[rapBin-1][ptBin-1], onia::pTRange[rapBin-1][ptBin]);
	TH1D* pT_NP   = new TH1D( "pTNP", "pTNP", nBins, onia::pTRange[rapBin-1][ptBin-1], onia::pTRange[rapBin-1][ptBin]);
	TH1D* pT_PSR   = new TH1D( "pTPSR", "pTPSR", nBins, onia::pTRange[rapBin-1][ptBin-1], onia::pTRange[rapBin-1][ptBin]);
	
	TH1D* rap_L   = new TH1D( "rapLSB", "rapLSB", nBins, onia::rapForPTRange[rapBin-1], onia::rapForPTRange[rapBin]);
	TH1D* rap_R   = new TH1D( "rapRSB", "rapRSB", nBins, onia::rapForPTRange[rapBin-1], onia::rapForPTRange[rapBin]);
	TH1D* rap_highct_L   = new TH1D( "rap_highct_LSB", "rap_highct_LSB", nBins, onia::rapForPTRange[rapBin-1], onia::rapForPTRange[rapBin]);
	TH1D* rap_highct_R   = new TH1D( "rap_highct_RSB", "rap_highct_RSB", nBins, onia::rapForPTRange[rapBin-1], onia::rapForPTRange[rapBin]);
	TH1D* rap_NP   = new TH1D( "rapNP", "rapNP",nBins, onia::rapForPTRange[rapBin-1], onia::rapForPTRange[rapBin]);
	TH1D* rap_PSR   = new TH1D( "rapPSR", "rapPSR", nBins, onia::rapForPTRange[rapBin-1], onia::rapForPTRange[rapBin]);

	//------------------------------------------------------------------------------------------------
	// store pT and y borders
	TVectorD* cpmBorder = new TVectorD(1, 2, onia::cpmRange[cpmBin-1], onia::cpmRange[cpmBin], "END");
	TVectorD* pTBorder = new TVectorD(1, 2, onia::pTRange[rapBin-1][ptBin-1], onia::pTRange[rapBin-1][ptBin], "END");
	TVectorD* yBorder = new TVectorD(1, 2, onia::rapForPTRange[rapBin-1], onia::rapForPTRange[rapBin], "END");
	cpmBorder->Print();
	pTBorder->Print();
	yBorder->Print();
	output->cd();
	cpmBorder->Write();
	pTBorder->Write();
	yBorder->Write();

	// set branches
	intree->SetBranchAddress("lepP", &lepP);
	intree->SetBranchAddress("lepN", &lepN);
	intree->SetBranchAddress("JpsiP", &jpsi);
	intree->SetBranchAddress("CPM",&cpmval);
	double jpsict = 0;
	intree->SetBranchAddress("Jpsict", &jpsict);

	//---------------------------------------------------------------------------------------------
	// INPUT FROM FIT
	// MASS FIT
	// reading fitting parameters
	std::stringstream mfitresult;
	mfitresult << "m_fitresult_rap" << rapBin << "_pt" << ptBin << "_cpm" << cpmBin;
	RooFitResult *mresult = (RooFitResult*)ws->genobj(mfitresult.str().c_str());
	RooArgList mvarlist = mresult->floatParsFinal();

	// parameters
	RooRealVar *CBmass=(RooRealVar*)mvarlist.find("CBmass");
	RooRealVar *CBsigma1=(RooRealVar*)mvarlist.find("CBsigma");
	RooRealVar *CBsigma2=(RooRealVar*)mvarlist.find("CBsigma2");
	RooRealVar *fracCB1=(RooRealVar*)mvarlist.find("fracCB1");
	double mass = CBmass->getVal();
	double Sigma1 = CBsigma1->getVal();
	double Sigma2 = CBsigma2->getVal();
	double fCB1 = fracCB1->getVal();
	double sigma = sqrt( pow(Sigma1,2)*fCB1 + pow(Sigma2,2)*(1-fCB1) );

	// variables
	RooRealVar *m = ws->var("JpsiMass");
	RooRealVar *bkgLambda = ws->var("bkgLambda");
	RooRealVar *CBalpha = ws->var("CBalpha");
	RooRealVar *CBn = ws->var("CBn");
	RooRealVar *ct = ws->var("Jpsict");
	RooRealVar *cterr = ws->var("JpsictErr");

	// pdf
	RooAbsPdf *bkgMass = (RooAbsPdf*)ws->pdf("bkgMassShape");
	RooAbsPdf *signalMass = (RooAbsPdf*)ws->pdf("sigMassShape"); //massFull
	RooAbsPdf *PRpdf = (RooAbsPdf*)ws->pdf("TotalPromptLifetime");
	RooAbsPdf *NPpdf = (RooAbsPdf*)ws->pdf("nonPromptSSD");
	RooAbsPdf *BGpdf = (RooAbsPdf*)ws->pdf("backgroundlifetime");
	RooAbsPdf *BGpdfLSB = (RooAbsPdf*)ws->pdf("backgroundlifetimeLpre");
	RooAbsPdf *LSBpdf = (RooAbsPdf*)ws->pdf("backgroundlifetimeL");
	RooAbsPdf *BGpdfRSB = (RooAbsPdf*)ws->pdf("backgroundlifetimeRpre");
	RooAbsPdf *RSBpdf = (RooAbsPdf*)ws->pdf("backgroundlifetimeR");
	RooAbsPdf *lifetime = (RooAbsPdf*)ws->pdf("fulllifetime");

	// load snapshot with all results
	std::stringstream masssnapshotname;
	masssnapshotname << "m_snapshot_rap" << rapBin << "_pt" << ptBin << "_cpm" << cpmBin;
	ws->loadSnapshot(masssnapshotname.str().c_str());

	TF1* funcBG = (TF1*)bkgMass->asTF(*m, *bkgLambda, *m);
	TF1* funcSig = (TF1*)signalMass->asTF(*m, RooArgList(*CBalpha, *fracCB1, *CBn, *CBsigma1, *CBsigma2, *CBmass), *m);

	//-------------------------------------------------------------------------------------------------
	// mass edges for left and right sideband
	double massMinL = onia::massMin;
	double massMaxR = onia::massMax;
	double massMaxL = 0,
				 massMinR = 0,
				 massMinSR = 0,
				 massMaxSR = 0;
	if(MC){
		//massMaxL = mass - sigma;
		//massMinR = mass + sigma;
		massMaxL = mass - 10*sigma;
		massMinR = mass + 10*sigma;
		massMinSR = massMaxL;
		massMaxSR = massMinR;
		std::cout << "Using 10 sigma for MC" << std::endl;
	}
	else{
		massMaxL = mass - onia::nSigBkgLow  * sigma;
		massMinR = mass + onia::nSigBkgHigh * sigma;
		massMinSR = mass - onia::nSigMass * sigma;
		massMaxSR = mass + onia::nSigMass * sigma;
	}
	std::cout << "-------------------------------------------------------------\n" <<
		"left  sideband: mass window " << massMinL  << " < M < " << massMaxL  << " GeV\n" <<
		"right sideband: mass window " << massMinR  << " < M < " << massMaxR  << " GeV\n" <<
		"signal  region: mass window " << massMinSR  << " < M < " << massMaxSR  << " GeV\n" <<
		"-------------------------------------------------------------\n" << std::endl;

	// LIFETIME FIT
	// get dataset and number of events in dataset
	std::stringstream dataBin, data, cutSR;
	dataBin << "data_rap" << rapBin << "_pt" << ptBin << "_cpm" << cpmBin << "_SR";
	data    << "data_rap" << rapBin << "_pt" << ptBin << "_cpm" << cpmBin;
	cutSR   << "JpsiMass > " << massMinSR << " && JpsiMass < " << massMaxSR;

	RooAbsData* Data    = ws->data(data.str().c_str());
	RooDataSet *binData = (RooDataSet*)Data->reduce(
			Cut(cutSR.str().c_str()),
			Name(dataBin.str().c_str()),
			Title("Data For Fitting"));
	int entries = binData->numEntries();

	// load snapshot with all results
	std::stringstream snapshotname;
	snapshotname << "l_snapshot_rap" << rapBin << "_pt" << ptBin << "_cpm" << cpmBin;
	ws->loadSnapshot(snapshotname.str().c_str());

	// get parameters for calculating fraction of left sideband
	RooRealVar *mean_ws = (RooRealVar*)ws->var("MeanSR");
	double mean = mean_ws->getVal();
	RooRealVar *mean_LSB_ws = (RooRealVar*)ws->var("MeanSBL");
	double mean_LSB = mean_LSB_ws->getVal();
	RooRealVar *mean_RSB_ws = (RooRealVar*)ws->var("MeanSBR");
	double mean_RSB = mean_RSB_ws->getVal();
	//double fracLSB = 1. - (mean - mean_LSB)/(mean_RSB - mean_LSB); // using the mean value
	double fracLSB = calcuFracL(ws, mass, sigma); // using the median value
	if(FracLSB!=-1) fracLSB = double(FracLSB)/100.;

	RooRealVar* fBkg_ws = (RooRealVar*)ws->var("FracBkg");
	double fBkg = fBkg_ws->getVal();

	// reading fitting parameters
	std::stringstream fitresult;
	fitresult << "l_fitresult_rap" << rapBin << "_pt" << ptBin << "_cpm" << cpmBin;
	RooFitResult *result = (RooFitResult*)ws->genobj(fitresult.str().c_str());
	RooArgList varlist = result->floatParsFinal();

	// get parameters
	RooRealVar* fP_ws = (RooRealVar*)varlist.find("fPrompt");
	double fP1 = fP_ws->getVal();
	double fNP1 = 1 - fBkg - fP1;
	// for signal contamination
	RooAbsReal* fPLSB_ws = (RooAbsReal*)ws->function("fPromptSBL");
	double fPLSB = fPLSB_ws->getVal();
	RooRealVar* fBkgLSB_ws = (RooRealVar*)varlist.find("fBkgSBL");
	double fBkgLSB = fBkgLSB_ws->getVal();
	double fNPLSB = 1. - fPLSB - fBkgLSB;
	RooAbsReal* fPRSB_ws = (RooAbsReal*)ws->function("fPromptSBR");
	double fPRSB = fPRSB_ws->getVal();
	RooRealVar* fBkgRSB_ws = (RooRealVar*)varlist.find("fBkgSBR");
	double fBkgRSB = fBkgRSB_ws->getVal();
	double fNPRSB = 1. - fPRSB - fBkgRSB;
	std::cout << "-------------------------------------------------------------\n"
		<< "fit parameters for signal contamination: \n" 
		<< "LSB: fP = " << fPLSB << ", fBkg = " << fBkgLSB << ", fNP = " << fNPLSB << "\n"
		<< "RSB: fP = " << fPRSB << ", fBkg = " << fBkgRSB << ", fNP = " << fNPRSB << "\n"
		<< "fP/fNP = " << fPLSB/fNPLSB << " (LSB) = " << fPRSB/fNPRSB << " (RSB) \n"
		<< "-------------------------------------------------------------" << std::endl;

	// calculate entries in signal region
	double nP = fP1*entries;
	double nNP = fNP1*entries;
	double nBG = fBkg*entries;

	std::cout << "-------------------------------------------------------------\n" <<
		"total number of events in signal region: " << entries << "\n" <<
		"prompt events in signal region: " << nP << "\n" <<
		"non prompt events in signal region: " << nNP << "\n" <<
		"background events in signal region: " << nBG << "\n" <<
		"fraction of left sideband: " << fracLSB << "\n" <<
		"-------------------------------------------------------------\n" << std::endl;

	int NumEvt = 0., maxEvt = 100000;
	if(Data->numEntries() < maxEvt)
		NumEvt = Data->numEntries();
	else
		NumEvt = maxEvt;

	RooDataSet *dataJpsictErr = (RooDataSet*)Data->reduce(SelectVars(RooArgSet(*cterr)),
			EventRange(0,NumEvt),Name("dataJpsictErr"));

	double meanPt = getMeanPt(rapBin,ptBin,cpmBin,infilename.c_str());

	// define sigma of prompt p.d.f., got from fit the trend
	// define function y = a + b * pT
	double a = 0.073, b = 0.0027;
	//proper decay length
	double L_decay = a + b * meanPt;
	//pseudo-proper decay length
	double l_pdecay = L_decay * mass / meanPt ;

	cout<<"l_pdecay: "<<l_pdecay<<endl;

	// used before new lifetime model
	//double l_pdecay = 0., scale = 1.;
	//if(nState==4) l_pdecay = L_decay * mass / meanPt ;
	//if(nState==5) {
	//    if(rapBin == 1) scale = 1.19; // values for new lifetime model; old value: 1.21 - first bin (7-10) was deleted
	//    if(rapBin == 2) scale = 1.27; // old value 1.29
	//    if(rapBin == 3) scale = 1.28; // old value: 1.31
	//    l_pdecay = scale * L_decay * 3.092 / meanPt ;
	//}

	double nSigma = 0.;

	if(ctauScen==0){ //default values
		if(nState==4) nSigma = 2.5;
		if(nState==5) nSigma = 2.0;
		cout<<"nSigma: "<<nSigma<<endl;
	}
	else if(ctauScen==1){
		if(nState==4) nSigma = 3.5;
		if(nState==5) nSigma = 3.0;
		cout<<"nSigma: "<<nSigma<<endl;
	}
	else if(ctauScen==2){
		if(nState==4) nSigma = 1.5;
		if(nState==5) nSigma = 1.0;
		cout<<"nSigma: "<<nSigma<<endl;
	}
	else if(ctauScen==4){
		if(nState==4) nSigma = 1.0;
		if(nState==5) nSigma = 1.0;
		cout<<"nSigma: "<<nSigma<<endl;
	}
	else if(ctauScen==5){
		if(nState==4) nSigma = 3.0;
		if(nState==5) nSigma = 3.0;
		cout<<"nSigma: "<<nSigma<<endl;
	}
	else if(ctauScen==6){
		if(nState==4) nSigma = 2.0;
		if(nState==5) nSigma = 2.0;
		cout<<"nSigma: "<<nSigma<<endl;
	}
	else if(ctauScen==7){
		if(nState==4) nSigma = 3.0;
		if(nState==5) nSigma = 3.0;
		cout<<"nSigma: "<<nSigma<<endl;
	}
	else{ //default values
		if(nState==4) nSigma = 2.5;
		if(nState==5) nSigma = 2.0;
		cout<<"nSigma: "<<nSigma<<endl;
	} //more can be added if need

	double ctauCut = nSigma*l_pdecay;

	if(ctauScen==3){
		ctauCut = 0.1 ; // mm 
	}

	cout << "ctauCut: " << ctauCut << endl;

	double ctCutMinPR = -ctauCut, ctCutMaxPR = ctauCut;
	double ctCutMinNP =  ctauCut, ctCutMaxNP = 6.;
	if(ctauScen==7){ ctCutMinNP = 5. * l_pdecay; }


	// histogram method
	// calculate fractions in prompt signal region
	//vector<double> InteRltsPR = calculateInte(ws,dataJpsictErr,ctCutMinPR,ctCutMaxPR);
	//double fPinP  = InteRltsPR[0];
	//double fNPinP = InteRltsPR[1];
	//double fBGinP = InteRltsPR[2];

	//vector<double> InteRltsSB = calculateInteSB(ws,dataJpsictErr,ctCutMinPR,ctCutMaxPR);
  //double fBkgLSBtotInP = InteRltsSB[2];
	//double fBkgLSBinP    = InteRltsSB[0];
	//double fBkgRSBtotInP = InteRltsSB[3];
	//double fBkgRSBinP    = InteRltsSB[1];
	//double fSRinPLSB = 1 - (fBkgLSBinP * fBkgLSB)/fBkgLSBtotInP;
	//double fSRinPRSB = 1 - (fBkgRSBinP * fBkgRSB)/fBkgRSBtotInP;

	// integral method
	// calculate fractions in the prompt signal region
	ct->setRange("PR", ctCutMinPR, ctCutMaxPR);
	// prompt fraction in signal region
	RooAbsReal* i_fPR = PRpdf->createIntegral(RooArgSet(*ct), NormSet(RooArgSet(*ct)), Range("PR"));
	double fPinP = i_fPR->getVal();
	// non prompt fraction in signal region
	RooAbsReal* i_fNPR = NPpdf->createIntegral(RooArgSet(*ct), NormSet(RooArgSet(*ct)), Range("PR"));
	double fNPinP = i_fNPR->getVal();
	// background fraction in signal region
	RooAbsReal* i_fBkg = BGpdf->createIntegral(RooArgSet(*ct), NormSet(RooArgSet(*ct)), Range("PR"));
	double fBGinP = i_fBkg->getVal();

	// calculate signal contamination in prompt sidebands
	// LSB
	RooAbsReal* i_fBkgLSBtotInP = LSBpdf->createIntegral(RooArgSet(*ct), NormSet(RooArgSet(*ct)), Range("PR"));
	RooAbsReal* i_fBkgLSBinP = BGpdfLSB->createIntegral(RooArgSet(*ct), NormSet(RooArgSet(*ct)), Range("PR"));
	// fSRinPLSB = ftot_LSB_ctau - fBkg_LSB_ctau * fBkgLSB
	double fSRinPLSB = 1 - (i_fBkgLSBinP->getVal() * fBkgLSB)/i_fBkgLSBtotInP->getVal();
	// RSB
	RooAbsReal* i_fBkgRSBtotInP = RSBpdf->createIntegral(RooArgSet(*ct), NormSet(RooArgSet(*ct)), Range("PR"));
	RooAbsReal* i_fBkgRSBinP = BGpdfRSB->createIntegral(RooArgSet(*ct), NormSet(RooArgSet(*ct)), Range("PR"));
	double fSRinPRSB = 1 - (i_fBkgRSBinP->getVal() * fBkgRSB)/i_fBkgRSBtotInP->getVal();

	std::cout << "-------------------------------------------------------------\n" <<
		nSigma << " sigma used for ctau cut" << "\n" <<
		"pseudo proper decay length: " << l_pdecay << "\n" << 
		"ctau cut at " << ctauCut << "\n" <<
		"fraction of prompt events within ctau cut: " << fPinP << " \n" <<
		"fraction of non prompt events within ctau cut: " << fNPinP << " \n" <<
		"fraction of background events within ctau cut: " << fBGinP << " \n" <<
		"-------------------------------------------------------------\n" <<
		"signal contamination in \n" <<
		"LSB: " << fSRinPLSB << "\n" <<
		"RSB: " << fSRinPRSB << "\n" <<
		"-------------------------------------------------------------\n" << std::endl;

	//------------------------------------------------------------------------------------------------
	// histogram method
	// calculate fractions in non prompt signal region
	//vector<double> InteRltsNP = calculateInte(ws,dataJpsictErr,ctCutMinNP,ctCutMaxNP);
	//double fPbkg  = InteRltsNP[0];
	//double fNPbkg = InteRltsNP[1];
	//double fBGbkg = InteRltsNP[2];

	// integal method
	// integral of all PDFs over non prompt region
	ct->setRange("NPR", ctCutMinNP, ctCutMaxNP);
	// prompt fraction in non prompt region
	RooAbsReal* i_fPRinNP = PRpdf->createIntegral(RooArgSet(*ct), NormSet(RooArgSet(*ct)), Range("NPR"));
	double fPbkg = i_fPRinNP->getVal();
	// non prompt fraction in signal region
	RooAbsReal* i_fNPRinNP = NPpdf->createIntegral(RooArgSet(*ct), NormSet(RooArgSet(*ct)), Range("NPR"));
	double fNPbkg = i_fNPRinNP->getVal();
	// background fraction in signal region
	RooAbsReal* i_fBkgInNP = BGpdf->createIntegral(RooArgSet(*ct), NormSet(RooArgSet(*ct)), Range("NPR"));
	double fBGbkg = i_fBkgInNP->getVal();

	double fBGinNP = 0;  // comb. background fraction in high ctau signal region
	double fP      = 0;  // prompt fraction in prompt signal region
	double fNPB    = 0;  // non prompt background events in prompt region
	double fBGsig  = 0;  // background in prompt signal region
	double fTBGsig = 0;  // total background fraction

	// for MC set fractions to 0.001 except prompt fraction
	if(MC || PolLSB || PolRSB){
		fBGinNP = 0.001;
		fP      = 0.999;
		fNPB    = 0.001;
		fBGsig  = 0.001;
		fSRinPLSB = 0.000001;
		fSRinPRSB = 0.000001;
		if(PolLSB) fTBGsig = fSRinPLSB;
		else if(PolRSB) fTBGsig = fSRinPRSB;
		else fTBGsig = 0.001;
	}
	// use real fractions for data
	else{
		fBGinNP = fBGbkg*fBkg/(fBGbkg*fBkg + fNPbkg*fNP1 + fPbkg*fP1);
		fP      =   fPinP*fP1/(fNPinP*fNP1 + fPinP*fP1 + fBGinP*fBkg);
		fNPB    = fNPinP*fNP1/(fNPinP*fNP1 + fPinP*fP1 + fBGinP*fBkg);
		fBGsig  = fBGinP*fBkg/(fNPinP*fNP1 + fPinP*fP1 + fBGinP*fBkg);
		fTBGsig =  fNPB + fBGsig;
		// in case of measuring non prompt polarization
		if(PolNP){
			fTBGsig = fBGinNP;
		}
	}

	//------------------------------------------------------------------------------------------------
	// evaluate error on non-prompt fraction
	double fPerr = 0., fNPerr = 0., fBGerr = 0.;
	if(!MC && !PolLSB && !PolRSB && !PolNP){
		if(doCtauUncer){
			int nEvents = 50;
			double promptFrac = 0, nonpromptFrac = 0, bkgFrac = 0;

			RooRealVar* bkgTauDSD = (RooRealVar*)ws->var("bkgTauDSD");
			RooRealVar* bkgTauFD = (RooRealVar*)ws->var("bkgTauFD");
			RooRealVar* bkgTauSSD_SBL = (RooRealVar*)ws->var("bkgTauSSD_SBL");
			RooRealVar* bkgTauSSD_SBR = (RooRealVar*)ws->var("bkgTauSSD_SBR");
			RooRealVar* fBkgDSD_SBL = (RooRealVar*)ws->var("fBkgDSD_SBL");
			RooRealVar* fBkgDSD_SBR = (RooRealVar*)ws->var("fBkgDSD_SBR");
			RooRealVar* fBkgSBL = (RooRealVar*)ws->var("fBkgSBL");
			RooRealVar* fBkgSBR = (RooRealVar*)ws->var("fBkgSBR");
			RooRealVar* fBkgSSDR_SBL = (RooRealVar*)ws->var("fBkgSSDR_SBL");
			RooRealVar* fBkgSSDR_SBR = (RooRealVar*)ws->var("fBkgSSDR_SBR");
			RooRealVar* fPrompt = (RooRealVar*)ws->var("fPrompt");
			RooRealVar* fBKG = (RooRealVar*)ws->var("fBkg");
			RooRealVar* fracGauss2 = (RooRealVar*)ws->var("fracGauss2");
			RooRealVar* nonPromptTau = (RooRealVar*)ws->var("nonPromptTau");
			RooArgSet *paraVars = new RooArgSet(*bkgTauDSD,*bkgTauFD,*bkgTauSSD_SBL,*bkgTauSSD_SBR,
					*fBkgDSD_SBL,*fBkgDSD_SBR, *fBkgSSDR_SBL,*fBkgSSDR_SBR,
					*nonPromptTau);
			paraVars->add(RooArgSet(*fPrompt,*fBKG,*fracGauss2,*fBkgSBL,*fBkgSBR));

			// create Hesse pdf and generate dataset
			RooAbsPdf *multiVarPdf = (RooAbsPdf*)result->createHessePdf(*paraVars);
			RooDataSet *multiVarData = (RooDataSet*)multiVarPdf->generate(*paraVars,nEvents);

			TH1D* histPFracDist = new TH1D();
			TH1D* histNPFracDist = new TH1D();
			TH1D* histBGFracDist = new TH1D();

			double IntePR=0., InteNP=0., InteBG=0.;
			double fracP1=0., fracNP1=0., fracBkg=0.;
			for(int n = 0; n < nEvents; n++) {
				if(n%40==0) std::cout << (double)n/nEvents*100 << "%" << std::endl;
				RooArgSet* args = (RooArgSet*)multiVarData->get(n);

				double bkgTauFD_ = ((RooRealVar*)args->find("bkgTauFD"))->getVal();
				ws->var("bkgTauFD")->setVal(bkgTauFD_);
				double bkgTauSSD_SBL_ = ((RooRealVar*)args->find("bkgTauSSD_SBL"))->getVal();
				ws->var("bkgTauSSD_SBL")->setVal(bkgTauSSD_SBL_);
				double bkgTauSSD_SBR_ = ((RooRealVar*)args->find("bkgTauSSD_SBR"))->getVal();
				ws->var("bkgTauSSD_SBR")->setVal(bkgTauSSD_SBR_);
				double fBkgDSD_SBL_ = ((RooRealVar*)args->find("fBkgDSD_SBL"))->getVal();
				ws->var("fBkgDSD_SBL")->setVal(fBkgDSD_SBL_);
				double fBkgDSD_SBR_ = ((RooRealVar*)args->find("fBkgDSD_SBR"))->getVal();
				ws->var("fBkgDSD_SBR")->setVal(fBkgDSD_SBR_);
				double fBkgSBL_ = ((RooRealVar*)args->find("fBkgSBL"))->getVal();
				ws->var("fBkgSBL")->setVal(fBkgSBL_);
				double fBkgSBR_ = ((RooRealVar*)args->find("fBkgSBR"))->getVal();
				ws->var("fBkgSBR")->setVal(fBkgSBR_);
				if(!(nState==4 && ptBin > 7)){
					double fBkgSSDR_SBL_ = ((RooRealVar*)args->find("fBkgSSDR_SBL"))->getVal();
					ws->var("fBkgSSDR_SBL")->setVal(fBkgSSDR_SBL_);
					double fBkgSSDR_SBR_ = ((RooRealVar*)args->find("fBkgSSDR_SBR"))->getVal();
					ws->var("fBkgSSDR_SBR")->setVal(fBkgSSDR_SBR_);
				}
				double fPrompt_ = ((RooRealVar*)args->find("fPrompt"))->getVal();
				ws->var("fPrompt")->setVal(fPrompt_);
				double fBkg_ = ((RooRealVar*)args->find("fBkg"))->getVal();
				ws->var("fBkg")->setVal(fBkg_);
				double nonPromptTau_ = ((RooRealVar*)args->find("nonPromptTau"))->getVal();
				ws->var("nonPromptTau")->setVal(nonPromptTau_);
				double bkgTauDSD_ = ((RooRealVar*)args->find("bkgTauDSD"))->getVal();
				ws->var("bkgTauDSD")->setVal(bkgTauDSD_);
				double fracGauss2_ = ((RooRealVar*)args->find("fracGauss2"))->getVal();
				ws->var("fracGauss2")->setVal(fracGauss2_);

				fracP1 = fPrompt_; 
				fracBkg = fBkg_; 
				fracNP1 =  1.-fracP1-fracBkg;

				// histogram method
				//vector<double> InteRltsTemp = calculateInte(ws,dataJpsictErr,ctCutMinPR,ctCutMaxPR);
				//IntePR    = InteRltsTemp[0];
				//InteNP    = InteRltsTemp[1];
				//InteBG    = InteRltsTemp[2];
				//
				//promptFrac    = IntePR * fracP1;
				//nonpromptFrac = InteNP * fracNP1;
				//bkgFrac       = InteBG * fracBkg;

				// integral method
				// prompt fraction in signal region
				i_fPR = PRpdf->createIntegral(RooArgSet(*ct), NormSet(RooArgSet(*ct)), Range("PR"));
				double fPR = i_fPR->getVal();
				// non prompt fraction in signal region
				i_fNPR = NPpdf->createIntegral(RooArgSet(*ct), NormSet(RooArgSet(*ct)), Range("PR"));
				double fNPR = i_fNPR->getVal();
				// background fraction in signal region
				i_fBkg = BGpdf->createIntegral(RooArgSet(*ct), NormSet(RooArgSet(*ct)), Range("PR"));
				double fBkgInSR = i_fBkg->getVal();

				promptFrac = fPR * fracP1;
				nonpromptFrac = fNPR * fracNP1;
				bkgFrac = fBkgInSR * fracBkg;

				histPFracDist->Fill(promptFrac);
				histNPFracDist->Fill(nonpromptFrac);
				histBGFracDist->Fill(bkgFrac);
			}

			fPerr  = histPFracDist->GetRMS() /(fNPinP*fNP1 + fPinP*fP1 + fBGinP*fBkg);
			fNPerr = histNPFracDist->GetRMS()/(fNPinP*fNP1 + fPinP*fP1 + fBGinP*fBkg);
			fBGerr = histBGFracDist->GetRMS()/(fNPinP*fNP1 + fPinP*fP1 + fBGinP*fBkg);

		} // doCtauUncer
	} // if(!MC && ...)

	// calculate error on total background fraction
	double fTBGerr = TMath::Sqrt(TMath::Power(fNPerr,2) + TMath::Power(fBGerr,2));
	std::cout << "error on total background fraction: " << fTBGerr << std::endl;

	output->cd();

	//------------------------------------------------------------------------------------------------
	// fill histogram with combinatorial background fraction
	std::stringstream bkgname;
	bkgname << ";;fraction of comb. BG in " << onia::nSigMass << "  sigma window, prompt region";
	TH1D* hFracBG = new TH1D("comb_background_fraction", bkgname.str().c_str(), 1, 0., 1.);
	hFracBG->SetBinContent(1, fBGsig);
	hFracBG->SetBinError(1, fBGerr);
	hFracBG->Write();

	// fill histogram with non prompt background fraction
	std::stringstream NPbkgname;
	NPbkgname << ";;fraction of non prompt BG in " << onia::nSigMass << "  sigma window, prompt region";
	TH1D* hFracNPBG = new TH1D("nonprompt_background_fraction", NPbkgname.str().c_str(), 1, 0., 1.);
	hFracNPBG->SetBinContent(1, fNPB);
	hFracNPBG->SetBinError(1, fNPerr);
	hFracNPBG->Write();

	// fill histogram with total background fraction
	std::stringstream tbkgname;
	tbkgname << ";;fraction of total BG in " << onia::nSigMass << "  sigma window, prompt region";
	TH1D* hFracTBG = new TH1D("background_fraction", tbkgname.str().c_str(), 1, 0., 1.);
	if(scaleFracBg){
		//fTBGsig = 0.001;
		char Filename[200];
		sprintf(Filename,"%s/results_Psi%dS_rap%d_pT%d_cpm%d.root",polDataPath,nState-3,rapBin,ptBin,cpmBin);
		TFile *resultFile = new TFile(Filename,"R");
		TH1D* SubtractedBG_test=(TH1D*)resultFile->Get("SubtractedBG_test");
		double ratio = SubtractedBG_test -> GetMean();
		std::cout << "fBG_sub / fBG : " << ratio << std::endl;
		fTBGsig = fTBGsig / ratio ;
		resultFile->Close();
	}
	output->cd();
	hFracTBG->SetBinContent(1, fTBGsig);
	//hFracTBG->SetBinContent(1, fTBGsig*2.);
	hFracTBG->SetBinError(1, fTBGerr);
	hFracTBG->Write();

	// fill histogram with prompt fraction
	std::stringstream Pname;
	Pname << ";;fraction of prompt events in " << onia::nSigMass << "  sigma window, prompt region (1-fNP-fBkg)";
	TH1D* hFracP = new TH1D("prompt_fraction", Pname.str().c_str(), 1, 0., 1.);
	hFracP->SetBinContent(1, fP);
	hFracP->SetBinError(1, fPerr);
	hFracP->Write();

	// fill histogram with fraction of LSB
	TH1D* hFracLSB = new TH1D("fraction_LSB", ";;f_{LSB}", 1, 0., 1.);
	hFracLSB->SetBinContent(1, fracLSB);
	hFracLSB->Write();

	// fill histogram with events in signal region
	std::stringstream evtSRname;
	evtSRname << ";;events in " << onia::nSigMass << "  sigma window";
	TH1D* hEvtSR = new TH1D("events_SR", evtSRname.str().c_str(), 3, 0., 3.);
	// fill histogram only for data
	if(!MC){
		hEvtSR->SetBinContent(1, nP);
		hEvtSR->SetBinContent(2, nNP);
		hEvtSR->SetBinContent(3, nBG);
		hEvtSR->Write();
	}

	// fill histogram with events in prompt signal region
	std::stringstream evtPSRname;
	evtPSRname << ";;events in " << onia::nSigMass << "  sigma window and low ctau region";
	TH1D* hEvtPSR = new TH1D("events_promptSR", evtPSRname.str().c_str(), 3, 0., 3.);
	// fill histogram only for data
	if(!MC){
		hEvtPSR->SetBinContent(1, fPinP*nP);
		hEvtPSR->SetBinContent(2, fNPinP*nNP);
		hEvtPSR->SetBinContent(3, fBGinP*nBG);
		hEvtPSR->Write();
	}

	// fill histogram with weighted sigma
	TH1D* hsigma = new TH1D("weighted_sigma", ";;weighted #sigma", 1, 0, 1);
	hsigma->SetBinContent(1, sigma);
	hsigma->Write();

	//------------------------------------------------------------------------------------------------
	// build the 3D (pT, |y|, M) histos for the L and R mass sideband
	TH3D* hBG_pTRapMass_L = new TH3D("hBG_pTRapMass_L", ";p_{T} [GeV/c]; |y|; M [GeV]",
			7, onia::pTRange[rapBin-1][ptBin-1], onia::pTRange[rapBin-1][ptBin],
			2, onia::rapForPTRange[rapBin-1], onia::rapForPTRange[rapBin],
			7, massMinSR, massMaxSR); // signal mass window!
	hBG_pTRapMass_L->Sumw2();

	TH3D* hBG_pTRapMass_R = new TH3D("hBG_pTRapMass_R", ";p_{T} [GeV/c]; |y|; M [GeV]",
			7, onia::pTRange[rapBin-1][ptBin-1], onia::pTRange[rapBin-1][ptBin],
			2, onia::rapForPTRange[rapBin-1], onia::rapForPTRange[rapBin],
			7, massMinSR, massMaxSR); // signal mass window!
	hBG_pTRapMass_R->Sumw2();

	TH3D* hBG_pTRapMass_highct_L = new TH3D("hBG_pTRapMass_highct_L", ";p_{T} [GeV/c]; |y|; M [GeV]",
			7, onia::pTRange[rapBin-1][ptBin-1], onia::pTRange[rapBin-1][ptBin],
			2, onia::rapForPTRange[rapBin-1], onia::rapForPTRange[rapBin],
			7, massMinSR, massMaxSR);
	hBG_pTRapMass_highct_L->Sumw2();

	TH3D* hBG_pTRapMass_highct_R = new TH3D("hBG_pTRapMass_highct_R", ";p_{T} [GeV/c]; |y|; M [GeV]",
			7, onia::pTRange[rapBin-1][ptBin-1], onia::pTRange[rapBin-1][ptBin],
			2, onia::rapForPTRange[rapBin-1], onia::rapForPTRange[rapBin],
			7, massMinSR, massMaxSR);
	hBG_pTRapMass_highct_R->Sumw2();

	TH3D* hNP_pTRapMass = new TH3D("hNP_pTRapMass_NP", ";p_{T} [GeV/c]; |y|; M [GeV]",
			7, onia::pTRange[rapBin-1][ptBin-1], onia::pTRange[rapBin-1][ptBin],
			2, onia::rapForPTRange[rapBin-1], onia::rapForPTRange[rapBin],
			7, massMinSR, massMaxSR);
	hNP_pTRapMass->Sumw2();

	TH3D* hSR_pTRapMass = new TH3D("hSR_pTRapMass", ";p_{T} [GeV/c]; |y|; M [GeV]",
			7, onia::pTRange[rapBin-1][ptBin-1], onia::pTRange[rapBin-1][ptBin],
			2, onia::rapForPTRange[rapBin-1], onia::rapForPTRange[rapBin],
			7, massMinSR, massMaxSR);
	hSR_pTRapMass->Sumw2();

	//------------------------------------------------------------------------------------------------
	// loop through tree, fill background histos and save data in pt and y bins
	int index = -1;
	int n = intree->GetEntries();
	int MCevents = 0;
	int count=0;

	for(int i = 0; i < n; i++){

		long iEntry = intree->LoadTree(i);
		intree->GetEntry(iEntry);
		if(i % 100000 == 0) {std::cout << "entry " << i << " out of " << n << std::endl;}

		// ------------------------- TLorentzVecotrs -------------------------
		if(jpsi->Pt() >= onia::pTRange[rapBin-1][ptBin-1] &&
				jpsi->Pt() < onia::pTRange[rapBin-1][ptBin] &&
				TMath::Abs(jpsi->Rapidity()) >= onia::rapForPTRange[rapBin-1] &&
				TMath::Abs(jpsi->Rapidity()) < onia::rapForPTRange[rapBin]){

			//store TLorentzVectors of the two muons in the given pT and rap cell
			// left sideband
			if(PolLSB){
				if(jpsi->M() >= massMinL && jpsi->M() < massMaxL && TMath::Abs(jpsict) < ctauCut){
					outtree->Fill();
					pT_PSR->Fill(jpsi->Pt());
					cpm_PSR->Fill(cpmval);
					rap_PSR->Fill(TMath::Abs(jpsi->Rapidity()));
				}
			} // PolLSB
			// right sideband
			else if(PolRSB){
				if(jpsi->M() >= massMinR && jpsi->M() <= massMaxR && TMath::Abs(jpsict) < ctauCut){
					outtree->Fill();
					pT_PSR->Fill(jpsi->Pt());
					cpm_PSR->Fill(cpmval);
					rap_PSR->Fill(TMath::Abs(jpsi->Rapidity()));
				}
			} // PolRSB
			// for non prompt data
			else if(PolNP){
				if(jpsi->M() >= massMinSR && jpsi->M() < massMaxSR && TMath::Abs(jpsict) > ctauCut){
					outtree->Fill();
					pT_PSR->Fill(jpsi->Pt());
					cpm_PSR->Fill(cpmval);
					rap_PSR->Fill(TMath::Abs(jpsi->Rapidity()));
				}
			} // PolNP
			// for prompt data and MC
			//store only events from signal region
			else if(MC){
				if(jpsi->M() >= massMinSR && jpsi->M() < massMaxSR){
					outtree->Fill();
					pT_PSR->Fill(jpsi->Pt());
					cpm_PSR->Fill(cpmval);
					rap_PSR->Fill(TMath::Abs(jpsi->Rapidity()));
					if(MC) MCevents++;
				}
			}
			else{
				if(jpsi->M() >= massMinSR && jpsi->M() < massMaxSR && TMath::Abs(jpsict) < ctauCut){
					count++;
					outtree->Fill();
					pT_PSR->Fill(jpsi->Pt());
					cpm_PSR->Fill(cpmval);
					rap_PSR->Fill(TMath::Abs(jpsi->Rapidity()));
				}
			} // else

			//---------------------------- mass histograms for background model ------------------------------
			// for MC: fill mass histograms with random mass from signal region
			// fill rapidity and pT histograms with all events
			if(MC){
				pT_L->Fill(jpsi->Pt());
				pT_R->Fill(jpsi->Pt());
				pT_highct_L->Fill(jpsi->Pt());
				pT_highct_R->Fill(jpsi->Pt());
				pT_NP->Fill(jpsi->Pt());
				rap_L->Fill(TMath::Abs(jpsi->Rapidity()));
				rap_R->Fill(TMath::Abs(jpsi->Rapidity()));
				rap_highct_L->Fill(TMath::Abs(jpsi->Rapidity()));
				rap_highct_R->Fill(TMath::Abs(jpsi->Rapidity()));
				rap_NP->Fill(TMath::Abs(jpsi->Rapidity()));

				hNP_pTRapMass->Fill(jpsi->Pt(), TMath::Abs(jpsi->Rapidity()), gRandom->Uniform(massMinSR, massMaxSR));
				hSR_pTRapMass->Fill(jpsi->Pt(), TMath::Abs(jpsi->Rapidity()), gRandom->Uniform(massMinSR, massMaxSR));
				// split events up to fill in left and right background histograms
				// not all events are filled to be able to put together background histogram
				if(i%3==0){
					hBG_pTRapMass_L->Fill(jpsi->Pt(), TMath::Abs(jpsi->Rapidity()), gRandom->Uniform(massMinSR, massMaxSR));
					hBG_pTRapMass_highct_L->Fill(jpsi->Pt(), TMath::Abs(jpsi->Rapidity()), gRandom->Uniform(massMinSR, massMaxSR));
				}
				else if(i%5==0){
					hBG_pTRapMass_R->Fill(jpsi->Pt(), TMath::Abs(jpsi->Rapidity()), gRandom->Uniform(massMinSR, massMaxSR));
					hBG_pTRapMass_highct_R->Fill(jpsi->Pt(), TMath::Abs(jpsi->Rapidity()), gRandom->Uniform(massMinSR, massMaxSR));
				}
			} // if (MC)

			else{
				// store cosTheta and phi distributions of the background
				// events gets index 0 if it is in the left sideband and 1 if it is in the right one
				// events with index 2 are from the high ctau non prompt region
				// events with index 3 and 4 are from the high ctau region from left and right sideband
				// events with index 5 are from the prompt signal region
				if(jpsi->M() >= massMinL && jpsi->M() < massMaxL && TMath::Abs(jpsict) < ctauCut){
					index = 0;
					hBG_pTRapMass_L->Fill(jpsi->Pt(), TMath::Abs(jpsi->Rapidity()), funcBG->GetRandom(massMinSR, massMaxSR));
					pT_L->Fill(jpsi->Pt());
					cpm_L->Fill(cpmval);
					rap_L->Fill(TMath::Abs(jpsi->Rapidity()));
				}
				else if(jpsi->M() >= massMinR && jpsi->M() < massMaxR && TMath::Abs(jpsict) < ctauCut){
					index = 1;
					hBG_pTRapMass_R->Fill(jpsi->Pt(), TMath::Abs(jpsi->Rapidity()), funcBG->GetRandom(massMinSR, massMaxSR));
					pT_R->Fill(jpsi->Pt());
					cpm_R->Fill(cpmval);
					rap_R->Fill(TMath::Abs(jpsi->Rapidity()));
				}
				else if(jpsi->M() >= massMinSR && jpsi->M() < massMaxSR && jpsict >= ctauCut){
					index = 2;
					hNP_pTRapMass->Fill(jpsi->Pt(), TMath::Abs(jpsi->Rapidity()), funcSig->GetRandom(massMinSR, massMaxSR));
					pT_NP->Fill(jpsi->Pt());
					cpm_NP->Fill(cpmval);
					rap_NP->Fill(TMath::Abs(jpsi->Rapidity()));
				}
				else if(jpsi->M() >= massMinL && jpsi->M() < massMaxL && jpsict >= ctauCut){
					index = 3;
					if(PolNP)
						hBG_pTRapMass_highct_L->Fill(jpsi->Pt(), TMath::Abs(jpsi->Rapidity()), funcBG->GetRandom(massMinSR, massMaxSR));
					else
						hBG_pTRapMass_highct_L->Fill(jpsi->Pt(), TMath::Abs(jpsi->Rapidity()), funcSig->GetRandom(massMinSR, massMaxSR));
					pT_highct_L->Fill(jpsi->Pt());
					cpm_highct_L->Fill(cpmval);
					rap_highct_L->Fill(TMath::Abs(jpsi->Rapidity()));
				}
				else if(jpsi->M() >= massMinR && jpsi->M() < massMaxR && jpsict >= ctauCut){
					index = 4;
					if(PolNP)
						hBG_pTRapMass_highct_R->Fill(jpsi->Pt(), TMath::Abs(jpsi->Rapidity()), funcBG->GetRandom(massMinSR, massMaxSR));
					else
						hBG_pTRapMass_highct_R->Fill(jpsi->Pt(), TMath::Abs(jpsi->Rapidity()), funcSig->GetRandom(massMinSR, massMaxSR));
					pT_highct_R->Fill(jpsi->Pt());
					cpm_highct_R->Fill(cpmval);
					rap_highct_R->Fill(TMath::Abs(jpsi->Rapidity()));
				}
				else if(jpsi->M() >= massMinSR && jpsi->M() < massMaxSR && jpsict < ctauCut){
					index = 5;
					hSR_pTRapMass->Fill(jpsi->Pt(), TMath::Abs(jpsi->Rapidity()), funcBG->GetRandom(massMinSR, massMaxSR));
				}

				else continue;
			}// else (filling data histograms)

			///////////////////////
			calcPol(*lepP, *lepN);
			///////////////////////

			for(int iFrame = 0; iFrame < onia::kNbFrames; iFrame++){

				// folding in phi
				double phiFolded = thisPhi[iFrame];
				double thetaAdjusted = thisCosTh[iFrame];
				if(thisPhi[iFrame] > -90. && thisPhi[iFrame] < 0.) phiFolded *= -1;
				else if(thisPhi[iFrame] > 90 && thisPhi[iFrame] < 180){
					phiFolded = 180. - thisPhi[iFrame];
					thetaAdjusted *= -1;
				}
				else if(thisPhi[iFrame] > -180. && thisPhi[iFrame] < -90.){
					phiFolded = 180. + thisPhi[iFrame];
					thetaAdjusted *= -1;
				}

				// if bool folding is true, folding is applied to all background histograms
				if(folding){
					thisPhi[iFrame] = phiFolded;
					thisCosTh[iFrame] = thetaAdjusted;
				}

				// filling histograms
				if(MC){
					hNPBG_cosThetaPhi[iFrame]->Fill(thisCosTh[iFrame], thisPhi[iFrame]);
					hSR_cosThetaPhi[iFrame]->Fill(thisCosTh[iFrame], thisPhi[iFrame]);
					if(i%3==0){
						hBG_cosThetaPhiL[iFrame]->Fill(thisCosTh[iFrame], thisPhi[iFrame]);
						hBGinNP_cosThetaPhiL[iFrame]->Fill(thisCosTh[iFrame], thisPhi[iFrame]);
					}
					else if(i%5==0){
						hBG_cosThetaPhiR[iFrame]->Fill(thisCosTh[iFrame], thisPhi[iFrame]);
						hBGinNP_cosThetaPhiR[iFrame]->Fill(thisCosTh[iFrame], thisPhi[iFrame]);
					}
				}
				// for data: fill histograms according to the different regions
				else{
					if(index == 0) hBG_cosThetaPhiL[iFrame]->Fill(thisCosTh[iFrame], thisPhi[iFrame]);
					else if(index == 1) hBG_cosThetaPhiR[iFrame]->Fill(thisCosTh[iFrame], thisPhi[iFrame]);
					else if(index == 2) hNPBG_cosThetaPhi[iFrame]->Fill(thisCosTh[iFrame], thisPhi[iFrame]);
					else if(index == 3) hBGinNP_cosThetaPhiL[iFrame]->Fill(thisCosTh[iFrame], thisPhi[iFrame]);
					else if(index == 4) hBGinNP_cosThetaPhiR[iFrame]->Fill(thisCosTh[iFrame], thisPhi[iFrame]);
					else if(index == 5) hSR_cosThetaPhi[iFrame]->Fill(thisCosTh[iFrame], thisPhi[iFrame]);
				}

			} // iFrame

		} // if(onia...)
	} // i

	std::cout << "total events in PRSR: " << count << std::endl;

	// fill histograms with number of events for MC
	if(MC){
		// number of events in signal region
		hEvtSR->SetBinContent(1, fP*MCevents);
		// number of events in prompt signal region (same as in signal region)
		hEvtPSR->SetBinContent(1, fP*MCevents);
		hEvtSR->Write();
		hEvtPSR->Write();
	}

	//---------------- binning algorithm
	int nBinsPhi = 16, nBinsCosth = 160;
	int totalBins = 0, filledBins = 0;
	for(int binCosth = 0; binCosth < hBGinNP_cosThetaPhiR[2]->GetNbinsX(); binCosth++){
		for(int binPhi = 0; binPhi < hBGinNP_cosThetaPhiR[2]->GetNbinsY(); binPhi++){
			totalBins++;
			//use NP histo (same physical coverage, but better filled, no holes -> better estimate of coverage)
			int binContent = hNPBG_cosThetaPhi[2]->GetBinContent(binCosth+1,binPhi+1);
			if(binContent>0) filledBins++;
		}
	}

	double coverage = 2*(double)filledBins/(double)totalBins;
	nBinsCosth = 16*2/coverage;
	// find 2^n closest to nBinsCosth, but above the actual number
	nBinsCosth = findEvenNum((double)nBinsCosth);
	// set maximum binning to 64
	if(nBinsCosth > 64) nBinsCosth = 64;

	std::cout << "------------------------------------------------" << "\n"
		<< "Starting binning algorithm" << "\n"
		<< "filled bins: " << filledBins << "\n"
		<< "total bins: " << totalBins << "\n"
		<< "bin coverage: " << coverage << "\n"
		<< "starting point for binning in phi: " << nBinsPhi << "\n"
		<< "starting point for binning in cosTheta: " << nBinsCosth << "\n"
		<< "------------------------------------------------" <<std::endl;

	// calculate the integral of the lowstatBG histo (calculate all integrals of the 4 BG regions, and use the one with the smallest integral)
	int IntBG = hBGinNP_cosThetaPhiR[2]->Integral();
	if(IntBG > hBGinNP_cosThetaPhiL[2]->Integral()){
		IntBG = hBGinNP_cosThetaPhiL[2]->Integral();
		std::cout << "left high ct integral is smaller" << std::endl;
	}
	if(IntBG > hBG_cosThetaPhiL[2]->Integral()){
		IntBG = hBG_cosThetaPhiL[2]->Integral();
		std::cout << "left low ct integral is smaller" << std::endl;
	}
	if(IntBG > hBG_cosThetaPhiR[2]->Integral()){
		IntBG = hBG_cosThetaPhiR[2]->Integral();
		std::cout << "right low ct integral is smaller" << std::endl;
	}

	int IntNPBG = hNPBG_cosThetaPhi[2]->Integral();

	int nBinsPhiBG = nBinsPhi,
			nBinsCosthBG = nBinsCosth,
			nBinsPhiNPBG = nBinsPhi,
			nBinsCosthNPBG = nBinsCosth;

	// calculate average events per-bin cell for background histo
	double Naverage = (double)IntBG/((double)nBinsPhi*nBinsCosth*coverage/2.);
	std::cout << "average cell coverage: " << Naverage << std::endl;

	// if average events per bin is bigger than 10, no rebinning is needed
	if(Naverage > 10){
		std::cout << "Rebinning is not necessary in this case." << "\n"
			<< "Ending binning algorithm." << "\n"
			<< "------------------------------------------------" << std::endl;
	}
	// otherwise rebin
	else{
		std::cout << "------------------------------------------------" << "\n"
			<< "old cosTheta binning: " << nBinsCosth << "\n"
			<< "old phi binning: " << nBinsPhi << std::endl;

		//set nBinsPhi to the lowest 2^n, such that nBinsPhi > nBinsCosth*coverage/2
		nBinsPhi = findEvenNum(nBinsCosth*coverage/2.);

		std::cout << "closest 2^n number to cosTheta bins: " << nBinsCosth << "\n"
			<< "lowest 2^n number so that phi bins > cosTheta bins * coverage/2: " << nBinsPhi << "\n"
			<< "------------------------------------------------" << std::endl;

		// set minimum binning
		int nBinsPhiMin = 8,
				nBinsCosthMin = 8;
		if(folding) nBinsPhiMin = 16;

		//BG
		nBinsCosthBG = nBinsCosth;
		nBinsPhiBG = nBinsPhi;
		double NaverageBG = 0.;

		for(int i = 0; i < 500; i++){

			std::cout << "looping for correct binning in background histogram" << std::endl;
			// If the mimimum number of bins for both phi and costh are reached, stop the loop
			if(nBinsPhiBG/2 < nBinsPhiMin && nBinsCosthBG/2 < nBinsCosthMin) break;

			//Change the binning, first in phi, then in costh:
			if(nBinsPhiBG/2 >= nBinsPhiMin) nBinsPhiBG = nBinsPhiBG/2;  //This ensures a mimimum number of bins in phi, e.g. 4
			NaverageBG = (double)IntBG/((double)nBinsPhiBG*nBinsCosthBG*coverage/2.);
			std::cout << "average bin content per cell after " << i << " phi rebinning: " << NaverageBG << std::endl;
			if(NaverageBG > 10) break;

			if(nBinsCosthBG/2 >= nBinsCosthMin) nBinsCosthBG = nBinsCosthBG/2; //This ensures a mimimum number of bins in costh, e.g. 4
			NaverageBG = (double)IntBG/((double)nBinsPhiBG*nBinsCosthBG*coverage/2.);
			std::cout << "average bin content per cell after " << i << " cosTheta rebinning: " << NaverageBG << std::endl;
			if(NaverageBG > 10) break;
		}
		std::cout << "average bin content per cell exceeds 10: " << NaverageBG << "\n"
			<< "phi bins = " << nBinsPhiBG << ", cosTheta bins = " << nBinsCosthBG << "\n"
			<< "------------------------------------------------" << std::endl;

		//NPBG
		nBinsCosthNPBG=nBinsCosth;
		nBinsPhiNPBG=nBinsPhi;
		double NaverageNPBG=0.;

		for(int i = 0; i < 500; i++){

			std::cout << "looping for correct binning in non prompt histogram" << std::endl;
			// If the mimimum number of bins for both phi and costh are reached, stop the loop
			if(nBinsPhiNPBG/2 < nBinsPhiMin && nBinsCosthNPBG/2 < nBinsCosthMin) break;

			//Change the binning, first in phi, then in costh:
			if(nBinsPhiNPBG/2 >= nBinsPhiMin) nBinsPhiNPBG = nBinsPhiNPBG/2;  //This ensures a mimimum number of bins in phi, e.g. 4
			NaverageNPBG = (double)IntNPBG/((double)nBinsPhiNPBG*nBinsCosthNPBG*coverage/2.);
			std::cout << "average bin content per cell after " << i << " phi rebinning: " << NaverageNPBG << std::endl;
			if(NaverageNPBG > 10) break;

			if(nBinsCosthNPBG/2 >= nBinsCosthMin) nBinsCosthNPBG = nBinsCosthNPBG/2; //This ensures a mimimum number of bins in costh, e.g. 4
			NaverageNPBG = (double)IntNPBG/((double)nBinsPhiNPBG*nBinsCosthNPBG*coverage/2.);
			std::cout << "average bin content per cell after " << i << " cosTheta rebinning: " << NaverageNPBG << std::endl;
			if(NaverageNPBG > 10) break;
		}
		std::cout << "average bin content per cell exceeds 10: " << NaverageNPBG << "\n"
			<< "phi bins = " << nBinsPhiNPBG << ", cosTheta bins = " << nBinsCosthNPBG << "\n"
			<< "------------------------------------------------" << std::endl;

		// when forceBinning, set binning of Psi1S consistently to non prompt binning and Psi2S consistently to background binning
		if(forceBinning){

			if(nState==4){
				std::cout << "Force consistent binning equal to binning of non prompt histogram" << std::endl;
				nBinsPhiBG = nBinsPhiNPBG;
				nBinsCosthBG = nBinsCosthNPBG;
			}
			else if(nState==5){
				std::cout << "Force consistent binning equal to binning of background histogram" << std::endl;
				nBinsPhiNPBG = nBinsPhiBG;
				nBinsCosthNPBG = nBinsCosthBG;
			}

		} // if(forceBinning)

	} // else Naverage < 10

	std::cout << "final binning for background histogram: " << "\n"
		<< "phi bins: " << nBinsPhiBG << "\n"
		<< "cosTheta bins: " << nBinsCosthBG << "\n"
		<< "final binning for non prompt histogram" << "\n"
		<< "phi bins: " << nBinsPhiNPBG << "\n"
		<< "cosTheta bins: " << nBinsCosthNPBG << "\n"
		<< "------------------------------------------------" << std::endl;

	//loop again with new binning
	for(int iFrame = 0; iFrame < onia::kNbFrames; iFrame++){
		//book the 2D (cosTheta, phi) histos for the L and R mass sideband
		std::stringstream nameL, nameR, nameNP, nameBGinNPL, nameBGinNPR, nameSR, title;
		nameL << "hBG_cosThetaPhi_" << onia::frameLabel[iFrame] << "_L";
		nameR << "hBG_cosThetaPhi_" << onia::frameLabel[iFrame] << "_R";
		nameNP << "hNPBG_cosThetaPhi_" << onia::frameLabel[iFrame];
		nameBGinNPL << "hBGinNP_cosThetaPhi_" << onia::frameLabel[iFrame] << "_L";
		nameBGinNPR << "hBGinNP_cosThetaPhi_" << onia::frameLabel[iFrame] << "_R";
		nameSR << "hSR_cosThetaPhi_" << onia::frameLabel[iFrame];
		title << ";cos#theta_{"<< onia::frameLabel[iFrame] << "};#phi_{" << onia::frameLabel[iFrame] << "} [deg]";

		delete hBG_cosThetaPhiL[iFrame];
		hBG_cosThetaPhiL[iFrame] = new TH2D(nameL.str().c_str(), title.str().c_str(),
				nBinsCosthBG, onia::cosTMin, onia::cosTMax, nBinsPhiBG, onia::phiPolMin, onia::phiPolMax);
		hBG_cosThetaPhiL[iFrame]->Sumw2();
		delete hBG_cosThetaPhiR[iFrame];
		hBG_cosThetaPhiR[iFrame] = new TH2D(nameR.str().c_str(), title.str().c_str(),
				nBinsCosthBG, onia::cosTMin, onia::cosTMax, nBinsPhiBG, onia::phiPolMin, onia::phiPolMax);
		hBG_cosThetaPhiR[iFrame]->Sumw2();

		delete hNPBG_cosThetaPhi[iFrame];
		hNPBG_cosThetaPhi[iFrame] = new TH2D(nameNP.str().c_str(), title.str().c_str(),
				nBinsCosthNPBG, onia::cosTMin, onia::cosTMax, nBinsPhiNPBG, onia::phiPolMin, onia::phiPolMax);
		hNPBG_cosThetaPhi[iFrame]->Sumw2();

		delete hBGinNP_cosThetaPhiL[iFrame];
		hBGinNP_cosThetaPhiL[iFrame] = new TH2D(nameBGinNPL.str().c_str(), title.str().c_str(),
				nBinsCosthBG, onia::cosTMin, onia::cosTMax, nBinsPhiBG, onia::phiPolMin, onia::phiPolMax);
		hBGinNP_cosThetaPhiL[iFrame]->Sumw2();
		delete hBGinNP_cosThetaPhiR[iFrame];
		hBGinNP_cosThetaPhiR[iFrame] = new TH2D(nameBGinNPR.str().c_str(), title.str().c_str(),
				nBinsCosthBG, onia::cosTMin, onia::cosTMax, nBinsPhiBG, onia::phiPolMin, onia::phiPolMax);
		hBGinNP_cosThetaPhiR[iFrame]->Sumw2();

		delete hSR_cosThetaPhi[iFrame];
		hSR_cosThetaPhi[iFrame] = new TH2D(nameSR.str().c_str(), title.str().c_str(),
				nBinsCosthBG, onia::cosTMin, onia::cosTMax, nBinsPhiBG, onia::phiPolMin, onia::phiPolMax);
		hSR_cosThetaPhi[iFrame]->Sumw2();

	} // iFrame

	for(int i = 0; i < n; i++){

		long iEntry = intree->LoadTree(i);
		intree->GetEntry(iEntry);
		if(i % 100000 == 0) std::cout << "entry " << i << " out of " << n << std::endl;

		if(jpsi->Pt() >= onia::pTRange[rapBin-1][ptBin-1] &&
				jpsi->Pt() < onia::pTRange[rapBin-1][ptBin] &&
				TMath::Abs(jpsi->Rapidity()) >= onia::rapForPTRange[rapBin-1] &&
				TMath::Abs(jpsi->Rapidity()) < onia::rapForPTRange[rapBin]){

			if(!MC){
				if(jpsi->M() >= massMinL && jpsi->M() < massMaxL && TMath::Abs(jpsict) < ctauCut)
					index = 0;
				else if(jpsi->M() >= massMinR && jpsi->M() <= massMaxR && TMath::Abs(jpsict) < ctauCut)
					index = 1;
				else if(jpsi->M() >= massMinSR && jpsi->M() < massMaxSR && jpsict >= ctauCut)
					index = 2;
				else if(jpsi->M() >= massMinL && jpsi->M() < massMaxL && jpsict >= ctauCut)
					index = 3;
				else if(jpsi->M() >= massMinR && jpsi->M() <= massMaxR && jpsict >= ctauCut)
					index = 4;
				else if(jpsi->M() >= massMinSR && jpsi->M() < massMaxSR && jpsict < ctauCut)
					index = 5;
				else continue;
			}

			////////////////////////
			calcPol(*lepP, *lepN);
			////////////////////////

			for(int iFrame = 0; iFrame < onia::kNbFrames; iFrame++){

				// folding in phi
				double phiFolded = thisPhi[iFrame];
				double thetaAdjusted = thisCosTh[iFrame];
				if(thisPhi[iFrame] > -90. && thisPhi[iFrame] < 0.) phiFolded *= -1;
				else if(thisPhi[iFrame] > 90 && thisPhi[iFrame] < 180){
					phiFolded = 180. - thisPhi[iFrame];
					thetaAdjusted *= -1;
				}
				else if(thisPhi[iFrame] > -180. && thisPhi[iFrame] < -90.){
					phiFolded = 180. + thisPhi[iFrame];
					thetaAdjusted *= -1;
				}

				// if folding is true, apply folding in phi
				if(folding){
					thisPhi[iFrame] = phiFolded;
					thisCosTh[iFrame] = thetaAdjusted;
				}

				// filling histograms
				if(MC){
					hNPBG_cosThetaPhi[iFrame]->Fill(thisCosTh[iFrame], thisPhi[iFrame]);
					hSR_cosThetaPhi[iFrame]->Fill(thisCosTh[iFrame], thisPhi[iFrame]);
					if(i%3==0){
						hBG_cosThetaPhiL[iFrame]->Fill(thisCosTh[iFrame], thisPhi[iFrame]);
						hBGinNP_cosThetaPhiL[iFrame]->Fill(thisCosTh[iFrame], thisPhi[iFrame]);
					}
					if(i%5==0){
						hBG_cosThetaPhiR[iFrame]->Fill(thisCosTh[iFrame], thisPhi[iFrame]);
						hBGinNP_cosThetaPhiR[iFrame]->Fill(thisCosTh[iFrame], thisPhi[iFrame]);
					}
				}
				else{
					if(index == 0) hBG_cosThetaPhiL[iFrame]->Fill(thisCosTh[iFrame], thisPhi[iFrame]);
					if(index == 1) hBG_cosThetaPhiR[iFrame]->Fill(thisCosTh[iFrame], thisPhi[iFrame]);
					if(index == 2) hNPBG_cosThetaPhi[iFrame]->Fill(thisCosTh[iFrame], thisPhi[iFrame]);
					if(index == 3) hBGinNP_cosThetaPhiL[iFrame]->Fill(thisCosTh[iFrame], thisPhi[iFrame]);
					if(index == 4) hBGinNP_cosThetaPhiR[iFrame]->Fill(thisCosTh[iFrame], thisPhi[iFrame]);
					if(index == 5) hSR_cosThetaPhi[iFrame]->Fill(thisCosTh[iFrame], thisPhi[iFrame]);
				}
			} // iFrame
		} // if(onia...)
	} // i
	//------loop finished

	//---------------- end- binning algorithm

	//----------------------------------------------------------------------------------------------------
	// write background histos to file
	// 3D (pT, |y|, M) histos
	hBG_pTRapMass_L->Write();
	hBG_pTRapMass_R->Write();
	hSR_pTRapMass->Write();

	// subtract signal contamination from left sideband
	std::string namepTrapMassSRL = "SRinLSB_pTrapMass";
	TH3D* hSR_pTRapMass_L = (TH3D*) hSR_pTRapMass->Clone(namepTrapMassSRL.c_str());
	hSR_pTRapMass_L->Scale(fSRinPLSB/(1.*hSR_pTRapMass_L->Integral()));
	hBG_pTRapMass_L->Scale(1./(1.*hBG_pTRapMass_L->Integral()));
	hBG_pTRapMass_L = subtract3D(hBG_pTRapMass_L, hSR_pTRapMass_L);

	// subtract signal contamination from right sideband
	std::string namepTrapMassSRR = "SRinRSB_pTrapMass";
	TH3D* hSR_pTRapMass_R = (TH3D*) hSR_pTRapMass->Clone(namepTrapMassSRR.c_str());
	hSR_pTRapMass_R->Scale(fSRinPRSB/(1.*hSR_pTRapMass_R->Integral()));
	hBG_pTRapMass_R->Scale(1./(1.*hBG_pTRapMass_R->Integral()));
	hBG_pTRapMass_R = subtract3D(hBG_pTRapMass_R, hSR_pTRapMass_R);

	// add left and right sideband of low ctau region to combinatorial background histogram
	hBG_pTRapMass_L->Scale(fracLSB/(1.*hBG_pTRapMass_L->Integral()));
	hBG_pTRapMass_R->Scale((1.-fracLSB)/(1.*hBG_pTRapMass_R->Integral()));
	std::string namepTrapMasslowct = "comb_background_pTrapMass";
	TH3D* hBG_pTRapMass_lowct = (TH3D*) hBG_pTRapMass_L->Clone(namepTrapMasslowct.c_str());
	hBG_pTRapMass_lowct->Add(hBG_pTRapMass_R);
	hBG_pTRapMass_lowct->Write();

	// add left and right sideband of high ctau region to combinatorial background histogram in high ctau region
	hBG_pTRapMass_highct_L->Write();
	hBG_pTRapMass_highct_R->Write();
	hBG_pTRapMass_highct_L->Scale(fracLSB/(1.*hBG_pTRapMass_highct_L->Integral()));
	hBG_pTRapMass_highct_R->Scale((1.-fracLSB)/(1.*hBG_pTRapMass_highct_R->Integral()));
	std::string namepTrapMasshighct = "comb_background_highct_pTrapMass";
	TH3D* hBG_pTRapMass_highct = (TH3D*) hBG_pTRapMass_highct_L->Clone(namepTrapMasshighct.c_str());
	hBG_pTRapMass_highct->Add(hBG_pTRapMass_highct_R);
	hBG_pTRapMass_highct->Write();

	// create non prompt background histogram in signal region: (hNP)_norm - fBGinNP * (hBG_highct)_norm
	hNP_pTRapMass->Scale(1./(1.*hNP_pTRapMass->Integral()));
	hNP_pTRapMass->Write();
	hBG_pTRapMass_highct->Scale(fBGinNP/(1.*hBG_pTRapMass_highct->Integral()));
	std::string namepTrapMassNPS = "NPS_highct_pTrapMass";
	TH3D* hNPS_pTRapMass = (TH3D*) hNP_pTRapMass->Clone(namepTrapMassNPS.c_str());
	hNPS_pTRapMass = subtract3D(hNPS_pTRapMass, hBG_pTRapMass_highct);
	hNPS_pTRapMass->Write();

	// create total background
	std::string namepTrapMass = "background_pTrapMass";
	TH3D* hBG_pTRapMass = new TH3D();
	// for non prompt polarization: total background = high ct background
	if(PolLSB || PolRSB){
		hSR_pTRapMass->Scale(1./(1.*hSR_pTRapMass->Integral()));
		hBG_pTRapMass = (TH3D*) hSR_pTRapMass->Clone(namepTrapMass.c_str());
	}
	else if(PolNP)
		hBG_pTRapMass = (TH3D*) hBG_pTRapMass_highct->Clone(namepTrapMass.c_str());
	// add low ct background and non prompt background
	else{
		hNPS_pTRapMass->Scale(fNPB/(1.*hNPS_pTRapMass->Integral()));
		hBG_pTRapMass = (TH3D*) hNPS_pTRapMass->Clone(namepTrapMass.c_str());
		hBG_pTRapMass_lowct->Scale(fBGsig/(1.*hBG_pTRapMass_lowct->Integral()));
		hBG_pTRapMass->Add(hBG_pTRapMass_lowct);
	}
	hBG_pTRapMass->Write();

	double meancpm = 0;
	// mean cpm histos
	cpm_PSR->Scale(1./(1.*cpm_PSR->Integral()));
	TH1D* cpm_SRL = (TH1D*) cpm_PSR->Clone();
	cpm_SRL->Scale(fSRinPLSB/(1.*cpm_SRL->Integral()));
	cpm_L->Scale(1./(1.*cpm_L->Integral()));
	cpm_L->Add(cpm_SRL, -1.);
	if(PolLSB) meancpm = cpm_L->GetMean();

	TH1D* cpm_SRR = (TH1D*) cpm_PSR->Clone();
	cpm_SRR->Scale(fSRinPRSB/(1.*cpm_SRR->Integral()));
	cpm_R->Scale(1./(1.*cpm_R->Integral()));
	cpm_R->Add(cpm_SRR, -1.);
	if(PolRSB) meancpm = cpm_R->GetMean();

	cpm_L->Scale(fracLSB/(1.*cpm_L->Integral()));
	cpm_R->Scale((1.-fracLSB)/(1.*cpm_R->Integral()));
	cpm_L->Add(cpm_R);

	cpm_highct_L->Scale(fracLSB/(1.*cpm_highct_L->Integral()));
	cpm_highct_R->Scale((1.-fracLSB)/(1.*cpm_highct_R->Integral()));
	cpm_highct_L->Add(cpm_highct_R);

	cpm_NP->Scale(1./(1.*cpm_NP->Integral()));
	cpm_highct_L->Scale(fBGinNP/(1.*cpm_highct_L->Integral()));
	cpm_NP->Add(cpm_highct_L, -1.);
	if(PolNP) meancpm = cpm_NP->GetMean();

	cpm_L->Scale(fBGsig/(1.*cpm_L->Integral()));
	cpm_NP->Scale(fNPB/(1.*cpm_NP->Integral()));
	cpm_L->Add(cpm_NP);

	cpm_PSR->Add(cpm_L, -1.);
	if(!PolNP && !PolLSB && !PolRSB) meancpm = cpm_PSR->GetMean();

	std::stringstream meancpmname;
	meancpmname << ";;mean N_{ch}";
	TH1D* h_meancpm = new TH1D("mean_cpm", meancpmname.str().c_str(), 1, 0., 1.);
	h_meancpm->SetBinContent(1, meancpm);
	h_meancpm->Write();
	
	double meanPT = 0;
	// mean pT histos
	pT_PSR->Scale(1./(1.*pT_PSR->Integral()));
	TH1D* pT_SRL = (TH1D*) pT_PSR->Clone();
	pT_SRL->Scale(fSRinPLSB/(1.*pT_SRL->Integral()));
	pT_L->Scale(1./(1.*pT_L->Integral()));
	pT_L->Add(pT_SRL, -1.);
	if(PolLSB) meanPT = pT_L->GetMean();

	TH1D* pT_SRR = (TH1D*) pT_PSR->Clone();
	pT_SRR->Scale(fSRinPRSB/(1.*pT_SRR->Integral()));
	pT_R->Scale(1./(1.*pT_R->Integral()));
	pT_R->Add(pT_SRR, -1.);
	if(PolRSB) meanPT = pT_R->GetMean();

	pT_L->Scale(fracLSB/(1.*pT_L->Integral()));
	pT_R->Scale((1.-fracLSB)/(1.*pT_R->Integral()));
	pT_L->Add(pT_R);

	pT_highct_L->Scale(fracLSB/(1.*pT_highct_L->Integral()));
	pT_highct_R->Scale((1.-fracLSB)/(1.*pT_highct_R->Integral()));
	pT_highct_L->Add(pT_highct_R);

	pT_NP->Scale(1./(1.*pT_NP->Integral()));
	pT_highct_L->Scale(fBGinNP/(1.*pT_highct_L->Integral()));
	pT_NP->Add(pT_highct_L, -1.);
	if(PolNP) meanPT = pT_NP->GetMean();

	pT_L->Scale(fBGsig/(1.*pT_L->Integral()));
	pT_NP->Scale(fNPB/(1.*pT_NP->Integral()));
	pT_L->Add(pT_NP);

	pT_PSR->Add(pT_L, -1.);
	if(!PolNP && !PolLSB && !PolRSB) meanPT = pT_PSR->GetMean();

	std::stringstream meanPTname;
	meanPTname << ";;mean p_{T}";
	TH1D* h_meanPT = new TH1D("mean_pT", meanPTname.str().c_str(), 1, 0., 1.);
	h_meanPT->SetBinContent(1, meanPT);
	h_meanPT->Write();

	// mean y histos
	double meanY = 0;
	rap_PSR->Scale(1./(1.*rap_PSR->Integral()));
	TH1D* rap_SRL = (TH1D*) rap_PSR->Clone();
	rap_SRL->Scale(fSRinPLSB/(1.*rap_SRL->Integral()));
	rap_L->Scale(1./(1.*rap_L->Integral()));
	rap_L->Add(rap_SRL, -1.);
	if(PolLSB) meanY = rap_L->GetMean();

	TH1D* rap_SRR = (TH1D*) rap_PSR->Clone();
	rap_SRR->Scale(fSRinPRSB/(1.*rap_SRR->Integral()));
	rap_R->Scale(1./(1.*rap_R->Integral()));
	rap_R->Add(rap_SRR, -1.);
	if(PolRSB) meanY = rap_R->GetMean();

	rap_L->Scale(fracLSB/(1.*rap_L->Integral()));
	rap_R->Scale((1.-fracLSB)/(1.*rap_R->Integral()));
	rap_L->Add(rap_R);

	rap_highct_L->Scale(fracLSB/(1.*rap_highct_L->Integral()));
	rap_highct_R->Scale((1.-fracLSB)/(1.*rap_highct_R->Integral()));
	rap_highct_L->Add(rap_highct_R);

	rap_NP->Scale(1./(1.*rap_NP->Integral()));
	rap_highct_L->Scale(fBGinNP/(1.*rap_highct_L->Integral()));
	rap_NP->Add(rap_highct_L, -1.);
	if(PolNP) meanY = rap_NP->GetMean();

	rap_L->Scale(fBGsig/(1.*rap_L->Integral()));
	rap_NP->Scale(fNPB/(1.*rap_NP->Integral()));
	rap_L->Add(rap_NP);

	rap_PSR->Add(rap_L, -1.);
	if(!PolNP && !PolLSB && !PolRSB) meanY = rap_PSR->GetMean();

	std::stringstream meanYname;
	meanYname << ";;mean |y|";
	TH1D* h_meanY = new TH1D("mean_y", meanPTname.str().c_str(), 1, 0., 1.);
	h_meanY->SetBinContent(1, meanY);
	h_meanY->Write();

	//---rebin BG to be same as NPBG
	bool ResetBin=true;
	if(ResetBin){
		cout<<"Resetting bins...."<<endl;
		for(int iFrame = 0; iFrame < onia::kNbFrames; iFrame++){

			std::stringstream nameL, nameR, nameNP, nameBGinNPL, nameBGinNPR, title;
			nameL << "hBG_cosThetaPhi_" << onia::frameLabel[iFrame] << "_L";
			nameR << "hBG_cosThetaPhi_" << onia::frameLabel[iFrame] << "_R";
			nameNP << "hNPBG_cosThetaPhi_" << onia::frameLabel[iFrame];
			nameBGinNPL << "hBGinNP_cosThetaPhi_" << onia::frameLabel[iFrame] << "_L";
			nameBGinNPR << "hBGinNP_cosThetaPhi_" << onia::frameLabel[iFrame] << "_R";
			title << ";cos#theta_{"<< onia::frameLabel[iFrame] << "};#phi_{" << onia::frameLabel[iFrame] << "} [deg]";

			hBG_cosThetaPhiL[iFrame] = ReSetBin(hBG_cosThetaPhiL[iFrame], nBinsCosthNPBG, nBinsPhiNPBG, nameL, title);
			hBG_cosThetaPhiR[iFrame] = ReSetBin(hBG_cosThetaPhiR[iFrame], nBinsCosthNPBG, nBinsPhiNPBG, nameR, title);
			hBGinNP_cosThetaPhiL[iFrame] = ReSetBin(hBGinNP_cosThetaPhiL[iFrame], nBinsCosthNPBG, nBinsPhiNPBG, nameBGinNPL, title);
			hBGinNP_cosThetaPhiR[iFrame] = ReSetBin(hBGinNP_cosThetaPhiR[iFrame], nBinsCosthNPBG, nBinsPhiNPBG, nameBGinNPR, title);
		}//iFrame
	}
	//---rebin finished

	//======================================================
	for(int iFrame = 0; iFrame < onia::kNbFrames; iFrame++){
		hBG_cosThetaPhiL[iFrame]->Write();
		hBG_cosThetaPhiR[iFrame]->Write();
		hBGinNP_cosThetaPhiL[iFrame]->Write();
		hBGinNP_cosThetaPhiR[iFrame]->Write();
		hNPBG_cosThetaPhi[iFrame]->Write();
		hSR_cosThetaPhi[iFrame]->Write();

		// combinatorial background in signal region (prompt region)
		// subtract signal contamination in left and right sideband
		std::stringstream nameSRL, nameSRR;
		nameSRL << "hSR_cosThetaPhi_" << onia::frameLabel[iFrame] << "_L";
		nameSRR << "hSR_cosThetaPhi_" << onia::frameLabel[iFrame] << "_R";
		hSR_cosThetaPhiL[iFrame] = (TH2D*)hSR_cosThetaPhi[iFrame]->Clone(nameSRL.str().c_str());
		hSR_cosThetaPhiL[iFrame]->Scale(fSRinPLSB/(1.*hSR_cosThetaPhiL[iFrame]->Integral()));
		hBG_cosThetaPhiL[iFrame]->Scale(1./(1.*hBG_cosThetaPhiL[iFrame]->Integral()));
		hBG_cosThetaPhiL[iFrame] = subtract2D(hBG_cosThetaPhiL[iFrame], hSR_cosThetaPhiL[iFrame]);

		hSR_cosThetaPhiR[iFrame] = (TH2D*)hSR_cosThetaPhi[iFrame]->Clone(nameSRR.str().c_str());
		hSR_cosThetaPhiR[iFrame]->Scale(fSRinPRSB/(1.*hSR_cosThetaPhiR[iFrame]->Integral()));
		hBG_cosThetaPhiR[iFrame]->Scale(1./(1.*hBG_cosThetaPhiR[iFrame]->Integral()));
		hBG_cosThetaPhiR[iFrame] = subtract2D(hBG_cosThetaPhiR[iFrame], hSR_cosThetaPhiR[iFrame]);

		// combination of left and right sideband
		hBG_cosThetaPhiL[iFrame]->Scale(fracLSB/(1.*hBG_cosThetaPhiL[iFrame]->Integral()));
		hBG_cosThetaPhiR[iFrame]->Scale((1.-fracLSB)/(1.*hBG_cosThetaPhiR[iFrame]->Integral()));
		std::stringstream name;
		name << "comb_background_costhphi" << onia::frameLabel[iFrame];
		hBG_cosThetaPhi[iFrame] = (TH2D *) hBG_cosThetaPhiL[iFrame]->Clone(name.str().c_str());
		hBG_cosThetaPhi[iFrame]->Add(hBG_cosThetaPhiR[iFrame]);
		hBG_cosThetaPhi[iFrame]->Write();

		// combinatorial background in high ctau region
		// combination of left and right sideband in high ctau region
		hBGinNP_cosThetaPhiL[iFrame]->Scale(fracLSB/(1.*hBGinNP_cosThetaPhiL[iFrame]->Integral()));
		hBGinNP_cosThetaPhiR[iFrame]->Scale((1.-fracLSB)/(1.*hBGinNP_cosThetaPhiR[iFrame]->Integral()));
		std::stringstream nameBGinNP;
		nameBGinNP << "background_NPR_costhphi" << onia::frameLabel[iFrame];
		hBGinNP_cosThetaPhi[iFrame] = (TH2D *) hBGinNP_cosThetaPhiL[iFrame]->Clone(nameBGinNP.str().c_str());
		hBGinNP_cosThetaPhi[iFrame]->Add(hBGinNP_cosThetaPhiR[iFrame]);
		hBGinNP_cosThetaPhi[iFrame]->Write();

		// non prompt background in high ctau region
		// (hNPBG_cosThetaPhi + hBGinNP_cosThetaPhi)_norm - fBGinNP * (hBGinNP_cosThetaPhi)_norm
		hNPBG_cosThetaPhi[iFrame]->Scale(1./(1.*hNPBG_cosThetaPhi[iFrame]->Integral()));
		hBGinNP_cosThetaPhi[iFrame]->Scale(fBGinNP/(1.*hBGinNP_cosThetaPhi[iFrame]->Integral()));
		std::stringstream nameNPS;
		nameNPS << "background_NPSR_costhphi" << onia::frameLabel[iFrame];
		hNPS_cosThetaPhi[iFrame] = (TH2D *) hNPBG_cosThetaPhi[iFrame]->Clone(nameNPS.str().c_str());
		hNPS_cosThetaPhi[iFrame] = subtract2D(hNPS_cosThetaPhi[iFrame], hBGinNP_cosThetaPhi[iFrame]);
		hNPS_cosThetaPhi[iFrame]->Write();

		// total background
		std::stringstream nameTBG, nameTBGfolded, nameTBGunfolded;
		nameTBG << "background_costhphi" << onia::frameLabel[iFrame];
		nameTBGfolded << "background_folded_costhphi" << onia::frameLabel[iFrame];
		nameTBGunfolded << "background_unfolded_costhphi" << onia::frameLabel[iFrame];
		// for non promp polarization: only use high ctau background
		if(PolLSB || PolRSB){
			hSR_cosThetaPhi[iFrame]->Scale(1./(1.*hSR_cosThetaPhi[iFrame]->Integral()));
			hTBG_cosThetaPhi[iFrame] = (TH2D *) hSR_cosThetaPhi[iFrame]->Clone(nameTBGfolded.str().c_str());
		}
		else if(PolNP)
			hTBG_cosThetaPhi[iFrame] = (TH2D *) hBGinNP_cosThetaPhi[iFrame]->Clone(nameTBGfolded.str().c_str());
		else{
			// fNPBG * (hNPS_cosThetaPhi)_norm + fBGsig * (hBG_cosThetaPhi)_norm
			hBG_cosThetaPhi[iFrame]->Scale(fBGsig/(1.*hBG_cosThetaPhi[iFrame]->Integral()));
			hNPS_cosThetaPhi[iFrame]->Scale(fNPB/(1.*hNPS_cosThetaPhi[iFrame]->Integral()));
			hTBG_cosThetaPhi[iFrame] = (TH2D *) hNPS_cosThetaPhi[iFrame]->Clone(nameTBGfolded.str().c_str());
			hTBG_cosThetaPhi[iFrame]->Add(hBG_cosThetaPhi[iFrame]);
		}
		// write folded histogram to file
		hTBG_cosThetaPhi[iFrame]->Write();

		// get binning of total background histogram
		int nx = hTBG_cosThetaPhi[iFrame]->GetXaxis()->GetNbins();
		int ny = hTBG_cosThetaPhi[iFrame]->GetYaxis()->GetNbins();
		int yPhi = ny/4;
		int xCosTheta = nx/2;

		// unfold the total background histogram
		if(folding){

			if(iFrame == 0){
				std::cout << "---------------------------------------------------" << "\n"
					<< "Total background histogram" << "\n"
					<< "number of cosTheta bins: " << nx << "\n"
					<< "number of phi bins: " << ny << "\n"
					<< "phi bins " << 2*yPhi+1 << " to " << 3*yPhi << " are filled." << "\n"
					<< "---------------------------------------------------" << std::endl;
			}

			for (int j = 0; j <= nx; j++){
				for (int k = 2*yPhi+1; k <= 3*yPhi; k++){

					double c = hTBG_cosThetaPhi[iFrame]->GetBinContent(j,k);
					double e = hTBG_cosThetaPhi[iFrame]->GetBinError(j,k);

					// flip in cosTheta
					double l = nx + 1 - j;

					// set bin content and error of phiFolded in the other 3 (not yet filled) phi regions
					// 90 - 180: flip phi (upwards), flip cosTheta
					hTBG_cosThetaPhi[iFrame]->SetBinContent(l,6*yPhi+1-k,c);
					hTBG_cosThetaPhi[iFrame]->SetBinError(l,6*yPhi+1-k,e);
					// 0 - -90: flip phi (downwards)
					hTBG_cosThetaPhi[iFrame]->SetBinContent(j,ny+1-k,c);
					hTBG_cosThetaPhi[iFrame]->SetBinError(j,ny+1-k,e);
					// -90 - -180: flip cosTheta, shift phi
					hTBG_cosThetaPhi[iFrame]->SetBinContent(l,k-2*yPhi,c);
					hTBG_cosThetaPhi[iFrame]->SetBinError(l,k-2*yPhi,e);

				}
			}
		} // folding
		// write unfolded histogram to file
		hTBG_cosThetaPhi[iFrame]->SetName(nameTBGunfolded.str().c_str());
		hTBG_cosThetaPhi[iFrame]->Write();

		if(!normApproach){
			hTBG_cosThetaPhi[iFrame]->SetName(nameTBG.str().c_str());
			hTBG_cosThetaPhi[iFrame]->Write();
		}

		std::stringstream title, nameSR;
		title << ";cos#theta_{"<< onia::frameLabel[iFrame] << "};#phi_{" << onia::frameLabel[iFrame] << "} [deg]";
		nameSR << "hSR_rebinned_cosThetaPhi_" << onia::frameLabel[iFrame];
		// in case of normalization approach
		// set binning of background histogram to 16 x 16 if smaller than 16 x 16
		if(normApproach){
			if(nx < 16) nx = 16;
			if(ny < 16) ny = 16;
		}
		hTBG_cosThetaPhi[iFrame] = ReSetBin(hTBG_cosThetaPhi[iFrame], nx, ny, nameTBG, title);

		if(iFrame == 0){
			std::cout << "----------------------------------------------" << "\n"
				<< "Final binning of total background histogram: " << "\n"
				<< "cosTheta: " << nx << "\n"
				<< "phi: " << ny << "\n"
				<< "----------------------------------------------" << std::endl;
		}

		//------------ Normalization issue: set bins that are not filled in hSR to 0 in hTBG
		hSR_cosThetaPhi[iFrame] = (TH2D*)hTBG_cosThetaPhi[iFrame]->Clone(nameSR.str().c_str());

		// delete contents of hSR and fill it with events from prompt signal region
		for (int j = 0; j <= nx; j++){
			for (int k = 0; k <= ny; k++){
				hSR_cosThetaPhi[iFrame]->SetBinContent(j,k,0);
				hSR_cosThetaPhi[iFrame]->SetBinError(j,k,0);
			}
		}

	} // iFrame

	// loop through tree and fill hSR histogram
	std::cout << "Filling prompt signal region histogram" << std::endl;
	for(int i = 0; i < n; i++){

		long iEntry = intree->LoadTree(i);
		intree->GetEntry(iEntry);
		if(i % 100000 == 0) std::cout << "entry " << i << " out of " << n << std::endl;

		if(jpsi->Pt() >= onia::pTRange[rapBin-1][ptBin-1] &&
				jpsi->Pt() < onia::pTRange[rapBin-1][ptBin] &&
				TMath::Abs(jpsi->Rapidity()) >= onia::rapForPTRange[rapBin-1] &&
				TMath::Abs(jpsi->Rapidity()) < onia::rapForPTRange[rapBin]){

			if(PolLSB){
				if(jpsi->M() >= massMinL && jpsi->M() < massMaxL && jpsict < ctauCut){
					calcPol(*lepP, *lepN);
					for(int iFrame = 0; iFrame < onia::kNbFrames; iFrame++)
						hSR_cosThetaPhi[iFrame]->Fill(thisCosTh[iFrame], thisPhi[iFrame]);
				}
			}
			else if(PolRSB){
				if(jpsi->M() >= massMinR && jpsi->M() <= massMaxR && jpsict < ctauCut){
					calcPol(*lepP, *lepN);
					for(int iFrame = 0; iFrame < onia::kNbFrames; iFrame++)
						hSR_cosThetaPhi[iFrame]->Fill(thisCosTh[iFrame], thisPhi[iFrame]);
				}
			}
			else if(PolNP){
				if(jpsi->M() >= massMinSR && jpsi->M() <= massMaxSR && jpsict >= ctauCut){
					calcPol(*lepP, *lepN);
					for(int iFrame = 0; iFrame < onia::kNbFrames; iFrame++)
						hSR_cosThetaPhi[iFrame]->Fill(thisCosTh[iFrame], thisPhi[iFrame]);
				}
			}
			else{
				if(jpsi->M() >= massMinSR && jpsi->M() < massMaxSR && jpsict < ctauCut){
					calcPol(*lepP, *lepN);
					for(int iFrame = 0; iFrame < onia::kNbFrames; iFrame++)
						hSR_cosThetaPhi[iFrame]->Fill(thisCosTh[iFrame], thisPhi[iFrame]);
				}
			}

		} // if(onia...)
	} // for(int i ..)

	// set bins in hTBG_cosThetaPhi to 0 when bin is 0 in hSR_cosThetaPhi
	std::cout << "Setting bins in total background histogram to 0 when they are unfilled in prompt signal region histogram:" << std::endl;
	for(int iFrame = 0; iFrame < onia::kNbFrames; iFrame++){

		hSR_cosThetaPhi[iFrame]->Write();

		int nx = hTBG_cosThetaPhi[iFrame]->GetXaxis()->GetNbins();
		int ny = hTBG_cosThetaPhi[iFrame]->GetYaxis()->GetNbins();
		int zeroBins = 0;

		for (int j = 0; j <= nx; j++){
			for (int k = 0; k <= ny; k++){
				double c1 = hSR_cosThetaPhi[iFrame]->GetBinContent(j,k);
				double c2 = hTBG_cosThetaPhi[iFrame]->GetBinContent(j,k);
				double e2 = hTBG_cosThetaPhi[iFrame]->GetBinError(j,k);
				std::cout << c2 << " " << e2 << std::endl;
				if (c1 == 0 && c2 != 0){
					hTBG_cosThetaPhi[iFrame]->SetBinContent(j,k,0);
					hTBG_cosThetaPhi[iFrame]->SetBinError(j,k,0);
					zeroBins++;
					std::cout << "bin " << j << ", " << k << " was set to 0" << std::endl;
				}
			} // k
		} // j
		std::cout << "In the " << onia::frameLabel[iFrame] << " frame, " << zeroBins << " of " << nx*ny << " bins were set to 0." << std::endl;

		// in case of normalization approach:
		//write final background histogram with bins set to 0 to file
		if(normApproach){
			std::stringstream nameTBG;
			nameTBG << "background_costhphi" << onia::frameLabel[iFrame];
			hTBG_cosThetaPhi[iFrame]->SetName(nameTBG.str().c_str());
			hTBG_cosThetaPhi[iFrame]->Write();
		}

	}// iFrame

	outtree->Write();
	datafile->Close();
	fitfile->Close();

} // void

//-----------------------------------------------------------------------------------------------------------
vector<double> calculateInte(RooWorkspace *ws, RooDataSet *dataJpsictErr, double ctCutMin, double ctCutMax){

	double ctMin = -2., ctMax = 6.;
	int fineBins = 8000;

	RooAbsPdf *PRpdf = (RooAbsPdf*)ws->pdf("TotalPromptLifetime");
	RooAbsPdf *NPpdf = (RooAbsPdf*)ws->pdf("nonPromptSSD");
	RooAbsPdf *BGpdf = (RooAbsPdf*)ws->pdf("backgroundlifetime");

	RooRealVar *Jpsict = ws->var("Jpsict");
	RooRealVar *JpsictErr = ws->var("JpsictErr");
	Jpsict->setMin(ctMin);
	Jpsict->setMax(ctMax);

	RooDataSet *genDataPR = PRpdf->generate(*Jpsict,ProtoData(*dataJpsictErr));
	RooDataSet *genDataNP = NPpdf->generate(*Jpsict,ProtoData(*dataJpsictErr));
	RooDataSet *genDataBG = BGpdf->generate(*Jpsict,ProtoData(*dataJpsictErr));

	TH2F* histPR2D = (TH2F*)genDataPR->createHistogram("histPR2D",*Jpsict,Binning(fineBins),YVar(*JpsictErr,Binning(fineBins/10)));
	TH1F* histPR   = (TH1F*)histPR2D->ProjectionX();
	TH2F* histNP2D = (TH2F*)genDataNP->createHistogram("histNP2D",*Jpsict,Binning(fineBins),YVar(*JpsictErr,Binning(fineBins/10)));
	TH1F* histNP   = (TH1F*)histNP2D->ProjectionX();
	TH2F* histBG2D = (TH2F*)genDataBG->createHistogram("histBG2D",*Jpsict,Binning(fineBins),YVar(*JpsictErr,Binning(fineBins/10)));
	TH1F* histBG   = (TH1F*)histBG2D->ProjectionX();

	histPR->SetLineColor(kRed);
	histPR->SetMarkerColor(kRed);
	histNP->SetLineColor(kBlue);
	histNP->SetMarkerColor(kBlue);
	histBG->SetLineColor(kBlack);
	histBG->SetMarkerColor(kBlack);

	histPR->Scale(1./(1.*histPR->Integral()));
	histNP->Scale(1./(1.*histNP->Integral()));
	histBG->Scale(1./(1.*histBG->Integral()));

	int BinLow = 0, BinHigh = 0;
	bool getLow = false, getHigh = false;
	for(int bin = 0; bin < fineBins; bin++){
		if(ctCutMin > histPR->GetBinLowEdge(bin) && ctCutMin < histPR->GetBinLowEdge(bin+1)){
			BinLow = bin; getLow = true;}
		if(ctCutMax > histPR->GetBinLowEdge(bin) && ctCutMax < histPR->GetBinLowEdge(bin+1)){
			BinHigh = bin; getHigh = true;}
		if(getLow && getHigh) break;
	}

	double IntePR    = histPR->Integral(BinLow,BinHigh);
	double InteNP    = histNP->Integral(BinLow,BinHigh);
	double InteBG    = histBG->Integral(BinLow,BinHigh);

	vector<double> InteRlts;
	InteRlts.push_back(IntePR);
	InteRlts.push_back(InteNP);
	InteRlts.push_back(InteBG);

	delete histPR;
	delete histNP;
	delete histBG;
	delete histPR2D;
	delete histNP2D;
	delete histBG2D;

	return InteRlts;
}

vector<double> calculateInteSB(RooWorkspace *ws, RooDataSet *dataJpsictErr, double ctCutMin, double ctCutMax){

	double ctMin = -2., ctMax = 6.;
	int fineBins = 8000;

	RooAbsPdf *BGpdfLSB = (RooAbsPdf*)ws->pdf("backgroundlifetimeLpre");
	RooAbsPdf *LSBpdf = (RooAbsPdf*)ws->pdf("backgroundlifetimeL");
	RooAbsPdf *BGpdfRSB = (RooAbsPdf*)ws->pdf("backgroundlifetimeRpre");
	RooAbsPdf *RSBpdf = (RooAbsPdf*)ws->pdf("backgroundlifetimeR");

	RooRealVar *Jpsict = ws->var("Jpsict");
	RooRealVar *JpsictErr = ws->var("JpsictErr");
	Jpsict->setMin(ctMin);
	Jpsict->setMax(ctMax);

	RooDataSet *genDataBGLSB = BGpdfLSB->generate(*Jpsict,ProtoData(*dataJpsictErr));
	RooDataSet *genDataBGRSB = BGpdfRSB->generate(*Jpsict,ProtoData(*dataJpsictErr));
	RooDataSet *genDataLSB   = LSBpdf  ->generate(*Jpsict,ProtoData(*dataJpsictErr));
	RooDataSet *genDataRSB   = RSBpdf  ->generate(*Jpsict,ProtoData(*dataJpsictErr));
	
	TH2F* histBGLSB2D = (TH2F*)genDataBGLSB->createHistogram("histBGLSB2D",*Jpsict,Binning(fineBins),
			YVar(*JpsictErr,Binning(fineBins/10)));
	TH1F* histBGLSB   = (TH1F*)histBGLSB2D->ProjectionX();
	TH2F* histBGRSB2D = (TH2F*)genDataBGRSB->createHistogram("histBGRSB2D",*Jpsict,Binning(fineBins),
			YVar(*JpsictErr,Binning(fineBins/10)));
	TH1F* histBGRSB   = (TH1F*)histBGRSB2D->ProjectionX();
	TH2F* histLSB2D = (TH2F*)genDataLSB->createHistogram("histLSB2D",*Jpsict,Binning(fineBins),
			YVar(*JpsictErr,Binning(fineBins/10)));
	TH1F* histLSB   = (TH1F*)histLSB2D->ProjectionX();
	TH2F* histRSB2D = (TH2F*)genDataRSB->createHistogram("histRSB2D",*Jpsict,Binning(fineBins),
			YVar(*JpsictErr,Binning(fineBins/10)));
	TH1F* histRSB   = (TH1F*)histRSB2D->ProjectionX();

	histBGLSB->Scale(1./(1.*histBGLSB->Integral()));
	histBGRSB->Scale(1./(1.*histBGRSB->Integral()));
	histLSB->Scale(1./(1.*histLSB->Integral()));
	histRSB->Scale(1./(1.*histRSB->Integral()));

	int BinLow = 0, BinHigh = 0;
	bool getLow = false, getHigh = false;
	for(int bin = 0; bin < fineBins; bin++){
		if(ctCutMin > histBGLSB->GetBinLowEdge(bin) && ctCutMin < histBGLSB->GetBinLowEdge(bin+1)){
			BinLow = bin; getLow = true;}
		if(ctCutMax > histBGLSB->GetBinLowEdge(bin) && ctCutMax < histBGLSB->GetBinLowEdge(bin+1)){
			BinHigh = bin; getHigh = true;}
		if(getLow && getHigh) break;
	}

	double InteBGLSB  = histBGLSB->Integral(BinLow,BinHigh);
	double InteBGRSB  = histBGRSB->Integral(BinLow,BinHigh);
	double InteLSB    = histLSB  ->Integral(BinLow,BinHigh);
	double InteRSB    = histRSB  ->Integral(BinLow,BinHigh);

	vector<double> InteRlts;
	InteRlts.push_back(InteBGLSB);
	InteRlts.push_back(InteBGRSB);
	InteRlts.push_back(InteLSB);
	InteRlts.push_back(InteRSB);

	delete histBGLSB;
	delete histBGRSB;
	delete histLSB;
	delete histRSB;
	delete histBGLSB2D;
	delete histBGRSB2D;
	delete histLSB2D;
	delete histRSB2D;

	return InteRlts;
}

//=================================================
//=================================================
// manual subtraction
TH3D *subtract3D(TH3D* hist1, TH3D* hist2){

	int nx = hist1->GetXaxis()->GetNbins();
	int ny = hist1->GetYaxis()->GetNbins();
	int nz = hist1->GetZaxis()->GetNbins();

	for (int j = 0; j <= nx; j++){
		for (int k = 0; k <= ny; k++){
			for(int l = 0; l <= nz; l++){

				double c1 = hist1->GetBinContent(j,k,l);
				if (c1 > 0) {
					double c2 = hist2->GetBinContent(j,k,l);
					double c3 = c1 - 1.*c2;
					double e1 = hist1->GetBinError(j,k,l);  
					double e2 = hist2->GetBinError(j,k,l);  
					double e3 = TMath::Sqrt(e1*e1 + e2*e2); 
					if(c3 < 0){                             
						c3 = 0;                             
						e3 = 0;                             
					}
					hist1->SetBinContent(j,k,l,c3);
					hist1->SetBinError(j,k,l,e3);
				}

			} // j
		} // k
	} // l

	return hist1;
}

//=================================================
TH2D *subtract2D(TH2D* hist1, TH2D* hist2){

	int nx = hist1->GetXaxis()->GetNbins();
	int ny = hist1->GetYaxis()->GetNbins();

	for (int j = 0; j <= nx; j++){
		for (int k = 0; k <= ny; k++){
			double c1 = hist1->GetBinContent(j,k);
			if (c1 > 0) {
				double c2 = hist2->GetBinContent(j,k);
				double c3 = c1 - 1.*c2;
				double e1 = hist1->GetBinError(j,k);
				double e2 = hist2->GetBinError(j,k);
				double e3 = TMath::Sqrt(e1*e1 + e2*e2); 
				if(c3 < 0){
					c3 = 0;
					e3 = 0;
				}
				hist1->SetBinContent(j,k,c3);
				hist1->SetBinError(j,k, e3);
			}
		}
	}

	return hist1;
}


//=================================================
TH2D* ReSetBin(TH2D* hist, int nBinX, int nBinY, const std::stringstream& name, const std::stringstream& title){
	TH2D *tempHist = (TH2D*)hist->Clone("temp_BG_cosThetaPhiL");
	delete hist;
	hist = new TH2D(name.str().c_str(), title.str().c_str(),
			nBinX, onia::cosTMin, onia::cosTMax, nBinY, onia::phiPolMin, onia::phiPolMax);
	hist->Sumw2();
	TAxis *Xold = tempHist->GetXaxis();
	TAxis *Yold = tempHist->GetYaxis();
	TAxis *Xnew = hist->GetXaxis();
	TAxis *Ynew = hist->GetYaxis();
	for(int binX = 1; binX <= Xnew->GetNbins(); binX++){
		for(int binY = 1; binY <= Ynew->GetNbins(); binY++){
			double centerX = Xnew->GetBinCenter(binX);
			double centerY = Ynew->GetBinCenter(binY);

			//find the corresponding bin and bin error
			double binCont=0.,binErr=0.;
			bool findBin=false;
			for(int BinX = 1; BinX <= Xold->GetNbins(); BinX++){
				for(int BinY = 1; BinY <= Yold->GetNbins(); BinY++){
					double lowX = Xold->GetBinLowEdge(BinX);
					double upX  = Xold->GetBinUpEdge(BinX);
					double lowY = Yold->GetBinLowEdge(BinY);
					double upY  = Yold->GetBinUpEdge(BinY);
					if(centerX > lowX && centerX < upX && centerY > lowY && centerY < upY){
						binCont = tempHist->GetBinContent(BinX,BinY);
						binErr = tempHist->GetBinError(BinX,BinY);
						findBin=true;
					}
					if(findBin) break;
				}//BinY
			}//BinX
			//done
			hist->SetBinContent(binX,binY,binCont);
			hist->SetBinError(binX,binY,binErr);
		}//binY
	}//binX

	return hist;
}


//=================================================
int order(int n){
	int total=1;
	for(int i=0;i<n;i++)
		total=total*2;
	return total;
}

int findEvenNum(double number){
	int thisNum=0;
	for(int n=0;n<100;n++){
		if(number >= order(n) && number <= order(n+1)){
			thisNum=order(n+1);
			break;
		}
	}
	return thisNum;
}

double calcuFracL(RooWorkspace *ws, double mean, double sigma){
	RooRealVar JpsiMass(*ws->var("JpsiMass"));
	double sigMaxMass = mean+sigma*onia::nSigMass;
	double sigMinMass = mean-sigma*onia::nSigMass;
	double sbHighMass = mean+sigma*onia::nSigBkgHigh;
	double sbLowMass =  mean-sigma*onia::nSigBkgLow;
	JpsiMass.setRange("SR",sigMinMass,sigMaxMass);
	JpsiMass.setRange("LSB",JpsiMass.getMin(), sbLowMass);
	JpsiMass.setRange("RSB",sbHighMass, JpsiMass.getMax()); 

	RooAddPdf *bkgMassShape = (RooAddPdf*)ws->pdf("bkgMassShape");

	RooRealVar* InteLSB = (RooRealVar*)bkgMassShape->createIntegral(JpsiMass,NormSet(JpsiMass),Range("LSB"));
	RooRealVar* InteRSB = (RooRealVar*)bkgMassShape->createIntegral(JpsiMass,NormSet(JpsiMass),Range("RSB"));
	RooRealVar* InteSR  = (RooRealVar*)bkgMassShape->createIntegral(JpsiMass,NormSet(JpsiMass),Range("SR"));
	double inteLSB = InteLSB->getVal();
	double inteRSB = InteRSB->getVal();
	double inteSR  = InteSR->getVal();

	double mean_LSB, mean_RSB, mean_SR;

	double MassDist=0.001;
	for(int i=0; i<10000; i++){
		JpsiMass.setRange(Form("LSB_%d",i),JpsiMass.getMin(), JpsiMass.getMin()+i*MassDist);
		RooRealVar* tempLSB = (RooRealVar*)bkgMassShape->createIntegral(JpsiMass,NormSet(JpsiMass),Range(Form("LSB_%d",i)));
		if(tempLSB->getVal() > inteLSB/2.){
			mean_LSB = JpsiMass.getMin()+i*MassDist;
			break;
		}
	}
	for(int i=0; i<10000; i++){
		JpsiMass.setRange(Form("RSB_%d",i),sbHighMass,sbHighMass+i*MassDist);
		RooRealVar* tempRSB = (RooRealVar*)bkgMassShape->createIntegral(JpsiMass,NormSet(JpsiMass),Range(Form("RSB_%d",i)));
		if(tempRSB->getVal() > inteRSB/2.){
			mean_RSB = sbHighMass+i*MassDist;
			break;
		}
	}
	for(int i=0; i<10000; i++){
		JpsiMass.setRange(Form("SR_%d",i),sigMinMass,sigMinMass+i*MassDist);
		RooRealVar* tempSR = (RooRealVar*)bkgMassShape->createIntegral(JpsiMass,NormSet(JpsiMass),Range(Form("SR_%d",i)));
		if(tempSR->getVal() > inteSR/2.){
			mean_SR = sigMinMass+i*MassDist;
			break;
		}
	}
	double fracLCentral = (mean_RSB-mean_SR)/(mean_RSB-mean_LSB);
	cout<<"mean_LSB: "<<mean_LSB<<endl;
	cout<<"mean_RSB: "<<mean_RSB<<endl;
	cout<<"mean_SR: "<<mean_SR<<endl;
	cout<<"fracLCentral: "<<fracLCentral<<endl;

	return fracLCentral;
}

