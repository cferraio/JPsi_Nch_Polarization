#include <iostream>
#include <string>
#include <sstream>
#include "calculatePar.cc"
using namespace std;
using namespace RooFit;
using namespace onia;

double calcuFracL(RooWorkspace *ws, double mean, double sigma);
void PlotMassPar(int  nState=4);
void PlotLifePar(int  nState=4);
void PlotBFrac_1S(int  nState=4);
void PlotBFrac_2S(int  nState=5);

// evaluate Ctau cut to define PR and NP regions
double ctMin = -2., ctMax = 6.;
int fineBins=8000;
double legendsize=0.035;
vector<double> getSigma(RooWorkspace *ws, RooDataSet *dataJpsictErr, int rapBin, int ptBin);
void evaluate(double nSigma=2.5, int nState=4, int type=0, bool doCtauUncer=false); 
// type=0: PR(-ctau,ctau); type=1: NP(ctau,+infinity)
void plotEval(double nSigma=2.5, int nState=4, int type=0); 
void evaluateCtauCut(double nSigma=2.5, int nState=4, int type=0, bool doCtauUncer=false);
vector<double> calculateInte(RooWorkspace *ws, RooDataSet *dataJpsictErr, double ctCutMin, double ctCutMax);

//========================================================
// code to read input arguments
	template<typename T>
void fromSplit(const std::string& key, const std::string &arg, T& out)
{
	const char delim = '=';
	// Skip if key or delimiter not there
	if ((arg.find(key) == std::string::npos) ||
			(arg.find(delim) == std::string::npos))
		return;

	std::string skey, sval;
	std::stringstream sstr(arg);
	std::getline(sstr, skey, delim); // Dummy read to skip key
	std::getline(sstr, sval, delim); // Get value
	T tout;
	if (!(std::istringstream(sval) >> std::boolalpha >> tout))
		return;
	out = tout;
	std::cout << std::boolalpha << skey << ": "  << out << std::endl;
}

// Special version for string without the conversion 
	template<>
void fromSplit(const std::string& key, const std::string &arg, std::string &out)
{
	const char delim = '=';
	// Skip if key or delimiter not there
	if ((arg.find(key) == std::string::npos) ||
			(arg.find(delim) == std::string::npos))
		return;
	std::string skey, sval;
	std::stringstream sstr(arg);
	std::getline(sstr, skey, delim); // Dummy read to skip key
	std::getline(sstr, sval, delim); // Get value
	out = sval;
	std::cout << skey << ": "  << out << std::endl;
}

//=====================================================================
int main(int argc, char* argv[]){
	// set default values
	int nState = 999;
	bool doCtauUncer = false;

	// Loop over argument list                                                                                                                                                       
	for (int i=1; i < argc; i++){
		std::string arg = argv[i];
		fromSplit("nState", arg, nState);
		fromSplit("doCtauUncer", arg, doCtauUncer);
	}

	PlotMassPar(nState);
	PlotLifePar(nState);
	if(nState==4)
		PlotBFrac_1S(nState);
	if(nState==5)
		PlotBFrac_2S(nState);

	double nSigma=3.0;
	if(nState==4) nSigma=3.0; //2.5
	if(nState==5) nSigma=3.0; //2.0
	//nSigma = -1;
	//evaluateCtauCut(nSigma, nState, 0, doCtauUncer);
	evaluateCtauCut(nSigma, nState, 1, doCtauUncer);

	return 0;
}


//=============================================
void PlotMassPar(int  nState){
	int RapBins = onia::kNbRapForPTBins,
			PtBins  = onia::kNbPTMaxBins;

	std::stringstream savePath;
	savePath << "Fit/parameter/mass";
	gSystem->mkdir(savePath.str().c_str(),kTRUE);

	double pT[RapBins][PtBins];
	double mean[RapBins][PtBins];
	double sigmaWei[RapBins][PtBins];
	double bkgRatio3Sig[RapBins][PtBins];
	double sigma1[RapBins][PtBins];
	double sigma2[RapBins][PtBins];
	double evtBkgSB[RapBins][PtBins];
	double fracLSB[RapBins][PtBins];
	double alphaCB[RapBins][PtBins];
	double lambdaBG[RapBins][PtBins];
	double fracSigInLSB[RapBins][PtBins];
	double fracSigInRSB[RapBins][PtBins];
	double evtInLSB[RapBins][PtBins];
	double fracBkgInLSB[RapBins][PtBins];
	double fracBkgInRSB[RapBins][PtBins];

	double pTErr[RapBins][PtBins];
	double meanErr[RapBins][PtBins];
	double sigmaWeiErr[RapBins][PtBins];
	double bkgRatio3SigErr[RapBins][PtBins];
	double sigma1Err[RapBins][PtBins];
	double sigma2Err[RapBins][PtBins];
	double evtBkgSBErr[RapBins][PtBins];
	double fracLSBErr[RapBins][PtBins];
	double alphaCBErr[RapBins][PtBins];
	double lambdaBGErr[RapBins][PtBins];
	double fracSigInLSBErr[RapBins][PtBins];
	double fracSigInRSBErr[RapBins][PtBins];
	double evtInLSBErr[RapBins][PtBins];
	double fracBkgInLSBErr[RapBins][PtBins];
	double fracBkgInRSBErr[RapBins][PtBins];

	for(int rapBin = 1; rapBin < RapBins+1; rapBin++){
		for(int ptBin = 1; ptBin < PtBins+1; ptBin++){

			double ptMin=0.,ptMax=0.;
			ptMin = onia::pTRange[rapBin][ptBin-1];
			ptMax = onia::pTRange[rapBin][ptBin];
			pT[rapBin-1][ptBin-1] =  (ptMin + ptMax)/2.;
			pTErr[rapBin-1][ptBin-1] =  (ptMax - ptMin)/2.;

			mean[rapBin-1][ptBin-1]         = 0.;
			sigmaWei[rapBin-1][ptBin-1]     = 0.;
			bkgRatio3Sig[rapBin-1][ptBin-1] = 0.;
			sigma1[rapBin-1][ptBin-1]       = 0.;
			sigma2[rapBin-1][ptBin-1]       = 0.;
			evtBkgSB[rapBin-1][ptBin-1]     = 0.;
			fracLSB[rapBin-1][ptBin-1]      = 0.;
			fracSigInLSB[rapBin-1][ptBin-1]   = 0.;
			fracSigInRSB[rapBin-1][ptBin-1]   = 0.;
			evtInLSB[rapBin-1][ptBin-1]   = 0.;
			fracBkgInLSB[rapBin-1][ptBin-1]   = 0.;
			fracBkgInRSB[rapBin-1][ptBin-1]   = 0.;

			meanErr[rapBin-1][ptBin-1]         = 0.;
			sigmaWeiErr[rapBin-1][ptBin-1]     = 0.;
			bkgRatio3SigErr[rapBin-1][ptBin-1] = 0.;
			sigma1Err[rapBin-1][ptBin-1]       = 0.;
			sigma2Err[rapBin-1][ptBin-1]       = 0.;
			evtBkgSBErr[rapBin-1][ptBin-1]     = 0.;
			fracLSBErr[rapBin-1][ptBin-1]      = 0.;
			fracSigInLSBErr[rapBin-1][ptBin-1]   = 0.;
			fracSigInRSBErr[rapBin-1][ptBin-1]   = 0.;
			evtInLSBErr[rapBin-1][ptBin-1]   = 0.;
			fracBkgInLSBErr[rapBin-1][ptBin-1]   = 0.;
			fracBkgInRSBErr[rapBin-1][ptBin-1]   = 0.;
		}
	}

	TFile *inFile;
	RooWorkspace *ws;
	RooDataSet *data;
	char inName[200];

	for(int rapBin = 1; rapBin < RapBins+1; rapBin++){
		for(int ptBin = 1; ptBin < PtBins+1; ptBin++){

			sprintf(inName, Form("tmpFiles/fit_Psi%dS_rap%d_pt%d.root",nState-3,rapBin,ptBin));

			cout<<"inName: "<<inName<<endl;
			inFile=new TFile(inName,"R");
			if(!inFile) continue;
			ws=(RooWorkspace *)inFile->Get("ws_masslifetime");
			if(!ws){ cout<<">>=======Error: no workspace in root file=========="<<endl; continue; }
			RooRealVar JpsiMass(*ws->var("JpsiMass"));
			data=(RooDataSet *)ws->data(Form("data_rap%d_pt%d",rapBin,ptBin));
			cout<<"data is "<<data<<endl;

			RooRealVar *CBmass=(RooRealVar *)ws->var("CBmass");
			RooRealVar *CBsigma=(RooRealVar *)ws->var("CBsigma");
			RooRealVar *CBsigma2=(RooRealVar *)ws->var("CBsigma2");
			RooRealVar *CBalpha=(RooRealVar *)ws->var("CBalpha");
			RooRealVar *CBn=(RooRealVar *)ws->var("CBn");
			RooRealVar *bkgLambda=(RooRealVar *)ws->var("bkgLambda");
			RooRealVar *fracCB1_=(RooRealVar *)ws->var("fracCB1");
			RooRealVar *fracBkg_=(RooRealVar *)ws->var("fracBkg");

			double Mean = CBmass->getVal();
			double MeanErr = CBmass->getError();
			double Sigma = CBsigma->getVal();
			double SigmaErr = CBsigma->getError();
			double Sigma2 = CBsigma2->getVal();
			double Sigma2Err = CBsigma2->getError();
			double Alpha = CBalpha->getVal();
			double AlphaErr = CBalpha->getError();
			double cbN = CBn->getVal();
			double cbNErr = CBn->getError();
			double lambda = bkgLambda->getVal();
			double lambdaErr = bkgLambda->getError();
			double fracCB1 = fracCB1_->getVal();
			double fracCB1Err = fracCB1_->getError();
			double fracBkg = fracBkg_->getVal();
			double fracBkgErr = fracBkg_->getError();
			double SigmaWei = sqrt(pow(Sigma,2)*fracCB1+pow(Sigma2,2)*(1-fracCB1));
			double SigmaWeiErr =  (1./(2*SigmaWei))*
				sqrt(pow((pow(Sigma,2)-pow(Sigma2,2))*fracCB1Err,2) + pow(2*fracCB1*Sigma*SigmaErr,2) + pow(2*(1-fracCB1)*Sigma2*Sigma2Err,2) );

			if(Sigma > Sigma2){
				fracCB1 = 1-fracCB1;
				double temp = 0.;
				temp = Sigma; Sigma = Sigma2; Sigma2 = temp;
				temp = SigmaErr; SigmaErr = Sigma2Err; Sigma2Err = temp;
			}

			double sigMaxMass = Mean+SigmaWei*onia::nSigMass;
			double sigMinMass = Mean-SigmaWei*onia::nSigMass;
			double sbHighMass = Mean+SigmaWei*onia::nSigBkgHigh;
			double sbLowMass =  Mean-SigmaWei*onia::nSigBkgLow;
			cout<<"Mean: "<<Mean<<" SigmaWei: "<<SigmaWei<<endl;
			cout<<"Mean-5.5*SigmaWei: "<<Mean-SigmaWei*5.5<<endl;
			//double MassMin = onia::massMin; // = JpsiMass.getMin()
			//double MassMax = onia::massMax; // = JpsiMass.getMax()

			int nEntries = data->numEntries();
			JpsiMass.setRange("SR",sigMinMass,sigMaxMass);
			JpsiMass.setRange("SBL",JpsiMass.getMin(), sbLowMass);
			JpsiMass.setRange("SBR",sbHighMass, JpsiMass.getMax());

			// create datasets for LSB, RSB and SR
			std::stringstream cutSR, cutSBL, cutSBR;
			cutSR << "JpsiMass > " << sigMinMass << " && JpsiMass < " << sigMaxMass;
			cutSBL << "JpsiMass > " << onia::massMin << " && JpsiMass < " << sbLowMass;
			cutSBR << "JpsiMass > " << sbHighMass << " && JpsiMass < " << onia::massMax;

			std::stringstream binNameSR, binNameSBL, binNameSBR;
			binNameSR  << "data_rap" << rapBin << "_pt" << ptBin << "_SR";
			binNameSBL << "data_rap" << rapBin << "_pt" << ptBin << "_SBL";
			binNameSBR << "data_rap" << rapBin << "_pt" << ptBin << "_SBR";

			RooAbsData* dataSR, *dataSBL, *dataSBR;
			dataSR  = data->reduce(Cut(cutSR.str().c_str()));
			dataSBL = data->reduce(Cut(cutSBL.str().c_str()));
			dataSBR = data->reduce(Cut(cutSBR.str().c_str()));

			dataSR->SetNameTitle(binNameSR.str().c_str(), "data in signal region");
			dataSBL->SetNameTitle(binNameSBL.str().c_str(), "data in LSB");
			dataSBR->SetNameTitle(binNameSBR.str().c_str(), "data in RSB");

			TH1* histSR =  dataSR->createHistogram("histSR", JpsiMass,  Binning(120));
			TH1* histSBL = dataSBL->createHistogram("histSBL", JpsiMass, Binning(120));
			TH1* histSBR = dataSBR->createHistogram("histSBR", JpsiMass, Binning(120));

			double meanSR = histSR->GetMean();
			double meanSBL = histSBL->GetMean();
			double meanSBR = histSBR->GetMean();

			RooAddPdf *massPdf = (RooAddPdf*)ws->pdf("massModel");
			RooAbsPdf *bkgMassShape = (RooAbsPdf*)ws->pdf("bkgMassShape");
			RooAddPdf *sigMassShape = (RooAddPdf*)ws->pdf("sigMassShape");

			RooRealVar* fracFull3Sig=(RooRealVar*)massPdf->createIntegral(JpsiMass,NormSet(JpsiMass),Range("SR"));
			RooRealVar* fracBkg3Sig=(RooRealVar*)bkgMassShape->createIntegral(JpsiMass,NormSet(JpsiMass),Range("SR"));
			double FracFull3Sig = fracFull3Sig->getVal();
			double FracFull3SigErr = fracFull3Sig->getError();
			double FracBkg3Sig = fracBkg3Sig->getVal();
			double FracBkg3SigErr = fracBkg3Sig->getError();

			double evtFull3Sig = nEntries*FracFull3Sig;
			double evtBkg3Sig = nEntries*fracBkg*FracBkg3Sig;
			double BkgRatio3Sig = evtBkg3Sig/evtFull3Sig; //fracBkg*FracBkg3Sig/FracFull3Sig
			double BkgRatio3SigErr = BkgRatio3Sig*
				sqrt(pow(fracBkgErr/fracBkg,2)+pow(FracBkg3SigErr/FracBkg3Sig,2)+pow(FracFull3SigErr/FracFull3Sig,2));

			RooRealVar *fracBkgSBL = (RooRealVar*)bkgMassShape->createIntegral(JpsiMass,NormSet(JpsiMass),Range("SBL"));
			RooRealVar *fracBkgSBR = (RooRealVar*)bkgMassShape->createIntegral(JpsiMass,NormSet(JpsiMass),Range("SBR"));

			double FracBkgSBL = fracBkgSBL->getVal();
			double FracBkgSBR = fracBkgSBR->getVal();

			evtBkgSB[rapBin-1][ptBin-1]    = nEntries*fracBkg*(FracBkgSBL+FracBkgSBR);
			evtBkgSBErr[rapBin-1][ptBin-1] = nEntries*fracBkgErr*(FracBkgSBL+FracBkgSBR);

			//// fracSigInL(R)SB: signal fraction in L(R)SB
			RooRealVar *signalInSBL = (RooRealVar*)sigMassShape->createIntegral(JpsiMass,NormSet(JpsiMass),Range("SBL"));
			RooRealVar *signalInSBR = (RooRealVar*)sigMassShape->createIntegral(JpsiMass,NormSet(JpsiMass),Range("SBR"));
			RooRealVar *fullInSBL = (RooRealVar*)massPdf->createIntegral(JpsiMass,NormSet(JpsiMass),Range("SBL"));
			RooRealVar *fullInSBR = (RooRealVar*)massPdf->createIntegral(JpsiMass,NormSet(JpsiMass),Range("SBR"));

			double SignalInSBL = signalInSBL->getVal();
			double SignalInSBR = signalInSBR->getVal();
			double FullInSBL   = fullInSBL  ->getVal();
			double FullInSBR   = fullInSBR  ->getVal();

			//fracSigInLSB[rapBin-1][ptBin-1]    = SignalInSBL * (1-fracBkg) / FullInSBL ;
			//fracSigInRSB[rapBin-1][ptBin-1]    = SignalInSBR * (1-fracBkg) / FullInSBR ;
			//fracSigInLSBErr[rapBin-1][ptBin-1] = SignalInSBL * fracBkgErr  / FullInSBL ;
			//fracSigInRSBErr[rapBin-1][ptBin-1] = SignalInSBR *fracBkgErr   / FullInSBR ;

			evtInLSB[rapBin-1][ptBin-1] = nEntries * FracBkgSBL * fracBkg;

			fracBkgInLSB[rapBin-1][ptBin-1] = FracBkgSBL * fracBkg / FullInSBL ;
			fracBkgInRSB[rapBin-1][ptBin-1] = FracBkgSBR * fracBkg / FullInSBR ;
			fracBkgInLSBErr[rapBin-1][ptBin-1] = FracBkgSBL * fracBkgErr / FullInSBL ; 
			fracBkgInRSBErr[rapBin-1][ptBin-1] = FracBkgSBR * fracBkgErr / FullInSBR ; 

			fracSigInLSB[rapBin-1][ptBin-1]    = 1. - fracBkgInLSB[rapBin-1][ptBin-1];
			fracSigInRSB[rapBin-1][ptBin-1]    = 1. - fracBkgInRSB[rapBin-1][ptBin-1];
			fracSigInLSBErr[rapBin-1][ptBin-1] = fracBkgInLSBErr[rapBin-1][ptBin-1];
			fracSigInRSBErr[rapBin-1][ptBin-1] = fracBkgInRSBErr[rapBin-1][ptBin-1];

			//fracLSB[rapBin-1][ptBin-1]  = calcuFracL(ws, Mean, SigmaWei);
			fracLSB[rapBin-1][ptBin-1]  = 1. - (meanSR - meanSBL)/(meanSBR - meanSBL);

			mean[rapBin-1][ptBin-1] = Mean;
			sigmaWei[rapBin-1][ptBin-1] = SigmaWei*1000;
			bkgRatio3Sig[rapBin-1][ptBin-1] = BkgRatio3Sig;
			sigma1[rapBin-1][ptBin-1] = Sigma*1000;
			sigma2[rapBin-1][ptBin-1] = Sigma2*1000;

			alphaCB[rapBin-1][ptBin-1] = Alpha;
			lambdaBG[rapBin-1][ptBin-1] = -lambda;


			meanErr[rapBin-1][ptBin-1] = MeanErr;
			sigmaWeiErr[rapBin-1][ptBin-1] = SigmaWeiErr*1000;
			bkgRatio3SigErr[rapBin-1][ptBin-1] = BkgRatio3SigErr;
			sigma1Err[rapBin-1][ptBin-1] = SigmaErr*1000;
			sigma2Err[rapBin-1][ptBin-1] = Sigma2Err*1000;

			alphaCBErr[rapBin-1][ptBin-1] = AlphaErr;
			lambdaBGErr[rapBin-1][ptBin-1] = lambdaErr;

		}
	}

	for(int rapBin = 1; rapBin < RapBins+1; rapBin++){
		for(int ptBin = 1; ptBin < PtBins+1; ptBin++){
			cout<<"rap "<<rapBin<<"   pt "<<ptBin<<endl;
			cout<<"mean: "<<mean[rapBin-1][ptBin-1]<<"    meanErr: "<<meanErr[rapBin-1][ptBin-1]<<endl;
			cout<<"sigmaWei: "<<sigmaWei[rapBin-1][ptBin-1]<<"    sigmaWeiErr: "<<sigmaWeiErr[rapBin-1][ptBin-1]<<endl;
			cout<<"bkgRatio3Sig: "<<bkgRatio3Sig[rapBin-1][ptBin-1]<<"    bkgRatio3SigErr: "<<bkgRatio3SigErr[rapBin-1][ptBin-1]<<endl;
			cout<<"sigma1: "<<sigma1[rapBin-1][ptBin-1]<<"    sigma1Err: "<<sigma1Err[rapBin-1][ptBin-1]<<endl;
			cout<<"sigma2: "<<sigma2[rapBin-1][ptBin-1]<<"    sigma2Err: "<<sigma2Err[rapBin-1][ptBin-1]<<endl;
			cout<<"evtBkgSB: "<<evtBkgSB[rapBin-1][ptBin-1]<<"    evtBkgSBErr: "<<evtBkgSBErr[rapBin-1][ptBin-1]<<endl;
			cout<<"fracLSB: "<<fracLSB[rapBin-1][ptBin-1]<<"    fracLSBErr: "<<fracLSBErr[rapBin-1][ptBin-1]<<endl;
			cout<<"alphaCB: "<<alphaCB[rapBin-1][ptBin-1]<<"    alphaCBErr: "<<alphaCBErr[rapBin-1][ptBin-1]<<endl;
			cout<<"lambdaBG: "<<lambdaBG[rapBin-1][ptBin-1]<<"    lambdaBGErr: "<<lambdaBGErr[rapBin-1][ptBin-1]<<endl;
			cout<<"fracSigInLSB: "<<fracSigInLSB[rapBin-1][ptBin-1]<<"    fracSigInLSBErr: "<<fracSigInLSBErr[rapBin-1][ptBin-1]<<endl;
			cout<<"fracSigInRSB: "<<fracSigInRSB[rapBin-1][ptBin-1]<<"    fracSigInRSBErr: "<<fracSigInRSBErr[rapBin-1][ptBin-1]<<endl;
			cout<<"fracBkgInLSB: "<<fracBkgInLSB[rapBin-1][ptBin-1]<<"    fracBkgInLSBErr: "<<fracBkgInLSBErr[rapBin-1][ptBin-1]<<endl;
			cout<<"fracBkgInRSB: "<<fracBkgInRSB[rapBin-1][ptBin-1]<<"    fracBkgInRSBErr: "<<fracBkgInRSBErr[rapBin-1][ptBin-1]<<endl;
			cout<<"evtInLSB: "<<evtInLSB[rapBin-1][ptBin-1]<<"    evtInLSBErr: "<<evtInLSBErr[rapBin-1][ptBin-1]<<endl;
		}
	}

	RooRealVar JpsiMass(*ws->var("JpsiMass"));
	TGraphErrors *graph_mean[RapBins], *graph_sigmaWei[RapBins], *graph_bkgRatio3Sig[RapBins],
							 *graph_sigma1[RapBins], *graph_sigma2[RapBins], *graph_evtBkgSB[RapBins],
							 *graph_fracLSB[RapBins], *graph_alphaCB[RapBins], *graph_lambdaBG[RapBins], 
							 *graph_fracSigInLSB[RapBins], *graph_fracSigInRSB[RapBins], *graph_evtInLSB[RapBins],
							 *graph_fracBkgInLSB[RapBins], *graph_fracBkgInRSB[RapBins];

	for(int rapBin = 1; rapBin < RapBins+1; rapBin++){
		graph_mean[rapBin-1] = new TGraphErrors(PtBins,
				pT[rapBin-1], mean[rapBin-1], pTErr[rapBin-1], meanErr[rapBin-1]);
		graph_sigmaWei[rapBin-1] = new TGraphErrors(PtBins,
				pT[rapBin-1], sigmaWei[rapBin-1], pTErr[rapBin-1], sigmaWeiErr[rapBin-1]);
		graph_bkgRatio3Sig[rapBin-1] = new TGraphErrors(PtBins,
				pT[rapBin-1], bkgRatio3Sig[rapBin-1], pTErr[rapBin-1], bkgRatio3SigErr[rapBin-1]);
		graph_sigma1[rapBin-1] = new TGraphErrors(PtBins,
				pT[rapBin-1], sigma1[rapBin-1], pTErr[rapBin-1], sigma1Err[rapBin-1]);
		graph_sigma2[rapBin-1] = new TGraphErrors(PtBins,
				pT[rapBin-1], sigma2[rapBin-1], pTErr[rapBin-1], sigma2Err[rapBin-1]);
		graph_evtBkgSB[rapBin-1] = new TGraphErrors(PtBins,
				pT[rapBin-1], evtBkgSB[rapBin-1], pTErr[rapBin-1], evtBkgSBErr[rapBin-1]);
		graph_fracLSB[rapBin-1] = new TGraphErrors(PtBins,
				pT[rapBin-1], fracLSB[rapBin-1], pTErr[rapBin-1], fracLSBErr[rapBin-1]);
		graph_alphaCB[rapBin-1] = new TGraphErrors(PtBins,
				pT[rapBin-1], alphaCB[rapBin-1], pTErr[rapBin-1], alphaCBErr[rapBin-1]);
		graph_lambdaBG[rapBin-1] = new TGraphErrors(PtBins,
				pT[rapBin-1], lambdaBG[rapBin-1], pTErr[rapBin-1], lambdaBGErr[rapBin-1]);
		graph_fracSigInLSB[rapBin-1] = new TGraphErrors(PtBins,
				pT[rapBin-1], fracSigInLSB[rapBin-1], pTErr[rapBin-1], fracSigInLSBErr[rapBin-1]);
		graph_fracSigInRSB[rapBin-1] = new TGraphErrors(PtBins,
				pT[rapBin-1], fracSigInRSB[rapBin-1], pTErr[rapBin-1], fracSigInRSBErr[rapBin-1]);
		graph_evtInLSB[rapBin-1] = new TGraphErrors(PtBins,
				pT[rapBin-1], evtInLSB[rapBin-1], pTErr[rapBin-1], evtInLSBErr[rapBin-1]);
		graph_fracBkgInLSB[rapBin-1] = new TGraphErrors(PtBins,
				pT[rapBin-1], fracBkgInLSB[rapBin-1], pTErr[rapBin-1], fracBkgInLSBErr[rapBin-1]);
		graph_fracBkgInRSB[rapBin-1] = new TGraphErrors(PtBins,
				pT[rapBin-1], fracBkgInRSB[rapBin-1], pTErr[rapBin-1], fracBkgInRSBErr[rapBin-1]);

		//remove first pT bin for 2S
		if(nState==5){
			//graph_mean[rapBin-1]         -> RemovePoint(0);
			//graph_sigmaWei[rapBin-1]     -> RemovePoint(0);
			//graph_bkgRatio3Sig[rapBin-1] -> RemovePoint(0);
			//graph_sigma1[rapBin-1]       -> RemovePoint(0);
			//graph_sigma2[rapBin-1]       -> RemovePoint(0);
			//graph_evtBkgSB[rapBin-1]     -> RemovePoint(0);
			//graph_fracLSB[rapBin-1]      -> RemovePoint(0);
			//graph_alphaCB[rapBin-1]      -> RemovePoint(0);
			//graph_fracSigInLSB[rapBin-1] -> RemovePoint(0);
			//graph_fracSigInRSB[rapBin-1] -> RemovePoint(0);
			//graph_evtInLSB[rapBin-1]     -> RemovePoint(0);
			//graph_fracBkgInLSB[rapBin-1] -> RemovePoint(0);
			//graph_fracBkgInRSB[rapBin-1] -> RemovePoint(0);
		}
	}

	//double blX = 0.17-0.05, blY = 0.70+0.05, trX = 0.4-0.05, trY = 0.84+0.05;
	//double blX = 0.17, blY = 0.8, trX = 0.4, trY = 0.94;
	double blX = 0.12, blY = 0.8, trX = 0.4, trY = 0.96;
	TLegend* legend=new TLegend(blX,blY,trX,trY);
	legend->SetFillColor(kWhite);
	legend->SetTextFont(42);
	legend->SetTextSize(legendsize);
	legend->SetBorderSize(0.);
	legend->AddEntry(graph_mean[0],"|y| < 0.6","lp");
	legend->AddEntry(graph_mean[1],"0.6 < |y| < 1.2","lp");
	if(nState==5)
		legend->AddEntry(graph_mean[2],"1.2 < |y| < 1.5","lp");

	TLegend* legendFull=new TLegend(blX,blY,trX,trY);
	legendFull->SetFillColor(kWhite);
	legendFull->SetTextFont(42);
	legendFull->SetTextSize(legendsize);
	legendFull->SetBorderSize(0.);
	legendFull->AddEntry(graph_sigma1[0],"|y| < 0.6 #sigma_{1}","lp");
	legendFull->AddEntry(graph_sigma2[0],"|y| < 0.6 #sigma_{2}","lp");
	legendFull->AddEntry(graph_sigma1[1],"0.6 < |y| < 1.2 #sigma_{1}","lp");
	legendFull->AddEntry(graph_sigma2[1],"0.6 < |y| < 1.2 #sigma_{2}","lp");
	if(nState==5){
		legendFull->AddEntry(graph_sigma1[2],"1.2 < |y| < 1.5 #sigma_{1}","lp");
		legendFull->AddEntry(graph_sigma2[2],"1.2 < |y| < 1.5 #sigma_{2}","lp");
	}

	//0.75,0.75,0.82,0.84
	//double blX = 0.65, blY = 0.15, trX = 0.85, trY = 0.35;
	//blX = 0.625+0.05; blY = 0.70+0.05; trX = 0.825+0.05; trY = 0.84+0.05;
	//blX = 0.725; blY = 0.8; trX = 0.925; trY = 0.94;
	blX = 0.725; blY = 0.8; trX = 0.96; trY = 0.96;
	TLegend* legendRight=new TLegend(blX,blY,trX,trY);
	legendRight->SetFillColor(kWhite);
	legendRight->SetTextFont(42);
	legendRight->SetTextSize(legendsize);
	legendRight->SetBorderSize(0.);
	legendRight->AddEntry(graph_mean[0],"|y| < 0.6","lp");
	legendRight->AddEntry(graph_mean[1],"0.6 < |y| < 1.2","lp");
	if(nState==5)
		legendRight->AddEntry(graph_mean[2],"1.2 < |y| < 1.5","lp");

	//blX = 0.625+0.05; blY = 0.60+0.05; trX = 0.825+0.05; trY = 0.84+0.05;
	//blX = 0.725; blY = 0.7; trX = 0.925; trY = 0.94;
	blX = 0.725; blY = 0.7; trX = 0.96; trY = 0.96;
	TLegend* legendRightFull=new TLegend(blX,blY,trX,trY);
	legendRightFull->SetFillColor(kWhite);
	legendRightFull->SetTextFont(42);
	legendRightFull->SetTextSize(legendsize);
	legendRightFull->SetBorderSize(0.);
	legendRightFull->AddEntry(graph_sigma1[0],"|y| < 0.6 #sigma_{1}","lp");
	legendRightFull->AddEntry(graph_sigma2[0],"|y| < 0.6 #sigma_{2}","lp");
	legendRightFull->AddEntry(graph_sigma1[1],"0.6 < |y| < 1.2 #sigma_{1}","lp");
	legendRightFull->AddEntry(graph_sigma2[1],"0.6 < |y| < 1.2 #sigma_{2}","lp");
	if(nState==5){
		legendRightFull->AddEntry(graph_sigma1[2],"1.2 < |y| < 1.5 #sigma_{1}","lp");
		legendRightFull->AddEntry(graph_sigma2[2],"1.2 < |y| < 1.5 #sigma_{2}","lp");
	}

	TLegend* legendSigConRightFull=new TLegend(blX,blY,trX,trY);
	legendSigConRightFull->SetFillColor(kWhite);
	legendSigConRightFull->SetTextFont(42);
	legendSigConRightFull->SetTextSize(legendsize);
	legendSigConRightFull->SetBorderSize(0.);
	legendSigConRightFull->AddEntry(graph_fracSigInLSB[0],"|y| < 0.6 LSB","lp");
	legendSigConRightFull->AddEntry(graph_fracSigInRSB[0],"|y| < 0.6 RSB","lp");
	legendSigConRightFull->AddEntry(graph_fracSigInLSB[1],"0.6 < |y| < 1.2 LSB","lp");
	legendSigConRightFull->AddEntry(graph_fracSigInRSB[1],"0.6 < |y| < 1.2 RSB","lp");
	if(nState==5){
		legendSigConRightFull->AddEntry(graph_fracSigInLSB[2],"1.2 < |y| < 1.5 LSB","lp");
		legendSigConRightFull->AddEntry(graph_fracSigInRSB[2],"1.2 < |y| < 1.5 RSB","lp");
	}

	gStyle->SetPadBottomMargin(0.11); //0.12
	gStyle->SetPadLeftMargin(0.08); //0.12
	gStyle->SetPadRightMargin(0.02); //0.05
	gStyle->SetPadTopMargin(0.02); //0.05

	TCanvas *c1=new TCanvas("c1","");
	c1->SetTickx();
	c1->SetTicky();
	//c1->SetGridx();
	//c1->SetGridy();
	gStyle->SetTitleFillColor(10);
	gStyle->SetTitleBorderSize(1);
	gStyle->SetOptFit(1);            
	gStyle->SetOptStat(1);           
	gStyle->SetTitleFont(22);        
	gStyle->SetStatFont(22);         
	gStyle->SetStatColor(10);        
	gStyle->SetStatBorderSize(1);
	gStyle->SetLabelFont(22,"X");    
	gStyle->SetLabelFont(22,"Y");    
	gStyle->SetTitleXOffset(1.2);    
	gStyle->SetTitleYOffset(1.2);    
	gStyle->SetHistLineWidth(2);     
	gStyle->SetStatX(0.9);           
	gStyle->SetStatY(0.9);           
	gStyle->SetTitleX(0.15);         
	gStyle->SetTitleY(0.96);         

	double Xmin = 0.,  Xmax = 0.;
	double Ymin = 0.2, Ymax = 0.75;

	Xmin = 0.;    Xmax = 100.;
	Ymin = 3.085; Ymax = 3.095;
	if(nState==5) {Ymin = 3.674; Ymax = 3.685;}

	graph_mean[0]->SetTitle("");
	graph_mean[0]->GetXaxis()->SetTitle("p_{T} (GeV)");
	graph_mean[0]->GetYaxis()->SetTitle("mean (GeV)");
	graph_mean[0]->GetXaxis()->SetLimits(Xmin, Xmax);
	graph_mean[0]->GetYaxis()->SetRangeUser(Ymin, Ymax);
	graph_mean[0]->SetMarkerStyle(onia::marker_rapForPTBins[2]);
	graph_mean[0]->SetMarkerColor(onia::colour_rapForPTBins[2]);
	graph_mean[0]->SetLineColor(onia::colour_rapForPTBins[2]);
	graph_mean[1]->SetTitle("");
	graph_mean[1]->GetXaxis()->SetTitle("p_{T} (GeV)");
	graph_mean[1]->GetYaxis()->SetTitle("mean (GeV)");
	graph_mean[1]->GetXaxis()->SetLimits(Xmin, Xmax);
	graph_mean[1]->GetYaxis()->SetRangeUser(Ymin, Ymax);
	graph_mean[1]->SetMarkerStyle(onia::marker_rapForPTBins[3]);
	graph_mean[1]->SetMarkerColor(onia::colour_rapForPTBins[3]);
	graph_mean[1]->SetLineColor(onia::colour_rapForPTBins[3]);
	if(nState==5){
		graph_mean[2]->SetTitle("");
		graph_mean[2]->GetXaxis()->SetTitle("p_{T} (GeV)");
		graph_mean[2]->GetYaxis()->SetTitle("mean (GeV)");
		graph_mean[2]->GetXaxis()->SetLimits(Xmin, Xmax);
		graph_mean[2]->GetYaxis()->SetRangeUser(Ymin, Ymax);
		graph_mean[2]->SetMarkerStyle(onia::marker_rapForPTBins[4]);
		graph_mean[2]->SetMarkerColor(onia::colour_rapForPTBins[4]);
		graph_mean[2]->SetLineColor(onia::colour_rapForPTBins[4]);
	}

	Ymin = 0; Ymax = 90;

	graph_sigmaWei[0]->SetTitle("");
	graph_sigmaWei[0]->GetXaxis()->SetTitle("p_{T} (GeV)");
	graph_sigmaWei[0]->GetYaxis()->SetTitle("effective #sigma (MeV)");
	graph_sigmaWei[0]->GetXaxis()->SetLimits(Xmin, Xmax);
	graph_sigmaWei[0]->GetYaxis()->SetRangeUser(Ymin, Ymax);
	graph_sigmaWei[0]->SetMarkerStyle(onia::marker_rapForPTBins[2]);
	graph_sigmaWei[0]->SetMarkerColor(onia::colour_rapForPTBins[2]);
	graph_sigmaWei[0]->SetLineColor(onia::colour_rapForPTBins[2]);
	graph_sigmaWei[1]->SetTitle("");
	graph_sigmaWei[1]->GetXaxis()->SetTitle("p_{T} (GeV)");
	graph_sigmaWei[1]->GetYaxis()->SetTitle("effective #sigma (MeV)");
	graph_sigmaWei[1]->GetXaxis()->SetLimits(Xmin, Xmax);
	graph_sigmaWei[1]->GetYaxis()->SetRangeUser(Ymin, Ymax);
	graph_sigmaWei[1]->SetMarkerStyle(onia::marker_rapForPTBins[3]);
	graph_sigmaWei[1]->SetMarkerColor(onia::colour_rapForPTBins[3]);
	graph_sigmaWei[1]->SetLineColor(onia::colour_rapForPTBins[3]);
	if(nState==5){
		graph_sigmaWei[2]->SetTitle("");
		graph_sigmaWei[2]->GetXaxis()->SetTitle("p_{T} (GeV)");
		graph_sigmaWei[2]->GetYaxis()->SetTitle("effective #sigma (MeV)");
		graph_sigmaWei[2]->GetXaxis()->SetLimits(Xmin, Xmax);
		graph_sigmaWei[2]->GetYaxis()->SetRangeUser(Ymin, Ymax);
		graph_sigmaWei[2]->SetMarkerStyle(onia::marker_rapForPTBins[4]);
		graph_sigmaWei[2]->SetMarkerColor(onia::colour_rapForPTBins[4]);
		graph_sigmaWei[2]->SetLineColor(onia::colour_rapForPTBins[4]);
	}

	graph_sigma1[0]->SetTitle("");
	graph_sigma1[0]->GetXaxis()->SetTitle("p_{T} (GeV)");
	graph_sigma1[0]->GetYaxis()->SetTitle("#sigma_{1} (MeV)");
	graph_sigma1[0]->GetXaxis()->SetLimits(Xmin, Xmax);
	graph_sigma1[0]->GetYaxis()->SetRangeUser(Ymin, Ymax);
	graph_sigma1[0]->SetMarkerStyle(onia::marker_rapForPTBins[2]);
	graph_sigma1[0]->SetMarkerColor(onia::colour_rapForPTBins[2]);
	graph_sigma1[0]->SetLineColor(onia::colour_rapForPTBins[2]);
	graph_sigma1[1]->SetTitle("");
	graph_sigma1[1]->GetXaxis()->SetTitle("p_{T} (GeV)");
	graph_sigma1[1]->GetYaxis()->SetTitle("#sigma_{1} (MeV)");
	graph_sigma1[1]->GetXaxis()->SetLimits(Xmin, Xmax);
	graph_sigma1[1]->GetYaxis()->SetRangeUser(Ymin, Ymax);
	graph_sigma1[1]->SetMarkerStyle(onia::marker_rapForPTBins[3]);
	graph_sigma1[1]->SetMarkerColor(onia::colour_rapForPTBins[3]);
	graph_sigma1[1]->SetLineColor(onia::colour_rapForPTBins[3]);
	if(nState==5){
		graph_sigma1[2]->SetTitle("");
		graph_sigma1[2]->GetXaxis()->SetTitle("p_{T} (GeV)");
		graph_sigma1[2]->GetYaxis()->SetTitle("#sigma_{1} (MeV)");
		graph_sigma1[2]->GetXaxis()->SetLimits(Xmin, Xmax);
		graph_sigma1[2]->GetYaxis()->SetRangeUser(Ymin, Ymax);
		graph_sigma1[2]->SetMarkerStyle(onia::marker_rapForPTBins[4]);
		graph_sigma1[2]->SetMarkerColor(onia::colour_rapForPTBins[4]);
		graph_sigma1[2]->SetLineColor(onia::colour_rapForPTBins[4]);
	}

	Ymin = 0; Ymax = 90;

	graph_sigma2[0]->SetTitle("");
	graph_sigma2[0]->GetXaxis()->SetTitle("p_{T} (GeV)");
	graph_sigma2[0]->GetYaxis()->SetTitle("#sigma_{2} (MeV)");
	graph_sigma2[0]->GetXaxis()->SetLimits(Xmin, Xmax);
	graph_sigma2[0]->GetYaxis()->SetRangeUser(Ymin, Ymax);
	graph_sigma2[0]->SetMarkerStyle(onia::marker_rapForPTBins[2]);
	graph_sigma2[0]->SetMarkerColor(onia::colour_rapForPTBins[2]);
	graph_sigma2[0]->SetLineColor(onia::colour_rapForPTBins[2]);
	graph_sigma2[1]->SetTitle("");
	graph_sigma2[1]->GetXaxis()->SetTitle("p_{T} (GeV)");
	graph_sigma2[1]->GetYaxis()->SetTitle("#sigma_{2} (MeV)");
	graph_sigma2[1]->GetXaxis()->SetLimits(Xmin, Xmax);
	graph_sigma2[1]->GetYaxis()->SetRangeUser(Ymin, Ymax);
	graph_sigma2[1]->SetMarkerStyle(onia::marker_rapForPTBins[3]);
	graph_sigma2[1]->SetMarkerColor(onia::colour_rapForPTBins[3]);
	graph_sigma2[1]->SetLineColor(onia::colour_rapForPTBins[3]);
	if(nState==5){
		graph_sigma2[2]->SetTitle("");
		graph_sigma2[2]->GetXaxis()->SetTitle("p_{T} (GeV)");
		graph_sigma2[2]->GetYaxis()->SetTitle("#sigma_{2} (MeV)");
		graph_sigma2[2]->GetXaxis()->SetLimits(Xmin, Xmax);
		graph_sigma2[2]->GetYaxis()->SetRangeUser(Ymin, Ymax);
		graph_sigma2[2]->SetMarkerStyle(onia::marker_rapForPTBins[4]);
		graph_sigma2[2]->SetMarkerColor(onia::colour_rapForPTBins[4]);
		graph_sigma2[2]->SetLineColor(onia::colour_rapForPTBins[4]);
	}

	Ymin = 0.02; Ymax = 0.09;
	if(nState==5){Ymin = 0.1; Ymax = 0.6;}

	graph_bkgRatio3Sig[0]->SetTitle("");
	graph_bkgRatio3Sig[0]->GetXaxis()->SetTitle("p_{T} (GeV)");
	graph_bkgRatio3Sig[0]->GetYaxis()->SetTitle("Bg fraction(#pm3#sigma)");
	graph_bkgRatio3Sig[0]->GetXaxis()->SetLimits(Xmin, Xmax);
	graph_bkgRatio3Sig[0]->GetYaxis()->SetRangeUser(Ymin, Ymax);
	graph_bkgRatio3Sig[0]->SetMarkerStyle(onia::marker_rapForPTBins[2]);
	graph_bkgRatio3Sig[0]->SetMarkerColor(onia::colour_rapForPTBins[2]);
	graph_bkgRatio3Sig[0]->SetLineColor(onia::colour_rapForPTBins[2]);
	graph_bkgRatio3Sig[1]->SetTitle("");
	graph_bkgRatio3Sig[1]->GetXaxis()->SetTitle("p_{T} (GeV)");
	graph_bkgRatio3Sig[1]->GetYaxis()->SetTitle("Bg fraction(#pm3#sigma)");
	graph_bkgRatio3Sig[1]->GetXaxis()->SetLimits(Xmin, Xmax);
	graph_bkgRatio3Sig[1]->GetYaxis()->SetRangeUser(Ymin, Ymax);
	graph_bkgRatio3Sig[1]->SetMarkerStyle(onia::marker_rapForPTBins[3]);
	graph_bkgRatio3Sig[1]->SetMarkerColor(onia::colour_rapForPTBins[3]);
	graph_bkgRatio3Sig[1]->SetLineColor(onia::colour_rapForPTBins[3]);
	if(nState==5){
		graph_bkgRatio3Sig[2]->SetTitle("");
		graph_bkgRatio3Sig[2]->GetXaxis()->SetTitle("p_{T} (GeV)");
		graph_bkgRatio3Sig[2]->GetYaxis()->SetTitle("Bg fraction(#pm3#sigma)");
		graph_bkgRatio3Sig[2]->GetXaxis()->SetLimits(Xmin, Xmax);
		graph_bkgRatio3Sig[2]->GetYaxis()->SetRangeUser(Ymin, Ymax);
		graph_bkgRatio3Sig[2]->SetMarkerStyle(onia::marker_rapForPTBins[4]);
		graph_bkgRatio3Sig[2]->SetMarkerColor(onia::colour_rapForPTBins[4]);
		graph_bkgRatio3Sig[2]->SetLineColor(onia::colour_rapForPTBins[4]);
	}

	double x=0,y=0;
	graph_evtBkgSB[0]->GetPoint(0,x,y);
	//Ymin = 100; Ymax = y*3;
	//if(nState==5){Ymin = 1000; Ymax = y*10;}
	Ymin = 8; Ymax = 1.e6+100;
	if(nState==5){Ymin = 8; Ymax = 1.e6+100;}

	graph_evtBkgSB[0]->SetTitle("");
	graph_evtBkgSB[0]->GetXaxis()->SetTitle("p_{T} (GeV)");
	graph_evtBkgSB[0]->GetYaxis()->SetTitle("number of events in SB");
	graph_evtBkgSB[0]->GetXaxis()->SetLimits(Xmin, Xmax);
	graph_evtBkgSB[0]->GetYaxis()->SetRangeUser(Ymin, Ymax);
	graph_evtBkgSB[0]->SetMarkerStyle(onia::marker_rapForPTBins[2]);
	graph_evtBkgSB[0]->SetMarkerColor(onia::colour_rapForPTBins[2]);
	graph_evtBkgSB[0]->SetLineColor(onia::colour_rapForPTBins[2]);
	graph_evtBkgSB[1]->SetTitle("");
	graph_evtBkgSB[1]->GetXaxis()->SetTitle("p_{T} (GeV)");
	graph_evtBkgSB[1]->GetYaxis()->SetTitle("number of events in SB");
	graph_evtBkgSB[1]->GetXaxis()->SetLimits(Xmin, Xmax);
	graph_evtBkgSB[1]->GetYaxis()->SetRangeUser(Ymin, Ymax);
	graph_evtBkgSB[1]->SetMarkerStyle(onia::marker_rapForPTBins[3]);
	graph_evtBkgSB[1]->SetMarkerColor(onia::colour_rapForPTBins[3]);
	graph_evtBkgSB[1]->SetLineColor(onia::colour_rapForPTBins[3]);
	if(nState==5){
		graph_evtBkgSB[2]->SetTitle("");
		graph_evtBkgSB[2]->GetXaxis()->SetTitle("p_{T} (GeV)");
		graph_evtBkgSB[2]->GetYaxis()->SetTitle("number of events in SB");
		graph_evtBkgSB[2]->GetXaxis()->SetLimits(Xmin, Xmax);
		graph_evtBkgSB[2]->GetYaxis()->SetRangeUser(Ymin, Ymax);
		graph_evtBkgSB[2]->SetMarkerStyle(onia::marker_rapForPTBins[4]);
		graph_evtBkgSB[2]->SetMarkerColor(onia::colour_rapForPTBins[4]);
		graph_evtBkgSB[2]->SetLineColor(onia::colour_rapForPTBins[4]);
	}

	Ymin = 0.; Ymax = 1.;
	graph_fracLSB[0]->SetTitle("");
	graph_fracLSB[0]->GetXaxis()->SetTitle("p_{T} (GeV)");
	graph_fracLSB[0]->GetYaxis()->SetTitle("frac_{LSB}");
	graph_fracLSB[0]->GetXaxis()->SetLimits(Xmin, Xmax);
	graph_fracLSB[0]->GetYaxis()->SetRangeUser(Ymin, Ymax);
	graph_fracLSB[0]->SetMarkerStyle(onia::marker_rapForPTBins[2]);
	graph_fracLSB[0]->SetMarkerColor(onia::colour_rapForPTBins[2]);
	graph_fracLSB[0]->SetLineColor(onia::colour_rapForPTBins[2]);
	graph_fracLSB[1]->SetTitle("");
	graph_fracLSB[1]->GetXaxis()->SetTitle("p_{T} (GeV)");
	graph_fracLSB[1]->GetYaxis()->SetTitle("frac_{LSB}");
	graph_fracLSB[1]->GetXaxis()->SetLimits(Xmin, Xmax);
	graph_fracLSB[1]->GetYaxis()->SetRangeUser(Ymin, Ymax);
	graph_fracLSB[1]->SetMarkerStyle(onia::marker_rapForPTBins[3]);
	graph_fracLSB[1]->SetMarkerColor(onia::colour_rapForPTBins[3]);
	graph_fracLSB[1]->SetLineColor(onia::colour_rapForPTBins[3]);
	if(nState==5){
		graph_fracLSB[2]->SetTitle("");
		graph_fracLSB[2]->GetXaxis()->SetTitle("p_{T} (GeV)");
		graph_fracLSB[2]->GetYaxis()->SetTitle("frac_{LSB}");
		graph_fracLSB[2]->GetXaxis()->SetLimits(Xmin, Xmax);
		graph_fracLSB[2]->GetYaxis()->SetRangeUser(Ymin, Ymax);
		graph_fracLSB[2]->SetMarkerStyle(onia::marker_rapForPTBins[4]);
		graph_fracLSB[2]->SetMarkerColor(onia::colour_rapForPTBins[4]);
		graph_fracLSB[2]->SetLineColor(onia::colour_rapForPTBins[4]);
	}

	Ymin = 0.5; Ymax = 4.;
	graph_alphaCB[0]->SetTitle("");
	graph_alphaCB[0]->GetXaxis()->SetTitle("p_{T} (GeV)");
	graph_alphaCB[0]->GetYaxis()->SetTitle("#alpha_{CB}");
	graph_alphaCB[0]->GetXaxis()->SetLimits(Xmin, Xmax);
	graph_alphaCB[0]->GetYaxis()->SetRangeUser(Ymin, Ymax);
	graph_alphaCB[0]->SetMarkerStyle(onia::marker_rapForPTBins[2]);
	graph_alphaCB[0]->SetMarkerColor(onia::colour_rapForPTBins[2]);
	graph_alphaCB[0]->SetLineColor(onia::colour_rapForPTBins[2]);
	graph_alphaCB[1]->SetTitle("");
	graph_alphaCB[1]->GetXaxis()->SetTitle("p_{T} (GeV)");
	graph_alphaCB[1]->GetYaxis()->SetTitle("#alpha_{CB}");
	graph_alphaCB[1]->GetXaxis()->SetLimits(Xmin, Xmax);
	graph_alphaCB[1]->GetYaxis()->SetRangeUser(Ymin, Ymax);
	graph_alphaCB[1]->SetMarkerStyle(onia::marker_rapForPTBins[3]);
	graph_alphaCB[1]->SetMarkerColor(onia::colour_rapForPTBins[3]);
	graph_alphaCB[1]->SetLineColor(onia::colour_rapForPTBins[3]);
	if(nState==5){
		graph_alphaCB[2]->SetTitle("");
		graph_alphaCB[2]->GetXaxis()->SetTitle("p_{T} (GeV)");
		graph_alphaCB[2]->GetYaxis()->SetTitle("#alpha_{CB}");
		graph_alphaCB[2]->GetXaxis()->SetLimits(Xmin, Xmax);
		graph_alphaCB[2]->GetYaxis()->SetRangeUser(Ymin, Ymax);
		graph_alphaCB[2]->SetMarkerStyle(onia::marker_rapForPTBins[4]);
		graph_alphaCB[2]->SetMarkerColor(onia::colour_rapForPTBins[4]);
		graph_alphaCB[2]->SetLineColor(onia::colour_rapForPTBins[4]);
	}

	Ymin = 0.; Ymax = 4.;
	graph_lambdaBG[0]->SetTitle("");
	graph_lambdaBG[0]->GetXaxis()->SetTitle("p_{T} (GeV)");
	graph_lambdaBG[0]->GetYaxis()->SetTitle("#lambda_{BG}");
	graph_lambdaBG[0]->GetXaxis()->SetLimits(Xmin, Xmax);
	graph_lambdaBG[0]->GetYaxis()->SetRangeUser(Ymin, Ymax);
	graph_lambdaBG[0]->SetMarkerStyle(onia::marker_rapForPTBins[2]);
	graph_lambdaBG[0]->SetMarkerColor(onia::colour_rapForPTBins[2]);
	graph_lambdaBG[0]->SetLineColor(onia::colour_rapForPTBins[2]);
	graph_lambdaBG[1]->SetTitle("");
	graph_lambdaBG[1]->GetXaxis()->SetTitle("p_{T} (GeV)");
	graph_lambdaBG[1]->GetYaxis()->SetTitle("#lambda_{BG}");
	graph_lambdaBG[1]->GetXaxis()->SetLimits(Xmin, Xmax);
	graph_lambdaBG[1]->GetYaxis()->SetRangeUser(Ymin, Ymax);
	graph_lambdaBG[1]->SetMarkerStyle(onia::marker_rapForPTBins[3]);
	graph_lambdaBG[1]->SetMarkerColor(onia::colour_rapForPTBins[3]);
	graph_lambdaBG[1]->SetLineColor(onia::colour_rapForPTBins[3]);
	if(nState==5){
		graph_lambdaBG[2]->SetTitle("");
		graph_lambdaBG[2]->GetXaxis()->SetTitle("p_{T} (GeV)");
		graph_lambdaBG[2]->GetYaxis()->SetTitle("#lambda_{BG}");
		graph_lambdaBG[2]->GetXaxis()->SetLimits(Xmin, Xmax);
		graph_lambdaBG[2]->GetYaxis()->SetRangeUser(Ymin, Ymax);
		graph_lambdaBG[2]->SetMarkerStyle(onia::marker_rapForPTBins[4]);
		graph_lambdaBG[2]->SetMarkerColor(onia::colour_rapForPTBins[4]);
		graph_lambdaBG[2]->SetLineColor(onia::colour_rapForPTBins[4]);
	}

	Ymin = 0.0; Ymax = 0.4;
	graph_fracSigInLSB[0]->SetTitle("");
	graph_fracSigInLSB[0]->GetXaxis()->SetTitle("p_{T} (GeV)");
	graph_fracSigInLSB[0]->GetYaxis()->SetTitle("signal contamination in LSB");
	graph_fracSigInLSB[0]->GetXaxis()->SetLimits(Xmin, Xmax);
	graph_fracSigInLSB[0]->GetYaxis()->SetRangeUser(Ymin, Ymax);
	graph_fracSigInLSB[0]->SetMarkerStyle(onia::marker_rapForPTBins[2]);
	graph_fracSigInLSB[0]->SetMarkerColor(onia::colour_rapForPTBins[2]);
	graph_fracSigInLSB[0]->SetLineColor(onia::colour_rapForPTBins[2]);
	graph_fracSigInLSB[1]->SetTitle("");
	graph_fracSigInLSB[1]->GetXaxis()->SetTitle("p_{T} (GeV)");
	graph_fracSigInLSB[1]->GetYaxis()->SetTitle("signal contamination in LSB");
	graph_fracSigInLSB[1]->GetXaxis()->SetLimits(Xmin, Xmax);
	graph_fracSigInLSB[1]->GetYaxis()->SetRangeUser(Ymin, Ymax);
	graph_fracSigInLSB[1]->SetMarkerStyle(onia::marker_rapForPTBins[3]);
	graph_fracSigInLSB[1]->SetMarkerColor(onia::colour_rapForPTBins[3]);
	graph_fracSigInLSB[1]->SetLineColor(onia::colour_rapForPTBins[3]);
	if(nState==5){
		graph_fracSigInLSB[2]->SetTitle("");
		graph_fracSigInLSB[2]->GetXaxis()->SetTitle("p_{T} (GeV)");
		graph_fracSigInLSB[2]->GetYaxis()->SetTitle("signal contamination in LSB");
		graph_fracSigInLSB[2]->GetXaxis()->SetLimits(Xmin, Xmax);
		graph_fracSigInLSB[2]->GetYaxis()->SetRangeUser(Ymin, Ymax);
		graph_fracSigInLSB[2]->SetMarkerStyle(onia::marker_rapForPTBins[4]);
		graph_fracSigInLSB[2]->SetMarkerColor(onia::colour_rapForPTBins[4]);
		graph_fracSigInLSB[2]->SetLineColor(onia::colour_rapForPTBins[4]);
	}

	Ymin = 0.0; Ymax = 0.4;
	graph_fracSigInRSB[0]->SetTitle("");
	graph_fracSigInRSB[0]->GetXaxis()->SetTitle("p_{T} (GeV)");
	graph_fracSigInRSB[0]->GetYaxis()->SetTitle("signal contamination in RSB");
	graph_fracSigInRSB[0]->GetXaxis()->SetLimits(Xmin, Xmax);
	graph_fracSigInRSB[0]->GetYaxis()->SetRangeUser(Ymin, Ymax);
	graph_fracSigInRSB[0]->SetMarkerStyle(onia::marker_rapForPTBins[2]);
	graph_fracSigInRSB[0]->SetMarkerColor(onia::colour_rapForPTBins[2]);
	graph_fracSigInRSB[0]->SetLineColor(onia::colour_rapForPTBins[2]);
	graph_fracSigInRSB[1]->SetTitle("");
	graph_fracSigInRSB[1]->GetXaxis()->SetTitle("p_{T} (GeV)");
	graph_fracSigInRSB[1]->GetYaxis()->SetTitle("signal contamination in RSB");
	graph_fracSigInRSB[1]->GetXaxis()->SetLimits(Xmin, Xmax);
	graph_fracSigInRSB[1]->GetYaxis()->SetRangeUser(Ymin, Ymax);
	graph_fracSigInRSB[1]->SetMarkerStyle(onia::marker_rapForPTBins[3]);
	graph_fracSigInRSB[1]->SetMarkerColor(onia::colour_rapForPTBins[3]);
	graph_fracSigInRSB[1]->SetLineColor(onia::colour_rapForPTBins[3]);
	if(nState==5){
		graph_fracSigInRSB[2]->SetTitle("");
		graph_fracSigInRSB[2]->GetXaxis()->SetTitle("p_{T} (GeV)");
		graph_fracSigInRSB[2]->GetYaxis()->SetTitle("signal contamination in RSB");
		graph_fracSigInRSB[2]->GetXaxis()->SetLimits(Xmin, Xmax);
		graph_fracSigInRSB[2]->GetYaxis()->SetRangeUser(Ymin, Ymax);
		graph_fracSigInRSB[2]->SetMarkerStyle(onia::marker_rapForPTBins[4]);
		graph_fracSigInRSB[2]->SetMarkerColor(onia::colour_rapForPTBins[4]);
		graph_fracSigInRSB[2]->SetLineColor(onia::colour_rapForPTBins[4]);
	}

	Ymin = 8; Ymax = 1.e6+100;
	graph_evtInLSB[0]->SetTitle("");
	graph_evtInLSB[0]->GetXaxis()->SetTitle("p_{T} (GeV)");
	graph_evtInLSB[0]->GetYaxis()->SetTitle("number of events in LSB");
	graph_evtInLSB[0]->GetXaxis()->SetLimits(Xmin, Xmax);
	graph_evtInLSB[0]->GetYaxis()->SetRangeUser(Ymin, Ymax);
	graph_evtInLSB[0]->SetMarkerStyle(onia::marker_rapForPTBins[2]);
	graph_evtInLSB[0]->SetMarkerColor(onia::colour_rapForPTBins[2]);
	graph_evtInLSB[0]->SetLineColor(onia::colour_rapForPTBins[2]);
	graph_evtInLSB[1]->SetTitle("");
	graph_evtInLSB[1]->GetXaxis()->SetTitle("p_{T} (GeV)");
	graph_evtInLSB[1]->GetYaxis()->SetTitle("number of events in LSB");
	graph_evtInLSB[1]->GetXaxis()->SetLimits(Xmin, Xmax);
	graph_evtInLSB[1]->GetYaxis()->SetRangeUser(Ymin, Ymax);
	graph_evtInLSB[1]->SetMarkerStyle(onia::marker_rapForPTBins[3]);
	graph_evtInLSB[1]->SetMarkerColor(onia::colour_rapForPTBins[3]);
	graph_evtInLSB[1]->SetLineColor(onia::colour_rapForPTBins[3]);
	if(nState==5){
		graph_evtInLSB[2]->SetTitle("");
		graph_evtInLSB[2]->GetXaxis()->SetTitle("p_{T} (GeV)");
		graph_evtInLSB[2]->GetYaxis()->SetTitle("number of events in LSB");
		graph_evtInLSB[2]->GetXaxis()->SetLimits(Xmin, Xmax);
		graph_evtInLSB[2]->GetYaxis()->SetRangeUser(Ymin, Ymax);
		graph_evtInLSB[2]->SetMarkerStyle(onia::marker_rapForPTBins[4]);
		graph_evtInLSB[2]->SetMarkerColor(onia::colour_rapForPTBins[4]);
		graph_evtInLSB[2]->SetLineColor(onia::colour_rapForPTBins[4]);
	}

	Ymin = 0.5; Ymax = 1.2;
	graph_fracBkgInLSB[0]->SetTitle("");
	graph_fracBkgInLSB[0]->GetXaxis()->SetTitle("p_{T} (GeV)");
	graph_fracBkgInLSB[0]->GetYaxis()->SetTitle("Bg fraction in LSB");
	graph_fracBkgInLSB[0]->GetXaxis()->SetLimits(Xmin, Xmax);
	graph_fracBkgInLSB[0]->GetYaxis()->SetRangeUser(Ymin, Ymax);
	graph_fracBkgInLSB[0]->SetMarkerStyle(onia::marker_rapForPTBins[2]);
	graph_fracBkgInLSB[0]->SetMarkerColor(onia::colour_rapForPTBins[2]);
	graph_fracBkgInLSB[0]->SetLineColor(onia::colour_rapForPTBins[2]);
	graph_fracBkgInLSB[1]->SetTitle("");
	graph_fracBkgInLSB[1]->GetXaxis()->SetTitle("p_{T} (GeV)");
	graph_fracBkgInLSB[1]->GetYaxis()->SetTitle("Bg fraction in LSB");
	graph_fracBkgInLSB[1]->GetXaxis()->SetLimits(Xmin, Xmax);
	graph_fracBkgInLSB[1]->GetYaxis()->SetRangeUser(Ymin, Ymax);
	graph_fracBkgInLSB[1]->SetMarkerStyle(onia::marker_rapForPTBins[3]);
	graph_fracBkgInLSB[1]->SetMarkerColor(onia::colour_rapForPTBins[3]);
	graph_fracBkgInLSB[1]->SetLineColor(onia::colour_rapForPTBins[3]);
	if(nState==5){
		graph_fracBkgInLSB[2]->SetTitle("");
		graph_fracBkgInLSB[2]->GetXaxis()->SetTitle("p_{T} (GeV)");
		graph_fracBkgInLSB[2]->GetYaxis()->SetTitle("Bg fraction in LSB");
		graph_fracBkgInLSB[2]->GetXaxis()->SetLimits(Xmin, Xmax);
		graph_fracBkgInLSB[2]->GetYaxis()->SetRangeUser(Ymin, Ymax);
		graph_fracBkgInLSB[2]->SetMarkerStyle(onia::marker_rapForPTBins[4]);
		graph_fracBkgInLSB[2]->SetMarkerColor(onia::colour_rapForPTBins[4]);
		graph_fracBkgInLSB[2]->SetLineColor(onia::colour_rapForPTBins[4]);
	}

	Ymin = 0.5; Ymax = 1.2;
	graph_fracBkgInRSB[0]->SetTitle("");
	graph_fracBkgInRSB[0]->GetXaxis()->SetTitle("p_{T} (GeV)");
	graph_fracBkgInRSB[0]->GetYaxis()->SetTitle("Bg fraction in RSB");
	graph_fracBkgInRSB[0]->GetXaxis()->SetLimits(Xmin, Xmax);
	graph_fracBkgInRSB[0]->GetYaxis()->SetRangeUser(Ymin, Ymax);
	graph_fracBkgInRSB[0]->SetMarkerStyle(onia::marker_rapForPTBins[2]);
	graph_fracBkgInRSB[0]->SetMarkerColor(onia::colour_rapForPTBins[2]);
	graph_fracBkgInRSB[0]->SetLineColor(onia::colour_rapForPTBins[2]);
	graph_fracBkgInRSB[1]->SetTitle("");
	graph_fracBkgInRSB[1]->GetXaxis()->SetTitle("p_{T} (GeV)");
	graph_fracBkgInRSB[1]->GetYaxis()->SetTitle("Bg fraction in RSB");
	graph_fracBkgInRSB[1]->GetXaxis()->SetLimits(Xmin, Xmax);
	graph_fracBkgInRSB[1]->GetYaxis()->SetRangeUser(Ymin, Ymax);
	graph_fracBkgInRSB[1]->SetMarkerStyle(onia::marker_rapForPTBins[3]);
	graph_fracBkgInRSB[1]->SetMarkerColor(onia::colour_rapForPTBins[3]);
	graph_fracBkgInRSB[1]->SetLineColor(onia::colour_rapForPTBins[3]);
	if(nState==5){
		graph_fracBkgInRSB[2]->SetTitle("");
		graph_fracBkgInRSB[2]->GetXaxis()->SetTitle("p_{T} (GeV)");
		graph_fracBkgInRSB[2]->GetYaxis()->SetTitle("Bg fraction in RSB");
		graph_fracBkgInRSB[2]->GetXaxis()->SetLimits(Xmin, Xmax);
		graph_fracBkgInRSB[2]->GetYaxis()->SetRangeUser(Ymin, Ymax);
		graph_fracBkgInRSB[2]->SetMarkerStyle(onia::marker_rapForPTBins[4]);
		graph_fracBkgInRSB[2]->SetMarkerColor(onia::colour_rapForPTBins[4]);
		graph_fracBkgInRSB[2]->SetLineColor(onia::colour_rapForPTBins[4]);
	}

	//double left=0.43, top=0.8+0.05, textSize=0.055;
	//if(nState==5) left=0.41;
	double left=0.45, top=0.87, textSize=0.055;
	if(nState==5) left=0.43;
	TLatex *latex=new TLatex();
	latex->SetTextFont(42);
	latex->SetNDC(kTRUE);
	latex->SetTextSize(textSize);
	double step=textSize*1.3;
	///////////////////
	graph_mean[0]->Draw("AP");
	graph_mean[1]->Draw("P");
	if(nState==5)
		graph_mean[2]->Draw("P");
	legendRight->Draw();
	if(nState==4)
		latex->DrawLatex(left,top,"J/#psi");
	if(nState==5)
		latex->DrawLatex(left,top,"#psi(2S)");
	c1->SaveAs(Form("%s/mean.pdf",savePath.str().c_str()));

	graph_sigmaWei[0]->Draw("AP");
	graph_sigmaWei[1]->Draw("P");
	if(nState==5){
		graph_sigmaWei[2]->Draw("P");
		legendRight->Draw();
		latex->DrawLatex(left,top,"#psi(2S)");
	}
	else{
		legend->Draw();
		latex->DrawLatex(left,top,"J/#psi");
	}
	c1->SaveAs(Form("%s/sigmaEffective.pdf",savePath.str().c_str()));

	graph_sigma1[0]->Draw("AP");
	graph_sigma1[1]->Draw("P");
	if(nState==5)
		graph_sigma1[2]->Draw("P");
	legend->Draw();
	c1->SaveAs(Form("%s/sigma1.pdf",savePath.str().c_str()));

	graph_sigma2[0]->Draw("AP");
	graph_sigma2[1]->Draw("P");
	if(nState==5)
		graph_sigma2[2]->Draw("P");
	legend->Draw();
	c1->SaveAs(Form("%s/sigma2.pdf",savePath.str().c_str()));

	graph_bkgRatio3Sig[0]->Draw("AP");
	graph_bkgRatio3Sig[1]->Draw("P");
	if(nState==5)
		graph_bkgRatio3Sig[2]->Draw("P");
	legend->Draw();
	c1->SaveAs(Form("%s/bkgRatio3Sig.pdf",savePath.str().c_str()));

	c1->SetLogy();
	graph_evtBkgSB[0]->Draw("AP");
	graph_evtBkgSB[1]->Draw("P");
	if(nState==5)
		graph_evtBkgSB[2]->Draw("P");
	legendRight->Draw();
	if(nState==4)
		latex->DrawLatex(left,top,"J/#psi");
	if(nState==5)
		latex->DrawLatex(left,top,"#psi(2S)");
	c1->SaveAs(Form("%s/evtBkgSB.pdf",savePath.str().c_str()));
	c1->SetLogy(0);

	graph_fracLSB[0]->Draw("AP");
	graph_fracLSB[1]->Draw("P");
	if(nState==5)
		graph_fracLSB[2]->Draw("P");
	legendRight->Draw();
	if(nState==4)
		latex->DrawLatex(left,top,"J/#psi");
	if(nState==5)
		latex->DrawLatex(left,top,"#psi(2S)");
	c1->SaveAs(Form("%s/fracLSB.pdf",savePath.str().c_str()));

	graph_alphaCB[0]->Draw("AP");
	graph_alphaCB[1]->Draw("P");
	if(nState==5)
		graph_alphaCB[2]->Draw("P");
	legendRight->Draw();
	if(nState==4)
		latex->DrawLatex(left,top,"J/#psi");
	if(nState==5)
		latex->DrawLatex(left,top,"#psi(2S)");
	c1->SaveAs(Form("%s/alphaCB.pdf",savePath.str().c_str()));

	graph_lambdaBG[0]->Draw("AP");
	graph_lambdaBG[1]->Draw("P");
	if(nState==5)
		graph_lambdaBG[2]->Draw("P");
	legendRight->Draw();
	if(nState==4)
		latex->DrawLatex(left,top,"J/#psi");
	if(nState==5)
		latex->DrawLatex(left,top,"#psi(2S)");
	c1->SaveAs(Form("%s/lambdaBG.pdf",savePath.str().c_str()));

	graph_fracSigInLSB[0]->Draw("AP");
	graph_fracSigInLSB[1]->Draw("P");
	if(nState==5)
		graph_fracSigInLSB[2]->Draw("P");
	legendRight->Draw();
	if(nState==4)
		latex->DrawLatex(left,top,"J/#psi");
	if(nState==5)
		latex->DrawLatex(left,top,"#psi(2S)");
	c1->SaveAs(Form("%s/fracSigInLSB.pdf",savePath.str().c_str()));

	graph_fracSigInRSB[0]->Draw("AP");
	graph_fracSigInRSB[1]->Draw("P");
	if(nState==5)
		graph_fracSigInRSB[2]->Draw("P");
	legendRight->Draw();
	if(nState==4)
		latex->DrawLatex(left,top,"J/#psi");
	if(nState==5)
		latex->DrawLatex(left,top,"#psi(2S)");
	c1->SaveAs(Form("%s/fracSigInRSB.pdf",savePath.str().c_str()));

	c1->SetLogy();
	graph_evtInLSB[0]->Draw("AP");
	graph_evtInLSB[1]->Draw("P");
	if(nState==5)
		graph_evtInLSB[2]->Draw("P");
	legendRight->Draw();
	if(nState==4)
		latex->DrawLatex(left,top,"J/#psi");
	if(nState==5)
		latex->DrawLatex(left,top,"#psi(2S)");
	c1->SaveAs(Form("%s/evtInLSB.pdf",savePath.str().c_str()));
	c1->SetLogy(0);

	graph_fracBkgInLSB[0]->Draw("AP");
	graph_fracBkgInLSB[1]->Draw("P");
	if(nState==5)
		graph_fracBkgInLSB[2]->Draw("P");
	legendRight->Draw();
	if(nState==4)
		latex->DrawLatex(left,top,"J/#psi");
	if(nState==5)
		latex->DrawLatex(left,top,"#psi(2S)");
	c1->SaveAs(Form("%s/fracBkgInLSB.pdf",savePath.str().c_str()));

	graph_fracBkgInRSB[0]->Draw("AP");
	graph_fracBkgInRSB[1]->Draw("P");
	if(nState==5)
		graph_fracBkgInRSB[2]->Draw("P");
	legendRight->Draw();
	if(nState==4)
		latex->DrawLatex(left,top,"J/#psi");
	if(nState==5)
		latex->DrawLatex(left,top,"#psi(2S)");
	c1->SaveAs(Form("%s/fracBkgInRSB.pdf",savePath.str().c_str()));


	//ploting signal contamination in one single frame
	graph_fracSigInLSB[0]->SetMarkerStyle(20);
	graph_fracSigInLSB[1]->SetMarkerStyle(21);
	graph_fracSigInRSB[0]->SetMarkerStyle(24);
	graph_fracSigInRSB[1]->SetMarkerStyle(25);

	graph_fracSigInRSB[0]->SetMarkerSize(1.2);
	graph_fracSigInRSB[1]->SetMarkerSize(1.2);
	if(nState==5){
		graph_fracSigInLSB[2]->SetMarkerStyle(22);
		graph_fracSigInRSB[2]->SetMarkerStyle(26);
		graph_fracSigInRSB[2]->SetMarkerSize(1.2);
	}
	graph_fracSigInLSB[0]->GetYaxis()->SetTitle("signal contamination in L(R)SB");
	graph_fracSigInLSB[0]->Draw("AP");
	graph_fracSigInLSB[1]->Draw("P");
	graph_fracSigInRSB[0]->Draw("P");
	graph_fracSigInRSB[1]->Draw("P");
	if(nState==5){
		graph_fracSigInLSB[2]->Draw("P");
		graph_fracSigInRSB[2]->Draw("P");
		legendSigConRightFull->Draw();
		latex->DrawLatex(left,top,"#psi(2S)");
	}
	else{
		legendSigConRightFull->Draw();
		latex->DrawLatex(left,top,"J/#psi");
	}
	c1->SaveAs(Form("%s/fracSigInSB.pdf",savePath.str().c_str()));
	graph_fracSigInLSB[0]->GetYaxis()->SetTitle("signal contamination in LSB");


	////
	graph_sigma1[0]->SetMarkerStyle(20);
	graph_sigma1[1]->SetMarkerStyle(21);
	graph_sigma2[0]->SetMarkerStyle(24);
	graph_sigma2[1]->SetMarkerStyle(25);

	graph_sigma2[0]->SetMarkerSize(1.2);
	graph_sigma2[1]->SetMarkerSize(1.2);
	if(nState==5){
		graph_sigma1[2]->SetMarkerStyle(22);
		graph_sigma2[2]->SetMarkerStyle(26);
		graph_sigma2[2]->SetMarkerSize(1.2);
	}

	graph_sigma2[0]->GetYaxis()->SetTitle("#sigma_{1,2} (MeV)");
	graph_sigma2[0]->Draw("AP");
	graph_sigma2[1]->Draw("P");
	graph_sigma1[0]->Draw("P");
	graph_sigma1[1]->Draw("P");
	if(nState==5){
		graph_sigma2[2]->Draw("P");
		graph_sigma1[2]->Draw("P");
		legendRightFull->Draw();
		latex->DrawLatex(left,top,"#psi(2S)");
	}
	else{
		legendFull->Draw();
		latex->DrawLatex(left,top,"J/#psi");
	}
	c1->SaveAs(Form("%s/sigma1_2.pdf",savePath.str().c_str()));
	graph_sigma2[0]->GetYaxis()->SetTitle("sigma_{2} (MeV)");

	///
	TFile *outfile  = new TFile(Form("%s/MassPar_%dS.root",savePath.str().c_str(),nState-3),"RECREATE");
	for(int rapBin = 1; rapBin < RapBins+1; rapBin++){
		graph_mean[rapBin-1]->SetName(Form("graph_mean_%d",rapBin)); graph_mean[rapBin-1]->Write();
		graph_sigmaWei[rapBin-1]->SetName(Form("graph_sigmaWei_%d",rapBin)); graph_sigmaWei[rapBin-1]->Write();
		graph_sigma1[rapBin-1]->SetName(Form("graph_sigma1_%d",rapBin)); graph_sigma1[rapBin-1]->Write();
		graph_sigma2[rapBin-1]->SetName(Form("graph_sigma2_%d",rapBin)); graph_sigma2[rapBin-1]->Write();
		graph_bkgRatio3Sig[rapBin-1]->SetName(Form("graph_bkgRatio3Sig_%d",rapBin)); graph_bkgRatio3Sig[rapBin-1]->Write();
		graph_evtBkgSB[rapBin-1]->SetName(Form("graph_evtBkgSB_%d",rapBin)); graph_evtBkgSB[rapBin-1]->Write();
		graph_fracLSB[rapBin-1]->SetName(Form("graph_fracLSB_%d",rapBin)); graph_fracLSB[rapBin-1]->Write();
	}
	return;
}

//========================================
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

//==============================================
void PlotLifePar(int  nState) {
	int RapBins = onia::kNbRapForPTBins,
			PtBins  = onia::kNbPTMaxBins;
	int LR=2;

	std::stringstream savePath;
	savePath << "Fit/parameter/life";
	gSystem->mkdir(savePath.str().c_str(),kTRUE);

	double pT[RapBins][PtBins];
	double promptMean[LR][RapBins][PtBins];
	double promptCtRe[LR][RapBins][PtBins];
	double promptCtRe2[LR][RapBins][PtBins];
	double promptCtReWei[LR][RapBins][PtBins];
	double fracGauss2[LR][RapBins][PtBins];
	double bkgTauSSDR[LR][RapBins][PtBins];
	double bkgTauSSDL[LR][RapBins][PtBins];
	double bkgTauDSD[LR][RapBins][PtBins];
	double fBkgP[LR][RapBins][PtBins];
	double fBkgSSD[LR][RapBins][PtBins];
	double fBkgSSDL[LR][RapBins][PtBins];
	double fBkgSSDR[LR][RapBins][PtBins];
	double fBkgDSD[LR][RapBins][PtBins];

	double fP[LR][RapBins][PtBins];
	double fNP[LR][RapBins][PtBins];
	double fBkg[LR][RapBins][PtBins];

	double fracBkgSSDL_SBL[LR][RapBins][PtBins];
	double fracBkgSSDR_SBL[LR][RapBins][PtBins];
	double fracBkgDSD_SBL[LR][RapBins][PtBins];
	double fracBkgSSDL_SBR[LR][RapBins][PtBins];
	double fracBkgSSDR_SBR[LR][RapBins][PtBins];
	double fracBkgDSD_SBR[LR][RapBins][PtBins];
	double ratioLD_SBL[LR][RapBins][PtBins];
	double ratioLD_SBR[LR][RapBins][PtBins];


	double pTErr[RapBins][PtBins];
	double promptMeanErr[LR][RapBins][PtBins];
	double promptCtReErr[LR][RapBins][PtBins];
	double promptCtRe2Err[LR][RapBins][PtBins];
	double promptCtReWeiErr[LR][RapBins][PtBins];
	double fracGauss2Err[LR][RapBins][PtBins];
	double bkgTauSSDRErr[LR][RapBins][PtBins];
	double bkgTauSSDLErr[LR][RapBins][PtBins];
	double bkgTauDSDErr[LR][RapBins][PtBins];
	double fBkgPErr[LR][RapBins][PtBins];
	double fBkgSSDErr[LR][RapBins][PtBins];
	double fBkgSSDLErr[LR][RapBins][PtBins];
	double fBkgSSDRErr[LR][RapBins][PtBins];
	double fBkgDSDErr[LR][RapBins][PtBins];

	double fPErr[LR][RapBins][PtBins];
	double fNPErr[LR][RapBins][PtBins];
	double fBkgErr[LR][RapBins][PtBins];

	double fracBkgSSDLErr_SBL[LR][RapBins][PtBins];
	double fracBkgSSDRErr_SBL[LR][RapBins][PtBins];
	double fracBkgDSDErr_SBL[LR][RapBins][PtBins];
	double fracBkgSSDLErr_SBR[LR][RapBins][PtBins];
	double fracBkgSSDRErr_SBR[LR][RapBins][PtBins];
	double fracBkgDSDErr_SBR[LR][RapBins][PtBins];
	double ratioLDErr_SBL[LR][RapBins][PtBins];
	double ratioLDErr_SBR[LR][RapBins][PtBins];

	for(int LRBin = 1; LRBin < 2; LRBin++){
		for(int rapBin = 1; rapBin < RapBins+1; rapBin++){
			for(int ptBin = 1; ptBin < PtBins+1; ptBin++){

				double ptMin=0.,ptMax=0.;
				ptMin = onia::pTRange[rapBin][ptBin-1];
				ptMax = onia::pTRange[rapBin][ptBin];
				pT[rapBin-1][ptBin-1] =  (ptMin + ptMax)/2.;
				pTErr[rapBin-1][ptBin-1] =  (ptMax - ptMin)/2.;

				promptMean[0][rapBin-1][ptBin-1] =  0.;
				promptMean[1][rapBin-1][ptBin-1] =  0.;
				promptMeanErr[0][rapBin-1][ptBin-1] =  0.;
				promptMeanErr[1][rapBin-1][ptBin-1] =  0.;

				promptCtRe[0][rapBin-1][ptBin-1] = 0.;
				promptCtRe[1][rapBin-1][ptBin-1] = 0.;
				promptCtReErr[0][rapBin-1][ptBin-1] = 0.;
				promptCtReErr[1][rapBin-1][ptBin-1] = 0.;

				promptCtRe2[0][rapBin-1][ptBin-1] = 0;
				promptCtRe2[1][rapBin-1][ptBin-1] = 0;
				promptCtRe2Err[0][rapBin-1][ptBin-1] = 0;
				promptCtRe2Err[1][rapBin-1][ptBin-1] = 0;

				promptCtReWei[0][rapBin-1][ptBin-1] = 0;
				promptCtReWei[1][rapBin-1][ptBin-1] = 0;
				promptCtReWeiErr[0][rapBin-1][ptBin-1] = 0;
				promptCtReWeiErr[1][rapBin-1][ptBin-1] = 0;

				fracGauss2[0][rapBin-1][ptBin-1] = 0;
				fracGauss2[1][rapBin-1][ptBin-1] = 0;
				fracGauss2Err[0][rapBin-1][ptBin-1] = 0;
				fracGauss2Err[1][rapBin-1][ptBin-1] = 0;

				bkgTauSSDR[0][rapBin-1][ptBin-1] = 0;
				bkgTauSSDR[1][rapBin-1][ptBin-1] = 0;
				bkgTauSSDRErr[0][rapBin-1][ptBin-1] = 0;
				bkgTauSSDRErr[1][rapBin-1][ptBin-1] = 0;

				bkgTauSSDL[0][rapBin-1][ptBin-1] = 0;
				bkgTauSSDL[1][rapBin-1][ptBin-1] = 0;
				bkgTauSSDLErr[0][rapBin-1][ptBin-1] = 0;
				bkgTauSSDLErr[1][rapBin-1][ptBin-1] = 0;

				bkgTauDSD[0][rapBin-1][ptBin-1] = 0;
				bkgTauDSD[1][rapBin-1][ptBin-1] = 0;
				bkgTauDSDErr[0][rapBin-1][ptBin-1] = 0;
				bkgTauDSDErr[1][rapBin-1][ptBin-1] = 0;

				fBkgSSDL[0][rapBin-1][ptBin-1] = 0;
				fBkgSSDL[1][rapBin-1][ptBin-1] = 0;
				fBkgSSDLErr[0][rapBin-1][ptBin-1] = 0;
				fBkgSSDLErr[1][rapBin-1][ptBin-1] = 0;

				fBkgSSDR[0][rapBin-1][ptBin-1] = 0;
				fBkgSSDR[1][rapBin-1][ptBin-1] = 0;
				fBkgSSDRErr[0][rapBin-1][ptBin-1] = 0;
				fBkgSSDRErr[1][rapBin-1][ptBin-1] = 0;

				fBkgDSD[0][rapBin-1][ptBin-1] = 0;
				fBkgDSD[1][rapBin-1][ptBin-1] = 0;
				fBkgDSDErr[0][rapBin-1][ptBin-1] = 0;
				fBkgDSDErr[1][rapBin-1][ptBin-1] = 0;

				fracBkgSSDL_SBL[0][rapBin-1][ptBin-1] = 0;
				fracBkgSSDL_SBL[1][rapBin-1][ptBin-1] = 0;
				fracBkgSSDLErr_SBL[0][rapBin-1][ptBin-1] = 0;
				fracBkgSSDLErr_SBL[1][rapBin-1][ptBin-1] = 0;

				fracBkgSSDL_SBR[0][rapBin-1][ptBin-1] = 0;
				fracBkgSSDL_SBR[1][rapBin-1][ptBin-1] = 0;
				fracBkgSSDLErr_SBR[0][rapBin-1][ptBin-1] = 0;
				fracBkgSSDLErr_SBR[1][rapBin-1][ptBin-1] = 0;

				fracBkgSSDR_SBL[0][rapBin-1][ptBin-1] = 0;
				fracBkgSSDR_SBL[1][rapBin-1][ptBin-1] = 0;
				fracBkgSSDRErr_SBL[0][rapBin-1][ptBin-1] = 0;
				fracBkgSSDRErr_SBL[1][rapBin-1][ptBin-1] = 0;
				fracBkgSSDR_SBR[0][rapBin-1][ptBin-1] = 0;
				fracBkgSSDR_SBR[1][rapBin-1][ptBin-1] = 0;
				fracBkgSSDRErr_SBR[0][rapBin-1][ptBin-1] = 0;
				fracBkgSSDRErr_SBR[1][rapBin-1][ptBin-1] = 0;

				fracBkgDSD_SBL[0][rapBin-1][ptBin-1] = 0;
				fracBkgDSD_SBL[1][rapBin-1][ptBin-1] = 0;
				fracBkgDSDErr_SBL[0][rapBin-1][ptBin-1] = 0;
				fracBkgDSDErr_SBL[1][rapBin-1][ptBin-1] = 0;

				fracBkgDSD_SBR[0][rapBin-1][ptBin-1] = 0;
				fracBkgDSD_SBR[1][rapBin-1][ptBin-1] = 0;
				fracBkgDSDErr_SBR[0][rapBin-1][ptBin-1] = 0;
				fracBkgDSDErr_SBR[1][rapBin-1][ptBin-1] = 0;

				ratioLD_SBL[0][rapBin-1][ptBin-1] = 0;
				ratioLD_SBL[1][rapBin-1][ptBin-1] = 0;
				ratioLDErr_SBL[0][rapBin-1][ptBin-1] = 0;
				ratioLDErr_SBL[1][rapBin-1][ptBin-1] = 0;

				ratioLD_SBR[0][rapBin-1][ptBin-1] = 0;
				ratioLD_SBR[1][rapBin-1][ptBin-1] = 0;
				ratioLDErr_SBR[0][rapBin-1][ptBin-1] = 0;
				ratioLDErr_SBR[1][rapBin-1][ptBin-1] = 0;

				fP[0][rapBin-1][ptBin-1] = 0;
				fP[1][rapBin-1][ptBin-1] = 0;
				fPErr[0][rapBin-1][ptBin-1] = 0;
				fPErr[1][rapBin-1][ptBin-1] = 0;

				fNP[0][rapBin-1][ptBin-1] = 0;
				fNP[1][rapBin-1][ptBin-1] = 0;
				fNPErr[0][rapBin-1][ptBin-1] = 0;
				fNPErr[1][rapBin-1][ptBin-1] = 0;

				fBkg[0][rapBin-1][ptBin-1] = 0;
				fBkg[1][rapBin-1][ptBin-1] = 0;
				fBkgErr[0][rapBin-1][ptBin-1] = 0;
				fBkgErr[1][rapBin-1][ptBin-1] = 0;
			}
		}
	}

	TFile *inFile;
	RooWorkspace *ws;
	char inName[200];
	for(int LRBin = 1; LRBin < 2; LRBin++){
		for(int rapBin = 1; rapBin < RapBins+1; rapBin++){
			for(int ptBin = 1; ptBin < PtBins+1; ptBin++){

				sprintf(inName, Form("tmpFiles/fit_Psi%dS_rap%d_pt%d.root",nState-3,rapBin,ptBin));

				cout<<"inName: "<<inName<<endl;
				inFile=new TFile(inName,"R");
				if(!inFile) continue;
				ws=(RooWorkspace *)inFile->Get("ws_masslifetime");
				if(!ws){ cout<<">>=======Error: no workspace in root file=========="<<endl; continue; }
				RooRealVar *promptMean_=(RooRealVar*)ws->var("promptMean");
				if(!promptMean_) continue;

				RooRealVar *ctResolution_=(RooRealVar*)ws->var("ctResolution");
				RooRealVar *ctResolution2_=(RooRealVar*)ws->var("ctResolution2");
				RooRealVar *fracGauss2_=(RooRealVar*)ws->var("fracGauss2");

				RooRealVar *bkgTauSSD_SBL_=(RooRealVar*)ws->var("bkgTauSSD_SBL");
				RooRealVar *bkgTauFD_=(RooRealVar*)ws->var("bkgTauFD");
				RooRealVar *bkgTauDSD_=(RooRealVar*)ws->var("bkgTauDSD");

				RooRealVar *bkgTauSSD_SBR_=(RooRealVar*)ws->var("bkgTauSSD_SBR");

				RooRealVar *fBkgSSDR_SBL_ = (RooRealVar*)ws->var("fBkgSSDR_SBL");

				RooRealVar *fBkgSSDR_SBR_ = (RooRealVar*)ws->var("fBkgSSDR_SBR");

				RooRealVar *fBkgDSD_SBL_;
				RooRealVar *fBkgDSD_SBR_;
				fBkgDSD_SBL_ = (RooRealVar*)ws->var("fBkgDSD_SBL");
				fBkgDSD_SBR_ = (RooRealVar*)ws->var("fBkgDSD_SBR");

				double PromptMean = promptMean_->getVal();
				double PromptMeanErr = promptMean_->getError();
				double PromptCtRe = ctResolution_->getVal();
				double PromptCtReErr = ctResolution_->getError();
				double PromptCtRe2 = ctResolution2_->getVal();
				double PromptCtRe2Err = ctResolution2_->getError();
				double FracGauss2 = fracGauss2_->getVal();
				double FracGauss2Err = fracGauss2_->getError();
				cout<<"rap "<<rapBin<<" pt "<<ptBin<<endl;
				cout<<"FracGauss2: "<<FracGauss2<<endl;
				cout<<"PromptCtRe2: "<<PromptCtRe2<<endl<<endl;

				double bkgTauSSD_SBL = bkgTauSSD_SBL_->getVal();
				double bkgTauSSD_SBLErr = bkgTauSSD_SBL_->getError();
				double BkgTauFD = bkgTauFD_->getVal();
				double BkgTauFDErr = bkgTauFD_->getError();
				double BkgTauDSD = bkgTauDSD_->getVal();
				double BkgTauDSDErr = bkgTauDSD_->getError();

				double bkgTauSSD_SBR = bkgTauSSD_SBR_->getVal();
				double bkgTauSSD_SBRErr = bkgTauSSD_SBR_->getError();

				double fBkgSSDR_SBL = fBkgSSDR_SBL_->getVal();
				double fBkgSSDRErr_SBL = fBkgSSDR_SBL_->getError();

				double fBkgSSDR_SBR = fBkgSSDR_SBR_->getVal();
				double fBkgSSDRErr_SBR = fBkgSSDR_SBR_->getError();

				double fBkgDSD_SBL,fBkgDSDErr_SBL,fBkgSSDL_SBL,fBkgSSDLErr_SBL,
							 fBkgDSD_SBR,fBkgDSDErr_SBR,fBkgSSDL_SBR,fBkgSSDLErr_SBR;
				double RatioLD_SBL=0, RatioLDErr_SBL=0, RatioLD_SBR=0, RatioLDErr_SBR=0;

				fBkgDSD_SBL = fBkgDSD_SBL_->getVal();
				fBkgDSDErr_SBL = fBkgDSD_SBL_->getError();
				fBkgSSDL_SBL = 1-fBkgSSDR_SBL-fBkgDSD_SBL;
				fBkgSSDLErr_SBL = sqrt(pow(fBkgSSDRErr_SBL,2)+pow(fBkgDSDErr_SBL,2));

				fBkgDSD_SBR = fBkgDSD_SBR_->getVal();
				fBkgDSDErr_SBR = fBkgDSD_SBR_->getError();
				fBkgSSDL_SBR = 1-fBkgSSDR_SBR-fBkgDSD_SBR;
				fBkgSSDLErr_SBR = sqrt(pow(fBkgSSDRErr_SBR,2)+pow(fBkgDSDErr_SBR,2));

				if(fBkgDSD_SBL!=0.){
					RatioLD_SBL = fBkgSSDL_SBL/fBkgDSD_SBL;
					RatioLDErr_SBL = RatioLD_SBL*sqrt(pow(fBkgSSDLErr_SBL/fBkgSSDL_SBL,2)+pow(fBkgDSDErr_SBL/fBkgDSD_SBL,2));
				}
				if(fBkgDSD_SBR!=0.) {
					RatioLD_SBR = fBkgSSDL_SBR/fBkgDSD_SBR;
					RatioLDErr_SBR = RatioLD_SBR*sqrt(pow(fBkgSSDLErr_SBR/fBkgSSDL_SBR,2)+pow(fBkgDSDErr_SBR/fBkgDSD_SBR,2));
				}

				promptMean[0][rapBin-1][ptBin-1] =  PromptMean;
				promptMean[1][rapBin-1][ptBin-1] =  PromptMean;
				promptMeanErr[0][rapBin-1][ptBin-1] =  0.;
				promptMeanErr[1][rapBin-1][ptBin-1] =  0.;

				promptCtRe[0][rapBin-1][ptBin-1] = PromptCtRe;
				promptCtRe[1][rapBin-1][ptBin-1] = PromptCtRe;
				promptCtReErr[0][rapBin-1][ptBin-1] = 0.;
				promptCtReErr[1][rapBin-1][ptBin-1] = 0.;

				promptCtRe2[0][rapBin-1][ptBin-1] = PromptCtRe2;
				promptCtRe2[1][rapBin-1][ptBin-1] = PromptCtRe2;
				promptCtRe2Err[0][rapBin-1][ptBin-1] = PromptCtRe2Err;
				promptCtRe2Err[1][rapBin-1][ptBin-1] = PromptCtRe2Err;

				double FracGauss1 = 1- FracGauss2;
				double FracGauss1Err = FracGauss2Err;
				double PromptCtReWei = sqrt(pow(PromptCtRe,2)*FracGauss1+pow(PromptCtRe2,2)*(1-FracGauss1));
				double PromptCtReWeiErr = (1./(2*PromptCtReWei))*
					((pow(PromptCtRe,2)+pow(PromptCtRe2,2))*FracGauss1Err +
					 2*FracGauss1*PromptCtRe*PromptCtReErr + 2*(1-FracGauss1)*PromptCtRe2*PromptCtRe2Err);

				promptCtReWei[0][rapBin-1][ptBin-1] = PromptCtReWei;
				promptCtReWei[1][rapBin-1][ptBin-1] = PromptCtReWei;
				promptCtReWeiErr[0][rapBin-1][ptBin-1] = PromptCtReWeiErr;
				promptCtReWeiErr[1][rapBin-1][ptBin-1] = PromptCtReWeiErr;

				fracGauss2[0][rapBin-1][ptBin-1] = FracGauss2;
				fracGauss2[1][rapBin-1][ptBin-1] = FracGauss2;
				fracGauss2Err[0][rapBin-1][ptBin-1] = FracGauss2Err;
				fracGauss2Err[1][rapBin-1][ptBin-1] = FracGauss2Err;

				bkgTauSSDR[0][rapBin-1][ptBin-1] = bkgTauSSD_SBL;
				bkgTauSSDR[1][rapBin-1][ptBin-1] = bkgTauSSD_SBR;
				bkgTauSSDRErr[0][rapBin-1][ptBin-1] = bkgTauSSD_SBLErr;
				bkgTauSSDRErr[1][rapBin-1][ptBin-1] = bkgTauSSD_SBRErr;

				bkgTauSSDL[0][rapBin-1][ptBin-1] = BkgTauFD;
				bkgTauSSDL[1][rapBin-1][ptBin-1] = BkgTauFD;
				bkgTauSSDLErr[0][rapBin-1][ptBin-1] = BkgTauFDErr;
				bkgTauSSDLErr[1][rapBin-1][ptBin-1] = BkgTauFDErr;

				bkgTauDSD[0][rapBin-1][ptBin-1] = BkgTauDSD;
				bkgTauDSD[1][rapBin-1][ptBin-1] = BkgTauDSD;
				bkgTauDSDErr[0][rapBin-1][ptBin-1] = BkgTauDSDErr;
				bkgTauDSDErr[1][rapBin-1][ptBin-1] = BkgTauDSDErr;

				fracBkgSSDL_SBL[0][rapBin-1][ptBin-1] = fBkgSSDL_SBL;
				fracBkgSSDL_SBR[0][rapBin-1][ptBin-1] = fBkgSSDL_SBR;
				fracBkgSSDLErr_SBL[0][rapBin-1][ptBin-1] = fBkgSSDLErr_SBL;
				fracBkgSSDLErr_SBR[0][rapBin-1][ptBin-1] = fBkgSSDLErr_SBR;

				fracBkgSSDR_SBL[0][rapBin-1][ptBin-1] = fBkgSSDR_SBL;
				fracBkgSSDR_SBR[0][rapBin-1][ptBin-1] = fBkgSSDR_SBR;
				fracBkgSSDRErr_SBL[0][rapBin-1][ptBin-1] = fBkgSSDRErr_SBL;
				fracBkgSSDRErr_SBR[0][rapBin-1][ptBin-1] = fBkgSSDRErr_SBR;

				fracBkgDSD_SBL[0][rapBin-1][ptBin-1] = fBkgDSD_SBL;
				fracBkgDSD_SBR[0][rapBin-1][ptBin-1] = fBkgDSD_SBR;
				fracBkgDSDErr_SBL[0][rapBin-1][ptBin-1] = fBkgDSDErr_SBL;
				fracBkgDSDErr_SBR[0][rapBin-1][ptBin-1] = fBkgDSDErr_SBR;

				ratioLD_SBL[0][rapBin-1][ptBin-1] = RatioLD_SBL;
				ratioLD_SBR[0][rapBin-1][ptBin-1] = RatioLD_SBR;
				ratioLDErr_SBL[0][rapBin-1][ptBin-1] = RatioLDErr_SBL;
				ratioLDErr_SBR[0][rapBin-1][ptBin-1] = RatioLDErr_SBR;
			}
		}
	}

	RooRealVar JpsiMass(*ws->var("JpsiMass"));

	TGraphErrors *graph_promptMean[LR][RapBins], *graph_promptCtRe[LR][RapBins], *graph_promptCtRe2[LR][RapBins],
							 *graph_fracGauss2[LR][RapBins], *graph_promptCtReWei[LR][RapBins],
							 *graph_bkgTauSSDR[LR][RapBins], *graph_bkgTauSSDL[LR][RapBins], *graph_bkgTauDSD[LR][RapBins],
							 *graph_fracBkgSSDL_SBL[LR][RapBins], *graph_fracBkgSSDR_SBL[LR][RapBins], *graph_fracBkgDSD_SBL[LR][RapBins],
							 *graph_fracBkgSSDL_SBR[LR][RapBins], *graph_fracBkgSSDR_SBR[LR][RapBins], *graph_fracBkgDSD_SBR[LR][RapBins],
							 *graph_ratioLD_SBL[LR][RapBins], *graph_ratioLD_SBR[LR][RapBins];

	for(int rapBin = 1; rapBin < RapBins+1; rapBin++){
		graph_promptMean[0][rapBin-1] = new TGraphErrors(PtBins, pT[rapBin-1], promptMean[0][rapBin-1], pTErr[rapBin-1], promptMeanErr[0][rapBin-1]);
		graph_promptMean[1][rapBin-1] = new TGraphErrors(PtBins, pT[rapBin-1], promptMean[1][rapBin-1], pTErr[rapBin-1], promptMeanErr[1][rapBin-1]);

		graph_promptCtRe[0][rapBin-1] = new TGraphErrors(PtBins, pT[rapBin-1], promptCtRe[0][rapBin-1], pTErr[rapBin-1], promptCtReErr[0][rapBin-1]);
		graph_promptCtRe[1][rapBin-1] = new TGraphErrors(PtBins, pT[rapBin-1], promptCtRe[1][rapBin-1], pTErr[rapBin-1], promptCtReErr[1][rapBin-1]);

		graph_promptCtRe2[0][rapBin-1] = new TGraphErrors(PtBins, pT[rapBin-1], promptCtRe2[0][rapBin-1], pTErr[rapBin-1], promptCtRe2Err[0][rapBin-1]);
		graph_promptCtRe2[1][rapBin-1] = new TGraphErrors(PtBins, pT[rapBin-1], promptCtRe2[1][rapBin-1], pTErr[rapBin-1], promptCtRe2Err[1][rapBin-1]);

		graph_promptCtReWei[0][rapBin-1] = new TGraphErrors(PtBins, pT[rapBin-1], promptCtReWei[0][rapBin-1], pTErr[rapBin-1], promptCtReWeiErr[0][rapBin-1]);
		graph_promptCtReWei[1][rapBin-1] = new TGraphErrors(PtBins, pT[rapBin-1], promptCtReWei[1][rapBin-1], pTErr[rapBin-1], promptCtReWeiErr[1][rapBin-1]);

		graph_fracGauss2[0][rapBin-1] = new TGraphErrors(PtBins, pT[rapBin-1], fracGauss2[0][rapBin-1], pTErr[rapBin-1], fracGauss2Err[0][rapBin-1]);
		graph_fracGauss2[1][rapBin-1] = new TGraphErrors(PtBins, pT[rapBin-1], fracGauss2[1][rapBin-1], pTErr[rapBin-1], fracGauss2Err[1][rapBin-1]);

		graph_bkgTauSSDR[0][rapBin-1] = new TGraphErrors(PtBins, pT[rapBin-1], bkgTauSSDR[0][rapBin-1], pTErr[rapBin-1], bkgTauSSDRErr[0][rapBin-1]);
		graph_bkgTauSSDR[1][rapBin-1] = new TGraphErrors(PtBins, pT[rapBin-1], bkgTauSSDR[1][rapBin-1], pTErr[rapBin-1], bkgTauSSDRErr[1][rapBin-1]);

		graph_bkgTauSSDL[0][rapBin-1] = new TGraphErrors(PtBins, pT[rapBin-1], bkgTauSSDL[0][rapBin-1], pTErr[rapBin-1], bkgTauSSDLErr[0][rapBin-1]);
		graph_bkgTauSSDL[1][rapBin-1] = new TGraphErrors(PtBins, pT[rapBin-1], bkgTauSSDL[1][rapBin-1], pTErr[rapBin-1], bkgTauSSDLErr[1][rapBin-1]);

		graph_bkgTauDSD[0][rapBin-1] = new TGraphErrors(PtBins, pT[rapBin-1], bkgTauDSD[0][rapBin-1], pTErr[rapBin-1], bkgTauDSDErr[0][rapBin-1]);
		graph_bkgTauDSD[1][rapBin-1] = new TGraphErrors(PtBins, pT[rapBin-1], bkgTauDSD[1][rapBin-1], pTErr[rapBin-1], bkgTauDSDErr[1][rapBin-1]);

		graph_fracBkgSSDL_SBL[0][rapBin-1] = new TGraphErrors(PtBins, pT[rapBin-1], fracBkgSSDL_SBL[0][rapBin-1], pTErr[rapBin-1], fracBkgSSDLErr_SBL[0][rapBin-1]);
		graph_fracBkgSSDR_SBL[0][rapBin-1] = new TGraphErrors(PtBins, pT[rapBin-1], fracBkgSSDR_SBL[0][rapBin-1], pTErr[rapBin-1], fracBkgSSDRErr_SBL[0][rapBin-1]);
		graph_fracBkgDSD_SBL[0][rapBin-1] = new TGraphErrors(PtBins, pT[rapBin-1], fracBkgDSD_SBL[0][rapBin-1], pTErr[rapBin-1], fracBkgDSDErr_SBL[0][rapBin-1]);
		graph_fracBkgSSDL_SBR[0][rapBin-1] = new TGraphErrors(PtBins, pT[rapBin-1], fracBkgSSDL_SBR[0][rapBin-1], pTErr[rapBin-1], fracBkgSSDLErr_SBR[0][rapBin-1]);
		graph_fracBkgSSDR_SBR[0][rapBin-1] = new TGraphErrors(PtBins, pT[rapBin-1], fracBkgSSDR_SBR[0][rapBin-1], pTErr[rapBin-1], fracBkgSSDRErr_SBR[0][rapBin-1]);
		graph_fracBkgDSD_SBR[0][rapBin-1] = new TGraphErrors(PtBins, pT[rapBin-1], fracBkgDSD_SBR[0][rapBin-1], pTErr[rapBin-1], fracBkgDSDErr_SBR[0][rapBin-1]);

		graph_ratioLD_SBL[0][rapBin-1] = new TGraphErrors(PtBins, pT[rapBin-1], ratioLD_SBL[0][rapBin-1], pTErr[rapBin-1], ratioLDErr_SBL[0][rapBin-1]);
		graph_ratioLD_SBR[0][rapBin-1] = new TGraphErrors(PtBins, pT[rapBin-1], ratioLD_SBR[0][rapBin-1], pTErr[rapBin-1], ratioLDErr_SBR[0][rapBin-1]);

		//remove first pT bin for 2S
		if(nState==5){
			//graph_promptMean[0][rapBin-1] -> RemovePoint(0);
			//graph_promptMean[1][rapBin-1] -> RemovePoint(0);
			//graph_promptCtRe[0][rapBin-1] -> RemovePoint(0);
			//graph_promptCtRe[1][rapBin-1] -> RemovePoint(0);
			//graph_promptCtRe2[0][rapBin-1] -> RemovePoint(0);
			//graph_promptCtRe2[1][rapBin-1] -> RemovePoint(0);
			//graph_promptCtReWei[0][rapBin-1] -> RemovePoint(0);
			//graph_promptCtReWei[1][rapBin-1] -> RemovePoint(0);
			//graph_fracGauss2[0][rapBin-1] -> RemovePoint(0);
			//graph_fracGauss2[1][rapBin-1] -> RemovePoint(0);
			//graph_bkgTauSSDR[0][rapBin-1] -> RemovePoint(0);
			//graph_bkgTauSSDR[1][rapBin-1] -> RemovePoint(0);
			//graph_bkgTauSSDL[0][rapBin-1] -> RemovePoint(0);
			//graph_bkgTauSSDL[1][rapBin-1] -> RemovePoint(0);
			//graph_bkgTauDSD[0][rapBin-1] -> RemovePoint(0);
			//graph_bkgTauDSD[1][rapBin-1] -> RemovePoint(0);
			//graph_fracBkgSSDL_SBL[0][rapBin-1] -> RemovePoint(0);
			//graph_fracBkgSSDR_SBL[0][rapBin-1] -> RemovePoint(0);
			//graph_fracBkgDSD_SBL[0][rapBin-1] -> RemovePoint(0);
			//graph_fracBkgSSDL_SBR[0][rapBin-1] -> RemovePoint(0);
			//graph_fracBkgSSDR_SBR[0][rapBin-1] -> RemovePoint(0);
			//graph_fracBkgDSD_SBR[0][rapBin-1] -> RemovePoint(0);
			//graph_ratioLD_SBL[0][rapBin-1] -> RemovePoint(0);
			//graph_ratioLD_SBR[0][rapBin-1] -> RemovePoint(0);
		}
	}

	//TLegend* LifetimeLegend=new TLegend(0.75,0.75,0.88,0.88);
	TLegend* LifetimeLegend=new TLegend(0.75,0.75,0.84,0.84);
	LifetimeLegend->SetFillColor(kWhite);
	LifetimeLegend->SetTextFont(42);
	LifetimeLegend->SetTextSize(legendsize);
	LifetimeLegend->SetBorderSize(0.);
	LifetimeLegend->AddEntry(graph_promptMean[0][0],"SBL","lp");
	LifetimeLegend->AddEntry(graph_promptMean[1][0],"SBR","lp");

	TLegend* LifetimeLegendLR=new TLegend(0.75,0.75,0.84,0.84);
	LifetimeLegendLR->SetFillColor(kWhite);
	LifetimeLegendLR->SetTextFont(42);
	LifetimeLegendLR->SetTextSize(legendsize);
	LifetimeLegendLR->SetBorderSize(0.);
	LifetimeLegendLR->AddEntry(graph_promptMean[0][0],"SBL(R)","lp");

	gStyle->SetPadBottomMargin(0.11); //0.12
	gStyle->SetPadLeftMargin(0.08); //0.12
	gStyle->SetPadRightMargin(0.02); //0.05
	gStyle->SetPadTopMargin(0.02); //0.05

	TCanvas *c1=new TCanvas("c1","");
	c1->SetTickx();
	c1->SetTicky();
	//c1->SetGridx();
	//c1->SetGridy();
	gStyle->SetTitleFillColor(10);
	gStyle->SetTitleBorderSize(1);
	gStyle->SetOptFit(1);          
	gStyle->SetOptStat(1);         
	gStyle->SetTitleFont(22);      
	gStyle->SetStatFont(22);       
	gStyle->SetStatColor(10);      
	gStyle->SetStatBorderSize(1);
	gStyle->SetLabelFont(22,"X");  
	gStyle->SetLabelFont(22,"Y");  
	gStyle->SetTitleXOffset(1.2);  
	gStyle->SetTitleYOffset(1.2);  
	gStyle->SetHistLineWidth(2);   
	gStyle->SetStatX(0.9);         
	gStyle->SetStatY(0.9);         
	gStyle->SetTitleX(0.15);       
	gStyle->SetTitleY(0.96);       

	double Xmin = 0., Xmax = 100.;

	for(int rapBin = 1; rapBin < RapBins+1; rapBin++){

		//shape parameter plot
		c1->SetLogy(0);
		graph_promptMean[0][rapBin-1]->SetTitle("");
		graph_promptMean[0][rapBin-1]->GetXaxis()->SetTitle("p_{T} (GeV)");
		graph_promptMean[0][rapBin-1]->GetYaxis()->SetTitle("prompt Mean (mm)");
		graph_promptMean[0][rapBin-1]->GetYaxis()->SetRangeUser(-0.015, 0.015);
		graph_promptMean[0][rapBin-1]->GetXaxis()->SetLimits(Xmin, Xmax);
		graph_promptMean[0][rapBin-1]->SetMarkerStyle(onia::marker_rapForPTBins[2]);
		graph_promptMean[0][rapBin-1]->SetMarkerColor(onia::colour_rapForPTBins[2]);
		graph_promptMean[0][rapBin-1]->SetLineColor(onia::colour_rapForPTBins[2]);
		graph_promptMean[0][rapBin-1]->Draw("AP");
		graph_promptMean[1][rapBin-1]->SetMarkerStyle(onia::marker_rapForPTBins[3]);
		graph_promptMean[1][rapBin-1]->SetMarkerColor(onia::colour_rapForPTBins[3]);
		graph_promptMean[1][rapBin-1]->SetLineColor(onia::colour_rapForPTBins[3]);
		//graph_promptMean[1][rapBin-1]->Draw("P");
		LifetimeLegendLR->Draw();
		c1->SaveAs(Form("%s/promptMean_rap%d.pdf",savePath.str().c_str(),rapBin));

		graph_promptCtRe[0][rapBin-1]->SetTitle("");
		graph_promptCtRe[0][rapBin-1]->GetXaxis()->SetTitle("p_{T} (GeV)");
		graph_promptCtRe[0][rapBin-1]->GetYaxis()->SetTitle("prompt Resolution (mm)");
		graph_promptCtRe[0][rapBin-1]->GetYaxis()->SetRangeUser(0.0,2.0);
		graph_promptCtRe[0][rapBin-1]->GetXaxis()->SetLimits(Xmin, Xmax);
		graph_promptCtRe[0][rapBin-1]->SetMarkerStyle(onia::marker_rapForPTBins[2]);
		graph_promptCtRe[0][rapBin-1]->SetMarkerColor(onia::colour_rapForPTBins[2]);
		graph_promptCtRe[0][rapBin-1]->SetLineColor(onia::colour_rapForPTBins[2]);
		graph_promptCtRe[0][rapBin-1]->Draw("AP");
		graph_promptCtRe[1][rapBin-1]->SetMarkerStyle(onia::marker_rapForPTBins[3]);
		graph_promptCtRe[1][rapBin-1]->SetMarkerColor(onia::colour_rapForPTBins[3]);
		graph_promptCtRe[1][rapBin-1]->SetLineColor(onia::colour_rapForPTBins[3]);
		//graph_promptCtRe[1][rapBin-1]->Draw("P");
		LifetimeLegendLR->Draw();
		c1->SaveAs(Form("%s/promptCtRe_rap%d.pdf",savePath.str().c_str(),rapBin));

		graph_promptCtRe2[0][rapBin-1]->SetTitle("");
		graph_promptCtRe2[0][rapBin-1]->GetXaxis()->SetTitle("p_{T} (GeV)");
		graph_promptCtRe2[0][rapBin-1]->GetYaxis()->SetTitle("prompt Resolution2 (mm)");
		graph_promptCtRe2[0][rapBin-1]->GetYaxis()->SetRangeUser(-1.0,6.0);
		graph_promptCtRe2[0][rapBin-1]->GetXaxis()->SetLimits(Xmin, Xmax);
		graph_promptCtRe2[0][rapBin-1]->SetMarkerStyle(onia::marker_rapForPTBins[2]);
		graph_promptCtRe2[0][rapBin-1]->SetMarkerColor(onia::colour_rapForPTBins[2]);
		graph_promptCtRe2[0][rapBin-1]->SetLineColor(onia::colour_rapForPTBins[2]);
		graph_promptCtRe2[0][rapBin-1]->Draw("AP");
		graph_promptCtRe2[1][rapBin-1]->SetMarkerStyle(onia::marker_rapForPTBins[3]);
		graph_promptCtRe2[1][rapBin-1]->SetMarkerColor(onia::colour_rapForPTBins[3]);
		graph_promptCtRe2[1][rapBin-1]->SetLineColor(onia::colour_rapForPTBins[3]);
		//graph_promptCtRe2[1][rapBin-1]->Draw("P");
		LifetimeLegendLR->Draw();
		c1->SaveAs(Form("%s/promptCtRe2_rap%d.pdf",savePath.str().c_str(),rapBin));

		graph_promptCtReWei[0][rapBin-1]->SetTitle("");
		graph_promptCtReWei[0][rapBin-1]->GetXaxis()->SetTitle("p_{T} (GeV)");
		graph_promptCtReWei[0][rapBin-1]->GetYaxis()->SetTitle("prompt effective resolution (mm)");
		graph_promptCtReWei[0][rapBin-1]->GetYaxis()->SetRangeUser(-1.0,6.0);
		graph_promptCtReWei[0][rapBin-1]->GetXaxis()->SetLimits(Xmin, Xmax);
		graph_promptCtReWei[0][rapBin-1]->SetMarkerStyle(onia::marker_rapForPTBins[2]);
		graph_promptCtReWei[0][rapBin-1]->SetMarkerColor(onia::colour_rapForPTBins[2]);
		graph_promptCtReWei[0][rapBin-1]->SetLineColor(onia::colour_rapForPTBins[2]);
		graph_promptCtReWei[0][rapBin-1]->Draw("AP");
		graph_promptCtReWei[1][rapBin-1]->SetMarkerStyle(onia::marker_rapForPTBins[3]);
		graph_promptCtReWei[1][rapBin-1]->SetMarkerColor(onia::colour_rapForPTBins[3]);
		graph_promptCtReWei[1][rapBin-1]->SetLineColor(onia::colour_rapForPTBins[3]);
		//graph_promptCtReWei[1][rapBin-1]->Draw("P");
		LifetimeLegendLR->Draw();
		c1->SaveAs(Form("%s/promptCtReWei_rap%d.pdf",savePath.str().c_str(),rapBin));

		graph_fracGauss2[0][rapBin-1]->SetTitle("");
		graph_fracGauss2[0][rapBin-1]->GetXaxis()->SetTitle("p_{T} (GeV)");
		graph_fracGauss2[0][rapBin-1]->GetYaxis()->SetTitle("fraction of 2nd Gaussian");
		graph_fracGauss2[0][rapBin-1]->GetYaxis()->SetRangeUser(-0.1,0.6);
		graph_fracGauss2[0][rapBin-1]->GetXaxis()->SetLimits(Xmin, Xmax);
		graph_fracGauss2[0][rapBin-1]->SetMarkerStyle(onia::marker_rapForPTBins[2]);
		graph_fracGauss2[0][rapBin-1]->SetMarkerColor(onia::colour_rapForPTBins[2]);
		graph_fracGauss2[0][rapBin-1]->SetLineColor(onia::colour_rapForPTBins[2]);
		graph_fracGauss2[0][rapBin-1]->Draw("AP");
		graph_fracGauss2[1][rapBin-1]->SetMarkerStyle(onia::marker_rapForPTBins[3]);
		graph_fracGauss2[1][rapBin-1]->SetMarkerColor(onia::colour_rapForPTBins[3]);
		graph_fracGauss2[1][rapBin-1]->SetLineColor(onia::colour_rapForPTBins[3]);
		//graph_fracGauss2[1][rapBin-1]->Draw("P");
		LifetimeLegendLR->Draw();
		c1->SaveAs(Form("%s/fracGauss2_rap%d.pdf",savePath.str().c_str(),rapBin));

		//
		graph_bkgTauSSDR[0][rapBin-1]->SetTitle("");
		graph_bkgTauSSDR[0][rapBin-1]->GetXaxis()->SetTitle("p_{T} (GeV)");
		graph_bkgTauSSDR[0][rapBin-1]->GetYaxis()->SetTitle("Tau of SSR ");
		graph_bkgTauSSDR[0][rapBin-1]->GetYaxis()->SetRangeUser(0.1,0.8); //0.25,0.55
		graph_bkgTauSSDR[0][rapBin-1]->GetXaxis()->SetLimits(Xmin, Xmax);
		graph_bkgTauSSDR[0][rapBin-1]->SetMarkerStyle(onia::marker_rapForPTBins[2]);
		graph_bkgTauSSDR[0][rapBin-1]->SetMarkerColor(onia::colour_rapForPTBins[2]);
		graph_bkgTauSSDR[0][rapBin-1]->SetLineColor(onia::colour_rapForPTBins[2]);
		graph_bkgTauSSDR[0][rapBin-1]->Draw("AP");
		graph_bkgTauSSDR[1][rapBin-1]->SetMarkerStyle(onia::marker_rapForPTBins[3]);
		graph_bkgTauSSDR[1][rapBin-1]->SetMarkerColor(onia::colour_rapForPTBins[3]);
		graph_bkgTauSSDR[1][rapBin-1]->SetLineColor(onia::colour_rapForPTBins[3]);
		graph_bkgTauSSDR[1][rapBin-1]->Draw("P");
		LifetimeLegend->Draw();
		c1->SaveAs(Form("%s/bkgTauSSR_rap%d.pdf",savePath.str().c_str(),rapBin));

		graph_bkgTauSSDL[0][rapBin-1]->SetTitle("");
		graph_bkgTauSSDL[0][rapBin-1]->GetXaxis()->SetTitle("p_{T} (GeV)");
		graph_bkgTauSSDL[0][rapBin-1]->GetYaxis()->SetTitle("Tau of SSL");
		graph_bkgTauSSDL[0][rapBin-1]->GetYaxis()->SetRangeUser(-.05,.25);
		graph_bkgTauSSDL[0][rapBin-1]->GetXaxis()->SetLimits(Xmin, Xmax);
		graph_bkgTauSSDL[0][rapBin-1]->SetMarkerStyle(onia::marker_rapForPTBins[2]);
		graph_bkgTauSSDL[0][rapBin-1]->SetMarkerColor(onia::colour_rapForPTBins[2]);
		graph_bkgTauSSDL[0][rapBin-1]->SetLineColor(onia::colour_rapForPTBins[2]);
		graph_bkgTauSSDL[0][rapBin-1]->Draw("AP");
		graph_bkgTauSSDL[1][rapBin-1]->SetMarkerStyle(onia::marker_rapForPTBins[3]);
		graph_bkgTauSSDL[1][rapBin-1]->SetMarkerColor(onia::colour_rapForPTBins[3]);
		graph_bkgTauSSDL[1][rapBin-1]->SetLineColor(onia::colour_rapForPTBins[3]);
		//graph_bkgTauSSDL[1][rapBin-1]->Draw("P");
		LifetimeLegendLR->Draw();
		c1->SaveAs(Form("%s/bkgTauSSL_rap%d.pdf",savePath.str().c_str(),rapBin));

		graph_bkgTauDSD[0][rapBin-1]->SetTitle("");
		graph_bkgTauDSD[0][rapBin-1]->GetXaxis()->SetTitle("p_{T} (GeV)");
		graph_bkgTauDSD[0][rapBin-1]->GetYaxis()->SetTitle("Tau of DS");
		graph_bkgTauDSD[0][rapBin-1]->GetYaxis()->SetRangeUser(-0.01,0.04);
		graph_bkgTauDSD[0][rapBin-1]->GetXaxis()->SetLimits(Xmin, Xmax);
		graph_bkgTauDSD[0][rapBin-1]->SetMarkerStyle(onia::marker_rapForPTBins[2]);
		graph_bkgTauDSD[0][rapBin-1]->SetMarkerColor(onia::colour_rapForPTBins[2]);
		graph_bkgTauDSD[0][rapBin-1]->SetLineColor(onia::colour_rapForPTBins[2]);
		graph_bkgTauDSD[0][rapBin-1]->Draw("AP");
		graph_bkgTauDSD[1][rapBin-1]->SetMarkerStyle(onia::marker_rapForPTBins[3]);
		graph_bkgTauDSD[1][rapBin-1]->SetMarkerColor(onia::colour_rapForPTBins[3]);
		graph_bkgTauDSD[1][rapBin-1]->SetLineColor(onia::colour_rapForPTBins[3]);
		//graph_bkgTauDSD[1][rapBin-1]->Draw("P");
		LifetimeLegendLR->Draw();
		c1->SaveAs(Form("%s/bkgTauDS_rap%d.pdf",savePath.str().c_str(),rapBin));

		//
		graph_fracBkgSSDL_SBL[0][rapBin-1]->SetTitle("");
		graph_fracBkgSSDL_SBL[0][rapBin-1]->GetXaxis()->SetTitle("p_{T} (GeV)");
		graph_fracBkgSSDL_SBL[0][rapBin-1]->GetYaxis()->SetTitle("fraction of SSL");
		graph_fracBkgSSDL_SBL[0][rapBin-1]->GetYaxis()->SetRangeUser(-0.05,0.15);
		graph_fracBkgSSDL_SBL[0][rapBin-1]->GetXaxis()->SetLimits(Xmin, Xmax);
		graph_fracBkgSSDL_SBL[0][rapBin-1]->SetMarkerStyle(onia::marker_rapForPTBins[2]);
		graph_fracBkgSSDL_SBL[0][rapBin-1]->SetMarkerColor(onia::colour_rapForPTBins[2]);
		graph_fracBkgSSDL_SBL[0][rapBin-1]->SetLineColor(onia::colour_rapForPTBins[2]);
		graph_fracBkgSSDL_SBL[0][rapBin-1]->Draw("AP");
		graph_fracBkgSSDL_SBR[0][rapBin-1]->SetMarkerStyle(onia::marker_rapForPTBins[3]);
		graph_fracBkgSSDL_SBR[0][rapBin-1]->SetMarkerColor(onia::colour_rapForPTBins[3]);
		graph_fracBkgSSDL_SBR[0][rapBin-1]->SetLineColor(onia::colour_rapForPTBins[3]);
		graph_fracBkgSSDL_SBR[0][rapBin-1]->Draw("P");
		LifetimeLegend->Draw();
		c1->SaveAs(Form("%s/fBkgSSL_rap%d.pdf",savePath.str().c_str(),rapBin));

		graph_fracBkgSSDR_SBL[0][rapBin-1]->SetTitle("");
		graph_fracBkgSSDR_SBL[0][rapBin-1]->GetXaxis()->SetTitle("p_{T} (GeV)");
		graph_fracBkgSSDR_SBL[0][rapBin-1]->GetYaxis()->SetTitle("fraction of SSR");
		graph_fracBkgSSDR_SBL[0][rapBin-1]->GetXaxis()->SetLimits(Xmin, Xmax);
		if(nState==4)graph_fracBkgSSDR_SBL[0][rapBin-1]->GetYaxis()->SetRangeUser(0.4,0.9);
		if(nState==5)graph_fracBkgSSDR_SBL[0][rapBin-1]->GetYaxis()->SetRangeUser(0.1,1.0);
		graph_fracBkgSSDR_SBL[0][rapBin-1]->SetMarkerStyle(onia::marker_rapForPTBins[2]);
		graph_fracBkgSSDR_SBL[0][rapBin-1]->SetMarkerColor(onia::colour_rapForPTBins[2]);
		graph_fracBkgSSDR_SBL[0][rapBin-1]->SetLineColor(onia::colour_rapForPTBins[2]);
		graph_fracBkgSSDR_SBL[0][rapBin-1]->Draw("AP");
		graph_fracBkgSSDR_SBR[0][rapBin-1]->SetMarkerStyle(onia::marker_rapForPTBins[3]);
		graph_fracBkgSSDR_SBR[0][rapBin-1]->SetMarkerColor(onia::colour_rapForPTBins[3]);
		graph_fracBkgSSDR_SBR[0][rapBin-1]->SetLineColor(onia::colour_rapForPTBins[3]);
		graph_fracBkgSSDR_SBR[0][rapBin-1]->Draw("P");
		LifetimeLegend->Draw();
		c1->SaveAs(Form("%s/fBkgSSR_rap%d.pdf",savePath.str().c_str(),rapBin));

		graph_fracBkgDSD_SBL[0][rapBin-1]->SetTitle("");
		graph_fracBkgDSD_SBL[0][rapBin-1]->GetXaxis()->SetTitle("p_{T} (GeV)");
		graph_fracBkgDSD_SBL[0][rapBin-1]->GetYaxis()->SetTitle("fraction of DS");
		graph_fracBkgDSD_SBL[0][rapBin-1]->GetXaxis()->SetLimits(Xmin, Xmax);
		if(nState==4) graph_fracBkgDSD_SBL[0][rapBin-1]->GetYaxis()->SetRangeUser(0.1,0.7);
		if(nState==5) graph_fracBkgDSD_SBL[0][rapBin-1]->GetYaxis()->SetRangeUser(0.2,1.0);
		graph_fracBkgDSD_SBL[0][rapBin-1]->SetMarkerStyle(onia::marker_rapForPTBins[2]);
		graph_fracBkgDSD_SBL[0][rapBin-1]->SetMarkerColor(onia::colour_rapForPTBins[2]);
		graph_fracBkgDSD_SBL[0][rapBin-1]->SetLineColor(onia::colour_rapForPTBins[2]);
		graph_fracBkgDSD_SBL[0][rapBin-1]->Draw("AP");
		graph_fracBkgDSD_SBR[0][rapBin-1]->SetMarkerStyle(onia::marker_rapForPTBins[3]);
		graph_fracBkgDSD_SBR[0][rapBin-1]->SetMarkerColor(onia::colour_rapForPTBins[3]);
		graph_fracBkgDSD_SBR[0][rapBin-1]->SetLineColor(onia::colour_rapForPTBins[3]);
		graph_fracBkgDSD_SBR[0][rapBin-1]->Draw("P");
		LifetimeLegend->Draw();
		c1->SaveAs(Form("%s/fBkgDS_rap%d.pdf",savePath.str().c_str(),rapBin));

		graph_ratioLD_SBL[0][rapBin-1]->SetTitle("");
		graph_ratioLD_SBL[0][rapBin-1]->GetXaxis()->SetTitle("p_{T} (GeV)");
		graph_ratioLD_SBL[0][rapBin-1]->GetYaxis()->SetTitle("ratio fraction of SSL/DS");
		graph_ratioLD_SBL[0][rapBin-1]->GetYaxis()->SetRangeUser(-0.1,0.9);
		graph_ratioLD_SBL[0][rapBin-1]->GetXaxis()->SetLimits(Xmin, Xmax);
		graph_ratioLD_SBL[0][rapBin-1]->SetMarkerStyle(onia::marker_rapForPTBins[2]);
		graph_ratioLD_SBL[0][rapBin-1]->SetMarkerColor(onia::colour_rapForPTBins[2]);
		graph_ratioLD_SBL[0][rapBin-1]->SetLineColor(onia::colour_rapForPTBins[2]);
		graph_ratioLD_SBL[0][rapBin-1]->Draw("AP");
		graph_ratioLD_SBR[0][rapBin-1]->SetMarkerStyle(onia::marker_rapForPTBins[3]);
		graph_ratioLD_SBR[0][rapBin-1]->SetMarkerColor(onia::colour_rapForPTBins[3]);
		graph_ratioLD_SBR[0][rapBin-1]->SetLineColor(onia::colour_rapForPTBins[3]);
		graph_ratioLD_SBR[0][rapBin-1]->Draw("P");
		LifetimeLegend->Draw();
		c1->SaveAs(Form("%s/ratioLD_rap%d.pdf",savePath.str().c_str(),rapBin));
	}

	//double blX = 0.6, blY = 0.75, trX = 0.84, trY = 0.84;
	//double blX = 0.65, blY = 0.75, trX = 0.86, trY = 0.86;
	//if(nState==5) blY = 0.65;
	double blX = 0.7, blY = 0.8, trX = 0.93, trY = 0.93;
	if(nState==5) blY = 0.7;
	TLegend* LifetimeLegendLR_2rap=new TLegend(blX,blY,trX,trY);
	LifetimeLegendLR_2rap->SetFillColor(kWhite);
	LifetimeLegendLR_2rap->SetTextFont(42);
	LifetimeLegendLR_2rap->SetTextSize(legendsize);
	LifetimeLegendLR_2rap->SetBorderSize(0.);
	LifetimeLegendLR_2rap->AddEntry(graph_bkgTauSSDL[0][0],"|y| < 0.6 L(R)SB","lp");
	LifetimeLegendLR_2rap->AddEntry(graph_bkgTauSSDL[0][1],"0.6 < |y| < 1.2 L(R)SB","lp");
	if(nState==5)
		LifetimeLegendLR_2rap->AddEntry(graph_bkgTauSSDL[0][2],"1.2 < |y| < 1.5 L(R)SB","lp");

	TLegend* LifetimeLegend_2rap=new TLegend(blX,blY,trX,trY);
	LifetimeLegend_2rap->SetFillColor(kWhite);
	LifetimeLegend_2rap->SetTextFont(42);
	LifetimeLegend_2rap->SetTextSize(legendsize);
	LifetimeLegend_2rap->SetBorderSize(0.);
	LifetimeLegend_2rap->AddEntry(graph_bkgTauSSDR[0][0],"|y| < 0.6 LSB","lp");
	LifetimeLegend_2rap->AddEntry(graph_bkgTauSSDR[1][0],"|y| < 0.6 RSB","lp");
	LifetimeLegend_2rap->AddEntry(graph_bkgTauSSDR[0][1],"0.6 < |y| < 1.2 LSB","lp");
	LifetimeLegend_2rap->AddEntry(graph_bkgTauSSDR[1][1],"0.6 < |y| < 1.2 RSB","lp");
	if(nState==5){
		LifetimeLegend_2rap->AddEntry(graph_bkgTauSSDR[0][2],"1.2 < |y| < 1.5 LSB","lp");
		LifetimeLegend_2rap->AddEntry(graph_bkgTauSSDR[1][2],"1.2 < |y| < 1.5 RSB","lp");
	}

	TLegend* LifetimeLegend_SR_2rap=new TLegend(blX,blY,trX,trY);
	LifetimeLegend_SR_2rap->SetFillColor(kWhite);
	LifetimeLegend_SR_2rap->SetTextFont(42);
	LifetimeLegend_SR_2rap->SetTextSize(legendsize);
	LifetimeLegend_SR_2rap->SetBorderSize(0.);
	LifetimeLegend_SR_2rap->AddEntry(graph_fracGauss2[0][0],"|y| < 0.6","lp");
	LifetimeLegend_SR_2rap->AddEntry(graph_fracGauss2[0][1],"0.6 < |y| < 1.2","lp");
	if(nState==5)
		LifetimeLegend_SR_2rap->AddEntry(graph_fracGauss2[0][2],"1.2 < |y| < 1.5","lp");


	//
	graph_fracGauss2[0][0]->SetMarkerStyle(21);
	graph_fracGauss2[0][0]->SetMarkerColor(onia::colour_rapForPTBins[2]);
	graph_fracGauss2[0][0]->SetLineColor(onia::colour_rapForPTBins[2]);
	graph_fracGauss2[0][1]->SetMarkerStyle(20);
	graph_fracGauss2[0][1]->SetMarkerColor(onia::colour_rapForPTBins[3]);
	graph_fracGauss2[0][1]->SetLineColor(onia::colour_rapForPTBins[3]);
	if(nState==5){
		graph_fracGauss2[0][2]->SetMarkerStyle(22);
		graph_fracGauss2[0][2]->SetMarkerColor(onia::colour_rapForPTBins[4]);
		graph_fracGauss2[0][2]->SetLineColor(onia::colour_rapForPTBins[4]);
	}
	graph_fracGauss2[0][0]->Draw("AP");
	graph_fracGauss2[0][1]->Draw("P");
	if(nState==5)
		graph_fracGauss2[0][2]->Draw("P");
	LifetimeLegend_SR_2rap->Draw();
	c1->SaveAs(Form("%s/frac2ndGauss.pdf",savePath.str().c_str()));

	graph_bkgTauSSDR[0][0]->SetMarkerStyle(21);
	graph_bkgTauSSDR[0][0]->SetMarkerColor(onia::colour_rapForPTBins[2]);
	graph_bkgTauSSDR[0][0]->SetLineColor(onia::colour_rapForPTBins[2]);
	graph_bkgTauSSDR[0][1]->SetMarkerStyle(20);
	graph_bkgTauSSDR[0][1]->SetMarkerColor(onia::colour_rapForPTBins[3]);
	graph_bkgTauSSDR[0][1]->SetLineColor(onia::colour_rapForPTBins[3]);
	graph_bkgTauSSDR[1][0]->SetMarkerStyle(25);
	graph_bkgTauSSDR[1][0]->SetMarkerColor(onia::colour_rapForPTBins[2]);
	graph_bkgTauSSDR[1][0]->SetLineColor(onia::colour_rapForPTBins[2]);
	graph_bkgTauSSDR[1][1]->SetMarkerStyle(24);
	graph_bkgTauSSDR[1][1]->SetMarkerColor(onia::colour_rapForPTBins[3]);
	graph_bkgTauSSDR[1][1]->SetLineColor(onia::colour_rapForPTBins[3]);
	if(nState==5){
		graph_bkgTauSSDR[0][2]->SetMarkerStyle(22);
		graph_bkgTauSSDR[0][2]->SetMarkerColor(onia::colour_rapForPTBins[4]);
		graph_bkgTauSSDR[0][2]->SetLineColor(onia::colour_rapForPTBins[4]);
		graph_bkgTauSSDR[1][2]->SetMarkerStyle(26);
		graph_bkgTauSSDR[1][2]->SetMarkerColor(onia::colour_rapForPTBins[4]);
		graph_bkgTauSSDR[1][2]->SetLineColor(onia::colour_rapForPTBins[4]);
	}

	graph_bkgTauSSDR[0][0]->Draw("AP");
	graph_bkgTauSSDR[0][1]->Draw("P");
	graph_bkgTauSSDR[1][0]->Draw("P");
	graph_bkgTauSSDR[1][1]->Draw("P");
	if(nState==5){
		graph_bkgTauSSDR[0][2]->Draw("P");
		graph_bkgTauSSDR[1][2]->Draw("P");
	}
	LifetimeLegend_2rap->Draw();

	c1->SaveAs(Form("%s/bkgTauSSR.pdf",savePath.str().c_str()));


	graph_bkgTauSSDL[0][0]->SetMarkerStyle(onia::marker_rapForPTBins[2]);
	graph_bkgTauSSDL[0][0]->SetMarkerColor(onia::colour_rapForPTBins[2]);
	graph_bkgTauSSDL[0][0]->SetLineColor(onia::colour_rapForPTBins[2]);
	graph_bkgTauSSDL[0][1]->SetMarkerStyle(onia::marker_rapForPTBins[3]);
	graph_bkgTauSSDL[0][1]->SetMarkerColor(onia::colour_rapForPTBins[3]);
	graph_bkgTauSSDL[0][1]->SetLineColor(onia::colour_rapForPTBins[3]);
	if(nState==5){
		graph_bkgTauSSDL[0][2]->SetMarkerStyle(onia::marker_rapForPTBins[4]);
		graph_bkgTauSSDL[0][2]->SetMarkerColor(onia::colour_rapForPTBins[4]);
		graph_bkgTauSSDL[0][2]->SetLineColor(onia::colour_rapForPTBins[4]);
	}
	graph_bkgTauSSDL[0][0]->Draw("AP");
	graph_bkgTauSSDL[0][1]->Draw("P");
	if(nState==5)
		graph_bkgTauSSDL[0][2]->Draw("P");
	LifetimeLegendLR_2rap->Draw();
	c1->SaveAs(Form("%s/bkgTauSSL.pdf",savePath.str().c_str()));

	graph_bkgTauDSD[0][0]->SetMarkerStyle(onia::marker_rapForPTBins[2]);
	graph_bkgTauDSD[0][0]->SetMarkerColor(onia::colour_rapForPTBins[2]);
	graph_bkgTauDSD[0][0]->SetLineColor(onia::colour_rapForPTBins[2]);
	graph_bkgTauDSD[0][1]->SetMarkerStyle(onia::marker_rapForPTBins[3]);
	graph_bkgTauDSD[0][1]->SetMarkerColor(onia::colour_rapForPTBins[3]);
	graph_bkgTauDSD[0][1]->SetLineColor(onia::colour_rapForPTBins[3]);
	if(nState==5){
		graph_bkgTauDSD[0][2]->SetMarkerStyle(onia::marker_rapForPTBins[4]);
		graph_bkgTauDSD[0][2]->SetMarkerColor(onia::colour_rapForPTBins[4]);
		graph_bkgTauDSD[0][2]->SetLineColor(onia::colour_rapForPTBins[4]);
	}
	graph_bkgTauDSD[0][0]->Draw("AP");
	graph_bkgTauDSD[0][1]->Draw("P");
	if(nState==5)
		graph_bkgTauDSD[0][2]->Draw("P");
	LifetimeLegendLR_2rap->Draw();
	c1->SaveAs(Form("%s/bkgTauDS.pdf",savePath.str().c_str()));

	graph_fracBkgSSDL_SBL[0][0]->SetMarkerStyle(21);
	graph_fracBkgSSDL_SBL[0][0]->SetMarkerColor(onia::colour_rapForPTBins[2]);
	graph_fracBkgSSDL_SBL[0][0]->SetLineColor(onia::colour_rapForPTBins[2]);
	graph_fracBkgSSDL_SBL[0][1]->SetMarkerStyle(20);
	graph_fracBkgSSDL_SBL[0][1]->SetMarkerColor(onia::colour_rapForPTBins[3]);
	graph_fracBkgSSDL_SBL[0][1]->SetLineColor(onia::colour_rapForPTBins[3]);
	graph_fracBkgSSDL_SBR[0][0]->SetMarkerStyle(25);
	graph_fracBkgSSDL_SBR[0][0]->SetMarkerColor(onia::colour_rapForPTBins[2]);
	graph_fracBkgSSDL_SBR[0][0]->SetLineColor(onia::colour_rapForPTBins[2]);
	graph_fracBkgSSDL_SBR[0][1]->SetMarkerStyle(24);
	graph_fracBkgSSDL_SBR[0][1]->SetMarkerColor(onia::colour_rapForPTBins[3]);
	graph_fracBkgSSDL_SBR[0][1]->SetLineColor(onia::colour_rapForPTBins[3]);
	if(nState==5){
		graph_fracBkgSSDL_SBL[0][2]->SetMarkerStyle(22);
		graph_fracBkgSSDL_SBL[0][2]->SetMarkerColor(onia::colour_rapForPTBins[4]);
		graph_fracBkgSSDL_SBL[0][2]->SetLineColor(onia::colour_rapForPTBins[4]);
		graph_fracBkgSSDL_SBR[0][2]->SetMarkerStyle(26);
		graph_fracBkgSSDL_SBR[0][2]->SetMarkerColor(onia::colour_rapForPTBins[4]);
		graph_fracBkgSSDL_SBR[0][2]->SetLineColor(onia::colour_rapForPTBins[4]);
	}
	graph_fracBkgSSDL_SBL[0][0]->Draw("AP");
	graph_fracBkgSSDL_SBL[0][1]->Draw("P");
	graph_fracBkgSSDL_SBL[0][0]->Draw("AP");
	graph_fracBkgSSDL_SBL[0][1]->Draw("P");
	graph_fracBkgSSDL_SBR[0][0]->Draw("P");
	graph_fracBkgSSDL_SBR[0][1]->Draw("P");
	if(nState==5){
		graph_fracBkgSSDL_SBL[0][2]->Draw("P");
		graph_fracBkgSSDL_SBR[0][2]->Draw("P");
	}
	LifetimeLegend_2rap->Draw();
	c1->SaveAs(Form("%s/fBkgSSL.pdf",savePath.str().c_str()));

	graph_fracBkgSSDR_SBL[0][0]->SetMarkerStyle(21);
	graph_fracBkgSSDR_SBL[0][0]->SetMarkerColor(onia::colour_rapForPTBins[2]);
	graph_fracBkgSSDR_SBL[0][0]->SetLineColor(onia::colour_rapForPTBins[2]);
	graph_fracBkgSSDR_SBL[0][1]->SetMarkerStyle(20);
	graph_fracBkgSSDR_SBL[0][1]->SetMarkerColor(onia::colour_rapForPTBins[3]);
	graph_fracBkgSSDR_SBL[0][1]->SetLineColor(onia::colour_rapForPTBins[3]);
	graph_fracBkgSSDR_SBR[0][0]->SetMarkerStyle(25);
	graph_fracBkgSSDR_SBR[0][0]->SetMarkerColor(onia::colour_rapForPTBins[2]);
	graph_fracBkgSSDR_SBR[0][0]->SetLineColor(onia::colour_rapForPTBins[2]);
	graph_fracBkgSSDR_SBR[0][1]->SetMarkerStyle(24);
	graph_fracBkgSSDR_SBR[0][1]->SetMarkerColor(onia::colour_rapForPTBins[3]);
	graph_fracBkgSSDR_SBR[0][1]->SetLineColor(onia::colour_rapForPTBins[3]);
	if(nState==5){
		graph_fracBkgSSDR_SBL[0][2]->SetMarkerStyle(22);
		graph_fracBkgSSDR_SBL[0][2]->SetMarkerColor(onia::colour_rapForPTBins[4]);
		graph_fracBkgSSDR_SBL[0][2]->SetLineColor(onia::colour_rapForPTBins[4]);
		graph_fracBkgSSDR_SBR[0][2]->SetMarkerStyle(26);
		graph_fracBkgSSDR_SBR[0][2]->SetMarkerColor(onia::colour_rapForPTBins[4]);
		graph_fracBkgSSDR_SBR[0][2]->SetLineColor(onia::colour_rapForPTBins[4]);
	}
	graph_fracBkgSSDR_SBL[0][0]->Draw("AP");
	graph_fracBkgSSDR_SBL[0][1]->Draw("P");
	graph_fracBkgSSDR_SBR[0][0]->Draw("P");
	graph_fracBkgSSDR_SBR[0][1]->Draw("P");
	if(nState==5){
		graph_fracBkgSSDR_SBL[0][2]->Draw("P");
		graph_fracBkgSSDR_SBR[0][2]->Draw("P");
	}
	LifetimeLegend_2rap->Draw();
	c1->SaveAs(Form("%s/fBkgSSR.pdf",savePath.str().c_str()));

	graph_fracBkgDSD_SBL[0][0]->SetMarkerStyle(21);
	graph_fracBkgDSD_SBL[0][0]->SetMarkerColor(onia::colour_rapForPTBins[2]);
	graph_fracBkgDSD_SBL[0][0]->SetLineColor(onia::colour_rapForPTBins[2]);
	graph_fracBkgDSD_SBL[0][1]->SetMarkerStyle(20);
	graph_fracBkgDSD_SBL[0][1]->SetMarkerColor(onia::colour_rapForPTBins[3]);
	graph_fracBkgDSD_SBL[0][1]->SetLineColor(onia::colour_rapForPTBins[3]);
	graph_fracBkgDSD_SBR[0][0]->SetMarkerStyle(25);
	graph_fracBkgDSD_SBR[0][0]->SetMarkerColor(onia::colour_rapForPTBins[2]);
	graph_fracBkgDSD_SBR[0][0]->SetLineColor(onia::colour_rapForPTBins[2]);
	graph_fracBkgDSD_SBR[0][1]->SetMarkerStyle(24);
	graph_fracBkgDSD_SBR[0][1]->SetMarkerColor(onia::colour_rapForPTBins[3]);
	graph_fracBkgDSD_SBR[0][1]->SetLineColor(onia::colour_rapForPTBins[3]);
	if(nState==5){
		graph_fracBkgDSD_SBL[0][2]->SetMarkerStyle(22);
		graph_fracBkgDSD_SBL[0][2]->SetMarkerColor(onia::colour_rapForPTBins[4]);
		graph_fracBkgDSD_SBL[0][2]->SetLineColor(onia::colour_rapForPTBins[4]);
		graph_fracBkgDSD_SBR[0][2]->SetMarkerStyle(26);
		graph_fracBkgDSD_SBR[0][2]->SetMarkerColor(onia::colour_rapForPTBins[4]);
		graph_fracBkgDSD_SBR[0][2]->SetLineColor(onia::colour_rapForPTBins[4]);
	}
	graph_fracBkgDSD_SBL[0][0]->Draw("AP");
	graph_fracBkgDSD_SBL[0][1]->Draw("P");
	graph_fracBkgDSD_SBR[0][0]->Draw("P");
	graph_fracBkgDSD_SBR[0][1]->Draw("P");
	if(nState==5){
		graph_fracBkgDSD_SBL[0][2]->Draw("P");
		graph_fracBkgDSD_SBR[0][2]->Draw("P");
	}
	LifetimeLegend_2rap->Draw();
	c1->SaveAs(Form("%s/fBkgDS.pdf",savePath.str().c_str()));
}


//==================================
void PlotBFrac_1S(int nState){
	int RapBins = onia::kNbRapForPTBins,
			PtBins  = onia::kNbPTMaxBins;

	std::stringstream savePath;
	savePath << "Fit/parameter/bfrac";
	gSystem->mkdir(savePath.str().c_str(),kTRUE);

	int const FitMe= 1;
	double pT[RapBins][PtBins];
	double FracBkg[FitMe][RapBins][PtBins];
	double FracPrompt[FitMe][RapBins][PtBins];
	double FracNonPrompt[FitMe][RapBins][PtBins];
	double BFrac[FitMe][RapBins][PtBins];

	double pTErr[RapBins][PtBins];
	double pTErrLow[RapBins][PtBins];
	double pTErrHigh[RapBins][PtBins];
	double FracBkgErr[FitMe][RapBins][PtBins];
	double FracPromptErr[FitMe][RapBins][PtBins];
	double FracNonPromptErr[FitMe][RapBins][PtBins];
	double BFracErr[FitMe][RapBins][PtBins];

	cout<<"Initializing parameters....."<<endl;
	for(int fit = 1; fit < FitMe+1; fit++){
		for(int rapBin = 1; rapBin < RapBins+1; rapBin++){
			for(int ptBin = 1; ptBin < PtBins+1; ptBin++){

				pT[rapBin-1][ptBin-1] =  0.;
				pTErrLow[rapBin-1][ptBin-1] =  0.;
				pTErrHigh[rapBin-1][ptBin-1] =  0.;

				FracBkg[fit-1][rapBin-1][ptBin-1]       = 0.;
				FracPrompt[fit-1][rapBin-1][ptBin-1]    = 0.;
				FracNonPrompt[fit-1][rapBin-1][ptBin-1] = 0.;
				BFrac[fit-1][rapBin-1][ptBin-1]         = 0.;

				FracBkgErr[fit-1][rapBin-1][ptBin-1]       = 0.;
				FracPromptErr[fit-1][rapBin-1][ptBin-1]    = 0.;
				FracNonPromptErr[fit-1][rapBin-1][ptBin-1] = 0.;
				BFracErr[fit-1][rapBin-1][ptBin-1]         = 0.;
			}
		}
	}

	cout<<"done!"<<endl;

	TFile *inFile;
	RooWorkspace *ws;
	RooDataSet *data;
	char inName[200];

	for(int fit = 1; fit < FitMe+1; fit++){
		for(int rapBin = 1; rapBin < RapBins+1; rapBin++){
			for(int ptBin = 1; ptBin < PtBins+1; ptBin++){

				sprintf(inName, Form("tmpFiles/fit_Psi%dS_rap%d_pt%d.root",nState-3,rapBin,ptBin));
				cout<<"inName: "<<inName<<endl;

				inFile=new TFile(inName,"R");
				if(!inFile) continue;
				ws=(RooWorkspace *)inFile->Get("ws_masslifetime");
				if(!ws){ cout<<">>=======Error: no workspace in root file=========="<<endl; continue; }

				if(fit==1){
					double ptMin = onia::pTRange[rapBin][ptBin-1];
					double ptMax = onia::pTRange[rapBin][ptBin];
					pT[rapBin-1][ptBin-1] =  getMeanPt(rapBin,ptBin,inName);
					pTErrLow[rapBin-1][ptBin-1] =  pT[rapBin-1][ptBin-1]-ptMin;
					pTErrHigh[rapBin-1][ptBin-1] =  ptMax-pT[rapBin-1][ptBin-1];
					cout<<"meanPT: rap"<<rapBin<<" pt"<<ptBin<<": "<<pT[rapBin-1][ptBin-1]<<endl;
				}

				RooRealVar *fBkg_;
				RooRealVar *fPrompt_;
				fBkg_ = (RooRealVar*)ws->var("fBkg");
				if(!fBkg_) continue;
				fPrompt_ = (RooRealVar*)ws->var("fPrompt");

				double fBkg = fBkg_->getVal();
				double fPrompt = fPrompt_->getVal();
				double fNonPrompt =  1.-fBkg-fPrompt;
				double fBkgErr = fBkg_->getError();
				double fPromptErr = fPrompt_->getError();
				double fNonPromptErr = sqrt(pow(fBkgErr,2)+pow(fPromptErr,2));

				cout<<"fBkgErr: "<<fBkgErr<<endl;
				cout<<"fPromptErr: "<<fPromptErr<<endl;

				FracBkg[fit-1][rapBin-1][ptBin-1] = fBkg;
				FracPrompt[fit-1][rapBin-1][ptBin-1] = fPrompt;
				FracNonPrompt[fit-1][rapBin-1][ptBin-1] = fNonPrompt;
				BFrac[fit-1][rapBin-1][ptBin-1] = fNonPrompt/(fNonPrompt+fPrompt);

				FracBkgErr[fit-1][rapBin-1][ptBin-1] = fBkgErr;
				FracPromptErr[fit-1][rapBin-1][ptBin-1] = fPromptErr;
				FracNonPromptErr[fit-1][rapBin-1][ptBin-1] = fNonPromptErr;
				BFracErr[fit-1][rapBin-1][ptBin-1] = BFrac[fit-1][rapBin-1][ptBin-1]*
					sqrt(pow(fNonPromptErr/fNonPrompt,2)+pow(fPromptErr/fPrompt,2));

			}
		}
	}

	for(int fit = 1; fit < FitMe+1; fit++){
		for(int rapBin = 1; rapBin < RapBins+1; rapBin++){
			for(int ptBin = 1; ptBin < PtBins+1; ptBin++){
				cout<<"fit "<<fit<<"  rap "<<rapBin<<"   pt "<<ptBin<<endl;
				cout<<"BFrac: "<<BFrac[fit-1][rapBin-1][ptBin-1]<<endl;
				cout<<"BFracErr: "<<BFracErr[fit-1][rapBin-1][ptBin-1]<<endl;
			}
		}
	}

	RooRealVar JpsiMass(*ws->var("JpsiMass"));
	TGraphAsymmErrors *graph_BFrac[FitMe][RapBins], *graph_FracBkg[FitMe][RapBins], *graph_FracNonPrompt[FitMe][RapBins];

	for(int fit = 1; fit < FitMe+1; fit++){
		for(int rapBin = 1; rapBin < RapBins+1; rapBin++){
			graph_BFrac[fit-1][rapBin-1] = new TGraphAsymmErrors(PtBins, pT[rapBin-1], BFrac[fit-1][rapBin-1],
					pTErrLow[rapBin-1], pTErrHigh[rapBin-1], BFracErr[fit-1][rapBin-1], BFracErr[fit-1][rapBin-1]);

			graph_FracNonPrompt[fit-1][rapBin-1] = new TGraphAsymmErrors(PtBins, pT[rapBin-1], FracNonPrompt[fit-1][rapBin-1],
					pTErrLow[rapBin-1], pTErrHigh[rapBin-1], FracNonPromptErr[fit-1][rapBin-1], FracNonPromptErr[fit-1][rapBin-1]);

			graph_FracBkg[fit-1][rapBin-1] = new TGraphAsymmErrors(PtBins, pT[rapBin-1], FracBkg[fit-1][rapBin-1],
					pTErrLow[rapBin-1], pTErrHigh[rapBin-1], FracBkgErr[fit-1][rapBin-1], FracBkgErr[fit-1][rapBin-1]);
		}
	}


	//TLegend* LifetimeLegend2D=new TLegend(0.75,0.75,0.82,0.84);
	//double blX = 0.65, blY = 0.15, trX = 0.85, trY = 0.35;
	//double blX = 0.16-0.05, blY = 0.72+0.05, trX = 0.4-0.05, trY = 0.84+0.05;
	//double blX = 0.15, blY = 0.82, trX = 0.4, trY = 0.94;
	double blX = 0.12, blY = 0.82, trX = 0.4, trY = 0.96;

	TLegend* LifetimeLegend=new TLegend(blX,blY,trX,trY);
	LifetimeLegend->SetFillColor(kWhite);
	LifetimeLegend->SetTextFont(42);
	LifetimeLegend->SetTextSize(legendsize);
	LifetimeLegend->SetBorderSize(0.);
	LifetimeLegend->AddEntry(graph_BFrac[0][0],"|y| < 0.6","lp");
	LifetimeLegend->AddEntry(graph_BFrac[0][1],"0.6 < |y| < 1.2","lp");

	TLegend* LifetimeLegend2Steps=new TLegend(blX,blY,trX,trY);
	LifetimeLegend2Steps->SetFillColor(kWhite);
	LifetimeLegend2Steps->SetTextFont(42);
	LifetimeLegend2Steps->SetTextSize(legendsize);
	LifetimeLegend2Steps->SetBorderSize(0.);
	LifetimeLegend2Steps->AddEntry(graph_BFrac[0][0],"|y| < 0.6","lp");
	LifetimeLegend2Steps->AddEntry(graph_BFrac[0][1],"0.6 < |y| < 1.2","lp");

	gStyle->SetPadBottomMargin(0.11); //0.12
	gStyle->SetPadLeftMargin(0.08); //0.12
	gStyle->SetPadRightMargin(0.02); //0.05
	gStyle->SetPadTopMargin(0.02); //0.05

	TCanvas *c1=new TCanvas("c1","");
	c1->SetTickx();
	c1->SetTicky();
	//c1->SetGridx();
	//c1->SetGridy();
	gStyle->SetTitleFillColor(10);
	gStyle->SetTitleBorderSize(1);
	gStyle->SetOptFit(1);          
	gStyle->SetOptStat(1);         
	gStyle->SetTitleFont(22);      
	gStyle->SetStatFont(22);       
	gStyle->SetStatColor(10);      
	gStyle->SetStatBorderSize(1);
	gStyle->SetLabelFont(22,"X");  
	gStyle->SetLabelFont(22,"Y");  
	gStyle->SetTitleXOffset(1.2);  
	gStyle->SetTitleYOffset(1.0);  
	gStyle->SetHistLineWidth(2);   
	gStyle->SetStatX(0.9);         
	gStyle->SetStatY(0.9);         
	gStyle->SetTitleX(0.15);       
	gStyle->SetTitleY(0.96);       

	//double Ymin = 0.2, Ymax = .8;
	//double Xmin = 0., Xmax = 75.;
	double Ymin = 0., Ymax = 1.;
	double Xmin = 0., Xmax = 100.;

	//double lvalue = 0.615+0.02, tvalue = 0.4-0.01;
	double lvalue = 0.68, tvalue = 0.41;
	double left=lvalue, top=tvalue, textSize=0.035;
	TLatex *latex=new TLatex();
	latex->SetTextFont(42);
	latex->SetNDC(kTRUE);
	latex->SetTextSize(textSize);
	double step=textSize*1.3;

	for(int fit = 1; fit < FitMe+1; fit++){
		for(int rapBin = 1; rapBin < RapBins+1; rapBin++){

		}
	}
	graph_BFrac[0][0]->SetTitle("");
	graph_BFrac[0][0]->GetXaxis()->SetTitle("p_{T} (GeV)");
	graph_BFrac[0][0]->GetYaxis()->SetTitle("B fraction");
	graph_BFrac[0][0]->GetYaxis()->SetRangeUser(Ymin, Ymax);
	graph_BFrac[0][0]->GetXaxis()->SetLimits(Xmin, Xmax);
	graph_BFrac[0][0]->SetMarkerStyle(onia::marker_rapForPTBins[2]);
	graph_BFrac[0][0]->SetMarkerColor(onia::colour_rapForPTBins[2]);
	graph_BFrac[0][0]->SetLineColor(onia::colour_rapForPTBins[2]);
	graph_BFrac[0][1]->SetTitle("");
	graph_BFrac[0][1]->GetXaxis()->SetTitle("p_{T} (GeV)");
	graph_BFrac[0][1]->GetYaxis()->SetTitle("B fraction");
	graph_BFrac[0][1]->GetYaxis()->SetRangeUser(Ymin, Ymax);
	graph_BFrac[0][1]->GetXaxis()->SetLimits(Xmin, Xmax);
	graph_BFrac[0][1]->SetMarkerStyle(onia::marker_rapForPTBins[3]);
	graph_BFrac[0][1]->SetMarkerColor(onia::colour_rapForPTBins[3]);
	graph_BFrac[0][1]->SetLineColor(onia::colour_rapForPTBins[3]);

	graph_BFrac[0][0]->Draw("AP");
	graph_BFrac[0][1]->Draw("P");
	LifetimeLegend2Steps->Draw();
	left=lvalue; top=tvalue;
	latex->DrawLatex(left,top,"J/#psi");
	top -= step;
	latex->DrawLatex(left,top,"CMS Preliminary");
	top -= step;
	latex->DrawLatex(left,top,"pp  #sqrt{s} = 7 TeV");
	c1->SaveAs(Form("%s/BFrac.pdf",savePath.str().c_str()));

	//
	graph_FracNonPrompt[0][0]->SetTitle("");
	graph_FracNonPrompt[0][0]->GetXaxis()->SetTitle("p_{T} (GeV)");
	graph_FracNonPrompt[0][0]->GetYaxis()->SetTitle("NP fraction");
	graph_FracNonPrompt[0][0]->GetYaxis()->SetRangeUser(.1, .8);
	graph_FracNonPrompt[0][0]->GetXaxis()->SetLimits(Xmin, Xmax);
	graph_FracNonPrompt[0][0]->SetMarkerStyle(onia::marker_rapForPTBins[2]);
	graph_FracNonPrompt[0][0]->SetMarkerColor(onia::colour_rapForPTBins[2]);
	graph_FracNonPrompt[0][0]->SetLineColor(onia::colour_rapForPTBins[2]);
	graph_FracNonPrompt[0][1]->SetTitle("");
	graph_FracNonPrompt[0][1]->GetXaxis()->SetTitle("p_{T} (GeV)");
	graph_FracNonPrompt[0][1]->GetYaxis()->SetTitle("NP fraction");
	graph_FracNonPrompt[0][1]->GetYaxis()->SetRangeUser(.1, .8);
	graph_FracNonPrompt[0][1]->GetXaxis()->SetLimits(Xmin, Xmax);

	graph_FracNonPrompt[0][1]->SetMarkerStyle(onia::marker_rapForPTBins[3]);
	graph_FracNonPrompt[0][1]->SetMarkerColor(onia::colour_rapForPTBins[3]);
	graph_FracNonPrompt[0][1]->SetLineColor(onia::colour_rapForPTBins[3]);

	graph_FracNonPrompt[0][0]->Draw("AP");
	graph_FracNonPrompt[0][1]->Draw("P");
	LifetimeLegend2Steps->Draw();
	left=lvalue; top=tvalue;
	latex->DrawLatex(left,top,"J/#psi");
	top -= step;
	latex->DrawLatex(left,top,"CMS Preliminary");
	top -= step;
	latex->DrawLatex(left,top,"pp  #sqrt{s} = 7 TeV");
	c1->SaveAs(Form("%s/FracNonPrompt.pdf",savePath.str().c_str()));

	graph_FracBkg[0][0]->SetTitle("");
	graph_FracBkg[0][0]->GetXaxis()->SetTitle("p_{T} (GeV)");
	graph_FracBkg[0][0]->GetYaxis()->SetTitle("Bg fraction");
	graph_FracBkg[0][0]->GetYaxis()->SetRangeUser(0., .7);
	graph_FracBkg[0][0]->GetXaxis()->SetLimits(Xmin, Xmax);
	graph_FracBkg[0][0]->SetMarkerStyle(onia::marker_rapForPTBins[2]);
	graph_FracBkg[0][0]->SetMarkerColor(onia::colour_rapForPTBins[2]);
	graph_FracBkg[0][0]->SetLineColor(onia::colour_rapForPTBins[2]);
	graph_FracBkg[0][1]->SetTitle("");
	graph_FracBkg[0][1]->GetXaxis()->SetTitle("p_{T} (GeV)");
	graph_FracBkg[0][1]->GetYaxis()->SetTitle("Bg fraction");
	graph_FracBkg[0][1]->GetYaxis()->SetRangeUser(0., .7);
	graph_FracBkg[0][1]->GetXaxis()->SetLimits(Xmin, Xmax);
	graph_FracBkg[0][1]->SetMarkerStyle(onia::marker_rapForPTBins[3]);
	graph_FracBkg[0][1]->SetMarkerColor(onia::colour_rapForPTBins[3]);

	graph_FracBkg[0][1]->SetLineColor(onia::colour_rapForPTBins[3]);

	graph_FracBkg[0][0]->Draw("AP");
	graph_FracBkg[0][1]->Draw("P");
	LifetimeLegend2Steps->Draw();
	left=lvalue; top=tvalue;
	latex->DrawLatex(left,top,"J/#psi");
	top -= step;
	latex->DrawLatex(left,top,"CMS Preliminary");
	top -= step;
	latex->DrawLatex(left,top,"pp  #sqrt{s} = 7 TeV");
	c1->SaveAs(Form("%s/FracBkg.pdf",savePath.str().c_str()));

	// Jpsi cross-section results //
	const int PTBins1 = 10, PTBins2 = 6;
	double pTBinRap1[PTBins1+1] = {8., 9., 10., 11., 12., 13.5, 15., 18., 30., 45., 70.}; //rap 0.0 - 0.9
	double pTBinRap2[PTBins2+1] = {8., 9., 10., 12., 15., 30., 45.}; //rap 0.9 - 1.2

	double pTBinRap1_weight[PTBins1] = { 8.569, 9.515, 10.494, 11.485, 12.711, 14.206, 16.329, 21.77, 35.03, 52.7}; //rap 0.0 - 0.9 
	double pTBinRap2_weight[PTBins2] = {8.524, 9.491, 10.921, 13.312, 18.95, 34.7}; ////rap 0.9 - 1.2

	double BFracRap1[PTBins1] = { 0.271, 0.285, 0.320, 0.345, 0.373, 0.417, 0.454, 0.535, 0.633, 0.646 };
	double BFracErrRap1[PTBins1] = { 0.006, 0.005, 0.005, 0.005, 0.005, 0.006, 0.006, 0.006, 0.015, 0.038 };
	double BFracRap2[PTBins1] = { 0.265, 0.281, 0.324, 0.380, 0.481, 0.616 };
	double BFracErrRap2[PTBins1] = { 0.007, 0.007, 0.006, 0.006, 0.007, 0.029 };

	double pTRap1[PTBins1];
	double pTRap1ErrLow[PTBins1];
	double pTRap1ErrHigh[PTBins1];
	double pTRap2[PTBins2];

	double pTRap2ErrLow[PTBins2];
	double pTRap2ErrHigh[PTBins2];
	for(int pt=0; pt<PTBins1; pt++){
		pTRap1[pt] = pTBinRap1_weight[pt];
		pTRap1ErrLow[pt] = pTRap1[pt]-pTBinRap1[pt];
		pTRap1ErrHigh[pt] = pTBinRap1[pt+1]-pTRap1[pt];
	}
	for(int pt=0; pt<PTBins2; pt++){
		pTRap2[pt] = pTBinRap2_weight[pt];
		pTRap2ErrLow[pt] = pTRap2[pt]-pTBinRap2[pt];
		pTRap2ErrHigh[pt] = pTBinRap2[pt+1]-pTRap2[pt];
	}

	TGraphAsymmErrors *graph_BFracRap1 = new TGraphAsymmErrors(PTBins1,
			pTRap1, BFracRap1, pTRap1ErrLow, pTRap1ErrHigh, BFracErrRap1, BFracErrRap1);
	TGraphAsymmErrors *graph_BFracRap2 = new TGraphAsymmErrors(PTBins2,
			pTRap2, BFracRap2, pTRap2ErrLow, pTRap2ErrHigh, BFracErrRap2, BFracErrRap2);

	graph_BFracRap1->SetTitle("");
	graph_BFracRap1->GetXaxis()->SetTitle("p_{T} (GeV)");
	graph_BFracRap1->GetYaxis()->SetTitle("B fraction");
	graph_BFracRap1->GetYaxis()->SetRangeUser(Ymin, Ymax);
	graph_BFracRap1->GetXaxis()->SetLimits(Xmin, Xmax);
	graph_BFracRap1->SetMarkerStyle(onia::marker_rapForPTBins[2]);
	graph_BFracRap1->SetMarkerColor(onia::colour_rapForPTBins[2]);
	graph_BFracRap1->SetLineColor(onia::colour_rapForPTBins[2]);

	graph_BFracRap2->SetTitle("");
	graph_BFracRap2->GetXaxis()->SetTitle("p_{T} (GeV)");
	graph_BFracRap2->GetYaxis()->SetTitle("B fraction");

	graph_BFracRap2->GetYaxis()->SetRangeUser(Ymin, Ymax);
	graph_BFracRap2->GetXaxis()->SetLimits(Xmin, Xmax);
	graph_BFracRap2->SetMarkerStyle(onia::marker_rapForPTBins[3]);
	graph_BFracRap2->SetMarkerColor(onia::colour_rapForPTBins[3]);
	graph_BFracRap2->SetLineColor(onia::colour_rapForPTBins[3]);


	TLegend* LifetimeLegendCrossSection=new TLegend(blX,blY,trX,trY);
	LifetimeLegendCrossSection->SetFillColor(kWhite);
	LifetimeLegendCrossSection->SetTextFont(42);
	LifetimeLegendCrossSection->SetTextSize(legendsize);
	LifetimeLegendCrossSection->SetBorderSize(0.);
	LifetimeLegendCrossSection->AddEntry(graph_BFracRap1,"|y| < 0.9 (BPH-10-014)","lp");
	LifetimeLegendCrossSection->AddEntry(graph_BFracRap2,"0.9 < |y| < 1.2 (BPH-10-014)","lp");

	TLegend* LifetimeLegendCrossSection2Steps=new TLegend(blX,blY,trX,trY);
	LifetimeLegendCrossSection2Steps->SetFillColor(kWhite);
	LifetimeLegendCrossSection2Steps->SetTextFont(42);
	LifetimeLegendCrossSection2Steps->SetTextSize(legendsize);
	LifetimeLegendCrossSection2Steps->SetBorderSize(0.);
	LifetimeLegendCrossSection2Steps->AddEntry(graph_BFrac[0][0],"|y| < 0.6","lp");
	LifetimeLegendCrossSection2Steps->AddEntry(graph_BFrac[0][1],"0.6 < |y| < 1.2","lp");
	LifetimeLegendCrossSection2Steps->AddEntry(graph_BFracRap1,"|y| < 0.9 (BPH-10-014)","lp");
	LifetimeLegendCrossSection2Steps->AddEntry(graph_BFracRap2,"0.9 < |y| < 1.2 (BPH-10-014)","lp");

	TLegend* LifetimeLegendCrossSection2Steps_rap1=new TLegend(blX,blY,trX,trY);
	LifetimeLegendCrossSection2Steps_rap1->SetFillColor(kWhite);
	LifetimeLegendCrossSection2Steps_rap1->SetTextFont(42);
	LifetimeLegendCrossSection2Steps_rap1->SetTextSize(legendsize);
	LifetimeLegendCrossSection2Steps_rap1->SetBorderSize(0.);
	LifetimeLegendCrossSection2Steps_rap1->AddEntry(graph_BFrac[0][0],"|y| < 0.6","lp");

	LifetimeLegendCrossSection2Steps_rap1->AddEntry(graph_BFracRap1,"|y| < 0.9 (BPH-10-014)","lp");

	TLegend* LifetimeLegendCrossSection2Steps_rap2=new TLegend(blX,blY,trX,trY);
	LifetimeLegendCrossSection2Steps_rap2->SetFillColor(kWhite);
	LifetimeLegendCrossSection2Steps_rap2->SetTextFont(42);
	LifetimeLegendCrossSection2Steps_rap2->SetTextSize(legendsize);
	LifetimeLegendCrossSection2Steps_rap2->SetBorderSize(0.);
	LifetimeLegendCrossSection2Steps_rap2->AddEntry(graph_BFrac[0][1],"0.6 < |y| < 1.2","lp");
	LifetimeLegendCrossSection2Steps_rap2->AddEntry(graph_BFracRap2,"0.9 < |y| < 1.2 (BPH-10-014)","lp");

	graph_BFracRap1->Draw("AP");
	graph_BFracRap2->Draw("P");
	LifetimeLegendCrossSection->Draw();
	left=lvalue; top=tvalue;
	latex->DrawLatex(left,top,"J/#psi");
	top -= step;
	latex->DrawLatex(left,top,"CMS Preliminary");
	top -= step;
	latex->DrawLatex(left,top,"pp  #sqrt{s} = 7 TeV");
	c1->SaveAs(Form("%s/BFrac_BPH-10-014.pdf",savePath.str().c_str()));

	graph_BFrac[0][0]->SetMarkerStyle(onia::marker_rapForPTBins[2]);
	graph_BFrac[0][0]->SetMarkerColor(onia::colour_rapForPTBins[2]);
	graph_BFrac[0][0]->SetLineColor(onia::colour_rapForPTBins[2]);
	graph_BFrac[0][1]->SetMarkerStyle(onia::marker_rapForPTBins[3]);
	graph_BFrac[0][1]->SetMarkerColor(onia::colour_rapForPTBins[3]);
	graph_BFrac[0][1]->SetLineColor(onia::colour_rapForPTBins[3]);
	graph_BFracRap1->SetMarkerStyle(onia::marker_rapForPTBins[4]);
	graph_BFracRap1->SetMarkerColor(onia::colour_rapForPTBins[4]);
	graph_BFracRap1->SetLineColor(onia::colour_rapForPTBins[4]);
	graph_BFracRap2->SetMarkerStyle(onia::marker_rapForPTBins[5]);
	graph_BFracRap2->SetMarkerColor(onia::colour_rapForPTBins[5]);
	graph_BFracRap2->SetLineColor(onia::colour_rapForPTBins[5]);

	graph_BFrac[0][0]->Draw("AP");
	graph_BFrac[0][1]->Draw("P");
	graph_BFracRap1->Draw("P");
	graph_BFracRap2->Draw("P");
	LifetimeLegendCrossSection2Steps->Draw();
	left=lvalue; top=tvalue;
	latex->DrawLatex(left,top,"J/#psi");
	top -= step;
	latex->DrawLatex(left,top,"CMS Preliminary");
	top -= step;
	latex->DrawLatex(left,top,"pp  #sqrt{s} = 7 TeV");
	c1->SaveAs(Form("%s/BFrac_BPH-10-014_compare.pdf",savePath.str().c_str()));

	graph_BFrac[0][0]->SetMarkerStyle(onia::marker_rapForPTBins[2]);
	graph_BFrac[0][0]->SetMarkerColor(onia::colour_rapForPTBins[2]);
	graph_BFrac[0][0]->SetLineColor(onia::colour_rapForPTBins[2]);
	graph_BFrac[0][1]->SetMarkerStyle(onia::marker_rapForPTBins[2]);
	graph_BFrac[0][1]->SetMarkerColor(onia::colour_rapForPTBins[2]);
	graph_BFrac[0][1]->SetLineColor(onia::colour_rapForPTBins[2]);
	graph_BFracRap1->SetMarkerStyle(onia::marker_rapForPTBins[3]);
	graph_BFracRap1->SetMarkerColor(onia::colour_rapForPTBins[3]);
	graph_BFracRap1->SetLineColor(onia::colour_rapForPTBins[3]);

	graph_BFracRap2->SetMarkerStyle(onia::marker_rapForPTBins[3]);
	graph_BFracRap2->SetMarkerColor(onia::colour_rapForPTBins[3]);
	graph_BFracRap2->SetLineColor(onia::colour_rapForPTBins[3]);

	graph_BFrac[0][0]->Draw("AP");
	graph_BFracRap1->Draw("P");
	LifetimeLegendCrossSection2Steps_rap1->Draw();
	left=lvalue; top=tvalue;
	latex->DrawLatex(left,top,"J/#psi");
	top -= step;
	latex->DrawLatex(left,top,"CMS Preliminary");
	top -= step;
	latex->DrawLatex(left,top,"pp  #sqrt{s} = 7 TeV");
	c1->SaveAs(Form("%s/BFrac_BPH-10-014_compare_rap1.pdf",savePath.str().c_str()));

	graph_BFrac[0][1]->Draw("AP");
	graph_BFracRap2->Draw("P");
	LifetimeLegendCrossSection2Steps_rap2->Draw();
	left=lvalue; top=tvalue;
	latex->DrawLatex(left,top,"J/#psi");
	top -= step;
	latex->DrawLatex(left,top,"CMS Preliminary");
	top -= step;
	latex->DrawLatex(left,top,"pp  #sqrt{s} = 7 TeV");
	c1->SaveAs(Form("%s/BFrac_BPH-10-014_compare_rap2.pdf",savePath.str().c_str()));

	///
	TFile *outfile  = new TFile(Form("%s/BFrac_%dS.root",savePath.str().c_str(),nState-3),"RECREATE");

	graph_BFrac[0][0]->SetName("graph_BFrac_rap1");
	graph_BFrac[0][1]->SetName("graph_BFrac_rap2");
	graph_FracBkg[0][0]->SetName("graph_FracBkg_rap1");
	graph_FracBkg[0][1]->SetName("graph_FracBkg_rap2");
	graph_FracNonPrompt[0][0]->SetName("graph_FracNonPrompt_rap1");
	graph_FracNonPrompt[0][1]->SetName("graph_FracNonPrompt_rap2");
	graph_BFracRap1->SetName("graph_BPH_10_014_BFrac_rap1");
	graph_BFracRap2->SetName("graph_BPH_10_014_BFrac_rap2");

	graph_BFrac[0][0]->Write();
	graph_BFrac[0][1]->Write();
	graph_FracBkg[0][0]->Write();
	graph_FracBkg[0][1]->Write();
	graph_FracNonPrompt[0][0]->Write();
	graph_FracNonPrompt[0][1]->Write();
	graph_BFracRap1->Write();
	graph_BFracRap2->Write();
	outfile->Close();
}

//===========================
void PlotBFrac_2S(int nState){
	int RapBins = onia::kNbRapForPTBins,
			PtBins  = onia::kNbPTMaxBins;

	std::stringstream savePath;
	savePath << "Fit/parameter/bfrac";
	gSystem->mkdir(savePath.str().c_str(),kTRUE);

	int const FitMe= 1;

	double pT[RapBins][PtBins];
	double FracBkg[FitMe][RapBins][PtBins];
	double FracPrompt[FitMe][RapBins][PtBins];
	double FracNonPrompt[FitMe][RapBins][PtBins];
	double BFrac[FitMe][RapBins][PtBins];

	double pTErr[RapBins][PtBins];
	double pTErrLow[RapBins][PtBins];
	double pTErrHigh[RapBins][PtBins];
	double FracBkgErr[FitMe][RapBins][PtBins];
	double FracPromptErr[FitMe][RapBins][PtBins];
	double FracNonPromptErr[FitMe][RapBins][PtBins];
	double BFracErr[FitMe][RapBins][PtBins];

	cout<<"Initializing parameters....."<<endl;
	for(int fit = 1; fit < FitMe+1; fit++){
		for(int rapBin = 1; rapBin < RapBins+1; rapBin++){
			for(int ptBin = 1; ptBin < PtBins+1; ptBin++){

				pT[rapBin-1][ptBin-1] =  0.;
				pTErrLow[rapBin-1][ptBin-1] =  0.;
				pTErrHigh[rapBin-1][ptBin-1] =  0.;

				FracBkg[fit-1][rapBin-1][ptBin-1]       = 0.;
				FracPrompt[fit-1][rapBin-1][ptBin-1]    = 0.;
				FracNonPrompt[fit-1][rapBin-1][ptBin-1] = 0.;
				BFrac[fit-1][rapBin-1][ptBin-1]         = 0.;

				FracBkgErr[fit-1][rapBin-1][ptBin-1]       = 0.;
				FracPromptErr[fit-1][rapBin-1][ptBin-1]    = 0.;
				FracNonPromptErr[fit-1][rapBin-1][ptBin-1] = 0.;
				BFracErr[fit-1][rapBin-1][ptBin-1]         = 0.;
			}
		}
	}
	cout<<"done!"<<endl;
	TFile *inFile;
	RooWorkspace *ws;
	RooDataSet *data;
	char inName[200];

	for(int fit = 1; fit < FitMe+1; fit++){
		for(int rapBin = 1; rapBin < RapBins+1; rapBin++){
			for(int ptBin = 1; ptBin < PtBins+1; ptBin++){

				sprintf(inName, Form("tmpFiles/fit_Psi%dS_rap%d_pt%d.root",nState-3,rapBin,ptBin));

				cout<<"inName: "<<inName<<endl;
				inFile=new TFile(inName,"R");
				if(!inFile) continue;
				ws=(RooWorkspace *)inFile->Get("ws_masslifetime");
				if(!ws){ cout<<">>=======Error: no workspace in root file=========="<<endl; continue; }

				if(fit==1){
					double ptMin = onia::pTRange[rapBin][ptBin-1];
					double ptMax = onia::pTRange[rapBin][ptBin];
					pT[rapBin-1][ptBin-1] =  getMeanPt(rapBin,ptBin,inName);
					pTErrLow[rapBin-1][ptBin-1] =  pT[rapBin-1][ptBin-1]-ptMin;
					pTErrHigh[rapBin-1][ptBin-1] =  ptMax-pT[rapBin-1][ptBin-1];
				}

				RooRealVar *fBkg_;
				RooRealVar *fPrompt_;
				fBkg_ = (RooRealVar*)ws->var("fBkg");
				if(!fBkg_) continue;
				fPrompt_ = (RooRealVar*)ws->var("fPrompt");

				double fBkg = fBkg_->getVal();
				double fPrompt = fPrompt_->getVal();
				double fNonPrompt =  1.-fBkg-fPrompt;
				double fBkgErr = fBkg_->getError();
				double fPromptErr = fPrompt_->getError();
				double fNonPromptErr = sqrt(pow(fBkgErr,2)+pow(fPromptErr,2));
				cout<<"fBkgErr: "<<fBkgErr<<endl;
				cout<<"fPromptErr: "<<fPromptErr<<endl;

				FracBkg[fit-1][rapBin-1][ptBin-1] = fBkg;
				FracPrompt[fit-1][rapBin-1][ptBin-1] = fPrompt;
				FracNonPrompt[fit-1][rapBin-1][ptBin-1] = fNonPrompt;
				BFrac[fit-1][rapBin-1][ptBin-1] = fNonPrompt/(fNonPrompt+fPrompt);

				FracBkgErr[fit-1][rapBin-1][ptBin-1] = fBkgErr;
				FracPromptErr[fit-1][rapBin-1][ptBin-1] = fPromptErr;
				FracNonPromptErr[fit-1][rapBin-1][ptBin-1] = fNonPromptErr;
				BFracErr[fit-1][rapBin-1][ptBin-1] = BFrac[fit-1][rapBin-1][ptBin-1]*
					sqrt(pow(fNonPromptErr/fNonPrompt,2)+pow(fPromptErr/fPrompt,2));

			}
		}
	}

	for(int fit = 1; fit < FitMe+1; fit++){
		for(int rapBin = 1; rapBin < RapBins+1; rapBin++){
			for(int ptBin = 1; ptBin < PtBins+1; ptBin++){
				cout<<"fit "<<fit<<"  rap "<<rapBin<<"   pt "<<ptBin<<endl;
				cout<<"BFrac: "<<BFrac[fit-1][rapBin-1][ptBin-1]<<endl;
				cout<<"BFracErr: "<<BFracErr[fit-1][rapBin-1][ptBin-1]<<endl;
			}
		}
	}

	RooRealVar JpsiMass(*ws->var("JpsiMass"));
	TGraphAsymmErrors *graph_BFrac[FitMe][RapBins], *graph_FracBkg[FitMe][RapBins], *graph_FracNonPrompt[FitMe][RapBins];

	for(int fit = 1; fit < FitMe+1; fit++){
		for(int rapBin = 1; rapBin < RapBins+1; rapBin++){
			graph_BFrac[fit-1][rapBin-1] = new TGraphAsymmErrors(PtBins, pT[rapBin-1], BFrac[fit-1][rapBin-1],
					pTErrLow[rapBin-1], pTErrHigh[rapBin-1], BFracErr[fit-1][rapBin-1], BFracErr[fit-1][rapBin-1]);

			graph_FracNonPrompt[fit-1][rapBin-1] = new TGraphAsymmErrors(PtBins, pT[rapBin-1], FracNonPrompt[fit-1][rapBin-1],
					pTErrLow[rapBin-1], pTErrHigh[rapBin-1], FracNonPromptErr[fit-1][rapBin-1], FracNonPromptErr[fit-1][rapBin-1]);

			graph_FracBkg[fit-1][rapBin-1] = new TGraphAsymmErrors(PtBins, pT[rapBin-1], FracBkg[fit-1][rapBin-1],
					pTErrLow[rapBin-1], pTErrHigh[rapBin-1], FracBkgErr[fit-1][rapBin-1], FracBkgErr[fit-1][rapBin-1]);
		}
	}



	//double blX = 0.65, blY = 0.15, trX = 0.85, trY = 0.35;
	//double blX = 0.16-0.05, blY = 0.72+0.05, trX = 0.4-0.05, trY = 0.84+0.05;
	//double blX = 0.15, blY = 0.82, trX = 0.4, trY = 0.93;
	double blX = 0.12, blY = 0.82, trX = 0.4, trY = 0.96;

	TLegend* LifetimeLegend=new TLegend(blX,blY,trX,trY);
	LifetimeLegend->SetFillColor(kWhite);
	LifetimeLegend->SetTextFont(42);
	LifetimeLegend->SetTextSize(legendsize);
	LifetimeLegend->SetBorderSize(0.);
	LifetimeLegend->AddEntry(graph_BFrac[0][0],"|y| < 0.6","lp");
	LifetimeLegend->AddEntry(graph_BFrac[0][1],"0.6 < |y| < 1.2","lp");

	TLegend* LifetimeLegend2Steps=new TLegend(blX,blY,trX,trY);
	LifetimeLegend2Steps->SetFillColor(kWhite);
	LifetimeLegend2Steps->SetTextFont(42);
	LifetimeLegend2Steps->SetTextSize(legendsize);
	LifetimeLegend2Steps->SetBorderSize(0.);
	LifetimeLegend2Steps->AddEntry(graph_BFrac[0][0],"|y| < 0.6","lp");
	LifetimeLegend2Steps->AddEntry(graph_BFrac[0][1],"0.6 < |y| < 1.2","lp");
	LifetimeLegend2Steps->AddEntry(graph_BFrac[0][2],"1.2 < |y| < 1.5","lp");

	gStyle->SetPadBottomMargin(0.11); //0.12
	gStyle->SetPadLeftMargin(0.08); //0.12
	gStyle->SetPadRightMargin(0.02); //0.05
	gStyle->SetPadTopMargin(0.02); //0.05

	TCanvas *c1=new TCanvas("c1","");
	c1->SetTickx();
	c1->SetTicky();
	//c1->SetGridx();
	//c1->SetGridy();
	gStyle->SetTitleFillColor(10);
	gStyle->SetTitleBorderSize(1);
	gStyle->SetOptFit(1);          
	gStyle->SetOptStat(1);         
	gStyle->SetTitleFont(22);      
	gStyle->SetStatFont(22);       
	gStyle->SetStatColor(10);      
	gStyle->SetStatBorderSize(1);
	gStyle->SetLabelFont(22,"X");  
	gStyle->SetLabelFont(22,"Y");  
	gStyle->SetTitleXOffset(1.2);  
	gStyle->SetTitleYOffset(1.0);  
	gStyle->SetHistLineWidth(2);   
	gStyle->SetStatX(0.9);         
	gStyle->SetStatY(0.9);         
	gStyle->SetTitleX(0.15);       
	gStyle->SetTitleY(0.96);       

	//double Ymin = 0.2, Ymax = .8;
	//double Xmin = 0., Xmax = 75.;
	double Ymin = 0., Ymax = 1.;
	double Xmin = 0., Xmax = 100.;

	//double lvalue = 0.615+0.02, tvalue = 0.4-0.01;
	double lvalue = 0.68, tvalue = 0.41;
	double left=lvalue, top=tvalue, textSize=0.035;
	TLatex *latex=new TLatex();
	latex->SetTextFont(42);
	latex->SetNDC(kTRUE);
	latex->SetTextSize(textSize);
	double step=textSize*1.3;

	for(int fit = 1; fit < FitMe+1; fit++){
		for(int rapBin = 1; rapBin < RapBins+1; rapBin++){

		}
	}
	graph_BFrac[0][0]->SetTitle("");
	graph_BFrac[0][0]->GetXaxis()->SetTitle("p_{T} (GeV)");
	graph_BFrac[0][0]->GetYaxis()->SetTitle("B fraction");
	graph_BFrac[0][0]->GetYaxis()->SetRangeUser(Ymin, Ymax);
	graph_BFrac[0][0]->GetXaxis()->SetLimits(Xmin, Xmax);
	graph_BFrac[0][0]->SetMarkerStyle(onia::marker_rapForPTBins[2]);
	graph_BFrac[0][0]->SetMarkerColor(onia::colour_rapForPTBins[2]);
	graph_BFrac[0][0]->SetLineColor(onia::colour_rapForPTBins[2]);
	graph_BFrac[0][1]->SetTitle("");
	graph_BFrac[0][1]->GetXaxis()->SetTitle("p_{T} (GeV)");
	graph_BFrac[0][1]->GetYaxis()->SetTitle("B fraction");
	graph_BFrac[0][1]->GetYaxis()->SetRangeUser(Ymin, Ymax);
	graph_BFrac[0][1]->GetXaxis()->SetLimits(Xmin, Xmax);
	graph_BFrac[0][1]->SetMarkerStyle(onia::marker_rapForPTBins[3]);
	graph_BFrac[0][1]->SetMarkerColor(onia::colour_rapForPTBins[3]);
	graph_BFrac[0][1]->SetLineColor(onia::colour_rapForPTBins[3]);
	graph_BFrac[0][2]->SetTitle("");
	graph_BFrac[0][2]->GetXaxis()->SetTitle("p_{T} (GeV)");
	graph_BFrac[0][2]->GetYaxis()->SetTitle("B fraction");
	graph_BFrac[0][2]->GetYaxis()->SetRangeUser(Ymin, Ymax);
	graph_BFrac[0][2]->GetXaxis()->SetLimits(Xmin, Xmax);
	graph_BFrac[0][2]->SetMarkerStyle(onia::marker_rapForPTBins[4]);
	graph_BFrac[0][2]->SetMarkerColor(onia::colour_rapForPTBins[4]);
	graph_BFrac[0][2]->SetLineColor(onia::colour_rapForPTBins[4]);

	graph_BFrac[0][0]->Draw("AP");
	graph_BFrac[0][1]->Draw("P");
	graph_BFrac[0][2]->Draw("P");
	LifetimeLegend2Steps->Draw();
	left=lvalue; top=tvalue;
	latex->DrawLatex(left,top,"#psi(2S)");
	top -= step;
	latex->DrawLatex(left,top,"CMS Preliminary");
	top -= step;
	latex->DrawLatex(left,top,"pp  #sqrt{s} = 7 TeV");
	c1->SaveAs(Form("%s/BFrac.pdf",savePath.str().c_str()));

	graph_FracNonPrompt[0][0]->SetTitle("");
	graph_FracNonPrompt[0][0]->GetXaxis()->SetTitle("p_{T} (GeV)");
	graph_FracNonPrompt[0][0]->GetYaxis()->SetTitle("NP fraction");
	graph_FracNonPrompt[0][0]->GetYaxis()->SetRangeUser(.1, .8);
	graph_FracNonPrompt[0][0]->GetXaxis()->SetLimits(Xmin, Xmax);

	graph_FracNonPrompt[0][0]->SetMarkerStyle(onia::marker_rapForPTBins[2]);
	graph_FracNonPrompt[0][0]->SetMarkerColor(onia::colour_rapForPTBins[2]);
	graph_FracNonPrompt[0][0]->SetLineColor(onia::colour_rapForPTBins[2]);
	graph_FracNonPrompt[0][1]->SetTitle("");
	graph_FracNonPrompt[0][1]->GetXaxis()->SetTitle("p_{T} (GeV)");
	graph_FracNonPrompt[0][1]->GetYaxis()->SetTitle("NP fraction");
	graph_FracNonPrompt[0][1]->GetYaxis()->SetRangeUser(.1,.8 );
	graph_FracNonPrompt[0][1]->GetXaxis()->SetLimits(Xmin, Xmax);
	graph_FracNonPrompt[0][1]->SetMarkerStyle(onia::marker_rapForPTBins[3]);
	graph_FracNonPrompt[0][1]->SetMarkerColor(onia::colour_rapForPTBins[3]);
	graph_FracNonPrompt[0][1]->SetLineColor(onia::colour_rapForPTBins[3]);
	graph_FracNonPrompt[0][2]->SetTitle("");
	graph_FracNonPrompt[0][2]->GetXaxis()->SetTitle("p_{T} (GeV)");
	graph_FracNonPrompt[0][2]->GetYaxis()->SetTitle("NP fraction");
	graph_FracNonPrompt[0][2]->GetYaxis()->SetRangeUser(.1, .8);
	graph_FracNonPrompt[0][2]->GetXaxis()->SetLimits(Xmin, Xmax);
	graph_FracNonPrompt[0][2]->SetMarkerStyle(onia::marker_rapForPTBins[4]);
	graph_FracNonPrompt[0][2]->SetMarkerColor(onia::colour_rapForPTBins[4]);
	graph_FracNonPrompt[0][2]->SetLineColor(onia::colour_rapForPTBins[4]);

	graph_FracNonPrompt[0][0]->Draw("AP");
	graph_FracNonPrompt[0][1]->Draw("P");
	graph_FracNonPrompt[0][2]->Draw("P");
	LifetimeLegend2Steps->Draw();
	left=lvalue; top=tvalue;
	latex->DrawLatex(left,top,"#psi(2S)");
	top -= step;
	latex->DrawLatex(left,top,"CMS Preliminary");
	top -= step;
	latex->DrawLatex(left,top,"pp  #sqrt{s} = 7 TeV");

	c1->SaveAs(Form("%s/FracNonPrompt.pdf",savePath.str().c_str()));

	graph_FracBkg[0][0]->SetTitle("");
	graph_FracBkg[0][0]->GetXaxis()->SetTitle("p_{T} (GeV)");
	graph_FracBkg[0][0]->GetYaxis()->SetTitle("Bg fraction");
	graph_FracBkg[0][0]->GetYaxis()->SetRangeUser(0., .7);
	graph_FracBkg[0][0]->GetXaxis()->SetLimits(Xmin, Xmax);
	graph_FracBkg[0][0]->SetMarkerStyle(onia::marker_rapForPTBins[2]);
	graph_FracBkg[0][0]->SetMarkerColor(onia::colour_rapForPTBins[2]);
	graph_FracBkg[0][0]->SetLineColor(onia::colour_rapForPTBins[2]);
	graph_FracBkg[0][1]->SetTitle("");
	graph_FracBkg[0][1]->GetXaxis()->SetTitle("p_{T} (GeV)");
	graph_FracBkg[0][1]->GetYaxis()->SetTitle("Bg fraction");
	graph_FracBkg[0][1]->GetYaxis()->SetRangeUser(0., .7);
	graph_FracBkg[0][1]->GetXaxis()->SetLimits(Xmin, Xmax);
	graph_FracBkg[0][1]->SetMarkerStyle(onia::marker_rapForPTBins[3]);
	graph_FracBkg[0][1]->SetMarkerColor(onia::colour_rapForPTBins[3]);
	graph_FracBkg[0][1]->SetLineColor(onia::colour_rapForPTBins[3]);
	graph_FracBkg[0][2]->SetTitle("");
	graph_FracBkg[0][2]->GetXaxis()->SetTitle("p_{T} (GeV)");
	graph_FracBkg[0][2]->GetYaxis()->SetTitle("B fraction");
	graph_FracBkg[0][2]->GetYaxis()->SetRangeUser(0, .7);
	graph_FracBkg[0][2]->GetXaxis()->SetLimits(Xmin, Xmax);
	graph_FracBkg[0][2]->SetMarkerStyle(onia::marker_rapForPTBins[4]);
	graph_FracBkg[0][2]->SetMarkerColor(onia::colour_rapForPTBins[4]);
	graph_FracBkg[0][2]->SetLineColor(onia::colour_rapForPTBins[4]);

	graph_FracBkg[0][0]->Draw("AP");
	graph_FracBkg[0][1]->Draw("P");

	graph_FracBkg[0][2]->Draw("P");
	LifetimeLegend2Steps->Draw();
	left=lvalue; top=tvalue;
	latex->DrawLatex(left,top,"#psi(2S)");
	top -= step;
	latex->DrawLatex(left,top,"CMS Preliminary");
	top -= step;
	latex->DrawLatex(left,top,"pp  #sqrt{s} = 7 TeV");
	c1->SaveAs(Form("%s/FracBkg.pdf",savePath.str().c_str()));
	//// Jpsi cross-section results //
	//const int PTBins1 = 10, PTBins2 = 6;
	//double pTBinRap1[PTBins1+1] = {8., 9., 10., 11., 12., 13.5, 15., 18., 30., 45., 70.}; //rap 0.0 - 0.9
	//double pTBinRap2[PTBins2+1] = {8., 9., 10., 12., 15., 30., 45.}; //rap 0.9 - 1.2

	//double pTBinRap1_weight[PTBins1] = { 8.569, 9.515, 10.494, 11.485, 12.711, 14.206, 16.329, 21.77, 35.03, 52.7}; //rap 0.0 - 0.9 
	//double pTBinRap2_weight[PTBins2] = {8.524, 9.491, 10.921, 13.312, 18.95, 34.7}; ////rap 0.9 - 1.2

	//double BFracRap1[PTBins1] = { 0.271, 0.285, 0.320, 0.345, 0.373, 0.417, 0.454, 0.535, 0.633, 0.646 };
	//double BFracErrRap1[PTBins1] = { 0.006, 0.005, 0.005, 0.005, 0.005, 0.006, 0.006, 0.006, 0.015, 0.038 };
	//double BFracRap2[PTBins1] = { 0.265, 0.281, 0.324, 0.380, 0.481, 0.616 };
	//double BFracErrRap2[PTBins1] = { 0.007, 0.007, 0.006, 0.006, 0.007, 0.029 };

	//Psi' cross-section results //
	//const int PTBins1 = 9, PTBins2 = 7, PTBins3 = 7;
	const int PTBins1 = 9, PTBins2 = 6, PTBins3 = 7;
	double pTBinRap1[PTBins1+1] = {6.5, 8.0, 9.0, 10.0, 11.0, 12.0, 13.5, 15.0, 18.0, 30.0}; //rap 0.0 - 1.2
	//double pTBinRap2[PTBins2+1] = {5.5, 6.5, 8.0, 9.0, 10.0, 12.0, 15.0, 30.0}; //rap 1.2 - 1.6
	double pTBinRap2[PTBins2+1] = {6.5, 8.0, 9.0, 10.0, 12.0, 15.0, 30.0}; //rap 1.2 - 1.6
	double pTBinRap3[PTBins3+1] = {5.5, 6.5, 8.0, 9.0, 10.0, 12.0, 15.0, 30.0}; //rap 1.6 - 2.4

	double pTBinRap1_weight[PTBins1] = {7.60, 8.54, 9.50, 10.50, 11.48, 12.72, 14.23, 16.30, 21.92}; //rap 0.0 - 1.2
	//double pTBinRap2_weight[PTBins2] = {6.15, 7.36, 8.49, 9.50, 10.97, 13.32, 19.07}; //rap 1.2 - 1.6
	double pTBinRap2_weight[PTBins2] = {7.36, 8.49, 9.50, 10.97, 13.32, 19.07}; //rap 1.2 - 1.6
	double pTBinRap3_weight[PTBins3] = {6.20, 7.25, 8.48, 9.52, 10.88, 13.38, 19.01}; //rap 1.6 - 2.4

	double BFracRap1[PTBins1] = {0.326, 0.350, 0.335, 0.378, 0.394, 0.431, 0.456, 0.487, 0.572};
	double BFracErrRap1[PTBins1] = {0.041, 0.025, 0.021, 0.023, 0.024, 0.022, 0.025, 0.024, 0.025};
	//double BFracRap2[PTBins2] = {0.246, 0.308, 0.309, 0.32, 0.366, 0.375, 0.51};
	//double BFracErrRap2[PTBins2] = {0.093, 0.030, 0.035, 0.03, 0.029, 0.032, 0.03};
	double BFracRap2[PTBins2] = {0.308, 0.309, 0.32, 0.366, 0.375, 0.51};
	double BFracErrRap2[PTBins2] = {0.030, 0.035, 0.03, 0.029, 0.032, 0.03};
	double BFracRap3[PTBins3] = {0.21, 0.19, 0.30, 0.30, 0.37, 0.372, 0.45};
	double BFracErrRap3[PTBins3] = {0.06, 0.03, 0.04, 0.04, 0.03, 0.035, 0.04};

	double pTRap1[PTBins1];
	double pTRap1ErrLow[PTBins1];
	double pTRap1ErrHigh[PTBins1];
	double pTRap2[PTBins2];
	double pTRap2ErrLow[PTBins2];
	double pTRap2ErrHigh[PTBins2];
	double pTRap3[PTBins3];
	double pTRap3ErrLow[PTBins3];
	double pTRap3ErrHigh[PTBins3];

	for(int pt=0; pt<PTBins1; pt++){
		pTRap1[pt] = pTBinRap1_weight[pt];
		pTRap1ErrLow[pt] = pTRap1[pt]-pTBinRap1[pt];
		pTRap1ErrHigh[pt] = pTBinRap1[pt+1]-pTRap1[pt];
	}
	for(int pt=0; pt<PTBins2; pt++){
		pTRap2[pt] = pTBinRap2_weight[pt];
		pTRap2ErrLow[pt] = pTRap2[pt]-pTBinRap2[pt];
		pTRap2ErrHigh[pt] = pTBinRap2[pt+1]-pTRap2[pt];
	}
	for(int pt=0; pt<PTBins3; pt++){
		pTRap3[pt] = pTBinRap3_weight[pt];
		pTRap3ErrLow[pt] = pTRap3[pt]-pTBinRap3[pt];
		pTRap3ErrHigh[pt] = pTBinRap3[pt+1]-pTRap3[pt];
	}

	TGraphAsymmErrors *graph_BFracRap1 = new TGraphAsymmErrors(PTBins1,
			pTRap1, BFracRap1, pTRap1ErrLow, pTRap1ErrHigh, BFracErrRap1, BFracErrRap1);
	TGraphAsymmErrors *graph_BFracRap2 = new TGraphAsymmErrors(PTBins2,
			pTRap2, BFracRap2, pTRap2ErrLow, pTRap2ErrHigh, BFracErrRap2, BFracErrRap2);
	TGraphAsymmErrors *graph_BFracRap3 = new TGraphAsymmErrors(PTBins3,
			pTRap3, BFracRap3, pTRap3ErrLow, pTRap3ErrHigh, BFracErrRap3, BFracErrRap3);

	graph_BFracRap1->SetTitle("");
	graph_BFracRap1->GetXaxis()->SetTitle("p_{T} (GeV)");
	graph_BFracRap1->GetYaxis()->SetTitle("B fraction");

	graph_BFracRap1->GetYaxis()->SetRangeUser(Ymin, Ymax);
	graph_BFracRap1->GetXaxis()->SetLimits(Xmin, Xmax);
	graph_BFracRap1->SetMarkerStyle(onia::marker_rapForPTBins[2]);
	graph_BFracRap1->SetMarkerColor(onia::colour_rapForPTBins[2]);
	graph_BFracRap1->SetLineColor(onia::colour_rapForPTBins[2]);

	graph_BFracRap2->SetTitle("");
	graph_BFracRap2->GetXaxis()->SetTitle("p_{T} (GeV)");
	graph_BFracRap2->GetYaxis()->SetTitle("B fraction");
	graph_BFracRap2->GetYaxis()->SetRangeUser(Ymin, Ymax);
	graph_BFracRap2->GetXaxis()->SetLimits(Xmin, Xmax);
	graph_BFracRap2->SetMarkerStyle(onia::marker_rapForPTBins[3]);
	graph_BFracRap2->SetMarkerColor(onia::colour_rapForPTBins[3]);
	graph_BFracRap2->SetLineColor(onia::colour_rapForPTBins[3]);

	graph_BFracRap3->SetTitle("");
	graph_BFracRap3->GetXaxis()->SetTitle("p_{T} (GeV)");
	graph_BFracRap3->GetYaxis()->SetTitle("B fraction");
	graph_BFracRap3->GetYaxis()->SetRangeUser(Ymin, Ymax);
	graph_BFracRap3->GetXaxis()->SetLimits(Xmin, Xmax);
	graph_BFracRap3->SetMarkerStyle(onia::marker_rapForPTBins[4]);
	graph_BFracRap3->SetMarkerColor(onia::colour_rapForPTBins[4]);
	graph_BFracRap3->SetLineColor(onia::colour_rapForPTBins[4]);


	TLegend* LifetimeLegendCrossSection=new TLegend(blX,blY,trX,trY);
	LifetimeLegendCrossSection->SetFillColor(kWhite);
	LifetimeLegendCrossSection->SetTextFont(42);
	LifetimeLegendCrossSection->SetTextSize(legendsize);
	LifetimeLegendCrossSection->SetBorderSize(0.);
	LifetimeLegendCrossSection->AddEntry(graph_BFracRap1,"|y| < 1.2 (BPH-10-014)","lp");
	LifetimeLegendCrossSection->AddEntry(graph_BFracRap2,"1.2 < |y| < 1.6 (BPH-10-014)","lp");
	LifetimeLegendCrossSection->AddEntry(graph_BFracRap3,"1.6 < |y| < 2.4 (BPH-10-014)","lp");

	TLegend* LifetimeLegendCrossSection2Steps=new TLegend(blX,blY,trX,trY);
	LifetimeLegendCrossSection2Steps->SetFillColor(kWhite);
	LifetimeLegendCrossSection2Steps->SetTextFont(42);
	LifetimeLegendCrossSection2Steps->SetTextSize(legendsize);
	LifetimeLegendCrossSection2Steps->SetBorderSize(0.);
	LifetimeLegendCrossSection2Steps->AddEntry(graph_BFrac[0][0],"|y| < 0.6","lp");
	LifetimeLegendCrossSection2Steps->AddEntry(graph_BFrac[0][1],"0.6 < |y| < 1.2","lp");
	LifetimeLegendCrossSection2Steps->AddEntry(graph_BFrac[0][2],"1.2 < |y| < 1.5","lp");
	LifetimeLegendCrossSection2Steps->AddEntry(graph_BFracRap1,"|y| < 1.2 (BPH-10-014)","lp");
	LifetimeLegendCrossSection2Steps->AddEntry(graph_BFracRap2,"1.2 < |y| < 1.6 (BPH-10-014)","lp");
	//LifetimeLegendCrossSection2Steps->AddEntry(graph_BFracRap3,"1.6 < |y| < 2.4 (BPH-10-014)","lp");

	TLegend* LifetimeLegendCrossSection2Steps_rap1=new TLegend(blX,blY,trX,trY);
	LifetimeLegendCrossSection2Steps_rap1->SetFillColor(kWhite);
	LifetimeLegendCrossSection2Steps_rap1->SetTextFont(42);
	LifetimeLegendCrossSection2Steps_rap1->SetTextSize(legendsize);
	LifetimeLegendCrossSection2Steps_rap1->SetBorderSize(0.);
	LifetimeLegendCrossSection2Steps_rap1->AddEntry(graph_BFrac[0][0],"|y| < 0.6","lp");
	LifetimeLegendCrossSection2Steps_rap1->AddEntry(graph_BFrac[0][1],"0.6 < |y| < 1.2","lp");
	LifetimeLegendCrossSection2Steps_rap1->AddEntry(graph_BFracRap1,"|y| < 1.2 (BPH-10-014)","lp");


	TLegend* LifetimeLegendCrossSection2Steps_rap2=new TLegend(blX,blY,trX,trY);
	LifetimeLegendCrossSection2Steps_rap2->SetFillColor(kWhite);
	LifetimeLegendCrossSection2Steps_rap2->SetTextFont(42);
	LifetimeLegendCrossSection2Steps_rap2->SetTextSize(legendsize);
	LifetimeLegendCrossSection2Steps_rap2->SetBorderSize(0.);
	LifetimeLegendCrossSection2Steps_rap2->AddEntry(graph_BFrac[0][2],"1.2 < |y| < 1.5","lp");
	LifetimeLegendCrossSection2Steps_rap2->AddEntry(graph_BFracRap2,"1.2 < |y| < 1.8 (BPH-10-014)","lp");

	graph_BFracRap1->Draw("AP");
	graph_BFracRap2->Draw("P");
	graph_BFracRap3->Draw("P");
	LifetimeLegendCrossSection->Draw();
	left=lvalue; top=tvalue;
	latex->DrawLatex(left,top,"#psi(2S)");
	top -= step;
	latex->DrawLatex(left,top,"CMS Preliminary");
	top -= step;
	latex->DrawLatex(left,top,"pp  #sqrt{s} = 7 TeV");
	c1->SaveAs(Form("%s/BFrac_BPH-10-014.pdf",savePath.str().c_str()));

	graph_BFrac[0][0]->SetMarkerStyle(onia::marker_rapForPTBins[2]);
	graph_BFrac[0][0]->SetMarkerColor(onia::colour_rapForPTBins[2]);
	graph_BFrac[0][0]->SetLineColor(onia::colour_rapForPTBins[2]);
	graph_BFrac[0][1]->SetMarkerStyle(onia::marker_rapForPTBins[3]);
	graph_BFrac[0][1]->SetMarkerColor(onia::colour_rapForPTBins[3]);
	graph_BFrac[0][1]->SetLineColor(onia::colour_rapForPTBins[3]);
	graph_BFrac[0][2]->SetMarkerStyle(onia::marker_rapForPTBins[1]);
	graph_BFrac[0][2]->SetMarkerColor(onia::colour_rapForPTBins[1]);
	graph_BFrac[0][2]->SetLineColor(onia::colour_rapForPTBins[1]);

	graph_BFracRap1->SetMarkerStyle(onia::marker_rapForPTBins[4]);
	graph_BFracRap1->SetMarkerColor(onia::colour_rapForPTBins[4]);
	graph_BFracRap1->SetLineColor(onia::colour_rapForPTBins[4]);
	graph_BFracRap2->SetMarkerStyle(onia::marker_rapForPTBins[5]);
	graph_BFracRap2->SetMarkerColor(onia::colour_rapForPTBins[5]);
	graph_BFracRap2->SetLineColor(onia::colour_rapForPTBins[5]);

	graph_BFrac[0][0]->Draw("AP");
	graph_BFrac[0][1]->Draw("P");
	graph_BFrac[0][2]->Draw("P");
	graph_BFracRap1->Draw("P");
	graph_BFracRap2->Draw("P");
	LifetimeLegendCrossSection2Steps->Draw();
	left=lvalue; top=tvalue;
	latex->DrawLatex(left,top,"#psi(2S)");
	top -= step;
	latex->DrawLatex(left,top,"CMS Preliminary");
	top -= step;
	latex->DrawLatex(left,top,"pp  #sqrt{s} = 7 TeV");
	c1->SaveAs(Form("%s/BFrac_BPH-10-014_compare.pdf",savePath.str().c_str()));

	graph_BFrac[0][0]->SetMarkerStyle(onia::marker_rapForPTBins[2]);
	graph_BFrac[0][0]->SetMarkerColor(onia::colour_rapForPTBins[2]);
	graph_BFrac[0][0]->SetLineColor(onia::colour_rapForPTBins[2]);
	graph_BFrac[0][1]->SetMarkerStyle(onia::marker_rapForPTBins[3]);
	graph_BFrac[0][1]->SetMarkerColor(onia::colour_rapForPTBins[3]);
	graph_BFrac[0][1]->SetLineColor(onia::colour_rapForPTBins[3]);

	graph_BFrac[0][2]->SetMarkerStyle(onia::marker_rapForPTBins[3]);
	graph_BFrac[0][2]->SetMarkerColor(onia::colour_rapForPTBins[3]);

	graph_BFrac[0][2]->SetLineColor(onia::colour_rapForPTBins[3]);

	graph_BFracRap1->SetMarkerStyle(onia::marker_rapForPTBins[4]);
	graph_BFracRap1->SetMarkerColor(onia::colour_rapForPTBins[4]);
	graph_BFracRap1->SetLineColor(onia::colour_rapForPTBins[4]);

	graph_BFracRap2->SetMarkerStyle(onia::marker_rapForPTBins[4]);
	graph_BFracRap2->SetMarkerColor(onia::colour_rapForPTBins[4]);
	graph_BFracRap2->SetLineColor(onia::colour_rapForPTBins[4]);

	graph_BFrac[0][0]->Draw("AP");
	graph_BFrac[0][1]->Draw("P");
	graph_BFracRap1->Draw("P");
	LifetimeLegendCrossSection2Steps_rap1->Draw();
	left=lvalue; top=tvalue;
	latex->DrawLatex(left,top,"#psi(2S)");
	top -= step;
	latex->DrawLatex(left,top,"CMS Preliminary");
	top -= step;
	latex->DrawLatex(left,top,"pp  #sqrt{s} = 7 TeV");
	c1->SaveAs(Form("%s/BFrac_BPH-10-014_compare_rap1.pdf",savePath.str().c_str()));

	graph_BFrac[0][2]->Draw("AP");
	graph_BFracRap2->Draw("P");
	LifetimeLegendCrossSection2Steps_rap2->Draw();
	left=lvalue; top=tvalue;
	latex->DrawLatex(left,top,"#psi(2S)");
	top -= step;
	latex->DrawLatex(left,top,"CMS Preliminary");
	top -= step;

	latex->DrawLatex(left,top,"pp  #sqrt{s} = 7 TeV");
	c1->SaveAs(Form("%s/BFrac_BPH-10-014_compare_rap2.pdf",savePath.str().c_str()));

	///
	TFile *outfile  = new TFile(Form("%s/BFrac_%dS.root",savePath.str().c_str(),nState-3),"RECREATE");

	graph_BFrac[0][0]->SetName("graph_BFrac_rap1");
	graph_BFrac[0][1]->SetName("graph_BFrac_rap2");
	graph_BFrac[0][2]->SetName("graph_BFrac_rap3");
	graph_FracNonPrompt[0][0]->SetName("graph_FracNonPrompt_rap1");
	graph_FracNonPrompt[0][1]->SetName("graph_FracNonPrompt_rap2");
	graph_FracNonPrompt[0][2]->SetName("graph_FracNonPrompt_rap3");
	graph_FracBkg[0][0]->SetName("graph_FracBkg_rap1");
	graph_FracBkg[0][1]->SetName("graph_FracBkg_rap2");
	graph_FracBkg[0][2]->SetName("graph_FracBkg_rap3");
	graph_BFracRap1->SetName("graph_BPH_10_014_BFrac_rap1");
	graph_BFracRap2->SetName("graph_BPH_10_014_BFrac_rap2");
	graph_BFracRap3->SetName("graph_BPH_10_014_BFrac_rap3");

	graph_BFrac[0][0]->Write();
	graph_BFrac[0][1]->Write();
	graph_BFrac[0][2]->Write();
	graph_FracNonPrompt[0][0]->Write();
	graph_FracNonPrompt[0][1]->Write();
	graph_FracNonPrompt[0][2]->Write();
	graph_FracBkg[0][0]->Write();
	graph_FracBkg[0][1]->Write();
	graph_FracBkg[0][2]->Write();
	graph_BFracRap1->Write();
	graph_BFracRap2->Write();

	graph_BFracRap3->Write();
	outfile->Close();
}

//===============================================
void evaluateCtauCut(double nSigma, int nState, int type, bool doCtauUncer){ // type=0: PR(-ctau,ctau); type=1: NP(ctau,+infinity)
	evaluate(nSigma, nState, type, doCtauUncer); 
	plotEval(nSigma, nState, type);
}

//===============================================
void evaluate(double nSigma, int nState, int type, bool doCtauUncer){ // type=0: PR(-ctau,ctau); type=1: NP(ctau,+infinity)
	cout<<"nSigma: "<<nSigma<<endl;
	cout<<"Psi: "<<nState<<"S"<<endl;
	cout<<"type: "<<type<<endl;

	int RapBins = onia::kNbRapForPTBins,
			PtBins  = onia::kNbPTMaxBins;

	std::stringstream savePath;
	savePath << "Fit/parameter/evaluateCtau";
	if(type==0)
		savePath << "/PR_"<<nSigma<<"sigma";;
	if(type==1)
		savePath << "/NP_"<<nSigma<<"sigma";;
	gSystem->mkdir(savePath.str().c_str(),kTRUE);

	double pT[RapBins][PtBins];
	double FracPR[RapBins][PtBins];
	double FracPRRelativeErr[RapBins][PtBins];
	double FracBG[RapBins][PtBins];
	double FracBGRelativeErr[RapBins][PtBins];
	double FracNP[RapBins][PtBins];
	double FracNPRelativeErr[RapBins][PtBins];
	double FracBGNP[RapBins][PtBins];
	double CtauCut[RapBins][PtBins];
	double PRprob[RapBins][PtBins];
	double evtPR[RapBins][PtBins];
	double evtNP[RapBins][PtBins];
	double evtBG[RapBins][PtBins];
	double sigmaP[RapBins][PtBins];
	double sigmaP_L[RapBins][PtBins];


	double pTErrLow[RapBins][PtBins];
	double pTErrHigh[RapBins][PtBins];
	double FracPRErr[RapBins][PtBins];
	double FracPRErrErr[RapBins][PtBins];
	double FracPRRelativeErrErr[RapBins][PtBins];
	double FracBGErr[RapBins][PtBins];
	double FracBGErrErr[RapBins][PtBins];
	double FracBGRelativeErrErr[RapBins][PtBins];
	double FracNPErr[RapBins][PtBins];
	double FracNPErrErr[RapBins][PtBins];
	double FracNPRelativeErrErr[RapBins][PtBins];
	double FracBGNPErr[RapBins][PtBins];
	double CtauCutErr[RapBins][PtBins];
	double PRprobErr[RapBins][PtBins];

	double evtPRErr[RapBins][PtBins];
	double evtNPErr[RapBins][PtBins];
	double evtBGErr[RapBins][PtBins];
	double sigmaPErr[RapBins][PtBins];
	double sigmaP_LErr[RapBins][PtBins];

	TH1F *histPRFracDist[RapBins][PtBins];
	TH1F *histNPFracDist[RapBins][PtBins];
	TH1F *histBGFracDist[RapBins][PtBins];
	TH1F *histCtFracDist[RapBins][PtBins];

	double PRinNumerator[RapBins][PtBins];
	double NPinNumerator[RapBins][PtBins];
	double BGinNumerator[RapBins][PtBins];
	double PRinNumeratorErr[RapBins][PtBins];
	double NPinNumeratorErr[RapBins][PtBins];
	double BGinNumeratorErr[RapBins][PtBins];

	cout<<"Initializing parameters....."<<endl;
	for(int rapBin = 1; rapBin < RapBins+1; rapBin++){
		for(int ptBin = 1; ptBin < PtBins+1; ptBin++){

			pT[rapBin-1][ptBin-1]                   = 0.;
			pTErrLow[rapBin-1][ptBin-1]             = 0.;
			pTErrHigh[rapBin-1][ptBin-1]            = 0.;

			FracPR[rapBin-1][ptBin-1]               = 0.;
			FracPRRelativeErr[rapBin-1][ptBin-1]    = 0.;
			FracBG[rapBin-1][ptBin-1]               = 0.;
			FracBGRelativeErr[rapBin-1][ptBin-1]    = 0.;
			FracNP[rapBin-1][ptBin-1]               = 0.;
			FracNPRelativeErr[rapBin-1][ptBin-1]    = 0.;
			FracBGNP[rapBin-1][ptBin-1]             = 0.;
			CtauCut[rapBin-1][ptBin-1]              = 0.;
			PRprob[rapBin-1][ptBin-1]               = 0.;
			evtPR[rapBin-1][ptBin-1]                = 0.;
			evtNP[rapBin-1][ptBin-1]                = 0.;
			evtBG[rapBin-1][ptBin-1]                = 0.;
			sigmaP[rapBin-1][ptBin-1]               = 0.;
			sigmaP_L[rapBin-1][ptBin-1]             = 0.;

			FracPRErr[rapBin-1][ptBin-1]            = 0.;

			FracPRErrErr[rapBin-1][ptBin-1]         = 0.;
			FracPRRelativeErrErr[rapBin-1][ptBin-1] = 0.;
			FracBGErr[rapBin-1][ptBin-1]            = 0.;
			FracBGErrErr[rapBin-1][ptBin-1]         = 0.;
			FracBGRelativeErrErr[rapBin-1][ptBin-1] = 0.;
			FracNPErr[rapBin-1][ptBin-1]            = 0.;
			FracNPErrErr[rapBin-1][ptBin-1]         = 0.;
			FracNPRelativeErrErr[rapBin-1][ptBin-1] = 0.;
			FracBGNPErr[rapBin-1][ptBin-1]          = 0.;
			CtauCutErr[rapBin-1][ptBin-1]           = 0.;
			PRprobErr[rapBin-1][ptBin-1]            = 0.;
			evtPRErr[rapBin-1][ptBin-1]             = 0.;
			evtNPErr[rapBin-1][ptBin-1]             = 0.;
			evtBGErr[rapBin-1][ptBin-1]             = 0.;
			sigmaPErr[rapBin-1][ptBin-1]            = 0.;
			sigmaP_LErr[rapBin-1][ptBin-1]          = 0.;

			PRinNumerator[rapBin-1][ptBin-1]        = 0.;
			NPinNumerator[rapBin-1][ptBin-1]        = 0.;
			BGinNumerator[rapBin-1][ptBin-1]        = 0.;
			PRinNumeratorErr[rapBin-1][ptBin-1]     = 0.;
			NPinNumeratorErr[rapBin-1][ptBin-1]     = 0.;
			BGinNumeratorErr[rapBin-1][ptBin-1]     = 0.;
		}
	}
	cout<<"done!"<<endl;

	TFile *inFile;
	RooWorkspace *ws;
	char inName[200];
	RooDataSet *data, *dataSR;
	for(int rapBin = 1; rapBin < RapBins+1; rapBin++){
		for(int ptBin = 1; ptBin < PtBins+1; ptBin++){

			sprintf(inName, Form("tmpFiles/fit_Psi%dS_rap%d_pt%d.root",nState-3,rapBin,ptBin));
			cout<<"inName: "<<inName<<endl;
			inFile=new TFile(inName,"R");
			if(!inFile) continue;
			ws=(RooWorkspace *)inFile->Get("ws_masslifetime");
			if(!ws) { cout<<">>=======Error: no workspace in root file=========="<<endl; continue; }

			double ptMin = onia::pTRange[rapBin][ptBin-1];
			double ptMax = onia::pTRange[rapBin][ptBin];
			pT[rapBin-1][ptBin-1] =  getMeanPt(rapBin,ptBin,inName);
			pTErrLow[rapBin-1][ptBin-1] =  pT[rapBin-1][ptBin-1]-ptMin;
			pTErrHigh[rapBin-1][ptBin-1] =  ptMax-pT[rapBin-1][ptBin-1];

			//////////////////////////////////////////Get Mean  and SigmaWei Value from MassFit////////////////////
			RooRealVar *CBmass=(RooRealVar *)ws->var("CBmass");
			RooRealVar *CBsigma=(RooRealVar *)ws->var("CBsigma");
			RooRealVar *CBsigma2=(RooRealVar *)ws->var("CBsigma2");
			RooRealVar *fracCB1_=(RooRealVar *)ws->var("fracCB1");
			double Mean = CBmass->getVal();
			double MeanErr = CBmass->getError();
			double Sigma = CBsigma->getVal();
			double Sigma2 = CBsigma2->getVal();
			double fracCB1 = fracCB1_->getVal();
			double SigmaWei = sqrt(pow(Sigma,2)*fracCB1+pow(Sigma2,2)*(1-fracCB1));
			double sigMaxMass = Mean+SigmaWei*onia::nSigMass;
			double sigMinMass = Mean-SigmaWei*onia::nSigMass;
			double sbHighMass = Mean+SigmaWei*onia::nSigBkgHigh;
			double sbLowMass =  Mean-SigmaWei*onia::nSigBkgLow;
			stringstream binName, binNameSR, cutSR;
			binName    << "data_rap"<<rapBin<<"_pt"<<ptBin;
			binNameSR  << "data_rap"<<rapBin<<"_pt"<<ptBin<<"_SR";
			cutSR      << "JpsiMass > "<<sigMinMass<<" && JpsiMass < "<<sigMaxMass;
			cout<<"cutSR: "<<cutSR.str().c_str()<<endl;

			RooDataSet *data   = (RooDataSet*)ws->data(binName.str().c_str());
			RooDataSet *dataSR = (RooDataSet*)data->reduce(
					Cut(cutSR.str().c_str()),
					Name(binNameSR.str().c_str()),
					Title("Data For Fitting"));
			int EventsSR = dataSR->numEntries();
			///////////////////////////////////////////////DONE/////////////////////////////////////////////////////

			RooRealVar *fBkg__ = (RooRealVar*)ws->var("fBkg");
			if(!fBkg__) continue;
			RooRealVar *fPrompt__ = (RooRealVar*)ws->var("fPrompt");
			double fracBG = fBkg__->getVal();
			double fracBGErr = fBkg__->getError();
			double fracPR = fPrompt__->getVal();
			double fracPRErr = fPrompt__->getError();
			double fracNP =  1.-fracBG-fracPR;
			double fracNPErr = sqrt(pow(fracBGErr,2)+pow(fracPRErr,2));

			RooRealVar Jpsict(*ws->var("Jpsict"));
			RooRealVar JpsictErr(*ws->var("JpsictErr"));
			Jpsict.setMin(ctMin);   Jpsict.setMax(ctMax);
			cout<<"Jpsict.getMin(): "<<Jpsict.getMin()<<endl;
			cout<<"Jpsict.getMax(): "<<Jpsict.getMax()<<endl;

			int NumEvt = 0., maxEvt = 100000;
			if(dataSR->numEntries() < maxEvt) NumEvt = dataSR->numEntries();
			else NumEvt = maxEvt;
			//Jpsict error distributions, asumming that it is same for BG, P, and NP
			RooDataSet *dataJpsictErr = (RooDataSet*)dataSR->reduce(SelectVars(RooArgSet(JpsictErr)),
					EventRange(0,NumEvt),Name("dataJpsictErr"));

			////////////////calculate sigma of prompt p.d.f
			//RooDataSet *dataSRFit=(RooDataSet *)ws->data(Form("data_rap%d_pt%d_SR",rapBin,ptBin));
			//int NumEvt = 0., maxEvt = 100000;
			//if(dataSRFit->numEntries() < maxEvt) NumEvt = dataSRFit->numEntries();
			//else NumEvt = maxEvt;
			//RooDataSet *dataJpsictErrFit = (RooDataSet*)dataSRFit->reduce(SelectVars(RooArgSet(JpsictErr)),
			//		EventRange(0,NumEvt),Name("dataJpsictErr"));

			vector<double> SigmaP;
			SigmaP = getSigma(ws, dataJpsictErr, rapBin, ptBin);

			sigmaP[rapBin-1][ptBin-1]    = SigmaP[0];
			sigmaPErr[rapBin-1][ptBin-1] = SigmaP[1];

			//sigma of pseudo-proper decay length
			sigmaP_L[rapBin-1][ptBin-1]    = SigmaP[0] * pT[rapBin-1][ptBin-1] / Mean;
			sigmaP_LErr[rapBin-1][ptBin-1] = sigmaP_L[rapBin-1][ptBin-1] *
				sqrt(pow(sigmaPErr[rapBin-1][ptBin-1]/sigmaP[rapBin-1][ptBin-1],2) + pow(MeanErr/Mean,2));

			//define function y = a + b * pT
			double a = 0.073, b = 0.0027;
			//proper decay length
			double L_decay = a + b * pT[rapBin-1][ptBin-1];
			//pseudo-proper decay length
			double l_pdecay =0., scale = 1.;
			if(nState==4) l_pdecay = L_decay * Mean / pT[rapBin-1][ptBin-1];
			else if(nState==5){
				//if(rapBin == 1) scale = 1.23; //1.21;
				//if(rapBin == 2) scale = 1.31; //1.29;
				//if(rapBin == 3) scale = 1.32; //1.31; 
				l_pdecay = scale * L_decay * 3.092 / pT[rapBin-1][ptBin-1];
			}
			cout<<"l_pdecay: "<<l_pdecay<<endl;
			//sigmaP[rapBin-1][ptBin-1] = l_pdecay;
			//sigmaPErr[rapBin-1][ptBin-1] = 0;
			////////////////////////////////////////////////

			//ctCut define as nSigma * sigma
			double ctCut = nSigma*l_pdecay;
			if(nSigma==-1) ctCut = 0.1; //mm

			double ctCutMin=0., ctCutMax=0.;
			if(type==0) {ctCutMin = -ctCut; ctCutMax = ctCut;}
			if(type==1) {ctCutMin = ctCut; ctCutMax = 6.;}
			cout<<"ctCutMin: "<<ctCutMin<<endl;
			cout<<"ctCutMax: "<<ctCutMax<<endl;

			CtauCut[rapBin-1][ptBin-1] = ctCut;
			CtauCutErr[rapBin-1][ptBin-1] = 0.;

			vector<double> InteRlts = calculateInte(ws,dataJpsictErr,ctCutMin,ctCutMax);
			double IntePR    = InteRlts[0];
			double InteNP    = InteRlts[1];
			double InteBG    = InteRlts[2];

			double numeratorPR = IntePR * fracPR,
						 numeratorNP = InteNP * fracNP,
						 numeratorBG = InteBG * fracBG;
			double denominator = IntePR * fracPR + InteNP * fracNP + InteBG * fracBG;

			FracPR[rapBin-1][ptBin-1]     = numeratorPR  / denominator ;
			FracNP[rapBin-1][ptBin-1]     = numeratorNP  / denominator ;
			FracBG[rapBin-1][ptBin-1]     = numeratorBG  / denominator ;
			FracBGNP[rapBin-1][ptBin-1]   = ( numeratorNP + numeratorBG ) / denominator ;
			//error on "denominator" is neglected, or there should be big bias in the fit.
			//we just need to evaluate uncertainties of numeratorPR, numeratorNP, numeratorBG

			//cout<<"fracPR: "<<fracPR<<endl;
			//cout<<"fracNP: "<<fracNP<<endl;
			//cout<<"fracBG: "<<fracBG<<endl;
			//cout<<"IntePR: "<<IntePR<<endl;
			//cout<<"InteNP: "<<InteNP<<endl;
			//cout<<"InteBG: "<<InteBG<<endl;
			//cout<<"FracPR: "<<FracPR[rapBin-1][ptBin-1]<<endl;
			//cout<<"FracNP: "<<FracNP[rapBin-1][ptBin-1]<<endl;
			//cout<<"FracBG: "<<FracBG[rapBin-1][ptBin-1]<<endl;

			PRinNumerator[rapBin-1][ptBin-1]    = numeratorPR;
			NPinNumerator[rapBin-1][ptBin-1]    = numeratorNP;
			BGinNumerator[rapBin-1][ptBin-1]    = numeratorBG;

			evtPR[rapBin-1][ptBin-1]     = numeratorPR*EventsSR;
			evtNP[rapBin-1][ptBin-1]     = numeratorNP*EventsSR;
			evtBG[rapBin-1][ptBin-1]     = numeratorBG*EventsSR;

			PRprob[rapBin-1][ptBin-1] = IntePR;
			PRprobErr[rapBin-1][ptBin-1] = 0.;

			if(doCtauUncer){
				//////***************evaluate Uncertainties**********************
				cout<<">>rap "<<rapBin<<" pt "<<ptBin<<" : evaluating Uncertainties........."<<endl;

				double rangeWindowPR = 0.02;
				double rangeWindowNP = 0.01;
				double rangeWindowBG = 0.02;
				double rangeWindowCt = 0.05;
				histPRFracDist[rapBin-1][ptBin-1] = new TH1F(Form("histPRFracDist_rap%d_pt%d",rapBin,ptBin),"",
						50,numeratorPR-rangeWindowPR,numeratorPR+rangeWindowPR);
				histNPFracDist[rapBin-1][ptBin-1] = new TH1F(Form("histNPFracDist_rap%d_pt%d",rapBin,ptBin),"",
						50,numeratorNP-rangeWindowNP, numeratorNP+rangeWindowNP);
				histBGFracDist[rapBin-1][ptBin-1] = new TH1F(Form("histBGFracDist_rap%d_pt%d",rapBin,ptBin),"",
						50,numeratorBG-rangeWindowBG, numeratorBG+rangeWindowBG);
				histCtFracDist[rapBin-1][ptBin-1] = new TH1F(Form("histCtFracDist_rap%d_pt%d",rapBin,ptBin),"",
						50,sigmaP[rapBin-1][ptBin-1]-rangeWindowCt, sigmaP[rapBin-1][ptBin-1]+rangeWindowCt);
				histPRFracDist[rapBin-1][ptBin-1]->SetXTitle("Prompt fraction");
				histNPFracDist[rapBin-1][ptBin-1]->SetXTitle("Non-prompt fraction");
				histCtFracDist[rapBin-1][ptBin-1]->SetXTitle("Ct RMS");

				RooFitResult* fitRlt = (RooFitResult*)ws->obj(Form("l_fitresult_rap%d_pt%d",rapBin,ptBin));
				if(!fitRlt) { cout<<">>=======Error: no Fit result object in workspace=========="<<endl; continue; }
				fitRlt->Print();
				int nEvents = 50; //200
				//if(nState==4) nEvents = 50; //100
				//-------------parameters---------------------
				//  bkgTauDSD        3.8546e-03 +/-  9.93e-04
				//  bkgTauFD         1.0575e-01 +/-  1.80e-02
				//  bkgTauSSD_SBL    3.9861e-01 +/-  1.53e-02
				//  bkgTauSSD_SBR    3.5718e-01 +/-  2.32e-02
				//  fBkg             3.5902e-01 +/-  7.51e-03
				//  fBkgDSD_SBL      4.1656e-01 +/-  1.95e-02
				//  fBkgDSD_SBR      6.7681e-01 +/-  2.14e-02
				//  fBkgSBL          9.6721e-01 +/-  5.53e-03 //new
				//  fBkgSBR          9.9624e-01 +/-  5.67e-03 //new
				//  fBkgSSDR_SBL     5.6511e-01 +/-  1.81e-02
				//  fBkgSSDR_SBR     2.9558e-01 +/-  1.83e-02
				//  fPrompt          2.2219e-01 +/-  9.67e-03
				//  fracGauss2       1.0377e-01 +/-  2.09e-02
				//  nonPromptTau     4.0834e-01 +/-  1.22e-02

				RooRealVar bkgTauDSD(*ws->var("bkgTauDSD"));
				RooRealVar bkgTauFD(*ws->var("bkgTauFD"));
				RooRealVar bkgTauSSD_SBL(*ws->var("bkgTauSSD_SBL"));
				RooRealVar bkgTauSSD_SBR(*ws->var("bkgTauSSD_SBR"));
				RooRealVar fBkgDSD_SBL(*ws->var("fBkgDSD_SBL"));
				RooRealVar fBkgDSD_SBR(*ws->var("fBkgDSD_SBR"));
				RooRealVar fBkgSBL(*ws->var("fBkgSBL"));
				RooRealVar fBkgSBR(*ws->var("fBkgSBR"));
				RooRealVar fBkgSSDR_SBL(*ws->var("fBkgSSDR_SBL")); // fixed when pT bin > 7 for 1S
				RooRealVar fBkgSSDR_SBR(*ws->var("fBkgSSDR_SBR")); // fixed when pT bin > 7 for 1S
				RooRealVar fPrompt(*ws->var("fPrompt"));
				RooRealVar fBkg(*ws->var("fBkg"));
				RooRealVar fracGauss2(*ws->var("fracGauss2"));
				RooRealVar nonPromptTau(*ws->var("nonPromptTau"));

				RooArgSet *paraVars = new RooArgSet(bkgTauDSD,bkgTauFD,bkgTauSSD_SBL,bkgTauSSD_SBR,
						fBkgDSD_SBL,fBkgDSD_SBR, fBkgSSDR_SBL,fBkgSSDR_SBR,
						nonPromptTau);
				paraVars->add(RooArgSet(fPrompt,fBkg,fracGauss2,fBkgSBL,fBkgSBR));

				RooAbsPdf *multiVarPdf = (RooAbsPdf*)fitRlt->createHessePdf(*paraVars);
				RooDataSet *multiVarData = (RooDataSet*)multiVarPdf->generate(*paraVars,nEvents);
				for(int n = 0; n < nEvents; n++) {
					if(n%5==0) cout<<(double)n/nEvents*100<<"%"<<endl;
					RooArgSet* args = (RooArgSet*)(multiVarData->get(n));

					double bkgTauFD_ = ((RooRealVar*)args->find("bkgTauFD"))->getVal();
					ws->var("bkgTauFD")->setVal(bkgTauFD_);
					double bkgTauSSD_SBL_ = ((RooRealVar*)args->find("bkgTauSSD_SBL"))->getVal();
					ws->var("bkgTauSSD_SBL")->setVal(bkgTauSSD_SBL_);
					double bkgTauSSD_SBR_ = ((RooRealVar*)args->find("bkgTauSSD_SBR"))->getVal();
					ws->var("bkgTauSSD_SBR")->setVal(bkgTauSSD_SBR_);
					double fBkgDSD_SBL_ = ((RooRealVar*)args->find("fBkgDSD_SBL"))->getVal();
					ws->var("fBkgDSD_SBL")->setVal(fBkgDSD_SBL_);
					double fBkgSBL_ = ((RooRealVar*)args->find("fBkgSBL"))->getVal();
					ws->var("fBkgSBL")->setVal(fBkgSBL_);
					double fBkgSBR_ = ((RooRealVar*)args->find("fBkgSBR"))->getVal();
					ws->var("fBkgSBR")->setVal(fBkgSBR_);
					double fBkgDSD_SBR_ = ((RooRealVar*)args->find("fBkgDSD_SBR"))->getVal();
					ws->var("fBkgDSD_SBR")->setVal(fBkgDSD_SBR_);
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

					fracPR = fPrompt_; fracBG = fBkg_; fracNP =  1.-fracBG-fracPR;
					vector<double> InteRltsTemp = calculateInte(ws,dataJpsictErr,ctCutMin,ctCutMax);
					IntePR    = InteRltsTemp[0];
					InteNP    = InteRltsTemp[1];
					InteBG    = InteRltsTemp[2];
					numeratorPR = IntePR * fracPR;
					numeratorNP = InteNP * fracNP;
					numeratorBG = InteBG * fracBG;

					histPRFracDist[rapBin-1][ptBin-1] -> Fill(numeratorPR);
					histNPFracDist[rapBin-1][ptBin-1] -> Fill(numeratorNP);
					histBGFracDist[rapBin-1][ptBin-1] -> Fill(numeratorBG);

					vector<double> SigmaP_Temp;
					SigmaP_Temp = getSigma(ws, dataJpsictErr, rapBin, ptBin);
					histCtFracDist[rapBin-1][ptBin-1] -> Fill(SigmaP_Temp[0]);
				}

				PRinNumeratorErr[rapBin-1][ptBin-1] =  histPRFracDist[rapBin-1][ptBin-1] -> GetRMS();
				NPinNumeratorErr[rapBin-1][ptBin-1] =  histNPFracDist[rapBin-1][ptBin-1] -> GetRMS();
				BGinNumeratorErr[rapBin-1][ptBin-1] =  histBGFracDist[rapBin-1][ptBin-1] -> GetRMS();

				FracPRErr[rapBin-1][ptBin-1]       =    PRinNumeratorErr[rapBin-1][ptBin-1] / denominator;
				FracNPErr[rapBin-1][ptBin-1]       =    NPinNumeratorErr[rapBin-1][ptBin-1] / denominator;
				FracBGErr[rapBin-1][ptBin-1]       =    BGinNumeratorErr[rapBin-1][ptBin-1] / denominator;
				FracBGNPErr[rapBin-1][ptBin-1]     = 0.;

				FracPRRelativeErr[rapBin-1][ptBin-1] = FracPRErr[rapBin-1][ptBin-1] / FracPR[rapBin-1][ptBin-1] *100;
				FracNPRelativeErr[rapBin-1][ptBin-1] = FracNPErr[rapBin-1][ptBin-1] / FracNP[rapBin-1][ptBin-1] *100;
				FracBGRelativeErr[rapBin-1][ptBin-1] = FracBGErr[rapBin-1][ptBin-1] / FracBG[rapBin-1][ptBin-1] *100;

				sigmaPErr[rapBin-1][ptBin-1] = histCtFracDist[rapBin-1][ptBin-1] -> GetRMS();
				////done!!!!
			}

			sigmaP_L[rapBin-1][ptBin-1]    = SigmaP[0] * pT[rapBin-1][ptBin-1] / Mean;
			sigmaP_LErr[rapBin-1][ptBin-1] = sigmaP_L[rapBin-1][ptBin-1] *
				sqrt(pow(sigmaPErr[rapBin-1][ptBin-1]/sigmaP[rapBin-1][ptBin-1],2) + pow(MeanErr/Mean,2));

		}
	}


	for(int rapBin = 1; rapBin < RapBins+1; rapBin++){
		for(int ptBin = 1; ptBin < PtBins+1; ptBin++){
			cout<<"rap "<<rapBin<<"   pt "<<ptBin<<endl;
			cout<<"------------------------------------"<<endl;
			cout<<"PRinNumerator: "<<PRinNumerator[rapBin-1][ptBin-1]<<"   PRinNumeratorErr: "<<PRinNumeratorErr[rapBin-1][ptBin-1]<<endl;
			cout<<"NPinNumerator: "<<NPinNumerator[rapBin-1][ptBin-1]<<"   NPinNumeratorErr: "<<NPinNumeratorErr[rapBin-1][ptBin-1]<<endl;
			cout<<"BGinNumerator: "<<BGinNumerator[rapBin-1][ptBin-1]<<"   BGinNumeratorErr: "<<BGinNumeratorErr[rapBin-1][ptBin-1]<<endl;
			cout<<"------------------------------------"<<endl;
			cout<<"FracPR:   "<<FracPR[rapBin-1][ptBin-1]<<  "  FracPRErr: "<<FracPRErr[rapBin-1][ptBin-1]<<endl;
			cout<<"FracBG:   "<<FracBG[rapBin-1][ptBin-1]<<  "  FracBGErr: "<<FracBGErr[rapBin-1][ptBin-1]<<endl;
			cout<<"FracNP:   "<<FracNP[rapBin-1][ptBin-1]<<  "  FracNPErr: "<<FracNPErr[rapBin-1][ptBin-1]<<endl;
			cout<<"FracBGNP: "<<FracBGNP[rapBin-1][ptBin-1]<<"  FracBGNPErr: "<<FracBGNPErr[rapBin-1][ptBin-1]<<endl;
			cout<<"CtauCut:  "<<CtauCut[rapBin-1][ptBin-1]<< "  CtauCutErr: "<<CtauCutErr[rapBin-1][ptBin-1]<<endl;
			cout<<"PRprob:   "<<PRprob[rapBin-1][ptBin-1]<<  "  PRprobErr: "<<PRprobErr[rapBin-1][ptBin-1]<<endl;
			cout<<"===================================="<<endl;
		}
	}

	double Ymin = 0., Ymax = 1.2;
	double Xmin = 0., Xmax = 100.;

	TGraphAsymmErrors *graph_FracPR[RapBins], *graph_FracBG[RapBins],
										*graph_FracNP[RapBins], *graph_FracBGNP[RapBins],
										*graph_FracPRErr[RapBins], *graph_FracBGErr[RapBins],
										*graph_FracNPErr[RapBins],
										*graph_FracPRRelativeErr[RapBins], *graph_FracBGRelativeErr[RapBins],
										*graph_FracNPRelativeErr[RapBins],
										*graph_evtPR[RapBins], *graph_evtNP[RapBins], *graph_evtBG[RapBins],
										*graph_sigmaP[RapBins], *graph_sigmaP_L[RapBins],
										*graph_CtauCut[RapBins], *graph_PRprob[RapBins];
	for(int rapBin = 1; rapBin < RapBins+1; rapBin++){
		graph_FracPR[rapBin-1] = new TGraphAsymmErrors(PtBins, pT[rapBin-1], FracPR[rapBin-1],
				pTErrLow[rapBin-1], pTErrHigh[rapBin-1], FracPRErr[rapBin-1], FracPRErr[rapBin-1]);
		graph_FracPRErr[rapBin-1] = new TGraphAsymmErrors(PtBins, pT[rapBin-1], FracPRErr[rapBin-1],
				pTErrLow[rapBin-1], pTErrHigh[rapBin-1], FracPRErrErr[rapBin-1], FracPRErrErr[rapBin-1]);
		graph_FracPRRelativeErr[rapBin-1] = new TGraphAsymmErrors(PtBins, pT[rapBin-1], FracPRRelativeErr[rapBin-1],
				pTErrLow[rapBin-1], pTErrHigh[rapBin-1], FracPRRelativeErrErr[rapBin-1], FracPRRelativeErrErr[rapBin-1]);
		graph_evtPR[rapBin-1] = new TGraphAsymmErrors(PtBins, pT[rapBin-1], evtPR[rapBin-1],
				pTErrLow[rapBin-1], pTErrHigh[rapBin-1], evtPRErr[rapBin-1], evtPRErr[rapBin-1]);

		graph_FracBG[rapBin-1] = new TGraphAsymmErrors(PtBins, pT[rapBin-1], FracBG[rapBin-1],
				pTErrLow[rapBin-1], pTErrHigh[rapBin-1], FracBGErr[rapBin-1], FracBGErr[rapBin-1]);
		graph_FracBGErr[rapBin-1] = new TGraphAsymmErrors(PtBins, pT[rapBin-1], FracBGErr[rapBin-1],
				pTErrLow[rapBin-1], pTErrHigh[rapBin-1], FracBGErrErr[rapBin-1], FracBGErrErr[rapBin-1]);
		graph_FracBGRelativeErr[rapBin-1] = new TGraphAsymmErrors(PtBins, pT[rapBin-1], FracBGRelativeErr[rapBin-1],
				pTErrLow[rapBin-1], pTErrHigh[rapBin-1], FracBGRelativeErrErr[rapBin-1], FracBGRelativeErrErr[rapBin-1]);
		graph_evtBG[rapBin-1] = new TGraphAsymmErrors(PtBins, pT[rapBin-1], evtBG[rapBin-1],
				pTErrLow[rapBin-1], pTErrHigh[rapBin-1], evtBGErr[rapBin-1], evtBGErr[rapBin-1]);

		graph_sigmaP[rapBin-1] = new TGraphAsymmErrors(PtBins, pT[rapBin-1], sigmaP[rapBin-1],
				pTErrLow[rapBin-1], pTErrHigh[rapBin-1], sigmaPErr[rapBin-1], sigmaPErr[rapBin-1]);

		graph_sigmaP_L[rapBin-1] = new TGraphAsymmErrors(PtBins, pT[rapBin-1], sigmaP_L[rapBin-1],
				pTErrLow[rapBin-1], pTErrHigh[rapBin-1], sigmaP_LErr[rapBin-1], sigmaP_LErr[rapBin-1]);

		graph_FracNP[rapBin-1] = new TGraphAsymmErrors(PtBins, pT[rapBin-1], FracNP[rapBin-1],
				pTErrLow[rapBin-1], pTErrHigh[rapBin-1], FracNPErr[rapBin-1], FracNPErr[rapBin-1]);
		graph_FracNPErr[rapBin-1] = new TGraphAsymmErrors(PtBins, pT[rapBin-1], FracNPErr[rapBin-1],
				pTErrLow[rapBin-1], pTErrHigh[rapBin-1], FracNPErrErr[rapBin-1], FracNPErrErr[rapBin-1]);
		graph_FracNPRelativeErr[rapBin-1] = new TGraphAsymmErrors(PtBins, pT[rapBin-1], FracNPRelativeErr[rapBin-1],
				pTErrLow[rapBin-1], pTErrHigh[rapBin-1], FracNPRelativeErrErr[rapBin-1], FracNPRelativeErrErr[rapBin-1]);
		graph_evtNP[rapBin-1] = new TGraphAsymmErrors(PtBins, pT[rapBin-1], evtNP[rapBin-1],
				pTErrLow[rapBin-1], pTErrHigh[rapBin-1], evtNPErr[rapBin-1], evtNPErr[rapBin-1]);

		graph_FracBGNP[rapBin-1] = new TGraphAsymmErrors(PtBins, pT[rapBin-1], FracBGNP[rapBin-1],
				pTErrLow[rapBin-1], pTErrHigh[rapBin-1], FracBGNPErr[rapBin-1], FracBGNPErr[rapBin-1]);
		graph_CtauCut[rapBin-1] = new TGraphAsymmErrors(PtBins, pT[rapBin-1], CtauCut[rapBin-1],
				pTErrLow[rapBin-1], pTErrHigh[rapBin-1], CtauCutErr[rapBin-1], CtauCutErr[rapBin-1]);
		graph_PRprob[rapBin-1] = new TGraphAsymmErrors(PtBins, pT[rapBin-1], PRprob[rapBin-1],
				pTErrLow[rapBin-1], pTErrHigh[rapBin-1], PRprobErr[rapBin-1], PRprobErr[rapBin-1]);

		graph_FracPR[rapBin-1]->SetName(Form("graph_FracPR_%d",rapBin-1));
		graph_FracBG[rapBin-1]->SetName(Form("graph_FracBG_%d",rapBin-1));
		graph_FracNP[rapBin-1]->SetName(Form("graph_FracNP_%d",rapBin-1));
		graph_FracPRErr[rapBin-1]->SetName(Form("graph_FracPRErr_%d",rapBin-1));
		graph_FracPRRelativeErr[rapBin-1]->SetName(Form("graph_FracPRRelativeErr_%d",rapBin-1));
		graph_FracBGErr[rapBin-1]->SetName(Form("graph_FracBGErr_%d",rapBin-1));
		graph_FracBGRelativeErr[rapBin-1]->SetName(Form("graph_FracBGRelativeErr_%d",rapBin-1));

		graph_FracNPErr[rapBin-1]->SetName(Form("graph_FracNPErr_%d",rapBin-1));
		graph_FracNPRelativeErr[rapBin-1]->SetName(Form("graph_FracNPRelativeErr_%d",rapBin-1));
		graph_FracBGNP[rapBin-1]->SetName(Form("graph_FracBGNP_%d",rapBin-1));
		graph_evtPR[rapBin-1]->SetName(Form("graph_evtPR_%d",rapBin-1));
		graph_evtNP[rapBin-1]->SetName(Form("graph_evtNP_%d",rapBin-1));
		graph_evtBG[rapBin-1]->SetName(Form("graph_evtBG_%d",rapBin-1));
		graph_sigmaP[rapBin-1]->SetName(Form("graph_sigmaP_%d",rapBin-1));
		graph_sigmaP_L[rapBin-1]->SetName(Form("graph_sigmaP_L_%d",rapBin-1));

		graph_CtauCut[rapBin-1]->SetName(Form("graph_CtauCut_%d",rapBin-1));
		graph_PRprob[rapBin-1]->SetName(Form("graph_PRprob_%d",rapBin-1));

		graph_FracPR[rapBin-1]->SetTitle("");
		graph_FracPR[rapBin-1]->GetXaxis()->SetTitle("p_{T} (GeV)");
		graph_FracPR[rapBin-1]->GetYaxis()->SetTitle("fraction");
		graph_FracPR[rapBin-1]->GetYaxis()->SetRangeUser(Ymin, Ymax);
		graph_FracPR[rapBin-1]->GetXaxis()->SetLimits(Xmin, Xmax);
		graph_FracPR[rapBin-1]->SetMarkerStyle(onia::marker_rapForPTBins[0]);
		graph_FracPR[rapBin-1]->SetMarkerColor(onia::colour_rapForPTBins[0]);
		graph_FracPR[rapBin-1]->SetLineColor(onia::colour_rapForPTBins[0]);

		graph_FracPRErr[rapBin-1]->SetTitle("");
		graph_FracPRErr[rapBin-1]->GetXaxis()->SetTitle("p_{T} (GeV)");
		graph_FracPRErr[rapBin-1]->GetYaxis()->SetTitle("fraction Error");
		graph_FracPRErr[rapBin-1]->GetYaxis()->SetRangeUser(0., 0.03);
		graph_FracPRErr[rapBin-1]->GetXaxis()->SetLimits(Xmin, Xmax);
		graph_FracPRErr[rapBin-1]->SetMarkerStyle(onia::marker_rapForPTBins[0]);
		graph_FracPRErr[rapBin-1]->SetMarkerColor(onia::colour_rapForPTBins[0]);
		graph_FracPRErr[rapBin-1]->SetLineColor(onia::colour_rapForPTBins[0]);

		graph_FracPRRelativeErr[rapBin-1]->SetTitle("");
		graph_FracPRRelativeErr[rapBin-1]->GetXaxis()->SetTitle("p_{T} (GeV)");
		graph_FracPRRelativeErr[rapBin-1]->GetYaxis()->SetTitle("fraction RelativeError (%)");
		graph_FracPRRelativeErr[rapBin-1]->GetYaxis()->SetRangeUser(0., 50.);
		graph_FracPRRelativeErr[rapBin-1]->GetXaxis()->SetLimits(Xmin, Xmax);
		graph_FracPRRelativeErr[rapBin-1]->SetMarkerStyle(onia::marker_rapForPTBins[0]);
		graph_FracPRRelativeErr[rapBin-1]->SetMarkerColor(onia::colour_rapForPTBins[0]);
		graph_FracPRRelativeErr[rapBin-1]->SetLineColor(onia::colour_rapForPTBins[0]);

		graph_FracNP[rapBin-1]->SetTitle("");
		graph_FracNP[rapBin-1]->GetXaxis()->SetTitle("p_{T} (GeV)");
		graph_FracNP[rapBin-1]->GetYaxis()->SetTitle("fraction");
		graph_FracNP[rapBin-1]->GetYaxis()->SetRangeUser(Ymin, Ymax);
		graph_FracNP[rapBin-1]->GetXaxis()->SetLimits(Xmin, Xmax);
		graph_FracNP[rapBin-1]->SetMarkerStyle(onia::marker_rapForPTBins[1]);
		graph_FracNP[rapBin-1]->SetMarkerColor(onia::colour_rapForPTBins[1]);
		graph_FracNP[rapBin-1]->SetLineColor(onia::colour_rapForPTBins[1]);

		graph_FracNPErr[rapBin-1]->SetTitle("");
		graph_FracNPErr[rapBin-1]->GetXaxis()->SetTitle("p_{T} (GeV)");
		graph_FracNPErr[rapBin-1]->GetYaxis()->SetTitle("fraction Error");
		graph_FracNPErr[rapBin-1]->GetYaxis()->SetRangeUser(0., 0.01);
		graph_FracNPErr[rapBin-1]->GetXaxis()->SetLimits(Xmin, Xmax);
		graph_FracNPErr[rapBin-1]->SetMarkerStyle(onia::marker_rapForPTBins[1]);
		graph_FracNPErr[rapBin-1]->SetMarkerColor(onia::colour_rapForPTBins[1]);
		graph_FracNPErr[rapBin-1]->SetLineColor(onia::colour_rapForPTBins[1]);

		graph_FracNPRelativeErr[rapBin-1]->SetTitle("");
		graph_FracNPRelativeErr[rapBin-1]->GetXaxis()->SetTitle("p_{T} (GeV)");

		graph_FracNPRelativeErr[rapBin-1]->GetYaxis()->SetTitle("fraction RelativeError (%)");
		graph_FracNPRelativeErr[rapBin-1]->GetYaxis()->SetRangeUser(0., 50.);
		graph_FracNPRelativeErr[rapBin-1]->GetXaxis()->SetLimits(Xmin, Xmax);
		graph_FracNPRelativeErr[rapBin-1]->SetMarkerStyle(onia::marker_rapForPTBins[1]);
		graph_FracNPRelativeErr[rapBin-1]->SetMarkerColor(onia::colour_rapForPTBins[1]);
		graph_FracNPRelativeErr[rapBin-1]->SetLineColor(onia::colour_rapForPTBins[1]);

		graph_FracBG[rapBin-1]->SetTitle("");
		graph_FracBG[rapBin-1]->GetXaxis()->SetTitle("p_{T} (GeV)");
		graph_FracBG[rapBin-1]->GetYaxis()->SetTitle("fraction");
		graph_FracBG[rapBin-1]->GetYaxis()->SetRangeUser(Ymin, Ymax);
		graph_FracBG[rapBin-1]->GetXaxis()->SetLimits(Xmin, Xmax);
		graph_FracBG[rapBin-1]->SetMarkerStyle(onia::marker_rapForPTBins[2]);
		graph_FracBG[rapBin-1]->SetMarkerColor(onia::colour_rapForPTBins[2]);
		graph_FracBG[rapBin-1]->SetLineColor(onia::colour_rapForPTBins[2]);

		graph_FracBGErr[rapBin-1]->SetTitle("");
		graph_FracBGErr[rapBin-1]->GetXaxis()->SetTitle("p_{T} (GeV)");
		graph_FracBGErr[rapBin-1]->GetYaxis()->SetTitle("fraction Error");
		graph_FracBGErr[rapBin-1]->GetYaxis()->SetRangeUser(0., 0.01);
		graph_FracBGErr[rapBin-1]->GetXaxis()->SetLimits(Xmin, Xmax);
		graph_FracBGErr[rapBin-1]->SetMarkerStyle(onia::marker_rapForPTBins[2]);
		graph_FracBGErr[rapBin-1]->SetMarkerColor(onia::colour_rapForPTBins[2]);
		graph_FracBGErr[rapBin-1]->SetLineColor(onia::colour_rapForPTBins[2]);

		graph_FracBGRelativeErr[rapBin-1]->SetTitle("");
		graph_FracBGRelativeErr[rapBin-1]->GetXaxis()->SetTitle("p_{T} (GeV)");
		graph_FracBGRelativeErr[rapBin-1]->GetYaxis()->SetTitle("fraction RelativeError (%)");
		graph_FracBGRelativeErr[rapBin-1]->GetYaxis()->SetRangeUser(0., 50.);
		graph_FracBGRelativeErr[rapBin-1]->GetXaxis()->SetLimits(Xmin, Xmax);
		graph_FracBGRelativeErr[rapBin-1]->SetMarkerStyle(onia::marker_rapForPTBins[2]);
		graph_FracBGRelativeErr[rapBin-1]->SetMarkerColor(onia::colour_rapForPTBins[2]);
		graph_FracBGRelativeErr[rapBin-1]->SetLineColor(onia::colour_rapForPTBins[2]);

		graph_FracBGNP[rapBin-1]->SetTitle("");
		graph_FracBGNP[rapBin-1]->GetXaxis()->SetTitle("p_{T} (GeV)");
		graph_FracBGNP[rapBin-1]->GetYaxis()->SetTitle("fraction");
		graph_FracBGNP[rapBin-1]->GetYaxis()->SetRangeUser(Ymin, Ymax);
		graph_FracBGNP[rapBin-1]->GetXaxis()->SetLimits(Xmin, Xmax);
		graph_FracBGNP[rapBin-1]->SetMarkerStyle(onia::marker_rapForPTBins[3]);
		graph_FracBGNP[rapBin-1]->SetMarkerColor(onia::colour_rapForPTBins[3]);
		graph_FracBGNP[rapBin-1]->SetLineColor(onia::colour_rapForPTBins[3]);

		graph_evtPR[rapBin-1]->SetTitle("");
		graph_evtPR[rapBin-1]->GetXaxis()->SetTitle("p_{T} (GeV)");
		graph_evtPR[rapBin-1]->GetYaxis()->SetTitle("number of events P in PR");
		graph_evtPR[rapBin-1]->GetYaxis()->SetRangeUser(8., 1.e6+100);
		graph_evtPR[rapBin-1]->GetXaxis()->SetLimits(Xmin, Xmax);
		graph_evtPR[rapBin-1]->SetMarkerStyle(onia::marker_rapForPTBins[rapBin+1]);
		graph_evtPR[rapBin-1]->SetMarkerColor(onia::colour_rapForPTBins[rapBin+1]);
		graph_evtPR[rapBin-1]->SetLineColor(onia::colour_rapForPTBins[rapBin+1]);

		graph_evtNP[rapBin-1]->SetTitle("");
		graph_evtNP[rapBin-1]->GetXaxis()->SetTitle("p_{T} (GeV)");
		graph_evtNP[rapBin-1]->GetYaxis()->SetTitle("number of events NP in PR");
		graph_evtNP[rapBin-1]->GetYaxis()->SetRangeUser(8., 1.e6+100);
		graph_evtNP[rapBin-1]->GetXaxis()->SetLimits(Xmin, Xmax);
		graph_evtNP[rapBin-1]->SetMarkerStyle(onia::marker_rapForPTBins[rapBin+1]);
		graph_evtNP[rapBin-1]->SetMarkerColor(onia::colour_rapForPTBins[rapBin+1]);
		graph_evtNP[rapBin-1]->SetLineColor(onia::colour_rapForPTBins[rapBin+1]);

		graph_evtBG[rapBin-1]->SetTitle("");
		graph_evtBG[rapBin-1]->GetXaxis()->SetTitle("p_{T} (GeV)");
		graph_evtBG[rapBin-1]->GetYaxis()->SetTitle("number of events BG in PR");
		graph_evtBG[rapBin-1]->GetYaxis()->SetRangeUser(8., 1.e6+100);
		graph_evtBG[rapBin-1]->GetXaxis()->SetLimits(Xmin, Xmax);
		graph_evtBG[rapBin-1]->SetMarkerStyle(onia::marker_rapForPTBins[rapBin+1]);
		graph_evtBG[rapBin-1]->SetMarkerColor(onia::colour_rapForPTBins[rapBin+1]);
		graph_evtBG[rapBin-1]->SetLineColor(onia::colour_rapForPTBins[rapBin+1]);

		graph_sigmaP[rapBin-1]->SetTitle("");
		graph_sigmaP[rapBin-1]->GetXaxis()->SetTitle("p_{T} (GeV)");
		graph_sigmaP[rapBin-1]->GetYaxis()->SetTitle("RMS of Prompt (mm)");
		graph_sigmaP[rapBin-1]->GetYaxis()->SetRangeUser(0., 0.06);
		graph_sigmaP[rapBin-1]->GetXaxis()->SetLimits(Xmin, Xmax);
		graph_sigmaP[rapBin-1]->SetMarkerStyle(onia::marker_rapForPTBins[rapBin+1]);
		graph_sigmaP[rapBin-1]->SetMarkerColor(onia::colour_rapForPTBins[rapBin+1]);
		graph_sigmaP[rapBin-1]->SetLineColor(onia::colour_rapForPTBins[rapBin+1]);

		graph_sigmaP_L[rapBin-1]->SetTitle("");
		graph_sigmaP_L[rapBin-1]->GetXaxis()->SetTitle("p_{T} (GeV)");
		graph_sigmaP_L[rapBin-1]->GetYaxis()->SetTitle("RMS * pT / M (mm)");
		graph_sigmaP_L[rapBin-1]->GetYaxis()->SetRangeUser(0., 0.4);
		graph_sigmaP_L[rapBin-1]->GetXaxis()->SetLimits(Xmin, Xmax);
		graph_sigmaP_L[rapBin-1]->SetMarkerStyle(onia::marker_rapForPTBins[rapBin+1]);
		graph_sigmaP_L[rapBin-1]->SetMarkerColor(onia::colour_rapForPTBins[rapBin+1]);
		graph_sigmaP_L[rapBin-1]->SetLineColor(onia::colour_rapForPTBins[rapBin+1]);

		char title[100];
		sprintf(title, "ctau_{#psi}^{%.1f#sigma} (mm)",nSigma);

		graph_CtauCut[rapBin-1]->SetTitle("");
		graph_CtauCut[rapBin-1]->GetXaxis()->SetTitle("p_{T} (GeV)");
		graph_CtauCut[rapBin-1]->GetYaxis()->SetTitle(title);
		graph_CtauCut[rapBin-1]->GetYaxis()->SetRangeUser(0., 0.1);
		graph_CtauCut[rapBin-1]->GetXaxis()->SetLimits(Xmin, Xmax);
		graph_CtauCut[rapBin-1]->SetMarkerStyle(onia::marker_rapForPTBins[rapBin+1]);
		graph_CtauCut[rapBin-1]->SetMarkerColor(onia::colour_rapForPTBins[rapBin+1]);
		graph_CtauCut[rapBin-1]->SetLineColor(onia::colour_rapForPTBins[rapBin+1]);

		graph_PRprob[rapBin-1]->SetTitle("");
		graph_PRprob[rapBin-1]->GetXaxis()->SetTitle("p_{T} (GeV)");
		graph_PRprob[rapBin-1]->GetYaxis()->SetTitle("Prompt Prob.");
		graph_PRprob[rapBin-1]->GetYaxis()->SetRangeUser(0., 1.1);
		graph_PRprob[rapBin-1]->GetXaxis()->SetLimits(Xmin, Xmax);
		graph_PRprob[rapBin-1]->SetMarkerStyle(onia::marker_rapForPTBins[rapBin+1]);
		graph_PRprob[rapBin-1]->SetMarkerColor(onia::colour_rapForPTBins[rapBin+1]);
		graph_PRprob[rapBin-1]->SetLineColor(onia::colour_rapForPTBins[rapBin+1]);
	}

	TFile *outfile  = new TFile(Form("%s/graphFrac_%.1fsigma.root",savePath.str().c_str(),nSigma),"RECREATE");
	for(int rapBin = 1; rapBin < RapBins+1; rapBin++){
		graph_FracPR[rapBin-1]->Write();
		graph_FracNP[rapBin-1]->Write();
		graph_FracBG[rapBin-1]->Write();
		graph_FracBGNP[rapBin-1]->Write();
		graph_CtauCut[rapBin-1]->Write();
		graph_PRprob[rapBin-1]->Write();

		graph_FracPRErr[rapBin-1]->Write();
		graph_FracNPErr[rapBin-1]->Write();

		graph_FracBGErr[rapBin-1]->Write();

		graph_FracPRRelativeErr[rapBin-1]->Write();
		graph_FracNPRelativeErr[rapBin-1]->Write();
		graph_FracBGRelativeErr[rapBin-1]->Write();

		graph_evtPR[rapBin-1]->Write();
		graph_evtNP[rapBin-1]->Write();
		graph_evtBG[rapBin-1]->Write();
		graph_sigmaP[rapBin-1]->Write();
		graph_sigmaP_L[rapBin-1]->Write();
	}
	outfile->Close();

	return;
}


//========================================
void plotEval(double nSigma, int nState, int type){ 
	cout<<"nSigma: "<<nSigma<<endl;
	cout<<"Psi: "<<nState<<"S"<<endl;
	cout<<"type: "<<type<<endl;

	int RapBins = onia::kNbRapForPTBins,
			PtBins  = onia::kNbPTMaxBins;

	std::stringstream savePath;
	savePath << "Fit/parameter/evaluateCtau";
	if(type==0)
		savePath << "/PR_"<<nSigma<<"sigma";
	if(type==1)
		savePath << "/NP_"<<nSigma<<"sigma";
	gSystem->mkdir(savePath.str().c_str(),kTRUE);

	TFile *infile = new TFile(Form("%s/graphFrac_%.1fsigma.root",savePath.str().c_str(),nSigma),"R");
	if(!infile){cout<<"no input file"<<endl; return;}

	TGraphAsymmErrors *graph_FracPR[RapBins], *graph_FracBG[RapBins],
										*graph_FracNP[RapBins], *graph_FracBGNP[RapBins],
										*graph_FracPRErr[RapBins], *graph_FracBGErr[RapBins],
										*graph_FracNPErr[RapBins],
										*graph_FracPRRelativeErr[RapBins], *graph_FracBGRelativeErr[RapBins],
										*graph_FracNPRelativeErr[RapBins],
										*graph_evtPR[RapBins], *graph_evtNP[RapBins], *graph_evtBG[RapBins],
										*graph_sigmaP[RapBins], *graph_sigmaP_L[RapBins],
										*graph_CtauCut[RapBins], *graph_PRprob[RapBins];
	for(int rapBin = 1; rapBin < RapBins+1; rapBin++){
		graph_FracPR[rapBin-1] = (TGraphAsymmErrors*)infile->Get(Form("graph_FracPR_%d",rapBin-1));
		graph_FracNP[rapBin-1] = (TGraphAsymmErrors*)infile->Get(Form("graph_FracNP_%d",rapBin-1));
		graph_FracBG[rapBin-1] = (TGraphAsymmErrors*)infile->Get(Form("graph_FracBG_%d",rapBin-1));
		graph_FracBGNP[rapBin-1] = (TGraphAsymmErrors*)infile->Get(Form("graph_FracBGNP_%d",rapBin-1));
		graph_CtauCut[rapBin-1] = (TGraphAsymmErrors*)infile->Get(Form("graph_CtauCut_%d",rapBin-1));
		graph_PRprob[rapBin-1] = (TGraphAsymmErrors*)infile->Get(Form("graph_PRprob_%d",rapBin-1));

		graph_FracPRErr[rapBin-1] = (TGraphAsymmErrors*)infile->Get(Form("graph_FracPRErr_%d",rapBin-1));
		graph_FracNPErr[rapBin-1] = (TGraphAsymmErrors*)infile->Get(Form("graph_FracNPErr_%d",rapBin-1));
		graph_FracBGErr[rapBin-1] = (TGraphAsymmErrors*)infile->Get(Form("graph_FracBGErr_%d",rapBin-1));

		graph_FracPRRelativeErr[rapBin-1] = (TGraphAsymmErrors*)infile->Get(Form("graph_FracPRRelativeErr_%d",rapBin-1));
		graph_FracNPRelativeErr[rapBin-1] = (TGraphAsymmErrors*)infile->Get(Form("graph_FracNPRelativeErr_%d",rapBin-1));
		graph_FracBGRelativeErr[rapBin-1] = (TGraphAsymmErrors*)infile->Get(Form("graph_FracBGRelativeErr_%d",rapBin-1));

		graph_evtPR[rapBin-1] = (TGraphAsymmErrors*)infile->Get(Form("graph_evtPR_%d",rapBin-1));
		graph_evtNP[rapBin-1] = (TGraphAsymmErrors*)infile->Get(Form("graph_evtNP_%d",rapBin-1));
		graph_evtBG[rapBin-1] = (TGraphAsymmErrors*)infile->Get(Form("graph_evtBG_%d",rapBin-1));

		graph_sigmaP[rapBin-1] = (TGraphAsymmErrors*)infile->Get(Form("graph_sigmaP_%d",rapBin-1));
		graph_sigmaP[rapBin-1] -> GetXaxis()->SetLimits(0.,100.);

		graph_sigmaP_L[rapBin-1] = (TGraphAsymmErrors*)infile->Get(Form("graph_sigmaP_L_%d",rapBin-1));
		graph_sigmaP_L[rapBin-1] -> GetXaxis()->SetLimits(0.,100.);

		if(nState==5){ // for Psi2S, remove first pT bin
			//graph_FracPR[rapBin-1] -> RemovePoint(0);
			//graph_FracNP[rapBin-1] -> RemovePoint(0);
			//graph_FracBG[rapBin-1] -> RemovePoint(0);
			//graph_FracPRErr[rapBin-1] -> RemovePoint(0);
			//graph_FracNPErr[rapBin-1] -> RemovePoint(0);
			//graph_FracBGErr[rapBin-1] -> RemovePoint(0);
			//graph_FracPRRelativeErr[rapBin-1] -> RemovePoint(0);
			//graph_FracNPRelativeErr[rapBin-1] -> RemovePoint(0);
			//graph_FracBGRelativeErr[rapBin-1] -> RemovePoint(0);

			//graph_CtauCut[rapBin-1] -> RemovePoint(0);
			//graph_PRprob[rapBin-1] -> RemovePoint(0);
			//graph_sigmaP[rapBin-1] -> RemovePoint(0);
			//graph_sigmaP_L[rapBin-1] -> RemovePoint(0);

			//graph_evtPR[rapBin-1] -> RemovePoint(0);  
			//graph_evtNP[rapBin-1] -> RemovePoint(0);
			//graph_evtBG[rapBin-1] -> RemovePoint(0);
		}

		//graph_FracNP[rapBin-1] -> SetMarkerColor(kRed+1);

		graph_FracPRRelativeErr[rapBin-1] ->GetYaxis() -> SetRangeUser(0.,50.);
		graph_FracNPRelativeErr[rapBin-1] ->GetYaxis() -> SetRangeUser(0.,50.);
		graph_FracBGRelativeErr[rapBin-1] ->GetYaxis() -> SetRangeUser(0.,50.);
	}

	gStyle->SetPadBottomMargin(0.11); //0.12
	gStyle->SetPadLeftMargin(0.08); //0.12
	gStyle->SetPadRightMargin(0.02); //0.05
	gStyle->SetPadTopMargin(0.05); //0.05

	TCanvas *c1=new TCanvas("c1","");
	c1->SetTickx();
	c1->SetTicky();
	//c1->SetGridx();
	//c1->SetGridy();
	gStyle->SetTitleFillColor(10);
	gStyle->SetTitleBorderSize(1);
	gStyle->SetOptFit(1);           
	gStyle->SetOptStat(1);          
	gStyle->SetTitleFont(22);       
	gStyle->SetStatFont(22);        
	gStyle->SetStatColor(10);       
	gStyle->SetStatBorderSize(1);
	gStyle->SetLabelFont(22,"X");   
	gStyle->SetLabelFont(22,"Y");   
	gStyle->SetTitleXOffset(1.2);   
	gStyle->SetTitleYOffset(1.0);   
	gStyle->SetHistLineWidth(2);    
	gStyle->SetStatX(0.9);          
	gStyle->SetStatY(0.9);          
	gStyle->SetTitleX(0.15);        
	gStyle->SetTitleY(0.96);        


	//double left=0.65, top=0.65, textSize=0.032;
	//double left=0.17, top=0.83, textSize=0.035;
	double leftVal=0.17, topVal=0.88;
	double left=leftVal, top=topVal, textSize=0.035;
	TLatex *latex=new TLatex();
	latex->SetTextFont(42);
	latex->SetNDC(kTRUE);
	latex->SetTextSize(textSize);
	double stepLatex=textSize*1.3;

	//double blX = 0.16, blY = 0.70, trX = 0.4, trY = 0.84;
	//double blX = 0.65, blY = 0.15, trX = 0.85, trY = 0.35;
	//double blX = 0.7, blY = 0.73, trX = 0.83, trY = 0.84;
	//double blX = 0.72, blY = 0.75, trX = 0.93, trY = 0.93;
	//double blX = 0.72, blY = 0.75, trX = 0.93, trY = 0.93;
	double blX = 0.72, blY = 0.77, trX = 0.93, trY = 0.93;
	TLegend* legend_rap1=new TLegend(blX,blY,trX,trY);
	legend_rap1->SetFillColor(kWhite);
	legend_rap1->SetTextFont(42);
	legend_rap1->SetTextSize(legendsize);
	legend_rap1->SetBorderSize(0.);
	legend_rap1->AddEntry(graph_FracPR[0],"Prompt","p");
	legend_rap1->AddEntry(graph_FracNP[0],"Nonprompt","p");
	legend_rap1->AddEntry(graph_FracBG[0],"Background","p");
	//legend_rap1->AddEntry(graph_FracBGNP[0],"Bg+NP","lp");

	TLegend* legend_rap2=new TLegend(blX,blY,trX,trY);
	legend_rap2->SetFillColor(kWhite);

	legend_rap2->SetTextFont(42);
	legend_rap2->SetTextSize(legendsize);
	legend_rap2->SetBorderSize(0.);
	legend_rap2->AddEntry(graph_FracPR[1],"Prompt","p");
	legend_rap2->AddEntry(graph_FracNP[1],"Nonprompt","p");
	legend_rap2->AddEntry(graph_FracBG[1],"Nonprompt","p");
	//legend_rap2->AddEntry(graph_FracBGNP[1],"Bg+NP","lp");

	TLegend* legend_rap3=new TLegend(blX,blY,trX,trY);
	legend_rap3->SetFillColor(kWhite);
	legend_rap3->SetTextFont(42);
	legend_rap3->SetTextSize(legendsize);
	legend_rap3->SetBorderSize(0.);
	if(nState==5){
		legend_rap3->AddEntry(graph_FracPR[2],"Prompt","p");
		legend_rap3->AddEntry(graph_FracNP[2],"Nonprompt","p");
		legend_rap3->AddEntry(graph_FracBG[2],"Nonprompt","p");
		//legend_rap3->AddEntry(graph_FracBGNP[2],"Bg+NP","lp");
	}

	graph_FracPR[0]->Draw("AP");
	graph_FracNP[0]->Draw("P");
	graph_FracBG[0]->Draw("P");
	//graph_FracBGNP[0]->Draw("P");
	legend_rap1->Draw();
	top=topVal;
	if(nState==4) latex->DrawLatex(left,top, "J/#psi");
	if(nState==5) latex->DrawLatex(left,top, "#psi(2S)");
	top -= stepLatex;
	latex->DrawLatex(left,top, "|y| < 0.6");
	c1->SaveAs(Form("%s/fraction_rap1.pdf",savePath.str().c_str()));

	graph_FracPRErr[0]->Draw("AP");
	graph_FracNPErr[0]->Draw("P");
	graph_FracBGErr[0]->Draw("P");
	legend_rap1->Draw();
	top=topVal;
	if(nState==4) latex->DrawLatex(left,top, "J/#psi");
	if(nState==5) latex->DrawLatex(left,top, "#psi(2S)");
	top -= stepLatex;
	latex->DrawLatex(left,top, "|y| < 0.6");
	c1->SaveAs(Form("%s/fractionErr_rap1.pdf",savePath.str().c_str()));

	graph_FracPRRelativeErr[0]->Draw("AP");
	graph_FracNPRelativeErr[0]->Draw("P");
	graph_FracBGRelativeErr[0]->Draw("P");
	legend_rap1->Draw();
	top=topVal;
	if(nState==4) latex->DrawLatex(left,top, "J/#psi");
	if(nState==5) latex->DrawLatex(left,top, "#psi(2S)");

	top -= stepLatex;
	latex->DrawLatex(left,top, "|y| < 0.6");
	c1->SaveAs(Form("%s/fractionRelativeErr_rap1.pdf",savePath.str().c_str()));

	graph_FracPR[1]->Draw("AP");
	graph_FracNP[1]->Draw("P");
	graph_FracBG[1]->Draw("P");
	//graph_FracBGNP[1]->Draw("P");
	legend_rap2->Draw();
	top=topVal;
	if(nState==4) latex->DrawLatex(left,top, "J/#psi");
	if(nState==5) latex->DrawLatex(left,top, "#psi(2S)");
	top -= stepLatex;
	latex->DrawLatex(left,top, "0.6 < |y| < 1.2");
	c1->SaveAs(Form("%s/fraction_rap2.pdf",savePath.str().c_str()));

	graph_FracPRErr[1]->Draw("AP");
	graph_FracNPErr[1]->Draw("P");
	graph_FracBGErr[1]->Draw("P");
	legend_rap2->Draw();
	top=topVal;
	if(nState==4) latex->DrawLatex(left,top, "J/#psi");
	if(nState==5) latex->DrawLatex(left,top, "#psi(2S)");
	top -= stepLatex;
	latex->DrawLatex(left,top, "0.6 < |y| < 1.2");
	c1->SaveAs(Form("%s/fractionErr_rap2.pdf",savePath.str().c_str()));

	graph_FracPRRelativeErr[1]->Draw("AP");
	graph_FracNPRelativeErr[1]->Draw("P");
	graph_FracBGRelativeErr[1]->Draw("P");
	legend_rap2->Draw();
	top=topVal;
	if(nState==4) latex->DrawLatex(left,top, "J/#psi");
	if(nState==5) latex->DrawLatex(left,top, "#psi(2S)");
	top -= stepLatex;
	latex->DrawLatex(left,top, "0.6 < |y| < 1.2");
	c1->SaveAs(Form("%s/fractionRelativeErr_rap2.pdf",savePath.str().c_str()));

	if(nState==5){
		graph_FracPR[2]->Draw("AP");
		graph_FracNP[2]->Draw("P");
		graph_FracBG[2]->Draw("P");
		legend_rap3->Draw();
		top=topVal;
		latex->DrawLatex(left,top, "#psi(2S)");
		top -= stepLatex;
		latex->DrawLatex(left,top, "1.2 < |y| < 1.5");
		c1->SaveAs(Form("%s/fraction_rap3.pdf",savePath.str().c_str()));

		graph_FracPRErr[2]->Draw("AP");
		graph_FracNPErr[2]->Draw("P");
		graph_FracBGErr[2]->Draw("P");
		legend_rap3->Draw();
		top=topVal;
		latex->DrawLatex(left,top, "#psi(2S)");
		top -= stepLatex;
		latex->DrawLatex(left,top, "1.2 < |y| < 1.5");
		c1->SaveAs(Form("%s/fractionErr_rap3.pdf",savePath.str().c_str()));

		graph_FracPRRelativeErr[2]->Draw("AP");
		graph_FracNPRelativeErr[2]->Draw("P");
		graph_FracBGRelativeErr[2]->Draw("P");
		legend_rap3->Draw();
		top=topVal;
		latex->DrawLatex(left,top, "#psi(2S)");
		top -= stepLatex;
		latex->DrawLatex(left,top, "1.2 < |y| < 1.5");
		c1->SaveAs(Form("%s/fractionRelativeErr_rap3.pdf",savePath.str().c_str()));
	}

	TLegend* legend_ctauCut=new TLegend(blX,blY,trX,trY);
	legend_ctauCut->SetFillColor(kWhite);
	legend_ctauCut->SetTextFont(42);
	legend_ctauCut->SetTextSize(legendsize);
	legend_ctauCut->SetBorderSize(0.);
	legend_ctauCut->AddEntry(graph_CtauCut[0],"|y| < 0.6","lp");
	legend_ctauCut->AddEntry(graph_CtauCut[1],"0.6 < |y| < 1.2","lp");
	if(nState==5)
		legend_ctauCut->AddEntry(graph_CtauCut[2],"1.2 < |y| < 1.5","lp");

	TLegend* legend_promptProb=new TLegend(blX,blY,trX,trY);
	legend_promptProb->SetFillColor(kWhite);
	legend_promptProb->SetTextFont(42);
	legend_promptProb->SetTextSize(legendsize);
	legend_promptProb->SetBorderSize(0.);
	legend_promptProb->AddEntry(graph_PRprob[0],"|y| < 0.6","lp");
	legend_promptProb->AddEntry(graph_PRprob[1],"0.6 < |y| < 1.2","lp");
	if(nState==5)
		legend_promptProb->AddEntry(graph_PRprob[2],"1.2 < |y| < 1.5","lp");

	TLegend* legend_evt=new TLegend(blX,blY,trX,trY);
	legend_evt->SetFillColor(kWhite);
	legend_evt->SetTextFont(42);
	legend_evt->SetTextSize(legendsize);
	legend_evt->SetBorderSize(0.);
	legend_evt->AddEntry(graph_evtPR[0],"|y| < 0.6","lp");
	legend_evt->AddEntry(graph_evtPR[1],"0.6 < |y| < 1.2","lp");
	if(nState==5)
		legend_evt->AddEntry(graph_evtPR[2],"1.2 < |y| < 1.5","lp");

	graph_CtauCut[0]->Draw("AP");
	graph_CtauCut[1]->Draw("P");
	if(nState==5)
		graph_CtauCut[2]->Draw("P");
	legend_ctauCut->Draw();

	//left=0.55;top=0.83;
	leftVal=0.6; topVal=0.88;
	left=leftVal;top=topVal;
	if(nState==4) latex->DrawLatex(left,top, "J/#psi");
	if(nState==5) latex->DrawLatex(left,top, "#psi(2S)");
	c1->SaveAs(Form("%s/CtauCut.pdf",savePath.str().c_str()));

	graph_PRprob[0]->Draw("AP");
	graph_PRprob[1]->Draw("P");
	if(nState==5)
		graph_PRprob[2]->Draw("P");
	legend_promptProb->Draw();
	left=leftVal;top=topVal;
	if(nState==4) latex->DrawLatex(left,top, "J/#psi");
	if(nState==5) latex->DrawLatex(left,top, "#psi(2S)");
	c1->SaveAs(Form("%s/PRprob.pdf",savePath.str().c_str()));


	c1->SetLogy();
	graph_evtPR[0]->Draw("AP");
	graph_evtPR[1]->Draw("P");
	if(nState==5)
		graph_evtPR[2]->Draw("P");
	legend_evt->Draw();
	left=leftVal;top=topVal;
	if(nState==4) latex->DrawLatex(left,top, "J/#psi");
	if(nState==5) latex->DrawLatex(left,top, "#psi(2S)");
	c1->SaveAs(Form("%s/evtP.pdf",savePath.str().c_str()));

	graph_evtNP[0]->Draw("AP");
	graph_evtNP[1]->Draw("P");
	if(nState==5)
		graph_evtNP[2]->Draw("P");
	legend_evt->Draw();
	left=leftVal;top=topVal;
	if(nState==4) latex->DrawLatex(left,top, "J/#psi");
	if(nState==5) latex->DrawLatex(left,top, "#psi(2S)");
	c1->SaveAs(Form("%s/evtNP.pdf",savePath.str().c_str()));

	graph_evtBG[0]->Draw("AP");
	graph_evtBG[1]->Draw("P");
	if(nState==5)
		graph_evtBG[2]->Draw("P");
	legend_evt->Draw();
	left=leftVal;top=topVal;
	if(nState==4) latex->DrawLatex(left,top, "J/#psi");
	if(nState==5) latex->DrawLatex(left,top, "#psi(2S)");
	c1->SaveAs(Form("%s/evtBG.pdf",savePath.str().c_str()));

	TLegend* legend_sigmaP=new TLegend(blX,blY,trX,trY);
	legend_sigmaP->SetFillColor(kWhite);
	legend_sigmaP->SetTextFont(42);
	legend_sigmaP->SetTextSize(legendsize);
	legend_sigmaP->SetBorderSize(0.);
	legend_sigmaP->AddEntry(graph_sigmaP[0],"|y| < 0.6","lp");
	legend_sigmaP->AddEntry(graph_sigmaP[1],"0.6 < |y| < 1.2","lp");
	if(nState==5)
		legend_sigmaP->AddEntry(graph_sigmaP[2],"1.2 < |y| < 1.5","lp");

	c1->SetLogy(0);
	graph_sigmaP[0]->Draw("AP");
	graph_sigmaP[1]->Draw("P");
	if(nState==5)
		graph_sigmaP[2]->Draw("P");
	legend_sigmaP->Draw();
	left=leftVal;top=topVal;
	if(nState==4) latex->DrawLatex(left,top, "J/#psi");
	if(nState==5) latex->DrawLatex(left,top, "#psi(2S)");
	c1->SaveAs(Form("%s/rms.pdf",savePath.str().c_str()));

	graph_sigmaP_L[0]->Draw("AP");
	graph_sigmaP_L[1]->Draw("P");
	if(nState==5)
		graph_sigmaP_L[2]->Draw("P");
	legend_sigmaP->Draw();
	left=leftVal;top=topVal;
	if(nState==4) latex->DrawLatex(left,top, "J/#psi");
	if(nState==5) latex->DrawLatex(left,top, "#psi(2S)");
	c1->SaveAs(Form("%s/rmsL.pdf",savePath.str().c_str()));


	//get average sigmaP in rapidity bins:
	//2S
	if(nState==5){
		int NumP =  graph_sigmaP[0]->GetN();
		double sigmaRap0=0, sigmaRap1=0, sigmaRap2=0;
		cout<<"NumP: "<<NumP<<endl;
		for(int i=0; i<NumP; i++){
			double x0,y0,x1,y1,x2,y2;
			graph_sigmaP[0]->GetPoint(i,x0,y0);
			graph_sigmaP[1]->GetPoint(i,x1,y1);
			graph_sigmaP[2]->GetPoint(i,x2,y2);
			double xNew = (x0+x1+x2)/3.;
			double yNew = (y0+y1+y2)/3.;

			if(i>0){
				cout<<"i: "<<i<<endl;
				sigmaRap0 += y0;
				sigmaRap1 += y1;
				sigmaRap2 += y2;
			}
		}
		sigmaRap0 = sigmaRap0/(NumP-1.);
		sigmaRap1 = sigmaRap1/(NumP-1.);
		sigmaRap2 = sigmaRap2/(NumP-1.);
		cout<<"psi 2S:"<<endl;
		cout<<"sigmaRap0: "<<sigmaRap0<<endl;
		cout<<"sigmaRap1: "<<sigmaRap1<<endl;
		cout<<"sigmaRap2: "<<sigmaRap2<<endl;
	}

	if(nState==4){
		int NumP =  graph_sigmaP[0]->GetN();
		double sigmaRap0=0, sigmaRap1=0;
		cout<<"NumP: "<<NumP<<endl;
		for(int i=0; i<NumP; i++){
			double x0,y0,x1,y1;
			graph_sigmaP[0]->GetPoint(i,x0,y0);
			graph_sigmaP[1]->GetPoint(i,x1,y1);
			double xNew = (x0+x1)/3.;
			double yNew = (y0+y1)/3.;

			if(i<11){
				cout<<"i: "<<i<<endl;
				sigmaRap0 += y0;
				sigmaRap1 += y1;
			}
		}
		sigmaRap0 = sigmaRap0/(NumP-1.);
		sigmaRap1 = sigmaRap1/(NumP-1.);
		cout<<"psi 1S:"<<endl;
		cout<<"sigmaRap0: "<<sigmaRap0<<endl;
		cout<<"sigmaRap1: "<<sigmaRap1<<endl;
	}
}

//===================================
vector<double> calculateInte(RooWorkspace *ws, RooDataSet *dataJpsictErr, double ctCutMin, double ctCutMax){
	RooAbsPdf *PRpdf = (RooAbsPdf*)ws->pdf("TotalPromptLifetime");
	RooAbsPdf *NPpdf = (RooAbsPdf*)ws->pdf("nonPromptSSD");
	RooAbsPdf *BGpdf = (RooAbsPdf*)ws->pdf("backgroundlifetime");

	RooRealVar Jpsict(*ws->var("Jpsict"));
	RooRealVar JpsictErr(*ws->var("JpsictErr"));
	Jpsict.setMin(ctMin);   Jpsict.setMax(ctMax);

	RooDataSet *genDataPR = PRpdf->generate(Jpsict,ProtoData(*dataJpsictErr));
	RooDataSet *genDataNP = NPpdf->generate(Jpsict,ProtoData(*dataJpsictErr));
	RooDataSet *genDataBG = BGpdf->generate(Jpsict,ProtoData(*dataJpsictErr));

	TH2F* histPR2D = (TH2F*)genDataPR->createHistogram("histPR2D",Jpsict,Binning(fineBins),YVar(JpsictErr,Binning(fineBins/10)));
	TH1F* histPR   = (TH1F*)histPR2D->ProjectionX();
	TH2F* histNP2D = (TH2F*)genDataNP->createHistogram("histNP2D",Jpsict,Binning(fineBins),YVar(JpsictErr,Binning(fineBins/10)));
	TH1F* histNP   = (TH1F*)histNP2D->ProjectionX();
	TH2F* histBG2D = (TH2F*)genDataBG->createHistogram("histBG2D",Jpsict,Binning(fineBins),YVar(JpsictErr,Binning(fineBins/10)));
	TH1F* histBG   = (TH1F*)histBG2D->ProjectionX();

	histPR->SetLineColor(kRed);
	histPR->SetMarkerColor(kRed);
	histNP->SetLineColor(kBlue);
	histNP->SetMarkerColor(kBlue);
	histBG->SetLineColor(kBlack);
	histBG->SetMarkerColor(kBlack);

	histPR->Scale(1./histPR->Integral());
	histNP->Scale(1./histNP->Integral());
	histBG->Scale(1./histBG->Integral());

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

	//delete PRpdf;
	//delete NPpdf;
	//delete BGpdf;
	delete histPR;
	delete histNP;
	delete histBG;
	delete histPR2D;
	delete histNP2D;
	delete histBG2D;
	delete genDataPR;
	delete genDataNP;
	delete genDataBG;

	return InteRlts;
}

//===================================
vector<double> getSigma(RooWorkspace *ws, RooDataSet *dataJpsictErr, int rapBin, int ptBin){
	int nbins=200;
	TGaxis::SetMaxDigits(3);
	RooRealVar Jpsict(*ws->var("Jpsict"));
	RooRealVar JpsictErr(*ws->var("JpsictErr"));
	Jpsict.setMin(-.2); Jpsict.setMax(.2);

	//RooDataSet *dataSR=(RooDataSet *)ws->data(Form("data_rap%d_pt%d_SR",rapBin,ptBin));
	//int NumEvt = 0., maxEvt = 100000;
	//if(dataSR->numEntries() < maxEvt) NumEvt = dataSR->numEntries();
	//else NumEvt = maxEvt;
	//RooDataSet *dataJpsictErr = (RooDataSet*)dataSR->reduce(SelectVars(RooArgSet(JpsictErr)),
	//		EventRange(0,NumEvt),Name("dataJpsictErr"));

	RooAddModel *Prompt = (RooAddModel*)ws->pdf("TotalPromptLifetime");
	RooDataSet *data = Prompt->generate(Jpsict,ProtoData(*dataJpsictErr));

	TH2F* hist2D = (TH2F*)data->createHistogram("hist2D",Jpsict,Binning(nbins),YVar(JpsictErr,Binning(nbins)));
	TH1F* hist = (TH1F*)hist2D->ProjectionX();

	hist->Scale(1./hist->Integral());
	hist->SetLineColor(kRed);
	hist->SetMarkerColor(kRed);
	hist->GetXaxis()->SetLimits(-.2,.2);

	double rms = hist->GetRMS();
	double rmsE = hist->GetRMSError();
	//cout<<"rms: "<<rms<<endl;
	//cout<<"rmsE: "<<rmsE<<endl;

	vector<double> sigma;
	sigma.push_back(rms);
	sigma.push_back(rmsE);

	//delete Prompt;
	delete hist2D;
	delete hist;
	delete data;

	return sigma;
}
