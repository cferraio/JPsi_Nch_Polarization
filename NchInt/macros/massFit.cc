#include "rootIncludes.inc"
#include "commonVar.h"

using namespace RooFit;

void massFit(const std::string &infilename, int rapBin, int ptBin, int nState, bool fitMassPR, bool fitMassNP){
	TFile *infile = new TFile(infilename.c_str(), "UPDATE");
	if(!infile){
		std::cout << "Error: failed to open file with dataset" << std::endl;
	}
	RooWorkspace *ws=(RooWorkspace *)infile->Get("ws_masslifetime");
	if(!ws){
		std::cout << "Error: failed to open workspace " << std::endl;
	}

	//------------------------------------------------------------------------------------------------------------------
	// mass pdf
	std::cout << "Building mass pdf" << std::endl;

	//define mass shape
	//signal (prompt+noprompt) 
	ws->factory("RooCBShape::massCBShape(JpsiMass,CBmass[3.1,3.05,3.15],CBsigma[0.02,0.0001,0.1],CBalpha[1,.0001,6],CBn[10,.0001,60])");
	ws->factory("RooCBShape::massCBShape2(JpsiMass,CBmass,CBsigma2[0.02,0.0001,0.1],CBalpha,CBn)");
	ws->factory("SUM::sigMassShape(fracCB1[0.5,0.,1.]*massCBShape, massCBShape2)");
	//backgrond
	ws->factory("RooExponential::bkgMassShape(JpsiMass,bkgLambda[0,-5,5])");
	//full Mass shape
	ws->factory("SUM::massModel(fracBkg[0.1,0.,1.]*bkgMassShape, sigMassShape)");

	ws->Print("v");

	//------------------------------------------------------------------------------------------------------------------
	// do fit
	std::cout << "Fitting" << std::endl;

	RooRealVar JpsiMass(*ws->var("JpsiMass"));

	if(nState==4){
		ws->var("CBmass")->setMax(3.15);
		ws->var("CBmass")->setMin(3.05);
		ws->var("CBmass")->setVal(3.1);

		ws->var("CBalpha")->setVal(2.);
		ws->var("CBn")->setVal(2.5);
		ws->var("CBn")->setConstant(kTRUE);

		ws->var("CBsigma")->setMax(0.07);
		ws->var("CBsigma2")->setMax(0.07);
		ws->var("CBsigma")->setVal(.02);
		ws->var("CBsigma2")->setVal(.02);
		ws->var("fracCB1")->setVal(0.5);
		ws->var("bkgLambda")->setVal(-1.5);

		ws->var("fracBkg")->setVal(0.1);
	}
	else if(nState==5){
		ws->var("CBmass")->setMax(3.75);
		ws->var("CBmass")->setMin(3.65);
		ws->var("CBmass")->setVal(3.69);
		ws->var("CBalpha")->setVal(2.);
		ws->var("CBn")->setVal(2.5);
		ws->var("CBn")->setConstant(kTRUE);

		//if(rapBin == 1 && ptBin == 1){
		//  ws->var("CBsigma")->setMax(0.035);
		//  ws->var("CBsigma2")->setMax(0.035);
		//}

		ws->var("CBsigma")->setMax(0.07);
		ws->var("CBsigma2")->setMax(0.07);
		if(rapBin == 3){
			ws->var("CBsigma")->setMin(0.03);
			ws->var("CBsigma2")->setMin(0.03);
		}
		if(rapBin == 3 && ptBin > 3)
			ws->var("fracCB1")->setMin(0.1);

		ws->var("CBsigma")->setVal(.03);
		ws->var("CBsigma2")->setVal(.04);
		if(rapBin == 1 && ptBin == 1){
			ws->var("fracCB1")->setVal(0.5);
			ws->var("bkgLambda")->setVal(-1.5);
			ws->var("fracBkg")->setVal(0.1);
		}
	}

	RooAddPdf *MPdf = (RooAddPdf*)ws->pdf("massModel");

	std::stringstream binName;
	binName << "data_rap" << rapBin << "_pt" << ptBin;
	RooDataSet *data = (RooDataSet*)ws->data(binName.str().c_str());

	data->Print("v");

	////////////////////////////////////////////////////////////////
	RooRealVar *JpsiPt = ws->var("JpsiPt");
	TH1* histPt = data->createHistogram("histPt", *JpsiPt, Binning(120));
	double meanPt = histPt->GetMean();

	//// define sigma of prompt p.d.f., got from fit the trend
	//// define function y = a + b * pT
	double a = 0.073, b = 0.0027;
	//proper decay length
	double L_decay = a + b * meanPt;
	//pseudo-proper decay length
	double l_pdecay = L_decay * onia::MpsiPDG / meanPt ;
	double nSigma = 0.;
	if(nState==4) nSigma = 3.0;
	if(nState==5) nSigma = 3.0;
	double ctauCut = nSigma*l_pdecay;
	double ctCutMinPR = -ctauCut, ctCutMaxPR = ctauCut;
	double ctCutMinNP =  ctauCut, ctCutMaxNP = 6.;

	std::stringstream cutlife;
	if(fitMassPR){
		cout << " fitting mass in prompt lifetime region " << endl;
		cutlife << "Jpsict > " << ctCutMinPR << " && Jpsict < " << ctCutMaxPR ;
		data = (RooDataSet*)data -> reduce(cutlife.str().c_str());
	}
	else if(fitMassNP){
		cout << " fitting mass in non-prompt lifetime region " << endl;
		cutlife << "Jpsict > " << ctCutMinNP << " && Jpsict < " << ctCutMaxNP ;
		data = (RooDataSet*)data -> reduce(cutlife.str().c_str());
	}
	////////////////////////////////////////////////////////////////

	data->Print("v");

	RooArgSet *NLLs = new RooArgSet();
	RooAbsReal *MassNLL = NULL; 

	MassNLL = (RooAbsReal *)MPdf->createNLL(*data, NumCPU(4));

	NLLs->add(*MassNLL);

	RooAddition *simNLL = new RooAddition("add","add",*NLLs);
	RooMinuit *mMinuit = new RooMinuit(*simNLL);

	mMinuit->setStrategy(2);
	mMinuit->setPrintEvalErrors(-1);
	mMinuit->setEvalErrorWall(false);
	mMinuit->setVerbose(false);

	mMinuit->simplex();
	mMinuit->migrad();
	mMinuit->migrad();
	mMinuit->hesse();

	std::stringstream rltName, snapshotName;
	rltName<< "m_fitresult_rap"<<rapBin<<"_pt"<<ptBin;
	snapshotName<< "m_snapshot_rap"<<rapBin<<"_pt"<<ptBin;
	RooFitResult *fitresult = (RooFitResult*)mMinuit->save(rltName.str().c_str());
	ws->import(*fitresult,kTRUE);
	fitresult->Print();

	ws->saveSnapshot(snapshotName.str().c_str(), ws->allVars());

	ws->Write();
	infile->Close();
}