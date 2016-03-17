#include "rootIncludes.inc"
#include "commonVar.h"
#include "RooHistPdf.h"
#include "RooDataHist.h"

using namespace RooFit;

void massFit(const std::string &infilename, int rapBin, int ptBin, int cpmBin, int nState, bool fitMassPR, bool fitMassNP){
	TFile *infile = new TFile(infilename.c_str(), "UPDATE");
	if(!infile){
		std::cout << "Error: failed to open file with dataset" << std::endl;
	}
	cout<<"ptBin is "<<ptBin<<" and rap bin is "<<rapBin<<endl;
	cout<<infilename<<endl;
	RooWorkspace *ws=(RooWorkspace *)infile->Get("ws_masslifetime");
	if(!ws){
		std::cout << "Error: failed to open workspace " << std::endl;
	}

	//------------------------------------------------------------------------------------------------------------------
	// mass pdf
	std::cout << "Building mass pdf" << std::endl;

	//define mass shape
	//signal (prompt+noprompt) 
	
	RooRealVar CBsigma_p0("CBsigma_p0","CBsigma_p0",0.02,0.015,0.025);
	RooRealVar CBsigma_p1("CBsigma_p1","CBsigma_p1",0.,-0.01,0.01);
	RooRealVar CBsigma_p2("CBsigma_p2","CBsigma_p2",0.01,0.005,0.025);
	ws->import(RooArgList(CBsigma_p0,CBsigma_p1,CBsigma_p2));
	RooFormulaVar CBsigma("CBsigma","@0+@1*abs(@3)+@2*abs(@3)*abs(@3)",RooArgList(*ws->var("CBsigma_p0"), *ws->var("CBsigma_p1"), *ws->var("CBsigma_p2"), *ws->var("JpsiRap")));
	ws->import(CBsigma);
	
	RooRealVar CBalpha_p0("CBalpha_p0","CBalpha_p0",1.729,1.2,2.5);
	RooRealVar CBalpha_p1("CBalpha_p1","CBalpha_p1",0.191,0.,0.5);
	ws->import(RooArgList(CBalpha_p0,CBalpha_p1));
	RooFormulaVar CBalpha("CBalpha","@0+@1*abs(@2)",RooArgList(*ws->var("CBalpha_p0"), *ws->var("CBalpha_p1"), *ws->var("JpsiRap")));
	ws->import(CBalpha);
	
	RooRealVar CBmass_p0("CBmass_p0","CBmass_p0",3.094,3.086,3.098);
	RooRealVar CBmass_p1("CBmass_p1","CBmass_p1",0.001,-0.002,0.002);
	RooRealVar CBmass_p2("CBmass_p2","CBmass_p2",-0.003,-0.005,0.001);
	ws->import(RooArgList(CBmass_p0,CBmass_p1,CBmass_p2));
	RooFormulaVar CBmass("CBmass","@0+@1*abs(@3)+@2*abs(@3)*abs(@3)",RooArgList(*ws->var("CBmass_p0"), *ws->var("CBmass_p1"), *ws->var("CBmass_p2"), *ws->var("JpsiRap")));
	ws->import(CBmass);
	
	ws->factory("RooCBShape::sigMassShape(JpsiMass,CBmass,CBsigma,CBalpha,CBn[2.5])");
	ws->factory("Gaussian::gaussMassShape(JpsiMass,CBmass,CBsigma)");

    RooRealVar fracCB1("fracCB1","fracCB1",0.5,0.,1.);
    ws->import(fracCB1);
    RooRealVar CBsigma2("CBsigma2","CBsigma2",0.5,0.,1.);
    ws->import(CBsigma2);

	//backgrond
	ws->factory("RooExponential::bkgMassShape(JpsiMass,bkgLambda[0,-5,5])");
	
	//full Mass shape
	ws->factory("SUM::massModel(fracBkg[0.1,0.,1.]*bkgMassShape, sigMassShape)");

	
	ws->Print("v");

	//------------------------------------------------------------------------------------------------------------------
	// do fit
	std::cout << "Fitting" << std::endl;

	

	if(nState==4){
//		ws->var("CBmass")->setMax(3.15);
//		ws->var("CBmass")->setMin(3.05);
//		ws->var("CBmass")->setVal(3.1);

//		ws->var("CBalpha")->setVal(2.);
		ws->var("CBn")->setVal(2.5);
		ws->var("CBn")->setConstant(kTRUE);



        ws->var("CBsigma_p1")->setVal(0.);
        ws->var("CBsigma_p1")->setConstant(kTRUE);

		ws->var("CBmass_p1")->setVal(0.);
		ws->var("CBmass_p1")->setConstant(kTRUE);
		ws->var("CBmass_p2")->setVal(0.);
		ws->var("CBmass_p2")->setConstant(kTRUE);

		ws->var("CBalpha_p1")->setVal(0.);
		ws->var("CBalpha_p1")->setConstant(kTRUE);

        ws->var("CBsigma_p2")->setVal(0.0125);
		ws->var("CBsigma_p2")->setConstant(kTRUE);
/////////////		

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
	binName << "data_rap" << rapBin << "_pt" << ptBin << "_cpm" << cpmBin;;
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

	MassNLL = (RooAbsReal *)ws->pdf("massModel")->createNLL(*data, ConditionalObservables(*ws->var("JpsiRap")), NumCPU(4));

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
	rltName<< "m_fitresult_rap"<<rapBin<<"_pt"<<ptBin<<"_cpm"<<cpmBin;
	snapshotName<< "m_snapshot_rap"<<rapBin<<"_pt"<<ptBin<<"_cpm"<<cpmBin;
	RooFitResult *fitresult = (RooFitResult*)mMinuit->save(rltName.str().c_str());
	ws->import(*fitresult,kTRUE);
	fitresult->Print();
	
	RooRealVar *JpsiMass = ws->var("JpsiMass");
	assert( 0 != JpsiMass );
	RooRealVar *JpsiRap = ws->var("JpsiRap");
	assert( 0 != JpsiRap );
	RooAbsPdf *sigMassShape = ws->pdf("sigMassShape");
	assert ( 0 != sigMassShape );
	
	RooDataSet *dataJpsiRap = (RooDataSet*)data->reduce(SelectVars(RooArgSet(*JpsiRap)),Name("dataJpsiRap"));

	int nevt=1000000;
	cout<<"generating datasets"<<endl;
	RooDataSet *SignalPseudoData = sigMassShape->generate(*JpsiMass,ProtoData(*dataJpsiRap),NumEvents(nevt));
	cout<<"finished generating datasets"<<endl;
	
	int nbinsHists=100;
	TH1F* sigMassShape_asHist = new TH1F("sigMassShape_asHist","sigMassShape_asHist", nbinsHists, onia::massMin, onia::massMax);
	SignalPseudoData->fillHistogram(sigMassShape_asHist,RooArgList(*JpsiMass));
	RooDataHist* signalMassShape_asRooDataHist = new RooDataHist("signalMassShape_asRooDataHist","signalMassShape_asRooDataHist", RooArgList(*JpsiMass), sigMassShape_asHist);
	RooHistPdf* signalMassShape_asHistPdf = new RooHistPdf("signalMassShape_asHistPdf","signalMassShape_asHistPdf", RooArgSet(*JpsiMass), *signalMassShape_asRooDataHist, 3);
	
	ws->import(*signalMassShape_asHistPdf);
	
	ws->factory("SUM::FullmassPdf(fracBkg*bkgMassShape, signalMassShape_asHistPdf)");

	ws->saveSnapshot(snapshotName.str().c_str(), ws->allVars());

	ws->Write();
	infile->Close();
}