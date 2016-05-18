#include "rootIncludes.inc"
#include "commonVar.h"
#include "RooUtils.h"
#include "RooAbsPdf.h"
#include "RooDataSet.h"
#include "RooRealVar.h"
#include "RooAbsReal.h"
#include "RooDataHist.h"
#include "RooHistPdf.h"
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

//	RooRealVar JpsiMass(*ws->var("JpsiMass"));

	if(nState==4){
//		ws->var("CBmass")->setMax(3.15);
//		ws->var("CBmass")->setMin(3.05);
//		ws->var("CBmass")->setVal(3.1);

//		ws->var("CBalpha")->setVal(2.);

		ws->var("CBn")->setVal(2.5);
		ws->var("CBn")->setConstant(kTRUE);


		ws->var("CBmass_p1")->setVal(0.);
		ws->var("CBmass_p1")->setConstant(kTRUE);

		ws->var("CBmass_p2")->setVal(0.);
		ws->var("CBmass_p2")->setConstant(kTRUE);

		ws->var("CBalpha_p1")->setVal(0.);
		ws->var("CBalpha_p1")->setConstant(kTRUE);

        ws->var("CBsigma_p2")->setVal(0.0125);
		ws->var("CBsigma_p2")->setConstant(kTRUE);
		
		ws->var("CBsigma_p1")->setVal(0.);
		ws->var("CBsigma_p1")->setConstant(kTRUE);

		

		ws->var("fracCB1")->setVal(0.5);
		ws->var("bkgLambda")->setVal(-1.5);

		ws->var("fracBkg")->setVal(0.1);
	}
/*	else if(nState==5){
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
	} */

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

	ws->saveSnapshot(snapshotName.str().c_str(), ws->allVars());
	
	////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////
	
	RooRealVar *JpsiMass = ws->var("JpsiMass");
	assert( 0 != JpsiMass );
	
	RooRealVar *JpsiRap = ws->var("JpsiRap");
	assert( 0 != JpsiRap );

	RooAbsPdf *massPdf = ws->pdf("massModel");
	assert ( 0 != massPdf );
	RooAbsPdf *bkgMassShape = ws->pdf("bkgMassShape");
	assert ( 0 != bkgMassShape );
	RooAbsPdf *sigMassShape = ws->pdf("sigMassShape");
	assert ( 0 != sigMassShape );
	RooAbsPdf *gaussMassShape = ws->pdf("gaussMassShape");
	assert ( 0 != gaussMassShape );

	int nEntries = data->numEntries();
	
	RooDataSet *dataJpsiRap = (RooDataSet*)data->reduce(SelectVars(RooArgSet(*JpsiRap)),Name("dataJpsiRap"));
	
	int nevt=1000000;
	cout<<"generating datasets"<<endl;
	RooDataSet *GaussPseudoData = gaussMassShape->generate(*JpsiMass,ProtoData(*dataJpsiRap),NumEvents(nevt/100.));
	RooDataSet *BackgroundPseudoData = bkgMassShape->generate(*JpsiMass,ProtoData(*dataJpsiRap),NumEvents(nevt/100.));
	RooDataSet *SignalPseudoData = sigMassShape->generate(*JpsiMass,ProtoData(*dataJpsiRap),NumEvents(nevt));
	cout<<"finished generating datasets"<<endl;
	
	int nbinsHists=100;
	TH1F* sigMassShape_asHist = new TH1F("sigMassShape_asHist","sigMassShape_asHist", nbinsHists, onia::massMin, onia::massMax);
	SignalPseudoData->fillHistogram(sigMassShape_asHist,RooArgList(*JpsiMass));
	RooDataHist* signalMassShape_asRooDataHist = new RooDataHist("signalMassShape_asRooDataHist","signalMassShape_asRooDataHist", RooArgList(*JpsiMass), sigMassShape_asHist);
	RooHistPdf* signalMassShape_asHistPdf = new RooHistPdf("signalMassShape_asHistPdf","signalMassShape_asHistPdf", RooArgSet(*JpsiMass), *signalMassShape_asRooDataHist, 3);


	
	ws->import(*signalMassShape_asHistPdf);

	ws->factory("SUM::FullmassPdf(fracBkg*bkgMassShape, signalMassShape_asHistPdf)");


    std::stringstream cutSR;
    std::stringstream cutLSB;
    std::stringstream cutRSB;
    cutSR << "TMath::Abs(JpsiMass-" << ws->var("CBmass_p0")->getVal() << ") < " << onia::nSigMass << "*(" << ws->var("CBsigma_p0")->getVal() << "+"<< ws->var("CBsigma_p1")->getVal() << "*TMath::Abs(JpsiRap) + "<< ws->var("CBsigma_p2")->getVal() << "*TMath::Abs(JpsiRap)*TMath::Abs(JpsiRap))";
    cutLSB << "JpsiMass < " << ws->var("CBmass_p0")->getVal() << " - " << onia::nSigBkgLow << "*(" << ws->var("CBsigma_p0")->getVal() << "+"<< ws->var("CBsigma_p1")->getVal() << "*TMath::Abs(JpsiRap) + "<< ws->var("CBsigma_p2")->getVal() << "*TMath::Abs(JpsiRap)*TMath::Abs(JpsiRap))";
    cutRSB << "JpsiMass > " << ws->var("CBmass_p0")->getVal() << " + " << onia::nSigBkgHigh << "*(" << ws->var("CBsigma_p0")->getVal() << "+"<< ws->var("CBsigma_p1")->getVal() << "*TMath::Abs(JpsiRap) + "<< ws->var("CBsigma_p2")->getVal() << "*TMath::Abs(JpsiRap)*TMath::Abs(JpsiRap))";

    cout<<cutSR.str().c_str()<<endl;
    cout<<cutLSB.str().c_str()<<endl;
    cout<<cutRSB.str().c_str()<<endl;

	std::stringstream binNameSR, binNameLSB, binNameRSB;
	binNameSR  << "data_rap" << rapBin << "_pt" << ptBin << "_cpm" << cpmBin << "_SR";
	binNameLSB << "data_rap" << rapBin << "_pt" << ptBin << "_cpm" << cpmBin << "_SBL";
	binNameRSB << "data_rap" << rapBin << "_pt" << ptBin << "_cpm" << cpmBin << "_SBR";

	RooAbsData* dataSR, *dataLSB, *dataRSB;
	int events=0;

	dataSR	= data->reduce(Cut(cutSR.str().c_str()));
	dataLSB = data->reduce(Cut(cutLSB.str().c_str()));
	dataRSB = data->reduce(Cut(cutRSB.str().c_str()));

	dataSR->SetNameTitle(binNameSR.str().c_str(), "data in signal region");
	dataLSB->SetNameTitle(binNameLSB.str().c_str(), "data in LSB");
	dataRSB->SetNameTitle(binNameRSB.str().c_str(), "data in RSB");

//	ws->import(*dataSR);
//	ws->import(*dataLSB);
//	ws->import(*dataRSB);


    RooAbsData* SignalPseudoDataSR = SignalPseudoData->reduce(Cut(cutSR.str().c_str()));
    RooAbsData* SignalPseudoDataLSB = SignalPseudoData->reduce(Cut(cutLSB.str().c_str()));
    RooAbsData* SignalPseudoDataRSB = SignalPseudoData->reduce(Cut(cutRSB.str().c_str()));

    SignalPseudoDataSR->Print();

	double fracSigEventsInSR= double(SignalPseudoDataSR->numEntries()) / double(SignalPseudoData->numEntries());
	cout<<"fracSigEventsInSR = "<<fracSigEventsInSR<<endl;
	double fracSigEventsInLSB= double(SignalPseudoDataLSB->numEntries()) / double(SignalPseudoData->numEntries());
	cout<<"fracSigEventsInLSB = "<<fracSigEventsInLSB<<endl;
	double fracSigEventsInRSB= double(SignalPseudoDataRSB->numEntries()) / double(SignalPseudoData->numEntries());
	cout<<"fracSigEventsInRSB = "<<fracSigEventsInRSB<<endl;

    RooAbsData* BackgroundPseudoDataSR = BackgroundPseudoData->reduce(Cut(cutSR.str().c_str()));
    RooAbsData* BackgroundPseudoDataLSB = BackgroundPseudoData->reduce(Cut(cutLSB.str().c_str()));
    RooAbsData* BackgroundPseudoDataRSB = BackgroundPseudoData->reduce(Cut(cutRSB.str().c_str()));

	double fracBGEventsInSR= double(BackgroundPseudoDataSR->numEntries()) / double(BackgroundPseudoData->numEntries());
	cout<<"fracBGEventsInSR = "<<fracBGEventsInSR<<endl;
	double fracBGEventsInLSB= double(BackgroundPseudoDataLSB->numEntries()) / double(BackgroundPseudoData->numEntries());
	cout<<"fracBGEventsInLSB = "<<fracBGEventsInLSB<<endl;
	double fracBGEventsInRSB= double(BackgroundPseudoDataRSB->numEntries()) / double(BackgroundPseudoData->numEntries());
	cout<<"fracBGEventsInRSB = "<<fracBGEventsInRSB<<endl;


	int nbinsSigmaDef=200;
	TH1F* hist1D = (TH1F*)GaussPseudoData->createHistogram("hist1D",*JpsiMass,Binning(nbinsSigmaDef));
	hist1D->Scale(1./hist1D->Integral());
	hist1D->SetLineColor(kRed);
	hist1D->SetMarkerColor(kRed);
	hist1D->GetXaxis()->SetLimits(-.2,.2);

	double massres = hist1D->GetRMS();
	double err_massres = hist1D->GetRMSError();



	double fracBG = ws->var("fracBkg")->getVal();
	double fracBGErr = ws->var("fracBkg")->getError();

	double relativeErr_fracBG=fracBGErr/fracBG;
	double relativeErr_fracSig=fracBGErr/(1.-fracBG);

	double n_SigInSR = nEntries*(1.-fracBG)*fracSigEventsInSR;
	double n_BGInSR = nEntries*(fracBG)*fracBGEventsInSR;

	double n_SigInLSB = nEntries*(1.-fracBG)*fracSigEventsInLSB;
	double n_BGInLSB = nEntries*(fracBG)*fracBGEventsInLSB;

	double n_SigInRSB = nEntries*(1.-fracBG)*fracSigEventsInRSB;
	double n_BGInRSB = nEntries*(fracBG)*fracBGEventsInRSB;

	double frac_SigInSR = n_SigInSR / (n_SigInSR+n_BGInSR);
	double frac_BGInSR = n_BGInSR / (n_SigInSR+n_BGInSR);

	double frac_SigInLSB = n_SigInLSB / (n_SigInLSB+n_BGInLSB);
	double frac_BGInLSB = n_BGInLSB / (n_SigInLSB+n_BGInLSB);

	double frac_SigInRSB = n_SigInRSB / (n_SigInRSB+n_BGInRSB);
	double frac_BGInRSB = n_BGInRSB / (n_SigInRSB+n_BGInRSB);

	cout<<"frac_SigInSR = "<<frac_SigInSR<<endl;
	cout<<"frac_BGInSR = "<<frac_BGInSR<<endl;

	cout<<"frac_SigInLSB = "<<frac_SigInLSB<<endl;
	cout<<"frac_BGInLSB = "<<frac_BGInLSB<<endl;

	cout<<"frac_SigInRSB = "<<frac_SigInRSB<<endl;
	cout<<"frac_BGInRSB = "<<frac_BGInRSB<<endl;
	
	 RooRealVar var_frac_SigInSR("var_frac_SigInSR","var_frac_SigInSR",frac_SigInSR);  
	 var_frac_SigInSR.setError(frac_SigInSR*relativeErr_fracSig); 
	 if(!ws->var("var_frac_SigInSR")) ws->import(var_frac_SigInSR); 
	 else {
	 	ws->var("var_frac_SigInSR")->setVal(frac_SigInSR); 
	 	ws->var("var_frac_SigInSR")->setError(frac_SigInSR*relativeErr_fracSig);
	 }
    RooRealVar FracBkg("FracBkg","FracBkg",frac_BGInSR);  FracBkg.setError(frac_BGInSR*relativeErr_fracBG); if(!ws->var("FracBkg")) ws->import(FracBkg); else {ws->var("FracBkg")->setVal(frac_BGInSR); ws->var("FracBkg")->setError(frac_BGInSR*relativeErr_fracBG);}
    RooRealVar FracBkgErr("FracBkgErr","FracBkgErr",ws->var("FracBkg")->getError());
    ws->import(FracBkgErr);
    
    RooRealVar var_frac_SigInLSB("var_frac_SigInLSB","var_frac_SigInLSB",frac_SigInLSB);  var_frac_SigInLSB.setError(frac_SigInLSB*relativeErr_fracSig); if(!ws->var("var_frac_SigInLSB")) ws->import(var_frac_SigInLSB); else {ws->var("var_frac_SigInLSB")->setVal(frac_SigInLSB); ws->var("var_frac_SigInLSB")->setError(frac_SigInLSB*relativeErr_fracSig);}
    

    RooRealVar FracBkgSBL("FracBkgSBL","FracBkgSBL",frac_BGInLSB);
    FracBkgSBL.setError(frac_BGInLSB*relativeErr_fracBG);
    if(!ws->var("FracBkgSBL")) ws->import(FracBkgSBL);
    else {
    	ws->var("FracBkgSBL")->setVal(frac_BGInLSB);
	    ws->var("FracBkgSBL")->setError(frac_BGInLSB*relativeErr_fracBG);
    }
    RooRealVar FracBkgSBLErr("FracBkgSBLErr","FracBkgSBLErr",ws->var("FracBkgSBL")->getError());
    ws->import(FracBkgSBLErr);
    
    RooRealVar FracBkgSBR("FracBkgSBR","FracBkgSBR",frac_BGInRSB);  FracBkgSBR.setError(frac_BGInRSB*relativeErr_fracBG);
    if(!ws->var("FracBkgSBR")) ws->import(FracBkgSBR);
    else {ws->var("FracBkgSBR")->setVal(frac_BGInRSB); ws->var("FracBkgSBR")->setError(frac_BGInRSB*relativeErr_fracBG);
    }
    RooRealVar FracBkgSBRErr("FracBkgSBRErr","FracBkgSBRErr",ws->var("FracBkgSBR")->getError());
    ws->import(FracBkgSBRErr);
    
    RooRealVar var_frac_SigInRSB("var_frac_SigInRSB","var_frac_SigInRSB",frac_SigInRSB);  var_frac_SigInRSB.setError(frac_SigInRSB*relativeErr_fracSig); if(!ws->var("var_frac_SigInRSB")) ws->import(var_frac_SigInRSB); else {ws->var("var_frac_SigInRSB")->setVal(frac_SigInRSB); ws->var("var_frac_SigInRSB")->setError(frac_SigInRSB*relativeErr_fracSig);}
    
    


	////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////


	ws->Write();
	infile->Close();
}