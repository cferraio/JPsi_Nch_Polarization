#include "calculatePar.cc"

void buildLifetimePDF(RooWorkspace *ws);
void doFit(RooWorkspace *ws, int nState, double BkgRatio3Sig, double fracBkgInSBL, double fracBkgInSBR, int rapBin, int ptBin, int cpmBin);

//=================================================================================
void lifetimeFit(const std::string &infilename, int rapBin, int ptBin, int cpmBin, int nState){

	TFile* infile = new TFile(infilename.c_str(), "UPDATE");
	if(infile->IsZombie()){
		std::cout << "Error: failed to open mass root file" << std::endl;
		return;
	}
	std::string workspacename = "ws_masslifetime";
	RooWorkspace *ws=dynamic_cast<RooWorkspace*>(infile->Get(workspacename.c_str()));
	if(ws == 0){
		std::cout << "Error: failed to get workspace " << workspacename << " from file " << infilename << std::endl;
		infile->ls();
		return;
	}

	std::stringstream dataname;
	dataname << "data_rap" << rapBin << "_pt" << ptBin << "_cpm" << cpmBin;
	RooAbsData* data = ws->data(dataname.str().c_str());
	cout<<"dataname is "<<dataname<<endl;

	// calculate weighted sigma from mass fit
	RooRealVar *JpsiMass = ws->var("JpsiMass");
	RooRealVar *CBmass =  ws->var("CBmass");
	RooRealVar *CBsigma = ws->var("CBsigma");
	RooRealVar *CBsigma2 = ws->var("CBsigma2");
	RooRealVar *fracCB1_ = ws->var("fracCB1");
	RooRealVar *fracBkg_=(RooRealVar *)ws->var("fracBkg");
	double Mean = CBmass->getVal();
	double Sigma = CBsigma->getVal();
	double Sigma2 = CBsigma2->getVal();
	double fracCB1 = fracCB1_->getVal();
	double SigmaWei = sqrt(pow(Sigma,2)*fracCB1+pow(Sigma2,2)*(1-fracCB1));
	double fracBkg = fracBkg_->getVal();
	double fracBkgErr = fracBkg_->getError();

	// calculate signal mass range
	double sigMaxMass = Mean+SigmaWei*onia::nSigMass;
	double sigMinMass = Mean-SigmaWei*onia::nSigMass;
	double sbHighMass = Mean+SigmaWei*onia::nSigBkgHigh;
	double sbLowMass =  Mean-SigmaWei*onia::nSigBkgLow;

	std::cout << "-------------- mass range --------------\n"
		<< "signal range: " << sigMinMass << " - " << sigMaxMass << "\n"
		<< "sbHighMass: " << sbHighMass << "\n"
		<< "sbLowMass: " << sbLowMass << "\n"
		<< "----------------------------------------" << std::endl;
	// create datasets for LSB, RSB and SR
	std::stringstream cutSR, cutSBL, cutSBR;
	cutSR << "JpsiMass > " << sigMinMass << " && JpsiMass < " << sigMaxMass; 
	cutSBL << "JpsiMass > " << onia::massMin << " && JpsiMass < " << sbLowMass; 
	cutSBR << "JpsiMass > " << sbHighMass << " && JpsiMass < " << onia::massMax; 
	std::stringstream binNameSR, binNameSBL, binNameSBR;
	binNameSR  << "data_rap" << rapBin << "_pt" << ptBin << "_cpm" << cpmBin << "_SR";
	binNameSBL << "data_rap" << rapBin << "_pt" << ptBin << "_cpm" << cpmBin << "_SBL";
	binNameSBR << "data_rap" << rapBin << "_pt" << ptBin << "_cpm" << cpmBin << "_SBR";
	RooAbsData* dataSR, *dataSBL, *dataSBR; 
	int events=0;
	if(nState==4 && (ptBin==1 && cpmBin>2 && cpmBin<7)){
		// for 1S, rap1,pt 1 2 3 4; rap2,pt 1 2 3 4, not use the full statistics
		//1S, [1,4], SR:231029, SBL:16500, SBR:5521
		//1S, [2,5]: SR:144617, SBL:7672,  SBR:2483
		events = data->numEntries()/(6);
		dataSR	= data->reduce(Cut(cutSR.str().c_str()),
				EventRange(0,events));
		dataSBL = data->reduce(Cut(cutSBL.str().c_str()),
				EventRange(0,events));
		dataSBR = data->reduce(Cut(cutSBR.str().c_str()),
				EventRange(0,events));
	}	
	//else if(nState==5 && rapBin==2 && ptBin==1){
	//	events = data->numEntries()/5;
	//	dataSR  = data->reduce(Cut(cutSR.str().c_str()),
	//			    EventRange(0,events));
	//	dataSBL = data->reduce(Cut(cutSBL.str().c_str()),
	//			    EventRange(0,events));
	//	dataSBR = data->reduce(Cut(cutSBR.str().c_str()),
	//			    EventRange(0,events));
	//}
	else{
		dataSR	= data->reduce(Cut(cutSR.str().c_str()));
		dataSBL = data->reduce(Cut(cutSBL.str().c_str()));
		dataSBR = data->reduce(Cut(cutSBR.str().c_str()));
	}
	//***
	if(nState==4 && ptBin==2){
//		cout<<"**** Using less statisticals for lifetime fit"<<endl;
		events = data->numEntries();
//		if(ptBin==2) events = data->numEntries()/(11-cpmBin);
		dataSR  = data->reduce(Cut(cutSR.str().c_str()),
				EventRange(0,events));
		dataSBL = data->reduce(Cut(cutSBL.str().c_str()),
				EventRange(0,events));
		dataSBR = data->reduce(Cut(cutSBR.str().c_str()),
				EventRange(0,events));
	} 
	//***
	dataSR->SetNameTitle(binNameSR.str().c_str(), "data in signal region");
	dataSBL->SetNameTitle(binNameSBL.str().c_str(), "data in LSB");
	dataSBR->SetNameTitle(binNameSBR.str().c_str(), "data in RSB");

	ws->import(*dataSR);
	ws->import(*dataSBL);
	ws->import(*dataSBR);

	std:: cout << "----------------------------" << "\n"
		<< "events in SR: " << dataSR->numEntries() << "\n"
		<< "events in LSB: " << dataSBL->numEntries() << "\n"
		<< "events in RSB: " << dataSBR->numEntries() << "\n"
		<< "----------------------------" << std::endl;

	// caculating median of different regions by filling events into histogram and getting the mean
	TH1* histSR =  dataSR->createHistogram("histSR", *JpsiMass,  Binning(120));
	TH1* histSBL = dataSBL->createHistogram("histSBL", *JpsiMass, Binning(120));
	TH1* histSBR = dataSBR->createHistogram("histSBR", *JpsiMass, Binning(120));

	double meanSR = histSR->GetMean();
	double meanSBL = histSBL->GetMean();
	double meanSBR = histSBR->GetMean();

	RooRealVar* MeanSR = new RooRealVar("MeanSR","MeanSR",meanSR);
	RooRealVar* MeanSBL = new RooRealVar("MeanSBL","MeanSBL",meanSBL);
	RooRealVar* MeanSBR = new RooRealVar("MeanSBR","MeanSBR",meanSBR);
	ws->import(RooArgList(*MeanSR, *MeanSBL, *MeanSBR));

	std::cout << "----------------------------" << "\n"
		<< "meanSR: " << meanSR << "\n"
		<< "meanSBL: " << meanSBL << "\n"
		<< "meanSBR: " << meanSBR 
		<< "----------------------------" << std::endl;

	// calculate background fraction in signal region
	double BkgRatio3Sig = getFracBkgIn3Sigma(rapBin, ptBin, cpmBin, infilename.c_str(), onia::nSigMass);

	double BkgRatio3SigErr = getFracBkgErrIn3Sigma(rapBin, ptBin, cpmBin, infilename.c_str(), onia::nSigMass);

	std::cout << "------------ background fraction ----------------" << "\n"
		<< "background fraction in "<<onia::nSigMass<<" sigma mass window: " << BkgRatio3Sig << "\n"
		<< "error on background fraction: " << BkgRatio3SigErr << "\n"
		<< "-------------------------------------------------" << std::endl;

	RooRealVar* FracBkg = new RooRealVar("FracBkg","FracBkg",BkgRatio3Sig);
	RooRealVar* FracBkgErr = new RooRealVar("FracBkgErr","FracBkgErr",BkgRatio3SigErr);
	ws->import(RooArgList(*FracBkg, *FracBkgErr));

	//**** added 2013-3-4
	RooAbsPdf *bkgMassShape = (RooAbsPdf*)ws->pdf("bkgMassShape");
	RooAddPdf *sigMassShape = (RooAddPdf*)ws->pdf("sigMassShape");
	RooAddPdf *massPdf = (RooAddPdf*)ws->pdf("massModel");
	JpsiMass->setRange("SR",sigMinMass,sigMaxMass);
	JpsiMass->setRange("SBL",onia::massMin, sbLowMass);
	JpsiMass->setRange("SBR",sbHighMass, onia::massMax);
	RooRealVar *fracBkgSBL = (RooRealVar*)bkgMassShape->createIntegral(*JpsiMass,NormSet(*JpsiMass),Range("SBL"));
	RooRealVar *fracBkgSBR = (RooRealVar*)bkgMassShape->createIntegral(*JpsiMass,NormSet(*JpsiMass),Range("SBR"));
	RooRealVar *fullInSBL  = (RooRealVar*)massPdf->createIntegral(*JpsiMass,NormSet(*JpsiMass),Range("SBL"));
	RooRealVar *fullInSBR  = (RooRealVar*)massPdf->createIntegral(*JpsiMass,NormSet(*JpsiMass),Range("SBR"));
	double fracBkgSBL_val = fracBkgSBL->getVal();
	double fracBkgSBR_val = fracBkgSBR->getVal();
	double fullInSBL_val  = fullInSBL->getVal();
	double fullInSBR_val  = fullInSBR->getVal();
	double fracBkgInSBL = fracBkgSBL_val * fracBkg / fullInSBL_val ; 
	double fracBkgInSBR = fracBkgSBR_val * fracBkg / fullInSBR_val ;
	double fracBkgInSBLErr = fracBkgSBL_val / fullInSBL_val * fracBkgErr ;
	double fracBkgInSBRErr = fracBkgSBR_val / fullInSBR_val * fracBkgErr ;

	RooRealVar* FracBkgSBL = new RooRealVar("FracBkgSBL","FracBkgSBL",fracBkgInSBL);
	RooRealVar* FracBkgSBLErr = new RooRealVar("FracBkgSBLErr","FracBkgSBLErr",fracBkgInSBLErr);
	ws->import(RooArgList(*FracBkgSBL, *FracBkgSBLErr));
	RooRealVar* FracBkgSBR = new RooRealVar("FracBkgSBR","FracBkgSBR",fracBkgInSBR);
	RooRealVar* FracBkgSBRErr = new RooRealVar("FracBkgSBRErr","FracBkgSBRErr",fracBkgInSBRErr);
	ws->import(RooArgList(*FracBkgSBR, *FracBkgSBRErr));
	cout<<"fracBkgInSBL: "<<fracBkgInSBL<<" fracBkgInSBLErr: "<<fracBkgInSBLErr<<endl;
	cout<<"fracBkgInSBR: "<<fracBkgInSBR<<" fracBkgInSBRErr: "<<fracBkgInSBRErr<<endl;
	//****

	// building lifetime pdf
	std::cout << ">>>Building Mass and LifeTime PDF" << std::endl;
	buildLifetimePDF(ws);

	// fitting
	std::cout << ">>>Fitting" << std::endl;
	doFit(ws, nState, BkgRatio3Sig, fracBkgInSBL, fracBkgInSBR, rapBin, ptBin, cpmBin);

	std::cout << ">>>Writing results to root file" << std::endl;
	ws->Write();
	infile->Close();
}

//=================================================================================

void buildLifetimePDF(RooWorkspace *ws){

	//----prompt
	//resolution function
	ws->factory("RooGaussModel::promptLifetime(Jpsict,promptMean[0,-.01,.01],ctResolution[1,.001,2], 1, JpsictErr)");
	((RooGaussModel*)ws->pdf("promptLifetime"))->advertiseFlatScaleFactorIntegral(true);

	ws->factory("RooGaussModel::promptLifetime2(Jpsict,promptMean, ctResolution2[1,.001,6], 1, JpsictErr)");
	((RooGaussModel*)ws->pdf("promptLifetime2"))->advertiseFlatScaleFactorIntegral(true);

	RooGaussModel* promptLifetime = (RooGaussModel*)ws->pdf("promptLifetime");
	RooGaussModel* promptLifetime2 = (RooGaussModel*)ws->pdf("promptLifetime2");
	RooRealVar fracGauss2("fracGauss2","fracGauss2",.01,.0,.45);
	RooAddModel TotalPromptLifetime("TotalPromptLifetime","TotalPromptLifetime",
			RooArgSet(*promptLifetime2,*promptLifetime),fracGauss2);
	ws->import(fracGauss2);
	ws->import(TotalPromptLifetime);
	TotalPromptLifetime.Print();

	//----noprompt
	ws->factory("RooDecay::nonPromptSSD(Jpsict,nonPromptTau[.3,.01,3],TotalPromptLifetime,RooDecay::SingleSided)");

	//----background
	//SBL
	ws->factory("RooDecay::backgroundSSD_SBL(Jpsict,bkgTauSSD_SBL[.4,0,3],TotalPromptLifetime,RooDecay::SingleSided)");
	
	ws->factory("RooDecay::backgroundFD(Jpsict,bkgTauFD[.2,.00001,3],TotalPromptLifetime,RooDecay::Flipped)");
	
	
	ws->factory("RooDecay::backgroundDSD(Jpsict,bkgTauDSD[.05,0,3],TotalPromptLifetime,RooDecay::DoubleSided)");
	ws->factory("SUM::backgroundlifetimeLpre(fBkgSSDR_SBL[.4,0,1.]*backgroundSSD_SBL,fBkgDSD_SBL[.2,0,1.]*backgroundDSD,backgroundFD)");
	////SBR
	ws->factory("RooDecay::backgroundSSD_SBR(Jpsict,bkgTauSSD_SBR[.4,0,3],TotalPromptLifetime,RooDecay::SingleSided)");
	ws->factory("SUM::backgroundlifetimeRpre(fBkgSSDR_SBR[.4,0,1.]*backgroundSSD_SBR,fBkgDSD_SBR[.2,0,1.]*backgroundDSD,backgroundFD)");
	//Signal region
	//interpolation
	ws->factory("expr::fBkgSSDR('@0+(@2-@3)*(@1-@0)/(@4-@3)',fBkgSSDR_SBL,fBkgSSDR_SBR,MeanSR,MeanSBL,MeanSBR)");
	ws->factory("expr::fBkgDSD('@0+(@2-@3)*(@1-@0)/(@4-@3)',fBkgDSD_SBL,fBkgDSD_SBR,MeanSR,MeanSBL,MeanSBR)");
	ws->factory("expr::bkgTauSSD('@0+(@2-@3)*(@1-@0)/(@4-@3)',bkgTauSSD_SBL,bkgTauSSD_SBR,MeanSR,MeanSBL,MeanSBR)");
	ws->factory("RooDecay::backgroundSSD(Jpsict,bkgTauSSD,TotalPromptLifetime,RooDecay::SingleSided)");
	ws->factory("SUM::backgroundlifetime(fBkgSSDR*backgroundSSD,fBkgDSD*backgroundDSD,backgroundFD)");
	//---final pdf
	//signal region
	ws->factory("SUM::fulllifetime(fBkg[0.5,0.,1.]*backgroundlifetime,fPrompt[0.5,0.,1.]*TotalPromptLifetime,nonPromptSSD)");
	ws->factory("RooGaussian::fBkgConstraint(fBkg, FracBkg, FracBkgErr)"); 
	ws->factory("PROD::lifetimeConstraint(fulllifetime, fBkgConstraint)"); 

	RooRealVar* fBkgSBL = new RooRealVar("fBkgSBL","fBkgSBL",0.8,0.,1.);
	RooRealVar* fBkgSBR = new RooRealVar("fBkgSBR","fBkgSBR",0.8,0.,1.);
	ws->import(*fBkgSBL); ws->import(*fBkgSBR);

	//f_P / f_NP should be same in LSB,SR,RSB, to interpolate f_P in L(R)SB from SR
	ws->factory("expr::fPromptSBL('@0*(1.-@1)/(1.-@2)',fPrompt,fBkgSBL,fBkg)");
	ws->factory("expr::fPromptSBR('@0*(1.-@1)/(1.-@2)',fPrompt,fBkgSBR,fBkg)");

	//SBL
	ws->factory("SUM::backgroundlifetimeLnoC(fBkgSBL*backgroundlifetimeLpre,fPromptSBL*TotalPromptLifetime,nonPromptSSD)");
	ws->factory("RooGaussian::fBkgSBLConstraint(fBkgSBL, FracBkgSBL, FracBkgSBLErr)");
	ws->factory("PROD::backgroundlifetimeL(backgroundlifetimeLnoC, fBkgSBLConstraint)"); 

	//SBR
	ws->factory("SUM::backgroundlifetimeRnoC(fBkgSBR*backgroundlifetimeRpre,fPromptSBR*TotalPromptLifetime,nonPromptSSD)");
	ws->factory("RooGaussian::fBkgSBRConstraint(fBkgSBR, FracBkgSBR, FracBkgSBRErr)");
	ws->factory("PROD::backgroundlifetimeR(backgroundlifetimeRnoC, fBkgSBRConstraint)"); 

	ws->Print("v");
}


//=================================================================================
void doFit(RooWorkspace *ws, int nState, double BkgRatio3Sig, double fracBkgInSBL, double fracBkgInSBR, int rapBin, int ptBin, int cpmBin){
	RooRealVar JpsiMass(*ws->var("JpsiMass"));
	RooRealVar Jpsict(*ws->var("Jpsict"));
	stringstream binNameSR, binNameSBL, binNameSBR;
	binNameSR  << "data_rap" << rapBin << "_pt" << ptBin << "_cpm" << cpmBin <<"_SR";
	binNameSBL << "data_rap" << rapBin << "_pt" << ptBin << "_cpm" << cpmBin <<"_SBL";
	binNameSBR << "data_rap" << rapBin << "_pt" << ptBin << "_cpm" << cpmBin <<"_SBR";
	RooDataSet *dataSR = (RooDataSet*)ws->data(binNameSR.str().c_str());
	RooDataSet *dataSBL = (RooDataSet*)ws->data(binNameSBL.str().c_str());
	RooDataSet *dataSBR = (RooDataSet*)ws->data(binNameSBR.str().c_str());
	cout<<"bin name is "<<binNameSR<<endl;

	if(nState==4){
		ws->var("fBkgSSDR_SBL")->setVal(.7); //0.4
		ws->var("fBkgSSDR_SBR")->setVal(.7); //0.4
		ws->var("fBkgDSD_SBL")->setVal(.2);
		ws->var("fBkgDSD_SBR")->setVal(.2);
		if(ptBin>7){
			//// old
			//ws->var("fBkgSSDR_SBL")->setVal(.77);
			//ws->var("fBkgSSDR_SBR")->setVal(.77);
			//ws->var("fBkgSSDR_SBL")->setConstant(kTRUE);
			//ws->var("fBkgSSDR_SBR")->setConstant(kTRUE);
			//ws->var("fBkgDSD_SBL")->setMax(.23);
			//ws->var("fBkgDSD_SBR")->setMax(.23);
			////new model
			ws->var("fBkgSSDR_SBL")->setVal(.785);
			ws->var("fBkgSSDR_SBR")->setVal(.785);
			ws->var("fBkgSSDR_SBL")->setConstant(kTRUE);
			ws->var("fBkgSSDR_SBR")->setConstant(kTRUE);
			ws->var("fBkgDSD_SBL")->setMax(.215);
			ws->var("fBkgDSD_SBR")->setMax(.215);
		}
		ws->var("bkgTauSSD_SBL")->setVal(.4);
		ws->var("bkgTauSSD_SBR")->setVal(.4);
		ws->var("bkgTauFD")->setVal(.1);
		ws->var("bkgTauDSD")->setVal(.01);

		ws->var("fBkg")->setVal(BkgRatio3Sig);
		ws->var("fPrompt")->setVal(.3);
		ws->var("fBkgSBL")->setVal(fracBkgInSBL);
		ws->var("fBkgSBR")->setVal(fracBkgInSBR);
		ws->var("nonPromptTau")->setVal(.4);

		//** set limit
		ws->var("bkgTauSSD_SBL")->setMax(1.);
		ws->var("bkgTauSSD_SBR")->setMax(1.);
		ws->var("bkgTauFD")->setMax(0.3);
		//ws->var("bkgTauDSD")->setMax(0.04); //0.05
		ws->var("nonPromptTau")->setMax(1.);
		//ws->var("fBkgSBL")->setMin(0.75);
		//ws->var("fBkgSBR")->setMin(0.84);
		if(ptBin>10){
			ws->var("bkgTauDSD")->setMax(0.005);
			ws->var("bkgTauDSD")->setVal(0.);
		}
		//**
		ws->var("promptMean")->setVal(0.);
		ws->var("ctResolution")->setVal(.9);
		ws->var("promptMean")->setConstant(kTRUE);
		ws->var("ctResolution")->setConstant(kTRUE);
		
		if(ptBin==2){
		ws->var("bkgTauSSD_SBR")->setVal(.0); 
/*		if(cpmBin==2){
			ws->var("fBkgSSDR_SBR")->setVal(.4);
			}
*/		 if(cpmBin==1){ws->var("bkgTauSSD_SBR")->setVal(.05);}
		if(cpmBin==3){ws->var("bkgTauSSD_SBR")->setVal(.05);}
		
		if(cpmBin>4){
		ws->var("bkgTauSSD_SBR")->setVal(.4);
		}
		}

		if(rapBin == 1)
			ws->var("ctResolution2")->setVal(1.1);
		else if(rapBin == 2)
			ws->var("ctResolution2")->setVal(1.5);
		ws->var("ctResolution2")->setConstant(kTRUE);
		ws->var("fracGauss2")->setVal(0.2);
		ws->var("fracGauss2")->setMin(0.04);
	}
	else if(nState==5){
		ws->var("fBkgSSDR_SBL")->setVal(.6);
		ws->var("fBkgSSDR_SBR")->setVal(.3);
		ws->var("fBkgDSD_SBL")->setVal(.4);
		ws->var("fBkgDSD_SBR")->setVal(.6);

		ws->var("bkgTauSSD_SBL")->setVal(.4);
		ws->var("bkgTauSSD_SBR")->setVal(.4);
		ws->var("bkgTauFD")->setVal(.1);
		ws->var("bkgTauDSD")->setVal(.01);

		ws->var("fBkg")->setVal(BkgRatio3Sig);
		ws->var("fPrompt")->setVal(.3);
		ws->var("fBkgSBL")->setVal(fracBkgInSBL);
		ws->var("fBkgSBR")->setVal(fracBkgInSBR);

		ws->var("nonPromptTau")->setVal(.4);

		//** set limit
		ws->var("bkgTauSSD_SBL")->setMax(1.);
		ws->var("bkgTauSSD_SBR")->setMax(1.);
		ws->var("bkgTauFD")->setMax(0.3);
		ws->var("bkgTauDSD")->setMax(0.05);
		ws->var("nonPromptTau")->setMax(1.);
		ws->var("fBkgDSD_SBL")->setMin(0.2);
		ws->var("fBkgDSD_SBR")->setMin(0.2);
		//**
		

		ws->var("promptMean")->setVal(0.);
		ws->var("ctResolution")->setVal(.9);
		ws->var("promptMean")->setConstant(kTRUE);
		ws->var("ctResolution")->setConstant(kTRUE);
		ws->var("ctResolution2")->setVal(3.);
		ws->var("ctResolution2")->setConstant(kTRUE);
		ws->var("fracGauss2")->setVal(0.01);
	}
	RooAbsPdf *ModelLifeSR = (RooAbsPdf*)ws->pdf("lifetimeConstraint");
	RooAbsPdf *ModelLifeSBL = (RooAbsPdf*)ws->pdf("backgroundlifetimeL");
	RooAbsPdf *ModelLifeSBR = (RooAbsPdf*)ws->pdf("backgroundlifetimeR");

	RooArgSet *NLLs = new RooArgSet();
	RooAbsReal *MLNLLSR = NULL, *MLNLLSBL = NULL, *MLNLLSBR = NULL;
	MLNLLSR = (RooAbsReal *)ModelLifeSR->createNLL(*dataSR,
			ConditionalObservables(RooArgSet(*ws->var("JpsictErr"))),
			Constrain(RooArgSet(*ws->var("fBkg"))),
			Extended(kFALSE),
			NumCPU(6));
	MLNLLSBL = (RooAbsReal *)ModelLifeSBL->createNLL(*dataSBL,
			ConditionalObservables(RooArgSet(*ws->var("JpsictErr"))),
			Constrain(RooArgSet(*ws->var("fBkgSBL"))),
			Extended(kFALSE),
			NumCPU(6));
	MLNLLSBR = (RooAbsReal *)ModelLifeSBR->createNLL(*dataSBR,
			ConditionalObservables(RooArgSet(*ws->var("JpsictErr"))),
			Constrain(RooArgSet(*ws->var("fBkgSBR"))),
			Extended(kFALSE),
			NumCPU(6));
	NLLs->add(*MLNLLSR);
	NLLs->add(*MLNLLSBL);
	NLLs->add(*MLNLLSBR);

	RooAddition *simNLL = new RooAddition("add","add",*NLLs);
	RooMinuit *lMinuit = new RooMinuit(*simNLL);
	lMinuit->setStrategy(1); 
	if(nState==4 && ptBin > 7) lMinuit->setStrategy(2); 
	lMinuit->setPrintEvalErrors(-1);
	lMinuit->setEvalErrorWall(false);
	lMinuit->setVerbose(false);
	lMinuit->setPrintLevel(-1);
	stringstream rltName, snapshotName;
	rltName<< "l_fitresult_rap"<<rapBin<<"_pt"<<ptBin<<"_cpm"<<cpmBin;
	snapshotName<< "l_snapshot_rap"<<rapBin<<"_pt"<<ptBin<<"_cpm"<<cpmBin;
	RooFitResult *fitresult;
	lMinuit->simplex();
	lMinuit->migrad();
	fitresult = (RooFitResult*)lMinuit->save(rltName.str().c_str());
	cout<<"fitresult->covQual(): "<<fitresult->covQual()<<endl;
	lMinuit->migrad();
	fitresult = (RooFitResult*)lMinuit->save(rltName.str().c_str());
	cout<<"fitresult->covQual(): "<<fitresult->covQual()<<endl;
	lMinuit->hesse();
	fitresult = (RooFitResult*)lMinuit->save(rltName.str().c_str());
	cout<<"fitresult->covQual(): "<<fitresult->covQual()<<endl;

	double BkgTauFD = ((RooRealVar*)ws->var("bkgTauFD"))->getVal();
	double BkgTauDSD = ((RooRealVar*)ws->var("bkgTauDSD"))->getVal();
	double CtResolution = ((RooRealVar*)ws->var("ctResolution"))->getVal();
	double CtResolution2 = ((RooRealVar*)ws->var("ctResolution2"))->getVal();
	int refitCount = 0;
	while(BkgTauFD < BkgTauDSD){
		if(refitCount>3) break;
		cout<<"bkgTauFD < bkgTauDSD"<<endl;
		ws->var("bkgTauFD")->setVal(BkgTauDSD);
		ws->var("bkgTauDSD")->setVal(BkgTauFD);
		lMinuit->migrad();
		fitresult = (RooFitResult*)lMinuit->save(rltName.str().c_str());
		cout<<"fitresult->covQual(): "<<fitresult->covQual()<<endl;
		lMinuit->migrad();
		fitresult = (RooFitResult*)lMinuit->save(rltName.str().c_str());
		cout<<"fitresult->covQual(): "<<fitresult->covQual()<<endl;
		lMinuit->hesse();
		fitresult = (RooFitResult*)lMinuit->save(rltName.str().c_str());
		cout<<"fitresult->covQual(): "<<fitresult->covQual()<<endl;

		BkgTauFD = ((RooRealVar*)ws->var("bkgTauFD"))->getVal();
		BkgTauDSD = ((RooRealVar*)ws->var("bkgTauDSD"))->getVal();
		refitCount++;
	}
	while(CtResolution > CtResolution2){
		if(refitCount>3) break;
		cout<<"CtResolution > CtResolution2"<<endl;
		ws->var("ctResolution")->setVal(CtResolution2);
		ws->var("ctResolution2")->setVal(CtResolution);
		lMinuit->migrad();
		fitresult = (RooFitResult*)lMinuit->save(rltName.str().c_str());
		cout<<"fitresult->covQual(): "<<fitresult->covQual()<<endl;
		lMinuit->migrad();
		fitresult = (RooFitResult*)lMinuit->save(rltName.str().c_str());
		cout<<"fitresult->covQual(): "<<fitresult->covQual()<<endl;
		lMinuit->hesse();
		fitresult = (RooFitResult*)lMinuit->save(rltName.str().c_str());
		cout<<"fitresult->covQual(): "<<fitresult->covQual()<<endl;

		CtResolution = ((RooRealVar*)ws->var("ctResolution"))->getVal();
		CtResolution2 = ((RooRealVar*)ws->var("ctResolution2"))->getVal();

		refitCount++;
	}

	cout<<">>>>>>Refitted "<<refitCount<<" times"<<endl;
	if(refitCount>3) cout<<">>>>>>Fit not converged"<<endl;

	fitresult->Print();
	ws->import(*fitresult,kTRUE);

	ws->saveSnapshot(snapshotName.str().c_str(), ws->allVars());
}

