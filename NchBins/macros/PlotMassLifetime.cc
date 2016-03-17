#include "rootIncludes.inc"
#include "commonVar.h"
#include "RooUtils.h"

using namespace RooFit;

void plotMass(RooWorkspace *ws, int rapBin, int ptBin, int cpmBin, int nState);
void plotMassLog(RooWorkspace *ws, int rapBin, int ptBin, int cpmBin, int nState);
void plotLifeBg(RooWorkspace *ws, int rapBin, int ptBin, int cpmBin, int nState);
void plotLifeSig(RooWorkspace *ws, int rapBin, int ptBin, int cpmBin, int nState);
void plotLifeSig_linear(RooWorkspace *ws, int rapBin, int ptBin, int cpmBin, int nState);
//separate pull and distribution
void plotLifeBgIndividual(RooWorkspace *ws, int rapBin, int ptBin, int cpmBin, int nState);
void plotLifeSigIndividual(RooWorkspace *ws, int rapBin, int ptBin, int cpmBin, int nState);
//==============================================

void PlotMassLifetime(const std::string &infilename, int rapBin, int ptBin, int cpmBin, int nState, int Plotting){
	RooWorkspace* ws = getFromTFile<RooWorkspace>(infilename, "ws_masslifetime");

	switch (Plotting) {
		case 1:
			std::cout << ">>>>Plotting mass" << std::endl;
			plotMass(ws, rapBin, ptBin, cpmBin, nState);
			plotMassLog(ws, rapBin, ptBin, cpmBin, nState);
			std::cout << ">>>>Plotting lifetime sidebands" << std::endl;
			plotLifeBg(ws, rapBin, ptBin, cpmBin, nState);
			std::cout << ">>>>Plotting lifetime signal region" << std::endl;
			plotLifeSig(ws, rapBin, ptBin, cpmBin, nState);
			plotLifeSig_linear(ws, rapBin, ptBin, cpmBin, nState);
			break;
		case 2:
			std::cout << ">>>>Plotting mass" << std::endl;
			plotMass(ws, rapBin, ptBin, cpmBin, nState);
			plotMassLog(ws, rapBin, ptBin, cpmBin, nState);
			break;
		case 3:
			std::cout << ">>>>Plotting lifetime sidebands" << std::endl;
			plotLifeBg(ws, rapBin, ptBin, cpmBin, nState);
			break;
		case 4:
			std::cout << ">>>>Plotting lifetime signal region" << std::endl;
			plotLifeSig(ws, rapBin, ptBin, cpmBin, nState);
			plotLifeSig_linear(ws, rapBin, ptBin, cpmBin, nState);
			break;
		case 5:
			std::cout << ">>>>Plotting lifetime sidebands, separate pull and distribution" << std::endl;
			plotLifeBgIndividual(ws, rapBin, ptBin, cpmBin, nState);
			break;
		case 6:
			std::cout << ">>>>Plotting lifetime signal region, separate pull and distribution" << std::endl;
			plotLifeSigIndividual(ws, rapBin, ptBin, cpmBin, nState);
			break;
		default:
			std::cerr << "I do not know what do do with this value of Plotting" << std::endl;
	}
	delete ws;
}

//==============================================
void plotMass(RooWorkspace *ws, int rapBin, int ptBin, int cpmBin, int nState){
	int  nbins=90; //0.005 bin size
	TGaxis::SetMaxDigits(3);

	if(nState == 4) nbins=90;
	if(nState == 5) nbins=120;

	RooRealVar *JpsiMass = ws->var("JpsiMass");
	assert( 0 != JpsiMass );

	RooPlot *massFrame = JpsiMass->frame(Bins(nbins));
	assert ( 0 != massFrame );
	massFrame->SetName(Form("mass_plot_rap%d_pt%d_cpm%d",rapBin,ptBin,cpmBin));
	massFrame->SetTitle("");
	massFrame->GetYaxis()->SetTitle("Events / 5 MeV");
	massFrame->GetYaxis()->SetTitleOffset(1.3);

	RooPlot *massFramePull = JpsiMass->frame(Bins(nbins));
	assert ( 0 != massFramePull );
	massFramePull->SetName(Form("pullmass_plot_rap%d_pt%d_cpm%d",rapBin,ptBin,cpmBin));
	massFramePull->SetTitle("");
	massFramePull->GetYaxis()->SetTitle("pull");
	massFramePull->GetXaxis()->SetTitleSize(0.08);
	massFramePull->GetYaxis()->SetTitleSize(0.08);
	massFramePull->GetXaxis()->SetLabelSize(0.08);
	massFramePull->GetYaxis()->SetLabelSize(0.08);
	massFramePull->GetYaxis()->SetTitleOffset(0.4);
	massFramePull->GetYaxis()->SetRangeUser(-5.5,5.5);

	RooAbsData *data= ws->data(Form("data_rap%d_pt%d_cpm%d",rapBin,ptBin,cpmBin));
	assert ( 0 != data );
	RooFitResult* fitRlt = dynamic_cast<RooFitResult*>(ws->obj(Form("m_fitresult_rap%d_pt%d_cpm%d",rapBin,ptBin,cpmBin)));
	assert ( 0 != fitRlt);

	double Mean = -1.0, MeanErr = -1.0;
	getVarFromWorkspace(ws, "CBmass", Mean, MeanErr);
	double Sigma = -1.0, SigmaErr = -1.0;
	getVarFromWorkspace(ws, "CBsigma", Sigma, SigmaErr);
	double Sigma2 = -1.0, Sigma2Err = -1.0;
	getVarFromWorkspace(ws, "CBsigma2", Sigma2, Sigma2Err);
	double Alpha = -1.0, AlphaErr = -1.0;
	getVarFromWorkspace(ws, "CBalpha", Alpha, AlphaErr);
	double cbN = -1.0, cbNErr = -1.0;
	getVarFromWorkspace(ws, "CBn", cbN, cbNErr);
	double lambda = -1.0, lambdaErr = -1.0;
	getVarFromWorkspace(ws, "bkgLambda", lambda, lambdaErr);
	double fracCB1 = -1.0, fracCB1Err = -1.0;
	getVarFromWorkspace(ws, "fracCB1", fracCB1, fracCB1Err);
	double fracBkg = -1.0, fracBkgErr = -1.0;
	getVarFromWorkspace(ws, "fracBkg", fracBkg, fracBkgErr);

	double SigmaWei = sqrt(pow(Sigma,2)*fracCB1+pow(Sigma2,2)*(1-fracCB1));
	double SigmaWeiErr =  (1./(2*SigmaWei))*
		sqrt(pow((pow(Sigma,2)-pow(Sigma2,2))*fracCB1Err,2) + pow(2*fracCB1*Sigma*SigmaErr,2) + pow(2*(1-fracCB1)*Sigma2*Sigma2Err,2));

	RooAbsPdf *massPdf = ws->pdf("massModel");
	assert ( 0 != massPdf );
	RooAbsPdf *bkgMassShape = ws->pdf("bkgMassShape");
	assert ( 0 != bkgMassShape );
	double sigMaxMass = Mean+SigmaWei*onia::nSigMass;
	double sigMinMass = Mean-SigmaWei*onia::nSigMass;
	double sbHighMass = Mean+SigmaWei*onia::nSigBkgHigh;
	double sbLowMass =  Mean-SigmaWei*onia::nSigBkgLow;
	cout<<"sigMaxMass: "<<sigMaxMass<<endl;
	cout<<"sigMinMass: "<<sigMinMass<<endl;
	cout<<"sbHighMass: "<<sbHighMass<<endl;
	cout<<"sbLowMass: "<<sbLowMass<<endl;

	int nEntries = data->numEntries();
	JpsiMass->setRange("SigRegion",sigMinMass,sigMaxMass);

	RooRealVar* fracFull3Sig=(RooRealVar*)massPdf->createIntegral(*JpsiMass,NormSet(*JpsiMass),Range("SigRegion"));
	RooRealVar* fracBkg3Sig=(RooRealVar*)bkgMassShape->createIntegral(*JpsiMass,NormSet(*JpsiMass),Range("SigRegion"));
	double evtFull3Sig = nEntries*fracFull3Sig->getVal();
	double evtBkg3Sig = nEntries*fracBkg*fracBkg3Sig->getVal();
	double BkgRatio3Sig = evtBkg3Sig/evtFull3Sig;

	cout<<"fracFull3Sig: "<<fracFull3Sig->getVal()<<endl;
	cout<<"fracBkg3Sig: "<<fracBkg3Sig->getVal()<<endl;
	cout<<"evtFull3Sig: "<<evtFull3Sig<<endl;
	cout<<"evtBkg3Sig: "<<evtBkg3Sig<<endl;
	cout<<"BkgRatio3Sig: "<<BkgRatio3Sig<<endl;

	data->plotOn(massFrame,MarkerSize(0.8));
	massPdf->plotOn(massFrame,
			LineWidth(2),
			ProjWData(*data));
	//------get chi2------------SHOULD DONE after PLOTTING------
	int parsFit=(fitRlt->floatParsFinal()).getSize(); //this used the full p.d.f.
	int nBins_Mass=massFrame->GetNbinsX();
	double chi2Pre_Mass=massFrame->chiSquare(parsFit);  //reduced chi-squared = chi2/ndof
	int ndof_Mass=nBins_Mass-parsFit;  //num of degree of freedom
	double chi2_Mass=chi2Pre_Mass*ndof_Mass;

	RooHist* hpull_mass = massFrame->pullHist() ;
	hpull_mass->SetMarkerSize(0.8);
	for(int i=0;i<hpull_mass->GetN();i++){
		hpull_mass->SetPointEYlow(i,0.);
		hpull_mass->SetPointEYhigh(i,0.);
	}
	massFramePull->addPlotable(hpull_mass,"P");

	massPdf->plotOn(massFrame,
			Components("bkgMassShape"),
			LineStyle(7),
			LineColor(kPink+3),
			LineWidth(2),
			ProjWData(*data));

	if(Sigma > Sigma2){

		fracCB1 = 1-fracCB1;
		double temp = 0.;
		temp = Sigma; Sigma = Sigma2; Sigma2 = temp;
		temp = SigmaErr; SigmaErr = Sigma2Err; Sigma2Err = temp;
		massPdf->plotOn(massFrame,
				Components("massCBShape"),
				LineStyle(2),
				LineColor(kGreen),
				LineWidth(2),
				ProjWData(*data));
		massPdf->plotOn(massFrame,
				Components("massCBShape2"),
				LineStyle(2),
				LineColor(kRed),
				LineWidth(2),
				ProjWData(*data));

	}  else{
		massPdf->plotOn(massFrame,
				Components("massCBShape"),
				LineStyle(2),
				LineColor(kRed),
				LineWidth(2),
				ProjWData(*data));
		massPdf->plotOn(massFrame,
				Components("massCBShape2"),
				LineStyle(2),
				LineColor(kGreen),
				LineWidth(2),
				ProjWData(*data));
	}


	double minY = 0.;

	double maxY = 0.;
	if(nState == 4) maxY = massFrame->GetMaximum()*0.3;
	if(nState == 5) maxY = massFrame->GetMaximum()*0.4;
	double lineWidth = 2.0;
	TLine *lineSBLow = new TLine(sbLowMass, minY, sbLowMass, maxY);
	TLine *lineSBHigh = new TLine(sbHighMass, minY, sbHighMass, maxY);
	TLine *lineSigLow = new TLine(sigMinMass, minY, sigMinMass, maxY);
	TLine *lineSigHigh = new TLine(sigMaxMass, minY, sigMaxMass, maxY);
	lineSBLow->SetLineWidth(lineWidth);lineSBHigh->SetLineWidth(lineWidth);
	lineSigLow->SetLineWidth(lineWidth);lineSigHigh->SetLineWidth(lineWidth);
	lineSBLow->SetLineColor(kBlue);lineSBHigh->SetLineColor(kBlue);
	lineSigLow->SetLineColor(kRed);lineSigHigh->SetLineColor(kRed);
	lineSBLow->SetLineStyle(7);lineSBHigh->SetLineStyle(7);
	lineSigLow->SetLineStyle(5);lineSigHigh->SetLineStyle(5);

	TH1* legendBlue = data->createHistogram("legendBlue",*JpsiMass,Binning(50)) ; legendBlue->SetLineColor(kBlue) ; legendBlue->SetLineStyle(kSolid) ; legendBlue->SetLineWidth(2) ;
	TH1* legendBlueDash = data->createHistogram("legendBlueDash",*JpsiMass,Binning(50)) ; legendBlueDash->SetLineColor(kBlue) ; legendBlueDash->SetLineStyle(5) ; legendBlueDash->SetLineWidth(2) ;
	TH1* legendRed = data->createHistogram("legendRed",*JpsiMass,Binning(50)) ; legendRed->SetLineColor(kRed) ; legendRed->SetLineStyle(2) ; legendRed->SetLineWidth(2) ;
	TH1* legendBlack = data->createHistogram("legendBlack",*JpsiMass,Binning(50)) ; legendBlack->SetLineColor(kBlack) ; legendBlack->SetLineStyle(kSolid) ; legendBlack->SetLineWidth(2) ;
	TH1* legendGreen = data->createHistogram("legendGreen",*JpsiMass,Binning(50)) ; legendGreen->SetLineColor(kGreen) ; legendGreen->SetLineStyle(2) ; legendGreen->SetLineWidth(2) ;
	TH1* legendGreenDash = data->createHistogram("legendGreenDash",*JpsiMass,Binning(50)) ; legendGreenDash->SetLineColor(kGreen) ; legendGreenDash->SetLineStyle(2) ; legendGreenDash->SetLineWidth(2) ;
	TH1* legendPink = data->createHistogram("legendPink",*JpsiMass,Binning(50)) ; legendPink->SetLineColor(kPink+3) ; legendPink->SetLineStyle(7) ; legendPink->SetLineWidth(2) ;

	TLegend* MassLegend=new TLegend(0.7,0.55,0.88,0.75);
	MassLegend->SetFillColor(kWhite);
	MassLegend->SetTextFont(42);
	MassLegend->SetTextSize(0.035);
	MassLegend->SetBorderSize(0.);
	MassLegend->AddEntry(legendBlue,"sum","l");
	MassLegend->AddEntry(legendRed,"signal CB_{1}","l");
	MassLegend->AddEntry(legendGreen,"signal CB_{2}","l");
	MassLegend->AddEntry(legendPink,"background","l");

	double left=0.7, top=0.9, textSize=0.025;
	TLatex *latex=new TLatex();
	latex->SetTextFont(42);
	latex->SetNDC(kTRUE);
	latex->SetTextSize(textSize);
	double step=textSize*1.3;

	gStyle->SetPadBottomMargin(0.08); //0.12
	gStyle->SetPadLeftMargin(0.09); //0.12
	gStyle->SetPadRightMargin(0.02); //0.05
	gStyle->SetPadTopMargin(0.05); //0.05
	
	TCanvas *c1=new TCanvas("c1","",800,700);

	massFrame->Draw(); MassLegend->Draw();
	lineSBLow->Draw("same"); lineSBHigh->Draw("same"); lineSigLow->Draw("same"); lineSigHigh->Draw("same");
	top=0.90; textSize=0.030; latex->SetTextSize(textSize);
	if(nState == 4)
		latex->DrawLatex(left,top,"J/#psi");
	if(nState == 5)
		latex->DrawLatex(left,top,"#psi(2S)");
	top-=step;
	textSize=0.025; latex->SetTextSize(textSize);


	if(rapBin==0) latex->DrawLatex(left,top,Form("%.1f < |y| < %.1f",onia::rapForPTRange[rapBin],onia::rapForPTRange[2]));
	else if(rapBin==1) latex->DrawLatex(left,top,Form("|y| < %.1f",onia::rapForPTRange[rapBin]));
	else latex->DrawLatex(left,top,Form("%.1f < |y| < %.1f",onia::rapForPTRange[rapBin-1],onia::rapForPTRange[rapBin]));
	top-=step;
	if(ptBin==0)
		latex->DrawLatex(left,top,Form("%.1f < p_{T} < %.1f GeV",onia::pTRange[rapBin][ptBin],onia::pTRange[rapBin][11]));
	else
		latex->DrawLatex(left,top,Form("%.1f < p_{T} < %.1f GeV",onia::pTRange[rapBin][ptBin-1],onia::pTRange[rapBin][ptBin]));
	if(cpmBin==0)
		latex->DrawLatex(left,top-step,Form("%.1f < N_{ch} < %.1f",onia::cpmRange[cpmBin],onia::cpmRange[onia::NchBins-1]));
	else
		latex->DrawLatex(left,top-step,Form("%.1f < N_{ch} < %.1f",onia::cpmRange[cpmBin-1],onia::cpmRange[cpmBin]));		

	top-=step;
	top-=0.5*step;
	latex->DrawLatex(left,top-step,Form("n_{CB}  =  %.1f",cbN));

	left=0.15; top=0.90; textSize=0.020;
	latex->SetTextSize(textSize);
	latex->DrawLatex(left,top,Form("#chi^{2}/ndof = %.2f / %d", chi2_Mass, ndof_Mass));
	top-=step;
	latex->DrawLatex(left,top,Form("mean   =  %.3f #pm %.3f MeV",Mean*1000, MeanErr*1000));
	top-=step;
	latex->DrawLatex(left,top,Form("#sigma_{1}  =  %.3f #pm %.3f MeV",Sigma*1000, SigmaErr*1000));
	top-=step;
	latex->DrawLatex(left,top,Form("#sigma_{2}  =  %.3f #pm %.3f MeV",Sigma2*1000, Sigma2Err*1000));
	top-=step;
	latex->DrawLatex(left,top,Form("#alpha  =  %.3f #pm %.3f",Alpha, AlphaErr));
	top-=step;
	latex->DrawLatex(left,top,Form("#lambda =  %.3f #pm %.3f",-lambda, lambdaErr)); // change to has positive value
	top-=step;
	latex->DrawLatex(left,top,Form("frac_{CB_{1}}  =  %.3f #pm %.3f",fracCB1, fracCB1Err));
	top-=step;
	latex->DrawLatex(left,top,Form("frac_{Bg}  =  %.3f #pm %.3f",fracBkg, fracBkgErr));
	top-=step;
	top-=step;
	top-=step;
	latex->DrawLatex(left,top,Form("effective #sigma =  %.3f MeV",SigmaWei*1000));
	top-=step;
	top-=0.5*step;
	latex->SetTextSize(textSize);
	latex->DrawLatex(left,top,Form("#frac{B}{B+S} (#pm3#sigma)  =  %.3f",BkgRatio3Sig));

	std::stringstream saveMass;
	saveMass << "Fit/mass_rap" << rapBin << "_pt" << ptBin << "_cpm" << cpmBin << ".pdf";
	c1->SaveAs(saveMass.str().c_str());

	delete c1;
	delete legendBlue;
	delete legendBlueDash;
	delete legendRed;
	delete legendBlack;
	delete legendGreen;
	delete legendGreenDash;
	delete legendPink;
	return;
}

//==============================================
void plotMassLog(RooWorkspace *ws, int rapBin, int ptBin, int cpmBin, int nState){
	int  nbins=90; //0.005 bin size
	TGaxis::SetMaxDigits(3);

	if(nState == 4) nbins=90;
	if(nState == 5) nbins=120;
	RooRealVar JpsiMass(*ws->var("JpsiMass"));

	RooPlot *massFrame=((RooRealVar*)ws->var("JpsiMass"))->frame(Bins(nbins));
	massFrame->SetName(Form("mass_plot_rap%d_pt%d_cpm%d",rapBin,ptBin,cpmBin));
	massFrame->SetTitle("");
	massFrame->GetYaxis()->SetTitle("Events / 5 MeV");
	massFrame->GetYaxis()->SetTitleOffset(1.3);

	RooPlot *massFramePull=((RooRealVar*)ws->var("JpsiMass"))->frame(Bins(nbins));
	massFramePull->SetName(Form("pullmass_plot_rap%d_pt%d_cpm%d",rapBin,ptBin,cpmBin));
	massFramePull->SetTitle("");
	massFramePull->GetYaxis()->SetTitle("pull");
	massFramePull->GetXaxis()->SetTitleSize(0.08);
	massFramePull->GetYaxis()->SetTitleSize(0.08);
	massFramePull->GetXaxis()->SetLabelSize(0.08);
	massFramePull->GetYaxis()->SetLabelSize(0.08);
	massFramePull->GetYaxis()->SetTitleOffset(0.4);
	massFramePull->GetYaxis()->SetRangeUser(-5.5,5.5);

	RooDataSet *data=(RooDataSet *)ws->data(Form("data_rap%d_pt%d_cpm%d",rapBin,ptBin,cpmBin));
	RooFitResult* fitRlt = (RooFitResult*)ws->obj(Form("m_fitresult_rap%d_pt%d_cpm%d",rapBin,ptBin,cpmBin));

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
		sqrt(pow((pow(Sigma,2)-pow(Sigma2,2))*fracCB1Err,2) + pow(2*fracCB1*Sigma*SigmaErr,2) + pow(2*(1-fracCB1)*Sigma2*Sigma2Err,2));
	cout<<">>=====Mean: "<<Mean<<endl;
	cout<<">>=====Sigma: "<<Sigma<<endl;
	cout<<">>=====Sigma2: "<<Sigma2<<endl;
	cout<<">>=====Alpha: "<<Alpha<<endl;
	cout<<">>=====cbN: "<<cbN<<endl;
	cout<<">>=====lambda: "<<lambda<<endl;
	cout<<">>=====fracCB1: "<<fracCB1<<endl;
	cout<<">>=====fracBkg: "<<fracBkg<<endl;
	cout<<">>=====SigmaWei: "<<SigmaWei<<endl;

	RooAddPdf *massPdf = (RooAddPdf*)ws->pdf("massModel");
	RooAddPdf *bkgMassShape = (RooAddPdf*)ws->pdf("bkgMassShape");

	double sigMaxMass = Mean+SigmaWei*onia::nSigMass;
	double sigMinMass = Mean-SigmaWei*onia::nSigMass;
	double sbHighMass = Mean+SigmaWei*onia::nSigBkgHigh;
	double sbLowMass =  Mean-SigmaWei*onia::nSigBkgLow;
	cout<<"sigMaxMass: "<<sigMaxMass<<endl;
	cout<<"sigMinMass: "<<sigMinMass<<endl;
	cout<<"sbHighMass: "<<sbHighMass<<endl;
	cout<<"sbLowMass: "<<sbLowMass<<endl;

	int nEntries = data->numEntries();
	JpsiMass.setRange("SigRegion",sigMinMass,sigMaxMass);

	RooRealVar* fracFull3Sig=(RooRealVar*)massPdf->createIntegral(JpsiMass,NormSet(JpsiMass),Range("SigRegion"));
	RooRealVar* fracBkg3Sig=(RooRealVar*)bkgMassShape->createIntegral(JpsiMass,NormSet(JpsiMass),Range("SigRegion"));
	double evtFull3Sig = nEntries*fracFull3Sig->getVal();
	double evtBkg3Sig = nEntries*fracBkg*fracBkg3Sig->getVal();
	double BkgRatio3Sig = evtBkg3Sig/evtFull3Sig;
	cout<<"fracFull3Sig: "<<fracFull3Sig->getVal()<<endl;
	cout<<"fracBkg3Sig: "<<fracBkg3Sig->getVal()<<endl;
	cout<<"evtFull3Sig: "<<evtFull3Sig<<endl;
	cout<<"evtBkg3Sig: "<<evtBkg3Sig<<endl;
	cout<<"BkgRatio3Sig: "<<BkgRatio3Sig<<endl;
	data->plotOn(massFrame,MarkerSize(0.8));
	massPdf->plotOn(massFrame,
			LineWidth(2),
			ProjWData(*data));
	//------get chi2------------SHOULD DONE after PLOTTING------
	int parsFit=(fitRlt->floatParsFinal()).getSize(); //this used the full p.d.f.
	int nBins_Mass=massFrame->GetNbinsX();
	double chi2Pre_Mass=massFrame->chiSquare(parsFit);  //reduced chi-squared = chi2/ndof
	int ndof_Mass=nBins_Mass-parsFit;  //num of degree of freedom
	double chi2_Mass=chi2Pre_Mass*ndof_Mass;

	TH1F* pull = new TH1F("pull","pull distribution", 100,-10.,10.);
	gSystem->mkdir("Fit/pull",kTRUE);
	TFile *pullFile = new TFile(Form("Fit/pull/pull_rap%d_pt%d_cpm%d_Mass.root",rapBin,ptBin,cpmBin),"RECREATE");

	RooHist* hpull_mass = massFrame->pullHist() ;
	hpull_mass->SetMarkerSize(0.8);
	for(int i=0;i<hpull_mass->GetN();i++){
		hpull_mass->SetPointEYlow(i,0.);
		hpull_mass->SetPointEYhigh(i,0.);
		double x,y;
		hpull_mass->GetPoint(i,x,y);
		pull->Fill(y);
	}
	pullFile->cd();
	pull->Write();
	pullFile->Close();

	massFramePull->addPlotable(hpull_mass,"P");

	massPdf->plotOn(massFrame,
			Components("bkgMassShape"),
			LineStyle(7),
			LineColor(kPink+3),
			LineWidth(2),
			ProjWData(*data));

	if(Sigma > Sigma2){

		fracCB1 = 1-fracCB1;
		double temp = 0.;
		temp = Sigma; Sigma = Sigma2; Sigma2 = temp;
		temp = SigmaErr; SigmaErr = Sigma2Err; Sigma2Err = temp;
		massPdf->plotOn(massFrame,
				Components("massCBShape"),
				LineStyle(2),
				LineColor(kGreen),
				LineWidth(2),
				ProjWData(*data));
		massPdf->plotOn(massFrame,
				Components("massCBShape2"),
				LineStyle(2),
				LineColor(kRed),
				LineWidth(2),
				ProjWData(*data));

	}  else{
		massPdf->plotOn(massFrame,
				Components("massCBShape"),
				LineStyle(2),
				LineColor(kRed),
				LineWidth(2),
				ProjWData(*data));
		massPdf->plotOn(massFrame,
				Components("massCBShape2"),
				LineStyle(2),
				LineColor(kGreen),
				LineWidth(2),
				ProjWData(*data));
	}
	double minY = 0.;
	double maxY = 0.;
	if(nState == 4) maxY = massFrame->GetMaximum()*0.1;
	if(nState == 5) maxY = massFrame->GetMaximum()*0.3;
	double lineWidth = 2.0;
	TLine *lineSBLow = new TLine(sbLowMass, minY, sbLowMass, maxY);
	TLine *lineSBHigh = new TLine(sbHighMass, minY, sbHighMass, maxY);
	TLine *lineSigLow = new TLine(sigMinMass, minY, sigMinMass, maxY);
	TLine *lineSigHigh = new TLine(sigMaxMass, minY, sigMaxMass, maxY);
	lineSBLow->SetLineWidth(lineWidth);lineSBHigh->SetLineWidth(lineWidth);
	lineSigLow->SetLineWidth(lineWidth);lineSigHigh->SetLineWidth(lineWidth);
	lineSBLow->SetLineColor(kBlue);lineSBHigh->SetLineColor(kBlue);
	lineSigLow->SetLineColor(kRed);lineSigHigh->SetLineColor(kRed);
	lineSBLow->SetLineStyle(7);lineSBHigh->SetLineStyle(7);
	lineSigLow->SetLineStyle(5);lineSigHigh->SetLineStyle(5);

	double Ymax = massFrame->GetMaximum();
	cout<<"Ymax: "<<Ymax<<endl;

	if(nState == 4) massFrame->SetMaximum(90000.);
	if(nState == 5) massFrame->SetMaximum(7000.);
	massFrame->SetMinimum(2.);
	TH1* legendBlue = data->createHistogram("legendBlue",JpsiMass,Binning(50)) ; legendBlue->SetLineColor(kBlue) ; legendBlue->SetLineStyle(kSolid) ; legendBlue->SetLineWidth(2) ;
	TH1* legendBlueDash = data->createHistogram("legendBlueDash",JpsiMass,Binning(50)) ; legendBlueDash->SetLineColor(kBlue) ; legendBlueDash->SetLineStyle(5) ; legendBlueDash->SetLineWidth(2) ;
	TH1* legendRed = data->createHistogram("legendRed",JpsiMass,Binning(50)) ; legendRed->SetLineColor(kRed) ; legendRed->SetLineStyle(2) ; legendRed->SetLineWidth(2) ;
	TH1* legendBlack = data->createHistogram("legendBlack",JpsiMass,Binning(50)) ; legendBlack->SetLineColor(kBlack) ; legendBlack->SetLineStyle(kSolid) ; legendBlack->SetLineWidth(2) ;
	TH1* legendGreen = data->createHistogram("legendGreen",JpsiMass,Binning(50)) ; legendGreen->SetLineColor(kGreen) ; legendGreen->SetLineStyle(2) ; legendGreen->SetLineWidth(2) ;
	TH1* legendPink = data->createHistogram("legendPink",JpsiMass,Binning(50)) ; legendPink->SetLineColor(kPink+3) ; legendPink->SetLineStyle(7) ; legendPink->SetLineWidth(2) ;

	TLegend* MassLegend=new TLegend(0.7,0.8,0.93,0.93);
	MassLegend->SetFillColor(kWhite);
	MassLegend->SetTextFont(42);
	MassLegend->SetTextSize(0.035);
	MassLegend->SetBorderSize(0.);
	MassLegend->AddEntry(legendBlue,"sum","l");
	MassLegend->AddEntry(legendRed,"signal CB_{1}","l");
	MassLegend->AddEntry(legendGreen,"signal CB_{2}","l");
	MassLegend->AddEntry(legendPink,"background","l");

	double left=0.15, top=0.90, textSize=0.025;
	TLatex *latex=new TLatex();
	latex->SetTextFont(42);
	latex->SetNDC(kTRUE);
	latex->SetTextSize(textSize);
	double step=textSize*1.3;

	gStyle->SetPadBottomMargin(0.08); //0.12
	gStyle->SetPadLeftMargin(0.09); //0.12
	gStyle->SetPadRightMargin(0.02); //0.05
	gStyle->SetPadTopMargin(0.02); //0.05

	TCanvas *c1=new TCanvas("c1","",800,700);

	c1->cd();

	c1->SetLogy(1);
	massFrame->Draw(); MassLegend->Draw();
	lineSBLow->Draw("same"); lineSBHigh->Draw("same"); lineSigLow->Draw("same"); lineSigHigh->Draw("same");
	textSize=0.030; latex->SetTextSize(textSize);
	if(nState == 4)
		latex->DrawLatex(left,top,"J/#psi");
	if(nState == 5)
		latex->DrawLatex(left,top,"#psi(2S)");
	top-=step;
	textSize=0.025; latex->SetTextSize(textSize);

	if(rapBin==0) latex->DrawLatex(left,top,Form("%.1f < |y| < %.1f",onia::rapForPTRange[rapBin],onia::rapForPTRange[2]));
	else if(rapBin==1) latex->DrawLatex(left,top,Form("|y| < %.1f",onia::rapForPTRange[rapBin]));
	else latex->DrawLatex(left,top,Form("%.1f < |y| < %.1f",onia::rapForPTRange[rapBin-1],onia::rapForPTRange[rapBin]));
	top-=step;
	if(ptBin==0)
		latex->DrawLatex(left,top,Form("%.1f < p_{T} < %.1f GeV",onia::pTRange[rapBin][ptBin],onia::pTRange[rapBin][11]));
	else
		latex->DrawLatex(left,top,Form("%.1f < p_{T} < %.1f GeV",onia::pTRange[rapBin][ptBin-1],onia::pTRange[rapBin][ptBin]));
	if(cpmBin==0)
		latex->DrawLatex(left,top-step,Form("%.1f < N_{ch} < %.1f",onia::cpmRange[cpmBin],onia::cpmRange[onia::NchBins-1]));
	else
		latex->DrawLatex(left,top-step,Form("%.1f < N_{ch} < %.1f",onia::cpmRange[cpmBin-1],onia::cpmRange[cpmBin]));			

	top-=step;
	top-=step;

	std::stringstream saveMasslog;
	saveMasslog << "Fit/mass_log_rap" << rapBin << "_pt" << ptBin << "_cpm" << cpmBin << ".pdf";
	c1->SaveAs(saveMasslog.str().c_str());

	delete c1;
	delete legendBlue;
	delete legendBlueDash;
	delete legendRed;
	delete legendBlack;
	delete legendGreen;
	delete legendPink;
	return;
}

//==============================================
void plotLifeBg(RooWorkspace *ws, int rapBin, int ptBin, int cpmBin, int nState){
	//int nbins=140; //bin size 0.025 mm
	double binSize = 0.025; // mm 
	double xaxisMin = -0.5, xaxisMax = 2.5;
	int nbins = (int)((xaxisMax-xaxisMin)/binSize);


	RooRealVar JpsiMass(*ws->var("JpsiMass"));
	RooRealVar *JpsictErr = ws->var("JpsictErr");

	RooPlot *ctauFrameBkgSBL=((RooRealVar*)ws->var("Jpsict"))->frame(Bins(nbins), Range(xaxisMin, xaxisMax));
	ctauFrameBkgSBL->SetName(Form("ctaubkg_plot_rap%d_pt%d_cpm%d",rapBin,ptBin,cpmBin));
	ctauFrameBkgSBL->SetYTitle("Events per 25 micron");
	ctauFrameBkgSBL->SetTitle("");
	ctauFrameBkgSBL->GetXaxis()->SetTitle("lifetime [mm]");

	RooPlot *ctauFrameBkgSBLPull=((RooRealVar*)ws->var("Jpsict"))->frame(Bins(nbins), Range(xaxisMin, xaxisMax));
	ctauFrameBkgSBLPull->SetName(Form("pullctaubkg_plot_rap%d_pt%d_cpm%d",rapBin,ptBin,cpmBin));
	ctauFrameBkgSBLPull->SetTitle("");
	ctauFrameBkgSBLPull->GetXaxis()->SetTitle("lifetime [mm]");
	ctauFrameBkgSBLPull->GetYaxis()->SetTitle("pull");
	ctauFrameBkgSBLPull->GetXaxis()->SetTitleSize(0.08);
	ctauFrameBkgSBLPull->GetYaxis()->SetTitleSize(0.08);
	ctauFrameBkgSBLPull->GetXaxis()->SetLabelSize(0.08);
	ctauFrameBkgSBLPull->GetYaxis()->SetLabelSize(0.08);
	ctauFrameBkgSBLPull->GetYaxis()->SetTitleOffset(0.4);
	ctauFrameBkgSBLPull->GetYaxis()->SetRangeUser(-5.5,5.5);

	RooPlot *ctauFrameBkgSBR=((RooRealVar*)ws->var("Jpsict"))->frame(Bins(nbins), Range(xaxisMin, xaxisMax));
	ctauFrameBkgSBR->SetName(Form("ctaubkg_plot_rap%d_pt%d_cpm%d",rapBin,ptBin,cpmBin));
	ctauFrameBkgSBR->SetYTitle("Events per 25 micron");
	ctauFrameBkgSBR->SetTitle("");
	ctauFrameBkgSBR->GetXaxis()->SetTitle("lifetime [mm]");

	RooPlot *ctauFrameBkgSBRPull=((RooRealVar*)ws->var("Jpsict"))->frame(Bins(nbins), Range(xaxisMin, xaxisMax));
	ctauFrameBkgSBRPull->SetName(Form("pullctaubkg_plot_rap%d_pt%d_cpm%d",rapBin,ptBin,cpmBin));
	ctauFrameBkgSBRPull->SetTitle("");
	ctauFrameBkgSBRPull->GetXaxis()->SetTitle("lifetime [mm]");
	ctauFrameBkgSBRPull->GetYaxis()->SetTitle("pull");
	ctauFrameBkgSBRPull->GetXaxis()->SetTitleSize(0.08);
	ctauFrameBkgSBRPull->GetYaxis()->SetTitleSize(0.08);
	ctauFrameBkgSBRPull->GetXaxis()->SetLabelSize(0.08);
	ctauFrameBkgSBRPull->GetYaxis()->SetLabelSize(0.08);
	ctauFrameBkgSBRPull->GetYaxis()->SetTitleOffset(0.4);
	ctauFrameBkgSBRPull->GetYaxis()->SetRangeUser(-5.5,5.5);

	RooDataSet *dataSBL=(RooDataSet *)ws->data(Form("data_rap%d_pt%d_cpm%d_SBL",rapBin,ptBin,cpmBin));
	RooDataSet *dataSBR=(RooDataSet *)ws->data(Form("data_rap%d_pt%d_cpm%d_SBR",rapBin,ptBin,cpmBin));

	RooFitResult* fitRlt = (RooFitResult*)ws->obj(Form("l_fitresult_rap%d_pt%d_cpm%d",rapBin,ptBin,cpmBin));
	if(!fitRlt) { cout<<">>=======Error: no Fit result object in workspace=========="<<endl; return; }
	fitRlt->Print();

	RooRealVar *promptMean_=(RooRealVar*)ws->var("promptMean");
	RooRealVar *ctResolution_=(RooRealVar*)ws->var("ctResolution");
	RooRealVar *ctResolution2_=(RooRealVar*)ws->var("ctResolution2");
	RooRealVar *fracGauss2_=(RooRealVar*)ws->var("fracGauss2");

	RooRealVar *bkgTauSSD_SBL_=(RooRealVar*)ws->var("bkgTauSSD_SBL");
	RooRealVar *bkgTauFD_=(RooRealVar*)ws->var("bkgTauFD");
	RooRealVar *bkgTauDSD_=(RooRealVar*)ws->var("bkgTauDSD");

	RooRealVar *bkgTauSSD_SBR_=(RooRealVar*)ws->var("bkgTauSSD_SBR");

	RooRealVar *fBkgSSDR_SBL_ = (RooRealVar*)ws->var("fBkgSSDR_SBL");
	RooRealVar *fBkgDSD_SBL_ = (RooRealVar*)ws->var("fBkgDSD_SBL");
	RooRealVar *fBkgSSDR_SBR_ = (RooRealVar*)ws->var("fBkgSSDR_SBR");
	RooRealVar *fBkgDSD_SBR_ = (RooRealVar*)ws->var("fBkgDSD_SBR");

	RooRealVar *fBkg_ = (RooRealVar*)ws->var("fBkg");

	RooRealVar *fBkgSBL_ = (RooRealVar*)ws->var("fBkgSBL");
	RooRealVar *fBkgSBR_ = (RooRealVar*)ws->var("fBkgSBR");

	RooRealVar *fPrompt_ = (RooRealVar*)ws->var("fPrompt");

	double promptMean = promptMean_->getVal();
	double promptMeanErr = promptMean_->getError();
	double promptCtRe = ctResolution_->getVal();
	double promptCtReErr = ctResolution_->getError();
	double promptCtRe2 = ctResolution2_->getVal();
	double promptCtRe2Err = ctResolution2_->getError();
	double fracGauss2 = fracGauss2_->getVal();
	double fracGauss2Err = fracGauss2_->getError();

	double bkgTauSSD_SBL = bkgTauSSD_SBL_->getVal();
	double bkgTauSSD_SBLErr = bkgTauSSD_SBL_->getError();
	double bkgTauFD = bkgTauFD_->getVal();
	double bkgTauFDErr = bkgTauFD_->getError();
	double bkgTauDSD = bkgTauDSD_->getVal();
	double bkgTauDSDErr = bkgTauDSD_->getError();

	double bkgTauSSD_SBR = bkgTauSSD_SBR_->getVal();
	double bkgTauSSD_SBRErr = bkgTauSSD_SBR_->getError();


	double fBkg = fBkg_->getVal();
	double fBkgErr = fBkg_->getError();

	double fBkgSBL = fBkgSBL_->getVal();
	double fBkgSBR = fBkgSBR_->getVal();
	double fBkgSBLErr = fBkgSBL_->getError();
	double fBkgSBRErr = fBkgSBR_->getError();

	double fPrompt = fPrompt_->getVal();
	double fPromptErr = fPrompt_->getError();
	double fNonPrompt =  1.-fBkg-fPrompt;
	double fNonPromptErr = sqrt(pow(fBkgErr,2)+pow(fPromptErr,2));

	//f_P / f_NP should be same in LSB,SR,RSB, to interpolate f_P in L(R)SB from SR
	double fPromptSBL = fPrompt * ( 1. - fBkgSBL ) / ( 1. - fBkg ) ;
	double fPromptSBR = fPrompt * ( 1. - fBkgSBR ) / ( 1. - fBkg ) ;
	double fPromptSBLErr = sqrt( 
			pow(fPromptErr*(1.-fBkgSBL )/(1.-fBkg),2) + 
			pow(fBkgSBLErr*fPrompt/(1.-fBkg),2) + 
			pow(fBkgErr*fPrompt*(1.-fBkgSBL )/pow(1.-fBkg,2),2) );
	double fPromptSBRErr = sqrt( 
			pow(fPromptErr*(1.-fBkgSBR )/(1.-fBkg),2) + 
			pow(fBkgSBRErr*fPrompt/(1.-fBkg),2) + 
			pow(fBkgErr*fPrompt*(1.-fBkgSBR )/pow(1.-fBkg,2),2) );

	double fNonPromptSBL = 1. - fBkgSBL - fPromptSBL ;
	double fNonPromptSBR = 1. - fBkgSBR - fPromptSBR ;
	double fNonPromptSBLErr = sqrt(
			pow(fPromptErr*(1.-fBkgSBL )/(1.-fBkg),2) +
			pow(fBkgSBLErr*(1.-fPrompt/(1.-fBkg)),2) +
			pow(fBkgErr*fPrompt*(1.-fBkgSBL )/pow(1.-fBkg,2),2) );
	double fNonPromptSBRErr = sqrt(
			pow(fPromptErr*(1.-fBkgSBR )/(1.-fBkg),2) +
			pow(fBkgSBRErr*(1.-fPrompt/(1.-fBkg)),2) +
			pow(fBkgErr*fPrompt*(1.-fBkgSBR )/pow(1.-fBkg,2),2) );

	cout<<"fBkgSBL: "<<fBkgSBL<<" fBkgSBR: "<<fBkgSBR<<endl;
	cout<<"fPromptSBL: "<<fPromptSBL<<" fPromptSBR: "<<fPromptSBR <<endl;
	cout<<"fNonPromptSBL: "<<fNonPromptSBL<<" fNonPromptSBR: "<<fNonPromptSBL<<endl;
	cout<<"fPromptSBL/fNonPromptSBL: "<<fPromptSBL/fNonPromptSBL<<endl;
	cout<<"fPromptSBR/fNonPromptSBR: "<<fPromptSBR/fNonPromptSBR<<endl;
	cout<<"fPrompt/fNonPrompt: "<<fPrompt/fNonPrompt<<endl;

	double fBkgSSDR_SBL = fBkgSSDR_SBL_->getVal();
	double fBkgSSDRErr_SBL = fBkgSSDR_SBL_->getError();
	double fBkgDSD_SBL = fBkgDSD_SBL_->getVal();
	double fBkgDSDErr_SBL = fBkgDSD_SBL_->getError();
	double fBkgSSDL_SBL = 1-fBkgSSDR_SBL-fBkgDSD_SBL;

	double fBkgSSDR_SBR = fBkgSSDR_SBR_->getVal();
	double fBkgSSDRErr_SBR = fBkgSSDR_SBR_->getError();
	double fBkgDSD_SBR = fBkgDSD_SBR_->getVal();
	double fBkgDSDErr_SBR = fBkgDSD_SBR_->getError();
	double fBkgSSDL_SBR = 1-fBkgSSDR_SBR-fBkgDSD_SBR;

	cout<<"promptMean: "<<promptMean<<endl;
	cout<<"promptCtRe: "<<promptCtRe<<endl;
	cout<<"fBkgSSDR_SBL: "<<fBkgSSDR_SBL<<endl;
	cout<<"fBkgSSDL_SBL: "<<fBkgSSDL_SBL<<endl;
	cout<<"fBkgSSDR_SBR: "<<fBkgSSDR_SBR<<endl;
	cout<<"fBkgSSDL_SBR: "<<fBkgSSDL_SBR<<endl;

	RooAddPdf *ModelLifeSBL = (RooAddPdf*)ws->pdf("backgroundlifetimeL");
	RooAddPdf *ModelLifeSBR = (RooAddPdf*)ws->pdf("backgroundlifetimeR");
	RooAddPdf *backgroundSSD_SBL = (RooAddPdf*)ws->pdf("backgroundSSD_SBL");
	RooAddPdf *backgroundSSD_SBR = (RooAddPdf*)ws->pdf("backgroundSSD_SBR");
	RooAddPdf *backgroundDSD = (RooAddPdf*)ws->pdf("backgroundDSD");
	RooAddPdf *backgroundFD = (RooAddPdf*)ws->pdf("backgroundFD");

	RooAddModel *Prompt = (RooAddModel*)ws->pdf("TotalPromptLifetime");
	RooAbsPdf *nonPromptSSD = (RooAbsPdf*)ws->pdf("nonPromptSSD");


	//relative fraction of Prompt: f_P/f_NP, same in LSB,SR,RSB
	RooRealVar *fracPrompt = new RooRealVar("fracPrompt","fracPrompt",fPromptSBL/(fPromptSBL+fNonPromptSBL));
	RooAddPdf *SignalPdf = new RooAddPdf("SignalPdf","SignalPdf",RooArgList(*Prompt,*nonPromptSSD),*fracPrompt);


	int parsFit;

	//ploting background SBL
	dataSBL->plotOn(ctauFrameBkgSBL,MarkerSize(0.8));
	ModelLifeSBL->plotOn(ctauFrameBkgSBL,
			ProjWData(*JpsictErr, *dataSBL),
			LineWidth(2),
			NumCPU(1));

	parsFit=(fitRlt->floatParsFinal()).getSize(); //this used the full p.d.f.
	int nBins_LSBL=ctauFrameBkgSBL->GetNbinsX();
	double chi2Pre_LSBL=ctauFrameBkgSBL->chiSquare(parsFit);  //reduced chi-squared = chi2/ndof
	int ndof_LSBL=nBins_LSBL-parsFit;  //num of degree of freedom
	double chi2_LSBL=chi2Pre_LSBL*ndof_LSBL;

	TH1F* pull = new TH1F("pull","pull distribution", 100,-10.,10.);
	gSystem->mkdir("Fit/pull",kTRUE);
	TFile *pullFile_LSB = new TFile(Form("Fit/pull/pull_rap%d_pt%d_cpm%d_LSB.root",rapBin,ptBin,cpmBin),"RECREATE");

	RooHist* hpull_ctauBkgSBL = ctauFrameBkgSBL->pullHist() ;
	hpull_ctauBkgSBL->SetMarkerSize(0.8);
	for(int i=0;i<hpull_ctauBkgSBL->GetN();i++){
		hpull_ctauBkgSBL->SetPointEYlow(i,0.);
		hpull_ctauBkgSBL->SetPointEYhigh(i,0.);
		double x,y;
		hpull_ctauBkgSBL->GetPoint(i,x,y);
		pull->Fill(y);
	}
	cout<<"hpull_ctauBkgSBL->GetN(): "<<hpull_ctauBkgSBL->GetN()<<endl;
	pullFile_LSB->cd();
	pull->Write();
	pullFile_LSB->Close();

	ctauFrameBkgSBLPull->addPlotable(hpull_ctauBkgSBL,"P");

	backgroundSSD_SBL->plotOn(ctauFrameBkgSBL,
			ProjWData(*JpsictErr, *dataSBL),
			Normalization(fBkgSBL*fBkgSSDR_SBL),
			LineStyle(2),
			LineColor(kRed),
			LineWidth(2),
			NumCPU(1));

	backgroundFD->plotOn(ctauFrameBkgSBL,
			ProjWData(*JpsictErr, *dataSBL),
			Normalization(fBkgSBL*(1.-fBkgDSD_SBL-fBkgSSDR_SBL)),
			LineStyle(2),
			LineColor(kGreen),
			LineWidth(2),
			NumCPU(1));

	backgroundDSD->plotOn(ctauFrameBkgSBL,
			ProjWData(*JpsictErr, *dataSBL),
			Normalization(fBkgSBL*fBkgDSD_SBL),
			LineStyle(1),
			LineColor(kBlack),
			LineWidth(2),
			NumCPU(1));

	SignalPdf->plotOn(ctauFrameBkgSBL,
			ProjWData(*JpsictErr, *dataSBL),
			Normalization(1.-fBkgSBL),
			LineStyle(7),
			LineColor(kPink+3),
			LineWidth(2));

	double Ymax = ctauFrameBkgSBL->GetMaximum();
	cout<<"Ymax: "<<Ymax<<endl;

	ctauFrameBkgSBL->SetMaximum(5*Ymax);
	ctauFrameBkgSBL->SetMinimum(1.1);

	dataSBR->plotOn(ctauFrameBkgSBR,MarkerSize(0.8));
	ModelLifeSBR->plotOn(ctauFrameBkgSBR,
			ProjWData(*JpsictErr, *dataSBR),
			LineWidth(2),
			NumCPU(1));

	parsFit=(fitRlt->floatParsFinal()).getSize(); //this used the full p.d.f.
	int nBins_LSBR=ctauFrameBkgSBR->GetNbinsX();
	double chi2Pre_LSBR=ctauFrameBkgSBR->chiSquare(parsFit);  //reduced chi-squared = chi2/ndof
	int ndof_LSBR=nBins_LSBR-parsFit;  //num of degree of freedom
	double chi2_LSBR=chi2Pre_LSBR*ndof_LSBR;

	pull = new TH1F("pull","pull distribution", 100,-10.,10.);
	gSystem->mkdir("Fit/pull",kTRUE);
	TFile *pullFile_RSB = new TFile(Form("Fit/pull/pull_rap%d_pt%d_cpm%d_RSB.root",rapBin,ptBin,cpmBin),"RECREATE");

	RooHist* hpull_ctauBkgSBR = ctauFrameBkgSBR->pullHist() ;
	hpull_ctauBkgSBR->SetMarkerSize(0.8);
	for(int i=0;i<hpull_ctauBkgSBR->GetN();i++){
		hpull_ctauBkgSBR->SetPointEYlow(i,0.);
		hpull_ctauBkgSBR->SetPointEYhigh(i,0.);
		double x,y;
		hpull_ctauBkgSBR->GetPoint(i,x,y);
		pull->Fill(y);
	}

	pullFile_RSB->cd();
	pull->Write();
	pullFile_RSB->Close();

	ctauFrameBkgSBRPull->addPlotable(hpull_ctauBkgSBR,"P");

	backgroundSSD_SBR->plotOn(ctauFrameBkgSBR,
			ProjWData(*JpsictErr, *dataSBR),
			Normalization(fBkgSBR*fBkgSSDR_SBR),
			LineStyle(2),
			LineColor(kRed),
			LineWidth(2),
			NumCPU(1));

	backgroundFD->plotOn(ctauFrameBkgSBR,
			ProjWData(*JpsictErr, *dataSBR),
			Normalization(fBkgSBR*(1.-fBkgSSDR_SBR-fBkgDSD_SBR)),
			LineStyle(2),
			LineColor(kGreen),
			LineWidth(2),
			NumCPU(1));

	backgroundDSD->plotOn(ctauFrameBkgSBR,
			ProjWData(*JpsictErr, *dataSBR),
			Normalization(fBkgSBR*fBkgDSD_SBR),
			LineStyle(1),
			LineColor(kBlack),
			LineWidth(2),
			NumCPU(1));

	SignalPdf->plotOn(ctauFrameBkgSBR,
			ProjWData(*JpsictErr, *dataSBR),
			Normalization(1.-fBkgSBR),
			LineStyle(7),
			LineColor(kPink+3),
			LineWidth(2));

	Ymax = ctauFrameBkgSBR->GetMaximum();
	cout<<"Ymax: "<<Ymax<<endl;
	ctauFrameBkgSBR->SetMaximum(5*Ymax);
	ctauFrameBkgSBR->SetMinimum(1.1);

	TH1* legendBlue = dataSBL->createHistogram("legendBlue",JpsiMass,Binning(50)) ; legendBlue->SetLineColor(kBlue) ; legendBlue->SetLineStyle(kSolid) ; legendBlue->SetLineWidth(2) ;
	TH1* legendBlueDash = dataSBL->createHistogram("legendBlueDash",JpsiMass,Binning(50)) ; legendBlueDash->SetLineColor(kBlue) ; legendBlueDash->SetLineStyle(5) ; legendBlueDash->SetLineWidth(2) ;
	TH1* legendRed = dataSBL->createHistogram("legendRed",JpsiMass,Binning(50)) ; legendRed->SetLineColor(kRed) ; legendRed->SetLineStyle(2) ; legendRed->SetLineWidth(2) ;
	TH1* legendBlack = dataSBL->createHistogram("legendBlack",JpsiMass,Binning(50)) ; legendBlack->SetLineColor(kBlack) ; legendBlack->SetLineStyle(kSolid) ; legendBlack->SetLineWidth(2) ;
	TH1* legendGreen = dataSBL->createHistogram("legendGreen",JpsiMass,Binning(50)) ; legendGreen->SetLineColor(kGreen) ; legendGreen->SetLineStyle(kSolid) ; legendGreen->SetLineWidth(2) ;
	TH1* legendGreenDash = dataSBL->createHistogram("legendGreenDash",JpsiMass,Binning(50)) ; legendGreenDash->SetLineColor(kGreen) ; legendGreenDash->SetLineStyle(2) ; legendGreenDash->SetLineWidth(2) ;
	TH1* legendPink = dataSBL->createHistogram("legendPink",JpsiMass,Binning(50)) ; legendPink->SetLineColor(kPink+3) ; legendPink->SetLineStyle(7) ; legendPink->SetLineWidth(2) ;

	TLegend* LifetimeLegendBkgSBL=new TLegend(0.35,0.7,0.5,0.93);
	LifetimeLegendBkgSBL->SetFillColor(kWhite);
	LifetimeLegendBkgSBL->SetTextFont(42);
	LifetimeLegendBkgSBL->SetTextSize(0.035);
	LifetimeLegendBkgSBL->SetBorderSize(0.);
	LifetimeLegendBkgSBL->AddEntry(legendBlue,"sum","l");
	LifetimeLegendBkgSBL->AddEntry(legendGreenDash,"SSL","l");
	LifetimeLegendBkgSBL->AddEntry(legendRed,"SSR ","l");
	LifetimeLegendBkgSBL->AddEntry(legendBlack,"DS","l");
	LifetimeLegendBkgSBL->AddEntry(legendPink,"signal","l");

	TLegend* LifetimeLegendBkgSBR=new TLegend(0.35,0.7,0.5,0.93);
	LifetimeLegendBkgSBR->SetFillColor(kWhite);
	LifetimeLegendBkgSBR->SetTextFont(42);
	LifetimeLegendBkgSBR->SetTextSize(0.035);
	LifetimeLegendBkgSBR->SetBorderSize(0.);
	LifetimeLegendBkgSBR->AddEntry(legendBlue,"sum","l");
	LifetimeLegendBkgSBR->AddEntry(legendGreenDash,"SSL","l");
	LifetimeLegendBkgSBR->AddEntry(legendRed,"SSR ","l");
	LifetimeLegendBkgSBR->AddEntry(legendBlack,"DS","l");
	LifetimeLegendBkgSBR->AddEntry(legendPink,"signal","l");

	double left=0.7, top=0.9, textSize=0.025;
	TLatex *latex=new TLatex();
	latex->SetTextFont(42);
	latex->SetNDC(kTRUE);
	latex->SetTextSize(textSize);
	double step=textSize*1.3;

	gStyle->SetPadBottomMargin(0.11); //0.12
	gStyle->SetPadLeftMargin(0.08); //0.12
	gStyle->SetPadRightMargin(0.02); //0.05
	gStyle->SetPadTopMargin(0.02); //0.05

	TCanvas *c1=new TCanvas("c1","",800,900);

	c1->cd();
	TPad *pad1 = new TPad("pad1","pad1",0.,0.,1.,0.3);
	pad1->SetGridy();
	pad1->SetBottomMargin(0.2);
	pad1->Draw();
	c1->cd();
	TPad *pad2 = new TPad("pad2","pad2",0.,0.3,1.,1.);
	pad2->Draw();

	//Bkg SBL
	pad2->cd(0);pad2->SetLogy(1);
	ctauFrameBkgSBL->Draw();
	LifetimeLegendBkgSBL->Draw();
	left=0.5; top=0.9; textSize=0.030; latex->SetTextSize(textSize);
	if(nState == 4) latex->DrawLatex(left,top,"J/#psi");
	if(nState == 5) latex->DrawLatex(left,top,"#psi(2S)");
	top-=step;
	textSize=0.025; latex->SetTextSize(textSize);

	if(rapBin==0) latex->DrawLatex(left,top,Form("%.1f < |y| < %.1f",onia::rapForPTRange[rapBin],onia::rapForPTRange[2]));
	else if(rapBin==1) latex->DrawLatex(left,top,Form("|y| < %.1f",onia::rapForPTRange[rapBin]));
	else latex->DrawLatex(left,top,Form("%.1f < |y| < %.1f",onia::rapForPTRange[rapBin-1],onia::rapForPTRange[rapBin]));
	top-=step;
	if(ptBin==0)
		latex->DrawLatex(left,top,Form("%.1f < p_{T} < %.1f GeV",onia::pTRange[rapBin][ptBin],onia::pTRange[rapBin][11]));
	else
		latex->DrawLatex(left,top,Form("%.1f < p_{T} < %.1f GeV",onia::pTRange[rapBin][ptBin-1],onia::pTRange[rapBin][ptBin]));
	if(cpmBin==0)
		latex->DrawLatex(left,top-step,Form("%.1f < N_{ch} < %.1f",onia::cpmRange[cpmBin],onia::cpmRange[onia::NchBins-1]));
	else
		latex->DrawLatex(left,top-step,Form("%.1f < N_{ch} < %.1f",onia::cpmRange[cpmBin-1],onia::cpmRange[cpmBin]));		

	top-=step;
	latex->SetTextColor(kRed);
	latex->DrawLatex(left,top-step,"Left Mass SideBand");
	latex->SetTextColor(kBlack);
	top-=step;
	latex->DrawLatex(left,top-step,Form("mean  =  %.1f",promptMean));
	top-=step;
	latex->DrawLatex(left,top-step,Form("resolution1 =  %.1f",promptCtRe));
	top-=step;
	latex->DrawLatex(left,top-step,Form("resolution2  =  %.1f",promptCtRe2));
	top-=step;
	if(ptBin>7)
		latex->DrawLatex(left,top,Form("frac_{SSR}   =  %.3f",fBkgSSDR_SBL));


	left=0.7; top=0.9;
	latex->DrawLatex(left,top,Form("#chi^{2}/ndof   =  %.2f/%d",chi2_LSBL,ndof_LSBL));
	top-=step;
	latex->DrawLatex(left,top,Form("frac_{DS}   =  %.3f #pm %.3f",fBkgDSD_SBL, fBkgDSDErr_SBL));
	top-=step;
	if(ptBin<8){
		latex->DrawLatex(left,top,Form("frac_{SSR}   =  %.3f #pm %.3f",fBkgSSDR_SBL, fBkgSSDRErr_SBL));
		top-=step;
	}
	latex->DrawLatex(left,top,Form("#tau_{SSR}   =  %.3f #pm %.3f",bkgTauSSD_SBL, bkgTauSSD_SBLErr));
	top-=step;
	latex->DrawLatex(left,top,Form("frac_{Bg}   =  %.3f #pm %.3f",fBkgSBL,fBkgSBLErr));
	top-=step;
	top-=step;

	latex->DrawLatex(left,top,Form("frac_{Gauss2}   =  %.3f #pm %.3f",fracGauss2, fracGauss2Err));
	top-=step;
	latex->DrawLatex(left,top,Form("#tau_{DS}   =  %.3f #pm %.3f",bkgTauDSD, bkgTauDSDErr));
	top-=step;
	latex->DrawLatex(left,top,Form("#tau_{SSL}  =  %.3f #pm %.3f",bkgTauFD, bkgTauFDErr));

	pad1->cd(0); pad1->SetLogy(0);
	ctauFrameBkgSBLPull->Draw();

	std::stringstream savectLSB;
	savectLSB << "Fit/ct_LSB_rap" << rapBin << "_pt" << ptBin << "_cpm" <<cpmBin << "_withPull.pdf";
	c1->SaveAs(savectLSB.str().c_str());

	//Bkg SBR
	pad2->cd(0);pad2->SetLogy(1);
	ctauFrameBkgSBR->Draw();
	LifetimeLegendBkgSBR->Draw();
	left=0.5; top=0.9; textSize=0.030; latex->SetTextSize(textSize);
	if(nState == 4) latex->DrawLatex(left,top,"J/#psi");
	if(nState == 5) latex->DrawLatex(left,top,"#psi(2S)");
	top-=step;
	textSize=0.025; latex->SetTextSize(textSize);

	if(rapBin==0) latex->DrawLatex(left,top,Form("%.1f < |y| < %.1f",onia::rapForPTRange[rapBin],onia::rapForPTRange[2]));
	else if(rapBin==1) latex->DrawLatex(left,top,Form("|y| < %.1f",onia::rapForPTRange[rapBin]));
	else latex->DrawLatex(left,top,Form("%.1f < |y| < %.1f",onia::rapForPTRange[rapBin-1],onia::rapForPTRange[rapBin]));
	top-=step;
	if(ptBin==0)
		latex->DrawLatex(left,top,Form("%.1f < p_{T} < %.1f GeV",onia::pTRange[rapBin][ptBin],onia::pTRange[rapBin][11]));
	else
		latex->DrawLatex(left,top,Form("%.1f < p_{T} < %.1f GeV",onia::pTRange[rapBin][ptBin-1],onia::pTRange[rapBin][ptBin]));
	if(cpmBin==0)
		latex->DrawLatex(left,top-step,Form("%.1f < N_{ch} < %.1f",onia::cpmRange[cpmBin],onia::cpmRange[onia::NchBins-1]));
	else
		latex->DrawLatex(left,top-step,Form("%.1f < N_{ch} < %.1f",onia::cpmRange[cpmBin-1],onia::cpmRange[cpmBin]));			

	top-=step;
	latex->SetTextColor(kRed);
	latex->DrawLatex(left,top-step,"Right Mass SideBand");
	latex->SetTextColor(kBlack);
	top-=step;

	left=0.7; top=0.9;
	latex->DrawLatex(left,top,Form("#chi^{2}/ndof   =  %.2f/%d",chi2_LSBR,ndof_LSBR));
	top-=step;
	latex->DrawLatex(left,top,Form("frac_{DS}   =  %.3f #pm %.3f",fBkgDSD_SBR, fBkgDSDErr_SBR));
	top-=step;
	if(ptBin<8){
		latex->DrawLatex(left,top,Form("frac_{SSR}   =  %.3f #pm %.3f",fBkgSSDR_SBR, fBkgSSDRErr_SBR));
		top-=step;
	}
	latex->DrawLatex(left,top,Form("#tau_{SSR}   =  %.3f #pm %.3f",bkgTauSSD_SBR, bkgTauSSD_SBRErr));
	top-=step;
	latex->DrawLatex(left,top,Form("frac_{Bg}   =  %.3f #pm %.3f",fBkgSBR,fBkgSBRErr));
	top-=step;
	top-=step;

	pad1->cd(0); pad1->SetLogy(0);
	ctauFrameBkgSBRPull->Draw();

	std::stringstream savectRSB;
	savectRSB << "Fit/ct_RSB_rap" << rapBin << "_pt" << ptBin << "_cpm" << cpmBin << "_withPull.pdf";
	c1->SaveAs(savectRSB.str().c_str());

	delete c1;
	delete legendBlue;
	delete legendBlueDash;
	delete legendRed;
	delete legendBlack;
	delete legendGreen;
	delete legendGreenDash;
	delete legendPink;
	return;
}

//==============================================
void plotLifeSig(RooWorkspace *ws, int rapBin, int ptBin, int cpmBin, int nState){
	//int nbins=230; //bin size 0.01 mm
	double binSize = 0.01; // mm 
	double xaxisMin = -0.5, xaxisMax = 2.5;
	int nbins = (int)((xaxisMax-xaxisMin)/binSize);


	RooRealVar JpsiMass(*ws->var("JpsiMass"));
	RooRealVar Jpsict(*ws->var("Jpsict"));
	RooRealVar *JpsictErr = ws->var("JpsictErr");

	RooPlot *ctauFrameSig=((RooRealVar*)ws->var("Jpsict"))->frame(Bins(nbins), Range(xaxisMin, xaxisMax));
	ctauFrameSig->SetName(Form("ctausig_plot_rap%d_pt%d_cpm%d",rapBin,ptBin,cpmBin));
	ctauFrameSig->SetYTitle("Events per 10 micron");
	ctauFrameSig->SetTitle("");
	ctauFrameSig->SetXTitle("lifetime [mm]");

	RooPlot *ctauFrameSigPull=((RooRealVar*)ws->var("Jpsict"))->frame(Bins(nbins), Range(xaxisMin, xaxisMax));
	ctauFrameSigPull->SetName(Form("pullctausig_plot_rap%d_pt%d_cpm%d",rapBin,ptBin,cpmBin));
	ctauFrameSigPull->SetTitle("");
	ctauFrameSigPull->SetXTitle("lifetime [mm]");
	ctauFrameSigPull->GetYaxis()->SetTitle("pull");
	ctauFrameSigPull->GetXaxis()->SetTitleSize(0.08);
	ctauFrameSigPull->GetYaxis()->SetTitleSize(0.08);
	ctauFrameSigPull->GetXaxis()->SetLabelSize(0.08);
	ctauFrameSigPull->GetYaxis()->SetLabelSize(0.08);
	ctauFrameSigPull->GetYaxis()->SetTitleOffset(0.4);
	ctauFrameSigPull->GetYaxis()->SetRangeUser(-5.5,5.5);


	RooDataSet *dataSR=(RooDataSet *)ws->data(Form("data_rap%d_pt%d_cpm%d_SR",rapBin,ptBin,cpmBin));

	RooFitResult* fitRlt = (RooFitResult*)ws->obj(Form("l_fitresult_rap%d_pt%d_cpm%d",rapBin,ptBin,cpmBin));
	if(!fitRlt) { cout<<">>=======Error: no Fit result object in workspace=========="<<endl; return; }
	fitRlt->Print();
	cout<<"dataSR->numEntries: "<<dataSR->numEntries()<<endl;

	RooRealVar *promptMean_=(RooRealVar*)ws->var("promptMean");
	RooRealVar *ctResolution_=(RooRealVar*)ws->var("ctResolution");
	RooRealVar *ctResolution2_=(RooRealVar*)ws->var("ctResolution2");
	RooRealVar *fracGauss2_=(RooRealVar*)ws->var("fracGauss2");

	RooRealVar *nonPromptTau_=(RooRealVar*)ws->var("nonPromptTau");
	RooRealVar *fBkg_ = (RooRealVar*)ws->var("fBkg");
	RooRealVar *fPrompt_ = (RooRealVar*)ws->var("fPrompt");

	double promptMean = promptMean_->getVal();
	double promptMeanErr = promptMean_->getError();
	double promptCtRe = ctResolution_->getVal();
	double promptCtReErr = ctResolution_->getError();
	double promptCtRe2 = ctResolution2_->getVal();
	double promptCtRe2Err = ctResolution2_->getError();
	double fracGauss2 = fracGauss2_->getVal();
	double fracGauss2Err = fracGauss2_->getError();
	double nonPromptTau = nonPromptTau_->getVal();
	double nonPromptTauErr = nonPromptTau_->getError();

	double fBkg = fBkg_->getVal();
	double fBkgErr = fBkg_->getError();
	double fPrompt = fPrompt_->getVal();
	double fPromptErr = fPrompt_->getError();
	double fNonPrompt =  1.-fBkg-fPrompt;

	cout<<"promptMean: "<<promptMean<<endl;
	cout<<"promptCtRe: "<<promptCtRe<<endl;
	cout<<"nonPromptTau: "<<nonPromptTau<<endl;
	cout<<"fBkg: "<<fBkg<<endl;
	cout<<"fBkgErr: "<<fBkgErr<<endl;
	cout<<"fPrompt: "<<fPrompt<<endl;
	cout<<"fNonPrompt: "<<fNonPrompt<<endl;

	RooAddPdf *ModelLife = (RooAddPdf*)ws->pdf("fulllifetime");

	RooAddModel *Prompt = (RooAddModel*)ws->pdf("TotalPromptLifetime");
	RooAbsPdf *nonPromptSSD = (RooAbsPdf*)ws->pdf("nonPromptSSD");
	RooAbsPdf *backgroundlifetime = (RooAbsPdf*)ws->pdf("backgroundlifetime");

	int parsFit;
	//ploting signal region
	dataSR->plotOn(ctauFrameSig,MarkerSize(0.8), Name("myHist"));
	ModelLife->plotOn(ctauFrameSig,
			ProjWData(*JpsictErr, *dataSR),
			LineWidth(2), NumCPU(1),
			Name("myCurve"));

	//------get chi2------------SHOULD DONE after PLOTTING------
	parsFit=(fitRlt->floatParsFinal()).getSize(); //this used the full p.d.f.

	int nBins_LSig=ctauFrameSig->GetNbinsX();
	double chi2Pre_LSig=ctauFrameSig->chiSquare(parsFit);  //reduced chi-squared = chi2/ndof
	int ndof_LSig=nBins_LSig-parsFit;  //num of degree of freedom
	double chi2_LSig=chi2Pre_LSig*ndof_LSig;

	RooHist* hpull_ctauSig = ctauFrameSig->pullHist() ;
	//RooHist* hpull_ctauSig = ctauFrameSig->pullHist("myHist","myCurve",kTRUE);

	//save pull distribution
	TH1F* pull = new TH1F("pull","pull distribution", 100,-10.,10.);
	gSystem->mkdir("Fit/pull",kTRUE);
	//TFile *pullFile = new TFile(Form("Fit/pull/pullAverage_rap%d_pt%d_cpm%d_SR.root",rapBin,ptBin,cpmBin),"RECREATE");
	TFile *pullFile = new TFile(Form("Fit/pull/pull_rap%d_pt%d_cpm%d_SR.root",rapBin,ptBin,cpmBin),"RECREATE");

	hpull_ctauSig->SetMarkerSize(0.8);
	for(int i=0;i<hpull_ctauSig->GetN();i++){
		hpull_ctauSig->SetPointEYlow(i,0.);
		hpull_ctauSig->SetPointEYhigh(i,0.);
		double x,y;
		hpull_ctauSig->GetPoint(i,x,y);
		pull->Fill(y);
	}
	pullFile->cd();
	pull->Write();
	pullFile->Close();

	ctauFrameSigPull->addPlotable(hpull_ctauSig,"P");

	Prompt->plotOn(ctauFrameSig,
			ProjWData(*JpsictErr, *dataSR),
			Normalization(fPrompt),
			LineStyle(5),
			LineColor(kBlue),
			LineWidth(2), NumCPU(1));

	nonPromptSSD->plotOn(ctauFrameSig,
			ProjWData(*JpsictErr, *dataSR),
			Normalization(fNonPrompt),
			LineStyle(2),
			LineColor(kRed),
			LineWidth(2), NumCPU(1));

	backgroundlifetime->plotOn(ctauFrameSig,
			ProjWData(*JpsictErr, *dataSR),
			Normalization(fBkg),
			LineStyle(7),
			LineColor(kPink+3),
			LineWidth(2), NumCPU(1));

	double Ymax = ctauFrameSig->GetMaximum();
	cout<<"Ymax: "<<Ymax<<endl;
	ctauFrameSig->SetMaximum(3*Ymax);
	ctauFrameSig->SetMinimum(1.1);


	TH1* legendBlue = dataSR->createHistogram("legendBlue",JpsiMass,Binning(50)) ; legendBlue->SetLineColor(kBlue) ; legendBlue->SetLineStyle(kSolid) ; legendBlue->SetLineWidth(2) ;
	TH1* legendBlueDash = dataSR->createHistogram("legendBlueDash",JpsiMass,Binning(50)) ; legendBlueDash->SetLineColor(kBlue) ; legendBlueDash->SetLineStyle(5) ; legendBlueDash->SetLineWidth(2) ;
	TH1* legendRed = dataSR->createHistogram("legendRed",JpsiMass,Binning(50)) ; legendRed->SetLineColor(kRed) ; legendRed->SetLineStyle(2) ; legendRed->SetLineWidth(2) ;
	TH1* legendBlack = dataSR->createHistogram("legendBlack",JpsiMass,Binning(50)) ; legendBlack->SetLineColor(kBlack) ; legendBlack->SetLineStyle(kSolid) ; legendBlack->SetLineWidth(2) ;
	TH1* legendGreen = dataSR->createHistogram("legendGreen",JpsiMass,Binning(50)) ; legendGreen->SetLineColor(kGreen) ; legendGreen->SetLineStyle(kSolid) ; legendGreen->SetLineWidth(2) ;
	TH1* legendGreenDash = dataSR->createHistogram("legendGreenDash",JpsiMass,Binning(50)) ; legendGreenDash->SetLineColor(kGreen) ; legendGreenDash->SetLineStyle(2) ; legendGreenDash->SetLineWidth(2) ;
	TH1* legendPink = dataSR->createHistogram("legendPink",JpsiMass,Binning(50)) ; legendPink->SetLineColor(kPink+3) ; legendPink->SetLineStyle(7) ; legendPink->SetLineWidth(2) ;

	TLegend* LifetimeLegendSig=new TLegend(0.31,0.8,0.43,0.96);
	LifetimeLegendSig->SetFillColor(kWhite);
	LifetimeLegendSig->SetTextFont(42);
	LifetimeLegendSig->SetTextSize(0.035);
	LifetimeLegendSig->SetBorderSize(0.);
	LifetimeLegendSig->AddEntry(legendBlue,"sum","l");
	LifetimeLegendSig->AddEntry(legendBlueDash,"prompt","l");
	LifetimeLegendSig->AddEntry(legendRed,"non-prompt","l");
	LifetimeLegendSig->AddEntry(legendPink,"background","l");

	double left=0.5, top=0.9, textSize=0.025;
	TLatex *latex=new TLatex();
	latex->SetTextFont(42);
	latex->SetNDC(kTRUE);
	latex->SetTextSize(textSize);
	double step=textSize*1.3;

	gStyle->SetPadBottomMargin(0.11); //0.12
	gStyle->SetPadLeftMargin(0.08); //0.12
	gStyle->SetPadRightMargin(0.02); //0.05
	gStyle->SetPadTopMargin(0.02); //0.05

	TCanvas *c1=new TCanvas("c1","",800,900);

	c1->cd();
	TPad *pad1 = new TPad("pad1","pad1",0.,0.,1.,0.3);
	pad1->SetGridy();
	pad1->SetBottomMargin(0.2);
	pad1->Draw();
	c1->cd();
	TPad *pad2 = new TPad("pad2","pad2",0.,0.3,1.,1.);
	pad2->Draw();

	//Sig
	pad2->cd(0);pad2->SetLogy(1);
	ctauFrameSig->Draw();
	LifetimeLegendSig->Draw();

	left=0.5; top=0.9; textSize=0.030; latex->SetTextSize(textSize);
	if(nState == 4) latex->DrawLatex(left,top,"J/#psi");
	if(nState == 5) latex->DrawLatex(left,top,"#psi(2S)");
	top-=step;
	textSize=0.025; latex->SetTextSize(textSize);

	if(rapBin==0) latex->DrawLatex(left,top,Form("%.1f < |y| < %.1f",onia::rapForPTRange[rapBin],onia::rapForPTRange[2]));
	else if(rapBin==1) latex->DrawLatex(left,top,Form("|y| < %.1f",onia::rapForPTRange[rapBin]));
	else latex->DrawLatex(left,top,Form("%.1f < |y| < %.1f",onia::rapForPTRange[rapBin-1],onia::rapForPTRange[rapBin]));
	top-=step;
	if(ptBin==0)
		latex->DrawLatex(left,top,Form("%.1f < p_{T} < %.1f GeV",onia::pTRange[rapBin][ptBin],onia::pTRange[rapBin][11]));
	else
		latex->DrawLatex(left,top,Form("%.1f < p_{T} < %.1f GeV",onia::pTRange[rapBin][ptBin-1],onia::pTRange[rapBin][ptBin]));
	if(cpmBin==0)
		latex->DrawLatex(left,top-step,Form("%.1f < N_{ch} < %.1f",onia::cpmRange[cpmBin],onia::cpmRange[onia::NchBins-1]));
	else
		latex->DrawLatex(left,top-step,Form("%.1f < N_{ch} < %.1f",onia::cpmRange[cpmBin-1],onia::cpmRange[cpmBin]));			

	top-=step;
	latex->SetTextColor(kRed);
	latex->DrawLatex(left,top-step,"Signal Mass Region");
	latex->SetTextColor(kBlack);
	top-=step;
	latex->DrawLatex(left,top-step,Form("mean  =  %.1f",promptMean));
	top-=step;
	latex->DrawLatex(left,top-step,Form("resolution1 =  %.1f",promptCtRe));
	top-=step;
	latex->DrawLatex(left,top-step,Form("resolution2  =  %.1f",promptCtRe2));

	left=0.7, top=0.9;
	latex->DrawLatex(left,top,Form("#chi^{2}/ndof   =  %.2f/%d",chi2_LSig,ndof_LSig));
	top-=step;
	latex->DrawLatex(left,top,Form("frac_{Gauss2}   =  %.3f #pm %.3f",fracGauss2, fracGauss2Err));
	top-=step;
	latex->DrawLatex(left,top,Form("#tau_{NP}  =  %.3f #pm %.3f",nonPromptTau, nonPromptTauErr));
	top-=step;
	latex->DrawLatex(left,top,Form("frac_{P}   =  %.3f #pm %.3f",fPrompt, fPromptErr));
	top-=step;
	latex->DrawLatex(left,top,Form("frac_{Bg}   =  %.3f #pm %.3f",fBkg, fBkgErr));
	top-=step;
	top-=step;
	latex->DrawLatex(left,top,Form("BFrac   =  %.3f ",fNonPrompt/(fNonPrompt+fPrompt)));

	pad1->cd(0); pad1->SetLogy(0);
	ctauFrameSigPull->Draw();

	std::stringstream savectlog;
	//savectlog << "Fit/ct_SR_log_rap" << rapBin << "_pt" << ptBin << "_cpm" << cpmBin << "_pullAverage_withPull.pdf";
	savectlog << "Fit/ct_SR_log_rap" << rapBin << "_pt" << ptBin << "_cpm" << cpmBin << "_withPull.pdf";
	c1->SaveAs(savectlog.str().c_str());

	delete c1;
	delete legendBlue;
	delete legendBlueDash;
	delete legendRed;
	delete legendBlack;
	delete legendGreen;
	delete legendGreenDash;
	delete legendPink;
	return;
}

//==============================================
void plotLifeSig_linear(RooWorkspace *ws, int rapBin, int ptBin, int cpmBin, int nState){
	TGaxis::SetMaxDigits(3);
	//int nbins=100; //bin size 0.005 mm
	double binSize = 0.005; // mm 
	double xaxisMin = -0.2, xaxisMax = 0.3;
	int nbins = (int)((xaxisMax-xaxisMin)/binSize);

	RooRealVar JpsiMass(*ws->var("JpsiMass"));
	RooRealVar Jpsict(*ws->var("Jpsict"));
	RooRealVar *JpsictErr = ws->var("JpsictErr");

	RooPlot *ctauFrameSig=((RooRealVar*)ws->var("Jpsict"))->frame(Bins(nbins), Range(xaxisMin, xaxisMax));
	ctauFrameSig->SetName(Form("ctausig_plot_rap%d_pt%d_cpm%d",rapBin,ptBin,cpmBin));
	ctauFrameSig->SetYTitle("Events per 5 micron");
	ctauFrameSig->SetTitle("");
	ctauFrameSig->GetYaxis()->SetTitleOffset(1.3);
	ctauFrameSig->SetXTitle("lifetime [mm]");

	RooPlot *ctauFrameSigPull=((RooRealVar*)ws->var("Jpsict"))->frame(Bins(nbins), Range(xaxisMin, xaxisMax));
	ctauFrameSigPull->SetName(Form("pullctausig_plot_rap%d_pt%d_cpm%d",rapBin,ptBin,cpmBin));
	ctauFrameSigPull->SetTitle("");
	ctauFrameSigPull->SetXTitle("lifetime [mm]");
	ctauFrameSigPull->GetYaxis()->SetTitle("pull");
	ctauFrameSigPull->GetXaxis()->SetTitleSize(0.08);
	ctauFrameSigPull->GetYaxis()->SetTitleSize(0.08);
	ctauFrameSigPull->GetXaxis()->SetLabelSize(0.08);
	ctauFrameSigPull->GetYaxis()->SetLabelSize(0.08);
	ctauFrameSigPull->GetYaxis()->SetTitleOffset(0.4);
	ctauFrameSigPull->GetYaxis()->SetRangeUser(-5.5,5.5);


	RooDataSet *dataSR=(RooDataSet *)ws->data(Form("data_rap%d_pt%d_cpm%d_SR",rapBin,ptBin,cpmBin));

	RooFitResult* fitRlt = (RooFitResult*)ws->obj(Form("l_fitresult_rap%d_pt%d_cpm%d",rapBin,ptBin,cpmBin));
	if(!fitRlt) { cout<<">>=======Error: no Fit result object in workspace=========="<<endl; return; }
	fitRlt->Print();

	RooRealVar *promptMean_=(RooRealVar*)ws->var("promptMean");
	RooRealVar *ctResolution_=(RooRealVar*)ws->var("ctResolution");
	RooRealVar *ctResolution2_=(RooRealVar*)ws->var("ctResolution2");
	RooRealVar *fracGauss2_=(RooRealVar*)ws->var("fracGauss2");

	RooRealVar *nonPromptTau_=(RooRealVar*)ws->var("nonPromptTau");
	RooRealVar *fBkg_ = (RooRealVar*)ws->var("fBkg");
	RooRealVar *fPrompt_ = (RooRealVar*)ws->var("fPrompt");

	double promptMean = promptMean_->getVal();
	double promptMeanErr = promptMean_->getError();
	double promptCtRe = ctResolution_->getVal();
	double promptCtReErr = ctResolution_->getError();
	double promptCtRe2 = ctResolution2_->getVal();
	double promptCtRe2Err = ctResolution2_->getError();
	double fracGauss2 = fracGauss2_->getVal();
	double fracGauss2Err = fracGauss2_->getError();
	double nonPromptTau = nonPromptTau_->getVal();
	double nonPromptTauErr = nonPromptTau_->getError();
	double fBkg = fBkg_->getVal();
	double fBkgErr = fBkg_->getError();
	double fPrompt = fPrompt_->getVal();
	double fPromptErr = fPrompt_->getError();
	double fNonPrompt =  1.-fBkg-fPrompt;

	cout<<"promptMean: "<<promptMean<<endl;
	cout<<"promptCtRe: "<<promptCtRe<<endl;
	cout<<"nonPromptTau: "<<nonPromptTau<<endl;
	cout<<"fBkg: "<<fBkg<<endl;
	cout<<"fPrompt: "<<fPrompt<<endl;
	cout<<"fNonPrompt: "<<fNonPrompt<<endl;

	RooAddPdf *ModelLife = (RooAddPdf*)ws->pdf("fulllifetime");

	RooAddModel *Prompt = (RooAddModel*)ws->pdf("TotalPromptLifetime");
	RooAbsPdf *nonPromptSSD = (RooAbsPdf*)ws->pdf("nonPromptSSD");
	RooAbsPdf *backgroundlifetime = (RooAbsPdf*)ws->pdf("backgroundlifetime");

	int parsFit;
	//ploting signal region
	dataSR->plotOn(ctauFrameSig,MarkerSize(0.8));
	ModelLife->plotOn(ctauFrameSig,
			ProjWData(*JpsictErr, *dataSR),
			LineWidth(2), NumCPU(1));

	//------get chi2------------SHOULD DONE after PLOTTING------
	parsFit=(fitRlt->floatParsFinal()).getSize(); //this used the full p.d.f.

	int nBins_LSig=ctauFrameSig->GetNbinsX();
	double chi2Pre_LSig=ctauFrameSig->chiSquare(parsFit);  //reduced chi-squared = chi2/ndof
	int ndof_LSig=nBins_LSig-parsFit;  //num of degree of freedom
	double chi2_LSig=chi2Pre_LSig*ndof_LSig;

	RooHist* hpull_ctauSig = ctauFrameSig->pullHist() ;
	hpull_ctauSig->SetMarkerSize(0.8);
	for(int i=0;i<hpull_ctauSig->GetN();i++){
		hpull_ctauSig->SetPointEYlow(i,0.);
		hpull_ctauSig->SetPointEYhigh(i,0.);
	}
	ctauFrameSigPull->addPlotable(hpull_ctauSig,"P");

	Prompt->plotOn(ctauFrameSig,
			ProjWData(*JpsictErr, *dataSR),
			Normalization(fPrompt),
			LineStyle(5),
			LineColor(kBlue),
			LineWidth(2), NumCPU(1));

	nonPromptSSD->plotOn(ctauFrameSig,
			ProjWData(*JpsictErr, *dataSR),
			Normalization(fNonPrompt),
			LineStyle(2),
			LineColor(kRed),
			LineWidth(2), NumCPU(1));

	backgroundlifetime->plotOn(ctauFrameSig,
			ProjWData(*JpsictErr, *dataSR),
			Normalization(fBkg),
			LineStyle(7),
			LineColor(kPink+3),
			LineWidth(2), NumCPU(1));
	double Ymax = ctauFrameSig->GetMaximum();
	cout<<"Ymax: "<<Ymax<<endl;
	ctauFrameSig->SetMaximum(1.1*Ymax);

	TH1* legendBlue = dataSR->createHistogram("legendBlue",JpsiMass,Binning(50)) ; legendBlue->SetLineColor(kBlue) ; legendBlue->SetLineStyle(kSolid) ; legendBlue->SetLineWidth(2) ;
	TH1* legendBlueDash = dataSR->createHistogram("legendBlueDash",JpsiMass,Binning(50)) ; legendBlueDash->SetLineColor(kBlue) ; legendBlueDash->SetLineStyle(5) ; legendBlueDash->SetLineWidth(2) ;
	TH1* legendRed = dataSR->createHistogram("legendRed",JpsiMass,Binning(50)) ; legendRed->SetLineColor(kRed) ; legendRed->SetLineStyle(2) ; legendRed->SetLineWidth(2) ;
	TH1* legendBlack = dataSR->createHistogram("legendBlack",JpsiMass,Binning(50)) ; legendBlack->SetLineColor(kBlack) ; legendBlack->SetLineStyle(kSolid) ; legendBlack->SetLineWidth(2) ;
	TH1* legendGreen = dataSR->createHistogram("legendGreen",JpsiMass,Binning(50)) ; legendGreen->SetLineColor(kGreen) ; legendGreen->SetLineStyle(kSolid) ; legendGreen->SetLineWidth(2) ;
	TH1* legendGreenDash = dataSR->createHistogram("legendGreenDash",JpsiMass,Binning(50)) ; legendGreenDash->SetLineColor(kGreen) ; legendGreenDash->SetLineStyle(2) ; legendGreenDash->SetLineWidth(2) ;
	TH1* legendPink = dataSR->createHistogram("legendPink",JpsiMass,Binning(50)) ; legendPink->SetLineColor(kPink+3) ; legendPink->SetLineStyle(7) ; legendPink->SetLineWidth(2) ;

	TLegend* LifetimeLegendSig=new TLegend(0.12,0.8,0.25,0.93);
	LifetimeLegendSig->SetFillColor(kWhite);
	LifetimeLegendSig->SetTextFont(42);
	LifetimeLegendSig->SetTextSize(0.035);
	LifetimeLegendSig->SetBorderSize(0.);
	LifetimeLegendSig->AddEntry(legendBlue,"sum","l");
	LifetimeLegendSig->AddEntry(legendBlueDash,"prompt","l");
	LifetimeLegendSig->AddEntry(legendRed,"non-prompt","l");
	LifetimeLegendSig->AddEntry(legendPink,"background","l");

	double left=0.75, top=0.9, textSize=0.032;
	TLatex *latex=new TLatex();
	latex->SetTextFont(42);
	latex->SetNDC(kTRUE);
	latex->SetTextSize(textSize);
	double step=textSize*1.3;

	gStyle->SetPadBottomMargin(0.11); //0.12
	gStyle->SetPadLeftMargin(0.1); //0.12
	gStyle->SetPadRightMargin(0.02); //0.05
	gStyle->SetPadTopMargin(0.05); //0.05

	TCanvas *c1=new TCanvas("c1","",800,900);

	c1->cd();
	TPad *pad1 = new TPad("pad1","pad1",0.,0.,1.,0.3);
	pad1->SetGridy();
	pad1->SetBottomMargin(0.2);
	pad1->Draw();
	c1->cd();
	TPad *pad2 = new TPad("pad2","pad2",0.,0.3,1.,1.);
	pad2->Draw();

	//Sig
	pad2->cd(0);pad2->SetLogy(0);
	ctauFrameSig->Draw();
	LifetimeLegendSig->Draw();

	left=0.75; top=0.9; textSize=0.030; latex->SetTextSize(textSize);
	if(nState == 4) latex->DrawLatex(left,top,"J/#psi");
	if(nState == 5) latex->DrawLatex(left,top,"#psi(2S)");
	top-=step;
	textSize=0.025; latex->SetTextSize(textSize);

	if(rapBin==0) latex->DrawLatex(left,top,Form("%.1f < |y| < %.1f",onia::rapForPTRange[rapBin],onia::rapForPTRange[2]));
	else if(rapBin==1) latex->DrawLatex(left,top,Form("|y| < %.1f",onia::rapForPTRange[rapBin]));
	else latex->DrawLatex(left,top,Form("%.1f < |y| < %.1f",onia::rapForPTRange[rapBin-1],onia::rapForPTRange[rapBin]));
	top-=step;
	if(ptBin==0)
		latex->DrawLatex(left,top,Form("%.1f < p_{T} < %.1f GeV",onia::pTRange[rapBin][ptBin],onia::pTRange[rapBin][11]));
	else
		latex->DrawLatex(left,top,Form("%.1f < p_{T} < %.1f GeV",onia::pTRange[rapBin][ptBin-1],onia::pTRange[rapBin][ptBin]));
	if(cpmBin==0)
		latex->DrawLatex(left,top-step,Form("%.1f < N_{ch} < %.1f",onia::cpmRange[cpmBin],onia::cpmRange[onia::NchBins-1]));
	else
		latex->DrawLatex(left,top-step,Form("%.1f < N_{ch} < %.1f",onia::cpmRange[cpmBin-1],onia::cpmRange[cpmBin]));			

	top-=step;
	latex->SetTextColor(kRed);
	latex->DrawLatex(left,top-step,"Signal Mass Region");
	latex->SetTextColor(kBlack);
	top-=step;

	pad1->cd(0); pad1->SetLogy(0);
	ctauFrameSigPull->Draw();

	std::stringstream savect;
	savect << "Fit/ct_SR_rap" << rapBin << "_pt" << ptBin << "_cpm" << cpmBin << "_withPull.pdf";
	c1->SaveAs(savect.str().c_str());

	delete c1;
	delete legendBlue;
	delete legendBlueDash;
	delete legendRed;
	delete legendBlack;
	delete legendGreen;
	delete legendGreenDash;
	delete legendPink;
	return;
}

//==============================================
void plotLifeBgIndividual(RooWorkspace *ws, int rapBin, int ptBin, int cpmBin, int nState){
	//int nbins=140; //bin size 0.025 mm
	double binSize = 0.025; // mm 
	double xaxisMin = -0.5, xaxisMax = 2.5;
	int nbins = (int)((xaxisMax-xaxisMin)/binSize);


	RooRealVar JpsiMass(*ws->var("JpsiMass"));
	RooRealVar *JpsictErr = ws->var("JpsictErr");

	RooPlot *ctauFrameBkgSBL=((RooRealVar*)ws->var("Jpsict"))->frame(Bins(nbins), Range(xaxisMin, xaxisMax));
	ctauFrameBkgSBL->SetName(Form("ctaubkg_plot_rap%d_pt%d_cpm%d",rapBin,ptBin,cpmBin));
	ctauFrameBkgSBL->SetYTitle("Events per 25 micron");
	ctauFrameBkgSBL->SetTitle("");
	ctauFrameBkgSBL->GetXaxis()->SetTitle("lifetime [mm]");

	RooPlot *ctauFrameBkgSBLPull=((RooRealVar*)ws->var("Jpsict"))->frame(Bins(nbins), Range(xaxisMin, xaxisMax));
	ctauFrameBkgSBLPull->SetName(Form("pullctaubkg_plot_rap%d_pt%d_cpm%d",rapBin,ptBin,cpmBin));
	ctauFrameBkgSBLPull->SetTitle("");
	ctauFrameBkgSBLPull->GetXaxis()->SetTitle("lifetime [mm]");
	ctauFrameBkgSBLPull->GetYaxis()->SetTitle("pull");
	ctauFrameBkgSBLPull->GetXaxis()->SetTitleSize(0.08);
	ctauFrameBkgSBLPull->GetYaxis()->SetTitleSize(0.08);
	ctauFrameBkgSBLPull->GetXaxis()->SetLabelSize(0.08);
	ctauFrameBkgSBLPull->GetYaxis()->SetLabelSize(0.08);
	ctauFrameBkgSBLPull->GetYaxis()->SetTitleOffset(0.4);
	ctauFrameBkgSBLPull->GetYaxis()->SetRangeUser(-5.5,5.5);

	RooPlot *ctauFrameBkgSBR=((RooRealVar*)ws->var("Jpsict"))->frame(Bins(nbins), Range(xaxisMin, xaxisMax));
	ctauFrameBkgSBR->SetName(Form("ctaubkg_plot_rap%d_pt%d_cpm%d",rapBin,ptBin,cpmBin));
	ctauFrameBkgSBR->SetYTitle("Events per 25 micron");
	ctauFrameBkgSBR->SetTitle("");
	ctauFrameBkgSBR->GetXaxis()->SetTitle("lifetime [mm]");

	RooPlot *ctauFrameBkgSBRPull=((RooRealVar*)ws->var("Jpsict"))->frame(Bins(nbins), Range(xaxisMin, xaxisMax));
	ctauFrameBkgSBRPull->SetName(Form("pullctaubkg_plot_rap%d_pt%d_cpm%d",rapBin,ptBin,cpmBin));
	ctauFrameBkgSBRPull->SetTitle("");
	ctauFrameBkgSBRPull->GetXaxis()->SetTitle("lifetime [mm]");
	ctauFrameBkgSBRPull->GetYaxis()->SetTitle("pull");
	ctauFrameBkgSBRPull->GetXaxis()->SetTitleSize(0.08);
	ctauFrameBkgSBRPull->GetYaxis()->SetTitleSize(0.08);
	ctauFrameBkgSBRPull->GetXaxis()->SetLabelSize(0.08);
	ctauFrameBkgSBRPull->GetYaxis()->SetLabelSize(0.08);
	ctauFrameBkgSBRPull->GetYaxis()->SetTitleOffset(0.4);
	ctauFrameBkgSBRPull->GetYaxis()->SetRangeUser(-5.5,5.5);

	RooDataSet *dataSBL=(RooDataSet *)ws->data(Form("data_rap%d_pt%d_cpm%d_SBL",rapBin,ptBin,cpmBin));
	RooDataSet *dataSBR=(RooDataSet *)ws->data(Form("data_rap%d_pt%d_cpm%d_SBR",rapBin,ptBin,cpmBin));

	RooFitResult* fitRlt = (RooFitResult*)ws->obj(Form("l_fitresult_rap%d_pt%d_cpm%d",rapBin,ptBin,cpmBin));
	if(!fitRlt) { cout<<">>=======Error: no Fit result object in workspace=========="<<endl; return; }
	fitRlt->Print();

	RooRealVar *promptMean_=(RooRealVar*)ws->var("promptMean");
	RooRealVar *ctResolution_=(RooRealVar*)ws->var("ctResolution");
	RooRealVar *ctResolution2_=(RooRealVar*)ws->var("ctResolution2");
	RooRealVar *fracGauss2_=(RooRealVar*)ws->var("fracGauss2");

	RooRealVar *bkgTauSSD_SBL_=(RooRealVar*)ws->var("bkgTauSSD_SBL");
	RooRealVar *bkgTauFD_=(RooRealVar*)ws->var("bkgTauFD");
	RooRealVar *bkgTauDSD_=(RooRealVar*)ws->var("bkgTauDSD");

	RooRealVar *bkgTauSSD_SBR_=(RooRealVar*)ws->var("bkgTauSSD_SBR");

	RooRealVar *fBkgSSDR_SBL_ = (RooRealVar*)ws->var("fBkgSSDR_SBL");
	RooRealVar *fBkgDSD_SBL_ = (RooRealVar*)ws->var("fBkgDSD_SBL");
	RooRealVar *fBkgSSDR_SBR_ = (RooRealVar*)ws->var("fBkgSSDR_SBR");
	RooRealVar *fBkgDSD_SBR_ = (RooRealVar*)ws->var("fBkgDSD_SBR");

	RooRealVar *fBkg_ = (RooRealVar*)ws->var("fBkg");

	RooRealVar *fBkgSBL_ = (RooRealVar*)ws->var("fBkgSBL");
	RooRealVar *fBkgSBR_ = (RooRealVar*)ws->var("fBkgSBR");

	RooRealVar *fPrompt_ = (RooRealVar*)ws->var("fPrompt");

	double promptMean = promptMean_->getVal();
	double promptMeanErr = promptMean_->getError();
	double promptCtRe = ctResolution_->getVal();
	double promptCtReErr = ctResolution_->getError();
	double promptCtRe2 = ctResolution2_->getVal();
	double promptCtRe2Err = ctResolution2_->getError();
	double fracGauss2 = fracGauss2_->getVal();
	double fracGauss2Err = fracGauss2_->getError();

	double bkgTauSSD_SBL = bkgTauSSD_SBL_->getVal();
	double bkgTauSSD_SBLErr = bkgTauSSD_SBL_->getError();
	double bkgTauFD = bkgTauFD_->getVal();
	double bkgTauFDErr = bkgTauFD_->getError();
	double bkgTauDSD = bkgTauDSD_->getVal();
	double bkgTauDSDErr = bkgTauDSD_->getError();

	double bkgTauSSD_SBR = bkgTauSSD_SBR_->getVal();
	double bkgTauSSD_SBRErr = bkgTauSSD_SBR_->getError();


	double fBkg = fBkg_->getVal();
	double fBkgErr = fBkg_->getError();

	double fBkgSBL = fBkgSBL_->getVal();
	double fBkgSBR = fBkgSBR_->getVal();
	double fBkgSBLErr = fBkgSBL_->getError();
	double fBkgSBRErr = fBkgSBR_->getError();

	double fPrompt = fPrompt_->getVal();
	double fPromptErr = fPrompt_->getError();
	double fNonPrompt =  1.-fBkg-fPrompt;
	double fNonPromptErr = sqrt(pow(fBkgErr,2)+pow(fPromptErr,2));

	//f_P / f_NP should be same in LSB,SR,RSB, to interpolate f_P in L(R)SB from SR
	double fPromptSBL = fPrompt * ( 1. - fBkgSBL ) / ( 1. - fBkg ) ;
	double fPromptSBR = fPrompt * ( 1. - fBkgSBR ) / ( 1. - fBkg ) ;
	double fPromptSBLErr = sqrt( 
			pow(fPromptErr*(1.-fBkgSBL )/(1.-fBkg),2) + 
			pow(fBkgSBLErr*fPrompt/(1.-fBkg),2) + 
			pow(fBkgErr*fPrompt*(1.-fBkgSBL )/pow(1.-fBkg,2),2) );
	double fPromptSBRErr = sqrt( 
			pow(fPromptErr*(1.-fBkgSBR )/(1.-fBkg),2) + 
			pow(fBkgSBRErr*fPrompt/(1.-fBkg),2) + 
			pow(fBkgErr*fPrompt*(1.-fBkgSBR )/pow(1.-fBkg,2),2) );

	double fNonPromptSBL = 1. - fBkgSBL - fPromptSBL ;
	double fNonPromptSBR = 1. - fBkgSBR - fPromptSBR ;
	double fNonPromptSBLErr = sqrt(
			pow(fPromptErr*(1.-fBkgSBL )/(1.-fBkg),2) +
			pow(fBkgSBLErr*(1.-fPrompt/(1.-fBkg)),2) +
			pow(fBkgErr*fPrompt*(1.-fBkgSBL )/pow(1.-fBkg,2),2) );
	double fNonPromptSBRErr = sqrt(
			pow(fPromptErr*(1.-fBkgSBR )/(1.-fBkg),2) +
			pow(fBkgSBRErr*(1.-fPrompt/(1.-fBkg)),2) +
			pow(fBkgErr*fPrompt*(1.-fBkgSBR )/pow(1.-fBkg,2),2) );

	cout<<"fBkgSBL: "<<fBkgSBL<<" fBkgSBR: "<<fBkgSBR<<endl;
	cout<<"fPromptSBL: "<<fPromptSBL<<" fPromptSBR: "<<fPromptSBR <<endl;
	cout<<"fNonPromptSBL: "<<fNonPromptSBL<<" fNonPromptSBR: "<<fNonPromptSBL<<endl;
	cout<<"fPromptSBL/fNonPromptSBL: "<<fPromptSBL/fNonPromptSBL<<endl;
	cout<<"fPromptSBR/fNonPromptSBR: "<<fPromptSBR/fNonPromptSBR<<endl;
	cout<<"fPrompt/fNonPrompt: "<<fPrompt/fNonPrompt<<endl;

	double fBkgSSDR_SBL = fBkgSSDR_SBL_->getVal();
	double fBkgSSDRErr_SBL = fBkgSSDR_SBL_->getError();
	double fBkgDSD_SBL = fBkgDSD_SBL_->getVal();
	double fBkgDSDErr_SBL = fBkgDSD_SBL_->getError();
	double fBkgSSDL_SBL = 1-fBkgSSDR_SBL-fBkgDSD_SBL;

	double fBkgSSDR_SBR = fBkgSSDR_SBR_->getVal();
	double fBkgSSDRErr_SBR = fBkgSSDR_SBR_->getError();
	double fBkgDSD_SBR = fBkgDSD_SBR_->getVal();
	double fBkgDSDErr_SBR = fBkgDSD_SBR_->getError();
	double fBkgSSDL_SBR = 1-fBkgSSDR_SBR-fBkgDSD_SBR;

	cout<<"promptMean: "<<promptMean<<endl;
	cout<<"promptCtRe: "<<promptCtRe<<endl;
	cout<<"fBkgSSDR_SBL: "<<fBkgSSDR_SBL<<endl;
	cout<<"fBkgSSDL_SBL: "<<fBkgSSDL_SBL<<endl;
	cout<<"fBkgSSDR_SBR: "<<fBkgSSDR_SBR<<endl;
	cout<<"fBkgSSDL_SBR: "<<fBkgSSDL_SBR<<endl;

	RooAddPdf *ModelLifeSBL = (RooAddPdf*)ws->pdf("backgroundlifetimeL");
	RooAddPdf *ModelLifeSBR = (RooAddPdf*)ws->pdf("backgroundlifetimeR");
	RooAddPdf *backgroundSSD_SBL = (RooAddPdf*)ws->pdf("backgroundSSD_SBL");
	RooAddPdf *backgroundSSD_SBR = (RooAddPdf*)ws->pdf("backgroundSSD_SBR");
	RooAddPdf *backgroundDSD = (RooAddPdf*)ws->pdf("backgroundDSD");
	RooAddPdf *backgroundFD = (RooAddPdf*)ws->pdf("backgroundFD");

	RooAddModel *Prompt = (RooAddModel*)ws->pdf("TotalPromptLifetime");
	RooAbsPdf *nonPromptSSD = (RooAbsPdf*)ws->pdf("nonPromptSSD");


	//relative fraction of Prompt: f_P/f_NP, same in LSB,SR,RSB
	RooRealVar *fracPrompt = new RooRealVar("fracPrompt","fracPrompt",fPromptSBL/(fPromptSBL+fNonPromptSBL));
	RooAddPdf *SignalPdf = new RooAddPdf("SignalPdf","SignalPdf",RooArgList(*Prompt,*nonPromptSSD),*fracPrompt);


	int parsFit;

	//ploting background SBL
	dataSBL->plotOn(ctauFrameBkgSBL,MarkerSize(0.8));
	ModelLifeSBL->plotOn(ctauFrameBkgSBL,
			ProjWData(*JpsictErr, *dataSBL),
			LineWidth(2),
			NumCPU(1));

	parsFit=(fitRlt->floatParsFinal()).getSize(); //this used the full p.d.f.
	int nBins_LSBL=ctauFrameBkgSBL->GetNbinsX();
	double chi2Pre_LSBL=ctauFrameBkgSBL->chiSquare(parsFit);  //reduced chi-squared = chi2/ndof
	int ndof_LSBL=nBins_LSBL-parsFit;  //num of degree of freedom
	double chi2_LSBL=chi2Pre_LSBL*ndof_LSBL;

	RooHist* hpull_ctauBkgSBL = ctauFrameBkgSBL->pullHist() ;
	hpull_ctauBkgSBL->SetMarkerSize(0.8);
	for(int i=0;i<hpull_ctauBkgSBL->GetN();i++){
		hpull_ctauBkgSBL->SetPointEYlow(i,0.);
		hpull_ctauBkgSBL->SetPointEYhigh(i,0.);
	}
	cout<<"hpull_ctauBkgSBL->GetN(): "<<hpull_ctauBkgSBL->GetN()<<endl;

	ctauFrameBkgSBLPull->addPlotable(hpull_ctauBkgSBL,"P");

	backgroundSSD_SBL->plotOn(ctauFrameBkgSBL,
			ProjWData(*JpsictErr, *dataSBL),
			Normalization(fBkgSBL*fBkgSSDR_SBL),
			LineStyle(2),
			LineColor(kRed),
			LineWidth(2),
			NumCPU(1));

	backgroundFD->plotOn(ctauFrameBkgSBL,
			ProjWData(*JpsictErr, *dataSBL),
			Normalization(fBkgSBL*(1.-fBkgDSD_SBL-fBkgSSDR_SBL)),
			LineStyle(2),
			LineColor(kGreen),
			LineWidth(2),
			NumCPU(1));

	backgroundDSD->plotOn(ctauFrameBkgSBL,
			ProjWData(*JpsictErr, *dataSBL),
			Normalization(fBkgSBL*fBkgDSD_SBL),
			LineStyle(1),
			LineColor(kBlack),
			LineWidth(2),
			NumCPU(1));

	SignalPdf->plotOn(ctauFrameBkgSBL,
			ProjWData(*JpsictErr, *dataSBL),
			Normalization(1.-fBkgSBL),
			LineStyle(7),
			LineColor(kPink+3),
			LineWidth(2));

	double Ymax = ctauFrameBkgSBL->GetMaximum();
	cout<<"Ymax: "<<Ymax<<endl;

	ctauFrameBkgSBL->SetMaximum(5*Ymax);
	ctauFrameBkgSBL->SetMinimum(1.1);

	dataSBR->plotOn(ctauFrameBkgSBR,MarkerSize(0.8));
	ModelLifeSBR->plotOn(ctauFrameBkgSBR,
			ProjWData(*JpsictErr, *dataSBR),
			LineWidth(2),
			NumCPU(1));

	parsFit=(fitRlt->floatParsFinal()).getSize(); //this used the full p.d.f.
	int nBins_LSBR=ctauFrameBkgSBR->GetNbinsX();
	double chi2Pre_LSBR=ctauFrameBkgSBR->chiSquare(parsFit);  //reduced chi-squared = chi2/ndof
	int ndof_LSBR=nBins_LSBR-parsFit;  //num of degree of freedom
	double chi2_LSBR=chi2Pre_LSBR*ndof_LSBR;

	RooHist* hpull_ctauBkgSBR = ctauFrameBkgSBR->pullHist() ;
	hpull_ctauBkgSBR->SetMarkerSize(0.8);
	for(int i=0;i<hpull_ctauBkgSBR->GetN();i++){
		hpull_ctauBkgSBR->SetPointEYlow(i,0.);
		hpull_ctauBkgSBR->SetPointEYhigh(i,0.);
	}
	ctauFrameBkgSBRPull->addPlotable(hpull_ctauBkgSBR,"P");

	backgroundSSD_SBR->plotOn(ctauFrameBkgSBR,
			ProjWData(*JpsictErr, *dataSBR),
			Normalization(fBkgSBR*fBkgSSDR_SBR),
			LineStyle(2),
			LineColor(kRed),
			LineWidth(2),
			NumCPU(1));

	backgroundFD->plotOn(ctauFrameBkgSBR,
			ProjWData(*JpsictErr, *dataSBR),
			Normalization(fBkgSBR*(1.-fBkgSSDR_SBR-fBkgDSD_SBR)),
			LineStyle(2),
			LineColor(kGreen),
			LineWidth(2),
			NumCPU(1));

	backgroundDSD->plotOn(ctauFrameBkgSBR,
			ProjWData(*JpsictErr, *dataSBR),
			Normalization(fBkgSBR*fBkgDSD_SBR),
			LineStyle(1),
			LineColor(kBlack),
			LineWidth(2),
			NumCPU(1));

	SignalPdf->plotOn(ctauFrameBkgSBR,
			ProjWData(*JpsictErr, *dataSBR),
			Normalization(1.-fBkgSBR),
			LineStyle(7),
			LineColor(kPink+3),
			LineWidth(2));

	Ymax = ctauFrameBkgSBR->GetMaximum();
	cout<<"Ymax: "<<Ymax<<endl;
	ctauFrameBkgSBR->SetMaximum(5*Ymax);
	ctauFrameBkgSBR->SetMinimum(1.1);

	TH1* legendBlue = dataSBL->createHistogram("legendBlue",JpsiMass,Binning(50)) ; legendBlue->SetLineColor(kBlue) ; legendBlue->SetLineStyle(kSolid) ; legendBlue->SetLineWidth(2) ;
	TH1* legendBlueDash = dataSBL->createHistogram("legendBlueDash",JpsiMass,Binning(50)) ; legendBlueDash->SetLineColor(kBlue) ; legendBlueDash->SetLineStyle(5) ; legendBlueDash->SetLineWidth(2) ;
	TH1* legendRed = dataSBL->createHistogram("legendRed",JpsiMass,Binning(50)) ; legendRed->SetLineColor(kRed) ; legendRed->SetLineStyle(2) ; legendRed->SetLineWidth(2) ;
	TH1* legendBlack = dataSBL->createHistogram("legendBlack",JpsiMass,Binning(50)) ; legendBlack->SetLineColor(kBlack) ; legendBlack->SetLineStyle(kSolid) ; legendBlack->SetLineWidth(2) ;
	TH1* legendGreen = dataSBL->createHistogram("legendGreen",JpsiMass,Binning(50)) ; legendGreen->SetLineColor(kGreen) ; legendGreen->SetLineStyle(kSolid) ; legendGreen->SetLineWidth(2) ;
	TH1* legendGreenDash = dataSBL->createHistogram("legendGreenDash",JpsiMass,Binning(50)) ; legendGreenDash->SetLineColor(kGreen) ; legendGreenDash->SetLineStyle(2) ; legendGreenDash->SetLineWidth(2) ;
	TH1* legendPink = dataSBL->createHistogram("legendPink",JpsiMass,Binning(50)) ; legendPink->SetLineColor(kPink+3) ; legendPink->SetLineStyle(7) ; legendPink->SetLineWidth(2) ;

	TLegend* LifetimeLegendBkgSBL=new TLegend(0.35,0.7,0.5,0.93);
	LifetimeLegendBkgSBL->SetFillColor(kWhite);
	LifetimeLegendBkgSBL->SetTextFont(42);
	LifetimeLegendBkgSBL->SetTextSize(0.035);
	LifetimeLegendBkgSBL->SetBorderSize(0.);
	LifetimeLegendBkgSBL->AddEntry(legendBlue,"sum","l");
	LifetimeLegendBkgSBL->AddEntry(legendGreenDash,"SSL","l");
	LifetimeLegendBkgSBL->AddEntry(legendRed,"SSR ","l");
	LifetimeLegendBkgSBL->AddEntry(legendBlack,"DS","l");
	LifetimeLegendBkgSBL->AddEntry(legendPink,"signal","l");

	TLegend* LifetimeLegendBkgSBR=new TLegend(0.35,0.7,0.5,0.93);
	LifetimeLegendBkgSBR->SetFillColor(kWhite);
	LifetimeLegendBkgSBR->SetTextFont(42);
	LifetimeLegendBkgSBR->SetTextSize(0.035);
	LifetimeLegendBkgSBR->SetBorderSize(0.);
	LifetimeLegendBkgSBR->AddEntry(legendBlue,"sum","l");
	LifetimeLegendBkgSBR->AddEntry(legendGreenDash,"SSL","l");
	LifetimeLegendBkgSBR->AddEntry(legendRed,"SSR ","l");
	LifetimeLegendBkgSBR->AddEntry(legendBlack,"DS","l");
	LifetimeLegendBkgSBR->AddEntry(legendPink,"signal","l");

	double left=0.7, top=0.9, textSize=0.025;
	TLatex *latex=new TLatex();
	latex->SetTextFont(42);
	latex->SetNDC(kTRUE);
	latex->SetTextSize(textSize);
	double step=textSize*1.3;

	gStyle->SetPadBottomMargin(0.11); //0.12
	gStyle->SetPadLeftMargin(0.08); //0.12
	gStyle->SetPadRightMargin(0.02); //0.05
	gStyle->SetPadTopMargin(0.02); //0.05

	TCanvas *c1=new TCanvas("c1","c1",800,700);
	c1->SetBottomMargin(0.08);

	TCanvas *c2=new TCanvas("c2","c2",800,300);
	c2->SetBottomMargin(0.2);

	//Bkg SBL
	c1->cd(0);c1->SetLogy(1);
	ctauFrameBkgSBL->Draw();
	LifetimeLegendBkgSBL->Draw();
	left=0.5; top=0.9; textSize=0.030; latex->SetTextSize(textSize);
	if(nState == 4) latex->DrawLatex(left,top,"J/#psi");
	if(nState == 5) latex->DrawLatex(left,top,"#psi(2S)");
	top-=step;
	textSize=0.025; latex->SetTextSize(textSize);

	if(rapBin==0) latex->DrawLatex(left,top,Form("%.1f < |y| < %.1f",onia::rapForPTRange[rapBin],onia::rapForPTRange[2]));
	else if(rapBin==1) latex->DrawLatex(left,top,Form("|y| < %.1f",onia::rapForPTRange[rapBin]));
	else latex->DrawLatex(left,top,Form("%.1f < |y| < %.1f",onia::rapForPTRange[rapBin-1],onia::rapForPTRange[rapBin]));
	top-=step;
	if(ptBin==0)
		latex->DrawLatex(left,top,Form("%.1f < p_{T} < %.1f GeV",onia::pTRange[rapBin][ptBin],onia::pTRange[rapBin][11]));
	else
		latex->DrawLatex(left,top,Form("%.1f < p_{T} < %.1f GeV",onia::pTRange[rapBin][ptBin-1],onia::pTRange[rapBin][ptBin]));
	if(cpmBin==0)
		latex->DrawLatex(left,top-step,Form("%.1f < N_{ch} < %.1f",onia::cpmRange[cpmBin],onia::cpmRange[onia::NchBins-1]));
	else
		latex->DrawLatex(left,top-step,Form("%.1f < N_{ch} < %.1f",onia::cpmRange[cpmBin-1],onia::cpmRange[cpmBin]));		

	top-=step;
	latex->SetTextColor(kRed);
	latex->DrawLatex(left,top-step,"Left Mass SideBand");
	latex->SetTextColor(kBlack);
	top-=step;
	latex->DrawLatex(left,top-step,Form("mean  =  %.1f",promptMean));
	top-=step;
	latex->DrawLatex(left,top-step,Form("resolution1 =  %.1f",promptCtRe));
	top-=step;
	latex->DrawLatex(left,top-step,Form("resolution2  =  %.1f",promptCtRe2));
	top-=step;
	if(ptBin>7)
		latex->DrawLatex(left,top,Form("frac_{SSR}   =  %.3f",fBkgSSDR_SBL));

	left=0.7; top=0.9;
	latex->DrawLatex(left,top,Form("#chi^{2}/ndof   =  %.2f/%d",chi2_LSBL,ndof_LSBL));
	top-=step;
	latex->DrawLatex(left,top,Form("frac_{DS}   =  %.3f #pm %.3f",fBkgDSD_SBL, fBkgDSDErr_SBL));
	top-=step;
	if(ptBin<8){
		latex->DrawLatex(left,top,Form("frac_{SSR}   =  %.3f #pm %.3f",fBkgSSDR_SBL, fBkgSSDRErr_SBL));
		top-=step;
	}
	latex->DrawLatex(left,top,Form("#tau_{SSR}   =  %.3f #pm %.3f",bkgTauSSD_SBL, bkgTauSSD_SBLErr));
	top-=step;
	latex->DrawLatex(left,top,Form("frac_{Bg}   =  %.3f #pm %.3f",fBkgSBL,fBkgSBLErr));
	top-=step;
	top-=step;

	latex->DrawLatex(left,top,Form("frac_{Gauss2}   =  %.3f #pm %.3f",fracGauss2, fracGauss2Err));
	top-=step;
	latex->DrawLatex(left,top,Form("#tau_{DS}   =  %.3f #pm %.3f",bkgTauDSD, bkgTauDSDErr));
	top-=step;
	latex->DrawLatex(left,top,Form("#tau_{SSL}  =  %.3f #pm %.3f",bkgTauFD, bkgTauFDErr));

	std::stringstream savectLSB;
	savectLSB << "Fit/ct_LSB_rap" << rapBin << "_pt" << ptBin << "_cpm" << cpmBin << ".pdf";
	c1->SaveAs(savectLSB.str().c_str());

	// pull
	c2->cd(0); c2->SetLogy(0);
	ctauFrameBkgSBLPull->Draw();

	left=0.6; top=0.9; textSize=0.055; latex->SetTextSize(textSize);
	std::stringstream pullLatexLSB;
	if(nState == 4) pullLatexLSB << "J/#psi, ";
	if(nState == 5) pullLatexLSB << "#psi(2S), ";
	if(rapBin==0) pullLatexLSB << Form("%.1f < |y| < %.1f, ",onia::rapForPTRange[rapBin],onia::rapForPTRange[2]);
	else if(rapBin==1) pullLatexLSB << Form("|y| < %.1f, ",onia::rapForPTRange[rapBin]);
	else pullLatexLSB << Form("%.1f < |y| < %.1f, ",onia::rapForPTRange[rapBin-1],onia::rapForPTRange[rapBin]);
	if(ptBin==0) pullLatexLSB << Form("%.1f < p_{T} < %.1f GeV",onia::pTRange[rapBin][ptBin],onia::pTRange[rapBin][11]);
	else pullLatexLSB << Form("%.1f < N_{ch} < %.1f GeV",onia::pTRange[rapBin][ptBin-1],onia::pTRange[rapBin][ptBin]);
	if(cpmBin==0) pullLatexLSB << Form("%.1f < N_{ch} < %.1f",onia::cpmRange[cpmBin],onia::cpmRange[onia::NchBins-1]);
	else pullLatexLSB << Form("%.1f < n_{ch} < %.1f",onia::cpmRange[cpmBin-1],onia::cpmRange[cpmBin]);

	latex->DrawLatex(left,top,pullLatexLSB.str().c_str());

	std::stringstream savectLSBPull;
	savectLSBPull << "Fit/ct_LSB_rap" << rapBin << "_pt" << ptBin << "_cpm" << cpmBin << "_pull.pdf";
	c2->SaveAs(savectLSBPull.str().c_str());
	//---------------------------------------------------------------------------

	//Bkg SBR
	c1->cd(0);c1->SetLogy(1);
	ctauFrameBkgSBR->Draw();
	LifetimeLegendBkgSBR->Draw();
	left=0.5; top=0.9; textSize=0.030; latex->SetTextSize(textSize);
	if(nState == 4) latex->DrawLatex(left,top,"J/#psi");
	if(nState == 5) latex->DrawLatex(left,top,"#psi(2S)");
	top-=step;
	textSize=0.025; latex->SetTextSize(textSize);

	if(rapBin==0) latex->DrawLatex(left,top,Form("%.1f < |y| < %.1f",onia::rapForPTRange[rapBin],onia::rapForPTRange[2]));
	else if(rapBin==1) latex->DrawLatex(left,top,Form("|y| < %.1f",onia::rapForPTRange[rapBin]));
	else latex->DrawLatex(left,top,Form("%.1f < |y| < %.1f",onia::rapForPTRange[rapBin-1],onia::rapForPTRange[rapBin]));
	top-=step;
	if(ptBin==0)
		latex->DrawLatex(left,top,Form("%.1f < p_{T} < %.1f GeV",onia::pTRange[rapBin][ptBin],onia::pTRange[rapBin][11]));
	else
		latex->DrawLatex(left,top,Form("%.1f < p_{T} < %.1f GeV",onia::pTRange[rapBin][ptBin-1],onia::pTRange[rapBin][ptBin]));
	if(cpmBin==0)
		latex->DrawLatex(left,top-step,Form("%.1f < N_{ch} < %.1f",onia::cpmRange[cpmBin],onia::cpmRange[onia::NchBins-1]));
	else
		latex->DrawLatex(left,top-step,Form("%.1f < N_{ch} < %.1f",onia::cpmRange[cpmBin-1],onia::cpmRange[cpmBin]));		

	top-=step;
	latex->SetTextColor(kRed);
	latex->DrawLatex(left,top-step,"Right Mass SideBand");
	latex->SetTextColor(kBlack);
	top-=step;

	left=0.7; top=0.9;
	latex->DrawLatex(left,top,Form("#chi^{2}/ndof   =  %.2f/%d",chi2_LSBR,ndof_LSBR));
	top-=step;
	latex->DrawLatex(left,top,Form("frac_{DS}   =  %.3f #pm %.3f",fBkgDSD_SBR, fBkgDSDErr_SBR));
	top-=step;
	if(ptBin<8){
		latex->DrawLatex(left,top,Form("frac_{SSR}   =  %.3f #pm %.3f",fBkgSSDR_SBR, fBkgSSDRErr_SBR));
		top-=step;
	}
	latex->DrawLatex(left,top,Form("#tau_{SSR}   =  %.3f #pm %.3f",bkgTauSSD_SBR, bkgTauSSD_SBRErr));
	top-=step;
	latex->DrawLatex(left,top,Form("frac_{Bg}   =  %.3f #pm %.3f",fBkgSBR,fBkgSBRErr));
	top-=step;
	top-=step;

	std::stringstream savectRSB;
	savectRSB << "Fit/ct_RSB_rap" << rapBin << "_pt" << ptBin << "_cpm" << cpmBin << ".pdf";
	c1->SaveAs(savectRSB.str().c_str());

	c2->cd(0); c2->SetLogy(0);
	ctauFrameBkgSBRPull->Draw();

	left=0.6; top=0.9; textSize=0.055; latex->SetTextSize(textSize);
	std::stringstream pullLatexRSB;
	if(nState == 4) pullLatexRSB << "J/#psi, ";
	if(nState == 5) pullLatexRSB << "#psi(2S), ";
	if(rapBin==0) pullLatexRSB << Form("%.1f < |y| < %.1f, ",onia::rapForPTRange[rapBin],onia::rapForPTRange[2]);
	else if(rapBin==1) pullLatexRSB << Form("|y| < %.1f, ",onia::rapForPTRange[rapBin]);
	else pullLatexRSB << Form("%.1f < |y| < %.1f, ",onia::rapForPTRange[rapBin-1],onia::rapForPTRange[rapBin]);
	if(ptBin==0) pullLatexRSB << Form("%.1f < p_{T} < %.1f GeV",onia::pTRange[rapBin][ptBin],onia::pTRange[rapBin][11]);
	else pullLatexRSB << Form("%.1f < p_{T} < %.1f GeV",onia::pTRange[rapBin][ptBin-1],onia::pTRange[rapBin][ptBin]);
	if(cpmBin==0) pullLatexRSB << Form("%.1f < N_{ch} < %.1f",onia::cpmRange[cpmBin],onia::cpmRange[onia::NchBins-1]);
	else pullLatexRSB << Form("%.1f < n_{ch} < %.1f",onia::cpmRange[cpmBin-1],onia::cpmRange[cpmBin]);

	latex->DrawLatex(left,top,pullLatexRSB.str().c_str());

	std::stringstream savectRSBPull;
	savectRSBPull << "Fit/ct_RSB_rap" << rapBin << "_pt" << ptBin << "_cpm" << cpmBin << "_pull.pdf";
	c2->SaveAs(savectRSBPull.str().c_str());

	delete c1;
	delete c2;
	delete legendBlue;
	delete legendBlueDash;
	delete legendRed;
	delete legendBlack;
	delete legendGreen;
	delete legendGreenDash;
	delete legendPink;
	return;
}

//==============================================
void plotLifeSigIndividual(RooWorkspace *ws, int rapBin, int ptBin, int cpmBin, int nState){
	//int nbins=230; //bin size 0.01 mm
	double binSize = 0.01; // mm 
	double xaxisMin = -0.5, xaxisMax = 2.5;
	int nbins = (int)((xaxisMax-xaxisMin)/binSize);

	RooRealVar JpsiMass(*ws->var("JpsiMass"));
	RooRealVar Jpsict(*ws->var("Jpsict"));
	RooRealVar *JpsictErr = ws->var("JpsictErr");

	RooPlot *ctauFrameSig=((RooRealVar*)ws->var("Jpsict"))->frame(Bins(nbins), Range(xaxisMin, xaxisMax));
	ctauFrameSig->SetName(Form("ctausig_plot_rap%d_pt%d_cpm%d",rapBin,ptBin,cpmBin));
	ctauFrameSig->SetYTitle("Events per 10 micron");
	ctauFrameSig->SetTitle("");
	ctauFrameSig->SetXTitle("lifetime [mm]");

	RooPlot *ctauFrameSigPull=((RooRealVar*)ws->var("Jpsict"))->frame(Bins(nbins), Range(xaxisMin, xaxisMax));
	ctauFrameSigPull->SetName(Form("pullctausig_plot_rap%d_pt%d_cpm%d",rapBin,ptBin,cpmBin));
	ctauFrameSigPull->SetTitle("");
	ctauFrameSigPull->SetXTitle("lifetime [mm]");
	ctauFrameSigPull->GetYaxis()->SetTitle("pull");
	ctauFrameSigPull->GetXaxis()->SetTitleSize(0.08);
	ctauFrameSigPull->GetYaxis()->SetTitleSize(0.08);
	ctauFrameSigPull->GetXaxis()->SetLabelSize(0.08);
	ctauFrameSigPull->GetYaxis()->SetLabelSize(0.08);
	ctauFrameSigPull->GetYaxis()->SetTitleOffset(0.4);
	ctauFrameSigPull->GetYaxis()->SetRangeUser(-5.5,5.5);


	RooDataSet *dataSR=(RooDataSet *)ws->data(Form("data_rap%d_pt%d_cpm%d_SR",rapBin,ptBin,cpmBin));

	RooFitResult* fitRlt = (RooFitResult*)ws->obj(Form("l_fitresult_rap%d_pt%d_cpm%d",rapBin,ptBin,cpmBin));
	if(!fitRlt) { cout<<">>=======Error: no Fit result object in workspace=========="<<endl; return; }
	fitRlt->Print();
	cout<<"dataSR->numEntries: "<<dataSR->numEntries()<<endl;

	RooRealVar *promptMean_=(RooRealVar*)ws->var("promptMean");
	RooRealVar *ctResolution_=(RooRealVar*)ws->var("ctResolution");
	RooRealVar *ctResolution2_=(RooRealVar*)ws->var("ctResolution2");
	RooRealVar *fracGauss2_=(RooRealVar*)ws->var("fracGauss2");

	RooRealVar *nonPromptTau_=(RooRealVar*)ws->var("nonPromptTau");
	RooRealVar *fBkg_ = (RooRealVar*)ws->var("fBkg");
	RooRealVar *fPrompt_ = (RooRealVar*)ws->var("fPrompt");

	double promptMean = promptMean_->getVal();
	double promptMeanErr = promptMean_->getError();
	double promptCtRe = ctResolution_->getVal();
	double promptCtReErr = ctResolution_->getError();
	double promptCtRe2 = ctResolution2_->getVal();
	double promptCtRe2Err = ctResolution2_->getError();
	double fracGauss2 = fracGauss2_->getVal();
	double fracGauss2Err = fracGauss2_->getError();
	double nonPromptTau = nonPromptTau_->getVal();
	double nonPromptTauErr = nonPromptTau_->getError();

	double fBkg = fBkg_->getVal();
	double fBkgErr = fBkg_->getError();
	double fPrompt = fPrompt_->getVal();
	double fPromptErr = fPrompt_->getError();
	double fNonPrompt =  1.-fBkg-fPrompt;

	cout<<"promptMean: "<<promptMean<<endl;
	cout<<"promptCtRe: "<<promptCtRe<<endl;
	cout<<"nonPromptTau: "<<nonPromptTau<<endl;
	cout<<"fBkg: "<<fBkg<<endl;
	cout<<"fBkgErr: "<<fBkgErr<<endl;
	cout<<"fPrompt: "<<fPrompt<<endl;
	cout<<"fNonPrompt: "<<fNonPrompt<<endl;

	RooAddPdf *ModelLife = (RooAddPdf*)ws->pdf("fulllifetime");

	RooAddModel *Prompt = (RooAddModel*)ws->pdf("TotalPromptLifetime");
	RooAbsPdf *nonPromptSSD = (RooAbsPdf*)ws->pdf("nonPromptSSD");
	RooAbsPdf *backgroundlifetime = (RooAbsPdf*)ws->pdf("backgroundlifetime");

	int parsFit;
	//ploting signal region
	dataSR->plotOn(ctauFrameSig,MarkerSize(0.8), Name("myHist"));
	ModelLife->plotOn(ctauFrameSig,
			ProjWData(*JpsictErr, *dataSR),
			LineWidth(2), NumCPU(1),
			Name("myCurve"));

	//------get chi2------------SHOULD DONE after PLOTTING------
	parsFit=(fitRlt->floatParsFinal()).getSize(); //this used the full p.d.f.

	int nBins_LSig=ctauFrameSig->GetNbinsX();
	double chi2Pre_LSig=ctauFrameSig->chiSquare(parsFit);  //reduced chi-squared = chi2/ndof
	int ndof_LSig=nBins_LSig-parsFit;  //num of degree of freedom
	double chi2_LSig=chi2Pre_LSig*ndof_LSig;

	RooHist* hpull_ctauSig = ctauFrameSig->pullHist() ;
	//RooHist* hpull_ctauSig = ctauFrameSig->pullHist("myHist","myCurve",kTRUE);

	//save pull distribution
	TH1F* pull = new TH1F("pull","pull distribution", 100,-10.,10.);
	gSystem->mkdir("Fit/pull",kTRUE);
	//TFile *pullFile = new TFile(Form("Fit/pull/pullAverage_rap%d_pt%d_cpm%d_SR.root",rapBin,ptBincpmBin),"RECREATE");
	TFile *pullFile = new TFile(Form("Fit/pull/pull_rap%d_pt%d_cpm%d_SR.root",rapBin,ptBin,cpmBin),"RECREATE");

	hpull_ctauSig->SetMarkerSize(0.8);
	for(int i=0;i<hpull_ctauSig->GetN();i++){
		hpull_ctauSig->SetPointEYlow(i,0.);
		hpull_ctauSig->SetPointEYhigh(i,0.);
		double x,y;
		hpull_ctauSig->GetPoint(i,x,y);
		pull->Fill(y);
	}
	pullFile->cd();
	pull->Write();
	pullFile->Close();

	ctauFrameSigPull->addPlotable(hpull_ctauSig,"P");

	Prompt->plotOn(ctauFrameSig,
			ProjWData(*JpsictErr, *dataSR),
			Normalization(fPrompt),
			LineStyle(5),
			LineColor(kBlue),
			LineWidth(2), NumCPU(1));

	nonPromptSSD->plotOn(ctauFrameSig,
			ProjWData(*JpsictErr, *dataSR),
			Normalization(fNonPrompt),
			LineStyle(2),
			LineColor(kRed),
			LineWidth(2), NumCPU(1));

	backgroundlifetime->plotOn(ctauFrameSig,
			ProjWData(*JpsictErr, *dataSR),
			Normalization(fBkg),
			LineStyle(7),
			LineColor(kPink+3),
			LineWidth(2), NumCPU(1));

	double Ymax = ctauFrameSig->GetMaximum();
	cout<<"Ymax: "<<Ymax<<endl;
	ctauFrameSig->SetMaximum(3*Ymax);
	ctauFrameSig->SetMinimum(1.1);


	TH1* legendBlue = dataSR->createHistogram("legendBlue",JpsiMass,Binning(50)) ; legendBlue->SetLineColor(kBlue) ; legendBlue->SetLineStyle(kSolid) ; legendBlue->SetLineWidth(2) ;
	TH1* legendBlueDash = dataSR->createHistogram("legendBlueDash",JpsiMass,Binning(50)) ; legendBlueDash->SetLineColor(kBlue) ; legendBlueDash->SetLineStyle(5) ; legendBlueDash->SetLineWidth(2) ;
	TH1* legendRed = dataSR->createHistogram("legendRed",JpsiMass,Binning(50)) ; legendRed->SetLineColor(kRed) ; legendRed->SetLineStyle(2) ; legendRed->SetLineWidth(2) ;
	TH1* legendBlack = dataSR->createHistogram("legendBlack",JpsiMass,Binning(50)) ; legendBlack->SetLineColor(kBlack) ; legendBlack->SetLineStyle(kSolid) ; legendBlack->SetLineWidth(2) ;
	TH1* legendGreen = dataSR->createHistogram("legendGreen",JpsiMass,Binning(50)) ; legendGreen->SetLineColor(kGreen) ; legendGreen->SetLineStyle(kSolid) ; legendGreen->SetLineWidth(2) ;
	TH1* legendGreenDash = dataSR->createHistogram("legendGreenDash",JpsiMass,Binning(50)) ; legendGreenDash->SetLineColor(kGreen) ; legendGreenDash->SetLineStyle(2) ; legendGreenDash->SetLineWidth(2) ;
	TH1* legendPink = dataSR->createHistogram("legendPink",JpsiMass,Binning(50)) ; legendPink->SetLineColor(kPink+3) ; legendPink->SetLineStyle(7) ; legendPink->SetLineWidth(2) ;

	TLegend* LifetimeLegendSig=new TLegend(0.31,0.8,0.43,0.96);
	LifetimeLegendSig->SetFillColor(kWhite);
	LifetimeLegendSig->SetTextFont(42);
	LifetimeLegendSig->SetTextSize(0.035);
	LifetimeLegendSig->SetBorderSize(0.);
	LifetimeLegendSig->AddEntry(legendBlue,"sum","l");
	LifetimeLegendSig->AddEntry(legendBlueDash,"prompt","l");
	LifetimeLegendSig->AddEntry(legendRed,"non-prompt","l");
	LifetimeLegendSig->AddEntry(legendPink,"background","l");

	double left=0.5, top=0.9, textSize=0.025;
	TLatex *latex=new TLatex();
	latex->SetTextFont(42);
	latex->SetNDC(kTRUE);
	latex->SetTextSize(textSize);
	double step=textSize*1.3;

	gStyle->SetPadBottomMargin(0.11); //0.12
	gStyle->SetPadLeftMargin(0.08); //0.12
	gStyle->SetPadRightMargin(0.02); //0.05
	gStyle->SetPadTopMargin(0.02); //0.05

	TCanvas *c1=new TCanvas("c1","c1",800,700);
	c1->SetBottomMargin(0.08);

	TCanvas *c2=new TCanvas("c2","c2",800,300);
	c2->SetBottomMargin(0.2);

	//Sig
	c1->cd(0);c1->SetLogy(1);
	ctauFrameSig->Draw();
	LifetimeLegendSig->Draw();

	left=0.5; top=0.9; textSize=0.030; latex->SetTextSize(textSize);
	if(nState == 4) latex->DrawLatex(left,top,"J/#psi");
	if(nState == 5) latex->DrawLatex(left,top,"#psi(2S)");
	top-=step;
	textSize=0.025; latex->SetTextSize(textSize);

	if(rapBin==0) latex->DrawLatex(left,top,Form("%.1f < |y| < %.1f",onia::rapForPTRange[rapBin],onia::rapForPTRange[2]));
	else if(rapBin==1) latex->DrawLatex(left,top,Form("|y| < %.1f",onia::rapForPTRange[rapBin]));
	else latex->DrawLatex(left,top,Form("%.1f < |y| < %.1f",onia::rapForPTRange[rapBin-1],onia::rapForPTRange[rapBin]));
	top-=step;
	if(ptBin==0)
		latex->DrawLatex(left,top,Form("%.1f < p_{T} < %.1f GeV",onia::pTRange[rapBin][ptBin],onia::pTRange[rapBin][11]));
	else
		latex->DrawLatex(left,top,Form("%.1f < p_{T} < %.1f GeV",onia::pTRange[rapBin][ptBin-1],onia::pTRange[rapBin][ptBin]));
	if(cpmBin==0)
		latex->DrawLatex(left,top-step,Form("%.1f < N_{ch} < %.1f",onia::cpmRange[cpmBin],onia::cpmRange[onia::NchBins-1]));
	else
		latex->DrawLatex(left,top-step,Form("%.1f < N_{ch} < %.1f",onia::cpmRange[cpmBin-1],onia::cpmRange[cpmBin]));			

	top-=step;
	latex->SetTextColor(kRed);
	latex->DrawLatex(left,top-step,"Signal Mass Region");
	latex->SetTextColor(kBlack);
	top-=step;
	latex->DrawLatex(left,top-step,Form("mean  =  %.1f",promptMean));
	top-=step;
	latex->DrawLatex(left,top-step,Form("resolution1 =  %.1f",promptCtRe));
	top-=step;
	latex->DrawLatex(left,top-step,Form("resolution2  =  %.1f",promptCtRe2));

	left=0.7, top=0.9;
	latex->DrawLatex(left,top,Form("#chi^{2}/ndof   =  %.2f/%d",chi2_LSig,ndof_LSig));
	top-=step;
	latex->DrawLatex(left,top,Form("frac_{Gauss2}   =  %.3f #pm %.3f",fracGauss2, fracGauss2Err));
	top-=step;
	latex->DrawLatex(left,top,Form("#tau_{NP}  =  %.3f #pm %.3f",nonPromptTau, nonPromptTauErr));
	top-=step;
	latex->DrawLatex(left,top,Form("frac_{P}   =  %.3f #pm %.3f",fPrompt, fPromptErr));
	top-=step;
	latex->DrawLatex(left,top,Form("frac_{Bg}   =  %.3f #pm %.3f",fBkg, fBkgErr));
	top-=step;
	top-=step;
	latex->DrawLatex(left,top,Form("BFrac   =  %.3f ",fNonPrompt/(fNonPrompt+fPrompt)));

	std::stringstream savectlog;
	savectlog << "Fit/ct_SR_log_rap" << rapBin << "_pt" << ptBin << "_cpm" << cpmBin << ".pdf";
	c1->SaveAs(savectlog.str().c_str());

	c2->cd(0); c2->SetLogy(0);
	ctauFrameSigPull->Draw();

	left=0.6; top=0.9; textSize=0.055; latex->SetTextSize(textSize);
	std::stringstream pullLatexSR;
	if(nState == 4) pullLatexSR << "J/#psi, ";
	if(nState == 5) pullLatexSR << "#psi(2S), ";
	if(rapBin==0) pullLatexSR << Form("%.1f < |y| < %.1f, ",onia::rapForPTRange[rapBin],onia::rapForPTRange[2]);
	else if(rapBin==1) pullLatexSR << Form("|y| < %.1f, ",onia::rapForPTRange[rapBin]);
	else pullLatexSR << Form("%.1f < |y| < %.1f, ",onia::rapForPTRange[rapBin-1],onia::rapForPTRange[rapBin]);
	if(ptBin==0) pullLatexSR << Form("%.1f < p_{T} < %.1f GeV",onia::pTRange[rapBin][ptBin],onia::pTRange[rapBin][11]);
	else pullLatexSR << Form("%.1f < p_{T} < %.1f GeV",onia::pTRange[rapBin][ptBin-1],onia::pTRange[rapBin][ptBin]);
	if(cpmBin==0) pullLatexSR << Form("%.1f < N_{ch} < %.1f",onia::cpmRange[cpmBin],onia::cpmRange[onia::NchBins-1]);
	else pullLatexSR << Form("%.1f < n_{ch} < %.1f",onia::cpmRange[cpmBin-1],onia::cpmRange[cpmBin]);
	
	latex->DrawLatex(left,top,pullLatexSR.str().c_str());

	std::stringstream savectlogPull;
	savectlogPull << "Fit/ct_SR_log_rap" << rapBin << "_pt" << ptBin << "_cpm" << cpmBin << "_pull.pdf";
	c2->SaveAs(savectlogPull.str().c_str());

	delete c1;
	delete c2;
	delete legendBlue;
	delete legendBlueDash;
	delete legendRed;
	delete legendBlack;
	delete legendGreen;
	delete legendGreenDash;
	delete legendPink;
	return;
}
