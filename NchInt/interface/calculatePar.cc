#include "rootIncludes.inc"
#include "commonVar.h"
#include "RooUtils.h"

using namespace RooFit;

double getMeanPt(int rapBin, int ptBin, const std::string &massfile){
    RooWorkspace* ws = getFromTFile<RooWorkspace>(massfile, "ws_masslifetime");

    // Get stuff from dataset (make sure it is not 0)
    RooRealVar *JpsiMass = ws->var("JpsiMass");
    RooRealVar *JpsiPt = ws->var("JpsiPt");
    RooAbsData *data = ws->data(Form("data_rap%d_pt%d",rapBin,ptBin));
    assert( 0 != JpsiMass && 0 != JpsiPt && 0 != data);

    TH1* histPt = data->createHistogram("histPt", *JpsiPt, Binning(120));
    TH1* histMass = data->createHistogram("histMass", *JpsiMass, Binning(120));

    double meanPt = histPt->GetMean();
    double meanMass = histPt->GetMean();

		delete histPt;
		delete histMass;
    return meanPt;
}

double getFracBkg(int rapBin, int ptBin, const std::string &massfile){
    RooWorkspace* ws = getFromTFile<RooWorkspace>(massfile, "ws_masslifetime");

    // Get stuff from dataset (make sure it is not 0)
    RooRealVar *fBkg_ = ws->var("fracBkg");
    assert( 0 != fBkg_ );

    return fBkg_->getVal();
}

double getFracBkgIn3Sigma(int rapBin, int ptBin, const std::string &massfile, int nSigMass){
    RooWorkspace* ws = getFromTFile<RooWorkspace>(massfile, "ws_masslifetime");

    RooAbsData *data=  ws->data(Form("data_rap%d_pt%d",rapBin,ptBin));
    RooRealVar *JpsiMass = ws->var("JpsiMass");

    RooRealVar *CBmass =  ws->var("CBmass");
    RooRealVar *CBsigma=  ws->var("CBsigma");
    RooRealVar *CBsigma2= ws->var("CBsigma2");
    RooRealVar *fracCB1_= ws->var("fracCB1");
    RooRealVar *fracBkg_= ws->var("fracBkg");

    double Mean = CBmass->getVal();
    double Sigma = CBsigma->getVal();
    double SigmaErr = CBsigma->getError();
    double Sigma2 = CBsigma2->getVal();
    double Sigma2Err = CBsigma2->getError();
    double fracCB1 = fracCB1_->getVal();
    double fracCB1Err = fracCB1_->getError();
    double fracBkg = fracBkg_->getVal();
    double fracBkgErr = fracBkg_->getError();
    double SigmaWei = sqrt(pow(Sigma,2)*fracCB1+pow(Sigma2,2)*(1-fracCB1));
    double SigmaWeiErr = sqrt(pow((Sigma-Sigma2)*fracCB1Err,2)+pow(fracCB1*SigmaErr,2)+pow((1-fracCB1)*Sigma2Err,2)); //to correct

    double sigMaxMass = Mean+SigmaWei*nSigMass;
    double sigMinMass = Mean-SigmaWei*nSigMass;
    double sbHighMass = Mean+SigmaWei*onia::nSigBkgHigh;
    double sbLowMass =  Mean-SigmaWei*onia::nSigBkgLow;

    //calculate the background fraction in 3*sigma signal region
    int nEntries = data->numEntries();
    JpsiMass->setRange("SigRegion",sigMinMass,sigMaxMass);
    RooAbsPdf *massPdf = ws->pdf("massModel");
    RooAbsPdf *bkgMassShape = ws->pdf("bkgMassShape");
    RooRealVar *fracFull3Sig = (RooRealVar*)massPdf->createIntegral(*JpsiMass,NormSet(*JpsiMass),Range("SigRegion"));
    RooRealVar *fracBkg3Sig = (RooRealVar*)bkgMassShape->createIntegral(*JpsiMass,NormSet(*JpsiMass),Range("SigRegion"));
    double FracFull3Sig = fracFull3Sig->getVal();
    double FracFull3SigErr = fracFull3Sig->getError();
    double FracBkg3Sig = fracBkg3Sig->getVal();
    double FracBkg3SigErr = fracBkg3Sig->getError();

    double evtFull3Sig = nEntries*FracFull3Sig;
    double evtBkg3Sig = nEntries*fracBkg*FracBkg3Sig;
    double BkgRatio3Sig = evtBkg3Sig/evtFull3Sig; //fracBkg*FracBkg3Sig/FracFull3Sig
    double BkgRatio3SigErr = BkgRatio3Sig*
        sqrt(pow(fracBkgErr/fracBkg,2)+pow(FracBkg3SigErr/FracBkg3Sig,2)+pow(FracFull3SigErr/FracFull3Sig,2));

    return BkgRatio3Sig;
}

double getFracBkgErr(int rapBin, int ptBin, const std::string &massfile){
    TFile *inFile=new TFile(massfile.c_str(),"R");
    if(!inFile) {std::cout<<"Error: failed to opne root file"<<std::endl;}
    RooWorkspace *ws=(RooWorkspace *)inFile->Get("ws_masslifetime");
    if(!ws) {std::cout<<"Error: failed to opne workspace "<<std::endl;}

    RooRealVar *fBkg_;
    fBkg_ = (RooRealVar*)ws->var("fracBkg");
    if(!fBkg_){std::cout<<"not needed variable in the workspace"<<std::endl; return -1;}

    double FracBkgErr = fBkg_->getError();
    return FracBkgErr;
}

double getFracBkgErrIn3Sigma(int rapBin, int ptBin, const std::string &massfile, int nSigMass){
    TFile *inFile=new TFile(massfile.c_str(),"R");
    if(!inFile) {std::cout<<"Error: failed to opne root file"<<std::endl;}
    RooWorkspace *ws=(RooWorkspace *)inFile->Get("ws_masslifetime");
    if(!ws) {std::cout<<"Error: failed to opne workspace "<<std::endl;}
    RooDataSet *data=(RooDataSet *)ws->data(Form("data_rap%d_pt%d",rapBin,ptBin));
    RooRealVar JpsiMass(*ws->var("JpsiMass"));

    RooRealVar *CBmass=(RooRealVar *)ws->var("CBmass");
    RooRealVar *CBsigma=(RooRealVar *)ws->var("CBsigma");
    RooRealVar *CBsigma2=(RooRealVar *)ws->var("CBsigma2");
    RooRealVar *fracCB1_=(RooRealVar *)ws->var("fracCB1");
    RooRealVar *fracBkg_=(RooRealVar *)ws->var("fracBkg");
    double Mean = CBmass->getVal();
    double Sigma = CBsigma->getVal();
    double SigmaErr = CBsigma->getError();
    double Sigma2 = CBsigma2->getVal();
    double Sigma2Err = CBsigma2->getError();
    double fracCB1 = fracCB1_->getVal();
    double fracCB1Err = fracCB1_->getError();
    double fracBkg = fracBkg_->getVal();
    double fracBkgErr = fracBkg_->getError();
    double SigmaWei = sqrt(pow(Sigma,2)*fracCB1+pow(Sigma2,2)*(1-fracCB1));
    double SigmaWeiErr = sqrt(pow((Sigma-Sigma2)*fracCB1Err,2)+pow(fracCB1*SigmaErr,2)+pow((1-fracCB1)*Sigma2Err,2)); //to correct

    double sigMaxMass = Mean+SigmaWei*nSigMass;
    double sigMinMass = Mean-SigmaWei*nSigMass;
    double sbHighMass = Mean+SigmaWei*onia::nSigBkgHigh;
    double sbLowMass =  Mean-SigmaWei*onia::nSigBkgLow;

    //calculate the background fraction in 3*sigma signal region
    int nEntries = data->numEntries();
    JpsiMass.setRange("SigRegion",sigMinMass,sigMaxMass);
    RooAddPdf *massPdf = (RooAddPdf*)ws->pdf("massModel");
    RooAddPdf *bkgMassShape = (RooAddPdf*)ws->pdf("bkgMassShape");
    RooRealVar* fracFull3Sig=(RooRealVar*)massPdf->createIntegral(JpsiMass,NormSet(JpsiMass),Range("SigRegion"));
    RooRealVar* fracBkg3Sig=(RooRealVar*)bkgMassShape->createIntegral(JpsiMass,NormSet(JpsiMass),Range("SigRegion"));
    double FracFull3Sig = fracFull3Sig->getVal();
    double FracFull3SigErr = fracFull3Sig->getError();
    double FracBkg3Sig = fracBkg3Sig->getVal();
    double FracBkg3SigErr = fracBkg3Sig->getError();

    double evtFull3Sig = nEntries*FracFull3Sig;
    double evtBkg3Sig = nEntries*fracBkg*FracBkg3Sig;
    double BkgRatio3Sig = evtBkg3Sig/evtFull3Sig; //fracBkg*FracBkg3Sig/FracFull3Sig
    double BkgRatio3SigErr = BkgRatio3Sig*
        sqrt(pow(fracBkgErr/fracBkg,2)+pow(FracBkg3SigErr/FracBkg3Sig,2)+pow(FracFull3SigErr/FracFull3Sig,2));

    return BkgRatio3SigErr;
}
