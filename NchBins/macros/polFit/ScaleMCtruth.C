#include "Riostream.h"
#include "TSystem.h"
#include "TStyle.h"
#include "TLine.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TCanvas.h"
#include "TString.h"
#include "TFile.h"
#include "TROOT.h"
#include "TTree.h"
#include "TFrame.h"
#include "TMath.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "TRandom3.h"
#include "TF1.h"
#include "TF2.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TGraphAsymmErrors.h"
#include "TGraph2D.h"
#include "TGraph.h"
#include "TEfficiency.h"
#include "TFitterMinuit.h"
#include "TBox.h"
#include "TMath.h"
#include "TMatrixDSym.h"
#include "TMatrixDSymEigen.h"
#include "TMatrixD.h"
#include "TMatrix.h"

#include <vector>

#include "ScaleChi2Func.h"

void ScaleMCtruth(){

	gROOT->SetBatch();

	gStyle->SetPalette(1,0);
	gStyle->SetPadBottomMargin(0.12);
	gStyle->SetPadLeftMargin(0.13);
	gStyle->SetPadRightMargin(0.15);

	gStyle->SetTickLength(-0.02, "xyz");
	gStyle->SetLabelOffset(0.02, "x");
	gStyle->SetLabelOffset(0.02, "y");
	gStyle->SetTitleOffset(1.3, "x");
	gStyle->SetTitleOffset(1.4, "y");
	gStyle->SetTitleFillColor(kWhite);




	delete gRandom;
	gRandom = new TRandom3(0);  // better random generator

	char EffFile[200];
	char EffFileData[200];
	char graphName[200];

	sprintf(EffFile,"EffFiles/singleMuTruthEff_30April2013_40GeVrap1_startFromTracks_noKinemCuts_Dimuon0JpsiNoVertexing_GENvariables.root");
	//sprintf(EffFile,"EffFiles/singleMuTruthEff_7March2012_40GeVrap1_TnPBins.root");
	sprintf(EffFileData,"EffFiles/EfficiencyFactorized_Dimuon0Jpsi_combined_DATA_MC_Trk80Cuts_14Mar2012.root");
	//sprintf(EffFileData,"EffFiles/EfficiencyFactorized_Dimuon0Jpsi_combinedMC_DATA_run1_Trk80Cuts_scaled_14June2012.root");
	//sprintf(EffFileData,"EffFiles/EfficiencyFactorized_Dimuon0Jpsi_combinedMC_DATA_run1_Trk80Cuts_scaled_sanity_drM1_newFactor_21June2012.root");
	//sprintf(EffFileData,"EffFiles/EfficiencyFactorized_Dimuon0Jpsi_combinedMC_DATA_run1_Trk80Cuts_scaled_sanity_Mixed_24June2012.root");

	char Date[100];
	sprintf(Date,"May20");
	char JobID[200];
	sprintf(JobID,"May20_GENvariables");
	//sprintf(JobID,"June27_SOFT_Approval_changedScaling");
	char filename[200];
	sprintf(filename,"Scaling/%s",JobID);
	gSystem->mkdir("Scaling");
	gSystem->mkdir(filename);
	const int etaBins=8;
	double etaRange[9]={0,0.2,0.3,0.6,0.8,1,1.2,1.4,1.6};
	Int_t binY = 0;
	double edmEta[etaBins];

	int pTBinsNew = 1000;

	TFile *fInEff = new TFile(EffFile);
	TFile *fInEffData = new TFile(EffFileData);

	sprintf(filename,"Scaling/%s/Sigmas.txt",JobID);
	ofstream myfile;
	myfile.open (filename);





	/////////// Loop over eta bins ////////////////////////////////////


	for(int etaBin=0;etaBin<etaBins;etaBin++){
		int etaBinMin=0;
		int etaBinMax=7;
		if(etaBin<etaBinMin||etaBin>etaBinMax) continue;
		cout<<"working on bin "<<etaRange[etaBin]<<" < eta < "<<etaRange[etaBin+1]<<endl;

		///////// GET MC AND DATA EFFICIENCIES IN THE CORRECT FORMAT ////////////////////////

		sprintf(graphName,"totEff_MCTRUTH_PT_AETA%d",etaBin);
		TEfficiency* TEffMC=(TEfficiency*) fInEff->Get(graphName);

		TH1* hEffTOT=(TH1*)TEffMC->GetTotalHistogram();
		TH1* hEffPASS=(TH1*)TEffMC->GetPassedHistogram();
		hEffPASS->Divide(hEffTOT);

		hEffPASS->Print();
		TEffMC=(TEfficiency*) fInEff->Get(graphName);

		int nBinspT_=hEffPASS->GetNbinsX();
		int nBinsBelow20GeV=0;
		int nBinsAbove20GeV=0;
		int nBinsBelowInfo=0;
		double meanEffAbove20GeV=0;

		for(int pTBin=0;pTBin<nBinspT_;pTBin++){
			if(TEffMC->GetEfficiency(TEffMC->FindFixBin(hEffPASS->GetXaxis()->GetBinCenter(pTBin)))<0.000001) nBinsBelowInfo=pTBin+2;
			if(hEffPASS->GetXaxis()->GetBinCenter(pTBin)<20) {nBinsBelow20GeV=pTBin+1;}
			if(hEffPASS->GetXaxis()->GetBinCenter(pTBin)>20){
				nBinsAbove20GeV=pTBin+1-nBinsBelow20GeV;
				binY=TEffMC->FindFixBin(hEffPASS->GetXaxis()->GetBinCenter(pTBin));
				meanEffAbove20GeV+=TEffMC->GetEfficiency(binY);
			}
		}

		meanEffAbove20GeV=meanEffAbove20GeV/nBinsAbove20GeV;
		cout<<"meanEffAbove20GeV "<<meanEffAbove20GeV<<endl;
		cout<<"nBinscheck "<<nBinsBelow20GeV+nBinsAbove20GeV<<" vs real "<<nBinspT_<<endl;


		const int nBinspT=nBinspT_-nBinsAbove20GeV+1-nBinsBelowInfo;
		double ptCentre[nBinspT];
		double effMean[nBinspT];
		double err_effMean_low[nBinspT];
		double err_effMean_high[nBinspT];

		int histBinstart=nBinsBelowInfo;
		for(int pTBin=0;pTBin<nBinsBelow20GeV-nBinsBelowInfo;pTBin++){
			ptCentre[pTBin]=hEffPASS->GetXaxis()->GetBinCenter(histBinstart);
			binY=TEffMC->FindFixBin(ptCentre[pTBin]);
			effMean[pTBin] = TEffMC->GetEfficiency(binY);
			err_effMean_low[pTBin] = TEffMC->GetEfficiencyErrorLow(binY);
			err_effMean_high[pTBin] = TEffMC->GetEfficiencyErrorUp(binY);
			//				  cout<<"ptCentre[pTBin] "<<ptCentre[pTBin]<<endl;
			//				  cout<<"binY "<<binY<<endl;
			//				  cout<<"effMean "<<effMean[pTBin]<<endl;
			histBinstart++;
		}

		ptCentre[nBinspT-1]=100;
		//effMean[nBinspT-1] = meanEffAbove20GeV;
		effMean[nBinspT-1] = effMean[nBinspT-2];
		err_effMean_low[nBinspT-1] = 0;
		err_effMean_high[nBinspT-1] = 0;


		TGraphAsymmErrors *graphMC = new TGraphAsymmErrors(nBinspT,ptCentre,effMean,0,0,err_effMean_low,err_effMean_high);
		//graphMC->SaveAs("graphMC.root");

		bool removeStrangePoints=false;
		if(removeStrangePoints){
			if(etaBin==0) graphMC->RemovePoint(12);
			if(etaBin==0) graphMC->RemovePoint(17);
			if(etaBin==0) graphMC->RemovePoint(1);

			if(etaBin==1) graphMC->RemovePoint(17);

			if(etaBin==2) graphMC->RemovePoint(10);
			if(etaBin==2) graphMC->RemovePoint(11);
			if(etaBin==2) graphMC->RemovePoint(11);
			if(etaBin==2) graphMC->RemovePoint(12);

			if(etaBin==3) graphMC->RemovePoint(16);
			if(etaBin==3) graphMC->RemovePoint(16);

			if(etaBin==4) graphMC->RemovePoint(9);
			if(etaBin==4) graphMC->RemovePoint(11);
			if(etaBin==4) graphMC->RemovePoint(11);
			if(etaBin==4) graphMC->RemovePoint(14);
			if(etaBin==4) graphMC->RemovePoint(14);
			if(etaBin==4) graphMC->RemovePoint(15);

			if(etaBin==5) graphMC->RemovePoint(12);
			if(etaBin==5) graphMC->RemovePoint(14);
			if(etaBin==5) graphMC->RemovePoint(14);

			if(etaBin==6) graphMC->RemovePoint(15);
			if(etaBin==6) graphMC->RemovePoint(15);
			if(etaBin==6) graphMC->RemovePoint(15);
			if(etaBin==6) graphMC->RemovePoint(16);
			if(etaBin==6) graphMC->RemovePoint(17);
			if(etaBin==6) graphMC->RemovePoint(18);
			if(etaBin==6) graphMC->RemovePoint(18);
			if(etaBin==6) graphMC->RemovePoint(10);
			if(etaBin==6) graphMC->RemovePoint(20);

			if(etaBin==7) graphMC->RemovePoint(15);
			if(etaBin==7) graphMC->RemovePoint(15);
			if(etaBin==7) graphMC->RemovePoint(17);
			if(etaBin==7) graphMC->RemovePoint(19);
			if(etaBin==7) graphMC->RemovePoint(19);
			if(etaBin==7) graphMC->RemovePoint(20);
		}

		if(etaBin==6) graphMC->RemovePoint(4);
		if(etaBin==6) graphMC->RemovePoint(4);
		//	        if(etaBin==6) graphMC->RemovePoint(1);
		//	        if(etaBin==6) graphMC->RemovePoint(1);
		//	        if(etaBin==6) graphMC->RemovePoint(0);
		//	        if(etaBin==7) graphMC->RemovePoint(7);

		if(etaBin==7) graphMC->RemovePoint(5);
		if(etaBin==7) graphMC->RemovePoint(6);

		sprintf(graphName,"gEff_DATA_PT_AETA%d",etaBin);
		TGraphAsymmErrors *graphDATA = new TGraphAsymmErrors(*((TGraphAsymmErrors *) fInEffData->Get(graphName)));

		graphDATA->SaveAs("tmp/graphDataBefore.root");

		///// Alter according to MuonID, May19:
		/*
			 double SFerror;
			 double SF;

			 if(etaBin==0||etaBin==2||etaBin==3||etaBin==4||etaBin==5) { SF=1.019; SFerror=0.005; }//for the bin |eta| < 1.2 (excl. 0.2-0.3)
			 if(etaBin==1) { SF=1.032; SFerror=0.01; }//for the bin 0.2 < |eta| < 0.3
			 if(etaBin==6||etaBin==7) { SF=1.019; SFerror=0.01; }//for the bin 1.2 < |eta| < 1.6



		//	  SF=1;
		//	  SFerror=0;

		//	  SF = 1.019 +- 0.005 for the bin |eta| < 1.2 (excl. 0.2-0.3) (0.5% relative error)
		//	  SF = 1.019 +- 0.010 for the bin 1.2 < |eta| < 1.6 (1% relative error)
		//	  SF = 1.032 +- 0.010 for the bin 0.2 < |eta| < 0.3 (1% relative error)

		int nPtALTER_=graphDATA->GetN();
		const int nDimALTER=nPtALTER_;
		double effDATAALTER[nDimALTER];
		double pTDATAALTER[nDimALTER];
		double err_pTDATAALTER_low[nDimALTER];
		double err_pTDATAALTER_high[nDimALTER];
		double efferrALTER_low[nDimALTER];
		double efferrALTER_high[nDimALTER];


		for(int i=0;i<nPtALTER_;i++){

		graphDATA->GetPoint(i,pTDATAALTER[i],effDATAALTER[i]);
		graphDATA->SetPoint(i,pTDATAALTER[i],effDATAALTER[i]*SF);

		err_pTDATAALTER_low[i]=graphDATA->GetErrorXlow(i);
		err_pTDATAALTER_high[i]=graphDATA->GetErrorXhigh(i);
		efferrALTER_low[i]=graphDATA->GetErrorYlow(i);
		efferrALTER_high[i]=graphDATA->GetErrorYhigh(i);

		efferrALTER_low[i]=TMath::Sqrt(TMath::Power(efferrALTER_low[i],2)+TMath::Power(SFerror,2));
		efferrALTER_high[i]=TMath::Sqrt(TMath::Power(efferrALTER_high[i],2)+TMath::Power(SFerror,2));

		graphDATA->SetPointError(i,err_pTDATAALTER_low[i],err_pTDATAALTER_high[i],efferrALTER_low[i],efferrALTER_high[i]);


		}

		//	  graphDATA = new TGraphAsymmErrors(nPtALTER_,pTDATAALTER,effDATAALTER,err_pTDATAALTER_low,err_pTDATAALTER_high,efferrALTER_low,efferrALTER_high);

		graphDATA->SaveAs("tmp/graphDataAfter.root");
		*/
		/////// End Altering Data Graph

		double pTdist = 150./double(pTBinsNew);


		///////// 'Fit' Settings ////////////////////////

		double pTshiftEst;
		double pTscaleEst;
		double effscaleEst;
		double effscale2Est;
		double effshiftEst;


		///////////////////////////////////////////////////
		/////////////// DO 'FIT' //////////////////////////
		///////////////////////////////////////////////////


		//pars[0]=pTshift
		//pars[1]=pTscale
		//pars[2]=effshift
		//pars[3]=effscale


		double amin, edm, errdef;
		int nvpar, nparx;
		bool fit_effscale=false;
		bool fit_effscale2=false;
		bool fit_pTscale=true;

		//	   minuitx->SetParameter(0,"name",val,firststep,startposition,wtf?);

		FcnPVLogL* fcnx = new FcnPVLogL(graphMC,graphDATA,pTdist,pTBinsNew);
		TFitterMinuit* minuitx = new TFitterMinuit();
		minuitx->SetMinuitFCN(fcnx);
		/*	   minuitx->SetParameter(0,"pTshift",0,1,-1,1);
					 minuitx->SetParameter(1,"pTscale",1,0.4,0.6,1.4);
					 minuitx->SetParameter(2,"effshift",0,0.4,-0.4,0.4);
					 minuitx->SetParameter(3,"effscale",1,0.4,0.6,1.4);
					 */
		minuitx->SetParameter(0,"pTshift",0,1e-4,-3.5,3.5);
		minuitx->SetParameter(1,"pTscale",1,1e-2,0.4,1.6);
		minuitx->SetParameter(2,"effshift",0,1e-4,-0.1,0.1);
		if(fit_effscale||fit_effscale2)
			minuitx->SetParameter(2,"effshift",0,1e-4,-0.6,0.6);
		minuitx->SetParameter(3,"effscale",0,1e-4,-0.0035,0.0035);
		minuitx->SetParameter(4,"effscale2",1,1e-6,0.4,1.6);

		if(!fit_pTscale) minuitx->FixParameter(1);
		if(!fit_effscale) minuitx->FixParameter(3);
		if(!fit_effscale2) minuitx->FixParameter(4);

		//	   minuitx->FixParameter(0);
		//	   minuitx->FixParameter(1);
		//	   minuitx->FixParameter(2);
		//	   minuitx->FixParameter(3);
		//	   minuitx->FixParameter(4);
		//	   minuitx->SetStrategy(2);
		//	   minuitx->SetMaxIterations(1000);
		//	   minuitx->SetMinimumTolerance(1e2);
		minuitx->SetPrintLevel(3);

		minuitx->CreateMinimizer();
		int nFits=10;
		if(etaBin==6) nFits=1;
		if(etaBin==7) nFits=1;
		for(int iFits=0;iFits<nFits;iFits++){
			minuitx->Minimize();
			minuitx->PrintResults(3,1.);
			minuitx->GetStats(amin, edm, errdef, nvpar, nparx);
			cout<<"GetStats"<<endl;
			cout<<"amin = "<<amin<<endl;
			cout<<"edm = "<<edm<<endl;
			cout<<"errdef = "<<errdef<<endl;
			cout<<"nvpar = "<<nvpar<<endl;
			cout<<"nparx = "<<nparx<<endl;
			//		   graphDATA->RemovePoint(0);
			if(edm<1e-8) break;
		}

		std::vector<double> myPar(nparx);

		Double_t* arglist;


		//	   minuitx->SetParameter(0,"pTshift",-0.5995,1e-4,-1.5,1.5);
		//	   minuitx->FixParameter(0);
		//	   minuitx->FixParameter(2);

		bool exeMINOS=false;
		if(exeMINOS){
			minuitx->ExecuteCommand("MINOS",arglist,nvpar);
			minuitx->PrintResults(3,1.);
		}

		//	   minuitx->ReleaseParameter(0);

		/*	   minuitx->FixParameter(0);
					 minuitx->FixParameter(2);
					 minuitx->Minimize();
					 minuitx->ReleaseParameter(0);
					 minuitx->FixParameter(1);
					 minuitx->Minimize();
					 minuitx->ReleaseParameter(2);
					 minuitx->FixParameter(0);
					 minuitx->Minimize();
					 minuitx->ReleaseParameter(1);
					 minuitx->ReleaseParameter(0);
					 */
		edmEta[etaBin]=edm;


		pTshiftEst=minuitx->GetParameter(0);
		pTscaleEst=minuitx->GetParameter(1);
		effshiftEst=minuitx->GetParameter(2);
		effscaleEst=minuitx->GetParameter(3);
		effscale2Est=minuitx->GetParameter(4);


		double err_pTshiftEst=minuitx->GetParError(0);
		double err_pTscaleEst=minuitx->GetParError(1);
		double err_effshiftEst=minuitx->GetParError(2);
		double err_effscaleEst=minuitx->GetParError(3);
		double err_effscale2Est=minuitx->GetParError(4);

		double eplus_pTshiftEst, eminus_pTshiftEst, eparab_pTshiftEst, globcc_pTshiftEst;
		double eplus_pTscaleEst, eminus_pTscaleEst, eparab_pTscaleEst, globcc_pTscaleEst;
		double eplus_effshiftEst, eminus_effshiftEst, eparab_effshiftEst, globcc_effshiftEst;
		double eplus_effscaleEst, eminus_effscaleEst, eparab_effscaleEst, globcc_effscaleEst;
		minuitx->GetErrors(0, eplus_pTshiftEst, eminus_pTshiftEst, eparab_pTshiftEst, globcc_pTshiftEst);
		minuitx->GetErrors(1, eplus_pTscaleEst, eminus_pTscaleEst, eparab_pTscaleEst, globcc_pTscaleEst);
		minuitx->GetErrors(2, eplus_effshiftEst, eminus_effshiftEst, eparab_effshiftEst, globcc_effshiftEst);
		minuitx->GetErrors(3, eplus_effscaleEst, eminus_effscaleEst, eparab_effscaleEst, globcc_effscaleEst);
		//	   cout<<"pTshift errors : eplus = "<<eplus_pTshiftEst<<" , eminus = "<<eminus_pTshiftEst<<" , eparab = "<<eparab_pTshiftEst<<" , globcc = "<<globcc_pTshiftEst<<" , fiterr = "<<err_pTshiftEst<<endl;
		//	   cout<<"pTscale errors : eplus = "<<eplus_pTscaleEst<<" , eminus = "<<eminus_pTscaleEst<<" , eparab = "<<eparab_pTscaleEst<<" , globcc = "<<globcc_pTscaleEst<<" , fiterr = "<<err_pTscaleEst<<endl;
		//	   cout<<"effshift errors: eplus = "<<eplus_effshiftEst<<" , eminus = "<<eminus_effshiftEst<<" , eparab = "<<eparab_effshiftEst<<" , globcc = "<<globcc_effshiftEst<<" , fiterr = "<<err_effshiftEst<<endl;
		//	   cout<<"effscale errors: eplus = "<<eplus_effscaleEst<<" , eminus = "<<eminus_effscaleEst<<" , eparab = "<<eparab_effscaleEst<<" , globcc = "<<globcc_effscaleEst<<" , fiterr = "<<err_effscaleEst<<endl;

		// Calc Correlation variations

		bool FixedPar[nparx];
		for(int i=0;i<nparx;i++){
			FixedPar[i] = minuitx->IsFixed(i);
		}

		int jump_i=0;
		int jump_j=0;

		TMatrixDSym C = TMatrixDSym(nvpar);
		for(int i=0;i<nvpar;i++){
			for(int j=0;j<nvpar;j++){
				C[i][j] = minuitx->GetCovarianceMatrixElement(i,j);
			}
		}

		cout<< "Covariance Matrix:"<<endl;
		C.Print();

		TMatrixD mu = TMatrixD(nvpar,1);

		jump_i=0;
		for(int i=0;i<nvpar;i++){
			if(FixedPar[i]) jump_i++;
			mu[i][0]=minuitx->GetParameter(i+jump_i);
		}

		cout<< "Estimated parameter vector:"<<endl;
		mu.Print();

		TMatrixDSymEigen eigensys = TMatrixDSymEigen(C);
		TMatrixD O = TMatrixD(eigensys.GetEigenVectors());
		TMatrixD OT = TMatrixD(TMatrixD::kTransposed, O);

		cout<< "Matrix of eigenvectors (covariance matrix):"<<endl;
		O.Print();

		cout<< "Matrix of eigenvectors, transposed"<<endl;
		OT.Print();

		TMatrixDSym Lambda = TMatrixDSym(nvpar);
		for(int i=0;i<nvpar;i++){
			Lambda[i][i] = eigensys.GetEigenValues()[i];
		}

		cout<< "Eigenvalues:"<<endl;
		Lambda.Print();

		cout<< "Should be zero:"<<endl;
		(C - TMatrixD(TMatrixD(O, TMatrixD::kMult, Lambda), TMatrixD::kMult, OT)).Print();

		TMatrixD LM = TMatrixD(nvpar,1);//Shift in Eigenbasis
		TMatrixD shift = TMatrixD(nvpar,1);//shift in Parameterbasis
		TMatrixD varPar = TMatrixD(nparx,nvpar);//shift in Parameterbasis

		jump_j=0;
		for(int i=0;i<nvpar;i++){
			for(int k=0;k<nvpar;k++){
				LM[k][0]=0;
			}
			LM[i][0] = TMath::Sqrt(Lambda[i][i]);
			shift = TMatrixD(O,TMatrixD::kMult, LM);
			jump_j=0;
			for(int j=0;j<nparx;j++){
				if(FixedPar[j]) {jump_j++;varPar[j][i]=0;}
				else {varPar[j][i]=shift[j-jump_j][0];}
			}
			cout<< "Shift vector "<<i<<":"<<endl;
			shift.Print();

		}


		cout<< "shift Matrix (rows=free parameters):"<<endl;
		varPar.Print();


		// check and plot chi2 function

		int chi2CheckBins=200;
		int nSigma;
		bool plot1Dchi2=true;

		nSigma=2;
		double delta=nSigma*err_pTshiftEst;
		TH1D*  chi2_par0= new TH1D( "chi2_par0", "chi2_par0", chi2CheckBins,  pTshiftEst-delta, pTshiftEst+delta);;
		if(plot1Dchi2){
			for(int i=1;i<chi2CheckBins+1;i++){
				myPar[0] = pTshiftEst-delta+delta*i*2/chi2CheckBins;
				myPar[1] = pTscaleEst;
				myPar[2] = effshiftEst;
				myPar[3] = effscaleEst;
				myPar[4] = effscale2Est;
				double chi2 = (*fcnx)(myPar);
				chi2_par0->SetBinContent(i,chi2);
			}}

		delta=nSigma*err_pTscaleEst;
		TH1D*  chi2_par1= new TH1D( "chi2_par1", "chi2_par1", chi2CheckBins,  pTscaleEst-delta, pTscaleEst+delta);;
		if(plot1Dchi2&&fit_pTscale){
			for(int i=1;i<chi2CheckBins+1;i++){
				myPar[0] = pTshiftEst;
				myPar[1] = pTscaleEst-delta+delta/chi2CheckBins*i*2;
				myPar[2] = effshiftEst;
				myPar[3] = effscaleEst;
				myPar[4] = effscale2Est;
				double chi2 = (*fcnx)(myPar);
				chi2_par1->SetBinContent(i,chi2);
			}}

		TH1D*  chi2_par2= new TH1D( "chi2_par2", "chi2_par2", chi2CheckBins,  effshiftEst-nSigma*err_effshiftEst, effshiftEst+nSigma*err_effshiftEst);;
		if(plot1Dchi2){
			for(int i=1;i<chi2CheckBins+1;i++){
				myPar[0] = pTshiftEst;
				myPar[1] = pTscaleEst;
				myPar[2] = effshiftEst-nSigma*err_effshiftEst+err_effshiftEst/chi2CheckBins*i*nSigma*2;
				myPar[3] = effscaleEst;
				myPar[4] = effscale2Est;
				double chi2 = (*fcnx)(myPar);
				chi2_par2->SetBinContent(i,chi2);
			}}

		TH1D*  chi2_par3= new TH1D( "chi2_par3", "chi2_par3", chi2CheckBins,  effscaleEst-nSigma*err_effscaleEst, effscaleEst+nSigma*err_effscaleEst);;
		if(plot1Dchi2&&fit_effscale){
			for(int i=1;i<chi2CheckBins+1;i++){
				myPar[0] = pTshiftEst;
				myPar[1] = pTscaleEst;
				myPar[2] = effshiftEst;
				myPar[3] = effscaleEst-nSigma*err_effscaleEst+err_effscaleEst/chi2CheckBins*i*nSigma*2;
				myPar[4] = effscale2Est;
				double chi2 = (*fcnx)(myPar);
				chi2_par3->SetBinContent(i,chi2);
			}}

		TH1D*  chi2_par4= new TH1D( "chi2_par4", "chi2_par4", chi2CheckBins,  effscale2Est-nSigma*err_effscale2Est, effscale2Est+nSigma*err_effscale2Est);;
		if(plot1Dchi2&&fit_effscale2){
			for(int i=1;i<chi2CheckBins+1;i++){
				myPar[0] = pTshiftEst;
				myPar[1] = pTscaleEst;
				myPar[2] = effshiftEst;
				myPar[3] = effscaleEst;
				myPar[4] = effscale2Est-nSigma*err_effscale2Est+err_effscale2Est/chi2CheckBins*i*nSigma*2;
				double chi2 = (*fcnx)(myPar);
				chi2_par4->SetBinContent(i,chi2);
			}}


		///////////////////////////
		bool plot2Dchi2=false;

		int chi2CheckBins2D=100;
		nSigma=2;
		TH2D* chi2_par0par1_2D   = new TH2D( "chi2_par0par1_2D", "chi2_par0par1_2D", chi2CheckBins2D,  pTscaleEst-nSigma*err_pTscaleEst, pTscaleEst+nSigma*err_pTscaleEst, chi2CheckBins2D,  pTshiftEst-nSigma*err_pTshiftEst, pTshiftEst+nSigma*err_pTshiftEst);
		if(plot2Dchi2&&fit_pTscale){
			for(int i=1;i<chi2CheckBins2D+1;i++){
				for(int j=1;j<chi2CheckBins2D+1;j++){
					myPar[0] = pTshiftEst-nSigma*err_pTshiftEst+err_pTshiftEst/chi2CheckBins2D*i*nSigma*2;
					myPar[1] = pTscaleEst-nSigma*err_pTscaleEst+err_pTscaleEst/chi2CheckBins2D*j*nSigma*2;
					myPar[2] = effshiftEst;
					myPar[3] = effscaleEst;
					myPar[4] = effscale2Est;
					double chi2 = (*fcnx)(myPar);
					if(chi2<10*amin)chi2_par0par1_2D->SetBinContent(j,i,chi2);
				}}
			cout<<"chi2_par0par1_2D done..."<<endl;
		}

		TH2D* chi2_par0par2_2D   = new TH2D( "chi2_par0par2_2D", "chi2_par0par2_2D", chi2CheckBins2D,  effshiftEst-nSigma*err_effshiftEst, effshiftEst+nSigma*err_effshiftEst, chi2CheckBins2D,  pTshiftEst-nSigma*err_pTshiftEst, pTshiftEst+nSigma*err_pTshiftEst);
		if(plot2Dchi2&&fit_pTscale){
			for(int i=1;i<chi2CheckBins2D+1;i++){
				for(int j=1;j<chi2CheckBins2D+1;j++){
					myPar[0] = pTshiftEst-nSigma*err_pTshiftEst+err_pTshiftEst/chi2CheckBins2D*i*nSigma*2;
					myPar[1] = pTscaleEst;
					myPar[2] = effshiftEst-nSigma*err_effshiftEst+err_effshiftEst/chi2CheckBins2D*j*nSigma*2;
					myPar[3] = effscaleEst;
					myPar[4] = effscale2Est;
					double chi2 = (*fcnx)(myPar);
					if(chi2<10*amin)chi2_par0par2_2D->SetBinContent(j,i,chi2);
				}}
			cout<<"chi2_par0par2_2D done..."<<endl;
		}

		TH2D* chi2_par2par1_2D   = new TH2D( "chi2_par2par1_2D", "chi2_par2par1_2D", chi2CheckBins2D,  pTscaleEst-nSigma*err_pTscaleEst, pTscaleEst+nSigma*err_pTscaleEst, chi2CheckBins2D,  effshiftEst-nSigma*err_effshiftEst, effshiftEst+nSigma*err_effshiftEst);
		if(plot2Dchi2&&fit_pTscale){
			for(int i=1;i<chi2CheckBins2D+1;i++){
				for(int j=1;j<chi2CheckBins2D+1;j++){
					myPar[0] = pTshiftEst;
					myPar[1] = pTscaleEst-nSigma*err_pTscaleEst+err_pTscaleEst/chi2CheckBins2D*j*nSigma*2;
					myPar[2] = effshiftEst-nSigma*err_effshiftEst+err_effshiftEst/chi2CheckBins2D*i*nSigma*2;
					myPar[3] = effscaleEst;
					myPar[4] = effscale2Est;
					double chi2 = (*fcnx)(myPar);
					if(chi2<10*amin)chi2_par2par1_2D->SetBinContent(j,i,chi2);
				}}
			cout<<"chi2_par2par1_2D done..."<<endl;
		}

		TH2D* chi2_par2par4_2D   = new TH2D( "chi2_par2par4_2D", "chi2_par2par4_2D", chi2CheckBins2D,  effscale2Est-nSigma*err_effscale2Est, effscale2Est+nSigma*err_effscale2Est, chi2CheckBins2D,  effshiftEst-nSigma*err_effshiftEst, effshiftEst+nSigma*err_effshiftEst);
		if(plot2Dchi2&&fit_effscale2){
			for(int i=1;i<chi2CheckBins2D+1;i++){
				for(int j=1;j<chi2CheckBins2D+1;j++){
					myPar[0] = pTshiftEst;
					myPar[1] = pTscaleEst;
					myPar[2] = effshiftEst-nSigma*err_effshiftEst+err_effshiftEst/chi2CheckBins2D*i*nSigma*2;
					myPar[3] = effscaleEst;
					myPar[4] = effscale2Est-nSigma*err_effscale2Est+err_effscale2Est/chi2CheckBins2D*j*nSigma*2;
					double chi2 = (*fcnx)(myPar);
					if(chi2<10*amin)chi2_par2par4_2D->SetBinContent(j,i,chi2);
				}}
			cout<<"chi2_par2par4_2D done..."<<endl;
		}

		chi2_par0->SetStats(0);
		chi2_par1->SetStats(0);
		chi2_par2->SetStats(0);
		chi2_par3->SetStats(0);
		chi2_par4->SetStats(0);
		chi2_par0par1_2D->SetStats(0);
		chi2_par0par2_2D->SetStats(0);
		chi2_par2par1_2D->SetStats(0);
		chi2_par2par4_2D->SetStats(0);

		char savename[200];

		TLine* RefLinepTshift=new TLine( pTshiftEst, -1000, pTshiftEst, 1000 );
		RefLinepTshift->SetLineWidth( 2 );
		RefLinepTshift->SetLineStyle( 2 );
		RefLinepTshift->SetLineColor( kGreen+2 );
		TLine* RefLinepTscale=new TLine( pTscaleEst, -1000, pTscaleEst, 1000 );
		RefLinepTscale->SetLineWidth( 2 );
		RefLinepTscale->SetLineStyle( 2 );
		RefLinepTscale->SetLineColor( kGreen+2 );
		TLine* RefLineeffshift=new TLine( effshiftEst, -1000, effshiftEst, 1000 );
		RefLineeffshift->SetLineWidth( 2 );
		RefLineeffshift->SetLineStyle( 2 );
		RefLineeffshift->SetLineColor( kGreen+2 );
		TLine* RefLineeffscale=new TLine( effscaleEst, -1000, effscaleEst, 1000 );
		RefLineeffscale->SetLineWidth( 2 );
		RefLineeffscale->SetLineStyle( 2 );
		RefLineeffscale->SetLineColor( kGreen+2 );
		TLine* RefLineeffscale2=new TLine( effscale2Est, -1000, effscale2Est, 1000 );
		RefLineeffscale2->SetLineWidth( 2 );
		RefLineeffscale2->SetLineStyle( 2 );
		RefLineeffscale2->SetLineColor( kGreen+2 );


		TCanvas *chi2CheckCanvas = new TCanvas("chi2CheckCanvas","chi2CheckCanvas",1000,800);
		chi2CheckCanvas->SetFillColor(kWhite);
		chi2CheckCanvas->Divide(2,3);
		chi2CheckCanvas->GetFrame()->SetFillColor(kWhite);
		chi2CheckCanvas->GetFrame()->SetBorderSize(0);
		chi2CheckCanvas->cd(1); gPad->SetLeftMargin(0.2); gPad->SetFillColor(kWhite);chi2_par0->GetYaxis()->SetTitle("chi^{2}");	chi2_par0->GetXaxis()->SetTitle("pT_{shift}");	chi2_par0->Draw(); RefLinepTshift->Draw( "same" );
		chi2CheckCanvas->cd(2); gPad->SetLeftMargin(0.2); gPad->SetFillColor(kWhite);chi2_par1->GetYaxis()->SetTitle("chi^{2}");	chi2_par1->GetXaxis()->SetTitle("pT_{scale}");	chi2_par1->Draw(); RefLinepTscale->Draw( "same" );
		chi2CheckCanvas->cd(3); gPad->SetLeftMargin(0.2); gPad->SetFillColor(kWhite);chi2_par2->GetYaxis()->SetTitle("chi^{2}");	chi2_par2->GetXaxis()->SetTitle("eff_{shift}");	chi2_par2->Draw(); RefLineeffshift->Draw( "same" );
		if(fit_effscale) {chi2CheckCanvas->cd(4); gPad->SetLeftMargin(0.2); gPad->SetFillColor(kWhite);chi2_par3->GetYaxis()->SetTitle("chi^{2}");	chi2_par3->GetXaxis()->SetTitle("eff_{scale}");	chi2_par3->Draw(); RefLineeffscale->Draw( "same" );}
		if(fit_effscale2) {chi2CheckCanvas->cd(5); gPad->SetLeftMargin(0.2); gPad->SetFillColor(kWhite);chi2_par4->GetYaxis()->SetTitle("chi^{2}");	chi2_par4->GetXaxis()->SetTitle("#tilde{eff}_{scale}");	chi2_par4->Draw(); RefLineeffscale2->Draw( "same" );}
		sprintf(savename,"Scaling/%s/eta%d_chi2Check.pdf",JobID,etaBin);
		chi2CheckCanvas->SaveAs(savename);
		delete chi2CheckCanvas;

		TCanvas *chi2CheckCanvas2D = new TCanvas("chi2CheckCanvas2D","chi2CheckCanvas2D",1000,800);
		chi2CheckCanvas2D->SetFillColor(kWhite);
		chi2CheckCanvas2D->Divide(2,2);
		chi2CheckCanvas2D->GetFrame()->SetFillColor(kWhite);
		chi2CheckCanvas2D->GetFrame()->SetBorderSize(0);
		chi2CheckCanvas2D->cd(1); gPad->SetLeftMargin(0.2); gPad->SetFillColor(kWhite);chi2_par0par1_2D->GetYaxis()->SetTitle("pT_{shift}");	chi2_par0par1_2D->GetXaxis()->SetTitle("pT_{scale}");	chi2_par0par1_2D->Draw("colz");
		chi2CheckCanvas2D->cd(2); gPad->SetLeftMargin(0.2); gPad->SetFillColor(kWhite);chi2_par0par2_2D->GetYaxis()->SetTitle("pT_{shift}");	chi2_par0par2_2D->GetXaxis()->SetTitle("eff_{shift}");	chi2_par0par2_2D->Draw("colz");
		chi2CheckCanvas2D->cd(3); gPad->SetLeftMargin(0.2); gPad->SetFillColor(kWhite);chi2_par2par1_2D->GetYaxis()->SetTitle("eff_{shift}");	chi2_par2par1_2D->GetXaxis()->SetTitle("pT_{scale}");	chi2_par2par1_2D->Draw("colz");
		chi2CheckCanvas2D->cd(4); gPad->SetLeftMargin(0.2); gPad->SetFillColor(kWhite);chi2_par2par4_2D->GetYaxis()->SetTitle("eff_{shift}");	chi2_par2par4_2D->GetXaxis()->SetTitle("#tilde{eff}_{scale}");	chi2_par2par4_2D->Draw("colz");
		sprintf(savename,"Scaling/%s/eta%d_chi2Check2D.pdf",JobID,etaBin);
		chi2CheckCanvas2D->SaveAs(savename);
		delete chi2CheckCanvas2D;

		delete chi2_par0;
		delete chi2_par1;
		delete chi2_par2;
		delete chi2_par3;
		delete chi2_par4;
		delete chi2_par0par1_2D;
		delete chi2_par0par2_2D;
		delete chi2_par2par1_2D;
		delete chi2_par2par4_2D;


		////////////////////// PLOT RESIDUALS //////////////////////////////////////////////////////

		bool plotResiduals=true;

		int chi2CheckBinsResidual=500;
		nSigma=2;
		TH1D* residual_par0   = new TH1D( "residual_par0", "residual_par0", chi2CheckBinsResidual,  pTscaleEst-nSigma*err_pTscaleEst, pTscaleEst+nSigma*err_pTscaleEst);
		TH1D* residual_par0Single1    = new TH1D( "residual_par0Single1 ", "residual_par0Single1 ", chi2CheckBinsResidual,  pTscaleEst-nSigma*err_pTscaleEst, pTscaleEst+nSigma*err_pTscaleEst);
		TH1D* residual_par0Single2    = new TH1D( "residual_par0Single2 ", "residual_par0Single2 ", chi2CheckBinsResidual,  pTscaleEst-nSigma*err_pTscaleEst, pTscaleEst+nSigma*err_pTscaleEst);
		TH1D* residual_par0Single3    = new TH1D( "residual_par0Single3 ", "residual_par0Single3 ", chi2CheckBinsResidual,  pTscaleEst-nSigma*err_pTscaleEst, pTscaleEst+nSigma*err_pTscaleEst);
		TH1D* residual_par0Single4    = new TH1D( "residual_par0Single4 ", "residual_par0Single4 ", chi2CheckBinsResidual,  pTscaleEst-nSigma*err_pTscaleEst, pTscaleEst+nSigma*err_pTscaleEst);
		TH1D* residual_par0Single5    = new TH1D( "residual_par0Single5 ", "residual_par0Single5 ", chi2CheckBinsResidual,  pTscaleEst-nSigma*err_pTscaleEst, pTscaleEst+nSigma*err_pTscaleEst);
		TH1D* residual_par0Single6    = new TH1D( "residual_par0Single6 ", "residual_par0Single6 ", chi2CheckBinsResidual,  pTscaleEst-nSigma*err_pTscaleEst, pTscaleEst+nSigma*err_pTscaleEst);
		TH1D* residual_par0Single7    = new TH1D( "residual_par0Single7 ", "residual_par0Single7 ", chi2CheckBinsResidual,  pTscaleEst-nSigma*err_pTscaleEst, pTscaleEst+nSigma*err_pTscaleEst);
		TH1D* residual_par0Single8    = new TH1D( "residual_par0Single8 ", "residual_par0Single8 ", chi2CheckBinsResidual,  pTscaleEst-nSigma*err_pTscaleEst, pTscaleEst+nSigma*err_pTscaleEst);
		TH1D* residual_par0Single9    = new TH1D( "residual_par0Single9 ", "residual_par0Single9 ", chi2CheckBinsResidual,  pTscaleEst-nSigma*err_pTscaleEst, pTscaleEst+nSigma*err_pTscaleEst);
		TH1D* residual_par0Single10   = new TH1D( "residual_par0Single10", "residual_par0Single10", chi2CheckBinsResidual,  pTscaleEst-nSigma*err_pTscaleEst, pTscaleEst+nSigma*err_pTscaleEst);
		TH1D* residual_par0Single11   = new TH1D( "residual_par0Single11", "residual_par0Single11", chi2CheckBinsResidual,  pTscaleEst-nSigma*err_pTscaleEst, pTscaleEst+nSigma*err_pTscaleEst);
		TH1D* residual_par0Single12   = new TH1D( "residual_par0Single12", "residual_par0Single12", chi2CheckBinsResidual,  pTscaleEst-nSigma*err_pTscaleEst, pTscaleEst+nSigma*err_pTscaleEst);
		TH1D* residual_par0Single13   = new TH1D( "residual_par0Single13", "residual_par0Single13", chi2CheckBinsResidual,  pTscaleEst-nSigma*err_pTscaleEst, pTscaleEst+nSigma*err_pTscaleEst);

		if(plotResiduals){

			for(int j=1;j<chi2CheckBinsResidual+1;j++){

				myPar[0] = pTshiftEst;
				myPar[1] = pTscaleEst-nSigma*err_pTscaleEst+err_pTscaleEst/chi2CheckBinsResidual*j*nSigma*2;
				myPar[2] = effshiftEst;
				myPar[3] = effscaleEst;
				myPar[4] = effscale2Est;

				double ptCentre_res[pTBinsNew];
				double effMean_res[pTBinsNew];
				double residual=0;
				double residualSingle[13];

				for(int pTBin=0;pTBin<pTBinsNew;pTBin++){
					ptCentre_res[pTBin]=pTBin*pTdist;
					effMean_res[pTBin]=graphMC->Eval(pTBin*pTdist/myPar[1]-myPar[0])*(1+pTBin*pTdist*myPar[3])*myPar[4]+myPar[2];
				}

				TGraphAsymmErrors *Ngraph = new TGraphAsymmErrors(pTBinsNew,ptCentre_res,effMean_res,0,0,0,0);

				char FitOptions[200];
				sprintf(FitOptions,"EFNRQ");
				TF1* f1local;
				int FitRange=2;

				// Calculate chi2
				int nPt_=graphDATA->GetN();
				const int nDim=nPt_;
				double effDATA[nDim];
				double pTDATA[nDim];
				double efferr[nDim];
				double effModel[nDim];
				double chi2[nDim];
				double chiSq=0;

				for(int i=0;i<nPt_;i++){
					graphDATA->GetPoint(i,pTDATA[i],effDATA[i]);
					//      f1local = new TF1("f1local","pol1",pTDATA[i]-FitRange*graphDATA->GetErrorXlow(i),pTDATA[i]+FitRange*graphDATA->GetErrorXhigh(i));
					//      Ngraph->Fit("f1local",FitOptions);

					//	  effModel[i] = f1local->GetParameter(0) + f1local->GetParameter(1)*pTDATA[i];// + f1local->GetParameter(2)*pTDATA[i]*pTDATA[i];
					effModel[i] = Ngraph->Eval(pTDATA[i]);//(Ngraph->Eval(pTDATA[i])+Ngraph->Eval(pTDATA[i]+pTdist)+Ngraph->Eval(pTDATA[i]-pTdist))/3;
					if(effDATA[i]-effModel[i]>0)efferr[i]=graphDATA->GetErrorYlow(i);
					else efferr[i]=graphDATA->GetErrorYhigh(i);
					chi2[i]=(effModel[i]-effDATA[i])*(effModel[i]-effDATA[i])/(efferr[i]*efferr[i]);
					chiSq=chiSq+chi2[i];
					delete f1local;
					residual+=(effDATA[i]-effModel[i]);
					residualSingle[i]=(effDATA[i]-effModel[i]);///(efferr[i]*efferr[i]);
					if(i==0 )	residual_par0Single1 	->SetBinContent(j,residualSingle[i]);
					if(i==1 )	residual_par0Single2 	->SetBinContent(j,residualSingle[i]);
					if(i==2 )	residual_par0Single3 	->SetBinContent(j,residualSingle[i]);
					if(i==3 )	residual_par0Single4 	->SetBinContent(j,residualSingle[i]);
					if(i==4 )	residual_par0Single5 	->SetBinContent(j,residualSingle[i]);
					if(i==5 )	residual_par0Single6 	->SetBinContent(j,residualSingle[i]);
					if(i==6 )	residual_par0Single7 	->SetBinContent(j,residualSingle[i]);
					if(i==7 )	residual_par0Single8 	->SetBinContent(j,residualSingle[i]);
					if(i==8 )	residual_par0Single9 	->SetBinContent(j,residualSingle[i]);
					if(i==9 )	residual_par0Single10	->SetBinContent(j,residualSingle[i]);
					if(i==10)	residual_par0Single11	->SetBinContent(j,residualSingle[i]);
					if(i==11)	residual_par0Single12	->SetBinContent(j,residualSingle[i]);
					if(i==12)	residual_par0Single13	->SetBinContent(j,residualSingle[i]);
				}



				residual_par0->SetBinContent(j,residual);

			}

		}

		residual_par0->SetStats(0);
		residual_par0Single1 ->SetStats(0);
		residual_par0Single2 ->SetStats(0);
		residual_par0Single3 ->SetStats(0);
		residual_par0Single4 ->SetStats(0);
		residual_par0Single5 ->SetStats(0);
		residual_par0Single6 ->SetStats(0);
		residual_par0Single7 ->SetStats(0);
		residual_par0Single8 ->SetStats(0);
		residual_par0Single9 ->SetStats(0);
		residual_par0Single10->SetStats(0);
		residual_par0Single11->SetStats(0);
		residual_par0Single12->SetStats(0);
		residual_par0Single13->SetStats(0);

		TCanvas *residualCheckCanvas = new TCanvas("residualCheckCanvas","residualCheckCanvas",1000,800);
		residualCheckCanvas->SetFillColor(kWhite);
		residualCheckCanvas->Divide(5,4);
		residualCheckCanvas->GetFrame()->SetFillColor(kWhite);
		residualCheckCanvas->GetFrame()->SetBorderSize(0);
		residualCheckCanvas->cd(1 ); gPad->SetLeftMargin(0.2); gPad->SetFillColor(kWhite);residual_par0->GetYaxis()->SetTitle("residual");	residual_par0->GetXaxis()->SetTitle("pT_{shift}");	residual_par0->Draw();
		residualCheckCanvas->cd(2 ); gPad->SetLeftMargin(0.2); gPad->SetFillColor(kWhite);residual_par0Single1 ->GetYaxis()->SetTitle("residual1 ");	residual_par0Single1 ->GetXaxis()->SetTitle("pT_{scale}");	residual_par0Single1 ->Draw(); RefLinepTscale->Draw( "same" );
		residualCheckCanvas->cd(3 ); gPad->SetLeftMargin(0.2); gPad->SetFillColor(kWhite);residual_par0Single2 ->GetYaxis()->SetTitle("residual2 ");	residual_par0Single2 ->GetXaxis()->SetTitle("pT_{scale}");	residual_par0Single2 ->Draw(); RefLinepTscale->Draw( "same" );
		residualCheckCanvas->cd(4 ); gPad->SetLeftMargin(0.2); gPad->SetFillColor(kWhite);residual_par0Single3 ->GetYaxis()->SetTitle("residual3 ");	residual_par0Single3 ->GetXaxis()->SetTitle("pT_{scale}");	residual_par0Single3 ->Draw(); RefLinepTscale->Draw( "same" );
		residualCheckCanvas->cd(5 ); gPad->SetLeftMargin(0.2); gPad->SetFillColor(kWhite);residual_par0Single4 ->GetYaxis()->SetTitle("residual4 ");	residual_par0Single4 ->GetXaxis()->SetTitle("pT_{scale}");	residual_par0Single4 ->Draw(); RefLinepTscale->Draw( "same" );
		residualCheckCanvas->cd(6 ); gPad->SetLeftMargin(0.2); gPad->SetFillColor(kWhite);residual_par0Single5 ->GetYaxis()->SetTitle("residual5 ");	residual_par0Single5 ->GetXaxis()->SetTitle("pT_{scale}");	residual_par0Single5 ->Draw(); RefLinepTscale->Draw( "same" );
		residualCheckCanvas->cd(7 ); gPad->SetLeftMargin(0.2); gPad->SetFillColor(kWhite);residual_par0Single6 ->GetYaxis()->SetTitle("residual6 ");	residual_par0Single6 ->GetXaxis()->SetTitle("pT_{scale}");	residual_par0Single6 ->Draw(); RefLinepTscale->Draw( "same" );
		residualCheckCanvas->cd(8 ); gPad->SetLeftMargin(0.2); gPad->SetFillColor(kWhite);residual_par0Single7 ->GetYaxis()->SetTitle("residual7 ");	residual_par0Single7 ->GetXaxis()->SetTitle("pT_{scale}");	residual_par0Single7 ->Draw(); RefLinepTscale->Draw( "same" );
		residualCheckCanvas->cd(9 ); gPad->SetLeftMargin(0.2); gPad->SetFillColor(kWhite);residual_par0Single8 ->GetYaxis()->SetTitle("residual8 ");	residual_par0Single8 ->GetXaxis()->SetTitle("pT_{scale}");	residual_par0Single8 ->Draw(); RefLinepTscale->Draw( "same" );
		residualCheckCanvas->cd(10); gPad->SetLeftMargin(0.2); gPad->SetFillColor(kWhite);residual_par0Single9 ->GetYaxis()->SetTitle("residual9 ");	residual_par0Single9 ->GetXaxis()->SetTitle("pT_{scale}");	residual_par0Single9 ->Draw(); RefLinepTscale->Draw( "same" );
		residualCheckCanvas->cd(11); gPad->SetLeftMargin(0.2); gPad->SetFillColor(kWhite);residual_par0Single10->GetYaxis()->SetTitle("residual10");	residual_par0Single10->GetXaxis()->SetTitle("pT_{scale}");	residual_par0Single10->Draw(); RefLinepTscale->Draw( "same" );
		residualCheckCanvas->cd(12); gPad->SetLeftMargin(0.2); gPad->SetFillColor(kWhite);residual_par0Single11->GetYaxis()->SetTitle("residual11");	residual_par0Single11->GetXaxis()->SetTitle("pT_{scale}");	residual_par0Single11->Draw(); RefLinepTscale->Draw( "same" );
		residualCheckCanvas->cd(13); gPad->SetLeftMargin(0.2); gPad->SetFillColor(kWhite);residual_par0Single12->GetYaxis()->SetTitle("residual12");	residual_par0Single12->GetXaxis()->SetTitle("pT_{scale}");	residual_par0Single12->Draw(); RefLinepTscale->Draw( "same" );
		residualCheckCanvas->cd(14); gPad->SetLeftMargin(0.2); gPad->SetFillColor(kWhite);residual_par0Single13->GetYaxis()->SetTitle("residual13");	residual_par0Single13->GetXaxis()->SetTitle("pT_{scale}");	residual_par0Single13->Draw(); RefLinepTscale->Draw( "same" );
		sprintf(savename,"Scaling/%s/eta%d_residualCheck.pdf",JobID,etaBin);
		residualCheckCanvas->SaveAs(savename);
		delete residualCheckCanvas;
		delete residual_par0;
		delete residual_par0Single1 ;
		delete residual_par0Single2 ;
		delete residual_par0Single3 ;
		delete residual_par0Single4 ;
		delete residual_par0Single5 ;
		delete residual_par0Single6 ;
		delete residual_par0Single7 ;
		delete residual_par0Single8 ;
		delete residual_par0Single9 ;
		delete residual_par0Single10;
		delete residual_par0Single11;
		delete residual_par0Single12;
		delete residual_par0Single13;





		///////////////////////////////////////////////////
		/////////////// END OF DO 'FIT' ///////////////////
		///////////////////////////////////////////////////

		double plotEffmin=0.3;
		double plotpTmax=40;

		///////// DRAW stuff ////////////////////////


		TCanvas *plotCanvas = new TCanvas("plotCanvas","plotCanvas",1000,800);
		plotCanvas->SetFillColor(kWhite);
		plotCanvas->SetGrid();
		plotCanvas->GetFrame()->SetFillColor(kWhite);
		plotCanvas->GetFrame()->SetBorderSize(0);
		plotCanvas->SetRightMargin(0.05) ;

		TH1F *plotHisto = new TH1F;
		plotHisto = plotCanvas->DrawFrame(0,plotEffmin,plotpTmax,1.1);
		plotHisto->SetXTitle("p_{T}(#mu) [GeV/c]");
		plotHisto->GetYaxis()->SetTitleOffset(1.5);

		TLine* RefLine = new TLine( 0, 1,plotpTmax, 1 );
		RefLine->SetLineWidth( 2 );
		RefLine->SetLineStyle( 2 );
		RefLine->SetLineColor( kGreen+2 );
		RefLine->Draw( "same" );

		TLegend* plotLegend=new TLegend(0.55,0.15,0.9375,0.3);
		plotLegend->SetFillColor(kWhite);
		plotLegend->SetTextFont(72);
		plotLegend->SetTextSize(0.02);
		plotLegend->SetBorderSize(1);
		char legendentry[200];

		TBox* box = new TBox( 20.5,0.55 , 40, 0.75);
		box->SetLineColor(kBlack);
		box->SetLineWidth(2);
		box->SetLineStyle(1);
		box->SetFillColor(kWhite);
		//	   		box->Draw( "same" );


		//pars[0]=pTshift
		//pars[1]=pTscale
		//pars[2]=effshift
		//pars[3]=effscale

		double ptCentre_red[pTBinsNew];
		double effMean_red[pTBinsNew];

		double pTshift_Draw;
		double pTscale_Draw;
		double effshift_Draw;
		double effscale_Draw;
		double effscale2_Draw;

		///////////////////////////////////////////////////////////////////////////////////
		pTshift_Draw=minuitx->GetParameter(0);
		pTscale_Draw=minuitx->GetParameter(1);
		effshift_Draw=minuitx->GetParameter(2);
		effscale_Draw=minuitx->GetParameter(3);
		effscale2_Draw=minuitx->GetParameter(4);

		for(int pTBin=0;pTBin<pTBinsNew;pTBin++){
			ptCentre_red[pTBin]=pTBin*pTdist;
			effMean_red[pTBin]=graphMC->Eval(pTBin*pTdist/pTscale_Draw-pTshift_Draw)*(1+pTBin*pTdist*effscale_Draw)*effscale2_Draw+effshift_Draw;
		}

		TGraphAsymmErrors *Ngraph = new TGraphAsymmErrors(pTBinsNew,ptCentre_red,effMean_red,0,0,0,0);

		int nDim=13;
		int nParam=nvpar;
		double chisquare=amin/(nDim-nParam);
		double BIC=chisquare+nParam*TMath::Log(nDim);

		cout<<"chiSq/ndf = "<<chisquare<<endl;
		cout<<"BIC       = "<<BIC<<endl;

		int SystLineStyle=7;
		double SystLineWidth=1;
		double CentralLineWidth=1;

		Ngraph->SetLineColor(kRed);
		Ngraph->SetLineWidth(CentralLineWidth);
		Ngraph->Draw("same");
		graphDATA->SetMarkerSize(0.8);
		graphDATA->SetLineWidth(0.2);
		graphDATA->Draw("P");

		sprintf(legendentry,"Data TnP");
		plotLegend->AddEntry(graphDATA,legendentry,"ple");
		sprintf(legendentry,"Parametrization,  #chi^{2}/ndf = %1.3f, BIC = %1.3f", chisquare, BIC);
		plotLegend->AddEntry(Ngraph,legendentry,"l");
		plotLegend->Draw();

		char texTex[200];
		sprintf(texTex,"Factorized single muon efficiencies, %1.1f < |#eta| < %1.1f",etaRange[etaBin],etaRange[etaBin+1]);
		TLatex *text = new TLatex(4,1.035,texTex);
		text->SetTextSize(0.035);
		text->Draw( "same" );

		if(fit_pTscale) sprintf(texTex,"pT_{scale} = %1.3f#pm%1.3f, pT_{shift} = %1.3f#pm%1.3f GeV",pTscaleEst,err_pTscaleEst,pTshiftEst,err_pTshiftEst);
		if(!fit_pTscale) sprintf(texTex,"pT_{shift} = %1.3f#pm%1.3f GeV",pTshiftEst,err_pTshiftEst);
		TLatex *text2 = new TLatex(21,0.65,texTex);
		text2->SetTextSize(0.0225);
		text2->Draw( "same" );

		if(fit_effscale&&!fit_effscale2) sprintf(texTex,"eff_{scale} = %1.2f#pm%1.2f *10^{-3}, eff_{shift} = %1.3f#pm%1.3f",effscaleEst*1000,err_effscaleEst*1000,effshiftEst,err_effshiftEst);
		if(fit_effscale&&fit_effscale2) sprintf(texTex,"eff_{scale} = %1.2f#pm%1.2f *10^{-3}, #tilde{eff}_{scale} = %1.3f#pm%1.3f, eff_{shift} = %1.3f#pm%1.3f",effscaleEst*1000,err_effscaleEst*1000,effscale2Est,err_effscale2Est,effshiftEst,err_effshiftEst);
		if(fit_effscale2&&!fit_effscale) sprintf(texTex,"#tilde{eff}_{scale} = %1.3f#pm%1.3f, eff_{shift} = %1.3f#pm%1.3f",effscale2Est,err_effscale2Est,effshiftEst,err_effshiftEst);
		if(!fit_effscale&&!fit_effscale2) sprintf(texTex,"eff_{shift} = %1.3f#pm%1.3f",effshiftEst,err_effshiftEst);
		TLatex *text3 = new TLatex(21,0.55,texTex);
		text3->SetTextSize(0.0225);
		if(fit_effscale&&fit_effscale2) text3->SetTextSize(0.015);
		text3->Draw( "same" );

		sprintf(savename,"Scaling/%s/eta%d_noMC.pdf",JobID,etaBin);
		plotCanvas->SaveAs(savename);

		graphMC->SetMarkerSize(0.8);
		graphMC->SetMarkerStyle(20);
		graphMC->SetMarkerColor(kGreen+2);
		graphMC->SetLineWidth(0.2);
		graphMC->Draw("P");
		sprintf(legendentry,"MC Truth");
		plotLegend->AddEntry(graphMC,legendentry,"ple");

		sprintf(savename,"Scaling/%s/eta%d.pdf",JobID,etaBin);
		plotCanvas->SaveAs(savename);
		plotCanvas->Close();
		delete plotCanvas;


		int nPt_=graphDATA->GetN();
		nDim = nPt_;
		double effDATA[nDim];
		double pTDATA[nDim];
		double efferr[nDim];
		double effModel[nDim];
		double pull[nPt_];
		double errXlow[nDim];
		double errXhigh[nDim];
		double errY[nDim];


		TH1D* pull_histo=new TH1D("pull_histo","pull_histo",24,-3,3);

		for(int i=0;i<nPt_;i++){
			graphDATA->GetPoint(i,pTDATA[i],effDATA[i]);
			effModel[i] = Ngraph->Eval(pTDATA[i]);
			if(effDATA[i]-effModel[i]>0)efferr[i]=graphDATA->GetErrorYlow(i);
			else efferr[i]=graphDATA->GetErrorYhigh(i);
			pull[i]=(effModel[i]-effDATA[i])/(efferr[i])*(-1);
			errXlow[i]=graphDATA->GetErrorXlow(i);
			errXhigh[i]=graphDATA->GetErrorXhigh(i);
			errY[i]=0;
			pull_histo->Fill(pull[i]);
		}

		TGraphAsymmErrors *Pullgraph = new TGraphAsymmErrors(nDim,pTDATA,pull,errXlow,errXhigh,errY,errY);

		TCanvas *pullCanvas = new TCanvas("pullCanvas","pullCanvas",1000,800);
		pullCanvas->SetFillColor(kWhite);
		pullCanvas->SetGrid();
		pullCanvas->GetFrame()->SetFillColor(kWhite);
		pullCanvas->GetFrame()->SetBorderSize(0);
		pullCanvas->SetRightMargin(0.05) ;

		TH1F *pullHisto = new TH1F;
		pullHisto = pullCanvas->DrawFrame(0,-4,plotpTmax,4);
		pullHisto->SetXTitle("p_{T}(#mu) [GeV/c]");
		pullHisto->GetYaxis()->SetTitleOffset(1.5);
		pullHisto->SetYTitle("(#epsilon_{data}-#epsilon_{model})/#sigma_{data}");

		TLine* RefLinePull = new TLine( 0, 0,plotpTmax, 0 );
		RefLinePull->SetLineWidth( 2 );
		RefLinePull->SetLineStyle( 2 );
		RefLinePull->SetLineColor( kGreen+2 );
		RefLinePull->Draw( "same" );
		RefLinePull = new TLine( 0, 2,plotpTmax, 2 );
		RefLinePull->SetLineWidth( 2 );
		RefLinePull->SetLineStyle( 2 );
		RefLinePull->SetLineColor( kBlue-6 );
		RefLinePull->Draw( "same" );
		RefLinePull = new TLine( 0, -2,plotpTmax, -2 );
		RefLinePull->SetLineWidth( 2 );
		RefLinePull->SetLineStyle( 2 );
		RefLinePull->SetLineColor( kBlue-6 );
		RefLinePull->Draw( "same" );

		Pullgraph->SetMarkerStyle(20);
		Pullgraph->SetMarkerSize(0.8);
		Pullgraph->Draw("P");

		sprintf(savename,"Scaling/%s/eta%d_Pull.pdf",JobID,etaBin);
		pullCanvas->SaveAs(savename);
		pullCanvas->Close();
		delete pullCanvas;



		pull_histo->Print();


		TCanvas *pullHistoCanvas = new TCanvas("pullHistoCanvas","pullHistoCanvas",1000,800);
		pullHistoCanvas->SetFillColor(kWhite);
		pullHistoCanvas->SetGrid();
		pullHistoCanvas->GetFrame()->SetFillColor(kWhite);
		pullHistoCanvas->GetFrame()->SetBorderSize(0);
		pullHistoCanvas->SetRightMargin(0.05) ;

		double pullHistoPlotMax=1.5*pull_histo->GetMaximum();
		TH1F *pullHistoHisto = new TH1F;
		pullHistoHisto = pullHistoCanvas->DrawFrame(-4,0,4,pullHistoPlotMax);
		pullHistoHisto->SetYTitle("counts");
		pullHistoHisto->GetYaxis()->SetTitleOffset(1.5);
		pullHistoHisto->SetXTitle("(#epsilon_{data}-#epsilon_{model})/#sigma_{data}");

		TLine* RefLinePullHisto = new TLine( 0, 0,0, pullHistoPlotMax );
		RefLinePullHisto->SetLineWidth( 2 );
		RefLinePullHisto->SetLineStyle( 2 );
		RefLinePullHisto->SetLineColor( kGreen+2 );
		RefLinePullHisto->Draw( "same" );
		RefLinePullHisto = new TLine( -1, 0,-1, pullHistoPlotMax );
		RefLinePullHisto->SetLineWidth( 2 );
		RefLinePullHisto->SetLineStyle( 2 );
		RefLinePullHisto->SetLineColor( kBlue-6 );
		RefLinePullHisto->Draw( "same" );
		RefLinePullHisto = new TLine( 1, 0,1, pullHistoPlotMax );
		RefLinePullHisto->SetLineWidth( 2 );
		RefLinePullHisto->SetLineStyle( 2 );
		RefLinePullHisto->SetLineColor( kBlue-6 );
		RefLinePullHisto->Draw( "same" );

		pull_histo->SetMarkerStyle(20);
		pull_histo->SetMarkerSize(0.8);
		pull_histo->Draw("E,same");

		char texTexPull[200];
		sprintf(texTexPull,"#mu = %1.3f +- %1.3f",pull_histo->GetMean(),pull_histo->GetMeanError());
		TLatex *textPull = new TLatex(1.1,pullHistoPlotMax*0.95,texTexPull);
		textPull->SetTextSize(0.035);
		textPull->Draw( "same" );
		sprintf(texTexPull,"R.M.S. = %1.3f +- %1.3f",pull_histo->GetRMS(),pull_histo->GetRMSError());
		TLatex *textPull2 = new TLatex(1.1,pullHistoPlotMax*0.9,texTexPull);
		textPull2->SetTextSize(0.035);
		textPull2->Draw( "same" );

		sprintf(savename,"Scaling/%s/eta%d_PullHisto.pdf",JobID,etaBin);
		pullHistoCanvas->SaveAs(savename);
		pullHistoCanvas->Close();
		delete pullHistoCanvas;




		TCanvas *SystCanvas = new TCanvas("SystCanvas","SystCanvas",1000,800);
		SystCanvas->SetFillColor(kWhite);
		SystCanvas->SetGrid();
		SystCanvas->GetFrame()->SetFillColor(kWhite);
		SystCanvas->GetFrame()->SetBorderSize(0);
		SystCanvas->SetRightMargin(0.05) ;

		TH1F *SystHisto = new TH1F;
		SystHisto = SystCanvas->DrawFrame(0,plotEffmin,plotpTmax,1.1);
		SystHisto->SetXTitle("p_{T}(#mu) [GeV/c]");
		SystHisto->GetYaxis()->SetTitleOffset(1.5);

		RefLine->Draw( "same" );
		text->Draw( "same" );
		text2->Draw( "same" );
		text3->Draw( "same" );

		TLegend* SystLegend=new TLegend(0.55,0.15,0.9375,0.3);
		SystLegend->SetFillColor(kWhite);
		SystLegend->SetTextFont(72);
		SystLegend->SetTextSize(0.02);
		SystLegend->SetBorderSize(1);

		sprintf(legendentry,"Data TnP");
		SystLegend->AddEntry(graphDATA,legendentry,"ple");
		sprintf(legendentry,"Parametrization,  #chi^{2}/ndf = %1.3f, BIC = %1.3f", chisquare, BIC);
		SystLegend->AddEntry(Ngraph,legendentry,"l");

		Ngraph->SetName(graphName);
		Ngraph->Draw("same");



		bool DrawEasySystematicLines=false;

		///////// Create pTshift systematic lines
		if(DrawEasySystematicLines){
			pTshift_Draw=minuitx->GetParameter(0)+minuitx->GetParError(0);
			pTscale_Draw=minuitx->GetParameter(1);
			effshift_Draw=minuitx->GetParameter(2);
			effscale_Draw=minuitx->GetParameter(3);
			effscale2_Draw=minuitx->GetParameter(4);

			for(int pTBin=1;pTBin<pTBinsNew;pTBin++){
				ptCentre_red[pTBin]=pTBin*pTdist;
				effMean_red[pTBin]=graphMC->Eval(pTBin*pTdist/pTscale_Draw-pTshift_Draw)*(1+pTBin*pTdist*effscale_Draw)*effscale2_Draw+effshift_Draw;
			}

			TGraphAsymmErrors *Ngraph_pTshiftPLUS = new TGraphAsymmErrors(pTBinsNew,ptCentre_red,effMean_red,0,0,0,0);

			Ngraph_pTshiftPLUS->SetLineColor(kBlue);
			Ngraph_pTshiftPLUS->SetLineWidth(SystLineWidth);
			Ngraph_pTshiftPLUS->SetLineStyle(SystLineStyle);
			Ngraph_pTshiftPLUS->SetName(graphName);
			Ngraph_pTshiftPLUS->Draw("same");

			pTshift_Draw=minuitx->GetParameter(0)-minuitx->GetParError(0);
			pTscale_Draw=minuitx->GetParameter(1);
			effshift_Draw=minuitx->GetParameter(2);
			effscale_Draw=minuitx->GetParameter(3);
			effscale2_Draw=minuitx->GetParameter(4);

			for(int pTBin=1;pTBin<pTBinsNew;pTBin++){
				ptCentre_red[pTBin]=pTBin*pTdist;
				effMean_red[pTBin]=graphMC->Eval(pTBin*pTdist/pTscale_Draw-pTshift_Draw)*(1+pTBin*pTdist*effscale_Draw)*effscale2_Draw+effshift_Draw;
			}

			TGraphAsymmErrors *Ngraph_pTshiftMINUS = new TGraphAsymmErrors(pTBinsNew,ptCentre_red,effMean_red,0,0,0,0);

			Ngraph_pTshiftMINUS->SetLineColor(kBlue);
			Ngraph_pTshiftMINUS->SetLineWidth(SystLineWidth);
			Ngraph_pTshiftMINUS->SetLineStyle(SystLineStyle);
			Ngraph_pTshiftMINUS->SetName(graphName);
			Ngraph_pTshiftMINUS->Draw("same");

			sprintf(legendentry,"Parametrization, pT_{shift} #pm 1 #sigma");
			SystLegend->AddEntry(Ngraph_pTshiftMINUS,legendentry,"l");

			///////// Create pTscale systematic lines

			pTshift_Draw=minuitx->GetParameter(0);
			pTscale_Draw=minuitx->GetParameter(1)+minuitx->GetParError(1);
			effshift_Draw=minuitx->GetParameter(2);
			effscale_Draw=minuitx->GetParameter(3);
			effscale2_Draw=minuitx->GetParameter(4);

			for(int pTBin=1;pTBin<pTBinsNew;pTBin++){
				ptCentre_red[pTBin]=pTBin*pTdist;
				effMean_red[pTBin]=graphMC->Eval(pTBin*pTdist/pTscale_Draw-pTshift_Draw)*(1+pTBin*pTdist*effscale_Draw)*effscale2_Draw+effshift_Draw;
			}

			TGraphAsymmErrors *Ngraph_pTscalePLUS = new TGraphAsymmErrors(pTBinsNew,ptCentre_red,effMean_red,0,0,0,0);

			Ngraph_pTscalePLUS->SetLineColor(kBlue-6);
			Ngraph_pTscalePLUS->SetLineWidth(SystLineWidth);
			Ngraph_pTscalePLUS->SetLineStyle(SystLineStyle);
			Ngraph_pTscalePLUS->SetName(graphName);
			if(fit_pTscale) Ngraph_pTscalePLUS->Draw("same");

			pTshift_Draw=minuitx->GetParameter(0);
			pTscale_Draw=minuitx->GetParameter(1)-minuitx->GetParError(1);
			effshift_Draw=minuitx->GetParameter(2);
			effscale_Draw=minuitx->GetParameter(3);
			effscale2_Draw=minuitx->GetParameter(4);

			for(int pTBin=1;pTBin<pTBinsNew;pTBin++){
				ptCentre_red[pTBin]=pTBin*pTdist;
				effMean_red[pTBin]=graphMC->Eval(pTBin*pTdist/pTscale_Draw-pTshift_Draw)*(1+pTBin*pTdist*effscale_Draw)*effscale2_Draw+effshift_Draw;
			}

			TGraphAsymmErrors *Ngraph_pTscaleMINUS = new TGraphAsymmErrors(pTBinsNew,ptCentre_red,effMean_red,0,0,0,0);

			Ngraph_pTscaleMINUS->SetLineColor(kBlue-6);
			Ngraph_pTscaleMINUS->SetLineWidth(1);
			Ngraph_pTscaleMINUS->SetLineStyle(SystLineStyle);
			Ngraph_pTscaleMINUS->SetName(graphName);
			if(fit_pTscale) Ngraph_pTscaleMINUS->Draw("same");

			sprintf(legendentry,"Parametrization, pT_{scale} #pm 1 #sigma");
			if(fit_pTscale) SystLegend->AddEntry(Ngraph_pTscaleMINUS,legendentry,"l");



			///////// Create effshift systematic lines

			pTshift_Draw=minuitx->GetParameter(0);
			pTscale_Draw=minuitx->GetParameter(SystLineWidth);
			effshift_Draw=minuitx->GetParameter(2)+minuitx->GetParError(2);
			effscale_Draw=minuitx->GetParameter(3);
			effscale2_Draw=minuitx->GetParameter(4);

			for(int pTBin=1;pTBin<pTBinsNew;pTBin++){
				ptCentre_red[pTBin]=pTBin*pTdist;
				effMean_red[pTBin]=graphMC->Eval(pTBin*pTdist/pTscale_Draw-pTshift_Draw)*(1+pTBin*pTdist*effscale_Draw)*effscale2_Draw+effshift_Draw;
			}

			TGraphAsymmErrors *Ngraph_effshiftPLUS = new TGraphAsymmErrors(pTBinsNew,ptCentre_red,effMean_red,0,0,0,0);

			Ngraph_effshiftPLUS->SetLineColor(kGreen+3);
			Ngraph_effshiftPLUS->SetLineWidth(1);
			Ngraph_effshiftPLUS->SetLineStyle(SystLineStyle);
			Ngraph_effshiftPLUS->SetName(graphName);
			Ngraph_effshiftPLUS->Draw("same");

			pTshift_Draw=minuitx->GetParameter(0);
			pTscale_Draw=minuitx->GetParameter(1);
			effshift_Draw=minuitx->GetParameter(2)-minuitx->GetParError(2);
			effscale_Draw=minuitx->GetParameter(3);
			effscale2_Draw=minuitx->GetParameter(4);

			for(int pTBin=1;pTBin<pTBinsNew;pTBin++){
				ptCentre_red[pTBin]=pTBin*pTdist;
				effMean_red[pTBin]=graphMC->Eval(pTBin*pTdist/pTscale_Draw-pTshift_Draw)*(1+pTBin*pTdist*effscale_Draw)*effscale2_Draw+effshift_Draw;
			}

			TGraphAsymmErrors *Ngraph_effshiftMINUS = new TGraphAsymmErrors(pTBinsNew,ptCentre_red,effMean_red,0,0,0,0);

			Ngraph_effshiftMINUS->SetLineColor(kGreen+3);
			Ngraph_effshiftMINUS->SetLineWidth(SystLineWidth);
			Ngraph_effshiftMINUS->SetLineStyle(SystLineStyle);
			Ngraph_effshiftMINUS->SetName(graphName);
			Ngraph_effshiftMINUS->Draw("same");

			sprintf(legendentry,"Parametrization, eff_{shift} #pm 1 #sigma");
			SystLegend->AddEntry(Ngraph_effshiftMINUS,legendentry,"l");


			///////// Create effscale systematic lines

			pTshift_Draw=minuitx->GetParameter(0);
			pTscale_Draw=minuitx->GetParameter(1);
			effshift_Draw=minuitx->GetParameter(2);
			effscale_Draw=minuitx->GetParameter(3)+minuitx->GetParError(3);
			effscale2_Draw=minuitx->GetParameter(4);

			for(int pTBin=1;pTBin<pTBinsNew;pTBin++){
				ptCentre_red[pTBin]=pTBin*pTdist;
				effMean_red[pTBin]=graphMC->Eval(pTBin*pTdist/pTscale_Draw-pTshift_Draw)*(1+pTBin*pTdist*effscale_Draw)*effscale2_Draw+effshift_Draw;
			}

			TGraphAsymmErrors *Ngraph_effscalePLUS = new TGraphAsymmErrors(pTBinsNew,ptCentre_red,effMean_red,0,0,0,0);

			Ngraph_effscalePLUS->SetLineColor(kOrange);
			Ngraph_effscalePLUS->SetLineWidth(SystLineWidth);
			Ngraph_effscalePLUS->SetLineStyle(SystLineStyle);
			Ngraph_effscalePLUS->SetName(graphName);
			if(fit_effscale) Ngraph_effscalePLUS->Draw("same");

			pTshift_Draw=minuitx->GetParameter(0);
			pTscale_Draw=minuitx->GetParameter(1);
			effshift_Draw=minuitx->GetParameter(2);
			effscale_Draw=minuitx->GetParameter(3)-minuitx->GetParError(3);
			effscale2_Draw=minuitx->GetParameter(4);

			for(int pTBin=1;pTBin<pTBinsNew;pTBin++){
				ptCentre_red[pTBin]=pTBin*pTdist;
				effMean_red[pTBin]=graphMC->Eval(pTBin*pTdist/pTscale_Draw-pTshift_Draw)*(1+pTBin*pTdist*effscale_Draw)*effscale2_Draw+effshift_Draw;
			}

			TGraphAsymmErrors *Ngraph_effscaleMINUS = new TGraphAsymmErrors(pTBinsNew,ptCentre_red,effMean_red,0,0,0,0);

			Ngraph_effscaleMINUS->SetLineColor(kOrange);
			Ngraph_effscaleMINUS->SetLineWidth(SystLineWidth);
			Ngraph_effscaleMINUS->SetLineStyle(SystLineStyle);
			Ngraph_effscaleMINUS->SetName(graphName);
			if(fit_effscale) Ngraph_effscaleMINUS->Draw("same");

			sprintf(legendentry,"Parametrization, eff_{scale} #pm 1 #sigma");
			if(fit_effscale) SystLegend->AddEntry(Ngraph_effscaleMINUS,legendentry,"l");


			///////// Create effscale2 systematic lines

			pTshift_Draw=minuitx->GetParameter(0);
			pTscale_Draw=minuitx->GetParameter(1);
			effshift_Draw=minuitx->GetParameter(2);
			effscale2_Draw=minuitx->GetParameter(3);
			effscale2_Draw=minuitx->GetParameter(4)+minuitx->GetParError(4);

			for(int pTBin=1;pTBin<pTBinsNew;pTBin++){
				ptCentre_red[pTBin]=pTBin*pTdist;
				effMean_red[pTBin]=graphMC->Eval(pTBin*pTdist/pTscale_Draw-pTshift_Draw)*(1+pTBin*pTdist*effscale_Draw)*effscale2_Draw+effshift_Draw;
			}

			TGraphAsymmErrors *Ngraph_effscale2PLUS = new TGraphAsymmErrors(pTBinsNew,ptCentre_red,effMean_red,0,0,0,0);

			Ngraph_effscale2PLUS->SetLineColor(kOrange+7);
			Ngraph_effscale2PLUS->SetLineWidth(SystLineWidth);
			Ngraph_effscale2PLUS->SetLineStyle(SystLineStyle);
			Ngraph_effscale2PLUS->SetName(graphName);
			if(fit_effscale2) Ngraph_effscale2PLUS->Draw("same");

			pTshift_Draw=minuitx->GetParameter(0);
			pTscale_Draw=minuitx->GetParameter(1);
			effshift_Draw=minuitx->GetParameter(2);
			effscale2_Draw=minuitx->GetParameter(3);
			effscale2_Draw=minuitx->GetParameter(4)-minuitx->GetParError(4);

			for(int pTBin=1;pTBin<pTBinsNew;pTBin++){
				ptCentre_red[pTBin]=pTBin*pTdist;
				effMean_red[pTBin]=graphMC->Eval(pTBin*pTdist/pTscale_Draw-pTshift_Draw)*(1+pTBin*pTdist*effscale_Draw)*effscale2_Draw+effshift_Draw;
			}

			TGraphAsymmErrors *Ngraph_effscale2MINUS = new TGraphAsymmErrors(pTBinsNew,ptCentre_red,effMean_red,0,0,0,0);

			Ngraph_effscale2MINUS->SetLineColor(kOrange+7);
			Ngraph_effscale2MINUS->SetLineWidth(1);
			Ngraph_effscale2MINUS->SetLineStyle(SystLineStyle);
			Ngraph_effscale2MINUS->SetName(graphName);
			if(fit_effscale2) Ngraph_effscale2MINUS->Draw("same");

			sprintf(legendentry,"Parametrization, #tilde{eff}_{scale} #pm 1 #sigma");
			if(fit_effscale2) SystLegend->AddEntry(Ngraph_effscale2MINUS,legendentry,"l");

		}

		///////// Create Correlated lines
		TGraphAsymmErrors* Ngraph_corr[nvpar*2];

		if(!DrawEasySystematicLines){

			int j=0;
			int nCorrSigmas=1;
			for(int i=0;i<nvpar;i++){
				for(int signum=-1;signum<2;signum=signum+2){
					pTshift_Draw=minuitx->GetParameter(0)+signum*varPar[0][i]*nCorrSigmas;
					pTscale_Draw=minuitx->GetParameter(1)+signum*varPar[1][i]*nCorrSigmas;
					effshift_Draw=minuitx->GetParameter(2)+signum*varPar[2][i]*nCorrSigmas;
					effscale_Draw=minuitx->GetParameter(3)+signum*varPar[3][i]*nCorrSigmas;
					effscale2_Draw=minuitx->GetParameter(4)+signum*varPar[4][i]*nCorrSigmas;

					for(int pTBin=0;pTBin<pTBinsNew;pTBin++){
						ptCentre_red[pTBin]=pTBin*pTdist;
						effMean_red[pTBin]=graphMC->Eval(pTBin*pTdist/pTscale_Draw-pTshift_Draw)*(1+pTBin*pTdist*effscale_Draw)*effscale2_Draw+effshift_Draw;
					}

					Ngraph_corr[j] = new TGraphAsymmErrors(pTBinsNew,ptCentre_red,effMean_red,0,0,0,0);

					if(signum<0)Ngraph_corr[j]->SetLineColor(kBlue-6);
					if(signum>0)Ngraph_corr[j]->SetLineColor(kBlue-6);
					Ngraph_corr[j]->SetLineWidth(SystLineWidth);
					Ngraph_corr[j]->SetLineStyle(SystLineStyle);
					Ngraph_corr[j]->SetName(graphName);
					Ngraph_corr[j]->Draw("same");

					j++;
				}
			}
			sprintf(legendentry,"Error band");
			SystLegend->AddEntry(Ngraph_corr[0],legendentry,"l");

		}





		Ngraph->Draw("same");
		graphDATA->Draw("P");

		SystLegend->Draw();

		sprintf(savename,"Scaling/%s/eta%d_Syst.pdf",JobID,etaBin);
		SystCanvas->SaveAs(savename);
		SystCanvas->Close();

		delete SystCanvas;








		TCanvas *ApprovalCanvas = new TCanvas("ApprovalCanvas","ApprovalCanvas",1150,800);
		ApprovalCanvas->SetFillColor(kWhite);
		ApprovalCanvas->GetFrame()->SetFillColor(kWhite);
		ApprovalCanvas->GetFrame()->SetBorderSize(0);
		ApprovalCanvas->SetRightMargin(0.05) ;
		ApprovalCanvas->SetTopMargin(0.05) ;

		double plotEffminAPP=0.45;
		double plotpTmaxAPP=30;

		TH1F *ApprovalHisto = new TH1F;
		ApprovalHisto = ApprovalCanvas->DrawFrame(0,plotEffminAPP,plotpTmaxAPP,1.1);
		ApprovalHisto->SetXTitle("single muon transverse momentum [GeV]");
		ApprovalHisto->SetYTitle("single muon efficiency");
		ApprovalHisto->GetYaxis()->SetTitleSize(0.04);
		//		   		ApprovalHisto->GetYaxis()->SetTitleOffset(1.5);


		TLine* RefLineApproval = new TLine( 0, 1,plotpTmaxAPP, 1 );
		RefLineApproval->SetLineWidth( 2 );
		RefLineApproval->SetLineStyle( 3 );
		RefLineApproval->SetLineColor( kBlack );

		RefLineApproval->Draw( "same" );

		TLegend* ApprovalLegend=new TLegend(0.6375,0.15,0.925,0.35);
		ApprovalLegend->SetFillColor(kWhite);
		ApprovalLegend->SetTextFont(72);
		ApprovalLegend->SetTextSize(0.0375);
		ApprovalLegend->SetBorderSize(1);


		int j=0;
		int nCorrSigmas=1;
		for(int i=0;i<nvpar;i++){

			Ngraph_corr[j]->SetLineColor(kGreen+2);
			Ngraph_corr[j]->Draw("same");

			j++;
		}

		Ngraph->SetLineWidth(2.);
		graphDATA->SetMarkerSize(1.25);

		Ngraph->Draw("same");
		graphDATA->Draw("P");


		double TextAppMin;
		TextAppMin=12.5;
		if(etaBin==0) TextAppMin=TextAppMin+2.;

		char texTexAPP[200];
		sprintf(texTexAPP,"Muon efficiencies, %1.1f < |#eta(#mu)| < %1.1f",etaRange[etaBin],etaRange[etaBin+1]);
		if(etaBin==0) sprintf(texTexAPP,"Muon efficiencies, |#eta(#mu)| < %1.1f", etaRange[etaBin+1]);
		TLatex *textAPP = new TLatex(TextAppMin,1.035,texTexAPP);
		textAPP->SetTextSize(0.045);
		textAPP->Draw( "same" );

		sprintf(texTexAPP,"CMS preliminary");
		TLatex *textAPPPre = new TLatex(20.5,0.665,texTexAPP);
		textAPPPre->SetTextSize(0.035);
		textAPPPre->Draw( "same" );

		sprintf(legendentry,"Data T&P");
		ApprovalLegend->AddEntry(graphDATA,legendentry,"ple");
		sprintf(legendentry,"Parametrization");
		ApprovalLegend->AddEntry(Ngraph,legendentry,"l");
		sprintf(legendentry,"Error curves");
		ApprovalLegend->AddEntry(Ngraph_corr[0],legendentry,"l");
		ApprovalLegend->Draw();

		sprintf(savename,"Scaling/%s/eta%d_Approval.pdf",JobID,etaBin);
		ApprovalCanvas->SaveAs(savename);
		ApprovalCanvas->Close();

		delete ApprovalCanvas;






		myfile << pTscaleEst << "\n";

		bool SaveOutputRootFiles=true;

		if(SaveOutputRootFiles){

			char outfilename[200];
			sprintf(outfilename,"Scaling/%s/ParametrizedFactDataEff_%s_Central.root",JobID,Date);
			TFile *outfile;
			if(etaBin==0) outfile = new TFile(outfilename,"RECREATE");
			if(etaBin>0)  outfile = new TFile(outfilename,"UPDATE");
			outfile->cd();
			Ngraph->Draw("P");
			Ngraph->Write();

			outfile->Write();
			outfile->Close();
			delete outfile;
			outfile = NULL;

			if(!DrawEasySystematicLines){

				for(int i=0;i<nvpar*2;i++){
					if(i==0) sprintf(outfilename,"Scaling/%s/ParametrizedFactDataEff_%s_pTshift_minus.root",JobID,Date);
					if(i==1) sprintf(outfilename,"Scaling/%s/ParametrizedFactDataEff_%s_pTshift_plus.root",JobID,Date);
					if(i==2) sprintf(outfilename,"Scaling/%s/ParametrizedFactDataEff_%s_pTscale_minus.root",JobID,Date);
					if(i==3) sprintf(outfilename,"Scaling/%s/ParametrizedFactDataEff_%s_pTscale_plus.root",JobID,Date);
					if(i==4) sprintf(outfilename,"Scaling/%s/ParametrizedFactDataEff_%s_effshift_minus.root",JobID,Date);
					if(i==5) sprintf(outfilename,"Scaling/%s/ParametrizedFactDataEff_%s_effshift_plus.root",JobID,Date);
					if(etaBin==0) outfile = new TFile(outfilename,"RECREATE");
					if(etaBin>0)  outfile = new TFile(outfilename,"UPDATE");

					outfile->cd();
					Ngraph_corr[i]->Draw("P");
					Ngraph_corr[i]->Write();

					outfile->Write();
					outfile->Close();
					delete outfile;
					outfile = NULL;

				}
			}
		}
	}


	/////////// End of Loop over eta bins ////////////////////////////////////


	myfile.close();

	//		  sprintf(EffFile,"EffFiles/ScaledMCEff_MCtruth18Jan.root");
	//		  hEvalEff2D->SaveAs(EffFile);


	//			return 0;

	for(int etaBin=0;etaBin<etaBins;etaBin++){
		cout<<"etaBin "<<etaBin+1<<" edm = "<<edmEta[etaBin]<<endl;
	}


}
