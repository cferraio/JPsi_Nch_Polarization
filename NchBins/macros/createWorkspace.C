#include "rootIncludes.inc"
#include "commonVar.h"

#include <string>
#include <iostream>
#include <sstream>

using namespace RooFit;

void createWorkspace(const std::string &infilename, int nState, bool correctCtau, bool drawRapPt2D, bool drawPtCPM2D){
	gROOT->SetStyle("Plain");
	gStyle->SetTitleBorderSize(0);

	// Set some strings
	const std::string workspacename = "ws_masslifetime",
				treename = "selectedData";

	// Get the tree from the data file
	TFile *f = TFile::Open(infilename.c_str());
	TTree *tree = (TTree*)f->Get(treename.c_str());

	// Set branch addresses in tree to be able to import tree to roofit
	TLorentzVector* jpsi = new TLorentzVector;
	tree->SetBranchAddress("JpsiP",&jpsi);
	double CPMval = 0;
	tree->SetBranchAddress("CPM",&CPMval);
	double massErr = 0;
	tree->SetBranchAddress("JpsiMassErr",&massErr);
	double Vprob = 0;
	tree->SetBranchAddress("JpsiVprob",&Vprob);
	double lifetime = 0;
	tree->SetBranchAddress("Jpsict",&lifetime);
	double lifetimeErr = 0;
	tree->SetBranchAddress("JpsictErr",&lifetimeErr);

	// define variables necessary for J/Psi(Psi(2S)) mass,lifetime fit
	RooRealVar* JpsiMass =
		new RooRealVar("JpsiMass", "M [GeV]", onia::massMin, onia::massMax);
	RooRealVar* JpsiMassErr =
		new RooRealVar("JpsiMassErr", "#delta M [GeV]", 0, 5);
	RooRealVar* JpsiRap =
		new RooRealVar("JpsiRap", "y", -onia::rap, onia::rap);
	RooRealVar* JpsiPt =
		new RooRealVar("JpsiPt", "p_{T} [GeV]", 0. ,100.);
	RooRealVar* JpsiCPM =
		new RooRealVar("JpsiCPM", "N_{ch}", 0. ,100.);		
	RooRealVar* Jpsict =
		new RooRealVar("Jpsict", "lifetime [mm]", -1., 2.5);
	RooRealVar* JpsictErr =
		new RooRealVar("JpsictErr", "Error on lifetime [mm]", 0.0001, 1);
	RooRealVar* JpsiVprob =
		new RooRealVar("JpsiVprob", "", 0.01, 1.);

	// Set bins
	Jpsict->setBins(10000,"cache");
	Jpsict->setBins(100);
	JpsiMass->setBins(100);
	JpsictErr->setBins(100);

	// The list of data variables    
	RooArgList dataVars(*JpsiMass,*JpsiMassErr,*JpsiRap,*JpsiPt,*JpsiCPM,*Jpsict,*JpsictErr,*JpsiVprob);

	// construct dataset to contain events
	RooDataSet* fullData = new RooDataSet("fullData","The Full Data From the Input ROOT Trees",dataVars);

	int entries = tree->GetEntries();
	cout << "entries " << entries << endl;

	// loop through events in tree and save them to dataset
	for (int ientries = 0; ientries < entries; ientries++) {
	
		if (ientries%100000==0) std::cout << "event " << ientries << " of " << entries <<  std::endl;

		tree->GetEntry(ientries);

		double M =jpsi->M();
		double y=jpsi->Rapidity();
		double pt=jpsi->Pt();
		double cpm=CPMval;


		if (M > JpsiMass->getMin() && M < JpsiMass->getMax()
				&& massErr > JpsiMassErr->getMin() && massErr < JpsiMassErr->getMax()
				&& pt > JpsiPt->getMin() && pt < JpsiPt->getMax()
				&& cpm > JpsiCPM->getMin() && cpm < JpsiCPM->getMax()
				&& y > JpsiRap->getMin() && y < JpsiRap->getMax()
				&& lifetime > Jpsict->getMin() && lifetime < Jpsict->getMax()
				&& lifetimeErr > JpsictErr->getMin() && lifetimeErr < JpsictErr->getMax()
				&& Vprob > JpsiVprob->getMin() && Vprob < JpsiVprob->getMax()
			 ){

			JpsiPt      ->setVal(pt); 
			JpsiCPM		->setVal(cpm);
			JpsiRap     ->setVal(y); 
			JpsiMass    ->setVal(M);
			JpsiMassErr ->setVal(massErr);
			JpsiVprob   ->setVal(Vprob);

			//cout<<"before lifetime correction \n"
			//	<<"Jpsict: "<<lifetime<<" JpsictErr: "<<lifetimeErr<<endl;

			if(correctCtau){
				lifetime    = lifetime    * onia::MpsiPDG / M ;
				lifetimeErr = lifetimeErr * onia::MpsiPDG / M ;
				Jpsict    ->setVal(lifetime);
				JpsictErr ->setVal(lifetimeErr);
				//cout<<"MpsiPDG: "<<onia::MpsiPDG<<endl;
				//cout<<"after lifetime correction \n"
				//	<<"Jpsict: "<<lifetime<<" JpsictErr: "<<lifetimeErr<<endl;
			}
			else{
				Jpsict    ->setVal(lifetime);
				JpsictErr ->setVal(lifetimeErr);
			}

			fullData->add(dataVars);
		}
	}//ientries


	//------------------------------------------------------------------------------------------------------------------
	// Define workspace and import datasets

	////Get datasets binned in pT, cpm, and y

		for(int iRap = 1; iRap <= onia::kNbRapForPTBins; iRap++){ 

		Double_t yMin;
		Double_t yMax;
		if(iRap==0){
			yMin = onia::rapForPTRange[0];
			yMax = onia::rapForPTRange[onia::kNbRapForPTBins];
		} else{
			yMin = onia::rapForPTRange[iRap-1];
			yMax = onia::rapForPTRange[iRap];
		}

		for(int iPT = 1; iPT <= onia::kNbPTBins[iRap]; iPT++){
			//for(int iPT = 0; iPT <= 0; iPT++)

			Double_t ptMin;
			Double_t ptMax;
			if(iPT==0){
				ptMin = onia::pTRange[iRap][0];
				ptMax = onia::pTRange[iRap][onia::kNbPTBins[0]];
			} else{
				ptMin = onia::pTRange[iRap][iPT-1];
				ptMax = onia::pTRange[iRap][iPT];
			}
			
		  for(int iCPM = 1; iCPM <= onia::NchBins; iCPM++){
		  
		  
	  	    Double_t cpmMin;
			Double_t cpmMax;
			if(iCPM==0){
				cpmMin = onia::cpmRange[0];
				cpmMax = onia::cpmRange[onia::NchBins];
			} else{
				cpmMin = onia::cpmRange[iCPM-1];
				cpmMax = onia::cpmRange[iCPM];
			}

			// output file name and workspace
			std::stringstream outfilename;
			outfilename << "tmpFiles/backupWorkSpace/fit_Psi" << nState-3 << "S_rap" << iRap << "_pt" << iPT << "_cpm" << iCPM << ".root";
//			outfilename << "tmpFiles/fit_Psi" << nState-3 << "S_rap" << iRap << "_pt" << iPT << ".root";			
			RooWorkspace* ws = new RooWorkspace(workspacename.c_str());

			// define pt and y cuts on dataset
			std::stringstream cutString;
			cutString << "(JpsiCPM > " << cpmMin << " && JpsiCPM < "<< cpmMax << ") && " << "(JpsiPt >= " << ptMin << " && JpsiPt < "<< ptMax << ") && "
				<< "(TMath::Abs(JpsiRap) >= " << yMin << " && TMath::Abs(JpsiRap) < " << yMax << ")";

			cout << "cutString: " << cutString.str().c_str() << endl;

			// get the dataset for the fit
			RooDataSet* binData = (RooDataSet*)fullData->reduce(cutString.str().c_str());
			std::stringstream name;
			name << "data_rap" << iRap << "_pt" << iPT << "_cpm" << iCPM;;
			binData->SetNameTitle(name.str().c_str(), "Data For Fitting");    

			// Import variables to workspace
			ws->import(*binData);
			ws->writeToFile(outfilename.str().c_str());
			
		  }//iCPM
		}//iPT
	}//iRap

	////---------------------------------------------------------------
	////--Integrating rapidity and pt bins, in +/- 3*sigma mass window
	////---------------------------------------------------------------
	if(drawRapPt2D){
		double yMin = onia::rapForPTRange[0];
		double yMax = onia::rapForPTRange[onia::kNbRapForPTBins];
		double ptMin =  onia::pTRange[0][0];
		double ptMax =  onia::pTRange[0][onia::kNbPTBins[0]];
		double cpmMin =  onia::cpmRange[0];
		double cpmMax =  onia::cpmRange[onia::NchBins];		

		std::stringstream cutRapPt;
		cutRapPt << "(JpsiCPM > " << cpmMin << " && JpsiCPM < "<< cpmMax << ") && "
			<< "(JpsiPt > " << ptMin << " && JpsiPt < "<< ptMax << ") && "
			<< "(TMath::Abs(JpsiRap) > " << yMin << " && TMath::Abs(JpsiRap) < " << yMax << ")";
		cout<<"cutRapPt: "<<cutRapPt.str().c_str()<<endl;

		RooDataSet* rapPtData = (RooDataSet*)fullData->reduce(cutRapPt.str().c_str());
		std::stringstream nameRapPt;
		nameRapPt << "data_rap0_pt0_cpm0";
		rapPtData->SetNameTitle(nameRapPt.str().c_str(), "Data For full rap and pt");

		// output file name and workspace
		std::stringstream outfilename;
		outfilename << "tmpFiles/backupWorkSpace/fit_Psi" << nState-3 << "S_rap0_pt0_cpm0.root";
		RooWorkspace* ws_RapPt = new RooWorkspace(workspacename.c_str());
		//Import variables to workspace
		ws_RapPt->import(*rapPtData);
		ws_RapPt->writeToFile(outfilename.str().c_str());

		TH2D* rapPt;
		TH1D* rap1p2;
		double MassMin;
		double MassMax;

		rap1p2 = new TH1D("rap1p2","rap1p2",30,0, 1.8); 
		if(nState==4){
			rapPt = new TH2D( "rapPt", "rapPt", 52,-1.3,1.3,144,0,72);
			MassMin=3.011;//massPsi1S-onia::nSigMass*sigma1S;
			MassMax=3.174;//massPsi1S+onia::nSigMass*sigma1S;
			// sigma  27.2 MeV
			// mean 3.093 GeV
		}
		if(nState==5){
			rapPt = new TH2D( "rapPt", "rapPt", 64,-1.6,1.6,144,0,72); //  rap<1.5
			//rapPt = new TH2D( "rapPt", "rapPt", 52,-1.3,1.3,144,0,72); //  rap<1.2
			MassMin=3.576;//massPsi2S-onia::nSigMass*sigma2S;
			MassMax=3.786;//massPsi2S+onia::nSigMass*sigma2S;
			// sigma 34.9 MeV // pT > 7
			// sigma 34.3 MeV // pT > 10
			// mean 3.681 GeV
		}

		cout<<"Plotting rap-Pt for Psi"<<nState-3<<"S"<<endl;
		cout<<"MassMin for rap-Pt plot = "<<MassMin<<endl;
		cout<<"MassMax for rap-Pt plot = "<<MassMax<<endl;

		TTree *rapPtTree = (TTree*)rapPtData->tree();
		std::stringstream cutMass;
		cutMass<<"(JpsiMass > " << MassMin << " && JpsiMass < "<< MassMax << ")";
		//following two methods can only be used in root_v30, 34 does not work
		rapPtTree->Draw("JpsiPt:JpsiRap>>rapPt",cutMass.str().c_str(),"colz");
		cout<<"debug"<<endl;
		rapPtTree->Draw("TMath::Abs(JpsiRap)>>rap1p2",cutMass.str().c_str());

		TCanvas* c2 = new TCanvas("c2","c2",1200,1500);
		rapPt->SetYTitle("p_{T}(#mu#mu) [GeV]");
		rapPt->SetXTitle("y(#mu#mu)");
		gStyle->SetPalette(1);
		gPad->SetFillColor(kWhite);
		rapPt->SetTitle(0);
		rapPt->SetStats(0);
		gPad->SetLeftMargin(0.15);
		gPad->SetRightMargin(0.17);
		rapPt->GetYaxis()->SetTitleOffset(1.5);
		rapPt->Draw("colz");

		TLine* rapPtLine;

		for(int iRap=0;iRap<onia::kNbRapForPTBins+1;iRap++){
			rapPtLine= new TLine( -onia::rapForPTRange[iRap], onia::pTRange[0][0], -onia::rapForPTRange[iRap], onia::pTRange[0][onia::kNbPTBins[iRap]] );
			rapPtLine->SetLineWidth( 2 );
			rapPtLine->SetLineStyle( 1 );
			rapPtLine->SetLineColor( kWhite );
			rapPtLine->Draw();
			rapPtLine= new TLine( onia::rapForPTRange[iRap], onia::pTRange[0][0], onia::rapForPTRange[iRap], onia::pTRange[0][onia::kNbPTBins[iRap]] );
			rapPtLine->SetLineWidth( 2 );
			rapPtLine->SetLineStyle( 1 );
			rapPtLine->SetLineColor( kWhite );
			rapPtLine->Draw();
			int pTBegin = 0;
			if(nState==5) pTBegin = 1;
			for(int iPt=pTBegin;iPt<onia::kNbPTBins[iRap]+1;iPt++){
				rapPtLine= new TLine( -onia::rapForPTRange[onia::kNbRapForPTBins], onia::pTRange[0][iPt], onia::rapForPTRange[onia::kNbRapForPTBins], onia::pTRange[0][iPt] );
				rapPtLine->SetLineWidth( 2 );
				rapPtLine->SetLineStyle( 1 );
				rapPtLine->SetLineColor( kWhite );
				rapPtLine->Draw();
			}
		}

		char savename[200];
		sprintf(savename,"Fit/rapPt_Psi%dS.pdf",nState-3);
		c2->SaveAs(savename);
		TCanvas* c3 = new TCanvas("c3","c3",1500,1200);
		rap1p2->SetYTitle("Events");
		rap1p2->SetXTitle("y(#mu#mu)");
		rap1p2->SetTitle(0);
		rap1p2->SetStats(0);
		rap1p2->GetYaxis()->SetTitleOffset(1.2);
		rap1p2->Draw();
		sprintf(savename,"Fit/rapDimuon_1p2_Psi%dS.pdf",nState-3);
		c3->SaveAs(savename);
	}
	
	if(drawPtCPM2D){
		double yMin = onia::rapForPTRange[0];
		double yMax = onia::rapForPTRange[onia::kNbRapForPTBins];
		double ptMin =  onia::pTRange[0][0];
		double ptMax =  onia::pTRange[0][onia::kNbPTBins[0]];
		double cpmMin =  onia::cpmRange[0];
		double cpmMax =  onia::cpmRange[onia::NchBins];		

		std::stringstream cutRapPt;
		cutRapPt << "(JpsiCPM > " << cpmMin << " && JpsiCPM < "<< cpmMax << ") && "
			<< "(JpsiPt > " << ptMin << " && JpsiPt < "<< ptMax << ") && "
			<< "(TMath::Abs(JpsiRap) > " << yMin << " && TMath::Abs(JpsiRap) < " << yMax << ")";
		cout<<"cutRapPt: "<<cutRapPt.str().c_str()<<endl;

		RooDataSet* rapPtData = (RooDataSet*)fullData->reduce(cutRapPt.str().c_str());
		std::stringstream nameRapPt;
		nameRapPt << "data_rap0_pt0_cpm0";
		rapPtData->SetNameTitle(nameRapPt.str().c_str(), "Data For full rap and pt");

		// output file name and workspace
		std::stringstream outfilename;
		outfilename << "tmpFiles/backupWorkSpace/fit_Psi" << nState-3 << "S_rap0_pt0_cpm0.root";
		RooWorkspace* ws_RapPt = new RooWorkspace(workspacename.c_str());
		//Import variables to workspace
		ws_RapPt->import(*rapPtData);
		ws_RapPt->writeToFile(outfilename.str().c_str());

		TH2D* PtCPM;
		double MassMin;
		double MassMax;

		if(nState==4){
			PtCPM = new TH2D( "PtCPM", "PtCPM", 100,0,50,200,0,100);
			MassMin=3.011;//massPsi1S-onia::nSigMass*sigma1S;
			MassMax=3.174;//massPsi1S+onia::nSigMass*sigma1S;
			// sigma  27.2 MeV
			// mean 3.093 GeV
		}
		if(nState==5){
			PtCPM = new TH2D( "PtCPM", "PtCPM", 100,0,50,200,0,100); //  rap<1.5
			//rapPt = new TH2D( "rapPt", "rapPt", 52,-1.3,1.3,144,0,72); //  rap<1.2
			MassMin=3.576;//massPsi2S-onia::nSigMass*sigma2S;
			MassMax=3.786;//massPsi2S+onia::nSigMass*sigma2S;
			// sigma 34.9 MeV // pT > 7
			// sigma 34.3 MeV // pT > 10
			// mean 3.681 GeV
		}

		cout<<"Plotting Pt-CPM for Psi"<<nState-3<<"S"<<endl;
		cout<<"MassMin for Pt-CPM plot = "<<MassMin<<endl;
		cout<<"MassMax for Pt-CPM plot = "<<MassMax<<endl;

		TTree *rapPtTree = (TTree*)rapPtData->tree();
		std::stringstream cutMass;
		cutMass<<"(JpsiMass > " << MassMin << " && JpsiMass < "<< MassMax << ")";
		//following two methods can only be used in root_v30, 34 does not work
		rapPtTree->Draw("JpsiCPM:JpsiPt>>PtCPM",cutMass.str().c_str(),"colz");
		cout<<"debug"<<endl;

		TCanvas* c2 = new TCanvas("c2","c2",1200,1500);
		PtCPM->SetYTitle("N_{ch}");
		PtCPM->SetXTitle("p_{T}(#mu#mu) [GeV]");
		gStyle->SetPalette(1);
		gPad->SetFillColor(kWhite);
		PtCPM->SetTitle(0);
		PtCPM->SetStats(0);
		gPad->SetLeftMargin(0.15);
		gPad->SetRightMargin(0.17);
		PtCPM->GetYaxis()->SetTitleOffset(1.5);
		PtCPM->Draw("colz");

		TLine* PtCPMLine;

		int iRap=0;
		for(int iPt=0;iPt<onia::kNbPTMaxBins+1;iPt++){
		int cpmBegin = 0;
			if(nState==5) cpmBegin = 1;
			for(int icpm=cpmBegin;icpm<onia::NchBins+1;icpm++){
			
			PtCPMLine= new TLine( onia::pTRange[iRap][0], onia::cpmRange[icpm], onia::pTRange[iRap][onia::kNbPTMaxBins], onia::cpmRange[icpm] );
			PtCPMLine->SetLineWidth( 2 );
			PtCPMLine->SetLineStyle( 1 );
			PtCPMLine->SetLineColor( kWhite );
			PtCPMLine->Draw();
			PtCPMLine= new TLine( onia::pTRange[iRap][iPt], onia::cpmRange[0], onia::pTRange[iRap][iPt], onia::cpmRange[onia::NchBins] );
			PtCPMLine->SetLineWidth( 2 );
			PtCPMLine->SetLineStyle( 1 );
			PtCPMLine->SetLineColor( kWhite );
			PtCPMLine->Draw();

			
//				PtCPMLine= new TLine( onia::pTRange[0][onia::kNbPTMaxBins], onia::cpmRange[icpm], onia::pTRange[0][onia::kNbPTMaxBins], onia::cpmRange[icpm] );
//				PtCPMLine->SetLineWidth( 2 );
//				PtCPMLine->SetLineStyle( 1 );
//				PtCPMLine->SetLineColor( kWhite );
//				PtCPMLine->Draw();
			}
		}

		char savename[200];
		sprintf(savename,"Fit/PtCPM_Psi%dS.pdf",nState-3);
		c2->SaveAs(savename);
	}

	f->Close();
}
