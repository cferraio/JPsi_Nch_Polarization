/*
 * CalcPlotAngles.C
 *
 *  Created on: Feb 25, 2013
 *      Author: valentinknuenz
 */

#include "Riostream.h"
#include "TSystem.h"
#include "TString.h"
#include "TFile.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TTree.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "TRandom3.h"
#include "TF1.h"
#include "TF2.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TGraphAsymmErrors.h"

void CalcPlotAngles(){

	gROOT->Reset();
	gROOT->SetBatch();
	gStyle->SetPalette(1);


	char JobID[200];
	char BaseID[200];
	char DataID[200];
	char FileID[200];


	int StateMin=4;
	int StateMax=4;

	int pTBinMin=9;
	int pTBinMax=9;

	int RapBinMin=1;
	int RapBinMax=1;

	bool debug=false;
	bool plotToyMC=false;

	sprintf(BaseID,"/afs/ihep.ac.cn/users/z/zhangll/fs/work/polarization/PsiPol2011/macros/DataFiles");
	sprintf(DataID,"SetOfCuts11_ctauScen5_FracLSB-1_newMLfit_30Apr2013_correctfLSB_test");
	//sprintf(DataID,"SetOfCuts11_ctauScen5_MC");
	sprintf(JobID,"%s_June13",DataID);

	if(plotToyMC){
		sprintf(BaseID,"/afs/ihep.ac.cn/users/z/zhangll/fs/work/polarization/PsiPol2011/Psi/ToyMC");
		//sprintf(DataID,"ToyMC_Psi1S_22Feb2013_FiducialTest/Sig_frame3scen3_Bkg_frame3scen3");
	}

	//sprintf(FileID,"%s/%s", BaseID, DataID);


	const double pbeam_ = 7000.; // exact number irrelevant as long as pbeam >> Mprot
	const double Mprot_ = 0.9382720;
	const double Mlepton_ = 0.10566;  // (muon)
	const double gPI_ = TMath::Pi();
	const double Ebeam_ = sqrt( pbeam_*pbeam_ + Mprot_*Mprot_ );
	TLorentzVector beam1_LAB_( 0., 0., pbeam_, Ebeam_ );
	TLorentzVector beam2_LAB_( 0., 0., -pbeam_, Ebeam_ );

	char FileDir[200];
	sprintf(FileDir,"AnglePlots");
	gSystem->mkdir(FileDir);
	sprintf(FileDir,"AnglePlots/%s",JobID);
	gSystem->mkdir(FileDir);


	for( int nState  = StateMin; nState  < StateMax+1; nState++ ){
		for( int irap  = RapBinMin; irap  < RapBinMax+1; irap++ ){ 
			if(nState==4&&irap==3) continue;
			for( int ipT   = pTBinMin; ipT   < pTBinMax+1;  ipT++  ){ 
				if(nState==5&&(ipT>6||ipT<2)) continue;

				//sprintf(DataID,"ToyMC_Psi%dS_22Feb2013_FiducialTest/Sig_frame3scen3_Bkg_frame3scen3",nState-3);
	      sprintf(DataID,"SetOfCuts11_ctauScen5_FracLSB-1_newMLfit_30Apr2013_correctfLSB_test");
	      //sprintf(DataID,"SetOfCuts11_ctauScen5_MC");
				sprintf(FileID,"%s/%s", BaseID, DataID);

				////////// HISTOS ///////////


				double yOffset=1.6;
				int MarkStyle=20;

				int nBinsCosth2D=16;
				int nBinsPhi2D=16;
				int nBinsCosth1D=32;
				int nBinsPhi1D=32;

				TH2D* costhphiPX   = new TH2D( "costhphiPX", "costhphiPX", nBinsCosth2D, -1,1, nBinsPhi2D,  -180, 180);
				TH1D* costhPX   = new TH1D( "costhPX", "costhPX", nBinsCosth1D, -1,1);
				TH1D* phiPX   = new TH1D( "phiPX", "phiPX", nBinsPhi1D,  -180, 180);

				costhphiPX->GetYaxis()->SetTitleOffset(yOffset);
				costhphiPX->GetXaxis()->SetTitle("cos#vartheta");
				costhphiPX->GetYaxis()->SetTitle("#varphi [deg.]");
				costhphiPX->SetStats(0);
				costhphiPX->SetTitle(0);

				costhPX->GetXaxis()->SetTitle("cos#vartheta");
				costhPX->SetStats(0);
				costhPX->SetTitle(0);
				costhPX->SetMarkerStyle(MarkStyle);

				phiPX->GetXaxis()->SetTitle("#varphi [deg.]");
				phiPX->SetStats(0);
				phiPX->SetTitle(0);
				phiPX->SetMarkerStyle(MarkStyle);

				TH2D* costhphiHX   = new TH2D( "costhphiHX", "costhphiHX", nBinsCosth2D, -1,1, nBinsPhi2D,  -180, 180);
				TH1D* costhHX   = new TH1D( "costhHX", "costhHX", nBinsCosth1D, -1,1);
				TH1D* costhHX_1   = new TH1D( "costhHX_1", "costhHX_1", nBinsCosth1D, -1,1);
				TH1D* phiHX   = new TH1D( "phiHX", "phiHX", nBinsPhi1D,  -180, 180);
				//TH1D* phiHX   = new TH1D( "phiHX", "phiHX", nBinsPhi1D,  70, 100);

				costhphiHX->GetYaxis()->SetTitleOffset(yOffset);
				costhphiHX->GetXaxis()->SetTitle("cos#vartheta");
				costhphiHX->GetYaxis()->SetTitle("#varphi [deg.]");
				costhphiHX->SetStats(0);
				costhphiHX->SetTitle(0);

				costhHX->GetXaxis()->SetTitle("cos#vartheta");
				costhHX->SetStats(0);
				costhHX->SetTitle(0);
				costhHX->SetMarkerStyle(MarkStyle);

				costhHX_1->GetXaxis()->SetTitle("cos#vartheta");
				costhHX_1->SetStats(0);
				costhHX_1->SetTitle(0);
				costhHX_1->SetMarkerStyle(MarkStyle);
				costhHX_1->SetMarkerSize(0.8);
				costhHX_1->SetMarkerColor(kRed);
				costhHX_1->SetLineColor(kRed);

				phiHX->GetXaxis()->SetTitle("#varphi [deg.]");
				phiHX->SetStats(0);
				phiHX->SetTitle(0);
				phiHX->SetMarkerStyle(MarkStyle);

				TH2D* costhphiCS   = new TH2D( "costhphiCS", "costhphiCS", nBinsCosth2D, -1,1, nBinsPhi2D,  -180, 180);
				TH1D* costhCS   = new TH1D( "costhCS", "costhCS", nBinsCosth1D, -1,1);
				TH1D* phiCS   = new TH1D( "phiCS", "phiCS", nBinsPhi1D,  -180, 180);

				costhphiCS->GetYaxis()->SetTitleOffset(yOffset);
				costhphiCS->GetXaxis()->SetTitle("cos#vartheta");
				costhphiCS->GetYaxis()->SetTitle("#varphi [deg.]");
				costhphiCS->SetStats(0);
				costhphiCS->SetTitle(0);

				costhCS->GetXaxis()->SetTitle("cos#vartheta");
				costhCS->SetStats(0);
				costhCS->SetTitle(0);
				costhCS->SetMarkerStyle(MarkStyle);

				phiCS->GetXaxis()->SetTitle("#varphi [deg.]");
				phiCS->SetStats(0);
				phiCS->SetTitle(0);
				phiCS->SetMarkerStyle(MarkStyle);

				//TH2D* costhCSphiHX   = new TH2D( "costhCSphiHX", "costhCSphiHX", nBinsCosth2D, -1,1, nBinsPhi2D,  -180, 180);
				TH2D* costhCSphiHX   = new TH2D( "costhCSphiHX", "costhCSphiHX", nBinsCosth2D/0.5, -0.4,0.4, nBinsPhi2D*2,  70., 110.);
				costhCSphiHX->GetYaxis()->SetTitleOffset(yOffset);
				costhCSphiHX->GetXaxis()->SetTitle("cos#vartheta");
				costhCSphiHX->GetYaxis()->SetTitle("#varphi [deg.]");
				costhCSphiHX->SetStats(0);
				costhCSphiHX->SetTitle(0);

				int nBinsMass=50; double MassMin=2.9; double MassMax=3.3; if(nState==5) {MassMin=3.5; MassMax=3.9;}
				int nBinspT=300; double pTMin=10; double pTMax=70;
				int nBinsrap=50; double rapMin=-2; double rapMax=2;

				int nBinsmuonpT=300; double muonpTMin=0; double muonpTMax=50;


				TH1D* dimuonMass   = new TH1D( "dimuonMass", "dimuonMass", nBinsMass, MassMin, MassMax);
				TH1D* dimuonpT   = new TH1D( "dimuonpT", "dimuonpT", nBinspT, pTMin, pTMax);
				TH1D* dimuonrap   = new TH1D( "dimuonrap", "dimuonrap", nBinsrap, rapMin, rapMax);
				//TH1D* dimuonrap   = new TH1D( "dimuonrap", "dimuonrap", nBinsrap, -0.65,0.65);
				TH1D* muonpT   = new TH1D( "muonpT", "muonpT", nBinsmuonpT, muonpTMin, muonpTMax);
				TH1D* muoneta   = new TH1D( "muoneta", "muoneta", nBinsrap, rapMin, rapMax);


				dimuonMass->GetXaxis()->SetTitle("dimuon mass [GeV]");
				dimuonMass->SetStats(0);
				dimuonMass->SetTitle(0);
				dimuonMass->SetMarkerStyle(MarkStyle);

				dimuonpT->GetXaxis()->SetTitle("dimuon pT [GeV]");
				dimuonpT->SetStats(0);
				dimuonpT->SetTitle(0);
				dimuonpT->SetMarkerStyle(MarkStyle);

				dimuonrap->GetXaxis()->SetTitle("dimuon rapidity [GeV]");
				dimuonrap->SetStats(0);
				dimuonrap->SetTitle(0);
				//dimuonrap->SetTitle("phi_HX>75. && phi_HX<100. && costh_CS>-0.2 && costh_CS<0.2");
				dimuonrap->SetMarkerStyle(MarkStyle);

				muonpT->GetXaxis()->SetTitle("muon pT [GeV]");
				muonpT->SetStats(0);
				muonpT->SetTitle(0);
				muonpT->SetMarkerStyle(MarkStyle);

				muoneta->GetXaxis()->SetTitle("muon rapidity [GeV]");
				muoneta->SetStats(0);
				muoneta->SetTitle(0);
				muoneta->SetMarkerStyle(MarkStyle);


				////////////////////////////////

				char InFileName[200];
				sprintf(InFileName,"%s/Psi%dS/tmpFiles/data_Psi%dS_rap%d_pT%d.root",FileID, nState-3, nState-3, irap, ipT);
				//if(plotToyMC) sprintf(InFileName,"%s/rap%d_pT%d/Generation1/data.root",FileID, irap, ipT);
				if(plotToyMC) sprintf(InFileName,"%s/rap%d_pT%d/data.root",FileID, irap, ipT);


				TFile* inFile = new TFile(InFileName,"R"); //"UPDATE");

				TTree* data = (TTree*)inFile->Get("selectedData");

				char OutFileName[200];
				sprintf(OutFileName,"%s/Psi%dS/tmpFiles/data_Psi%dS_rap%d_pT%d_reduced.root",FileID, nState-3, nState-3, irap, ipT);
				TFile* outFile = new TFile(OutFileName,"RECREATE");
				TTree *reducedData = new TTree ("selectedData", "selected events");

				TLorentzVector* lepP  = new TLorentzVector();  data->SetBranchAddress( "lepP",   &lepP );  // lepton 4-vectors
				TLorentzVector* lepN  = new TLorentzVector();  data->SetBranchAddress( "lepN",   &lepN );
				TLorentzVector* JpsiP = new TLorentzVector();  data->SetBranchAddress( "JpsiP",  &JpsiP );
				double Jpsict = 0;                             data->SetBranchAddress( "Jpsict", &Jpsict );
				double JpsictErr = 0;                          data->SetBranchAddress( "JpsictErr", &JpsictErr );
				double JpsiMassErr = 0;                        data->SetBranchAddress( "JpsiMassErr", &JpsiMassErr );
				double JpsiVprob = 0;                          data->SetBranchAddress( "JpsiVprob", &JpsiVprob );

				reducedData->Branch("lepP", "TLorentzVector", &lepP);
				reducedData->Branch("lepN", "TLorentzVector", &lepN);
				reducedData->Branch("JpsiP", "TLorentzVector", &JpsiP);
				reducedData->Branch("Jpsict", &Jpsict, "Jpsict/D");
				reducedData->Branch("JpsictErr", &JpsictErr, "JpsictErr/D");
				reducedData->Branch("JpsiMassErr", &JpsiMassErr, "JpsiMassErr/D");
				reducedData->Branch("JpsiVprob", &JpsiVprob, "JpsiVprob/D");

				int n_events=data->GetEntries();

				int n_step = n_events/5;  // to visualize progress of the event scan (50 steps)
				int n_step_=1;


				for ( int i_event = 1; i_event <= n_events; i_event++ ) {
					if (i_event%n_step == 0) {cout << n_step_*20 <<" % "<<endl; n_step_++;}

					if(debug&&i_event>n_events/100) continue;

					data->GetEvent( i_event-1 );

					//if( JpsiP->Pt() < 30. || JpsiP->Pt() > 35. ) cout << "JpsiP->Pt() " << JpsiP->Pt() << endl;

					//if( !(JpsiP->Pt() >= 34. && JpsiP->Pt() <= 35.) ) continue;

					double lepP_pT  = lepP->Pt();
					double lepN_pT  = lepN->Pt();

					double lepP_eta = lepP->PseudoRapidity();
					double lepN_eta = lepN->PseudoRapidity();

					// dilepton 4-vector:

					TLorentzVector dilepton = *lepP + *lepN;
					double pT   = dilepton.Pt();
					double rap  = dilepton.Rapidity();
					double mass = dilepton.M();

					// calculation of decay angles in three polarization frames

					// reference directions to calculate angles:

					TVector3 lab_to_dilep = -dilepton.BoostVector();

					TLorentzVector beam1_DILEP = beam1_LAB_;
					beam1_DILEP.Boost(lab_to_dilep);         // beam1 in the dilepton rest frame
					TLorentzVector beam2_DILEP = beam2_LAB_;
					beam2_DILEP.Boost(lab_to_dilep);         // beam2 in the dilepton rest frame

					TVector3 beam1_direction     = beam1_DILEP.Vect().Unit();
					TVector3 beam2_direction     = beam2_DILEP.Vect().Unit();
					TVector3 dilep_direction     = dilepton.Vect().Unit();
					TVector3 beam1_beam2_bisect  = ( beam1_direction - beam2_direction ).Unit();


					// all polarization frames have the same Y axis = the normal to the plane formed by
					// the directions of the colliding hadrons:

					TVector3 Yaxis = ( beam1_direction.Cross( beam2_direction ) ).Unit();

					// flip of y axis with rapidity:

					if ( rap < 0. ) Yaxis = - Yaxis;

					TVector3 perpendicular_to_beam = ( beam1_beam2_bisect.Cross( Yaxis ) ).Unit();


					// positive lepton in the dilepton rest frame:

					TLorentzVector lepton_DILEP = *lepP;
					lepton_DILEP.Boost(lab_to_dilep);

					// CS frame angles:

					TVector3 newZaxis = beam1_beam2_bisect;
					TVector3 newYaxis = Yaxis;
					TVector3 newXaxis = newYaxis.Cross( newZaxis );

					TRotation rotation;
					rotation.SetToIdentity();
					rotation.RotateAxes( newXaxis, newYaxis, newZaxis );
					rotation.Invert();  // transforms coordinates from the "xyz" frame to the new frame
					TVector3 lepton_DILEP_rotated = lepton_DILEP.Vect();
					lepton_DILEP_rotated.Transform(rotation);

					double costh_CS = lepton_DILEP_rotated.CosTheta();
					double phi_CS   = lepton_DILEP_rotated.Phi() * 180. / gPI_;
					double phith_CS;
					if ( costh_CS < 0. ) phith_CS = phi_CS - 135.;
					if ( costh_CS > 0. ) phith_CS = phi_CS - 45.;
					if ( phith_CS < -180. ) phith_CS = 360. + phith_CS;


					// HELICITY frame angles:

					newZaxis = dilep_direction;
					newYaxis = Yaxis;
					newXaxis = newYaxis.Cross( newZaxis );

					rotation.SetToIdentity();
					rotation.RotateAxes( newXaxis, newYaxis, newZaxis );
					rotation.Invert();
					lepton_DILEP_rotated = lepton_DILEP.Vect();
					lepton_DILEP_rotated.Transform(rotation);

					double costh_HX = lepton_DILEP_rotated.CosTheta();
					double phi_HX   = lepton_DILEP_rotated.Phi() * 180. / gPI_;
					double phith_HX;
					if ( costh_HX < 0. ) phith_HX = phi_HX - 135.;
					if ( costh_HX > 0. ) phith_HX = phi_HX - 45.;
					if ( phith_HX < -180. ) phith_HX = 360. + phith_HX;


					// PERPENDICULAR HELICITY frame angles:

					newZaxis = perpendicular_to_beam;
					newYaxis = Yaxis;
					newXaxis = newYaxis.Cross( newZaxis );

					rotation.SetToIdentity();
					rotation.RotateAxes( newXaxis, newYaxis, newZaxis );
					rotation.Invert();
					lepton_DILEP_rotated = lepton_DILEP.Vect();
					lepton_DILEP_rotated.Transform(rotation);

					double costh_PX = lepton_DILEP_rotated.CosTheta();
					double phi_PX   = lepton_DILEP_rotated.Phi() * 180. / gPI_;
					double phith_PX;
					if ( costh_PX < 0. ) phith_PX = phi_PX - 135.;
					if ( costh_PX > 0. ) phith_PX = phi_PX - 45.;
					if ( phith_PX < -180. ) phith_PX = 360. + phith_PX;



					// fill histograms
					costhphiPX->Fill( costh_PX, phi_PX, 1. );
					costhPX->Fill( costh_PX, 1. );
					phiPX->Fill( phi_PX, 1. );
					costhphiHX->Fill( costh_HX, phi_HX, 1. );
					costhHX->Fill( costh_HX, 1. );
					if(!(phi_HX>80. && phi_HX<95.)) costhHX_1->Fill( costh_HX, 1. ); 
					//if(phi_HX>70. && phi_HX<100.) 
						phiHX->Fill( phi_HX, 1. );
					costhphiCS->Fill( costh_CS, phi_CS, 1. );
					//if(phi_HX<80. || phi_HX>91.) 
						costhCS->Fill( costh_CS, 1. );
					phiCS->Fill( phi_CS, 1. );

					if(costh_CS>-0.4 && costh_CS<0.4 && phi_HX>70. && phi_HX<110.) 
						costhCSphiHX->Fill( costh_CS, phi_HX, 1. );

					dimuonMass->Fill( mass, 1. );
					dimuonpT->Fill( pT, 1. );
					//if(phi_HX>75. && phi_HX<100. && costh_CS>-0.2 && costh_CS<0.2) 
						dimuonrap->Fill( rap, 1. );
					muonpT->Fill( lepP_pT, 1. );
					muoneta->Fill( lepP_eta, 1. );
					muonpT->Fill( lepN_pT, 1. );
					muoneta->Fill( lepN_eta, 1. );

					//if(phi_HX<80. || phi_HX>91.)  reducedData->Fill();
					if(phi_HX>85. && phi_HX<95.)  reducedData->Fill();

				}

				//inFile->cd();
				//reducedData->Write();
				//inFile->Close();

				outFile->cd();
				reducedData->Write();
				outFile->Close();


				//Plot


				char plotname[200];

				TCanvas *c3;


				c3 = new TCanvas(plotname, plotname, 500, 500);
				c3->SetFillColor(kWhite);
				c3->SetTopMargin(0.05);
				c3->SetLeftMargin(0.13);
				c3->SetRightMargin(0.13); //
				costhphiPX->Draw("colz");
				sprintf(plotname, "%s/Psi%dS_costhphiPX_rap%d_pT%d.pdf", FileDir, nState-3, irap, ipT);
				c3->SaveAs(plotname);

				c3 = new TCanvas(plotname, plotname, 700, 500);
				c3->SetFillColor(kWhite);
				c3->SetTopMargin(0.05);
				c3->SetLeftMargin(0.13);
				c3->SetRightMargin(0.05);
				costhPX->Draw("e");
				sprintf(plotname, "%s/Psi%dS_costhPX_rap%d_pT%d.pdf", FileDir, nState-3, irap, ipT);
				c3->SaveAs(plotname);

				c3 = new TCanvas(plotname, plotname, 700, 500);
				c3->SetFillColor(kWhite);
				c3->SetTopMargin(0.05);
				c3->SetLeftMargin(0.13);
				c3->SetRightMargin(0.05);
				phiPX->GetYaxis()->SetRangeUser(0., phiPX->GetMaximum()*2.);
				phiPX->Draw("e");
				sprintf(plotname, "%s/Psi%dS_phiPX_rap%d_pT%d.pdf", FileDir, nState-3, irap, ipT);
				c3->SaveAs(plotname);


				c3 = new TCanvas(plotname, plotname, 500, 500);
				c3->SetFillColor(kWhite);
				c3->SetTopMargin(0.05);
				c3->SetLeftMargin(0.13);
				c3->SetRightMargin(0.13); //
				costhphiHX->Draw("colz");
				sprintf(plotname, "%s/Psi%dS_costhphiHX_rap%d_pT%d.pdf", FileDir, nState-3, irap, ipT);
				c3->SaveAs(plotname);

				c3 = new TCanvas(plotname, plotname, 700, 500);
				c3->SetFillColor(kWhite);
				c3->SetTopMargin(0.05);
				c3->SetLeftMargin(0.13);
				c3->SetRightMargin(0.05);
				costhHX->Draw("e");
				costhHX_1->Draw("esame");
				sprintf(plotname, "%s/Psi%dS_costhHX_rap%d_pT%d.pdf", FileDir, nState-3, irap, ipT);
				c3->SaveAs(plotname);
				costhHX->Scale(1./costhHX->Integral());
				costhHX_1->Scale(1./costhHX_1->Integral());
				costhHX->Draw("e");
				costhHX_1->Draw("esame");
				sprintf(plotname, "%s/Psi%dS_costhHX_rap%d_pT%d_Norm.pdf", FileDir, nState-3, irap, ipT);
				c3->SaveAs(plotname);

				c3 = new TCanvas(plotname, plotname, 700, 500);
				c3->SetFillColor(kWhite);
				c3->SetTopMargin(0.05);
				c3->SetLeftMargin(0.13);
				c3->SetRightMargin(0.05);
				phiHX->GetYaxis()->SetRangeUser(0., phiHX->GetMaximum()*2.);
				phiHX->Draw("e");
				sprintf(plotname, "%s/Psi%dS_phiHX_rap%d_pT%d.pdf", FileDir, nState-3, irap, ipT);
				c3->SaveAs(plotname);


				c3 = new TCanvas(plotname, plotname, 500, 500);
				c3->SetFillColor(kWhite);
				c3->SetTopMargin(0.05);
				c3->SetLeftMargin(0.13);
				c3->SetRightMargin(0.13); //
				costhphiCS->Draw("colz");
				sprintf(plotname, "%s/Psi%dS_costhphiCS_rap%d_pT%d.pdf", FileDir, nState-3, irap, ipT);
				c3->SaveAs(plotname);

				c3 = new TCanvas(plotname, plotname, 700, 500);
				c3->SetFillColor(kWhite);
				c3->SetTopMargin(0.05);
				c3->SetLeftMargin(0.13);
				c3->SetRightMargin(0.05);
				costhCS->Draw("e");
				sprintf(plotname, "%s/Psi%dS_costhCS_rap%d_pT%d.pdf", FileDir, nState-3, irap, ipT);
				c3->SaveAs(plotname);

				c3 = new TCanvas(plotname, plotname, 700, 500);
				c3->SetFillColor(kWhite);
				c3->SetTopMargin(0.05);
				c3->SetLeftMargin(0.13);
				c3->SetRightMargin(0.05);
				//phiCS->GetYaxis()->SetRangeUser(0., phiCS->GetMaximum()*2.);
				phiCS->Draw("e");
				sprintf(plotname, "%s/Psi%dS_phiCS_rap%d_pT%d.pdf", FileDir, nState-3, irap, ipT);
				c3->SaveAs(plotname);

				//
				c3 = new TCanvas(plotname, plotname, 500, 500);
				c3->SetFillColor(kWhite);
				c3->SetTopMargin(0.05);
				c3->SetLeftMargin(0.13);
				c3->SetRightMargin(0.13); //
				costhCSphiHX->Draw("colz");
				sprintf(plotname, "%s/Psi%dS_costhCSphiHX_rap%d_pT%d.pdf", FileDir, nState-3, irap, ipT);
				c3->SaveAs(plotname);

				c3 = new TCanvas(plotname, plotname, 700, 500);
				c3->SetFillColor(kWhite);
				c3->SetTopMargin(0.05);
				c3->SetLeftMargin(0.13);
				c3->SetRightMargin(0.05);
				dimuonMass->Draw("e");
				sprintf(plotname, "%s/Psi%dS_dimuonMass_rap%d_pT%d.pdf", FileDir, nState-3, irap, ipT);
				c3->SaveAs(plotname);

				c3 = new TCanvas(plotname, plotname, 1700, 500);
				c3->SetFillColor(kWhite);
				c3->SetTopMargin(0.05);
				c3->SetLeftMargin(0.13);
				c3->SetRightMargin(0.05);
				dimuonpT->Draw("e");
				sprintf(plotname, "%s/Psi%dS_dimuonpT_rap%d_pT%d.pdf", FileDir, nState-3, irap, ipT);
				c3->SaveAs(plotname);

				c3 = new TCanvas(plotname, plotname, 700, 500);
				c3->SetFillColor(kWhite);
				c3->SetTopMargin(0.05);
				c3->SetLeftMargin(0.13);
				c3->SetRightMargin(0.05);
				dimuonrap->Draw("e");
				sprintf(plotname, "%s/Psi%dS_dimuonrap_rap%d_pT%d.pdf", FileDir, nState-3, irap, ipT);
				c3->SaveAs(plotname);

				c3 = new TCanvas(plotname, plotname, 1700, 500);
				c3->SetFillColor(kWhite);
				c3->SetTopMargin(0.05);
				c3->SetLeftMargin(0.13);
				c3->SetRightMargin(0.05);
				muonpT->Draw("e");
				sprintf(plotname, "%s/Psi%dS_muonpT_rap%d_pT%d.pdf", FileDir, nState-3, irap, ipT);
				c3->SaveAs(plotname);

				c3 = new TCanvas(plotname, plotname, 700, 500);
				c3->SetFillColor(kWhite);
				c3->SetTopMargin(0.05);
				c3->SetLeftMargin(0.13);
				c3->SetRightMargin(0.05);
				muoneta->Draw("e");
				sprintf(plotname, "%s/Psi%dS_muoneta_rap%d_pT%d.pdf", FileDir, nState-3, irap, ipT);
				c3->SaveAs(plotname);




				delete costhphiPX;
				delete costhPX;
				delete phiPX;

				delete inFile;
				delete c3;

			}
		}
	}

}
