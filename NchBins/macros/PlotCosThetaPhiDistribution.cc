/*
 * PlotCosThetaDistribution.cc
 *
 *  Created on: Dec 8, 2011
 *      Author: valentinknuenz
 *  Modified on: Jan 22, 2013 linlinzhang
 *
 */

#include <iostream>
#include <string>
#include <sstream>
using namespace std;

#include "rootIncludes.inc"
#include "commonVar.h"

// binning
int const pT_binsFOR1Dhists=onia::kNbPTMaxBins;
int const Rap_binsFOR1Dhists=onia::kNbRapForPTBins;
int nBins1DCosth=20;
int nBins1DPhi=18;
// extremes and binning of lambda_gen extraction histos
const double l_min = -1;
const double l_max =  1;
const double l_step_1D = 0.02;

// cosTheta - phi - Distributions
TH2D* CosThPhDist[onia::kNbFrames][onia::kNbRapForPTBins+1][onia::kNbPTMaxBins+1][onia::NchBins+1];
TH1D* CosThDist[onia::kNbFrames][onia::kNbRapForPTBins+1][onia::kNbPTMaxBins+1][onia::NchBins+1];
TH1D* PhiDist[onia::kNbFrames][onia::kNbRapForPTBins+1][onia::kNbPTMaxBins+1][onia::NchBins+1];

// extraction histos
TH1D* h_costh2[onia::kNbFrames][onia::kNbRapForPTBins+1][onia::kNbPTMaxBins+1][onia::NchBins+1];
TH1D* h_cos2ph[onia::kNbFrames][onia::kNbRapForPTBins+1][onia::kNbPTMaxBins+1][onia::NchBins+1];
TH1D* h_sin2thcosph[onia::kNbFrames][onia::kNbRapForPTBins+1][onia::kNbPTMaxBins+1][onia::NchBins+1];

// results
TGraphAsymmErrors* graph_lamth[onia::kNbFrames][onia::kNbRapForPTBins+1][onia::NchBins+1];
TGraphAsymmErrors* graph_lamph[onia::kNbFrames][onia::kNbRapForPTBins+1][onia::NchBins+1];
TGraphAsymmErrors* graph_lamtp[onia::kNbFrames][onia::kNbRapForPTBins+1][onia::NchBins+1];
TGraphAsymmErrors* graph_lamthstar[onia::kNbFrames][onia::kNbRapForPTBins+1][onia::NchBins+1];
TGraphAsymmErrors* graph_lamphstar[onia::kNbFrames][onia::kNbRapForPTBins+1][onia::NchBins+1];
TGraphAsymmErrors* graph_lamtilde[onia::kNbFrames][onia::kNbRapForPTBins+1][onia::NchBins+1];

void LoadHistos(int iRapBin, int iPTBin, int iCPMBin, int nState, Char_t *DataPath);
void PlotHistos(int iRapBin, int iPTBin, int iCPMBin, int iFrame, int nState, Char_t *DataPath);

//===========================

int main(int argc, char** argv){
	Char_t *DataPath = "DataPath_Default";
	int nState = 999;

	for( int i=0;i < argc; ++i ) {
		if(std::string(argv[i]).find("nState") != std::string::npos) {char* nStatechar = argv[i]; char* nStatechar2 = strtok (nStatechar, "p"); nState = atof(nStatechar2); cout<<"nState = "<<nState<<endl;}
		if(std::string(argv[i]).find("DataPath") != std::string::npos) {char* DataPathchar = argv[i]; char* DataPathchar2 = strtok (DataPathchar, "="); DataPath = DataPathchar2; cout<<"DataPath = "<<DataPath<<endl;}
	}

	gStyle->SetPalette(1);
	gStyle->SetFrameBorderMode(0);

	std::cout << DataPath << std::endl;

	int nBinsCT=25;
	int nBinsPH=25;

	// initialization of all graphs
	for(int iRap = 0; iRap < onia::kNbRapForPTBins; iRap++){
	  for(int iPT = 0; iPT < onia::kNbPTBins[iRap+1]; iPT++){
		for(int iFrame = 0; iFrame < onia::kNbFrames; iFrame++){

			graph_lamth[iFrame][iRap][iPT] = new TGraphAsymmErrors();
			graph_lamph[iFrame][iRap][iPT] = new TGraphAsymmErrors();
			graph_lamtp[iFrame][iRap][iPT] = new TGraphAsymmErrors();
			graph_lamthstar[iFrame][iRap][iPT] = new TGraphAsymmErrors();
			graph_lamphstar[iFrame][iRap][iPT] = new TGraphAsymmErrors();
			graph_lamtilde[iFrame][iRap][iPT] = new TGraphAsymmErrors();

		}

		for(int iCPM = 0; iCPM < onia::NchBins; iCPM++){
			for(int iFrame = 0; iFrame < onia::kNbFrames; iFrame++){

				Char_t name[100];
				sprintf(name, "%s/tmpFiles/data_Psi%dS_rap%d_pT%d_cpm%d.root", DataPath,nState-3, iRap+1, iPT+1, iCPM+1);
				TFile *fIn = new TFile(name);
				TH2D *hNPBG_cosThetaPhi_PHX = (TH2D*)fIn->Get("hNPBG_cosThetaPhi_PHX");
				nBinsCT = hNPBG_cosThetaPhi_PHX->GetNbinsX();
				nBinsPH = hNPBG_cosThetaPhi_PHX->GetNbinsY();
				nBins1DCosth = nBinsCT; nBins1DPhi = nBinsPH;
				CosThPhDist[iFrame][iRap][iPT][iCPM] = new TH2D("","",nBinsCT, onia::cosTMin, onia::cosTMax,
						nBinsPH, onia::phiPolMin, onia::phiPolMax);
				CosThDist[iFrame][iRap][iPT][iCPM] = new TH1D("","",nBins1DCosth, onia::cosTMin, onia::cosTMax);
				PhiDist[iFrame][iRap][iPT][iCPM] = new TH1D("","",nBins1DPhi, onia::phiPolMin, onia::phiPolMax);

				h_costh2[iFrame][iRap][iPT][iCPM] = new TH1D( "", "", int((l_max-l_min)/l_step_1D), l_min, l_max );
				h_cos2ph[iFrame][iRap][iPT][iCPM] = new TH1D( "", "", int((l_max-l_min)/l_step_1D), l_min, l_max );
				h_sin2thcosph[iFrame][iRap][iPT][iCPM] = new TH1D( "", "", int((l_max-l_min)/l_step_1D), l_min, l_max );

			} //iFrame
			LoadHistos(iRap, iPT, iCPM, nState,DataPath);
			for(int iFrame = 0; iFrame < 3; iFrame++){
				PlotHistos(iRap, iPT, iCPM, iFrame,nState,DataPath);
			}

		} //iCPM
	  } //iPT
	} //iRap

	char filename[200];
	sprintf(filename, "%s/Figures/TGraphResults_gen_Psi%dS.root", DataPath, nState-3);
	TFile* results = new TFile(filename, "RECREATE");

	sprintf(filename, "%s/Figures/AngDistHist.root", DataPath);
	TFile* AngDistHistFile = new TFile(filename, "RECREATE", "AngDistHistFile");

	for(int iRap = 0; iRap < onia::kNbRapForPTBins; iRap++){
	  for(int iPT = 0; iPT < onia::kNbPTBins[iRap+1]; iPT++){
		for(int iFrame = 0; iFrame < 3; iFrame++){
			for(int iCPM = 0; iCPM < onia::NchBins; iCPM++){

				AngDistHistFile->cd();
				char histName[200];               
				// write histograms to file
				sprintf(histName,"Proj_%s_costh_rap%d_pT%d_cpm%d",onia::frameLabel[iFrame],iRap+1,iPT+1,iCPM+1);
				CosThDist[iFrame][iRap][iPT][iCPM]->SetName(histName);
				CosThDist[iFrame][iRap][iPT][iCPM]->Write();
				sprintf(histName,"Proj_%s_phi_rap%d_pT%d_cpm%d",onia::frameLabel[iFrame],iRap+1,iPT+1,iCPM+1);
				PhiDist[iFrame][iRap][iPT][iCPM]->SetName(histName);
				PhiDist[iFrame][iRap][iPT][iCPM]->Write();
				sprintf(histName,"CosthPhi_%s_rap%d_pT%d_cpm%d",onia::frameLabel[iFrame],iRap+1,iPT+1,iCPM+1);
				CosThPhDist[iFrame][iRap][iPT][iCPM]->SetName(histName);
				CosThPhDist[iFrame][iRap][iPT][iCPM]->Write();

				sprintf(histName,"costh2_%s_rap%d_pT%d_cpm%d",onia::frameLabel[iFrame],iRap+1,iPT+1,iCPM+1);
				h_costh2[iFrame][iRap][iPT][iCPM]->SetName(histName);
				h_costh2[iFrame][iRap][iPT][iCPM]->Write();
				sprintf(histName,"cos2ph_%s_rap%d_pT%d_cpm%d",onia::frameLabel[iFrame],iRap+1,iPT+1,iCPM+1);
				h_cos2ph[iFrame][iRap][iPT][iCPM]->SetName(histName);
				h_cos2ph[iFrame][iRap][iPT][iCPM]->Write();
				sprintf(histName,"sin2thcosph_%s_rap%d_pT%d_cpm%d",onia::frameLabel[iFrame],iRap+1,iPT+1,iCPM+1);
				h_sin2thcosph[iFrame][iRap][iPT][iCPM]->SetName(histName);
				h_sin2thcosph[iFrame][iRap][iPT][iCPM]->Write();

				// calculate the polarization parameters
				double costh2 = h_costh2[iFrame][iRap][iPT][iCPM]->GetMean();
				double lamth = (1. - 3. * costh2 ) / ( costh2 - 3./5. );
				double cos2ph = h_cos2ph[iFrame][iRap][iPT][iCPM]->GetMean();
				double lamph = cos2ph * (3. + lamth);
				double sin2thcosph = h_sin2thcosph[iFrame][iRap][iPT][iCPM]->GetMean();
				double lamtp = sin2thcosph * 5./4. * (3. + lamth);
				double lamtilde = (lamth + 3*lamph)/(1-lamph);
				double lamthstar = -9;
				double lamphstar = -9;

				// calculate mean pt of bin
				double meanpt = onia::pTRange[iRap+1][iPT] + (onia::pTRange[iRap+1][iPT+1] - onia::pTRange[iRap+1][iPT])/2; 
				double meancpm = onia::cpmRange[iPT+1] + (onia::cpmRange[iCPM+1] - onia::cpmRange[iCPM])/2; 

				// save to TGraph
				graph_lamth[iFrame][iRap][iPT]->SetPoint(iCPM, meancpm, lamth);
				graph_lamph[iFrame][iRap][iPT]->SetPoint(iCPM, meancpm, lamph);
				graph_lamtp[iFrame][iRap][iPT]->SetPoint(iCPM, meancpm, lamtp);
				graph_lamtilde[iFrame][iRap][iPT]->SetPoint(iCPM, meancpm, lamtilde);
				graph_lamthstar[iFrame][iRap][iPT]->SetPoint(iCPM, meancpm, lamthstar);
				graph_lamphstar[iFrame][iRap][iPT]->SetPoint(iCPM, meancpm, lamphstar);

			} // iCPM

			results->cd();
			char graphName[200];
			sprintf(graphName,"lth_%s_rap%d_pt%d", onia::frameLabel[iFrame],iRap+1,iPT+1);
			if(iFrame==2) sprintf(graphName,"lth_PX_rap%d_pt%d", iRap+1,iPT+1);
			graph_lamth[iFrame][iRap][iPT]->SetName(graphName);            
			graph_lamth[iFrame][iRap][iPT]->Write();
			sprintf(graphName,"lph_%s_rap%d_pt%d", onia::frameLabel[iFrame],iRap+1,iPT+1);
			if(iFrame==2) sprintf(graphName,"lph_PX_rap%d_pt%d", iRap+1,iPT+1);
			graph_lamph[iFrame][iRap][iPT]->SetName(graphName);
			graph_lamph[iFrame][iRap][iPT]->Write();
			sprintf(graphName,"ltp_%s_rap%d_pt%d", onia::frameLabel[iFrame],iRap+1,iPT+1);
			if(iFrame==2) sprintf(graphName,"ltp_PX_rap%d_pt%d", iRap+1,iPT+1);
			graph_lamtp[iFrame][iRap][iPT]->SetName(graphName);
			graph_lamtp[iFrame][iRap][iPT]->Write();
			sprintf(graphName,"ltilde_%s_rap%d_pt%d", onia::frameLabel[iFrame],iRap+1,iPT+1);
			if(iFrame==2) sprintf(graphName,"ltilde_PX_rap%d_pt%d", iRap+1,iPT+1);
			graph_lamtilde[iFrame][iRap][iPT]->SetName(graphName);
			graph_lamtilde[iFrame][iRap][iPT]->Write();
			sprintf(graphName,"lthstar_%s_rap%d_pt%d", onia::frameLabel[iFrame],iRap+1,iPT+1);
			if(iFrame==2) sprintf(graphName,"lthstar_PX_rap%d_pt%d", iRap+1,iPT+1);
			graph_lamthstar[iFrame][iRap][iPT]->SetName(graphName);
			graph_lamthstar[iFrame][iRap][iPT]->Write();
			sprintf(graphName,"lphstar_%s_rap%d_pt%d", onia::frameLabel[iFrame],iRap+1,iPT+1);
			if(iFrame==2) sprintf(graphName,"lphstar_PX_rap%d_pt%d", iRap+1,iPT+1);
			graph_lamphstar[iFrame][iRap][iPT]->SetName(graphName);
			graph_lamphstar[iFrame][iRap][iPT]->Write();

		} // iFrame
	  }// iPT
	} // iRap

	results->Write();
	results->Close();
	AngDistHistFile->Write();
	AngDistHistFile->Close();

	return 0;
}
//===========================
void PlotHistos(int iRapBin, int iPTBin, int iCPMBin, int iFrame,int nState, Char_t *DataPath){

	TGaxis::SetMaxDigits(3);

	double lvalue = 0.22, tvalue = 0.92;
	double left=lvalue, top=tvalue, textSize=0.035;
	TLatex *latex=new TLatex();
	latex->SetTextFont(42);
	latex->SetNDC(kTRUE);
	latex->SetTextSize(textSize);
	double step=textSize*1.3;

	Char_t name[500];
	sprintf(name, "c1_%s_rap%d_pT%d_cpm%d", onia::frameLabel[iFrame], iRapBin+1, iPTBin+1, iCPMBin+1);
	TCanvas *c1 = new TCanvas(name, "", 500, 500);
	gStyle->SetPalette(1);
	gStyle->SetOptStat(0);
	gPad->SetFillColor(kWhite);
	gPad->SetLeftMargin(0.15);
	gPad->SetRightMargin(0.15);

	char XTitle[200];
	char YTitle[200];
	sprintf(XTitle,"cos#theta_{%s}",onia::frameLabel[iFrame]);
	sprintf(YTitle,"#phi_{%s} [deg]",onia::frameLabel[iFrame]);
	double yOffset=1.4;
	CosThPhDist[iFrame][iRapBin][iPTBin][iCPMBin]->GetYaxis()->SetTitleOffset(yOffset);
	//gPad->SetLeftMargin(0.125);
	CosThPhDist[iFrame][iRapBin][iPTBin][iCPMBin]->SetStats(0);
	CosThPhDist[iFrame][iRapBin][iPTBin][iCPMBin]->GetYaxis()->SetTitle(YTitle);
	CosThPhDist[iFrame][iRapBin][iPTBin][iCPMBin]->GetXaxis()->SetTitle(XTitle);

	CosThPhDist[iFrame][iRapBin][iPTBin][iCPMBin]->Draw("colz");
	sprintf(name, "%s/Figures/cosThetaPhi_%s_rap%d_pT%d_cpm%d.pdf", DataPath,onia::frameLabel[iFrame], iRapBin+1, iPTBin+1, iCPMBin+1);
	cout<<name<<endl;

	if(iRapBin==0)
		latex->DrawLatex(left,top,Form("|y| < %.1f, %.1f < p_{T} < %.1f GeV,\n %.1f < N_{ch} < %.1f",
					onia::rapForPTRange[iRapBin+1],
					onia::pTRange[iRapBin+1][iPTBin],onia::pTRange[iRapBin+1][iPTBin+1],
					onia::cpmRange[iCPMBin],onia::cpmRange[iCPMBin+1]));
	else
		latex->DrawLatex(left,top,Form("%.1f < |y| < %.1f, %.1f < p_{T} < %.1f GeV,\n %.1f < N_{ch} < %.1f",
					onia::rapForPTRange[iRapBin],onia::rapForPTRange[iRapBin+1],
					onia::pTRange[iRapBin+1][iPTBin],onia::pTRange[iRapBin+1][iPTBin+1],
					onia::cpmRange[iCPMBin],onia::cpmRange[iCPMBin+1]));
	c1->Print(name);
}

//===========================
void LoadHistos(int iRapBin, int iPTBin, int iCPMBin, int nState, Char_t *DataPath){
	Char_t name[100];
	sprintf(name, "%s/tmpFiles/data_Psi%dS_rap%d_pT%d_cpm%d.root", DataPath,nState-3, iRapBin+1, iPTBin+1, iCPMBin+1);

	TFile *fIn = new TFile(name);
	TTree* selectedData = (TTree*)fIn->Get("selectedData");

	TLorentzVector *lepP = new TLorentzVector();
	TLorentzVector *lepN = new TLorentzVector();
	selectedData->SetBranchAddress("lepP", &lepP);
	selectedData->SetBranchAddress("lepN", &lepN);

	const double pbeam_ = 7000.; // exact number irrelevant as long as pbeam >> Mprot
	const double Mprot_ = 0.9382720;
	const double Mlepton_ = 0.10566;  // (muon)
	const double gPI_ = TMath::Pi();
	const double Ebeam_ = sqrt( pbeam_*pbeam_ + Mprot_*Mprot_ );
	TLorentzVector beam1_LAB_( 0., 0., pbeam_, Ebeam_ );
	TLorentzVector beam2_LAB_( 0., 0., -pbeam_, Ebeam_ );

	for(int i=0;i<selectedData->GetEntries();i++){

		selectedData->GetEvent( i );

		double lepP_pT  = lepP->Pt();
		double lepN_pT  = lepN->Pt();

		double lepP_eta = lepP->PseudoRapidity();
		double lepN_eta = lepN->PseudoRapidity();

		// dilepton 4-vector:
		TLorentzVector dilepton = *lepP + *lepN;
		double pT   = dilepton.Pt();
		double rap  = dilepton.Rapidity();
		double mass = dilepton.M();

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

		double cosalpha = sqrt( 1. - pow(costh_PX, 2.) ) * sin( lepton_DILEP_rotated.Phi() );

		// Filling Histograms of costh2, cos2ph and sin2thcosph for the extraction of the actual generated polarization
		// CS frame:
		double costh2 = pow(costh_CS,2.);
		double Phi = phi_CS/180. * gPI_;
		double cos2ph = cos(2.*Phi);
		double sin2thcosph= sin(2.*acos(costh_CS))*cos(Phi);
		h_costh2[0][iRapBin][iPTBin][iCPMBin]->Fill( costh2 );
		h_cos2ph[0][iRapBin][iPTBin][iCPMBin]->Fill( cos2ph );
		h_sin2thcosph[0][iRapBin][iPTBin][iCPMBin]->Fill( sin2thcosph );

		// HX frame:
		costh2 = pow(costh_HX,2.);
		Phi = phi_HX/180. * gPI_;
		cos2ph = cos(2.*Phi);
		sin2thcosph= sin(2.*acos(costh_HX))*cos(Phi);
		h_costh2[1][iRapBin][iPTBin][iCPMBin]->Fill( costh2 );
		h_cos2ph[1][iRapBin][iPTBin][iCPMBin]->Fill( cos2ph );
		h_sin2thcosph[1][iRapBin][iPTBin][iCPMBin]->Fill( sin2thcosph );

		// PHX frame:
		costh2 = pow(costh_PX,2.);
		Phi = phi_PX/180. * gPI_;
		cos2ph = cos(2.*Phi);
		sin2thcosph= sin(2.*acos(costh_PX))*cos(Phi);
		h_costh2[2][iRapBin][iPTBin][iCPMBin]->Fill( costh2 );
		h_cos2ph[2][iRapBin][iPTBin][iCPMBin]->Fill( cos2ph );
		h_sin2thcosph[2][iRapBin][iPTBin][iCPMBin]->Fill( sin2thcosph );

		// fill cosTheta-phi-distributions in different frames
		CosThPhDist[0][iRapBin][iPTBin][iCPMBin]->Fill(costh_CS,phi_CS);
		CosThPhDist[1][iRapBin][iPTBin][iCPMBin]->Fill(costh_HX,phi_HX);
		CosThPhDist[2][iRapBin][iPTBin][iCPMBin]->Fill(costh_PX,phi_PX);

		CosThDist[0][iRapBin][iPTBin][iCPMBin]->Fill(costh_CS);
		PhiDist[0][iRapBin][iPTBin][iCPMBin]->Fill(phi_CS);
		CosThDist[1][iRapBin][iPTBin][iCPMBin]->Fill(costh_HX);
		PhiDist[1][iRapBin][iPTBin][iCPMBin]->Fill(phi_HX);
		CosThDist[2][iRapBin][iPTBin][iCPMBin]->Fill(costh_PX);
		PhiDist[2][iRapBin][iPTBin][iCPMBin]->Fill(phi_PX);

	}

	//  cout<<"loaded"<<endl;
}
