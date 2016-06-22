#define PolData_cxx
#include "PolData.h"
#include "commonVar.h"
#include "effsAndCuts.h"

#include "TLorentzVector.h"
#include <TH2.h>
#include <TCanvas.h>

#include <string>
#include <iostream>
#include <sstream>

TH1F *Reco_StatEv;
TH1F *Reco_Onia_mass[onia::kNbPTMaxBins+1][onia::kNbRapForPTBins+1];
TH2F *Reco_Onia_rap_pT;
TH1F *Reco_Onia_pt[onia::kNbRapForPTBins+1];
TH1F *Reco_Onia_rap[onia::kNbPTMaxBins+1];

TTree *treeOut;
TLorentzVector *lepP, *lepN, *jpsi;


void PolData::Loop(int nState, bool rejectCowboys, int FidCuts, bool MC, bool RequestTrigger, bool removeEta0p2_0p3, bool cutDeltaREllDpt) {

	if (fChain == 0) return;

	Long64_t nentries = fChain->GetEntries();
	Long64_t cutAtRecEvent = nentries;
	Long64_t count = 0;
	Long64_t nb = 0;

	std::cout << "number of entries = " << nentries << std::endl;

	// create branches for certain output variables
	double jpsict = 0;
	double jpsictErr = 0;
	double jpsiMassErr = 0;
	double jpsiVprob = 0;

	treeOut->Branch("Jpsict", &jpsict, "Jpsict/D");
	treeOut->Branch("JpsictErr", &jpsictErr, "JpsictErr/D");
	treeOut->Branch("JpsiMassErr", &jpsiMassErr, "JpsiMassErr/D");
	treeOut->Branch("JpsiVprob", &jpsiVprob, "JpsiVprob/D");

  //double rndNumber;
  //ifstream rndFile; rndFile.open("/afs/ihep.ac.cn/users/z/zhangll/fs/work/polarization/PsiPol2011/macros/random.txt");


//	nentries = 2500000;
	//loop over the events
	for (Long64_t jentry=0; jentry<nentries; jentry++) {

	//for (Long64_t jentry=0; jentry<nentries/2; jentry++) {  // first half
	//for (Long64_t jentry=nentries/2; jentry<nentries; jentry++) {  // second half
	
	//for (Long64_t jentry=0; jentry<nentries/4; jentry++) {  //  one quarter
	//for (Long64_t jentry=nentries/4; jentry<2*nentries/4; jentry++) {  //  two quarter
	//for (Long64_t jentry=2*nentries/4; jentry<3*nentries/4; jentry++) {  //  three quarter
	//for (Long64_t jentry=3*nentries/4; jentry<nentries; jentry++) {  //  four quarter


		if(jentry % 100000 == 0) std::cout << "event " << jentry << " of " << nentries << std::endl;

		Long64_t ientry = LoadTree(jentry);
		if (ientry < 0) break;
		nb = fChain->GetEntry(jentry);

		//if( !(onia->Rapidity() < 0.) ) continue;
		//if( !(onia->Rapidity() > 0.) ) continue;
		//if( !(onia->Rapidity() > -0.6 && onia->Rapidity() < -0.3) ) continue;
		//if( !(onia->Rapidity() > -0.3 && onia->Rapidity() < 0.) ) continue;
		//if( !(onia->Rapidity() > 0. && onia->Rapidity() < 0.3) ) continue;
		//if( !(onia->Rapidity() > 0.3 && onia->Rapidity() < 0.6) ) continue;
		
		//rndFile>>rndNumber;
		//cout << "rndNumber " << rndNumber << endl;
		//if( !(rndNumber > 0.5) )  continue; // rndUp0p5
		//if( !(rndNumber < 0.5) )  continue; // rndLow0p5
		

		//if we process MC, we must ensure that we only consider reconstructed events
		if(onia->Pt() > 990.) continue;
		if(JpsiVprob < 0.01) continue;

		// count all events
		Reco_StatEv->Fill(0.5);

		//check the trigger flag: 0... no trigger, 1 ... triggered+matched, 3 ... triggered (HLT_DoubleMu0)
		//for a full list of accessible triggers, check https://espace.cern.ch/cms-quarkonia/onia-polarization/L1%20%20HLT/unprescaledTriggersVsRun.aspx
		int trigDecision = -99;

		// different trigger for different particles
		// Jpsi trigger paths:
		if(nState == 4){
			if(
					//HLT_Dimuon0_Jpsi_v1 == 1 || //  1e33: 165088 - 166967 (prescale of 20) and 1.4E33: 167039 - 167043 (prescale of 60)
					//HLT_Dimuon0_Jpsi_v2 == 1 || //  1e33: 166346 (prescale of 40)
					//HLT_Dimuon0_Jpsi_v3 == 1 ||  //1.4e33: 167078 - 167913 (prescale of 120)
					//HLT_Dimuon0_Jpsi_v5 == 1 || //2E33 (no cowboys)
					//HLT_Dimuon0_Jpsi_v6 == 1 || //3E33 (prescale of 120) (L1_DoubleMu0_HighQ; becomes inactive for Linst >= 5E33)
					//HLT_Dimuon0_Jpsi_v9 == 1 || //5E33 (prescale of 200)
					HLT_Dimuon10_Jpsi_Barrel_v1 == 1 || //  1e33: 165088 - 166967 and 1.4E33: 167039 - 167043 
					HLT_Dimuon10_Jpsi_Barrel_v2 == 1 || //  1e33: 166346
					HLT_Dimuon10_Jpsi_Barrel_v3 == 1 ||  //1.4e33: 167078 - 167913 (prescale of 2)
					HLT_Dimuon10_Jpsi_Barrel_v5 == 1 || //2E33 (no cowboys)
					HLT_Dimuon10_Jpsi_Barrel_v6 == 1 || //3E33 (L1_DoubleMu0_HighQ; becomes inactive for Linst >= 5E33)
					HLT_Dimuon10_Jpsi_Barrel_v9 == 1 ||
					HLT_Dimuon8_Jpsi_v3 == 1  ||
					HLT_Dimuon8_Jpsi_v4 == 1  ||
					HLT_Dimuon8_Jpsi_v5 == 1  ||
					HLT_Dimuon8_Jpsi_v6 == 1  ||
					HLT_Dimuon8_Jpsi_v7 == 1  ||
					HLT_Dimuon10_Jpsi_v3 == 1 ||
					HLT_Dimuon10_Jpsi_v4 == 1 ||
					HLT_Dimuon10_Jpsi_v6 == 1 ||
					HLT_Mu15_tkMu5_Onia_v1 == 1)  
					//5E33
				//HLT_Dimuon13_Jpsi_Barrel_v1 == 1 || //3E33 (L1_DoubleMu0_HighQ; becomes inactive for Linst >= 5E33)
				//HLT_Dimuon13_Jpsi_Barrel_v4 == 1) //5E33
				trigDecision = 1;

			if(trigDecision != 1 && RequestTrigger) continue;

			if(
					(HLT_Dimuon10_Jpsi_Barrel_v1 == 1 && onia->Pt() < 10) ||
					(HLT_Dimuon10_Jpsi_Barrel_v2 == 1 && onia->Pt() < 10) ||
					(HLT_Dimuon10_Jpsi_Barrel_v3 == 1 && onia->Pt() < 10) ||
					(HLT_Dimuon10_Jpsi_Barrel_v5 == 1 && onia->Pt() < 10) ||
					(HLT_Dimuon10_Jpsi_Barrel_v6 == 1 && onia->Pt() < 10) ||
					(HLT_Dimuon10_Jpsi_Barrel_v9 == 1 && onia->Pt() < 10) ||
					(HLT_Dimuon8_Jpsi_v3 == 1 && onia->Pt() < 10) ||
					(HLT_Dimuon8_Jpsi_v4 == 1 && onia->Pt() < 10) ||
					(HLT_Dimuon8_Jpsi_v5 == 1 && onia->Pt() < 10) ||
					(HLT_Dimuon8_Jpsi_v6 == 1 && onia->Pt() < 10) ||
					(HLT_Dimuon8_Jpsi_v7 == 1 && onia->Pt() < 10) ||
					(HLT_Dimuon10_Jpsi_v3 == 1 && onia->Pt() < 10) ||
					(HLT_Dimuon10_Jpsi_v4 == 1 && onia->Pt() < 10) ||
					(HLT_Dimuon10_Jpsi_v6 == 1 && onia->Pt() < 10) ||
					(HLT_Mu15_tkMu5_Onia_v1 == 1 && onia->Pt() < 10) 
					) 
				//(HLT_Dimuon13_Jpsi_Barrel_v1 == 1 && onia->Pt() < 13) ||
				//(HLT_Dimuon13_Jpsi_Barrel_v4 == 1 && onia->Pt() < 13))
				continue;

		} // if Jpsi

		// PsiPrime trigger paths
		if(nState == 5){
			if(
					HLT_Dimuon7_PsiPrime_v1 == 1 || //  1e33: 165088 - 166967 and 1.4E33: 167039 - 167043
					HLT_Dimuon7_PsiPrime_v2 == 1 || //  1e33: 166346
					HLT_Dimuon7_PsiPrime_v3 == 1 ||  //1.4e33: 167078 - 167913 (prescale of 2)
					HLT_Dimuon7_PsiPrime_v5 == 1 || //2E33 (no cowboys)
					HLT_Dimuon9_PsiPrime_v1 == 1 || //3E33 (L1_DoubleMu0_HighQ; becomes inactive for Linst >= 5E33)
					HLT_Dimuon9_PsiPrime_v4 == 1 ||
					HLT_Dimuon5_PsiPrime_v3 == 1 ||
					HLT_Dimuon5_PsiPrime_v4 == 1 ||
					HLT_Dimuon5_PsiPrime_v5 == 1 ||
					HLT_Dimuon5_PsiPrime_v6 == 1 ||
					HLT_Dimuon9_PsiPrime_v9 == 1 ||
					HLT_Mu15_tkMu5_Onia_v1 == 1)  //5E33
				//HLT_Dimuon11_PsiPrime_v1 == 1 || //3E33 (L1_DoubleMu0_HighQ; becomes inactive for Linst >= 5E33)
				//HLT_Dimuon11_PsiPrime_v4 == 1) //5E33
				trigDecision = 1;

			if(trigDecision != 1 && RequestTrigger) continue;

			if(
					(HLT_Dimuon7_PsiPrime_v1 == 1 && onia->Pt() < 7) ||
					(HLT_Dimuon7_PsiPrime_v2 == 1 && onia->Pt() < 7) ||
					(HLT_Dimuon7_PsiPrime_v3 == 1 && onia->Pt() < 7) ||
					(HLT_Dimuon7_PsiPrime_v5 == 1 && onia->Pt() < 7) ||
					(HLT_Dimuon9_PsiPrime_v1 == 1 && onia->Pt() < 9) ||
					(HLT_Dimuon9_PsiPrime_v4 == 1 && onia->Pt() < 9) ||
					(HLT_Dimuon9_PsiPrime_v9 == 1 && onia->Pt() < 9) ||
					(HLT_Dimuon11_PsiPrime_v1 == 1 && onia->Pt() < 11) ||
					(HLT_Dimuon5_PsiPrime_v3 == 1 && onia->Pt() < 9) ||
					(HLT_Dimuon5_PsiPrime_v4 == 1 && onia->Pt() < 9) ||
					(HLT_Dimuon5_PsiPrime_v5 == 1 && onia->Pt() < 9) ||
					(HLT_Dimuon5_PsiPrime_v6 == 1 && onia->Pt() < 9) ||
					(HLT_Mu15_tkMu5_Onia_v1 == 1 && onia->Pt() < 9) ||
					(HLT_Dimuon11_PsiPrime_v4 == 1 && onia->Pt() < 11))
				continue;

		} // if PsiPrime

		// Upsilon trigger paths:
		if(nState < 4){
			if(
					HLT_Dimuon5_Upsilon_Barrel_v1 == 1 || //  1e33: 165088 - 166967 and 1.4E33: 167039 - 167043
					HLT_Dimuon5_Upsilon_Barrel_v2 == 1 || //  1e33: 166346
					HLT_Dimuon5_Upsilon_Barrel_v3 == 1 ||  //1.4e33: 167078 - 167913 (prescale of 2)
					HLT_Dimuon5_Upsilon_Barrel_v5 == 1 || //2E33 (no cowboys)
					HLT_Dimuon7_Upsilon_Barrel_v1 == 1 || //3E33 (L1_DoubleMu0_HighQ; becomes inactive for Linst >= 5E33)
					HLT_Dimuon9_Upsilon_Barrel_v1 == 1 || //3E33 (L1_DoubleMu0_HighQ)
					HLT_Dimuon7_Upsilon_Barrel_v4 == 1 || //5E33 (becomes inactive for Linst >= 5E33)
					HLT_Dimuon9_Upsilon_Barrel_v4 == 1) //5E33
				trigDecision = 1;

			if(trigDecision != 1 && RequestTrigger) continue;

			if(
					HLT_Dimuon5_Upsilon_Barrel_v1 == 1 && onia->Pt() < 5.5  || //  1e33: 165088 - 166967 and 1.4E33: 167039 - 167043
					HLT_Dimuon5_Upsilon_Barrel_v2 == 1 && onia->Pt() < 5.5  || //  1e33: 166346
					HLT_Dimuon5_Upsilon_Barrel_v3 == 1 && onia->Pt() < 5.5  ||  //1.4e33: 167078 - 167913 (prescale of 2)
					HLT_Dimuon5_Upsilon_Barrel_v5 == 1 && onia->Pt() < 5.5  || //2E33 (no cowboys)
					HLT_Dimuon7_Upsilon_Barrel_v1 == 1 && onia->Pt() < 7.5  || //3E33 (L1_DoubleMu0_HighQ; becomes inactive for Linst >= 5E33)
					HLT_Dimuon9_Upsilon_Barrel_v1 == 1 && onia->Pt() < 9.5  || //3E33 (L1_DoubleMu0_HighQ)
					HLT_Dimuon7_Upsilon_Barrel_v4 == 1 && onia->Pt() < 7.5  || //5E33 (becomes inactive for Linst >= 5E33)
					HLT_Dimuon9_Upsilon_Barrel_v4 == 1 && onia->Pt() < 9.5 ) //5E33
				continue;

		} // if Upsilon

		// count events after trigger
		Reco_StatEv->Fill(1.5); 

		double onia_mass = onia->M();
		double onia_pt = onia->Pt();
		double onia_P = onia->P();
		double onia_eta = onia->PseudoRapidity();
		double onia_rap = onia->Rapidity();
		double onia_phi = onia->Phi();
		double onia_mT = sqrt(onia_mass*onia_mass + onia_pt*onia_pt);

		// restrict to barrel for upsilon and jpsi
		if(TMath::Abs(onia_rap) > onia::rap) continue;
		// count events after cut in rapidity
		Reco_StatEv->Fill(2.5); 

		double etaMuPos = muPos->PseudoRapidity();
		double etaMuNeg = muNeg->PseudoRapidity();
		double pTMuPos = muPos->Pt();
		double pTMuNeg = muNeg->Pt();
		double pMuPos = muPos->P();
		double pMuNeg = muNeg->P();
		double deltaPhi = muNeg->Phi() - muPos->Phi();
		if(deltaPhi > TMath::Pi()) deltaPhi -= 2.*TMath::Pi();
		else if(deltaPhi < -TMath::Pi()) deltaPhi += 2.*TMath::Pi();

		if(rejectCowboys)
			if(deltaPhi < 0.)  continue;
		// count events after cowboy rejection
		Reco_StatEv->Fill(3.5);


		// select mass window
		double massMin = 0, massMax = 0;
		// different mass region for different states
		massMin = onia::massMin; 
		massMax = onia::massMax;
		if(onia_mass < massMin || onia_mass > massMax) continue;   
		// count events after restriction in mass
		Reco_StatEv->Fill(4.5);

		//apply fiducial cuts
		bool muonsInAcc = kFALSE;
		if(isMuonInAcceptance(FidCuts-1, pTMuPos, etaMuPos) && isMuonInAcceptance(FidCuts-1, pTMuNeg, etaMuNeg)){
			muonsInAcc = kTRUE;
		}

		if(!muonsInAcc) continue;
		// count events after fiducial cuts
		Reco_StatEv->Fill(5.5);

		if(removeEta0p2_0p3)
			if( (TMath::Abs(etaMuPos) > 0.2 && TMath::Abs(etaMuPos) < 0.3) || 
					(TMath::Abs(etaMuNeg) > 0.2 && TMath::Abs(etaMuNeg) < 0.3) )
				continue;

		Reco_StatEv->Fill(6.5);


		double deltaEta = etaMuNeg - etaMuPos;
		double deltaPT = pTMuNeg - pTMuPos;
		double deltaREll = TMath::Sqrt(deltaEta*deltaEta + TMath::Power(1.2*deltaPhi,2));
		double deltaREllDpt = TMath::Sqrt(deltaREll*deltaREll + TMath::Power(0.00157*deltaPT,2));

		if(cutDeltaREllDpt){
			if(onia_pt>35. && onia_pt<40. && deltaREllDpt<0.18) continue;
			if(onia_pt>40. && onia_pt<50. && deltaREllDpt<0.16) continue;
			if(onia_pt>50. && onia_pt<70. && deltaREllDpt<0.14) continue;
		}

		Reco_StatEv->Fill(7.5);

		// filling mass, pt and y histograms
		// indices [0] contain all events while [1], etc. show events according to the pt and y bin

		// all events
		Reco_Onia_rap_pT->Fill(onia_rap, onia_pt);
		Reco_Onia_mass[0][0]->Fill(onia_mass);
		Reco_Onia_rap[0]->Fill(onia_rap);
		Reco_Onia_pt[0]->Fill(onia_pt);

		for(int iRap = 0; iRap < onia::kNbRapForPTBins; iRap++){
			// events integrated in pt in different rapidity bins
			if(TMath::Abs(onia_rap) > onia::rapForPTRange[iRap] && TMath::Abs(onia_rap) < onia::rapForPTRange[iRap+1]){
				Reco_Onia_pt[iRap+1]->Fill(onia_pt);
				Reco_Onia_mass[0][iRap+1]->Fill(onia_mass);
			}//if

			for(int iPT = 0; iPT < onia::kNbPTMaxBins; iPT++){		
				// events for different pT and y bins
				if(onia_pt > onia::pTRange[0][iPT] && onia_pt < onia::pTRange[0][iPT+1] &&
						TMath::Abs(onia_rap) > onia::rapForPTRange[iRap] && TMath::Abs(onia_rap) < onia::rapForPTRange[iRap+1]){
					Reco_Onia_mass[iPT+1][iRap+1]->Fill(onia_mass);
				}//if
			} // iPT
		} // iRap

		for(int iPT = 0; iPT < onia::kNbPTMaxBins; iPT++){	
			// events integrated in rapidity for different pT bins
			if(onia_pt > onia::pTRange[0][iPT] && onia_pt < onia::pTRange[0][iPT+1]){
				Reco_Onia_rap[iPT+1]->Fill(onia_rap);
				Reco_Onia_mass[iPT+1][0]->Fill(onia_mass); 
			} // if
		} // iPT

		// fill tree
		lepP = muPos;
		lepN = muNeg;
		jpsi = onia;
		jpsict = Jpsict;
		jpsictErr = JpsictErr;
		jpsiMassErr = JpsiMassErr;
		jpsiVprob = JpsiVprob;
		treeOut->Fill();

		//remaining of the events will be used for the analysis
		count++;

	} // for loop over events

	std::cout << "number of reconstructed events: " << count << " of a total of " << nentries << " events" << std::endl;

} // void
