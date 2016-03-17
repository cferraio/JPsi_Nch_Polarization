#include "Riostream.h"
#include "TSystem.h"
#include "TString.h"
#include "TFile.h"
#include "TROOT.h"
#include "TTree.h"
#include "TMath.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "TRandom3.h"
#include "TF1.h"
#include "TF2.h"
#include "TH1.h"
#include "TH2.h"
#include "TSystem.h"

//#include "genDefs.h"

bool HX_is_natural_sig = false;  // put both to false to generate in the CS frame
bool PX_is_natural_sig = false;

bool HX_is_natural_bkg = false;
bool PX_is_natural_bkg = false;



const double pbeam = 7000;
const double Mprot = 0.9382720;
const double Mlepton = 0.10566;  // (muon)
double gPI = TMath::Pi();
double Ebeam = sqrt( pbeam*pbeam + Mprot*Mprot );

TLorentzVector beam1_LAB( 0., 0., pbeam, Ebeam );
TLorentzVector beam2_LAB( 0., 0., -pbeam, Ebeam );


double func_rap_gen(double* x, double* par) {
	return   1.;
}

double func_pT_gen(double* x, double* par) {
	double beta = 3.; 
	double pTsq = 49.;
	if(par[0]==1) {beta = 3.46; pTsq = 47.3;}  // Upsi(1S)
	if(par[0]==2) {beta = 3.27; pTsq = 65.7;}  // Upsi(2S)
	if(par[0]==3) {beta = 3.05; pTsq = 80.5;}  // Upsi(3S)
	if(par[0]==4) {beta = 3.69; pTsq = 12.0;}  // Psi(1S)
	if(par[0]==5) {beta = 3.71; pTsq = 19.5;}  // Psi(2S)

	return x[0] * pow( 1. + 1./(beta - 2.) * x[0]*x[0] / pTsq, -beta  );
}


void polGen(double rapdilepton_min = 1,
		double rapdilepton_max = 1,
		double pTdilepton_min = 1,
		double pTdilepton_max = 1,
		double mass_signal_peak  =  1,
		double mass_signal_sigma =  1,
		double n_sigmas_signal = 1,
		int n_events = 50000,
		double f_BG = 0.5,
		double lambda_theta_sig=1,
		double lambda_phi_sig=1,
		double lambda_thetaphi_sig=1,
		double lambda_theta_bkg=1,
		double lambda_phi_bkg=1,
		double lambda_thetaphi_bkg=1,
		int frameSig=1,//CS...1, HX...2, PX...3
		int frameBkg=1,//CS...1, HX...2, PX...3
		int nGen=1,
		int nState=1,
		Char_t *dirstruct = "ToyDirectory_Default"
		){

	cout<<"/////////////////////////////////"<<endl;
	cout<<"running polGen.C ........////////"<<endl;
	cout<<"/////////////////////////////////"<<endl;

	char frameSigName[200];
	if(frameSig==1)sprintf(frameSigName,"CS");
	if(frameSig==2)sprintf(frameSigName,"HX");
	if(frameSig==3)sprintf(frameSigName,"PX");
	char frameBkgName[200];
	if(frameBkg==1)sprintf(frameBkgName,"CS");
	if(frameBkg==2)sprintf(frameBkgName,"HX");
	if(frameBkg==3)sprintf(frameBkgName,"PX");


	if(frameSig==2) HX_is_natural_sig=true; if(frameSig==3) PX_is_natural_sig=true; //else CS is the natural frame by default
	if(frameBkg==2) HX_is_natural_bkg=true; if(frameBkg==3) PX_is_natural_bkg=true; //else CS is the natural frame by default

	cout<<"Number of Events to be generated ........... "<<n_events<<endl;
	cout<<"pT min ..................................... "<<pTdilepton_min<<endl;
	cout<<"pT max ..................................... "<<pTdilepton_max<<endl;
	cout<<"Rapidity min ............................... "<<rapdilepton_min<<endl;
	cout<<"Rapidity max ............................... "<<rapdilepton_max<<endl;
	cout<<"Background Fraction ........................ "<<f_BG<<endl;
	cout<<"Injected lambda_theta Signal ............... "<<lambda_theta_sig<<" , in the "<<frameSigName<<" frame"<<endl;
	cout<<"Injected lambda_phi Signal ................. "<<lambda_phi_sig<<" , in the "<<frameSigName<<" frame"<<endl;
	cout<<"Injected lambda_thetaphi Signal ............ "<<lambda_thetaphi_sig<<" , in the "<<frameSigName<<" frame"<<endl;
	cout<<"Injected lambda_theta Background............ "<<lambda_theta_bkg<<" , in the "<<frameBkgName<<" frame"<<endl;
	cout<<"Injected lambda_phi Background ............. "<<lambda_phi_bkg<<" , in the "<<frameBkgName<<" frame"<<endl;
	cout<<"Injected lambda_thetaphi Background ........ "<<lambda_thetaphi_bkg<<" , in the "<<frameBkgName<<" frame"<<endl;
	cout<<"Number of Generation ....................... "<<nGen<<endl;
	cout<<"nState ..................................... "<<nState<<endl;
	cout<<"Directory Structure of Output .............. "<<dirstruct<<endl;


	double mass_min = mass_signal_peak - n_sigmas_signal*mass_signal_sigma;
	double mass_max = mass_signal_peak + n_sigmas_signal*mass_signal_sigma;

	char outfilename [500];
	sprintf(outfilename,"%s/genData.root",dirstruct);

	gROOT->Reset();

	delete gRandom;
	gRandom = new TRandom3(0);

	TF1* pT_distr = new TF1("pT_distr",func_pT_gen,pTdilepton_min,pTdilepton_max,1);
	pT_distr->SetParameter(0,nState);
	TF1* rap_distr = new TF1("rap_distr",func_rap_gen,rapdilepton_min,rapdilepton_max,0);

	TFile* hfileout = new TFile(outfilename, "RECREATE", "genData");
	TTree* genData = new TTree("genData","genData");

	// Structure of output ntuple

	TLorentzVector* lepP = new TLorentzVector(0.,0.,0.,0.);
	genData->Branch("lepP","TLorentzVector",&lepP);
	TLorentzVector* lepN = new TLorentzVector(0.,0.,0.,0.);
	genData->Branch("lepN","TLorentzVector",&lepN);

	double costh_CS;  genData->Branch("costh_CS",     &costh_CS,     "costh_CS/D");
	double phi_CS;    genData->Branch("phi_CS",       &phi_CS,       "phi_CS/D"  );
	double phith_CS;  genData->Branch("phith_CS",     &phith_CS,     "phith_CS/D");

	double costh_HX;  genData->Branch("costh_HX",     &costh_HX,     "costh_HX/D");
	double phi_HX;    genData->Branch("phi_HX",       &phi_HX,       "phi_HX/D"  );
	double phith_HX;  genData->Branch("phith_HX",     &phith_HX,     "phith_HX/D");

	double costh_PX;  genData->Branch("costh_PX",     &costh_PX,     "costh_PX/D");
	double phi_PX;    genData->Branch("phi_PX",       &phi_PX,       "phi_PX/D"  );
	double phith_PX;  genData->Branch("phith_PX",     &phith_PX,     "phith_PX/D");

	double cosalpha;  genData->Branch("cosalpha",     &cosalpha,     "cosalpha/D");

	double pT;        genData->Branch("pT",           &pT,           "pT/D"  );
	double rap;       genData->Branch("rap",          &rap,          "rap/D" );
	double mass;      genData->Branch("mass",         &mass,         "mass/D");

	double deltaHXCS; genData->Branch("deltaHXCS",    &deltaHXCS,    "deltaHXCS/D");

	int isBG;         genData->Branch("isBG",         &isBG,         "isBG/I");


	// extremes and binning of lambda_gen extraction histos
	const double l_min = -1;
	const double l_max =  1;
	const double l_step_1D = 0.02;


	TH1D* h_costh2_CS = new TH1D( "h_costh2_CS", "", int((l_max-l_min)/l_step_1D), l_min, l_max );
	TH1D* h_cos2ph_CS = new TH1D( "h_cos2ph_CS", "", int((l_max-l_min)/l_step_1D), l_min, l_max );
	TH1D* h_sin2thcosph_CS = new TH1D( "h_sin2thcosph_CS", "", int((l_max-l_min)/l_step_1D), l_min, l_max );

	TH1D* h_costh2_HX = new TH1D( "h_costh2_HX", "", int((l_max-l_min)/l_step_1D), l_min, l_max );
	TH1D* h_cos2ph_HX = new TH1D( "h_cos2ph_HX", "", int((l_max-l_min)/l_step_1D), l_min, l_max );
	TH1D* h_sin2thcosph_HX = new TH1D( "h_sin2thcosph_HX", "", int((l_max-l_min)/l_step_1D), l_min, l_max );

	TH1D* h_costh2_PX = new TH1D( "h_costh2_PX", "", int((l_max-l_min)/l_step_1D), l_min, l_max );
	TH1D* h_cos2ph_PX = new TH1D( "h_cos2ph_PX", "", int((l_max-l_min)/l_step_1D), l_min, l_max );
	TH1D* h_sin2thcosph_PX = new TH1D( "h_sin2thcosph_PX", "", int((l_max-l_min)/l_step_1D), l_min, l_max );

	double costh2;
	double cos2ph;
	double sin2thcosph;
	double Phi;

	const int n_step = n_events/5;
	int n_step_=1;


	cout << endl;
	cout << "Generating " << n_events << " dilepton events"<< endl;
	cout << "------------------------------------------------------------" << endl;
	cout << "Progress: "<<endl;



	/////////////////// CYCLE OF EVENTS ////////////////////////
	for(int i_event = 1; i_event <= n_events; i_event++){

		if (i_event%n_step == 0) {cout << n_step_*20 <<" % "<<endl; n_step_++;}

		// generation of dilepton in the pp event in the pp CM

		// mass

		isBG = 0;

		if ( gRandom->Uniform() < f_BG ) { mass = gRandom->Uniform(mass_min, mass_max); isBG = 1; }
		else { do { mass = gRandom->Gaus(mass_signal_peak, mass_signal_sigma); }
			while ( mass < mass_min || mass > mass_max ); }

			// pT:

			pT = pT_distr->GetRandom();

			// pL:

			double rap_sign = gRandom->Uniform(-1., 1.); rap_sign /= TMath::Abs(rap_sign);
			rap = rap_distr->GetRandom() * rap_sign;
			double mT = sqrt( mass*mass + pT*pT );
			double pL1 = 0.5 *mT * exp(rap);
			double pL2 = - 0.5 *mT * exp(-rap);
			double pL = pL1 + pL2;

			// Phi:

			double Phi = 2. * gPI * gRandom->Uniform(1.);

			// 4-vector:

			TLorentzVector dilepton;
			dilepton.SetXYZM( pT * cos(Phi) , pT * sin(Phi), pL, mass );


			// generation of polarization (generic reference frame)

			double lambda_theta    = lambda_theta_sig;
			double lambda_phi      = lambda_phi_sig;
			double lambda_thetaphi = lambda_thetaphi_sig;
			bool HX_is_natural = HX_is_natural_sig;
			bool PX_is_natural = PX_is_natural_sig;

			if ( isBG ) {
				lambda_theta    = lambda_theta_bkg;
				lambda_phi      = lambda_phi_bkg;
				lambda_thetaphi = lambda_thetaphi_bkg;
				HX_is_natural = HX_is_natural_bkg;
				PX_is_natural = PX_is_natural_bkg;
			}


			double costhphidistr_max = 1. + TMath::Abs(lambda_phi) + TMath::Abs(lambda_thetaphi);
			double costhphidistr_rnd;
			double costhphidistr;
			double costh_gen;
			double sinth_gen;
			double phi_gen;

			if ( lambda_theta > 0. ) costhphidistr_max += lambda_theta;

			do { costh_gen = -1. + 2. * gRandom->Uniform(1.);
				phi_gen   = 2. * gPI * gRandom->Uniform(1.);
				sinth_gen = sqrt( 1. - costh_gen*costh_gen );
				costhphidistr_rnd = costhphidistr_max * gRandom->Uniform(1.);
				costhphidistr = 1. + lambda_theta    * costh_gen*costh_gen
					+ lambda_phi      * sinth_gen*sinth_gen * cos(2.*phi_gen)
					+ lambda_thetaphi * 2.* sinth_gen*costh_gen * cos(phi_gen);
			} while ( costhphidistr_rnd > costhphidistr );


			// lepton momentum in the dilepton rest frame:

			double p_lepton_DILEP = sqrt( 0.25*mass*mass - Mlepton*Mlepton );

			TLorentzVector lepton_DILEP;

			lepton_DILEP.SetXYZM( p_lepton_DILEP * sinth_gen * cos(phi_gen),
					p_lepton_DILEP * sinth_gen * sin(phi_gen),
					p_lepton_DILEP * costh_gen,
					Mlepton );


			// reference directions to calculate angles:

			TVector3 lab_to_dilep = -dilepton.BoostVector();

			TLorentzVector beam1_DILEP = beam1_LAB;
			beam1_DILEP.Boost(lab_to_dilep);         // beam1 in the dilepton rest frame
			TLorentzVector beam2_DILEP = beam2_LAB;
			beam2_DILEP.Boost(lab_to_dilep);         // beam2 in the dilepton rest frame

			TVector3 beam1_direction     = beam1_DILEP.Vect().Unit();
			TVector3 beam2_direction     = beam2_DILEP.Vect().Unit();
			TVector3 dilep_direction     = dilepton.Vect().Unit();
			TVector3 beam1_beam2_bisect  = ( beam1_direction - beam2_direction ).Unit();


			deltaHXCS = dilep_direction.Angle(beam1_beam2_bisect) * 180./gPI;

			// all polarization frames have the same Y axis = the normal to the plane formed by
			// the directions of the colliding hadrons

			TVector3 Yaxis = ( beam1_direction.Cross( beam2_direction ) ).Unit();

			// flip of y axis with rapidity

			if ( rap < 0 ) Yaxis = - Yaxis;

			TVector3 perpendicular_to_beam = ( beam1_beam2_bisect.Cross( Yaxis ) ).Unit();


			// step 1: transform (rotation) lepton momentum components from generation frame
			// to the frame with x,y,z axes as in the laboratory

			TVector3 oldZaxis = beam1_beam2_bisect;
			if ( HX_is_natural ) oldZaxis = dilep_direction;
			if ( PX_is_natural ) oldZaxis = perpendicular_to_beam;

			TVector3 oldYaxis = Yaxis;
			TVector3 oldXaxis = oldYaxis.Cross(oldZaxis);

			TRotation rotation;
			rotation.RotateAxes(oldXaxis, oldYaxis, oldZaxis);
			// transforms coordinates from the "old" frame to the "xyz" frame

			TLorentzVector lepton_DILEP_xyz = lepton_DILEP;

			lepton_DILEP_xyz.Transform(rotation);
			// lepton_DILEP_xyz is the lepton in the dilepton rest frame
			// wrt to the lab axes

			// lepton 4-vectors in the LAB frame:

			TVector3 dilep_to_lab = dilepton.BoostVector();

			*lepP = lepton_DILEP_xyz;
			lepP->Boost(dilep_to_lab);
			lepN->SetPxPyPzE(-lepton_DILEP_xyz.Px(),-lepton_DILEP_xyz.Py(),-lepton_DILEP_xyz.Pz(),lepton_DILEP_xyz.E());
			lepN->Boost(dilep_to_lab);


			/////////////////////////////////////////////////////////////////////
			// CS frame


			TVector3 newZaxis = beam1_beam2_bisect;
			TVector3 newYaxis = Yaxis;
			TVector3 newXaxis = newYaxis.Cross( newZaxis );

			rotation.SetToIdentity();
			rotation.RotateAxes( newXaxis, newYaxis, newZaxis );
			rotation.Invert();  // transforms coordinates from the "xyz" frame to the new frame

			TVector3 lepton_DILEP_rotated = lepton_DILEP_xyz.Vect();

			lepton_DILEP_rotated.Transform(rotation);

			costh_CS = lepton_DILEP_rotated.CosTheta();

			phi_CS = lepton_DILEP_rotated.Phi() * 180. / gPI;

			if ( costh_CS < 0. ) phith_CS = phi_CS - 135.;
			if ( costh_CS > 0. ) phith_CS = phi_CS - 45.;

			if ( phith_CS < -180. ) phith_CS = 360. + phith_CS;


			/////////////////////////////////////////////////////////////////////
			// HELICITY frame

			newZaxis = dilep_direction;
			newYaxis = Yaxis;
			newXaxis = newYaxis.Cross( newZaxis );

			rotation.SetToIdentity();
			rotation.RotateAxes( newXaxis, newYaxis, newZaxis );
			rotation.Invert();

			lepton_DILEP_rotated = lepton_DILEP_xyz.Vect();

			lepton_DILEP_rotated.Transform(rotation);

			costh_HX = lepton_DILEP_rotated.CosTheta();

			phi_HX = lepton_DILEP_rotated.Phi() * 180. / gPI;

			if ( costh_HX < 0. ) phith_HX = phi_HX - 135.;
			if ( costh_HX > 0. ) phith_HX = phi_HX - 45.;

			if ( phith_HX < -180. ) phith_HX = 360. + phith_HX;


			/////////////////////////////////////////////////////////////////////
			// PERPENDICULAR HELICITY frame

			newZaxis = perpendicular_to_beam;
			newYaxis = Yaxis;
			newXaxis = newYaxis.Cross( newZaxis );

			rotation.SetToIdentity();
			rotation.RotateAxes( newXaxis, newYaxis, newZaxis );
			rotation.Invert();

			lepton_DILEP_rotated = lepton_DILEP_xyz.Vect();

			lepton_DILEP_rotated.Transform(rotation);

			costh_PX = lepton_DILEP_rotated.CosTheta();

			phi_PX = lepton_DILEP_rotated.Phi() * 180. / gPI;

			if ( costh_PX < 0. ) phith_PX = phi_PX - 135.;
			if ( costh_PX > 0. ) phith_PX = phi_PX - 45.;

			if ( phith_PX < -180. ) phith_PX = 360. + phith_PX;


			/////////////////////////////////////////////////////////////////////
			// invariant polarization angle

			cosalpha = sqrt( 1. - pow(costh_PX, 2.) ) * sin( lepton_DILEP_rotated.Phi() );



			////// Filling Histograms of costh2, cos2ph and sin2thcosph for the extraction of the actual generated polarization

			if ( !isBG ){
				costh2=pow(costh_CS,2.);
				Phi = phi_CS/180. * gPI ;
				cos2ph = cos(2.*Phi);
				sin2thcosph= sin(2.*acos(costh_CS))*cos(Phi);
				h_costh2_CS->Fill( costh2 );
				h_cos2ph_CS->Fill( cos2ph );
				h_sin2thcosph_CS->Fill( sin2thcosph );

				costh2=pow(costh_HX,2.);
				Phi = phi_HX/180. * gPI ;
				cos2ph = cos(2.*Phi);
				sin2thcosph= sin(2.*acos(costh_HX))*cos(Phi);
				h_costh2_HX->Fill( costh2 );
				h_cos2ph_HX->Fill( cos2ph );
				h_sin2thcosph_HX->Fill( sin2thcosph );

				costh2=pow(costh_PX,2.);
				Phi = phi_PX/180. * gPI ;
				cos2ph = cos(2.*Phi);
				sin2thcosph= sin(2.*acos(costh_PX))*cos(Phi);
				h_costh2_PX->Fill( costh2 );
				h_cos2ph_PX->Fill( cos2ph );
				h_sin2thcosph_PX->Fill( sin2thcosph );
			}


			//  filling of the ntuple:

			genData->Fill();


	} // end of external loop (generated events)

	cout << endl;


	double lamth_CS;
	double lamph_CS;
	double lamtp_CS;
	costh2=h_costh2_CS->GetMean();
	lamth_CS = (1. - 3. * costh2 ) / ( costh2 - 3./5. );
	cos2ph=h_cos2ph_CS->GetMean();
	lamph_CS = cos2ph * (3. + lamth_CS);
	sin2thcosph=h_sin2thcosph_CS->GetMean();
	lamtp_CS = sin2thcosph * 5./4. * (3. + lamth_CS);

	double lamth_HX;
	double lamph_HX;
	double lamtp_HX;
	costh2=h_costh2_HX->GetMean();
	lamth_HX = (1. - 3. * costh2 ) / ( costh2 - 3./5. );
	cos2ph=h_cos2ph_HX->GetMean();
	lamph_HX = cos2ph * (3. + lamth_HX);
	sin2thcosph=h_sin2thcosph_HX->GetMean();
	lamtp_HX = sin2thcosph * 5./4. * (3. + lamth_HX);

	double lamth_PX;
	double lamph_PX;
	double lamtp_PX;
	costh2=h_costh2_PX->GetMean();
	lamth_PX = (1. - 3. * costh2 ) / ( costh2 - 3./5. );
	cos2ph=h_cos2ph_PX->GetMean();
	lamph_PX = cos2ph * (3. + lamth_PX);
	sin2thcosph=h_sin2thcosph_PX->GetMean();
	lamtp_PX = sin2thcosph * 5./4. * (3. + lamth_PX);


	char resfilename[200];
	sprintf(resfilename,"%s/GenResults.root",dirstruct);
	TFile* GenResultFile = new TFile(resfilename, "RECREATE", "GenResultFile");

	TTree* GenResults = new TTree("GenResults","GenResults");
	GenResults->Branch("lthCS",         &lamth_CS,         "lthCS/D");
	GenResults->Branch("lphCS",         &lamph_CS,         "lphCS/D");
	GenResults->Branch("ltpCS",         &lamtp_CS,         "ltpCS/D");
	GenResults->Branch("lthHX",         &lamth_HX,         "lthHX/D");
	GenResults->Branch("lphHX",         &lamph_HX,         "lphHX/D");
	GenResults->Branch("ltpHX",         &lamtp_HX,         "ltpHX/D");
	GenResults->Branch("lthPX",         &lamth_PX,         "lthPX/D");
	GenResults->Branch("lphPX",         &lamph_PX,         "lphPX/D");
	GenResults->Branch("ltpPX",         &lamtp_PX,         "ltpPX/D");

	GenResults->Fill();


	//GenResultFile->cd(); //Linlin
	//GenResults->Write(); //Linlin
	GenResultFile->Write();
	GenResultFile->Close();

	//hfileout->cd(); //Linlin
	//genData->Write(); //Linlin
	hfileout->Write();
	hfileout->Close();

} // end of main
