//#include "commonVar.h"
//#include "rootIncludes.inc"

Double_t thisCosTh[onia::kNbFrames], thisPhi[onia::kNbFrames];
Double_t thisPhi_rad[onia::kNbFrames], thisCosPhi[onia::kNbFrames];

//=========================================
// calculation of decay angular parameters
//=========================================
void calcPol(TLorentzVector muplus_LAB,
	     TLorentzVector muminus_LAB){
  
  TLorentzVector qqbar_LAB = muplus_LAB + muminus_LAB;
  Double_t rapidity = qqbar_LAB.Rapidity();

  // boost beams and positive muon into the q-qbar rest frame:
  TVector3 LAB_to_QQBAR = -qqbar_LAB.BoostVector();

  TLorentzVector beam1_QQBAR = onia::beam1_LAB;
  beam1_QQBAR.Boost( LAB_to_QQBAR );

  TLorentzVector beam2_QQBAR = onia::beam2_LAB;
  beam2_QQBAR.Boost( LAB_to_QQBAR );

  TLorentzVector muplus_QQBAR = muplus_LAB;
  muplus_QQBAR.Boost( LAB_to_QQBAR );

  // reference directions in the Jpsi rest frame:

  TVector3 beam1_direction     = beam1_QQBAR.Vect().Unit();
  TVector3 beam2_direction     = beam2_QQBAR.Vect().Unit();
  TVector3 qqbar_direction     = qqbar_LAB.Vect().Unit();
  TVector3 beam1_beam2_bisect  = ( beam1_direction - beam2_direction ).Unit();

  // all polarization frames have the same Y axis = the normal to the plane formed by
  // the directions of the colliding hadrons
  TVector3 Yaxis = ( beam1_direction.Cross( beam2_direction ) ).Unit();
  if ( rapidity < 0. ) Yaxis = -Yaxis; //H: added (5 Dec 2010)

  /////////////////////////////////////////////////////////////////////
  // CS frame

  TVector3 newZaxis = beam1_beam2_bisect;
  TVector3 newYaxis = Yaxis;
  TVector3 newXaxis = newYaxis.Cross( newZaxis );

  TRotation rotation;
  rotation.RotateAxes( newXaxis, newYaxis, newZaxis );
  rotation.Invert();   // transforms coordinates from the "xyz" system
  // to the "new" (rotated) system having the polarization axis
  // as z axis

  TVector3 muplus_QQBAR_rotated(muplus_QQBAR.Vect());
  
  muplus_QQBAR_rotated.Transform( rotation );

  thisCosTh[onia::CS] = muplus_QQBAR_rotated.CosTheta();

  thisPhi_rad[onia::CS] = muplus_QQBAR_rotated.Phi();
  thisPhi[onia::CS] = muplus_QQBAR_rotated.Phi() * 180. / TMath::Pi();
  //if ( thisPhi[onia::CS] < 0. ) thisPhi[onia::CS]= 360. + thisPhi[onia::CS];      // phi defined in degrees from 0 to 360
  //  thisPhi[onia::CS] += 180.; //H: don't add anything...

  /////////////////////////////////////////////////////////////////////
  // HELICITY frame

  newZaxis = qqbar_direction;
  newYaxis = Yaxis;
  newXaxis = newYaxis.Cross( newZaxis );

  rotation.SetToIdentity();
  rotation.RotateAxes( newXaxis, newYaxis, newZaxis );
  rotation.Invert();

  muplus_QQBAR_rotated = muplus_QQBAR.Vect();

  muplus_QQBAR_rotated.Transform( rotation );

  thisCosTh[onia::HX] = muplus_QQBAR_rotated.CosTheta();

  thisPhi_rad[onia::HX] = muplus_QQBAR_rotated.Phi();
  thisPhi[onia::HX] = muplus_QQBAR_rotated.Phi() * 180. / TMath::Pi();
  //if ( thisPhi[onia::HX] < 0. ) thisPhi[onia::HX] = 360. + thisPhi[onia::HX]; // phi defined in degrees from 0 to 360
  //thisPhi[onia::HX] += 180.;//H: don't add anything...
 
  /////////////////////////////////////////////////////////////////////
  // GJ1 frame

  newZaxis = beam1_direction;
  newYaxis = Yaxis;
  newXaxis = newYaxis.Cross( newZaxis );

  rotation.SetToIdentity();
  rotation.RotateAxes( newXaxis, newYaxis, newZaxis );
  rotation.Invert();

  muplus_QQBAR_rotated = muplus_QQBAR.Vect();

  muplus_QQBAR_rotated.Transform( rotation );

  thisCosTh[onia::GJ1] = muplus_QQBAR_rotated.CosTheta();

  thisPhi_rad[onia::GJ1] = muplus_QQBAR_rotated.Phi();
  thisPhi[onia::GJ1] = muplus_QQBAR_rotated.Phi() * 180. / TMath::Pi();
  //if ( thisPhi[onia::GJ1] < 0. ) thisPhi[onia::GJ1] = 360. + thisPhi[onia::GJ1]; // phi defined in degrees from 0 to 360
  //thisPhi[onia::GJ1] += 180.;//H: don't add anything...

  /////////////////////////////////////////////////////////////////////
  // GJ2 frame

  newZaxis = beam2_direction;
  newYaxis = Yaxis;
  newXaxis = newYaxis.Cross( newZaxis );

  rotation.SetToIdentity();
  rotation.RotateAxes( newXaxis, newYaxis, newZaxis );
  rotation.Invert();

  muplus_QQBAR_rotated = muplus_QQBAR.Vect();

  muplus_QQBAR_rotated.Transform( rotation );

  thisCosTh[onia::GJ2] = muplus_QQBAR_rotated.CosTheta();

  thisPhi_rad[onia::GJ2] = muplus_QQBAR_rotated.Phi();
  thisPhi[onia::GJ2] = muplus_QQBAR_rotated.Phi() * 180. / TMath::Pi();
  //if ( thisPhi[onia::GJ2] < 0. ) thisPhi[onia::GJ2] = 360. + thisPhi[onia::GJ2]; // phi defined in degrees from 0 to 360
  //thisPhi[onia::GJ2] += 180.;//H: don't add anything...
  /////////////////////////////////////////////////////////////////////
  // sGJ frame (symmetrized GJ)

  newZaxis = beam1_direction; if( rapidity < 0. ) newZaxis = beam2_direction;
  newYaxis = Yaxis;

  // try to swith the following line on or off
  //if( rapidity < 0. ) newYaxis = -Yaxis;

  newXaxis = newYaxis.Cross( newZaxis );

  rotation.SetToIdentity();
  rotation.RotateAxes( newXaxis, newYaxis, newZaxis );
  rotation.Invert();

  muplus_QQBAR_rotated = muplus_QQBAR.Vect();

  muplus_QQBAR_rotated.Transform( rotation );

  thisCosTh[onia::sGJ] = muplus_QQBAR_rotated.CosTheta();

  thisPhi_rad[onia::sGJ] = muplus_QQBAR_rotated.Phi();
  thisPhi[onia::sGJ] = muplus_QQBAR_rotated.Phi() * 180. / TMath::Pi();
  //if ( thisPhi[onia::sGJ] < 0. ) thisPhi[onia::sGJ] = 360. + thisPhi[onia::sGJ]; // phi defined in degrees from 0 to 360
  //thisPhi[onia::sGJ] += 180.;//H: don't add anything...

  /////////////////////////////////////////////////////////////////////
  // PHX frame ("perpendicular helicity frame" - z axis perpendicular
  // to the CS axis)

  newZaxis = newZaxis = ( beam1_beam2_bisect.Cross( Yaxis ) ).Unit();
  newYaxis = Yaxis;
  newXaxis = newYaxis.Cross( newZaxis );

  rotation.SetToIdentity();
  rotation.RotateAxes( newXaxis, newYaxis, newZaxis );
  rotation.Invert();

  muplus_QQBAR_rotated = muplus_QQBAR.Vect();

  muplus_QQBAR_rotated.Transform( rotation );

  thisCosTh[onia::PHX] = muplus_QQBAR_rotated.CosTheta();

  thisPhi_rad[onia::PHX] = muplus_QQBAR_rotated.Phi();
  thisPhi[onia::PHX] = muplus_QQBAR_rotated.Phi() * 180. / TMath::Pi();
  //if ( thisPhi[onia::PHX] < 0. ) thisPhi[onia::PHX] = 360. + thisPhi[onia::PHX]; // phi defined in degrees from 0 to 360
  //thisPhi[onia::PHX] += 180.;//H: don't add anything...
}
