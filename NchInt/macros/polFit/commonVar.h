#include "TLorentzVector.h"
#include "TMath.h"

namespace onia{

  // beam energy in GeV
  const double pbeam = 3500.;
  // masses
  const double Mprot = 0.9382720;
  const double muMass = 0.105658;
	const double MpsiPDG = 3.096916;
  const double Ebeam = sqrt( pbeam*pbeam + Mprot*Mprot );
  const TLorentzVector beam1_LAB( 0., 0., pbeam, Ebeam );
  const TLorentzVector beam2_LAB( 0., 0., -pbeam, Ebeam );
	
  // mass ranges 
  const double massMin = 2.85;
  const double massMax = 3.3;
  const double nSigMass = 3.0;
  const double nSigBkgLow = 4.0;
  const double nSigBkgHigh = 3.5;
	
  // Binning
  const int kNbRapForPTBins = 2;
  double rapForPTRange[kNbRapForPTBins+1] = {0., 0.6, 1.2};
  //study the negative and positive rapidity sides separately
  double rapRange[2*kNbRapForPTBins+1] = {-1.2, -0.6, -0., 0.6, 1.2};

  //Jpsi
  const int kNbPTMaxBins = 12;
  const int kNbPTBins[kNbRapForPTBins+1] = {kNbPTMaxBins, kNbPTMaxBins, kNbPTMaxBins};//all y, y1 
  double pTRange[kNbRapForPTBins+1][kNbPTMaxBins+1] = {
    {10., 12., 14., 16., 18., 20., 22., 25., 30., 35., 40., 50., 70.},//all rapidities
    {10., 12., 14., 16., 18., 20., 22., 25., 30., 35., 40., 50., 70.},//mid-rapidity
    {10., 12., 14., 16., 18., 20., 22., 25., 30., 35., 40., 50., 70.}};//forward rapidity

  //number of reference frames
  const int kNbFrames = 6;
  const char *frameLabel[kNbFrames] = {"CS", "HX", "PHX", "sGJ", "GJ1", "GJ2"};
  enum {CS, HX, PHX, sGJ, GJ1, GJ2};
  
  //polarization variables
  const int kNbPolVar = 2; //cosTheta, phi
  enum {cosThPol,phiPol};
  //cosTheta
  const int kNbBinsCosT = 16;
  double cosTMin = -1., cosTMax = 1.;
  //phi for pol. 
  const int kNbBinsPhiPol = 16;
  double phiPolMin = -180., phiPolMax = 180.;


  //phase space limiting cuts:
  double rap = 1.2;

  //some make up to use the same colour and marker for each pT and rapidity bin
  //in every plotting macro:
  int colour_pT[] = {1, 2, 3, 4, 6, 7, 8, 49, 38, 46, 12, 40, 50};
  int marker_pT[] = {20, 21, 25, 22, 23, 26, 27, 28, 29, 30, 20, 20, 50};
  int colour_rapForPTBins[] = {1, 30, 4, 2, 3, kMagenta+1};
  int marker_rapForPTBins[] = {20, 21, 25, 20, 22, 29};

}
