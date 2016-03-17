/*
 * MattsResults.C
 *
 *  Created on: Dec 4, 2011
 *      Author: valentinknuenz
 */


#include "Riostream.h"
#include "TSystem.h"
#include "TStyle.h"
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
#include "TH3.h"
#include "TGraphAsymmErrors.h"
#include "TGraph2D.h"
#include "TGraph.h"


int main(){

	char filename[200], GraphName[200];
	char tmpfilename[200];
	sprintf(tmpfilename,"TGraphResults_1SUps_MattResults.root");		gSystem->Unlink(tmpfilename);
	sprintf(tmpfilename,"TGraphResults_1SUps_MattSyst.root");			gSystem->Unlink(tmpfilename);


	for(int nState=1;nState<4;nState++){


	const int nBinsInPt=8;

	double ptCentre_err[nBinsInPt]={1.,1.,1.,1.,2.,2.5,3.,8.5};

	/// arrays hier platzieren, output file schauen, dass es in zwei geht...
/*
///// MIT results
	double ptCentre_[nBinsInPt]={2.5,7.5,12.5,17.5,22.5,27.5};
	double ptCentre_err[nBinsInPt]={2.5,2.5,2.5,2.5,2.5,2.5};
	double lmean[nBinsInPt];
	double StatError[nBinsInPt];
	double SystError[nBinsInPt];
	double ZeroError[nBinsInPt]={0,0,0,0,0,0};

	double lth_CS[2][nBinsInPt]={{0.5,0.2,0.09,0.2,0.04,-0.1},{-0.1,-0.2,-0.1,-0.09,-0.01,-0.04}};
	double lph_CS[2][nBinsInPt]={{-0.09,-0.08,-0.03,-0.06,-0.02,0.3},{-0.08,-0.009,-0.02,-0.02,-0.07,-0.03}};
	double ltp_CS[2][nBinsInPt]={{-0.05,-0.1,-0.01,-0.02,0.02,0.02},{-0.07,-0.01,-0.008,0.008,0.05,-0.01}};
	double lti_CS[2][nBinsInPt]={{0.2,-0.05,0.02,-0.02,-0.02,1},{-0.4,-0.2,-0.06,-0.2,-0.2,-0.1}};

	double lth_CS_stat[2][nBinsInPt]={{0.05,0.02,0.02,0.04,0.05,0.08},{0.04,0.01,0.02,0.03,0.05,0.07}};
	double lph_CS_stat[2][nBinsInPt]={{0.009,0.02,0.02,0.04,0.05,0.06},{0.009,0.02,0.02,0.03,0.05,0.07}};
	double ltp_CS_stat[2][nBinsInPt]={{0.02,0.02,0.01,0.02,0.03,0.04},{0.02,0.02,0.01,0.02,0.03,0.04}};
	double lti_CS_stat[2][nBinsInPt]={{0.05,0.04,0.05,0.08,0.1,0.3},{0.05,0.04,0.05,0.08,0.1,0.2}};

	double lth_CS_syst[2][nBinsInPt]={{0.2,0.2,0.05,0.09,0.2,0.2},{0.6,0.2,0.07,0.03,0.08,0.2}};
	double lph_CS_syst[2][nBinsInPt]={{0.03,0.1,0.06,0.03,0.05,0.02},{0.04,0.08,0.05,0.05,0.06,0.05}};
	double ltp_CS_syst[2][nBinsInPt]={{0.1,0.1,0.07,0.02,0.01,0.01},{0.07,0.08,0.04,0.02,0.01,0.06}};
	double lti_CS_syst[2][nBinsInPt]={{0.2,0.3,0.2,0.06,0.05,0.3},{0.7,0.2,0.1,0.1,0.09,0.06}};


	double lth_HX[2][nBinsInPt]={{-0.3,-0.3,-0.1,-0.1,-0.0008,0.3},{0.03,-0.07,0.05,0.02,-0.09,0.06}};
	double lph_HX[2][nBinsInPt]={{0.01,0.09,0.03,0.05,0.01,0.1},{-0.05,-0.04,-0.04,-0.05,-0.06,-0.005}};
	double ltp_HX[2][nBinsInPt]={{-0.09,0.1,-0.004,-0.01,-0.02,-0.03},{-0.1,0.07,0.07,0.02,-0.06,0.03}};
	double lti_HX[2][nBinsInPt]={{-0.3,-0.06,-0.04,-0.01,0.03,0.7},{-0.1,-0.2,-0.07,-0.1,-0.2,0.05}};

	double lth_HX_stat[2][nBinsInPt]={{0.02,0.02,0.04,0.05,0.08,0.1},{0.05,0.03,0.04,0.05,0.08,0.1}};
	double lph_HX_stat[2][nBinsInPt]={{0.01,0.008,0.007,0.01,0.02,0.03},{0.009,0.01,0.01,0.01,0.02,0.04}};
	double ltp_HX_stat[2][nBinsInPt]={{0.02,0.01,0.01,0.02,0.03,0.04},{0.02,0.01,0.02,0.02,0.03,0.05}};
	double lti_HX_stat[2][nBinsInPt]={{0.04,0.04,0.05,0.07,0.1,0.2},{0.05,0.04,0.05,0.07,0.1,0.2}};

	double lth_HX_syst[2][nBinsInPt]={{0.2,0.2,0.09,0.06,0.09,0.1},{0.6,0.1,0.05,0.06,0.08,0.1}};
	double lph_HX_syst[2][nBinsInPt]={{0.04,0.06,0.01,0.01,0.03,0.08},{0.04,0.09,0.05,0.03,0.01,0.04}};
	double ltp_HX_syst[2][nBinsInPt]={{0.06,0.06,0.04,0.02,0.02,0.01},{0.2,0.06,0.05,0.04,0.03,0.03}};
	double lti_HX_syst[2][nBinsInPt]={{0.2,0.3,0.1,0.05,0.05,0.2},{0.7,0.3,0.2,0.1,0.07,0.06}};
*/

///// CDF 2012 results
	double ptBorders_[nBinsInPt+1]={0,2,4,6,8,12,17,23,40};
	double ptCentre_[nBinsInPt]={1.15,2.925,4.9,7.,9.95,14.,19.95,26.8};
	double ptCentre_err_low[nBinsInPt];
	double ptCentre_err_high[nBinsInPt];
	double lmean[nBinsInPt];
	double StatError[nBinsInPt];
	double StatError_low[nBinsInPt];
	double StatError_high[nBinsInPt];
	double TotalError_low[nBinsInPt];
	double TotalError_high[nBinsInPt];
	double SystError[nBinsInPt];
	double ZeroError[nBinsInPt]={0,0,0,0,0,0,0,0};

	double XXX=0.;


	double lth_CS[3][2][nBinsInPt]={{{ 0.214,
			-0.002,
			-0.053,
			 0.030,
			 0.012,
			 0.029,
			 0.033,
			 0.079},{0.,0.,0.,0.,0.,0.,0.,0.}},{{-0.591,
					 -0.115,
					 -0.042,
					 -0.194,
					  0.007,
					 -0.041,
					 -0.055,
					  0.151},{0.,0.,0.,0.,0.,0.,0.,0.}},{{ 0.371,
							  -0.659,
							  -0.421,
							  -0.049,
							  -0.286,
							  -0.217,
							  -0.198,
							   0.021},{0.,0.,0.,0.,0.,0.,0.,0.}}};
	double lph_CS[3][2][nBinsInPt]={{{-0.021,
			-0.049,
			-0.069,
			-0.104,
			-0.067,
			-0.015,
			-0.172,
			-0.082},{0.,0.,0.,0.,0.,0.,0.,0.}},{{ 0.200,
					-0.120,
					-0.031,
					-0.038,
					-0.072,
					 0.041,
					 0.048,
					-0.175},{0.,0.,0.,0.,0.,0.,0.,0.}},{{ 0.435,
							-0.053,
							-0.114,
							-0.077,
							-0.072,
							 0.058,
							 0.085,
							-0.023},{0.,0.,0.,0.,0.,0.,0.,0.}}};
	double ltp_CS[3][2][nBinsInPt]={{{ 0.022,
			 0.021,
			 0.012,
			-0.077,
			 0.031,
			 0.094,
			-0.065,
			-0.011},{0.,0.,0.,0.,0.,0.,0.,0.}},{{ 0.143,
					 0.144,
					 0.051,
					 0.071,
					 0.012,
					 0.126,
					-0.018,
					-0.224},{0.,0.,0.,0.,0.,0.,0.,0.}},{{ 0.052,
							 0.049,
							 0.181,
							 0.289,
							 0.065,
							 0.014,
							-0.139,
							-0.216},{0.,0.,0.,0.,0.,0.,0.,0.}}};
	double lti_CS[3][2][nBinsInPt]={{{ 0.149,
			-0.143,
			-0.244,
			-0.255,
			-0.176,
			-0.016,
			-0.412,
			-0.153},{0.,0.,0.,0.,0.,0.,0.,0.}},{{ 0.012,
					-0.423,
					-0.132,
					-0.296,
					-0.196,
					 0.086,
					 0.093,
					-0.317},{0.,0.,0.,0.,0.,0.,0.,0.}},{{ 2.961,
							-0.776,
							-0.688,
							-0.260,
							-0.467,
							-0.047,
							 0.064,
							-0.046},{0.,0.,0.,0.,0.,0.,0.,0.}}};

	double lth_CS_stat_low[3][2][nBinsInPt]={{{0.139,
			0.073,
			0.068,
			0.072,
			0.051,
			0.055,
			0.087,
			0.180},{0.,0.,0.,0.,0.,0.,0.,0.}},{{0.248,
					0.181,
					0.163,
					0.145,
					0.106,
					0.094,
					0.126,
					0.246},{0.,0.,0.,0.,0.,0.,0.,0.}},{{0.458,
							0.288,
							0.284,
							0.253,
							0.135,
							0.111,
							0.145,
							0.257},{0.,0.,0.,0.,0.,0.,0.,0.}}};
	double lph_CS_stat_low[3][2][nBinsInPt]={{{0.039,
			0.015,
			0.016,
			0.029,
			0.031,
			0.046,
			0.088,
			0.128},{0.,0.,0.,0.,0.,0.,0.,0.}},{{0.071,
					0.042,
					0.035,
					0.048,
					0.058,
					0.075,
					0.127,
					0.317},{0.,0.,0.,0.,0.,0.,0.,0.}},{{0.108,
							0.069,
							0.069,
							0.077,
							0.080,
							0.094,
							0.136,
							0.270},{0.,0.,0.,0.,0.,0.,0.,0.}}};
	double ltp_CS_stat_low[3][2][nBinsInPt]={{{0.053,
			0.027,
			0.030,
			0.039,
			0.034,
			0.039,
			0.057,
			0.123},{0.,0.,0.,0.,0.,0.,0.,0.}},{{0.103,
					0.065,
					0.066,
					0.077,
					0.067,
					0.072,
					0.090,
					0.181},{0.,0.,0.,0.,0.,0.,0.,0.}},{{0.177,
							0.109,
							0.117,
							0.129,
							0.092,
							0.091,
							0.088,
							0.199},{0.,0.,0.,0.,0.,0.,0.,0.}}};
	double lti_CS_stat_low[3][2][nBinsInPt]={{{0.164,
			0.078,
			0.073,
			0.089,
			0.077,
			0.107,
			0.160,
			0.310},{0.,0.,0.,0.,0.,0.,0.,0.}},{{0.392,
					0.182,
					0.182,
					0.181,
					0.154,
					0.197,
					0.333,
					0.465},{0.,0.,0.,0.,0.,0.,0.,0.}},{{1.201,
							0.302,
							0.288,
							0.281,
							0.198,
							0.245,
							0.399,
							0.562},{0.,0.,0.,0.,0.,0.,0.,0.}}};

	double lth_CS_stat_high[3][2][nBinsInPt]={{{+0.144,
			+0.074,
			+0.069,
			+0.074,
			+0.053,
			+0.057,
			+0.099,
			+0.212},{0.,0.,0.,0.,0.,0.,0.,0.}},{{+0.262,
					+0.187,
					+0.169,
					+0.152,
					+0.112,
					+0.100,
					+0.156,
					+0.380},{0.,0.,0.,0.,0.,0.,0.,0.}},{{+0.501,
							+0.307,
							+0.297,
							+0.270,
							+0.144,
							+0.123,
							+0.182,
							+0.368},{0.,0.,0.,0.,0.,0.,0.,0.}}};
	double lph_CS_stat_high[3][2][nBinsInPt]={{{+0.036,
			+0.015,
			+0.015,
			+0.027,
			+0.030,
			+0.043,
			+0.076,
			+0.103},{0.,0.,0.,0.,0.,0.,0.,0.}},{{+0.064,
					+0.041,
					+0.034,
					+0.046,
					+0.054,
					+0.069,
					+0.101,
					+0.213},{0.,0.,0.,0.,0.,0.,0.,0.}},{{+0.093,
							+0.065,
							+0.064,
							+0.070,
							+0.073,
							+0.084,
							+0.095,
							+0.194},{0.,0.,0.,0.,0.,0.,0.,0.}}};
	double ltp_CS_stat_high[3][2][nBinsInPt]={{{+0.053,
			+0.027,
			+0.031,
			+0.040,
			+0.035,
			+0.039,
			+0.056,
			+0.120},{0.,0.,0.,0.,0.,0.,0.,0.}},{{+0.104,
					+0.065,
					+0.066,
					+0.078,
					+0.068,
					+0.073,
					+0.090,
					+0.147},{0.,0.,0.,0.,0.,0.,0.,0.}},{{+0.180,
							+0.110,
							+0.119,
							+0.133,
							+0.093,
							+0.092,
							+0.081,
							+0.182},{0.,0.,0.,0.,0.,0.,0.,0.}}};
	double lti_CS_stat_high[3][2][nBinsInPt]={{{+0.169,
			+0.079,
			+0.074,
			+0.090,
			+0.079,
			+0.113,
			+0.172,
			+0.335},{0.,0.,0.,0.,0.,0.,0.,0.}},{{+0.427,
					+0.189,
					+0.190,
					+0.190,
					+0.162,
					+0.214,
					+0.381,
					+0.575},{0.,0.,0.,0.,0.,0.,0.,0.}},{{+1.499,
							+0.326,
							+0.307,
							+0.304,
							+0.212,
							+0.270,
							+0.462,
							+0.788},{0.,0.,0.,0.,0.,0.,0.,0.}}};

	double lth_CS_syst[3][2][nBinsInPt]={{{0.052,
			0.039,
			0.036,
			0.022,
			0.020,
			0.019,
			0.047,
			0.020},{0.,0.,0.,0.,0.,0.,0.,0.}},{{0.070,
					0.028,
					0.062,
					0.051,
					0.029,
					0.028,
					0.068,
					0.213},{0.,0.,0.,0.,0.,0.,0.,0.}},{{0.109,
							0.135,
							0.140,
							0.044,
							0.036,
							0.022,
							0.056,
							0.124},{0.,0.,0.,0.,0.,0.,0.,0.}}};
	double lph_CS_syst[3][2][nBinsInPt]={{{ 0.064,
			 0.022,
			 0.023,
			 0.024,
			 0.018,
			 0.045,
			 0.063,
			 0.085},{0.,0.,0.,0.,0.,0.,0.,0.}},{{0.149,
					 0.021,
					 0.033,
					 0.068,
					 0.018,
					 0.087,
					 0.083,
					 0.274},{0.,0.,0.,0.,0.,0.,0.,0.}},{{0.084,
							 0.234,
							 0.166,
							 0.094,
							 0.022,
							 0.058,
							 0.114,
							 0.164},{0.,0.,0.,0.,0.,0.,0.,0.}}};
	double ltp_CS_syst[3][2][nBinsInPt]={{{0.045,
			0.022,
			0.016,
			0.016,
			0.010,
			0.024,
			0.025,
			0.052},{0.,0.,0.,0.,0.,0.,0.,0.}},{{0.069,
					0.050,
					0.023,
					0.039,
					0.017,
					0.058,
					0.045,
					0.166},{0.,0.,0.,0.,0.,0.,0.,0.}},{{0.038,
							0.122,
							0.099,
							0.086,
							0.025,
							0.029,
							0.088,
							0.125},{0.,0.,0.,0.,0.,0.,0.,0.}}};
	double lti_CS_syst[3][2][nBinsInPt]={{{ 0.054,
			 0.057,
			 0.062,
			 0.063,
			 0.039,
			 0.058,
			 0.072,
			 0.074},{0.,0.,0.,0.,0.,0.,0.,0.}},{{0.086,
					 0.047,
					 0.107,
					 0.120,
					 0.044,
					 0.108,
					 0.116,
					 0.111},{0.,0.,0.,0.,0.,0.,0.,0.}},{{0.113,
							 0.071,
							 0.207,
							 0.178,
							 0.058,
							 0.127,
							 0.149,
							 0.154},{0.,0.,0.,0.,0.,0.,0.,0.}}};


//////////////// HX ////////////////////////////


	double lth_HX[3][2][nBinsInPt]={{{ 0.008,
			-0.078,
			-0.111,
			-0.125,
			-0.154,
			-0.133,
			-0.229,
			-0.213},{0.,0.,0.,0.,0.,0.,0.,0.}},{{-0.069,
					-0.287,
					-0.108,
					-0.085,
					-0.179,
					-0.091,
					 0.057,
					 0.004},{0.,0.,0.,0.,0.,0.,0.,0.}},{{ 0.323,
							 -0.190,
							 -0.007,
							 -0.481,
							  0.067,
							  0.143,
							  0.142,
							  0.138},{0.,0.,0.,0.,0.,0.,0.,0.}}};
	double lph_HX[3][2][nBinsInPt]={{{-0.025,
			-0.044,
			-0.070,
			-0.059,
			-0.022,
			-0.004,
			-0.060,
			 0.014},{0.,0.,0.,0.,0.,0.,0.,0.}},{{ 0.194,
					 -0.063,
					 -0.038,
					 -0.013,
					 -0.017,
					 -0.039,
					 -0.026,
					  0.091},{0.,0.,0.,0.,0.,0.,0.,0.}},{{ 0.515,
							  0.002,
							 -0.107,
							  0.139,
							 -0.190,
							 -0.120,
							 -0.076,
							  0.043},{0.,0.,0.,0.,0.,0.,0.,0.}}};
	double ltp_HX[3][2][nBinsInPt]={{{ 0.054,
			-0.020,
			 0.009,
			-0.020,
			-0.060,
			-0.020,
			 0.054,
			-0.046},{0.,0.,0.,0.,0.,0.,0.,0.}},{{-0.261,
					 0.045,
					-0.059,
					-0.083,
					-0.008,
					-0.013,
					-0.148,
					 0.189},{0.,0.,0.,0.,0.,0.,0.,0.}},{{ 0.434,
							 -0.046,
							  0.077,
							  0.050,
							  0.081,
							 -0.004,
							  0.046,
							 -0.287},{0.,0.,0.,0.,0.,0.,0.,0.}}};
	double lti_HX[3][2][nBinsInPt]={{{-0.064,
			-0.201,
			-0.300,
			-0.284,
			-0.215,
			-0.145,
			-0.385,
			-0.173},{0.,0.,0.,0.,0.,0.,0.,0.}},{{ 0.634,
					-0.449,
					-0.214,
					-0.124,
					-0.227,
					-0.201,
					-0.020,
					 0.306},{0.,0.,0.,0.,0.,0.,0.,0.}},{{ 3.855,
							 -0.184,
							 -0.297,
							 -0.077,
							 -0.421,
							 -0.194,
							 -0.081,
							  0.278},{0.,0.,0.,0.,0.,0.,0.,0.}}};

	double lth_HX_stat_low[3][2][nBinsInPt]={{{0.100,
			0.040,
			0.043,
			0.058,
			0.055,
			0.080,
			0.126,
			0.219},{0.,0.,0.,0.,0.,0.,0.,0.}},{{0.205,
					0.092,
					0.093,
					0.113,
					0.103,
					0.141,
					0.236,
					0.334},{0.,0.,0.,0.,0.,0.,0.,0.}},{{0.342,
							0.168,
							0.183,
							0.129,
							0.172,
							0.215,
							0.333,
							0.420},{0.,0.,0.,0.,0.,0.,0.,0.}}};
	double lph_HX_stat_low[3][2][nBinsInPt]={{{0.040,
			0.020,
			0.028,
			0.036,
			0.025,
			0.025,
			0.033,
			0.078},{0.,0.,0.,0.,0.,0.,0.,0.}},{{0.079,
					0.051,
					0.059,
					0.075,
					0.057,
					0.052,
					0.062,
					0.114},{0.,0.,0.,0.,0.,0.,0.,0.}},{{0.099,
							0.092,
							0.128,
							0.093,
							0.098,
							0.079,
							0.088,
							0.128},{0.,0.,0.,0.,0.,0.,0.,0.}}};
	double ltp_HX_stat_low[3][2][nBinsInPt]={{{0.054,
			0.021,
			0.023,
			0.032,
			0.028,
			0.032,
			0.049,
			0.090},{0.,0.,0.,0.,0.,0.,0.,0.}},{{0.130,
					0.053,
					0.051,
					0.063,
					0.055,
					0.063,
					0.075,
					0.127},{0.,0.,0.,0.,0.,0.,0.,0.}},{{0.197,
							0.097,
							0.100,
							0.086,
							0.094,
							0.090,
							0.119,
							0.168},{0.,0.,0.,0.,0.,0.,0.,0.}}};
	double lti_HX_stat_low[3][2][nBinsInPt]={{{ 0.136,
			 0.053,
			 0.055,
			 0.082,
			 0.074,
			 0.098,
			 0.146,
			 0.302},{0.,0.,0.,0.,0.,0.,0.,0.}},{{0.390,
					 0.120,
					 0.125,
					 0.169,
					 0.151,
					 0.179,
					 0.302,
					 0.521},{0.,0.,0.,0.,0.,0.,0.,0.}},{{1.335,
							 0.236,
							 0.228,
							 0.259,
							 0.193,
							 0.231,
							 0.365,
							 0.616},{0.,0.,0.,0.,0.,0.,0.,0.}}};

	double lth_HX_stat_high[3][2][nBinsInPt]={{{+0.102,
			+0.041,
			+0.044,
			+0.060,
			+0.056,
			+0.082,
			+0.131,
			+0.248},{0.,0.,0.,0.,0.,0.,0.,0.}},{{+0.201,
					+0.096,
					+0.100,
					+0.123,
					+0.110,
					+0.147,
					+0.256,
					+0.393},{0.,0.,0.,0.,0.,0.,0.,0.}},{{+0.381,
							+0.183,
							+0.208,
							+0.139,
							+0.187,
							+0.230,
							+0.365,
							+0.503},{0.,0.,0.,0.,0.,0.,0.,0.}}};
	double lph_HX_stat_high[3][2][nBinsInPt]={{{+0.037,
			+0.020,
			+0.027,
			+0.035,
			+0.025,
			+0.024,
			+0.032,
			+0.071},{0.,0.,0.,0.,0.,0.,0.,0.}},{{+0.068,
					+0.048,
					+0.054,
					+0.067,
					+0.052,
					+0.051,
					+0.060,
					+0.112},{0.,0.,0.,0.,0.,0.,0.,0.}},{{+0.085,
							+0.083,
							+0.110,
							+0.083,
							+0.089,
							+0.075,
							+0.083,
							+0.128},{0.,0.,0.,0.,0.,0.,0.,0.}}};
	double ltp_HX_stat_high[3][2][nBinsInPt]={{{+0.054,
			+0.021,
			+0.023,
			+0.032,
			+0.028,
			+0.032,
			+0.049,
			+0.091},{0.,0.,0.,0.,0.,0.,0.,0.}},{{+0.127,
					+0.054,
					+0.051,
					+0.063,
					+0.055,
					+0.063,
					+0.074,
					+0.131},{0.,0.,0.,0.,0.,0.,0.,0.}},{{+0.205,
							+0.097,
							+0.100,
							+0.086,
							+0.095,
							+0.090,
							+0.120,
							+0.164},{0.,0.,0.,0.,0.,0.,0.,0.}}};
	double lti_HX_stat_high[3][2][nBinsInPt]={{{+0.140,
			+0.054,
			+0.056,
			+0.084,
			+0.077,
			+0.102,
			+0.157,
			+0.363},{0.,0.,0.,0.,0.,0.,0.,0.}},{{+0.424,
					+0.125,
					+0.131,
					+0.178,
					+0.156,
					+0.193,
					+0.344,
					+0.585},{0.,0.,0.,0.,0.,0.,0.,0.}},{{+1.753,
							+0.256,
							+0.245,
							+0.288,
							+0.206,
							+0.250,
							+0.426,
							+0.879},{0.,0.,0.,0.,0.,0.,0.,0.}}};

	double lth_HX_syst[3][2][nBinsInPt]={{{0.049,
			0.030,
			0.029,
			0.039,
			0.024,
			0.036,
			0.076,
			0.111},{0.,0.,0.,0.,0.,0.,0.,0.}},{{0.046,
					0.056,
					0.052,
					0.070,
					0.030,
					0.076,
					0.121,
					0.165},{0.,0.,0.,0.,0.,0.,0.,0.}},{{0.115,
							0.075,
							0.087,
							0.134,
							0.039,
							0.074,
							0.124,
							0.135},{0.,0.,0.,0.,0.,0.,0.,0.}}};
	double lph_HX_syst[3][2][nBinsInPt]={{{0.068,
			0.021,
			0.022,
			0.017,
			0.018,
			0.041,
			0.025,
			0.058},{0.,0.,0.,0.,0.,0.,0.,0.}},{{0.132,
					0.014,
					0.031,
					0.059,
					0.018,
					0.097,
					0.051,
					0.160},{0.,0.,0.,0.,0.,0.,0.,0.}},{{0.070,
							0.189,
							0.143,
							0.053,
							0.034,
							0.059,
							0.064,
							0.116},{0.,0.,0.,0.,0.,0.,0.,0.}}};
	double ltp_HX_syst[3][2][nBinsInPt]={{{0.021,
			0.011,
			0.008,
			0.010,
			0.012,
			0.002,
			0.016,
			0.052},{0.,0.,0.,0.,0.,0.,0.,0.}},{{0.015,
					0.014,
					0.014,
					0.018,
					0.020,
					0.013,
					0.057,
					0.091},{0.,0.,0.,0.,0.,0.,0.,0.}},{{0.041,
							0.019,
							0.014,
							0.007,
							0.021,
							0.018,
							0.039,
							0.093},{0.,0.,0.,0.,0.,0.,0.,0.}}};
	double lti_HX_syst[3][2][nBinsInPt]={{{0.057,
			0.036,
			0.044,
			0.060,
			0.042,
			0.059,
			0.084,
			0.117},{0.,0.,0.,0.,0.,0.,0.,0.}},{{0.118,
					0.048,
					0.076,
					0.101,
					0.056,
					0.109,
					0.084,
					0.093},{0.,0.,0.,0.,0.,0.,0.,0.}},{{0.571,
							0.093,
							0.131,
							0.174,
							0.070,
							0.128,
							0.142,
							0.117},{0.,0.,0.,0.,0.,0.,0.,0.}}};

	for(int iLam = 1; iLam<19; iLam++){

	for(int rapBin = 1; rapBin < 3; rapBin++){


	if(iLam==1)  sprintf(GraphName,"lth_CS_rap%d",rapBin);
	if(iLam==2)  sprintf(GraphName,"lph_CS_rap%d",rapBin);
	if(iLam==3)  sprintf(GraphName,"ltp_CS_rap%d",rapBin);
	if(iLam==4)  sprintf(GraphName,"lthstar_CS_rap%d",rapBin);
	if(iLam==5)  sprintf(GraphName,"lphstar_CS_rap%d",rapBin);
	if(iLam==6)  sprintf(GraphName,"ltilde_CS_rap%d",rapBin);

	if(iLam==7)  sprintf(GraphName,"lth_HX_rap%d",rapBin);
	if(iLam==8)  sprintf(GraphName,"lph_HX_rap%d",rapBin);
	if(iLam==9)  sprintf(GraphName,"ltp_HX_rap%d",rapBin);
	if(iLam==10) sprintf(GraphName,"lthstar_HX_rap%d",rapBin);
	if(iLam==11) sprintf(GraphName,"lphstar_HX_rap%d",rapBin);
	if(iLam==12) sprintf(GraphName,"ltilde_HX_rap%d",rapBin);

	if(iLam==13) sprintf(GraphName,"lth_PX_rap%d",rapBin);
	if(iLam==14) sprintf(GraphName,"lph_PX_rap%d",rapBin);
	if(iLam==15) sprintf(GraphName,"ltp_PX_rap%d",rapBin);
	if(iLam==16) sprintf(GraphName,"lthstar_PX_rap%d",rapBin);
	if(iLam==17) sprintf(GraphName,"lphstar_PX_rap%d",rapBin);
	if(iLam==18) sprintf(GraphName,"ltilde_PX_rap%d",rapBin);

	for(int iPt=0;iPt<nBinsInPt;iPt++){
		if(iLam==1) {lmean[iPt]=lth_CS[nState-1][rapBin-1][iPt]; StatError_low[iPt]=lth_CS_stat_low[nState-1][rapBin-1][iPt]; StatError_high[iPt]=lth_CS_stat_high[nState-1][rapBin-1][iPt]; SystError[iPt]=lth_CS_syst[nState-1][rapBin-1][iPt];}
		if(iLam==2) {lmean[iPt]=lph_CS[nState-1][rapBin-1][iPt]; StatError_low[iPt]=lph_CS_stat_low[nState-1][rapBin-1][iPt]; StatError_high[iPt]=lph_CS_stat_high[nState-1][rapBin-1][iPt]; SystError[iPt]=lph_CS_syst[nState-1][rapBin-1][iPt];}
		if(iLam==3) {lmean[iPt]=ltp_CS[nState-1][rapBin-1][iPt]; StatError_low[iPt]=ltp_CS_stat_low[nState-1][rapBin-1][iPt]; StatError_high[iPt]=ltp_CS_stat_high[nState-1][rapBin-1][iPt]; SystError[iPt]=ltp_CS_syst[nState-1][rapBin-1][iPt];}
		if(iLam==4) {lmean[iPt]=1000; StatError_low[iPt]=0; StatError_high[iPt]=0; SystError[iPt]=0;}
		if(iLam==5) {lmean[iPt]=1000; StatError_low[iPt]=0; StatError_high[iPt]=0; SystError[iPt]=0;}
		if(iLam==6) {lmean[iPt]=lti_CS[nState-1][rapBin-1][iPt]; StatError_low[iPt]=lti_CS_stat_low[nState-1][rapBin-1][iPt]; StatError_high[iPt]=lti_CS_stat_high[nState-1][rapBin-1][iPt]; SystError[iPt]=lti_CS_syst[nState-1][rapBin-1][iPt];/*cout<<GraphName<<endl;cout<<lmean[iPt]<<endl;*/}

		if(iLam==7) {lmean[iPt]=lth_HX[nState-1][rapBin-1][iPt]; StatError_low[iPt]=lth_HX_stat_low[nState-1][rapBin-1][iPt]; StatError_high[iPt]=lth_HX_stat_high[nState-1][rapBin-1][iPt]; SystError[iPt]=lth_HX_syst[nState-1][rapBin-1][iPt];}
		if(iLam==8) {lmean[iPt]=lph_HX[nState-1][rapBin-1][iPt]; StatError_low[iPt]=lph_HX_stat_low[nState-1][rapBin-1][iPt]; StatError_high[iPt]=lph_HX_stat_high[nState-1][rapBin-1][iPt]; SystError[iPt]=lph_HX_syst[nState-1][rapBin-1][iPt];}
		if(iLam==9) {lmean[iPt]=ltp_HX[nState-1][rapBin-1][iPt]; StatError_low[iPt]=ltp_HX_stat_low[nState-1][rapBin-1][iPt]; StatError_high[iPt]=ltp_HX_stat_high[nState-1][rapBin-1][iPt]; SystError[iPt]=ltp_HX_syst[nState-1][rapBin-1][iPt];}
		if(iLam==10) {lmean[iPt]=1000; StatError_low[iPt]=0; StatError_high[iPt]=0; SystError[iPt]=0;}
		if(iLam==11) {lmean[iPt]=1000; StatError_low[iPt]=0; StatError_high[iPt]=0; SystError[iPt]=0;}
		if(iLam==12) {lmean[iPt]=lti_HX[nState-1][rapBin-1][iPt]; StatError_low[iPt]=lti_HX_stat_low[nState-1][rapBin-1][iPt]; StatError_high[iPt]=lti_HX_stat_high[nState-1][rapBin-1][iPt]; SystError[iPt]=lti_HX_syst[nState-1][rapBin-1][iPt];}

		if(iLam==13) {lmean[iPt]=1000; StatError_low[iPt]=0; StatError_high[iPt]=0; SystError[iPt]=0;}
		if(iLam==14) {lmean[iPt]=1000; StatError_low[iPt]=0; StatError_high[iPt]=0; SystError[iPt]=0;}
		if(iLam==15) {lmean[iPt]=1000; StatError_low[iPt]=0; StatError_high[iPt]=0; SystError[iPt]=0;}
		if(iLam==16) {lmean[iPt]=1000; StatError_low[iPt]=0; StatError_high[iPt]=0; SystError[iPt]=0;}
		if(iLam==17) {lmean[iPt]=1000; StatError_low[iPt]=0; StatError_high[iPt]=0; SystError[iPt]=0;}
		if(iLam==18) {lmean[iPt]=1000; StatError_low[iPt]=0; StatError_high[iPt]=0; SystError[iPt]=0;}


		TotalError_low[iPt]=TMath::Sqrt(TMath::Power(StatError_low[iPt],2)+TMath::Power(SystError[iPt],2)); TotalError_high[iPt]=TMath::Sqrt(TMath::Power(StatError_high[iPt],2)+TMath::Power(SystError[iPt],2));

		ptCentre_err_low[iPt]=ptCentre_[iPt]-ptBorders_[iPt];
		ptCentre_err_high[iPt]=ptBorders_[iPt+1]-ptCentre_[iPt];
}

	sprintf(filename,"TGraphResults_%dSUps_CDFStat.root",nState);
	TFile *outfileResStat;
	if(iLam==1&&rapBin==1) outfileResStat = new TFile(filename,"RECREATE");
	else outfileResStat = new TFile(filename,"UPDATE");

	TGraphAsymmErrors *MattResults = new TGraphAsymmErrors(nBinsInPt,ptCentre_,lmean,ptCentre_err_low,ptCentre_err_high,StatError_low,StatError_high);
	MattResults->SetName(GraphName);

	outfileResStat->cd();
	MattResults->Draw("aP");
	MattResults->Write();

	outfileResStat->Write();
	outfileResStat->Close();
	delete outfileResStat;
	outfileResStat = NULL;


	sprintf(filename,"TGraphResults_%dSUps_CDFSyst.root",nState);
	TFile *outfileResSyst;
	if(iLam==1&&rapBin==1) outfileResSyst = new TFile(filename,"RECREATE");
	else outfileResSyst = new TFile(filename,"UPDATE");

	TGraphAsymmErrors *MattSyst = new TGraphAsymmErrors(nBinsInPt,ptCentre_,lmean,ptCentre_err_low,ptCentre_err_high,SystError,SystError);
	MattSyst->SetName(GraphName);

	outfileResSyst->cd();
	MattSyst->Draw("aP");
	MattSyst->Write();

	outfileResSyst->Write();
	outfileResSyst->Close();
	delete outfileResSyst;
	outfileResSyst = NULL;


	sprintf(filename,"TGraphResults_%dSUps_CDFTotal.root",nState);
	TFile *outfileResTotal;
	if(iLam==1&&rapBin==1) outfileResTotal = new TFile(filename,"RECREATE");
	else outfileResTotal = new TFile(filename,"UPDATE");

	TGraphAsymmErrors *MattTotal = new TGraphAsymmErrors(nBinsInPt,ptCentre_,lmean,ptCentre_err_low,ptCentre_err_high,TotalError_low,TotalError_high);
	MattTotal->SetName(GraphName);

	outfileResTotal->cd();
	MattTotal->Draw("aP");
	MattTotal->Write();

	outfileResTotal->Write();
	outfileResTotal->Close();
	delete outfileResTotal;
	outfileResTotal = NULL;

	}
	}
	}

	return 0;
}
