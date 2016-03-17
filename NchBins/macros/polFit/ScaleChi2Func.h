#ifndef MN_FcnPVLogL_H_
#define MN_FcnPVLogL_H_
#include "Minuit2/FCNBase.h"
#include "TTree.h"
#include "TH1.h"
#include "TF1.h"
#include "TGraphAsymmErrors.h"
#include <iostream>
#include <string>
#include <sstream>
using namespace std;

#if !defined(__CINT__) && !defined(__MAKECINT__)

#include <vector>
#endif

class FcnPVLogL : public ROOT::Minuit2::FCNBase {
public:
//   FcnPVLogL() {}
FcnPVLogL(TGraphAsymmErrors *graphMC, TGraphAsymmErrors *graphDATA, double pTdist, int pTBinsNew) :  myMCgraph(graphMC), myDATAgraph(graphDATA), mypTdist(pTdist), mypTBinsNew(pTBinsNew), theErrorDef(1.) {}
 ~FcnPVLogL() {delete myMCgraph;}
 double Up() const {return theErrorDef;}
 double operator()(const std::vector<double>&) const;
private:
 TGraphAsymmErrors *myMCgraph;
 TGraphAsymmErrors *myDATAgraph;
 double mypTdist;
 int mypTBinsNew;
 double theErrorDef;
};
#endif //MN_FcnPVLogL_H_

double
FcnPVLogL::operator() (const std::vector<double>& pars) const
{

		  //pars[0]=pTshift
		  //pars[1]=pTscale
		  //pars[2]=effshift
		  //pars[3]=effscale

			double ptCentre[mypTBinsNew];
			double effMean[mypTBinsNew];
			double errGraph[mypTBinsNew];

// Set N = The critical part
///////////////////////////////////////////////////////////////////////////////////
			  for(int pTBin=0;pTBin<mypTBinsNew;pTBin++){
					  ptCentre[pTBin]=pTBin*mypTdist;
					  effMean[pTBin]=myMCgraph->Eval(pTBin*mypTdist/pars[1]-pars[0])*(1+pTBin*mypTdist*pars[3])*pars[4]+pars[2];
					  errGraph[pTBin]=0.;
			  }

//\tilde{\epsilon}_{Data}[p_T]=\epsilon_{MC}[\frac{1}{pT_{scale}}\cdot p_T-pT_{shift}]*(1+eff_{scale}\cdot p_T)\cdot \tilde{eff}_{scale}+eff_{shift}
///////////////////////////////////////////////////////////////////////////////////

			 	TGraphAsymmErrors *Ngraph = new TGraphAsymmErrors(mypTBinsNew,ptCentre,effMean,errGraph,errGraph,errGraph,errGraph);

			      char FitOptions[200];
			      sprintf(FitOptions,"EFNRQ");
			      TF1* f1local;
			      int FitRange=2;

// Calculate chi2
			  int nPt_=myDATAgraph->GetN();
			  const int nDim=nPt_;
			  double effDATA[nDim];
			  double pTDATA[nDim];
			  double efferr[nDim];
			  double effModel[nDim];
			  double chi2[nDim];
			  double chiSq=0;

			  for(int i=0;i<nPt_;i++){
			  myDATAgraph->GetPoint(i,pTDATA[i],effDATA[i]);

//		      f1local = new TF1("f1local","pol1",pTDATA[i]-FitRange*myDATAgraph->GetErrorXlow(i),pTDATA[i]+FitRange*myDATAgraph->GetErrorXhigh(i));
//		      Ngraph->Fit("f1local",FitOptions);

//			  effModel[i] = f1local->GetParameter(0) + f1local->GetParameter(1)*pTDATA[i];// + f1local->GetParameter(2)*pTDATA[i]*pTDATA[i];
//			  if(i==5) cout<<effModel[i]<<endl;
			  effModel[i] = Ngraph->Eval(pTDATA[i]);//(Ngraph->Eval(pTDATA[i])+Ngraph->Eval(pTDATA[i]+pTdist)+Ngraph->Eval(pTDATA[i]-pTdist))/3;
			  if(effDATA[i]-effModel[i]>0)efferr[i]=myDATAgraph->GetErrorYlow(i);
			  else efferr[i]=myDATAgraph->GetErrorYhigh(i);
			  chi2[i]=(effModel[i]-effDATA[i])*(effModel[i]-effDATA[i])/(efferr[i]*efferr[i]);
			  chiSq=chiSq+chi2[i];
			  delete f1local;
			  }

//			  cout<<"pTscaleEst = "<<pars[1]<<", "<<"pTshiftEst = "<<pars[0]<<", "<<"effscaleEst = "<<pars[3]<<", "<<"effshiftEst = "<<pars[2]<<", chiSq = "<<chiSq<<endl;

 return chiSq;
}

