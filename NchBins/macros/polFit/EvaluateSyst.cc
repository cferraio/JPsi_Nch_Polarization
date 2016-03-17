#include <iostream>
#include <string>
#include <sstream>
using namespace std;

#include "rootIncludes.inc"
#include "commonVar.h"
#include "ToyMC.h"

	template<typename T>
void fromSplit(const std::string& key, const std::string &arg, T& out)
{
	const char delim = '=';
	// Skip if key or delimiter not there
	if ((arg.find(key) == std::string::npos) ||
			(arg.find(delim) == std::string::npos))
		return;

	std::string skey, sval;
	std::stringstream sstr(arg);
	std::getline(sstr, skey, delim); // Dummy read to skip key
	std::getline(sstr, sval, delim); // Get value
	T tout;
	if (!(std::istringstream(sval) >> std::boolalpha >> tout))
		return;
	out = tout;
	std::cout << std::boolalpha << skey << ": "  << out << std::endl;
}

// Special version for string without the conversion
	template<>
void fromSplit(const std::string& key, const std::string &arg, std::string &out)
{
	const char delim = '=';
	// Skip if key or delimiter not there
	if ((arg.find(key) == std::string::npos) ||
			(arg.find(delim) == std::string::npos))
		return;
	std::string skey, sval;
	std::stringstream sstr(arg);
	std::getline(sstr, skey, delim); // Dummy read to skip key
	std::getline(sstr, sval, delim); // Get value
	out = sval;
	std::cout << skey << ": "  << out << std::endl;
}


//--------------------------------------------------------------------
int main(int argc, char* argv[]) {

	// Set defaults
	std::string storagedir,
	 	basedir,
	 	JobID1,
	 	JobID2,
	 	SystDir;
	int ptBinMin = 999,
		 	ptBinMax = 999,
		 	rapBinMin = 999,
		 	rapBinMax = 999,
		 	nState = 999;
	bool statErrConsideration = false,
			 centralsPlusSyst = false,
			 differentErrors = false,
			 sqrt12 = false,
			 removeFirstPoint = false,
			 TU = false;

	// Loop over argument list
	for (int i=1; i < argc; i++)
	{
		std::string arg = argv[i];
		fromSplit("ptBinMin", arg, ptBinMin);
		fromSplit("ptBinMax", arg, ptBinMax);
		fromSplit("rapBinMin", arg, rapBinMin);
		fromSplit("rapBinMax", arg, rapBinMax);
		fromSplit("nState", arg, nState);
		fromSplit("statErrConsideration", arg, statErrConsideration);
		fromSplit("centralsPlusSyst", arg, centralsPlusSyst);
		fromSplit("differentErrors", arg, differentErrors);
		fromSplit("sqrt12", arg, sqrt12);
		fromSplit("removeFirstPoint", arg, removeFirstPoint);
		fromSplit("TU", arg, TU);
		fromSplit("storagedir", arg, storagedir);
		fromSplit("basedir", arg, basedir);
		fromSplit("JobID1", arg, JobID1);
		fromSplit("JobID2", arg, JobID2);
		fromSplit("SystDir", arg, SystDir);
	}

	// input files
	std::stringstream filename1, filename2;
	filename1 << storagedir.c_str() << "/" << JobID1.c_str() << "/TGraphResults_Psi" << nState-3 << "S.root";
	filename2 << storagedir.c_str() << "/" << JobID2.c_str() << "/TGraphResults_Psi" << nState-3 << "S.root";
	TFile *infile1 = new TFile(filename1.str().c_str(),"READ");
	TFile *infile2 = new TFile(filename2.str().c_str(),"READ");
	std::cout << "input file 1: " << filename1.str().c_str() << "\n"
		<< "input file 2: " << filename2.str().c_str() << std::endl;

	// output file
	// delete output file if already there and create new one
	std::stringstream outfilename;
	outfilename << basedir.c_str() << "/macros/polFit/" << SystDir.c_str() << "/TGraphResults_Psi" << nState-3 << "S.root";
	TFile *outfile = new TFile(outfilename.str().c_str(),"RECREATE");

	//high pT rho uncertainty ( pT > 35 GeV ): 
	/////select the largest difference of RellCut to Default in 3 pT and 2 y bin, for each parameter(6) and frame(3)
	double lambdaHighest[18];
	double lambdaHighestVal[18] = {
		0.0712449 , 0.0414536 , 0.0212121 , 0.473577 , 0.0698799 , 0.0878234 , 
		0.0748523 , 0.021685 ,  0.0225437 , 0.878433 , 0.0601865 , 0.0896022 , 
		0.076699 ,  0.0226398 , 0.0230012 , 0.876299 , 0.0595996 , 0.105166
	};

	//low pT rho uncertainty ( 1S: pT < 35 GeV, 2S: all pT ): 
	////select the 1*rms of histogram integrating 7 pt x 2 rap bins of 1S, 4 pt x 3 rap bins of 2S
	double lambdaLowPt[18] = {
		0.0293649 , 0.0252533  , 0.0139001 , 0.04012   , 0.00448195 , 0.0482119 ,
		0.0443217 , 0.0120432  , 0.0149671 , 0.0390512 , 0.00464494 , 0.0463249 ,
		0.0465078 , 0.00993849 , 0.0137137 , 0.0391363 , 0.00473978 , 0.0470386
	};


	for(int iLam = 1; iLam < 19; iLam++){

		for(int rapBin = rapBinMin; rapBin < rapBinMax+1; rapBin++){

			// get graphs from file
			std::stringstream graphName;
			if(iLam==1)  graphName << "lth_CS_rap" << rapBin;
			if(iLam==2)  graphName << "lph_CS_rap" << rapBin;
			if(iLam==3)  graphName << "ltp_CS_rap" << rapBin;
			if(iLam==4)  graphName << "lthstar_CS_rap" << rapBin;
			if(iLam==5)  graphName << "lphstar_CS_rap" << rapBin;
			if(iLam==6)  graphName << "ltilde_CS_rap" << rapBin;

			if(iLam==7)  graphName << "lth_HX_rap" << rapBin;
			if(iLam==8)  graphName << "lph_HX_rap" << rapBin;
			if(iLam==9)  graphName << "ltp_HX_rap" << rapBin;
			if(iLam==10) graphName << "lthstar_HX_rap" << rapBin;
			if(iLam==11) graphName << "lphstar_HX_rap" << rapBin;
			if(iLam==12) graphName << "ltilde_HX_rap" << rapBin;

			if(iLam==13) graphName << "lth_PX_rap" << rapBin;
			if(iLam==14) graphName << "lph_PX_rap" << rapBin;
			if(iLam==15) graphName << "ltp_PX_rap" << rapBin;
			if(iLam==16) graphName << "lthstar_PX_rap" << rapBin;
			if(iLam==17) graphName << "lphstar_PX_rap" << rapBin;
			if(iLam==18) graphName << "ltilde_PX_rap" << rapBin;

			std::cout << graphName.str().c_str() << std::endl;

			int nFrame = 0;
			if(iLam > 0 && iLam < 7) nFrame = 1;
			if(iLam > 6 && iLam < 13) nFrame = 2;
			if(iLam > 12 && iLam < 19) nFrame = 3;

			TGraphAsymmErrors* graph1 = (TGraphAsymmErrors*) infile1->Get(graphName.str().c_str());
			TGraphAsymmErrors* graph2 = (TGraphAsymmErrors*) infile2->Get(graphName.str().c_str());

			// define arrays for storing values from TGraphs
			int nBinspT = ptBinMax - ptBinMin + 1;
			double ptCentre_[nBinspT],
						 ptCentreErr_low[nBinspT],
						 ptCentreErr_high[nBinspT],
						 lmean[nBinspT],
						 ptCentre1_[nBinspT],
						 ptCentreErr1_low[nBinspT],
						 ptCentreErr1_high[nBinspT],
						 ptCentre2_[nBinspT],
						 ptCentreErr2_low[nBinspT],
						 ptCentreErr2_high[nBinspT],
						 lmean1[nBinspT],
						 lmean2[nBinspT],
						 lmeanErr1_low[nBinspT],
						 lmeanErr1_high[nBinspT],
						 lmeanErr2_low[nBinspT],
						 lmeanErr2_high[nBinspT],
						 value[nBinspT-1],
						 valueErr_low[nBinspT-1],
						 valueErr_high[nBinspT-1],
						 ptCentreVal_[nBinspT-1];

			int pt = 0;

			for(int ptBin = ptBinMin; ptBin < ptBinMax+1; ptBin++) {

				graph1->GetPoint(pt,ptCentre1_[pt],lmean1[pt]);
				ptCentreErr1_high[pt]=graph1->GetErrorXhigh(pt);
				ptCentreErr1_low[pt]=graph1->GetErrorXlow(pt);
				graph2->GetPoint(pt,ptCentre2_[pt],lmean2[pt]);
				ptCentreErr2_high[pt]=graph2->GetErrorXhigh(pt);
				ptCentreErr2_low[pt]=graph2->GetErrorXlow(pt);

				if(lmean1[pt]>998) lmean1[pt]=9999;

				ptCentre_[pt]=(ptCentre1_[pt]+ptCentre2_[pt])/2.;
				ptCentreErr_low[pt]=(ptCentreErr1_low[pt]+ptCentreErr2_low[pt])/2.;
				ptCentreErr_high[pt]=(ptCentreErr1_high[pt]+ptCentreErr2_high[pt])/2.;
				//cout<<"ptCentre1_: "<<ptCentre1_[pt]<<endl;
				//cout<<"ptCentre2_: "<<ptCentre2_[pt]<<endl;
				//cout<<"ptCentreErr1_low: "<<ptCentreErr1_low[pt]<<endl;
				//cout<<"ptCentreErr2_low: "<<ptCentreErr2_low[pt]<<endl;
				//cout<<"ptCentreErr1_high: "<<ptCentreErr1_high[pt]<<endl;
				//cout<<"ptCentreErr2_high: "<<ptCentreErr2_high[pt]<<endl;

				lmeanErr1_high[pt]=graph1->GetErrorYhigh(pt);
				lmeanErr1_low[pt]=graph1->GetErrorYlow(pt);
				lmeanErr2_high[pt]=graph2->GetErrorYhigh(pt);
				lmeanErr2_low[pt]=graph2->GetErrorYlow(pt);

				// default setting: build difference
				lmean[pt]=lmean1[pt]-lmean2[pt];
				//lmean[pt]=fabs(lmean1[pt]-lmean2[pt])/2.;
				//lmean[pt]=fabs(lmean1[pt]-lmean2[pt])/TMath::Sqrt(12);
				//lmean[pt]=fabs(lmean1[pt]-lmean2[pt]);
				//lmean[pt] = lmean[pt] * lmean[pt] ;

				//lmean[pt] = lambdaLowPt[iLam-1] ;
				if(pt>8){
					if(lambdaHighest[iLam-1] < lmean[pt]) lambdaHighest[iLam-1] = lmean[pt];
					//lmean[pt] = lambdaHighestVal[iLam-1];
				}

				//lmean[pt] = lmean[pt] * lmean[pt] ;

				if(TU){
					double error1 = (lmeanErr1_low[pt] + lmeanErr1_high[pt])/2;
					double error2 = (lmeanErr2_low[pt] + lmeanErr2_high[pt])/2;
					if(error1 < error2) lmean[pt] = 0;
					else lmean[pt] = TMath::Sqrt(TMath::Power(error1,2) - TMath::Power(error2,2));
				}
				if(sqrt12) lmean[pt]=lmean[pt]/TMath::Sqrt(12.);
				if(removeFirstPoint && pt > 0){
					value[pt-1] = lmean1[pt];
					valueErr_low[pt-1] = ptCentreErr1_low[pt];
					valueErr_high[pt-1] = ptCentreErr1_high[pt];
					ptCentreVal_[pt-1] = ptCentre1_[pt];

					std::cout << pt << " set to " << pt-1 << ": pt = " << ptCentreVal_[pt-1] << " lambda = " << value[pt-1] << std::endl;
				}
				std::cout << pt << ": pT = " << ptCentre_[pt] << ", lambda = " << lmean[pt] << std::endl;

				double SigComp=1;

				double lmeanBuff=lmean[pt];
				if(statErrConsideration){
					cout << "StatErrCheck" << endl;
					bool getSqrt=true;
					if(lmeanBuff*lmeanBuff-SigComp*SigComp*lmeanErr2_low[pt]*lmeanErr2_low[pt]-SigComp*SigComp*lmeanErr1_high[pt]*lmeanErr1_high[pt]<0 && lmeanBuff < 0) {lmean[pt]=0;getSqrt=false;}
					if(lmeanBuff*lmeanBuff-SigComp*SigComp*lmeanErr1_low[pt]*lmeanErr1_low[pt]-SigComp*SigComp*lmeanErr2_high[pt]*lmeanErr2_high[pt]<0 && lmeanBuff > 0) {lmean[pt]=0;getSqrt=false;}
					if(lmeanBuff==0) {lmean[pt]=0;getSqrt=false;}
					if(lmeanBuff < 0 && getSqrt) lmean[pt]=-TMath::Sqrt(lmeanBuff*lmeanBuff-SigComp*SigComp*lmeanErr2_low[pt]*lmeanErr2_low[pt]-SigComp*SigComp*lmeanErr1_high[pt]*lmeanErr1_high[pt]);
					if(lmeanBuff > 0 && getSqrt) lmean[pt]=TMath::Sqrt(lmeanBuff*lmeanBuff-SigComp*SigComp*lmeanErr1_low[pt]*lmeanErr1_low[pt]-SigComp*SigComp*lmeanErr2_high[pt]*lmeanErr2_high[pt]);
				}

				//if(lmeanBuff < 0) lmean[pt]=lmeanBuff/TMath::Sqrt(lmeanErr2_low[pt]*lmeanErr2_low[pt]+lmeanErr1_high[pt]*lmeanErr1_high[pt]); //ifPull
				//else lmean[pt]=lmeanBuff/TMath::Sqrt(lmeanErr1_low[pt]*lmeanErr1_low[pt]+lmeanErr2_high[pt]*lmeanErr2_high[pt]); //ifPull

				pt++;
			}// ptBin

			TGraphAsymmErrors *graphSyst;
			// centrals plus systematics
			if(centralsPlusSyst) {
				graphSyst = new TGraphAsymmErrors(nBinspT,ptCentre2_,lmean2,ptCentreErr2_low,ptCentreErr2_high,lmean1,lmean1);// ifCentralsWithTotalSyst
				//std::cout << "centralsPlusSyst option" << std::endl;
			}
			// if 'take central value from JobID2, take error from JobID1'
			else if(differentErrors){
				graphSyst = new TGraphAsymmErrors(nBinspT,ptCentre_,lmean2,ptCentreErr_low,ptCentreErr_high,lmeanErr1_low,lmeanErr1_high);
				//std::cout << "differentErrors option" << std::endl;
			}
			else if(removeFirstPoint){
				graphSyst = new TGraphAsymmErrors(nBinspT-1,ptCentreVal_,value,valueErr_low,valueErr_high,0,0);
			}
			// default setting: difference
			else{
				graphSyst = new TGraphAsymmErrors(nBinspT,ptCentre_,lmean,ptCentreErr_low,ptCentreErr_high,0,0);
				//std::cout << "Default settings" << std::endl;
			}

			graphSyst->SetMarkerColor(ToyMC::MarkerColor[nFrame]);
			graphSyst->SetLineColor(ToyMC::MarkerColor[nFrame]);
			graphSyst->SetMarkerStyle(ToyMC::MarkerStyle[nState][rapBin]);
			graphSyst->SetMarkerSize(ToyMC::MarkerSize[nState][rapBin]);
			graphSyst->SetName(graphName.str().c_str());

			outfile->cd();
			graphSyst->Draw("P");
			graphSyst->Write();

		}// iRap
	} // iLam

	for(int iLam=1; iLam<19; iLam++){
		cout<<"lambdaHighest["<<iLam<<"]: "<<lambdaHighest[iLam-1]<<endl;
	}
	outfile->Write();
	std::cout << "new TGraphResults written" << std::endl;
	outfile->Close();
	return 0;

}// main
