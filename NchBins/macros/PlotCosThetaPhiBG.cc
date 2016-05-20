#include <iostream>
#include <string>
#include <sstream>
using namespace std;

#include "rootIncludes.inc"
#include "commonVar.h"

#include "TH2F.h"

enum{L,R};
const char *bgLabel[2] = {"L", "R"};
TH2F *hCosThetaPhi[onia::kNbRapForPTBins][onia::kNbPTMaxBins][onia::NchBins][onia::kNbFrames][2];
TH2F *hCosThetaPhiHighct[onia::kNbRapForPTBins][onia::kNbPTMaxBins][onia::NchBins][onia::kNbFrames][2];
TH2F *hCosThetaPhiNPBG[onia::kNbRapForPTBins][onia::kNbPTMaxBins][onia::NchBins][onia::kNbFrames];
TH2F *hTCosThetaPhi[onia::kNbRapForPTBins][onia::kNbPTMaxBins][onia::NchBins][onia::kNbFrames];
TH1F *events_SR[onia::kNbRapForPTBins][onia::kNbPTMaxBins][onia::NchBins];
TH1F *events_PRSR[onia::kNbRapForPTBins][onia::kNbPTMaxBins][onia::NchBins];
TH1F *mean_pT[onia::kNbRapForPTBins][onia::kNbPTMaxBins][onia::NchBins];
TH1F *mean_y[onia::kNbRapForPTBins][onia::kNbPTMaxBins][onia::NchBins];
TH1F *mean_cpm[onia::kNbRapForPTBins][onia::kNbPTMaxBins][onia::NchBins];
TH2F *nbin_costhCS[onia::kNbRapForPTBins][onia::kNbPTMaxBins][onia::NchBins];
TH2F *nbin_costhPHX[onia::kNbRapForPTBins][onia::kNbPTMaxBins][onia::NchBins];
TH2F *nbin_costhHX[onia::kNbRapForPTBins][onia::kNbPTMaxBins][onia::NchBins];
//TH2F *nbin_phi[onia::kNbRapForPTBins][onia::kNbPTMaxBins][onia::NchBins];

int evtPinPRSR[onia::kNbRapForPTBins][onia::kNbPTMaxBins][onia::NchBins],
		evtNPinPRSR[onia::kNbRapForPTBins][onia::kNbPTMaxBins][onia::NchBins],
		evtBGinPRSR[onia::kNbRapForPTBins][onia::kNbPTMaxBins][onia::NchBins];
double fracBGinPRSR[onia::kNbRapForPTBins][onia::kNbPTMaxBins][onia::NchBins],
			 fracNPinPRSR[onia::kNbRapForPTBins][onia::kNbPTMaxBins][onia::NchBins];
double MeanPt[onia::kNbRapForPTBins][onia::kNbPTMaxBins][onia::NchBins],
			 MeanRap[onia::kNbRapForPTBins][onia::kNbPTMaxBins][onia::NchBins],
			 MeanCpm[onia::kNbRapForPTBins][onia::kNbPTMaxBins][onia::NchBins];
int nbincosthCS[onia::kNbRapForPTBins][onia::kNbPTMaxBins][onia::NchBins],nbinphiCS[onia::kNbRapForPTBins][onia::kNbPTMaxBins][onia::NchBins];
int nbincosthHX[onia::kNbRapForPTBins][onia::kNbPTMaxBins][onia::NchBins],nbinphiHX[onia::kNbRapForPTBins][onia::kNbPTMaxBins][onia::NchBins];
int nbincosthPHX[onia::kNbRapForPTBins][onia::kNbPTMaxBins][onia::NchBins],nbinphiPHX[onia::kNbRapForPTBins][onia::kNbPTMaxBins][onia::NchBins];

void LoadHistos(Int_t iRapBin, Int_t iPTBin, Int_t iCPMBin, Int_t nState);
void PlotHistos(Int_t iRapBin, Int_t iPTBin, Int_t iCPMBin, Int_t iFrame, Int_t iWindow);
//===========================

//========================================================
// code to read input arguments
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


int main(int argc, char* argv[]){

	// set default values
	int nState = 999;

	// Loop over argument list                                                          
	for (int i=1; i < argc; i++){
		std::string arg = argv[i];
		fromSplit("nState", arg, nState);
	}

	gStyle->SetPalette(1);

	for(int iRap = 0; iRap < onia::kNbRapForPTBins; iRap++){
		for(int iPT = 0; iPT < onia::kNbPTBins[iRap+1]; iPT++){
		for(int iCPM = 0; iCPM < onia::NchBins; iCPM++){
			LoadHistos(iRap, iPT, iCPM, nState);
			for(int iFrame = 0; iFrame < onia::kNbFrames; iFrame++){
				if(iFrame<3) PlotHistos(iRap, iPT, iCPM, iFrame, L);
				if(iFrame<3) PlotHistos(iRap, iPT, iCPM, iFrame, R);
				}
			}
		}
	}


	FILE *NumFile;

	bool SaveTables=true;
	if(SaveTables){
		//for(int iRap = 0; iRap < onia::kNbRapForPTBins; iRap++){
		//	for(int iPT = 0; iPT < onia::kNbPTBins[iRap+1]; iPT++){
		//		cout<<"evtPinPRSR  ["<<iRap<<"]["<<iPT<<"]: "<<evtPinPRSR[iRap][iPT]<<endl;
		//		cout<<"evtNPinPRSR ["<<iRap<<"]["<<iPT<<"]: "<<evtNPinPRSR[iRap][iPT]<<endl;
		//		cout<<"evtBGinPRSR ["<<iRap<<"]["<<iPT<<"]: "<<evtBGinPRSR[iRap][iPT]<<endl;
		//		cout<<"fracBGinPRSR["<<iRap<<"]["<<iPT<<"]: "<<fracBGinPRSR[iRap][iPT]<<endl<<endl;
		//	}
		//}

		char framerap[200];
		char NumFileName[200];

		/// estimated number signal events and Bg fractions 
		sprintf(NumFileName,"NumEventsFracBG_Psi%dS.tex",nState-3);
		NumFile = fopen(NumFileName,"w");
		fprintf(NumFile, "\n");
		fprintf(NumFile,"\\documentclass{article}\n\\usepackage[applemac]{inputenc}\n\\usepackage{amsmath}\n\\usepackage{textcomp}\n\\pagestyle{plain}\n\\usepackage{graphicx}\n\\usepackage{multicol}\n\\usepackage{geometry}\n\\usepackage{subfigure}\n\\usepackage{booktabs}\n\\usepackage{setspace}\n\n\n\n\\begin{document}\n");

		fprintf(NumFile, "\n\n\n\n");

		fprintf(NumFile, "\n\n\n\n");

		fprintf(NumFile, "\\begin{table}\n\\centering\n \\caption{Number of prompt signal events, non-prompt fraction and background fraction in the PRSR, as a function of $p_{T}$ and $y$. The fractions are in \\% . } \n");
		fprintf(NumFile, "\\begin{tabular}{|c|ccc|ccc|ccc|}\n\\hline\n");
		fprintf(NumFile, "$N_{ch}$  & $N_{PR}^{PRSR}$ & $f_{NP}^{PRSR} $ & $f_{BG}^{PRSR} $ & $N_{PR}^{PRSR}$ & $f_{NP}^{PRSR} $ & $f_{BG}^{PRSR} $ & $N_{PR}^{PRSR}$ & $f_{NP}^{PRSR} $ & $f_{BG}^{PRSR} $  \\\\\n");


		if(nState==4){
			sprintf(framerap,"\\hline \\multicolumn{1}{|c|}{} & \\multicolumn{3}{|c|}{$15 < p_{T} < 25$ GeV} & \\multicolumn{3}{|c|}{$25 < p_{T} < 50$ GeV} & \\multicolumn{3}{|c|}{ } \\\\\n \\hline \n");
			fprintf(NumFile,framerap);
			fprintf(NumFile, "\\multicolumn{10}{|c|}{$\\Psi(1S)$} \\\\\n \\hline \n \\rule{0pt}{4mm} \n");
			int pt=0;
			int cpm=0;
			for(int cpmBin = 1; cpmBin < onia::NchBins+1; cpmBin++) {			
				fprintf(NumFile, "%1.0f--%1.0f   &  $%d  $  & $%1.1f $  &  $%1.1f $ &  $%d $ &  $%1.1f $ &  $%1.1f $ & -- & -- & -- \\\\\n", 
						onia::cpmRange[cpmBin-1], onia::cpmRange[cpmBin],
						evtPinPRSR[0][0][cpmBin-1], 100.*fracNPinPRSR[0][0][cpmBin-1], 100*fracBGinPRSR[0][0][cpmBin-1],
						evtPinPRSR[0][1][cpmBin-1], 100.*fracNPinPRSR[0][1][cpmBin-1], 100*fracBGinPRSR[0][1][cpmBin-1]);
						cpm++;
						}
		}

		if(nState==5){
			sprintf(framerap,"\\hline \\multicolumn{1}{|c|}{} & \\multicolumn{3}{|c|}{$%1.1f < |y| < %1.1f$ } & \\multicolumn{3}{|c|}{$%1.1f < |y| < %1.1f$} & \\multicolumn{3}{|c|}{$%1.1f < |y| < %1.1f$} \\\\\n \\hline \n",onia::rapForPTRange[0],onia::rapForPTRange[1], onia::rapForPTRange[1], onia::rapForPTRange[2], onia::rapForPTRange[2], onia::rapForPTRange[3]); 
			fprintf(NumFile,framerap);
			fprintf(NumFile, "\\multicolumn{10}{|c|}{$\\Psi(2S)$} \\\\\n \\hline \n \\rule{0pt}{4mm} \n");

			int pt=0;
			int cpm=0;
			for(int cpmBin = 1; cpmBin < onia::NchBins+1; cpmBin++) {
				fprintf(NumFile, "%1.0f--%1.0f   &  $%d  $  & $%1.1f $  &  $%1.1f $ &  $%d $ &  $%1.1f $ &  $%1.1f $ &  $%d $ &  $%1.1f $ &  $%1.1f $ \\\\\n", 
						onia::cpmRange[cpmBin-1], onia::cpmRange[cpmBin],
						evtPinPRSR[0][0][cpmBin-1], 100.*fracNPinPRSR[0][0][cpmBin-1], 100*fracBGinPRSR[0][0][cpmBin-1],
						evtPinPRSR[1][0][cpmBin-1], 100.*fracNPinPRSR[1][0][cpmBin-1], 100*fracBGinPRSR[1][0][cpmBin-1], 
						evtPinPRSR[2][0][cpmBin-1], 100.*fracNPinPRSR[2][0][cpmBin-1], 100*fracBGinPRSR[2][0][cpmBin-1] );
				cpm++;
			}
		}
		


		fprintf(NumFile, "\\hline\n");
		fprintf(NumFile, "\\end{tabular}\n");
		fprintf(NumFile, "\\label{tab:NumEventsFracBG}\n");
		fprintf(NumFile, "\\end{table}\n");
		fprintf(NumFile, "\n");

		fprintf(NumFile, "\\end{document}");

		fclose(NumFile);

		//// estimated mean Nch, pT, and y 
		sprintf(NumFileName,"meanPt_Psi%dS.tex",nState-3);
		NumFile = fopen(NumFileName,"w");
		fprintf(NumFile, "\n");
		fprintf(NumFile,"\\documentclass{article}\n\\usepackage[applemac]{inputenc}\n\\usepackage{amsmath}\n\\usepackage{textcomp}\n\\pagestyle{plain}\n\\usepackage{graphicx}\n\\usepackage{multicol}\n\\usepackage{geometry}\n\\usepackage{subfigure}\n\\usepackage{booktabs}\n\\usepackage{setspace}\n\n\n\n\\begin{document}\n");

		fprintf(NumFile, "\n\n\n\n");

		fprintf(NumFile, "\n\n\n\n");

		fprintf(NumFile, "\\begin{table}\n\\centering\n \\caption{Estimated mean $p_{T}$ and $|y|$, for each $\\psi(nS)$ kinematic bin} \n");
		fprintf(NumFile, "\\begin{tabular}{|c|ccc|ccc|ccc|}\n\\hline\n");
		fprintf(NumFile, "$N_{ch}$ & $\\hat{N_{ch}}$ & $\\hat{p_{T}}$ [GeV] & $\\hat{|y|}$ & $\\hat{N_{ch}}$ & $\\hat{p_{T}}$ [GeV] & $\\hat{|y|}$ & $\\hat{N_{ch}}$ & $\\hat{p_{T}}$ [GeV] & $\\hat{|y|}$ \\\\\n");

		if(nState==4){
			sprintf(framerap,"\\hline \\multicolumn{1}{|c|}{} & \\multicolumn{3}{|c|}{$%d < p_T < %d$ } & \\multicolumn{3}{|c|}{$%d < p_T < %d$ } & \\multicolumn{3}{|c|}{$%d < p_T < %d$ } \\\\\n \\hline \n",14, 25, 25, 50, 0, 0);
			fprintf(NumFile,framerap);
			fprintf(NumFile, "\\multicolumn{10}{|c|}{$\\Psi(1S)$} \\\\\n \\hline \n \\rule{0pt}{4mm} \n");
			int pt=0;
			int cpm=0;
			for(int cpmBin = 1; cpmBin < onia::NchBins+1; cpmBin++) {

				fprintf(NumFile, "%1.0f--%1.0f & $%1.3f$ & $%1.3f$ & $%1.3f$ & $%1.3f$ & $%1.3f$ & $%1.3f$ & -- & -- & -- \\\\\n", 
						onia::cpmRange[cpmBin-1], onia::cpmRange[cpmBin],
						MeanCpm[0][0][cpmBin-1],MeanPt[0][0][cpmBin-1], MeanRap[0][0][cpmBin-1],
						MeanCpm[0][1][cpmBin-1],MeanPt[0][1][cpmBin-1], MeanRap[0][1][cpmBin-1]);
				cpm++;
				}
		}

		if(nState==5){
			sprintf(framerap,"\\hline \\multicolumn{1}{|c|}{} & \\multicolumn{2}{|c|}{$%1.1f < |y| < %1.1f$ } & \\multicolumn{2}{|c|}{$%1.1f < |y| < %1.1f$} & \\multicolumn{2}{|c|}{$%1.1f < |y| < %1.1f$} \\\\\n \\hline \n",onia::rapForPTRange[0],onia::rapForPTRange[1], onia::rapForPTRange[1], onia::rapForPTRange[2], onia::rapForPTRange[2], onia::rapForPTRange[3]); 
			fprintf(NumFile,framerap);
			fprintf(NumFile, "\\multicolumn{7}{|c|}{$\\Psi(2S)$} \\\\\n \\hline \n \\rule{0pt}{4mm} \n");

			int pt=0;
			int cpm=0;
			for(int ptBin = 1; ptBin < onia::kNbPTBins[1]+1; ptBin++) {
			for(int cpmBin = 1; cpmBin < onia::NchBins+1; cpmBin++) {

				fprintf(NumFile, "%1.0f--%1.0f & $%1.3f$ & $%1.3f$ & $%1.3f$ & $%1.3f$ & $%1.3f$ & $%1.3f$ \\\\\n", 
						onia::pTRange[1][ptBin-1], onia::pTRange[1][ptBin],
						MeanPt[0][ptBin-1][cpmBin+1], MeanRap[0][ptBin-1][cpmBin+1],
						MeanPt[1][ptBin-1][cpmBin+1], MeanRap[1][ptBin-1][cpmBin+1],
						MeanPt[2][ptBin-1][cpmBin+1], MeanRap[2][ptBin-1][cpmBin+1]);
				cpm++;
				}
				pt++;
			}
		}

		fprintf(NumFile, "\\hline\n");
		fprintf(NumFile, "\\end{tabular}\n");
		fprintf(NumFile, "\\label{tab:NumEventsFracBG}\n");
		fprintf(NumFile, "\\end{table}\n");
		fprintf(NumFile, "\n");

		fprintf(NumFile, "\\end{document}");

		fclose(NumFile);


	}
	
			if(nState==4){
		   
		   cout<<endl;
	  	   cout<<endl;
		   cout<<"double ptCentre[nRapBins][nPtBins]={{";
	   
		   for(int iRap = 0; iRap < onia::kNbRapForPTBins; iRap++){
			 Int_t max_pt = onia::kNbPTBins[iRap]-1;
			 for(int iPT = 0; iPT <= max_pt; iPT++){
			 Int_t max_cpm = onia::NchBins-1;
			 for(int icpm = 0; icpm <= max_cpm; icpm++){      
				 cout<<MeanPt[iRap][iPT][icpm];
				 if(icpm<max_cpm) {cout<<", "; }
				 else if(icpm==max_cpm && iPT < onia::kNbPTBins[iRap]-1) {cout<<"},{"; }
				 else {cout<<"}};"; }
			 }
			 }
		   }
		 cout<<endl;
		 cout<<endl;
		 
		 cout<<"double cpmCentre[nRapBins][nPtBins]={{";
	   
		   for(int iRap = 0; iRap < onia::kNbRapForPTBins; iRap++){
			 Int_t max_pt = onia::kNbPTBins[iRap]-1;
			 for(int iPT = 0; iPT <= max_pt; iPT++){
			 Int_t max_cpm = onia::NchBins-1;
			 for(int icpm = 0; icpm <= max_cpm; icpm++){     
				 cout<<MeanCpm[iRap][iPT][icpm];
				 if(icpm<max_cpm) {cout<<", "; }
				 else if(icpm==max_cpm && iPT < onia::kNbPTBins[iRap]-1) {cout<<"},{"; }
				 else {cout<<"}};"; }
			 }
			 }
		   }
		 cout<<endl;
		 cout<<endl;
		 
		 cout<<"double meanRap[nRapBins][nPtBins]={{";
	   
		   for(int iRap = 0; iRap < onia::kNbRapForPTBins; iRap++){
			 Int_t max_pt = onia::kNbPTBins[iRap]-1;
			 for(int iPT = 0; iPT <= max_pt; iPT++){
			 Int_t max_cpm = onia::NchBins-1;
			 for(int icpm = 0; icpm <= max_cpm; icpm++){      
				 cout<<MeanRap[iRap][iPT][icpm];
				 if(icpm<max_cpm) {cout<<", "; }
				 else if(icpm==max_cpm && iPT < onia::kNbPTBins[iRap]-1) {cout<<"},{"; }
				 else {cout<<"}};"; }
			 }
			 }
		   }
		 cout<<endl;
		 cout<<endl;
		 
		 
		 cout<<"double fracBackground[nRapBins][nPtBins]={{";
	   
		   for(int iRap = 0; iRap < onia::kNbRapForPTBins; iRap++){
			 Int_t max_pt = onia::kNbPTBins[iRap]-1;
			 for(int iPT = 0; iPT <= max_pt; iPT++){
			 Int_t max_cpm = onia::NchBins-1;
			 for(int icpm = 0; icpm <= max_cpm; icpm++){      
				 cout<<fracBGinPRSR[iRap][iPT][icpm];
				 if(icpm<max_cpm) {cout<<", "; }
				 else if(icpm==max_cpm && iPT < onia::kNbPTBins[iRap]-1) {cout<<"},{"; }
				 else {cout<<"}};"; }
			 }
			 }
		   }
		 cout<<endl;
		 cout<<endl;
		 
		 cout<<"double numEvents[nRapBins][nPtBins]={{";
	   
		   for(int iRap = 0; iRap < onia::kNbRapForPTBins; iRap++){
			 Int_t max_pt = onia::kNbPTBins[iRap]-1;
			 for(int iPT = 0; iPT <= max_pt; iPT++){
			 Int_t max_cpm = onia::NchBins-1;
			 for(int icpm = 0; icpm <= max_cpm; icpm++){      
				 cout<<evtPinPRSR[iRap][iPT][icpm];
				 if(icpm<max_cpm) {cout<<", "; }
				 else if(icpm==max_cpm && iPT < onia::kNbPTBins[iRap]-1) {cout<<"},{"; }
				 else {cout<<"}};"; }
			 }
			 }
		   }
		 cout<<endl;
		 cout<<endl;
		 
		 
		 cout<<"//CS double binCosth[nRapBins][nPtBins]={{";
	   
		   for(int iRap = 0; iRap < onia::kNbRapForPTBins; iRap++){
			 Int_t max_pt = onia::kNbPTBins[iRap]-1;
			 for(int iPT = 0; iPT <= max_pt; iPT++){
			 Int_t max_cpm = onia::NchBins-1;
			 for(int icpm = 0; icpm <= max_cpm; icpm++){      
				 cout<<nbincosthCS[iRap][iPT][icpm];
				 if(icpm<max_cpm) {cout<<", "; }
				 else if(icpm==max_cpm && iPT < onia::kNbPTBins[iRap]-1) {cout<<"},{"; }
				 else {cout<<"}};"; }
			 }
			 }
		   }
		 cout<<endl;
		 cout<<endl;
		 
		 		 
		 cout<<"//CS double binPhi[nRapBins][nPtBins]={{";
	   
		   for(int iRap = 0; iRap < onia::kNbRapForPTBins; iRap++){
			 Int_t max_pt = onia::kNbPTBins[iRap]-1;
			 for(int iPT = 0; iPT <= max_pt; iPT++){
			 Int_t max_cpm = onia::NchBins-1;
			 for(int icpm = 0; icpm <= max_cpm; icpm++){      
				 cout<<nbinphiCS[iRap][iPT][icpm];
				 if(icpm<max_cpm) {cout<<", "; }
				 else if(icpm==max_cpm && iPT < onia::kNbPTBins[iRap]-1) {cout<<"},{"; }
				 else {cout<<"}};"; }
			 }
			 }
		   }
		 cout<<endl;
		 cout<<endl;
		 
		 cout<<"//HX double binCosth[nRapBins][nPtBins]={{";
	   
		   for(int iRap = 0; iRap < onia::kNbRapForPTBins; iRap++){
			 Int_t max_pt = onia::kNbPTBins[iRap]-1;
			 for(int iPT = 0; iPT <= max_pt; iPT++){
			 Int_t max_cpm = onia::NchBins-1;
			 for(int icpm = 0; icpm <= max_cpm; icpm++){      
				 cout<<nbincosthHX[iRap][iPT][icpm];
				 if(icpm<max_cpm) {cout<<", "; }
				 else if(icpm==max_cpm && iPT < onia::kNbPTBins[iRap]-1) {cout<<"},{"; }
				 else {cout<<"}};"; }
			 }
			 }
		   }
		 cout<<endl;
		 cout<<endl;
		 
		 cout<<"//HX double binPhi[nRapBins][nPtBins]={{";
	   
		   for(int iRap = 0; iRap < onia::kNbRapForPTBins; iRap++){
			 Int_t max_pt = onia::kNbPTBins[iRap]-1;
			 for(int iPT = 0; iPT <= max_pt; iPT++){
			 Int_t max_cpm = onia::NchBins-1;
			 for(int icpm = 0; icpm <= max_cpm; icpm++){      
				 cout<<nbinphiHX[iRap][iPT][icpm];
				 if(icpm<max_cpm) {cout<<", "; }
				 else if(icpm==max_cpm && iPT < onia::kNbPTBins[iRap]-1) {cout<<"},{"; }
				 else {cout<<"}};"; }
			 }
			 }
		   }
		 cout<<endl;
		 cout<<endl;
		 
		 cout<<"//PHX double binPhi[nRapBins][nPtBins]={{";
	   
		   for(int iRap = 0; iRap < onia::kNbRapForPTBins; iRap++){
			 Int_t max_pt = onia::kNbPTBins[iRap]-1;
			 for(int iPT = 0; iPT <= max_pt; iPT++){
			 Int_t max_cpm = onia::NchBins-1;
			 for(int icpm = 0; icpm <= max_cpm; icpm++){      
				 cout<<nbinphiPHX[iRap][iPT][icpm];
				 if(icpm<max_cpm) {cout<<", "; }
				 else if(icpm==max_cpm && iPT < onia::kNbPTBins[iRap]-1) {cout<<"},{"; }
				 else {cout<<"}};"; }
			 }
			 }
		   }
		 cout<<endl;
		 cout<<endl;
		 
		 cout<<"//PHX double binCosth[nRapBins][nPtBins]={{";
	   
		   for(int iRap = 0; iRap < onia::kNbRapForPTBins; iRap++){
			 Int_t max_pt = onia::kNbPTBins[iRap]-1;
			 for(int iPT = 0; iPT <= max_pt; iPT++){
			 Int_t max_cpm = onia::NchBins-1;
			 for(int icpm = 0; icpm <= max_cpm; icpm++){      
				 cout<<nbincosthPHX[iRap][iPT][icpm];
				 if(icpm<max_cpm) {cout<<", "; }
				 else if(icpm==max_cpm && iPT < onia::kNbPTBins[iRap]-1) {cout<<"},{"; }
				 else {cout<<"}};"; }
			 }
			 }
		   }
		 cout<<endl;
		 cout<<endl;
		 
	  
//	  int nx = hTBG_cosThetaPhi[iFrame]->GetXaxis()->GetNbins();
//		int ny = hTBG_cosThetaPhi[iFrame]->GetYaxis()->GetNbins();
		
		
  		}


	
	return 1;
}

//===========================
void PlotHistos(Int_t iRapBin, Int_t iPTBin, Int_t iCPMBin, Int_t iFrame, Int_t iWindow){

	TGaxis::SetMaxDigits(3);

	double lvalue = 0.28, tvalue = 0.92;
	double left=lvalue, top=tvalue, textSize=0.035;
	TLatex *latex=new TLatex();
//	latex->SetTextFont(42);
	latex->SetNDC(kTRUE);
	latex->SetTextSize(textSize);
	double step=textSize*1.3;

	gStyle->SetPadRightMargin(0.12);
	gStyle->SetOptStat(0);
	gStyle->SetFrameBorderMode(0);

	Char_t name[500];
	sprintf(name, "c1_%s_rap%d_pT%d__cpm%d_%s", onia::frameLabel[iFrame], iRapBin+1, iPTBin+1, iCPMBin+1, bgLabel[iWindow]);
	TCanvas *c1 = new TCanvas(name, "", 500, 500);
	gStyle->SetPalette(1);
	gPad->SetFillColor(kWhite);
	gPad->SetLeftMargin(0.15);
	gPad->SetRightMargin(0.15);

	double yOffset=1.4;
	hCosThetaPhi[iRapBin][iPTBin][iCPMBin][iFrame][iWindow]->GetYaxis()->SetTitleOffset(yOffset);
	gPad->SetLeftMargin(0.125);

	hCosThetaPhi[iRapBin][iPTBin][iCPMBin][iFrame][iWindow]->Draw("colz");

	if(iRapBin==0) 
		latex->DrawLatex(left,top,Form("|y| < %.1f, %.1f < p_{T} < %.1f GeV",
					onia::rapForPTRange[iRapBin+1],
					onia::pTRange[iRapBin+1][iPTBin],onia::pTRange[iRapBin+1][iPTBin+1]));
	else
		latex->DrawLatex(left,top,Form("%.1f < |y| < %.1f, %.1f < p_{T} < %.1f GeV",
					onia::rapForPTRange[iRapBin],onia::rapForPTRange[iRapBin+1],
					onia::pTRange[iRapBin+1][iPTBin],onia::pTRange[iRapBin+1][iPTBin+1]));

	sprintf(name, "Figures/cosThetaPhi_%s_rap%d_pT%d_cpm%d_%s.pdf", onia::frameLabel[iFrame], iRapBin+1, iPTBin+1, iCPMBin+1, bgLabel[iWindow]);
	c1->Print(name);

	//sprintf(name, "c2_%s_rap%d_pT%d_%s", onia::frameLabel[iFrame], iRapBin+1, iPTBin+1, bgLabel[iWindow]);
	//TCanvas *c2 = new TCanvas(name, "", 500, 500);
	//c2->SetFillColor(kWhite);
	hCosThetaPhiHighct[iRapBin][iPTBin][iCPMBin][iFrame][iWindow]->GetYaxis()->SetTitleOffset(yOffset);
	hCosThetaPhiHighct[iRapBin][iPTBin][iCPMBin][iFrame][iWindow]->Draw("colz");

	if(iRapBin==0) 
		latex->DrawLatex(left,top,Form("|y| < %.1f, %.1f < p_{T} < %.1f GeV",
					onia::rapForPTRange[iRapBin+1],
					onia::pTRange[iRapBin+1][iPTBin],onia::pTRange[iRapBin+1][iPTBin+1]));
	else
		latex->DrawLatex(left,top,Form("%.1f < |y| < %.1f, %.1f < p_{T} < %.1f GeV",
					onia::rapForPTRange[iRapBin],onia::rapForPTRange[iRapBin+1],
					onia::pTRange[iRapBin+1][iPTBin],onia::pTRange[iRapBin+1][iPTBin+1]));

	sprintf(name, "Figures/cosThetaPhi_%s_rap%d_pT%d_cpm%d_highct_%s.pdf", onia::frameLabel[iFrame], iRapBin+1, iPTBin+1, iCPMBin+1, bgLabel[iWindow]);
	c1->Print(name);

	if(iWindow==0){
		//sprintf(name, "c3_%s_rap%d_pT%d", onia::frameLabel[iFrame], iRapBin+1, iPTBin+1);
		//TCanvas *c3 = new TCanvas(name, "", 500, 500);
		//c3->SetFillColor(kWhite);
		hCosThetaPhiNPBG[iRapBin][iPTBin][iCPMBin][iFrame]->GetYaxis()->SetTitleOffset(yOffset);
		hCosThetaPhiNPBG[iRapBin][iPTBin][iCPMBin][iFrame]->Draw("colz");

		if(iRapBin==0) 
			latex->DrawLatex(left,top,Form("|y| < %.1f, %.1f < p_{T} < %.1f GeV",
						onia::rapForPTRange[iRapBin+1],
						onia::pTRange[iRapBin+1][iPTBin],onia::pTRange[iRapBin+1][iPTBin+1]));
		else
			latex->DrawLatex(left,top,Form("%.1f < |y| < %.1f, %.1f < p_{T} < %.1f GeV",
						onia::rapForPTRange[iRapBin],onia::rapForPTRange[iRapBin+1],
						onia::pTRange[iRapBin+1][iPTBin],onia::pTRange[iRapBin+1][iPTBin+1]));

		sprintf(name, "Figures/cosThetaPhiNPBG_%s_rap%d_pT%d_cpm%d.pdf", onia::frameLabel[iFrame], iRapBin+1, iPTBin+1, iCPMBin+1);
		c1->Print(name);
		//delete c3;

		//sprintf(name, "c4_%s_rap%d_pT%d", onia::frameLabel[iFrame], iRapBin+1, iPTBin+1);
		//TCanvas *c4 = new TCanvas(name, "", 500, 500);
		//c4->SetFillColor(kWhite);
		hTCosThetaPhi[iRapBin][iPTBin][iCPMBin][iFrame]->GetYaxis()->SetTitleOffset(yOffset);
		hTCosThetaPhi[iRapBin][iPTBin][iCPMBin][iFrame]->Draw("colz");

		if(iRapBin==0) 
			latex->DrawLatex(left,top,Form("|y| < %.1f, %.1f < p_{T} < %.1f GeV",
						onia::rapForPTRange[iRapBin+1],
						onia::pTRange[iRapBin+1][iPTBin],onia::pTRange[iRapBin+1][iPTBin+1]));
		else
			latex->DrawLatex(left,top,Form("%.1f < |y| < %.1f, %.1f < p_{T} < %.1f GeV",
						onia::rapForPTRange[iRapBin],onia::rapForPTRange[iRapBin+1],
						onia::pTRange[iRapBin+1][iPTBin],onia::pTRange[iRapBin+1][iPTBin+1]));

		sprintf(name, "Figures/cosThetaPhi_%s_rap%d_pT%d_cpm%d_total.pdf", onia::frameLabel[iFrame], iRapBin+1, iPTBin+1, iCPMBin+1);
		c1->Print(name);
		//delete c4;

	}

	delete c1;
	//delete c2;
}

//===========================
void LoadHistos(Int_t iRapBin, Int_t iPTBin, Int_t iCPMBin, Int_t nState){
	cout<<"rap "<<iRapBin+1<<" pt "<<iPTBin+1<<" cpm "<<iCPMBin+1<<endl;

	Char_t name[200];
	sprintf(name, "tmpFiles/data_Psi%dS_rap%d_pT%d_cpm%d.root", nState-3, iRapBin+1, iPTBin+1, iCPMBin+1);

	TFile *fIn = new TFile(name);

	sprintf(name, "mean_pT"); mean_pT[iRapBin][iPTBin][iCPMBin] = (TH1F*) fIn->Get(name);
	sprintf(name, "mean_y");  mean_y [iRapBin][iPTBin][iCPMBin] = (TH1F*) fIn->Get(name);
	sprintf(name, "mean_cpm");  mean_cpm [iRapBin][iPTBin][iCPMBin] = (TH1F*) fIn->Get(name);
	MeanPt [iRapBin][iPTBin][iCPMBin] = mean_pT[iRapBin][iPTBin][iCPMBin]->GetBinContent(1);
	MeanRap[iRapBin][iPTBin][iCPMBin] = mean_y [iRapBin][iPTBin][iCPMBin]->GetBinContent(1);
	MeanCpm[iRapBin][iPTBin][iCPMBin] = mean_cpm[iRapBin][iPTBin][iCPMBin]->GetBinContent(1);
	
	sprintf(name, "background_costhphiCS"); nbin_costhCS[iRapBin][iPTBin][iCPMBin] = (TH2F*) fIn->Get(name);
	nbincosthCS[iRapBin][iPTBin][iCPMBin]=nbin_costhCS[iRapBin][iPTBin][iCPMBin]->GetXaxis()->GetNbins();
	nbinphiCS[iRapBin][iPTBin][iCPMBin]=nbin_costhCS[iRapBin][iPTBin][iCPMBin]->GetYaxis()->GetNbins();
	
	sprintf(name, "background_costhphiHX"); nbin_costhHX[iRapBin][iPTBin][iCPMBin] = (TH2F*) fIn->Get(name);
	nbincosthHX[iRapBin][iPTBin][iCPMBin]=nbin_costhHX[iRapBin][iPTBin][iCPMBin]->GetXaxis()->GetNbins();
	nbinphiHX[iRapBin][iPTBin][iCPMBin]=nbin_costhHX[iRapBin][iPTBin][iCPMBin]->GetYaxis()->GetNbins();
	
	sprintf(name, "background_costhphiPHX"); nbin_costhPHX[iRapBin][iPTBin][iCPMBin] = (TH2F*) fIn->Get(name);
	nbincosthPHX[iRapBin][iPTBin][iCPMBin]=nbin_costhPHX[iRapBin][iPTBin][iCPMBin]->GetXaxis()->GetNbins();
	nbinphiPHX[iRapBin][iPTBin][iCPMBin]=nbin_costhPHX[iRapBin][iPTBin][iCPMBin]->GetYaxis()->GetNbins();

	sprintf(name, "events_SR");
	events_SR[iRapBin][iPTBin][iCPMBin] = (TH1F*) fIn->Get(name);
	cout<<"rap"<<iRapBin+1<<" pt"<<iPTBin+1<<" cpm "<<iCPMBin+1<<" Prompt signal in PR  : "
		<<events_SR[iRapBin][iPTBin][iCPMBin]->GetBinContent(1)<<endl;

	sprintf(name, "events_promptSR");
	events_PRSR[iRapBin][iPTBin][iCPMBin] = (TH1F*) fIn->Get(name);
	cout<<"rap"<<iRapBin+1<<" pt"<<iPTBin+1<<" cpm "<<iCPMBin+1<<" Prompt signal in PRSR: "
		<<events_PRSR[iRapBin][iPTBin][iCPMBin]->GetBinContent(1)<<endl;

	evtPinPRSR [iRapBin][iPTBin][iCPMBin] = (int)events_PRSR[iRapBin][iPTBin][iCPMBin]->GetBinContent(1);
	evtNPinPRSR[iRapBin][iPTBin][iCPMBin] = (int)events_PRSR[iRapBin][iPTBin][iCPMBin]->GetBinContent(2);
	evtBGinPRSR[iRapBin][iPTBin][iCPMBin] = (int)events_PRSR[iRapBin][iPTBin][iCPMBin]->GetBinContent(3);
	fracBGinPRSR[iRapBin][iPTBin][iCPMBin] = (double)evtBGinPRSR[iRapBin][iPTBin][iCPMBin] / 
		( (double)evtPinPRSR [iRapBin][iPTBin][iCPMBin] + (double)evtNPinPRSR[iRapBin][iPTBin][iCPMBin] + (double)evtBGinPRSR[iRapBin][iPTBin][iCPMBin] ) ;
	fracNPinPRSR[iRapBin][iPTBin][iCPMBin] = (double)evtNPinPRSR[iRapBin][iPTBin][iCPMBin] /
		( (double)evtPinPRSR [iRapBin][iPTBin][iCPMBin] + (double)evtNPinPRSR[iRapBin][iPTBin][iCPMBin] + (double)evtBGinPRSR[iRapBin][iPTBin][iCPMBin] ) ;

	for(int iFrame = 0; iFrame < onia::kNbFrames; iFrame++){
		sprintf(name, "hBG_cosThetaPhi_%s_L", onia::frameLabel[iFrame]);
		cout<<"name: "<<name<<endl;
		hCosThetaPhi[iRapBin][iPTBin][iCPMBin][iFrame][L] = (TH2F *) fIn->Get(name);
		//hCosThetaPhi[iRapBin][iPTBin][iFrame][L]->Print();
		sprintf(name, "hCosThetaPhi_%s_cpm%d_pT%d_rap%d_L", onia::frameLabel[iFrame], iCPMBin+1, iRapBin+1, iPTBin+1);
		hCosThetaPhi[iRapBin][iPTBin][iCPMBin][iFrame][L]->SetName(name);

		sprintf(name, "hBG_cosThetaPhi_%s_R", onia::frameLabel[iFrame]);
		hCosThetaPhi[iRapBin][iPTBin][iCPMBin][iFrame][R] = (TH2F *) fIn->Get(name);
		sprintf(name, "hCosThetaPhi_%s_cpm%d_pT%d_rap%d_R", onia::frameLabel[iFrame], iCPMBin+1, iRapBin+1, iPTBin+1);
		hCosThetaPhi[iRapBin][iPTBin][iCPMBin][iFrame][R]->SetName(name);

		sprintf(name, "hBGinNP_cosThetaPhi_%s_L", onia::frameLabel[iFrame]);
		hCosThetaPhiHighct[iRapBin][iPTBin][iCPMBin][iFrame][L] = (TH2F *) fIn->Get(name);
		//hCosThetaPhiHighct[iRapBin][iPTBin][iFrame][L]->Print();
		sprintf(name, "hCosThetaPhi_%s_cpm%d_pT%d_rap%d_highct_L", onia::frameLabel[iFrame], iCPMBin+1, iRapBin+1, iPTBin+1);
		hCosThetaPhiHighct[iRapBin][iPTBin][iCPMBin][iFrame][L]->SetName(name);

		sprintf(name, "hBGinNP_cosThetaPhi_%s_R", onia::frameLabel[iFrame]);
		hCosThetaPhiHighct[iRapBin][iPTBin][iCPMBin][iFrame][R] = (TH2F *) fIn->Get(name);
		sprintf(name, "hCosThetaPhi_%s_cpm%d_pT%d_rap%d_highct_R", onia::frameLabel[iFrame], iCPMBin+1, iRapBin+1, iPTBin+1);
		hCosThetaPhiHighct[iRapBin][iPTBin][iCPMBin][iFrame][R]->SetName(name);

		sprintf(name, "hNPBG_cosThetaPhi_%s", onia::frameLabel[iFrame]);
		hCosThetaPhiNPBG[iRapBin][iPTBin][iCPMBin][iFrame] = (TH2F *) fIn->Get(name);
		//hCosThetaPhiNPBG[iRapBin][iPTBin][iFrame]->Print();
		sprintf(name, "hCosThetaPhiNPBG_%s_cpm%d_pT%d_rap%d", onia::frameLabel[iFrame], iCPMBin+1, iRapBin+1, iPTBin+1);
		hCosThetaPhiNPBG[iRapBin][iPTBin][iCPMBin][iFrame]->SetName(name);

		sprintf(name, "background_costhphi%s", onia::frameLabel[iFrame]);
		hTCosThetaPhi[iRapBin][iPTBin][iCPMBin][iFrame] = (TH2F *) fIn->Get(name);
		//hTCosThetaPhi[iRapBin][iPTBin][iFrame]->Print();
		sprintf(name, "hTCosThetaPhi_%s_cpm%d_pT%d_rap%d", onia::frameLabel[iFrame], iCPMBin+1, iRapBin+1, iPTBin+1);
		hCosThetaPhiNPBG[iRapBin][iPTBin][iCPMBin][iFrame]->SetName(name);

		if(iFrame==2){
			//std::cout<<"PrintBin rap "<<iRapBin+1<<" pt "<<iPTBin+1<<" nBinsCosthBG: "<<
			//	hCosThetaPhi[iRapBin][iPTBin][iFrame][L]->GetNbinsX()<<std::endl;
			//std::cout<<"PrintBin rap "<<iRapBin+1<<" pt "<<iPTBin+1<<" nBinsPhiBG:   "<<
			//	hCosThetaPhi[iRapBin][iPTBin][iFrame][L]->GetNbinsY()<<std::endl;
			//std::cout<<"PrintBin rap "<<iRapBin+1<<" pt "<<iPTBin+1<<" nBinsCosthNPBG: "<<
			//	hCosThetaPhiNPBG[iRapBin][iPTBin][iFrame]->GetNbinsX()<<std::endl;
			//std::cout<<"PrintBin rap "<<iRapBin+1<<" pt "<<iPTBin+1<<" nBinsPhiNPBG:   "<<
			//	hCosThetaPhiNPBG[iRapBin][iPTBin][iFrame]->GetNbinsY()<<std::endl;
			std::cout<<"PrintBin rap "<<iRapBin+1<<" pt "<<iPTBin+1<<" cpm "<<iCPMBin+1<<" nBinsCosthTBG: "<<
				hTCosThetaPhi[iRapBin][iPTBin][iCPMBin][iFrame]->GetNbinsX()<<std::endl;
			std::cout<<"PrintBin rap "<<iRapBin+1<<" pt "<<iPTBin+1<<" cpm "<<iCPMBin+1<<" nBinsPhiTBG:   "<<
				hTCosThetaPhi[iRapBin][iPTBin][iCPMBin][iFrame]->GetNbinsY()<<std::endl;
		}

	}

}
