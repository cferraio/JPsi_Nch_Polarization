
#include <iostream>
#include <string>
#include <sstream>
using namespace std;

#include "rootIncludes.inc"
#include "commonVar.h"
#include "ToyMC.h"

int main(int argc, char** argv) {

	Char_t *storagedir = "Default"; //Storage Directory
	Char_t *basedir = "Default"; //Storage Directory
	Char_t *SystID = "Default";
	Char_t *JobID1 = "Default";

	int ptBinMin=1;
	int ptBinMax=1;
	int nState=1;
	bool Rap2ToRap1=false;
	bool Rap3ToRap1=false;

	for( int i=0;i < argc; ++i ) {

		if(std::string(argv[i]).find("JobID1") != std::string::npos) {char* JobID1char = argv[i]; char* JobID1char2 = strtok (JobID1char, "="); JobID1 = JobID1char2; cout<<"JobID1 = "<<JobID1<<endl;}

		if(std::string(argv[i]).find("SystID") != std::string::npos) {char* SystIDchar = argv[i]; char* SystIDchar2 = strtok (SystIDchar, "="); SystID = SystIDchar2; cout<<"SystID = "<<SystID<<endl;}
		if(std::string(argv[i]).find("storagedir") != std::string::npos) {char* storagedirchar = argv[i]; char* storagedirchar2 = strtok (storagedirchar, "="); storagedir = storagedirchar2; cout<<"storagedir = "<<storagedir<<endl;}
		if(std::string(argv[i]).find("basedir") != std::string::npos) {char* basedirchar = argv[i]; char* basedirchar2 = strtok (basedirchar, "="); basedir = basedirchar2; cout<<"basedir = "<<basedir<<endl;}

		if(std::string(argv[i]).find("ptBinMin") != std::string::npos) {char* ptBinMinchar = argv[i]; char* ptBinMinchar2 = strtok (ptBinMinchar, "p"); ptBinMin = atof(ptBinMinchar2); cout<<"ptBinMin = "<<ptBinMin<<endl;}
		if(std::string(argv[i]).find("ptBinMax") != std::string::npos) {char* ptBinMaxchar = argv[i]; char* ptBinMaxchar2 = strtok (ptBinMaxchar, "p"); ptBinMax = atof(ptBinMaxchar2); cout<<"ptBinMax = "<<ptBinMax<<endl;}
		if(std::string(argv[i]).find("nState") != std::string::npos) {char* nStatechar = argv[i]; char* nStatechar2 = strtok (nStatechar, "p"); nState = atof(nStatechar2); cout<<"nState = "<<nState<<endl;}

		if(std::string(argv[i]).find("Rap2ToRap1=1") != std::string::npos) {Rap2ToRap1=true; cout<<"Rap2ToRap1"<<endl;}
		if(std::string(argv[i]).find("Rap3ToRap1=1") != std::string::npos) {Rap3ToRap1=true; cout<<"Rap3ToRap1"<<endl;}
	}

	char tmpfilename[200];
	sprintf(tmpfilename,"%s/macros/polFit/Systematics/%s/ChangedTGraph/TGraphResults_Psi%dS.root",basedir,SystID,nState-3);
	gSystem->Unlink(tmpfilename);

	char filename[200];
	sprintf(filename,"%s/TGraphResults_Psi%dS.root",JobID1,nState-3);
	TFile *infile = new TFile(filename,"READ");

	char GraphName[200];

	for(int iFrame = 1; iFrame<4; iFrame++){

		int nRapBins = 2;
		if(nState==5) nRapBins = 3;

		for(int rapBin = 1; rapBin < nRapBins+1; rapBin++){

			cout << "rapBin: " << rapBin << endl;
			TGraphAsymmErrors* graph;

			for(int iParam=1; iParam<7; iParam++){
				if(iParam==1&&iFrame==1) sprintf(GraphName,"lth_CS_rap%d",rapBin);
				if(iParam==2&&iFrame==1) sprintf(GraphName,"lph_CS_rap%d",rapBin);
				if(iParam==3&&iFrame==1) sprintf(GraphName,"ltp_CS_rap%d",rapBin);
				if(iParam==4&&iFrame==1) sprintf(GraphName,"lthstar_CS_rap%d",rapBin);
				if(iParam==5&&iFrame==1) sprintf(GraphName,"lphstar_CS_rap%d",rapBin);
				if(iParam==6&&iFrame==1) sprintf(GraphName,"ltilde_CS_rap%d",rapBin);
				if(iParam==1&&iFrame==2) sprintf(GraphName,"lth_HX_rap%d",rapBin);
				if(iParam==2&&iFrame==2) sprintf(GraphName,"lph_HX_rap%d",rapBin);
				if(iParam==3&&iFrame==2) sprintf(GraphName,"ltp_HX_rap%d",rapBin);
				if(iParam==4&&iFrame==2) sprintf(GraphName,"lthstar_HX_rap%d",rapBin);
				if(iParam==5&&iFrame==2) sprintf(GraphName,"lphstar_HX_rap%d",rapBin);
				if(iParam==6&&iFrame==2) sprintf(GraphName,"ltilde_HX_rap%d",rapBin);
				if(iParam==1&&iFrame==3) sprintf(GraphName,"lth_PX_rap%d",rapBin);
				if(iParam==2&&iFrame==3) sprintf(GraphName,"lph_PX_rap%d",rapBin);
				if(iParam==3&&iFrame==3) sprintf(GraphName,"ltp_PX_rap%d",rapBin);
				if(iParam==4&&iFrame==3) sprintf(GraphName,"lthstar_PX_rap%d",rapBin);
				if(iParam==5&&iFrame==3) sprintf(GraphName,"lphstar_PX_rap%d",rapBin);
				if(iParam==6&&iFrame==3) sprintf(GraphName,"ltilde_PX_rap%d",rapBin);
				graph = (TGraphAsymmErrors*)infile->Get(GraphName);

				TGraphAsymmErrors* ChangedGraph = (TGraphAsymmErrors*)graph->Clone();

				sprintf(filename,"%s/macros/polFit/Systematics/%s/ChangedTGraph/TGraphResults_Psi%dS.root",basedir,SystID,nState-3);
				TFile *outfile = new TFile(filename,"UPDATE");

				outfile->cd();
				ChangedGraph->SetName(GraphName);
				ChangedGraph->Write();

				if(Rap2ToRap1 && rapBin==2){
					if(iParam==1&&iFrame==1) sprintf(GraphName,"lth_CS_rap%d",rapBin-1);
					if(iParam==2&&iFrame==1) sprintf(GraphName,"lph_CS_rap%d",rapBin-1);
					if(iParam==3&&iFrame==1) sprintf(GraphName,"ltp_CS_rap%d",rapBin-1);
					if(iParam==4&&iFrame==1) sprintf(GraphName,"lthstar_CS_rap%d",rapBin-1);
					if(iParam==5&&iFrame==1) sprintf(GraphName,"lphstar_CS_rap%d",rapBin-1);
					if(iParam==6&&iFrame==1) sprintf(GraphName,"ltilde_CS_rap%d",rapBin-1);
					if(iParam==1&&iFrame==2) sprintf(GraphName,"lth_HX_rap%d",rapBin-1);
					if(iParam==2&&iFrame==2) sprintf(GraphName,"lph_HX_rap%d",rapBin-1);
					if(iParam==3&&iFrame==2) sprintf(GraphName,"ltp_HX_rap%d",rapBin-1);
					if(iParam==4&&iFrame==2) sprintf(GraphName,"lthstar_HX_rap%d",rapBin-1);
					if(iParam==5&&iFrame==2) sprintf(GraphName,"lphstar_HX_rap%d",rapBin-1);
					if(iParam==6&&iFrame==2) sprintf(GraphName,"ltilde_HX_rap%d",rapBin-1);
					if(iParam==1&&iFrame==3) sprintf(GraphName,"lth_PX_rap%d",rapBin-1);
					if(iParam==2&&iFrame==3) sprintf(GraphName,"lph_PX_rap%d",rapBin-1);
					if(iParam==3&&iFrame==3) sprintf(GraphName,"ltp_PX_rap%d",rapBin-1);
					if(iParam==4&&iFrame==3) sprintf(GraphName,"lthstar_PX_rap%d",rapBin-1);
					if(iParam==5&&iFrame==3) sprintf(GraphName,"lphstar_PX_rap%d",rapBin-1);
					if(iParam==6&&iFrame==3) sprintf(GraphName,"ltilde_PX_rap%d",rapBin-1);
					outfile->cd();
					ChangedGraph->SetName(GraphName);
					ChangedGraph->Write();
				}

				if(Rap3ToRap1 && rapBin==3){
					if(iParam==1&&iFrame==1) sprintf(GraphName,"lth_CS_rap%d",rapBin-2);
					if(iParam==2&&iFrame==1) sprintf(GraphName,"lph_CS_rap%d",rapBin-2);
					if(iParam==3&&iFrame==1) sprintf(GraphName,"ltp_CS_rap%d",rapBin-2);
					if(iParam==4&&iFrame==1) sprintf(GraphName,"lthstar_CS_rap%d",rapBin-2);
					if(iParam==5&&iFrame==1) sprintf(GraphName,"lphstar_CS_rap%d",rapBin-2);
					if(iParam==6&&iFrame==1) sprintf(GraphName,"ltilde_CS_rap%d",rapBin-2);
					if(iParam==1&&iFrame==2) sprintf(GraphName,"lth_HX_rap%d",rapBin-2);
					if(iParam==2&&iFrame==2) sprintf(GraphName,"lph_HX_rap%d",rapBin-2);
					if(iParam==3&&iFrame==2) sprintf(GraphName,"ltp_HX_rap%d",rapBin-2);
					if(iParam==4&&iFrame==2) sprintf(GraphName,"lthstar_HX_rap%d",rapBin-2);
					if(iParam==5&&iFrame==2) sprintf(GraphName,"lphstar_HX_rap%d",rapBin-2);
					if(iParam==6&&iFrame==2) sprintf(GraphName,"ltilde_HX_rap%d",rapBin-2);
					if(iParam==1&&iFrame==3) sprintf(GraphName,"lth_PX_rap%d",rapBin-2);
					if(iParam==2&&iFrame==3) sprintf(GraphName,"lph_PX_rap%d",rapBin-2);
					if(iParam==3&&iFrame==3) sprintf(GraphName,"ltp_PX_rap%d",rapBin-2);
					if(iParam==4&&iFrame==3) sprintf(GraphName,"lthstar_PX_rap%d",rapBin-2);
					if(iParam==5&&iFrame==3) sprintf(GraphName,"lphstar_PX_rap%d",rapBin-2);
					if(iParam==6&&iFrame==3) sprintf(GraphName,"ltilde_PX_rap%d",rapBin-2);
					outfile->cd();
					ChangedGraph->SetName(GraphName);
					ChangedGraph->Write();
				}


				outfile->Write();
				outfile->Close();
				delete outfile;
				outfile = NULL;

			} //iParam


		} //rapBin

		cout<<"Switching Rapidity"<<endl;

	} //iFrame
	cout<<"Switching frame"<<endl;


	return 0;
}
