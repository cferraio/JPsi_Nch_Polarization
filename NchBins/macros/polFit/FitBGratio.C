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
#include "TH3.h"
#include "TGraphAsymmErrors.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLatex.h"

void FitBGratio(){

	  gROOT->Reset();
	  gROOT->SetBatch();



	  char PlotDir[200];
	  sprintf(PlotDir,"June7_2DfitsAllCells_Comb_AbsCosth07");
	  char PlotDirFull[200];
	  sprintf(PlotDirFull,"FigBuffer/BGratioFits/%s",PlotDir);
	  gSystem->mkdir(PlotDirFull);
	  char SystDirFull[200];
	  sprintf(SystDirFull,"Systematics/TheGreatRun_BKGratio/%s",PlotDir);
	  gSystem->mkdir(SystDirFull);

	  int FitDimension=2;
	  double AbsCosthMax=0.7;

	  const int nStates=3;
	  int fLSB[nStates+1]={0,72,46,30};

	  for( int nState  = 1; nState  < nStates+1; nState++ ){



		  TFile* inFile;
		  TFile* inFile2;


  TF2 *fcosthphi = new TF2( "fcosthphi", "[0]*(1.+[1]*x[0]*x[0]+[2]*(1.-x[0]*x[0])*cos(2.*x[1]*0.0174532925)+[3]*2.*x[0]*sqrt(1.-x[0]*x[0])*cos(x[1]*0.0174532925))", -AbsCosthMax, AbsCosthMax, -180., 180. );
  fcosthphi->SetParameters(1., 0.0, 0.0, 0.0);

  TF1 *fcosthphi1D = new TF1( "fcosthphi1D", "[0]*(1.+[1]*x[0]*x[0])", -AbsCosthMax, AbsCosthMax);
  fcosthphi1D->SetParameters(1., 0.0);

  const int npT = 11;//11
  const int nrap = 3;//3
  const int nframe = 4;//4

  double  lth[npT][nrap][nframe],   dlth[npT][nrap][nframe];
  double  lph[npT][nrap][nframe],   dlph[npT][nrap][nframe];
  double  lthph[npT][nrap][nframe], dlthph[npT][nrap][nframe];

  char histoName[50];

  TH2D* histoBGmodel;
  TH2D* histoBGdata;
  TH2D* histoBGratio;

  TH1D* histoBGmodel1D;
  TH1D* histoBGdata1D;
  TH1D* histoBGratio1D;

  TH1D* histoBGratio1D_rap1_pT0;
  TH1D* histoBGratio1D_rap2_pT0;

  TF1* Clone_rap1_pT0;
  TF1* Clone_rap2_pT0;

  TLatex *fitBGratio1DtexSumm1;
  TLatex *fitBGratio1DtexSumm2;

  TCanvas *c3;

  for( int iframe  = 1; iframe  < nframe; iframe++ ){


  for( int irap  = 0; irap  < nrap; irap++ ){
  for( int ipT   = 0; ipT   < npT;  ipT++  ){

	  char inFileName[500];
//	  sprintf(inFileName,"/afs/hephy.at/scratch/k/knuenz/CMSSW_4_2_4_patch2/src/UpsilonPol/macros/DataFiles/SetOfCuts11_HighCtau3_BGratioTest/AllStates_1.00Sigma_FracLSB%dPercent/data_%dSUps_rap%d_pT%d.root",fLSB[nState],nState,irap,ipT);
	  sprintf(inFileName,"/afs/hephy.at/scratch/k/knuenz/CMSSW_4_2_4_patch2/src/UpsilonPol/macros/DataFiles/SetOfCuts11_HighCtau3_BGratioTestComb/AllStates_1.00Sigma_FracLSB56Percent/data_%dSUps_rap%d_pT%d.root",nState,irap,ipT);
	  inFile = new TFile(inFileName,"READ");

	  if(iframe==1)    sprintf ( histoName, "background_costhphiCS");
	  if(iframe==2)    sprintf ( histoName, "background_costhphiHX");
	  if(iframe==3)    sprintf ( histoName, "background_costhphiPHX");

    cout << histoName << endl;

    histoBGmodel  = (TH2D*)inFile->Get( histoName );
    histoBGmodel->Sumw2();

    for(int iX=0;iX<histoBGmodel->GetNbinsX()+1;iX++){
        for(int iY=0;iY<histoBGmodel->GetNbinsY()+1;iY++){
        	int globalBin=histoBGmodel->GetBin(iX);
        	if(TMath::Abs(histoBGmodel->GetBinCenter(globalBin))>AbsCosthMax) {
        		histoBGmodel->SetBinContent(iX,iY,0.);
        		histoBGmodel->SetBinError(iX,iY,0.);
        	}
        }
    }

    histoBGmodel1D=histoBGmodel->ProjectionX("histoBGmodel1D", 0, -1,"e");
    histoBGmodel1D->Sumw2();

    for(int iX=0;iX<histoBGmodel1D->GetNbinsX()+1;iX++){
        	int globalBin=histoBGmodel1D->GetBin(iX);
        	if(TMath::Abs(histoBGmodel1D->GetBinCenter(globalBin))>AbsCosthMax) {
        		histoBGmodel1D->SetBinContent(iX,0.);
        		histoBGmodel1D->SetBinError(iX,0.);
        	}
        }



//	  sprintf(inFileName,"/afs/hephy.at/scratch/k/knuenz/CMSSW_4_2_4_patch2/src/UpsilonPol/macros/DataFiles/SetOfCuts11_HighCtau3_BGratioTest/AllStates_1.00Sigma_FracLSB%dPercent/Figures/Ups%dS/AngDistHist.root",fLSB[nState],nState);
	  sprintf(inFileName,"/afs/hephy.at/scratch/k/knuenz/CMSSW_4_2_4_patch2/src/UpsilonPol/macros/DataFiles/SetOfCuts11_HighCtau3_BGratioTestComb/AllStates_1.00Sigma_FracLSB56Percent/Figures/Ups%dS/AngDistHist.root",nState);
	  inFile2 = new TFile(inFileName,"READ");

	  if(iframe==1)    sprintf ( histoName, "CosthPhi_CS_rap%d_pT%d", irap, ipT );
	  if(iframe==2)    sprintf ( histoName, "CosthPhi_HX_rap%d_pT%d", irap, ipT );
	  if(iframe==3)    sprintf ( histoName, "CosthPhi_PHX_rap%d_pT%d", irap, ipT );

	    histoBGdata  = (TH2D*)inFile2->Get( histoName );
	    histoBGdata->Sumw2();

	    for(int iX=0;iX<histoBGdata->GetNbinsX()+1;iX++){
	        for(int iY=0;iY<histoBGdata->GetNbinsY()+1;iY++){
	        	int globalBin=histoBGdata->GetBin(iX);
	        	if(TMath::Abs(histoBGdata->GetBinCenter(globalBin))>AbsCosthMax) {
	        		histoBGdata->SetBinContent(iX,iY,0.);
	        		histoBGdata->SetBinError(iX,iY,0.);
	        	}
	        }
	    }

	    histoBGdata1D=histoBGdata->ProjectionX("histoBGdata1D", 0, -1,"e");
	    histoBGdata1D->Sumw2();

	    for(int iX=0;iX<histoBGdata1D->GetNbinsX()+1;iX++){
	        	int globalBin=histoBGdata1D->GetBin(iX);
	        	if(TMath::Abs(histoBGdata1D->GetBinCenter(globalBin))>AbsCosthMax) {
	        		histoBGdata1D->SetBinContent(iX,0.);
	        		histoBGdata1D->SetBinError(iX,0.);
	        	}
	        }

	    char savenamePlot[200];
	    gStyle->SetPalette(1);


	    TCanvas *c2 = new TCanvas("c2","c2",1200,800);
	    gPad->SetFillColor(kWhite);

	    c2->cd();

	    if(FitDimension==2){
		cout<<"Mean "<<histoBGdata->GetMean()<<endl;
	    histoBGdata->SetTitle(0);
	    histoBGdata->SetStats(0);
	    histoBGdata->Draw("colz");
	    sprintf(savenamePlot,"FigBuffer/BGratioFits/%s/histoBGdata_%dSUps_rap%d_pT%d.pdf",PlotDir,nState,irap,ipT);
	    if(ipT>5||ipT==0) c2->SaveAs(savenamePlot);

//	    histoBGdata->Print("all");
	    histoBGmodel->SetTitle(0);
	    histoBGmodel->SetStats(0);
	    histoBGmodel->Draw("colz");
	    sprintf(savenamePlot,"FigBuffer/BGratioFits/%s/histoBGmodel_%dSUps_rap%d_pT%d.pdf",PlotDir,nState,irap,ipT);
	    if(ipT>5||ipT==0) c2->SaveAs(savenamePlot);

	    histoBGratio=(TH2D*)histoBGdata->Clone("histoBGratio");
	    histoBGratio->Divide(histoBGmodel);
	    histoBGratio->Scale(100./histoBGratio->Integral());
	    cout<<"Mean "<<histoBGratio->GetMean()<<endl;
	    cout<<"Mean "<<histoBGmodel->GetMean()<<endl;

	    histoBGratio->SetTitle(0);
	    histoBGratio->SetStats(0);
	    histoBGratio->Draw("colz");
	    sprintf(savenamePlot,"FigBuffer/BGratioFits/%s/histoBGratio_%dSUps_rap%d_pT%d.pdf",PlotDir,nState,irap,ipT);
	    if(ipT>5||ipT==0) c2->SaveAs(savenamePlot);
	    histoBGratio->GetXaxis()->SetTitleOffset(1.65);
	    histoBGratio->GetYaxis()->SetTitleOffset(2.);
	    histoBGratio->Draw("lego");
	    sprintf(savenamePlot,"FigBuffer/BGratioFits/%s/histoBGratio_%dSUps_rap%d_pT%d_LEGO.pdf",PlotDir,nState,irap,ipT);
	    if(ipT>5||ipT==0) c2->SaveAs(savenamePlot);

	    histoBGratio->Fit("fcosthphi","0");

	     lth[ipT][irap][iframe]   = fcosthphi->GetParameter(1);
	    dlth[ipT][irap][iframe]   = fcosthphi->GetParError(1);
	     lph[ipT][irap][iframe]   = fcosthphi->GetParameter(2);
	    dlph[ipT][irap][iframe]   = fcosthphi->GetParError(2);
	     lthph[ipT][irap][iframe] = fcosthphi->GetParameter(3);
	    dlthph[ipT][irap][iframe] = fcosthphi->GetParError(3);

	    }

	    c2->cd();

	    if(FitDimension==1){

		    cout<<"Mean "<<histoBGdata1D->GetMean()<<endl;
		    histoBGdata1D->SetTitle(0);
		    histoBGdata1D->SetStats(0);
		    histoBGdata1D->Draw("colz");
		    sprintf(savenamePlot,"FigBuffer/BGratioFits/%s/histoBGdata_%dSUps_rap%d_pT%d.pdf",PlotDir,nState,irap,ipT);
		    if(ipT>5||ipT==0) c2->SaveAs(savenamePlot);

		    histoBGmodel1D->SetTitle(0);
		    histoBGmodel1D->SetStats(0);
		    histoBGmodel1D->Draw("colz");
		    sprintf(savenamePlot,"FigBuffer/BGratioFits/%s/histoBGmodel_%dSUps_rap%d_pT%d.pdf",PlotDir,nState,irap,ipT);
		    if(ipT>5||ipT==0) c2->SaveAs(savenamePlot);

		    histoBGratio1D=(TH1D*)histoBGdata1D->Clone("histoBGratio1D");
		    histoBGratio1D->Divide(histoBGmodel1D);
		    histoBGratio1D->Scale(100./histoBGratio1D->Integral());
		    cout<<"Mean "<<histoBGratio1D->GetMean()<<endl;
		    cout<<"Mean "<<histoBGmodel1D->GetMean()<<endl;

		    histoBGratio1D->SetMarkerStyle(20);
		    histoBGratio1D->SetTitle(0);
		    histoBGratio1D->SetStats(0);
		    histoBGratio1D->Draw("colz");
		    sprintf(savenamePlot,"FigBuffer/BGratioFits/%s/histoBGratio_%dSUps_rap%d_pT%d.pdf",PlotDir,nState,irap,ipT);
		    if(ipT>5||ipT==0) c2->SaveAs(savenamePlot);

		    histoBGratio1D->Fit("fcosthphi1D","0");

		    fcosthphi1D->SetLineColor(kGreen+2);
		    fcosthphi1D->SetLineWidth(2.);
		    fcosthphi1D->SetLineStyle(1);
		    fcosthphi1D->Draw("same");


		     lth[ipT][irap][iframe]   = fcosthphi1D->GetParameter(1);
		    dlth[ipT][irap][iframe]   = fcosthphi1D->GetParError(1);
		     lph[ipT][irap][iframe]   = 0;
		    dlph[ipT][irap][iframe]   = 0;
		     lthph[ipT][irap][iframe] = 0;
		    dlthph[ipT][irap][iframe] = 0;

		    char text[200];
		    sprintf(text,"#lambda_{#theta} = %1.3f #pm %1.3f",lth[ipT][irap][iframe],dlth[ipT][irap][iframe]);
		      TLatex *fitBGratio1Dtex = new TLatex(-0.8,histoBGratio1D->GetMaximum()*0.2,text);
		      fitBGratio1Dtex->SetTextSize(0.05)                                                                                                                                                                                                                                             ;
		      fitBGratio1Dtex->Draw( "same" )                                                                                                                                                                                                                                                 ;

		    sprintf(savenamePlot,"FigBuffer/BGratioFits/%s/histoBGratioWithFit_%dSUps_rap%d_pT%d.pdf",PlotDir,nState,irap,ipT);
		    if(ipT>5||ipT==0) c2->SaveAs(savenamePlot);


/*		    if(irap==1&&ipT==0&&iframe==3){

		    	  c3 = new TCanvas("c3","c3",1200,800);
		    	  gPad->SetFillColor(kWhite);

		    	sprintf(text,"#color[2]{#lambda_{#theta} = %1.3f #pm %1.3f}",lth[ipT][irap][iframe],dlth[ipT][irap][iframe]);
		    	fitBGratio1DtexSumm1 = new TLatex(-0.8,histoBGratio1D->GetMaximum()*0.2,text);
		    	fitBGratio1DtexSumm1->SetTextSize(0.05)                                                                                                                                                                                                                                             ;

		    	histoBGratio1D->SetMarkerColor(kRed);
		    	histoBGratio1D->SetMarkerStyle(20);
		    	histoBGratio1D->Draw();
		    	fcosthphi1D->SetLineColor(kRed);
		    	fcosthphi1D->Draw("same");
			    fitBGratio1DtexSumm1->Draw( "same" )                                                                                                                                                                                                                                                 ;



		    }

		    if(irap==2&&ipT==0&&iframe==3){

			    c3->SaveAs("tmp/Canvas.C");

			    char savesummary[200];
			    sprintf(savesummary,"FigBuffer/BGratioFits/%s/histoBGratioWithFit_%dSUps_IntegratedPt.pdf",PlotDir,nState);
			    c3->SaveAs(savesummary);

		    	histoBGratio1D_rap2_pT0=(TH1D*)histoBGratio1D->Clone("histoBGratio1D_rap2_pT0");
		    	Clone_rap2_pT0=(TF1*)fcosthphi1D->Clone("Clone_rap2_pT0");
		        sprintf(text,"#color[2]{#lambda_{#theta} = %1.3f #pm %1.3f}",lth[ipT][irap][iframe],dlth[ipT][irap][iframe]);
		          fitBGratio1DtexSumm2 = new TLatex(0.15,histoBGratio1D_rap2_pT0->GetMaximum()*0.2,text);
		          fitBGratio1DtexSumm2->SetTextSize(0.05)                                                                                                                                                                                                                                             ;
		    }
*/

	    }


    inFile->Close();
    inFile2->Close();

  }
  }



  }


  cout << endl;

  for( int iframe  = 1; iframe  < nframe; iframe++ ){
  for( int irap  = 0; irap  < nrap; irap++ ){
  for( int ipT   = 0; ipT   < npT;  ipT++  ){

      cout << "pT" << ipT << ", rap" << irap  << " in frame "<<iframe << ":  "
           <<          lth[ipT][irap][iframe]   << " +- " << dlth[ipT][irap][iframe]
           << "   " << lph[ipT][irap][iframe]   << " +- " << dlph[ipT][irap][iframe]
           << "   " << lthph[ipT][irap][iframe] << " +- " << dlthph[ipT][irap][iframe]
           << endl;
  }
  }
  }






	char filename[200];
	sprintf(filename,"/afs/hephy.at/scratch/k/knuenz/CMSSW_4_2_4_patch2/src/UpsilonPol/macros/polFit/Systematics/TheGreatRun_BKGratio/%s/TGraphResults_%dSUps.root",PlotDir,nState);
	gSystem->Unlink(filename);
	cout<<filename<<endl;
	char GraphName[200];


	for(int iLam = 1; iLam<19; iLam++){

	for(int rapBin = 1; rapBin < nrap; rapBin++){

	TFile *outfile = new TFile(filename,"UPDATE");

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



	const int nBinspT=10;
	double ptCentre_[nBinspT]={5.747, 6.50359, 7.60069, 8.47588, 9.59466, 10.9576, 13.7561, 17.7305, 23.5741, 35.862};
	double pTrange[nBinspT+1]={5., 6., 7., 8., 9., 10., 12., 16., 20., 30., 50.};
	double ptCentreErr_low[nBinspT];
	double ptCentreErr_high[nBinspT];
	double lmean[nBinspT];
	double errlmean[nBinspT];

	int pt=0;
	for(int ptBin = 1; ptBin < npT; ptBin++) {


		ptCentreErr_low[pt]=ptCentre_[pt]-pTrange[pt];
		ptCentreErr_high[pt]=pTrange[pt+1]-ptCentre_[pt];

		if(iLam==1) {lmean[pt]=lth[ptBin][rapBin][1]; errlmean[pt]=dlth[ptBin][rapBin][1];}
		if(iLam==2) {lmean[pt]=lph[ptBin][rapBin][1];errlmean[pt]=dlph[ptBin][rapBin][1];}
		if(iLam==3) {lmean[pt]=lthph[ptBin][rapBin][1];errlmean[pt]=dlthph[ptBin][rapBin][1];}
		if(iLam==4 || iLam==5) {lmean[pt]=0;errlmean[pt]=0;}
		if(iLam==6) {lmean[pt]=(lth[ptBin][rapBin][1]+3*lph[ptBin][rapBin][1])/(1-lph[ptBin][rapBin][1]);}

		if(iLam==7) {lmean[pt]=lth[ptBin][rapBin][2]; errlmean[pt]=dlth[ptBin][rapBin][2];}
		if(iLam==8) {lmean[pt]=lph[ptBin][rapBin][2];errlmean[pt]=dlph[ptBin][rapBin][2];}
		if(iLam==9) {lmean[pt]=lthph[ptBin][rapBin][2];errlmean[pt]=dlthph[ptBin][rapBin][2];}
		if(iLam==10 || iLam==11) {lmean[pt]=0;errlmean[pt]=0;}
		if(iLam==12) {lmean[pt]=(lth[ptBin][rapBin][2]+3*lph[ptBin][rapBin][2])/(1-lph[ptBin][rapBin][2]);}

		if(iLam==13) {lmean[pt]=lth[ptBin][rapBin][3]; errlmean[pt]=dlth[ptBin][rapBin][3];}
		if(iLam==14) {lmean[pt]=lph[ptBin][rapBin][3];errlmean[pt]=dlph[ptBin][rapBin][3];}
		if(iLam==15) {lmean[pt]=lthph[ptBin][rapBin][3];errlmean[pt]=dlthph[ptBin][rapBin][3];}
		if(iLam==16 || iLam==17) {lmean[pt]=0;errlmean[pt]=0;}
		if(iLam==18) {lmean[pt]=(lth[ptBin][rapBin][3]+3*lph[ptBin][rapBin][3])/(1-lph[ptBin][rapBin][3]);}

	pt++;
	}


	TGraphAsymmErrors *graphSyst = new TGraphAsymmErrors(nBinspT,ptCentre_,lmean,ptCentreErr_low,ptCentreErr_high,errlmean,errlmean);
	graphSyst->SetMarkerColor(kRed);
	graphSyst->SetLineColor(kRed);
	graphSyst->SetMarkerSize(2);
	graphSyst->SetName(GraphName);

	outfile->cd();
	graphSyst->Draw("P");
	graphSyst->Write();

	outfile->Write();
	outfile->Close();
	delete outfile;
	outfile = NULL;


	}


	}



	  }







}
