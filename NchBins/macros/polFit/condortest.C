{

char rootin[200];
char list[200];
char command[200];
TFile *f1[2][10][50];

for(int pt = 1; pt < 3; pt++)
{
  for(int cpm = 1; cpm < 11; cpm++)
  {
    for(int fit = 1; fit < 26; fit++)
    { 		  sprintf(rootin,"/data/users/ferraioc/Polarization/JPsi/NchBins/19May16_MassUpdateFixedErrBars_FracL25/results_Fit%d_Psi1S_rap1_pT%d_cpm%d.root",fit,pt,cpm);
        
    TFile* f1[pt-1][cpm-1][fit-1] = new TFile(rootin);
	if(f1[pt-1][cpm-1][fit-1]->IsZombie()){
		sprintf(list,"Error: failed to open pt%d_cpm%d_fit%d",pt,cpm,fit);
//		std::cout << "Error: failed to open pt"<<pt<<"_cpm"<<cpm<<"_fit"<<fit << std::endl;
//		cout<<list >> Condortest.txt;
		std::cout <<"/home/ferraioc/PolNew/CMSSW_5_3_20/src/JPsi_Nch_Polarization/NchBins/macros/polFit "<<fit-1<<" condor-simple.py "<<pt<<" "<<cpm<<endl;
		continue;
	}
    
//    cout<<rootin<<endl;
  //  cout<<h_lth_CS->GetEntries()<<endl;
    
//    f1[pt-1][cpm-1][fit-1]->Close();
    }
  }
}

}