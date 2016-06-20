{

char rootin[200];
char list[200];
char command[200];
TFile *f1[2][10][50];

//framework II
//Sig_frame3scen1_Bkg_frame1scen3/ Sig_frame3scen4_Bkg_frame1scen3/ 
//Sig_frame3scen2_Bkg_frame1scen3/ Sig_frame3scen5_Bkg_frame1scen3/ 

//framework 
//Sig_frame1scen3_Bkg_frame1scen3

for(int pt = 1; pt < 3; pt++)
{
  for(int cpm = 1; cpm < 11; cpm++)
  {
    for(int gen = 1; gen < 50; gen++)
    { 		  sprintf(rootin,"/data/users/ferraioc/Polarization/JPsi/NchBins/ToyMC/FrameworkII_19May2016/Sig_frame3scen5_Bkg_frame1scen3/rap1_pT%d_cpm%d/Generation%d/results.root ",pt,cpm,gen);
        
    TFile* f1[pt-1][cpm-1][gen-1] = new TFile(rootin);
	if(f1[pt-1][cpm-1][gen-1]->IsZombie()){
		sprintf(list,"Error: failed to open pt%d_cpm%d_gen%d",pt,cpm,gen);
		cout<<pt<<" "<<cpm<<" "<<gen-1<<" "<<"5"<<endl;
		continue;
	}
    
    }
  }
}

}