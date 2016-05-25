{

char rootin[200];
char list[200];
char command[200];
TFile *f1[2][10][50];

for(int pt = 1; pt < 3; pt++)
{
  for(int cpm = 1; cpm < 11; cpm++)
  {
    for(int gen = 1; gen < 26; gen++)
    { 		  sprintf(rootin,"/data/users/ferraioc/Polarization/JPsi/NchBins/ToyMC/FrameworkIII_19May2016/Sig_frame1scen3_Bkg_frame1scen4/rap1_pT%d_cpm%d/Generation%d/results.root ",pt,cpm,gen);
        
    TFile* f1[pt-1][cpm-1][gen-1] = new TFile(rootin);
	if(f1[pt-1][cpm-1][gen-1]->IsZombie()){
		sprintf(list,"Error: failed to open pt%d_cpm%d_gen%d",pt,cpm,gen);
		cout<<list<<endl;
		continue;
	}
    
    }
  }
}

}