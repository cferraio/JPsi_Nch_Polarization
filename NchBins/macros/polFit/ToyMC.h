namespace ToyMC{

double ScenarioSig [3][9]={{-1,-0.5,0,0.5,1,0,0,0.5,1},{0,0,0,0,0,0.5,-0.5,-0.75,-1},{0,0,0,0,0,0,0,0,0}};//lamth_Signal,lamph_Signal,lamthph_Signal

double ScenarioBkg [3][9]={{-1,-0.5,0,-0.25,1,0,0,2,4},{0,0,0,0.25,0,0.5,0.8,0.4,-0.4},{0,0,0,0,0,0,0,0,0}};////lamth_Bkg,lamph_Bkg,lamthph_Bkg


int MarkerStyle[6][4]={{0,0,0,0},{0,33,27,34},{0,20,24,29},{0,21,25,22},
	{0,33,27,34},{0,20,24,29}}; // for each state, rapBin (1= closed, 2=open)
int MarkerColor[6] = {0,1,1,1,1,1};//{0,600,632,418}; // for each frame
double MarkerSize[6][4]={{0,0,0},{0,2.75,2.75,2.75},{0,1.65,1.65,1.65},{0,1.65,1.65,1.65},
	{0,2.75,2.75,2.75},{0,1.65,1.65,1.65}};// for each state, rapBin

const int ncpmBins=2;
const int nPtBins=5;
const int nRapBins=2;


// costh and phi bins from data BG histogram 
int binCosth[nRapBins][nPtBins]={{35, 27, 32, 32, 16},{30, 24, 32, 32, 16}};
int binPhi  [nRapBins][nPtBins]={{16, 16, 16, 16, 16},{16, 16, 16, 16, 16}};

//// ToyBackground
double fracBackground[nRapBins][nPtBins]={{0.01,0.01,0.01,0.01,0.01},{0.01,0.01,0.01,0.01,0.01}};
//// from real data
//double fracBackground[nRapBins][nPtBins]={{0.183, 0.190, 0.222, 0.274, 0.370},{0.258, 0.263, 0.285, 0.357, 0.445},{0.303, 0.327, 0.362, 0.420, 0.524}};

//Psi(2): TPV cuts, prompt, Soft approval dataset after cowboy-fix
double ptCentre[nRapBins][nPtBins][ncpmBins]={{{11.893,11.8943},{15.7294,15.7294},{19.7371,19.7371},{25.1272,25.1272},{36.1102,36.1102}},{{11.8943,11.8943},{15.7294,15.7294},{19.7371,19.7371},{25.1272,25.1272},{36.1102,36.1102}}};
double cpmCentre[nRapBins][nPtBins][ncpmBins]={{{20,50},{20,50},{20,50},{20,50},{20,50}},{{20,50},{20,50},{20,50},{20,50},{20,50}}};
int numEvents[nRapBins][nPtBins]={{48068, 18033, 6802, 3944, 1288},{19481, 7892, 2808, 1777, 609}};//Psi(2S)
double meanRap[nRapBins][nPtBins]={{0.881283, 0.875584, 0.866698, 0.862725, 0.859332},{0,0,0,0,0}};
//Total Number of Signal Events in safe region: 251987

//double fracSignal[nRapBins][nPtBins]={{0.0358924, 0.0334104, 0.0277951, 0.026459, 0.0288096, 0.0335685, 0.041832, 0.0521845, 0.0490234, 0.0390959},{0.0140043, 0.0165193, 0.0167022, 0.0172273, 0.0181807, 0.0198244, 0.0239186, 0.0287584, 0.0290342, 0.0270056}};
double fracSignal[nRapBins][nPtBins]={{4.29656e-14, 2.75335e-14, 2.45359e-14, 2.55351e-14, 4.62963e-14},{0, 0, 0, 0, 0}};

const int nEffs=3;
const int FidCuts=3;

}
