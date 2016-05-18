#ifndef PolData_h
#define PolData_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include "TLorentzVector.h"

class PolData {
public :
   TTree *fChain;   //!pointer to the analyzed TTree or TChain
   int fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   TLorentzVector  *onia, *muPos, *muNeg;
   int eventNb;
   int runNb;
   int lumiBlock;
   int nPriVtx;
   
   double Jpsict;
   double vertexWeight;
   double JpsictErr;
   double JpsiMassErr;
   double JpsiVprob;
   double JpsiDistM1;
   double JpsiDphiM1;
   double JpsiDrM1;
   double JpsiDistM2;
   double JpsiDphiM2;
   double JpsiDrM2;

   // trigger
   // Jpsi
   int HLT_Dimuon0_Jpsi_v1;
   int HLT_Dimuon0_Jpsi_v2;
   int HLT_Dimuon0_Jpsi_v3;
   int HLT_Dimuon0_Jpsi_v5;
   int HLT_Dimuon0_Jpsi_v6;
   int HLT_Dimuon0_Jpsi_v9;
   int HLT_Dimuon8_Jpsi_v3;
   int HLT_Dimuon8_Jpsi_v4;
   int HLT_Dimuon8_Jpsi_v5;
   int HLT_Dimuon10_Jpsi_v3;
   int HLT_Dimuon10_Jpsi_v4;
   int HLT_Dimuon10_Jpsi_Barrel_v1;
   int HLT_Dimuon10_Jpsi_Barrel_v2;
   int HLT_Dimuon10_Jpsi_Barrel_v3;
   int HLT_Dimuon10_Jpsi_Barrel_v5;
   int HLT_Dimuon10_Jpsi_Barrel_v6;
   int HLT_Dimuon10_Jpsi_Barrel_v9;
   int HLT_Dimuon13_Jpsi_Barrel_v1;
   int HLT_Dimuon13_Jpsi_Barrel_v4;
   
   
   //PsiPrime
   int HLT_Dimuon6p5_Barrel_PsiPrime_v1;
   int HLT_Dimuon7_PsiPrime_v1;
   int HLT_Dimuon7_PsiPrime_v2;
   int HLT_Dimuon7_PsiPrime_v3;
   int HLT_Dimuon7_PsiPrime_v5;
   int HLT_Dimuon9_PsiPrime_v1;
   int HLT_Dimuon9_PsiPrime_v4;
   int HLT_Dimuon11_PsiPrime_v1;
   int HLT_Dimuon11_PsiPrime_v4;
   
   // Upsilon
   int HLT_Dimuon0_Upsilon_v1;
   int HLT_Dimuon0_Upsilon_v2;
   int HLT_Dimuon0_Upsilon_v3;
   int HLT_Dimuon0_Upsilon_Barrel_v1;
   int HLT_Dimuon5_Upsilon_Barrel_v1;
   int HLT_Dimuon5_Upsilon_Barrel_v2;
   int HLT_Dimuon5_Upsilon_Barrel_v3;
   int HLT_Dimuon5_Upsilon_Barrel_v5;
   int HLT_Dimuon7_Upsilon_Barrel_v1;
   int HLT_Dimuon7_Upsilon_Barrel_v4;
   int HLT_Dimuon9_Upsilon_Barrel_v1;
   int HLT_Dimuon9_Upsilon_Barrel_v4;

   PolData(TTree *tree=0);
   virtual ~PolData();
   virtual int Cut(Long64_t entry);
   virtual int GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void Init(TTree *tree);
   virtual void Loop(int selDimuType, bool rejectCowboys, int FidCuts, bool MC, bool RequestTrigger, bool removeEta0p2_0p3, bool cutDeltaREllDpt);
   virtual bool Notify();
   virtual void Show(Long64_t entry = -1);
};

#endif

#ifdef PolData_cxx
PolData::PolData(TTree *tree)
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("TTree_Onia2MuMu_V8_PromptReco_v4.root");
      
      if (!f) {
         f = new TFile("TTree_Onia2MuMu_V8_PromptReco_v4.root");
      }
      
      tree = (TTree*)gDirectory->Get("data");
   }
   
   Init(tree);
}

PolData::~PolData()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

int PolData::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}

Long64_t PolData::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (!fChain->InheritsFrom(TChain::Class()))  return centry;
   TChain *chain = (TChain*)fChain;
   if (chain->GetTreeNumber() != fCurrent) {
      fCurrent = chain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void PolData::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   
   onia = 0;
   muNeg = 0;
   muPos = 0;

   fChain->SetBranchAddress("JpsiP", &onia);
   fChain->SetBranchAddress("muNegP", &muNeg);
   fChain->SetBranchAddress("muPosP", &muPos);
   
   fChain->SetBranchAddress("eventNb", &eventNb);
   fChain->SetBranchAddress("runNb", &runNb);
   fChain->SetBranchAddress("lumiBlock", &lumiBlock);
   fChain->SetBranchAddress("nPriVtx", &nPriVtx);
   
   fChain->SetBranchAddress("Jpsict", &Jpsict);
   fChain->SetBranchAddress("vertexWeight", &vertexWeight);
   fChain->SetBranchAddress("JpsictErr", &JpsictErr);
   fChain->SetBranchAddress("JpsiMassErr", &JpsiMassErr);
   fChain->SetBranchAddress("JpsiVprob", &JpsiVprob);
   fChain->SetBranchAddress("JpsiDistM1", &JpsiDistM1);
   fChain->SetBranchAddress("JpsiDphiM1", &JpsiDphiM1);
   fChain->SetBranchAddress("JpsiDrM1", &JpsiDrM1);
   fChain->SetBranchAddress("JpsiDistM2", &JpsiDistM2);
   fChain->SetBranchAddress("JpsiDphiM2", &JpsiDphiM2);
   fChain->SetBranchAddress("JpsiDrM2", &JpsiDrM2);
 
   fChain->SetBranchAddress("HLT_Dimuon8_Jpsi_v3", &HLT_Dimuon8_Jpsi_v3);
   fChain->SetBranchAddress("HLT_Dimuon8_Jpsi_v4", &HLT_Dimuon8_Jpsi_v4);
   fChain->SetBranchAddress("HLT_Dimuon8_Jpsi_v5", &HLT_Dimuon8_Jpsi_v5);
   fChain->SetBranchAddress("HLT_Dimuon10_Jpsi_v3", &HLT_Dimuon10_Jpsi_v3);
   fChain->SetBranchAddress("HLT_Dimuon10_Jpsi_v4", &HLT_Dimuon10_Jpsi_v4);
   
   
   fChain->SetBranchAddress("HLT_Dimuon0_Jpsi_v1", &HLT_Dimuon0_Jpsi_v1);
   fChain->SetBranchAddress("HLT_Dimuon0_Jpsi_v2", &HLT_Dimuon0_Jpsi_v2);
   fChain->SetBranchAddress("HLT_Dimuon0_Jpsi_v3", &HLT_Dimuon0_Jpsi_v3);
   fChain->SetBranchAddress("HLT_Dimuon0_Jpsi_v5", &HLT_Dimuon0_Jpsi_v5);
   fChain->SetBranchAddress("HLT_Dimuon0_Jpsi_v6", &HLT_Dimuon0_Jpsi_v6);
   fChain->SetBranchAddress("HLT_Dimuon0_Jpsi_v9", &HLT_Dimuon0_Jpsi_v9);
   fChain->SetBranchAddress("HLT_Dimuon10_Jpsi_Barrel_v1", &HLT_Dimuon10_Jpsi_Barrel_v1);
   fChain->SetBranchAddress("HLT_Dimuon10_Jpsi_Barrel_v2", &HLT_Dimuon10_Jpsi_Barrel_v2);
   fChain->SetBranchAddress("HLT_Dimuon10_Jpsi_Barrel_v3", &HLT_Dimuon10_Jpsi_Barrel_v3);
   fChain->SetBranchAddress("HLT_Dimuon10_Jpsi_Barrel_v5", &HLT_Dimuon10_Jpsi_Barrel_v5);
   fChain->SetBranchAddress("HLT_Dimuon10_Jpsi_Barrel_v6", &HLT_Dimuon10_Jpsi_Barrel_v6);
   fChain->SetBranchAddress("HLT_Dimuon10_Jpsi_Barrel_v9", &HLT_Dimuon10_Jpsi_Barrel_v9);
   fChain->SetBranchAddress("HLT_Dimuon13_Jpsi_Barrel_v1", &HLT_Dimuon13_Jpsi_Barrel_v1);
   fChain->SetBranchAddress("HLT_Dimuon13_Jpsi_Barrel_v4", &HLT_Dimuon13_Jpsi_Barrel_v4);

   fChain->SetBranchAddress("HLT_Dimuon6p5_Barrel_PsiPrime_v1", &HLT_Dimuon6p5_Barrel_PsiPrime_v1);
   fChain->SetBranchAddress("HLT_Dimuon7_PsiPrime_v1", &HLT_Dimuon7_PsiPrime_v1);
   fChain->SetBranchAddress("HLT_Dimuon7_PsiPrime_v2", &HLT_Dimuon7_PsiPrime_v2);
   fChain->SetBranchAddress("HLT_Dimuon7_PsiPrime_v3", &HLT_Dimuon7_PsiPrime_v3);
   fChain->SetBranchAddress("HLT_Dimuon7_PsiPrime_v5", &HLT_Dimuon7_PsiPrime_v5);
   fChain->SetBranchAddress("HLT_Dimuon9_PsiPrime_v1", &HLT_Dimuon9_PsiPrime_v1);
   fChain->SetBranchAddress("HLT_Dimuon9_PsiPrime_v4", &HLT_Dimuon9_PsiPrime_v4);
   fChain->SetBranchAddress("HLT_Dimuon11_PsiPrime_v1", &HLT_Dimuon11_PsiPrime_v1);
   fChain->SetBranchAddress("HLT_Dimuon11_PsiPrime_v4", &HLT_Dimuon11_PsiPrime_v4);

   fChain->SetBranchAddress("HLT_Dimuon0_Upsilon_v1", &HLT_Dimuon0_Upsilon_v1);
   fChain->SetBranchAddress("HLT_Dimuon0_Upsilon_v2", &HLT_Dimuon0_Upsilon_v2);
   fChain->SetBranchAddress("HLT_Dimuon0_Upsilon_v3", &HLT_Dimuon0_Upsilon_v3);
   fChain->SetBranchAddress("HLT_Dimuon0_Upsilon_Barrel_v1", &HLT_Dimuon0_Upsilon_Barrel_v1);
   fChain->SetBranchAddress("HLT_Dimuon5_Upsilon_Barrel_v1", &HLT_Dimuon5_Upsilon_Barrel_v1);
   fChain->SetBranchAddress("HLT_Dimuon5_Upsilon_Barrel_v2", &HLT_Dimuon5_Upsilon_Barrel_v2);
   fChain->SetBranchAddress("HLT_Dimuon5_Upsilon_Barrel_v3", &HLT_Dimuon5_Upsilon_Barrel_v3);
   fChain->SetBranchAddress("HLT_Dimuon5_Upsilon_Barrel_v5", &HLT_Dimuon5_Upsilon_Barrel_v5);
   fChain->SetBranchAddress("HLT_Dimuon7_Upsilon_Barrel_v1", &HLT_Dimuon7_Upsilon_Barrel_v1);
   fChain->SetBranchAddress("HLT_Dimuon7_Upsilon_Barrel_v4", &HLT_Dimuon7_Upsilon_Barrel_v4);
   fChain->SetBranchAddress("HLT_Dimuon9_Upsilon_Barrel_v1", &HLT_Dimuon9_Upsilon_Barrel_v1);
   fChain->SetBranchAddress("HLT_Dimuon9_Upsilon_Barrel_v4", &HLT_Dimuon9_Upsilon_Barrel_v4);

   Notify();
}

Bool_t PolData::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void PolData::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
int PolData::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef PolData_cxx


