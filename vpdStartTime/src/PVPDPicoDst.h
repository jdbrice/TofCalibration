//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Jul  2 16:57:39 2008 by ROOT version 5.12/00h-rc2
// from TTree pvpd/tofr timing
// found on file: ../ntupleAll_5x.root
//////////////////////////////////////////////////////////

#ifndef PVPDPicoDst_h
#define PVPDPicoDst_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

class PVPDPicoDst {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leave types
   Int_t           run;
   Int_t           evt;
   Int_t           trgword;
   Float_t         vertexX;
   Float_t         vertexY;
   Float_t         vertexZ;
   Int_t           vpdEast;
   Int_t           vpdWest;
   Int_t           numberOfVpdEast;
   Int_t           numberOfVpdWest;
   Float_t         tDiff;
   Double_t        pvpdLeadingEdgeTimeEast[19];
   Double_t        pvpdLeadingEdgeTimeWest[19];
   Double_t        pvpdTotEast[19];
   Double_t        pvpdTotWest[19];

   // List of branches
   TBranch        *b_run;   //!
   TBranch        *b_evt;   //!
   TBranch        *b_trgword;   //!
   TBranch        *b_vertexX;   //!
   TBranch        *b_vertexY;   //!
   TBranch        *b_vertexZ;   //!
   TBranch        *b_vpdEast;   //!
   TBranch        *b_vpdWest;   //!
   TBranch        *b_numberOfVpdEast;   //!
   TBranch        *b_numberOfVpdWest;   //!
   TBranch        *b_tDiff;   //!
   TBranch        *b_pvpdLeadingEdgeTimeEast;   //!
   TBranch        *b_pvpdLeadingEdgeTimeWest;   //!
   TBranch        *b_pvpdTotEast;   //!
   TBranch        *b_pvpdTotWest;   //!

   PVPDPicoDst(TTree *tree=0);
   virtual ~PVPDPicoDst();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef PVPDPicoDst_cxx
PVPDPicoDst::PVPDPicoDst(TTree *tree)
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("../ntupleAll_5x.root");
      if (!f) {
         f = new TFile("../ntupleAll_5x.root");
      }
      tree = (TTree*)gDirectory->Get("pvpd");

   }
   Init(tree);
}

PVPDPicoDst::~PVPDPicoDst()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t PVPDPicoDst::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t PVPDPicoDst::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->IsA() != TChain::Class()) return centry;
   TChain *chain = (TChain*)fChain;
   if (chain->GetTreeNumber() != fCurrent) {
      fCurrent = chain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void PVPDPicoDst::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normaly not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("run", &run, &b_run);
   fChain->SetBranchAddress("evt", &evt, &b_evt);
   fChain->SetBranchAddress("trgword", &trgword, &b_trgword);
   fChain->SetBranchAddress("vertexX", &vertexX, &b_vertexX);
   fChain->SetBranchAddress("vertexY", &vertexY, &b_vertexY);
   fChain->SetBranchAddress("vertexZ", &vertexZ, &b_vertexZ);
   fChain->SetBranchAddress("vpdEast", &vpdEast, &b_vpdEast);
   fChain->SetBranchAddress("vpdWest", &vpdWest, &b_vpdWest);
   fChain->SetBranchAddress("numberOfVpdEast", &numberOfVpdEast, &b_numberOfVpdEast);
   fChain->SetBranchAddress("numberOfVpdWest", &numberOfVpdWest, &b_numberOfVpdWest);
   fChain->SetBranchAddress("tDiff", &tDiff, &b_tDiff);
   fChain->SetBranchAddress("pvpdLeadingEdgeTimeEast", pvpdLeadingEdgeTimeEast, &b_pvpdLeadingEdgeTimeEast);
   fChain->SetBranchAddress("pvpdLeadingEdgeTimeWest", pvpdLeadingEdgeTimeWest, &b_pvpdLeadingEdgeTimeWest);
   fChain->SetBranchAddress("pvpdTotEast", pvpdTotEast, &b_pvpdTotEast);
   fChain->SetBranchAddress("pvpdTotWest", pvpdTotWest, &b_pvpdTotWest);
   Notify();
}

Bool_t PVPDPicoDst::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normaly not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void PVPDPicoDst::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t PVPDPicoDst::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef PVPDPicoDst_cxx
