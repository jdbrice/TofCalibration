//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Fri Apr 24 11:49:14 2009 by ROOT version 5.12/00h-rc2
// from TTree StartPicoDst/StartPicoDst
// found on file: pvpdCali.root
//////////////////////////////////////////////////////////

#ifndef StartPicoDst_h
#define StartPicoDst_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

class StartPicoDst {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leave types
   Int_t           Ieast;
   Int_t           Iwest;
   Double_t        sumEast;
   Double_t        sumWest;
   Double_t        T0;
   Double_t        vzVpd;

   // List of branches
   TBranch        *b_Ieast;   //!
   TBranch        *b_Iwest;   //!
   TBranch        *b_sumEast;   //!
   TBranch        *b_sumWest;   //!
   TBranch        *b_T0;   //!
   TBranch        *b_vzVpd;   //!

   StartPicoDst(TTree *tree=0);
   virtual ~StartPicoDst();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef StartPicoDst_cxx
StartPicoDst::StartPicoDst(TTree *tree)
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("pvpdCali.root");
      if (!f) {
         f = new TFile("pvpdCali.root");
      }
      tree = (TTree*)gDirectory->Get("StartPicoDst");

   }
   Init(tree);
}

StartPicoDst::~StartPicoDst()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t StartPicoDst::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t StartPicoDst::LoadTree(Long64_t entry)
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

void StartPicoDst::Init(TTree *tree)
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

   fChain->SetBranchAddress("Ieast", &Ieast, &b_Ieast);
   fChain->SetBranchAddress("Iwest", &Iwest, &b_Iwest);
   fChain->SetBranchAddress("sumEast", &sumEast, &b_sumEast);
   fChain->SetBranchAddress("sumWest", &sumWest, &b_sumWest);
   fChain->SetBranchAddress("T0", &T0, &b_T0);
   fChain->SetBranchAddress("vzVpd", &vzVpd, &b_vzVpd);
   Notify();
}

Bool_t StartPicoDst::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normaly not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void StartPicoDst::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t StartPicoDst::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef StartPicoDst_cxx
