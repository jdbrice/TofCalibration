//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu Apr 23 21:43:10 2009 by ROOT version 5.12/00h-rc2
// from TTree TOTPicoDst/TOTPicoDst
// found on file: totCali_1.root
//////////////////////////////////////////////////////////

#ifndef TOTPicoDst_h
#define TOTPicoDst_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

class TOTPicoDst {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leave types
   Int_t           hits;
   Float_t         tof[8000];   //[hits]
   Float_t         vzVpd;
   Float_t         tag[8000];   //[hits]

   // List of branches
   TBranch        *b_hits;   //!
   TBranch        *b_tof;   //!
   TBranch        *b_vzVpd;   //!
   TBranch        *b_tag;   //!

   TOTPicoDst(TTree *tree=0);
   virtual ~TOTPicoDst();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef TOTPicoDst_cxx
TOTPicoDst::TOTPicoDst(TTree *tree)
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("totCali_1.root");
      if (!f) {
         f = new TFile("totCali_1.root");
      }
      tree = (TTree*)gDirectory->Get("TOTPicoDst");

   }
   Init(tree);
}

TOTPicoDst::~TOTPicoDst()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t TOTPicoDst::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t TOTPicoDst::LoadTree(Long64_t entry)
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

void TOTPicoDst::Init(TTree *tree)
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

   fChain->SetBranchAddress("hits", &hits, &b_hits);
   fChain->SetBranchAddress("tof", tof, &b_tof);
   fChain->SetBranchAddress("vzVpd", &vzVpd, &b_vzVpd);
   fChain->SetBranchAddress("tag", tag, &b_tag);
   Notify();
}

Bool_t TOTPicoDst::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normaly not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void TOTPicoDst::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t TOTPicoDst::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef TOTPicoDst_cxx
