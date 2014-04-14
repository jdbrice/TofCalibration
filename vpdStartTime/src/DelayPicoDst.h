//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu Apr 23 21:42:15 2009 by ROOT version 5.12/00h-rc2
// from TTree DelayPicoDst/DelayPicoDst
// found on file: delay_0.root
//////////////////////////////////////////////////////////

#ifndef DelayPicoDst_h
#define DelayPicoDst_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

class DelayPicoDst {
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

   DelayPicoDst(TTree *tree=0);
   virtual ~DelayPicoDst();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef DelayPicoDst_cxx
DelayPicoDst::DelayPicoDst(TTree *tree)
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("delay_0.root");
      if (!f) {
         f = new TFile("delay_0.root");
      }
      tree = (TTree*)gDirectory->Get("DelayPicoDst");

   }
   Init(tree);
}

DelayPicoDst::~DelayPicoDst()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t DelayPicoDst::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t DelayPicoDst::LoadTree(Long64_t entry)
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

void DelayPicoDst::Init(TTree *tree)
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

Bool_t DelayPicoDst::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normaly not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void DelayPicoDst::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t DelayPicoDst::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef DelayPicoDst_cxx
