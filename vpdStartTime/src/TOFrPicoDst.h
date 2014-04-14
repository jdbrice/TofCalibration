//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Apr 22 16:44:52 2009 by ROOT version 5.12/00h-rc2
// from TTree tof/BTof cell data
// found on file: ntupleAllDays.root
//////////////////////////////////////////////////////////

#ifndef TOFrPicoDst_h
#define TOFrPicoDst_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

class TOFrPicoDst {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leave types
   Int_t           run;
   Int_t           evt;
   Float_t         vertexX;
   Float_t         vertexY;
   Float_t         vertexZ;
   Int_t           vpdEast;
   Int_t           vpdWest;
   Int_t           numberOfVpdEast;
   Int_t           numberOfVpdWest;
   Float_t         tDiff;
   Double_t        tStart;
   Float_t         vpdVz;
   Double_t        vpdLeEast[19];
   Double_t        vpdLeWest[19];
   Double_t        vpdTotEast[19];
   Double_t        vpdTotWest[19];
   Int_t           nTofHits;
   Int_t           tray[8000];   //[nTofHits]
   Int_t           module[8000];   //[nTofHits]
   Int_t           cell[8000];   //[nTofHits]
   Double_t        leTime[8000];   //[nTofHits]
   Double_t        tot[8000];   //[nTofHits]
   Float_t         yLocal[8000];   //[nTofHits]
   Float_t         zLocal[8000];   //[nTofHits]
   Float_t         thetaLocal[8000];   //[nTofHits]
   Float_t         xGlobal[8000];   //[nTofHits]
   Float_t         yGlobal[8000];   //[nTofHits]
   Float_t         zGlobal[8000];   //[nTofHits]
   Int_t           trackId[8000];   //[nTofHits]
   Int_t           charge[8000];   //[nTofHits]
   Float_t         pt[8000];   //[nTofHits]
   Float_t         eta[8000];   //[nTofHits]
   Float_t         phi[8000];   //[nTofHits]
   Float_t         dcaX[8000];   //[nTofHits]
   Float_t         dcaY[8000];   //[nTofHits]
   Float_t         dcaZ[8000];   //[nTofHits]
   Float_t         length[8000];   //[nTofHits]
   Int_t           nHits[8000];   //[nTofHits]
   Int_t           nHitsFit[8000];   //[nTofHits]
   Int_t           nHitsDedx[8000];   //[nTofHits]
   Float_t         dedx[8000];   //[nTofHits]
   Float_t         nSigE[8000];   //[nTofHits]
   Float_t         nSigPi[8000];   //[nTofHits]
   Float_t         nSigK[8000];   //[nTofHits]
   Float_t         nSigP[8000];   //[nTofHits]
   Float_t         tofCorr[8000];   //[nTofHits]
   Float_t         beta[8000];   //[nTofHits]
   Double_t        T0;//added
   Double_t        vzVpd;//added
   Int_t           Iwest;//added
   Int_t           Ieast;//added

   // List of branches
   TBranch        *b_run;   //!
   TBranch        *b_evt;   //!
   TBranch        *b_vertexX;   //!
   TBranch        *b_vertexY;   //!
   TBranch        *b_vertexZ;   //!
   TBranch        *b_vpdEast;   //!
   TBranch        *b_vpdWest;   //!
   TBranch        *b_numberOfVpdEast;   //!
   TBranch        *b_numberOfVpdWest;   //!
   TBranch        *b_tDiff;   //!
   TBranch        *b_tStart;   //!
   TBranch        *b_vpdVz;   //!
   TBranch        *b_vpdLeEast;   //!
   TBranch        *b_vpdLeWest;   //!
   TBranch        *b_vpdTotEast;   //!
   TBranch        *b_vpdTotWest;   //!
   TBranch        *b_nTofHits;   //!
   TBranch        *b_tray;   //!
   TBranch        *b_module;   //!
   TBranch        *b_cell;   //!
   TBranch        *b_leTime;   //!
   TBranch        *b_tot;   //!
   TBranch        *b_yLocal;   //!
   TBranch        *b_zLocal;   //!
   TBranch        *b_thetaLocal;   //!
   TBranch        *b_xGlobal;   //!
   TBranch        *b_yGlobal;   //!
   TBranch        *b_zGlobal;   //!
   TBranch        *b_trackId;   //!
   TBranch        *b_charge;   //!
   TBranch        *b_pt;   //!
   TBranch        *b_eta;   //!
   TBranch        *b_phi;   //!
   TBranch        *b_dcaX;   //!
   TBranch        *b_dcaY;   //!
   TBranch        *b_dcaZ;   //!
   TBranch        *b_length;   //!
   TBranch        *b_nHits;   //!
   TBranch        *b_nHitsFit;   //!
   TBranch        *b_nHitsDedx;   //!
   TBranch        *b_dedx;   //!
   TBranch        *b_nSigE;   //!
   TBranch        *b_nSigPi;   //!
   TBranch        *b_nSigK;   //!
   TBranch        *b_nSigP;   //!
   TBranch        *b_tofCorr;   //!
   TBranch        *b_beta;   //!
   TBranch        *b_T0;
   TBranch        *b_vzVpd;
   TBranch        *b_Ieast;
   TBranch        *b_Iwest;


   TOFrPicoDst(TTree *tree=0);
   virtual ~TOFrPicoDst();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef TOFrPicoDst_cxx
TOFrPicoDst::TOFrPicoDst(TTree *tree)
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("ntupleAllDays.root");
      if (!f) {
         f = new TFile("ntupleAllDays.root");
      }
      tree = (TTree*)gDirectory->Get("tof");

   }
   Init(tree);
}

TOFrPicoDst::~TOFrPicoDst()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t TOFrPicoDst::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t TOFrPicoDst::LoadTree(Long64_t entry)
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

void TOFrPicoDst::Init(TTree *tree)
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
   fChain->SetBranchAddress("vertexX", &vertexX, &b_vertexX);
   fChain->SetBranchAddress("vertexY", &vertexY, &b_vertexY);
   fChain->SetBranchAddress("vertexZ", &vertexZ, &b_vertexZ);
   fChain->SetBranchAddress("vpdEast", &vpdEast, &b_vpdEast);
   fChain->SetBranchAddress("vpdWest", &vpdWest, &b_vpdWest);
   fChain->SetBranchAddress("numberOfVpdEast", &numberOfVpdEast, &b_numberOfVpdEast);
   fChain->SetBranchAddress("numberOfVpdWest", &numberOfVpdWest, &b_numberOfVpdWest);
   fChain->SetBranchAddress("tDiff", &tDiff, &b_tDiff);
   fChain->SetBranchAddress("tStart", &tStart, &b_tStart);
   fChain->SetBranchAddress("vpdVz", &vpdVz, &b_vpdVz);
   fChain->SetBranchAddress("vpdLeEast", vpdLeEast, &b_vpdLeEast);
   fChain->SetBranchAddress("vpdLeWest", vpdLeWest, &b_vpdLeWest);
   fChain->SetBranchAddress("vpdTotEast", vpdTotEast, &b_vpdTotEast);
   fChain->SetBranchAddress("vpdTotWest", vpdTotWest, &b_vpdTotWest);
   fChain->SetBranchAddress("nTofHits", &nTofHits, &b_nTofHits);
   fChain->SetBranchAddress("tray", tray, &b_tray);
   fChain->SetBranchAddress("module", module, &b_module);
   fChain->SetBranchAddress("cell", cell, &b_cell);
   fChain->SetBranchAddress("leTime", leTime, &b_leTime);
   fChain->SetBranchAddress("tot", tot, &b_tot);
   fChain->SetBranchAddress("yLocal", yLocal, &b_yLocal);
   fChain->SetBranchAddress("zLocal", zLocal, &b_zLocal);
   fChain->SetBranchAddress("thetaLocal", thetaLocal, &b_thetaLocal);
   fChain->SetBranchAddress("xGlobal", xGlobal, &b_xGlobal);
   fChain->SetBranchAddress("yGlobal", yGlobal, &b_yGlobal);
   fChain->SetBranchAddress("zGlobal", zGlobal, &b_zGlobal);
   fChain->SetBranchAddress("trackId", trackId, &b_trackId);
   fChain->SetBranchAddress("charge", charge, &b_charge);
   fChain->SetBranchAddress("pt", pt, &b_pt);
   fChain->SetBranchAddress("eta", eta, &b_eta);
   fChain->SetBranchAddress("phi", phi, &b_phi);
   fChain->SetBranchAddress("dcaX", dcaX, &b_dcaX);
   fChain->SetBranchAddress("dcaY", dcaY, &b_dcaY);
   fChain->SetBranchAddress("dcaZ", dcaZ, &b_dcaZ);
   fChain->SetBranchAddress("length", length, &b_length);
   fChain->SetBranchAddress("nHits", nHits, &b_nHits);
   fChain->SetBranchAddress("nHitsFit", nHitsFit, &b_nHitsFit);
   fChain->SetBranchAddress("nHitsDedx", nHitsDedx, &b_nHitsDedx);
   fChain->SetBranchAddress("dedx", dedx, &b_dedx);
   fChain->SetBranchAddress("nSigE", nSigE, &b_nSigE);
   fChain->SetBranchAddress("nSigPi", nSigPi, &b_nSigPi);
   fChain->SetBranchAddress("nSigK", nSigK, &b_nSigK);
   fChain->SetBranchAddress("nSigP", nSigP, &b_nSigP);
   fChain->SetBranchAddress("tofCorr", tofCorr, &b_tofCorr);
   fChain->SetBranchAddress("beta", beta, &b_beta);
   //fChain->SetBranchAddress("T0", &T0, &b_T0);
   //fChain->SetBranchAddress("vzVpd", &vzVpd, &b_vzVpd);
   //fChain->SetBranchAddress("Iwest", &Iwest, &b_Iwest);
   //fChain->SetBranchAddress("Ieast", &Ieast, &b_Ieast);

   Notify();
}

Bool_t TOFrPicoDst::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normaly not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void TOFrPicoDst::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t TOFrPicoDst::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef TOFrPicoDst_cxx
