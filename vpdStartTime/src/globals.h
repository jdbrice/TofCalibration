

#include <stdlib.h>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <vector>

#include <sys/types.h>
#include "math.h"
#include "string.h"

#ifndef __CINT__
#include "TROOT.h"
#include "Riostream.h"
#include "TFile.h"
#include "TChain.h"
#include "TH1.h"
#include "TH2.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TProfile.h"
#include "TTree.h"
#include "TNtuple.h"  
#include "TRandom.h"
#include "TF1.h"
#include "TMath.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TPostScript.h"
#include "TString.h"
 #include "TLeaf.h"

using namespace std;
using std::cout;   
using std::endl;  

#endif



const Int_t maxTOFrHits = 8000;
const Int_t minpVPDHits = 3;
const Double_t piMass = 0.13957;
const Double_t pMass = 0.938272;
const Double_t kMass = 0.493677;
const Double_t c_light = 29.9792458; 
const Float_t phaseDiff = 0;

const Int_t nPVPDChannel = 38; 
const Int_t nTray = 1; 
const Int_t nModule = 4; 

const Int_t nCell = 6; 
const Int_t nBoard = 1; 
const Int_t ntCell = 24;

const Double_t minTOT_pVPD = 10;
const Double_t maxTOT_pVPD = 40;


const Double_t minTOT_TOFr = 5; 
const Double_t maxTOT_TOFr = 40;
const Double_t minZ_TOFr = -3.05; 
const Double_t maxZ_TOFr = 3.05; 
const Double_t rCut = 1.; 
const Double_t vzCut = 6.; 
const Int_t nHitsFitCut = 25;
const Double_t massCut = 0.4; 
const Double_t tdiffCut = 0.8;

const Int_t nIter = 5; 
const Int_t maxpvpdNBin = 50; 
const Int_t maxNBin = 40; 
const Int_t maxNZBin = 40; 
const Bool_t fixBinWidth = kFALSE; 
const Bool_t fixNBin = !(fixBinWidth); 
const Double_t  mrmscut = 0.3;

const Int_t ntdcbin = 800;
const Double_t tlimit = 40.; 
const Int_t ntoftdcbin = 100;
const Double_t toftlimit = 2.;
const Double_t smallvalue = 0.0005;
const Double_t bigtime = 200.0; 
const Int_t maxiteration = 10; 
const Double_t pCut = 0.20; 



const Char_t *pvpdOutput = "pvpdCali";
const Char_t *t0Output = "t0Cali";
const Char_t *totOutput = "totCali";
const Char_t *zOutput = "zCali";
const Char_t *delayOutput = "delay";
const Char_t *resoOutput = "reso";
Int_t trayStart=0;
Int_t TDIGStart=0;
Double_t dedxcut = 2.8;


TF1 *piBetaCut;

double vpdDelay[38] = {11.852, 11.928, 11.764, 11.899, 11.915, 23.927, 11.735, 11.925, -0.014, 0.038,
  11.920, 0.143, -0.075, 24.103, 24.046, 0.0619, 0.0657, 0.0325, -0.026,         
  11.915, 10.695, 11.913, 11.915, 11.814, 24.026, 11.934, 11.800, 0.067, -0.500,
  11.876, 0.000, 0.0220, 24.016, 24.284, 0.046, 0.055, 0.019, -0.030             
}; 

const Int_t tsfM = 5; 


Int_t nBinPVPD[nPVPDChannel];
Double_t tofrTOTBins[nTray][nBoard][ntCell][maxNBin+1];
Double_t tofrTOTX[nTray][nBoard][ntCell][maxNBin];    
Double_t tofrZBins[nTray][nBoard][ntCell][maxNZBin+1];
Double_t pvpdTOTBins[nPVPDChannel][maxpvpdNBin+1];
Double_t pvpdTOTX[nPVPDChannel][maxpvpdNBin];   
Double_t pvpdTOTcorr[nPVPDChannel][maxpvpdNBin];  
Double_t tofrDelay[nTray][nBoard][nCell];
Double_t tofrRes[nTray][nBoard][nCell];

Double_t tofrTOT[nTray][nBoard][ntCell][maxNBin];
Double_t tofrZ[nTray][nBoard][ntCell][maxNZBin];
Double_t pvpdTOT[nPVPDChannel][maxpvpdNBin];

Double_t tofrTOTbad[nTray][nBoard][ntCell];
Double_t tofrZbad[nTray][nBoard][ntCell];
Double_t tofrTOTbadcorr[nTray][nBoard][ntCell][maxNBin];
Double_t tofrZbadcorr[nTray][nBoard][ntCell][maxNZBin];


vector<double> tofrtots[nTray][nBoard][ntCell];
vector<double> tofrzs[nTray][nBoard][ntCell];
vector<double> pvpdtots[nPVPDChannel];

TSplineFit *tsfPVPD[nPVPDChannel]; 
TSplineFit *tsfTOT[nTray][nBoard][ntCell]; 
TSplineFit *tsfZ[nTray][nBoard][ntCell]; 

TH1D *hPVPDMeanvsIt[nPVPDChannel]; 
TH1D *hPVPDFitMeanvsIt[nPVPDChannel]; 
TH1D *hPVPDFitSigmavsIt[nPVPDChannel]; 
TH1D *hPVPDResovsIt[nPVPDChannel]; 

TH2D *hPVPD[nPVPDChannel];
TH2D *hPVPDCorr[nPVPDChannel]; 
TH2D *hTOT[nTray][nBoard][ntCell];
TH2D *hZ[nTray][nBoard][ntCell];
TH1D *hpPVPD[nPVPDChannel];
TH1D *htPVPD1st[nPVPDChannel];
TH1D *htPVPD[nPVPDChannel];
TH1D *hpTOT[nTray][nBoard][ntCell];
TH1D *htTOT[nTray][nBoard][ntCell];
TH1D *hpZ[nTray][nBoard][ntCell];
TH1D *hDelay[nTray][nBoard][nCell];
TH1D *hDelay2[nTray][nBoard][nCell];
TH1D *htotPVPD[nPVPDChannel];


TH2D *hVzCorr;
TH2D *hVzCorr2;


TH2D *hVzCorrvzdiff;
TH2D *hVzCorr2vzdiff;
TH2D *hvxvvy;
TH2D *hvxvvyvzdiff;

TH1D *hResoPVPD11;
TH1D *hResoPVPD22;
TH1D *hResoPVPD33;
 
TH2D *hResoTOT;
TH2D *hResoZ;
TH2D *hResoD;
TH1D *hResoTOT1D;
TH1D *hResoZ1D;
TH1D *hResoD1D;

TH2D *hResoTOT_Pos;
TH2D *hResoZ_Pos;
TH2D *hResoD_Pos;
TH1D *hResoTOT1D_Pos;
TH1D *hResoZ1D_Pos;
TH1D *hResoD1D_Pos;

TH2D *hResoTOT_Neg;
TH2D *hResoZ_Neg;
TH2D *hResoD_Neg;
TH1D *hResoTOT1D_Neg;
TH1D *hResoZ1D_Neg;
TH1D *hResoD1D_Neg;

TH2D *hBetaP_TOT = new TH2D("hBetaP_TOT","1/beta vs p after slewing correction", 1000, 0, 5, 1000, 0, 5);
TH2D *hMassP_TOT = new TH2D("hMassP_TOT","M^2 vs p after slewing correction",1000, 0, 5, 2000, -1, 4);
TH2D *hBetaP_Z = new TH2D("hBetaP_Z","1/beta vs p after Z correction", 1000, 0, 5, 1000, 0, 5);
TH2D *hMassP_Z = new TH2D("hMassP_Z","M^2 vs p after Z correction",1000, 0, 5, 2000, -1, 4);
TH2D *hBetaP_D = new TH2D("hBetaP_D","1/beta vs p after T0 correction", 1000, 0, 5, 1000, 0, 5);
TH2D *hMassP_D = new TH2D("hMassP_D","M^2 vs p after T0 correction",1000, 0, 5, 2000, -1, 4);

TH2D *hBetaP_TOT_Pos = new TH2D("hBetaP_TOT_Pos","1/beta vs p after slewing correction positive", 1000, 0, 5, 1000, 0, 5); 
TH2D *hMassP_TOT_Pos = new TH2D("hMassP_TOT_Pos","M^2 vs p after slewing correction positive",1000, 0, 5, 2000, -1, 4); 
TH2D *hBetaP_Z_Pos = new TH2D("hBetaP_Z_Pos","1/beta vs p after Z correction positive", 1000, 0, 5, 1000, 0, 5);
TH2D *hMassP_Z_Pos = new TH2D("hMassP_Z_Pos","M^2 vs p after Z correction positive",1000, 0, 5, 2000, -1, 4);
TH2D *hBetaP_D_Pos = new TH2D("hBetaP_D_Pos","1/beta vs p after T0 correction positive", 1000, 0, 5, 1000, 0, 5);
TH2D *hMassP_D_Pos = new TH2D("hMassP_D_Pos","M^2 vs p after T0 correction positive",1000, 0, 5, 2000, -1, 4);

TH2D *hBetaP_TOT_Neg = new TH2D("hBetaP_TOT_Neg","1/beta vs p after slewing correction negative", 1000, 0, 5, 1000, 0, 5);
TH2D *hMassP_TOT_Neg = new TH2D("hMassP_TOT_Neg","M^2 vs p after slewing correction negative",1000, 0, 5, 2000, -1, 4);
TH2D *hBetaP_Z_Neg = new TH2D("hBetaP_Z_Neg","1/beta vs p after Z correction negative", 1000, 0, 5, 1000, 0, 5);
TH2D *hMassP_Z_Neg = new TH2D("hMassP_Z_Neg","M^2 vs p after Z correction negative",1000, 0, 5, 2000, -1, 4); 
TH2D *hBetaP_D_Neg = new TH2D("hBetaP_D_Neg","1/beta vs p after T0 correction negative", 1000, 0, 5, 1000, 0, 5);
TH2D *hMassP_D_Neg = new TH2D("hMassP_D_Neg","M^2 vs p after T0 correction negative",1000, 0, 5, 2000, -1, 4);
TH1D *hnull;


TH1D *hplus = new TH1D("hplus","dt positive", 200,-2.5,2.5);
TH1D *hminus = new TH1D("hminus","dt negtive", 200, -2.5, 2.5);

TH1D *hdeltatwest = new TH1D("hdeltatwest", "delta time between two channels: west", 20000, -100, 100);
TH1D *hdeltateast = new TH1D("hdeltateast", "delta time between two channels: east", 20000, -100, 100);


TFile *forwrite;

Int_t doTOT(TChain* chain, Int_t iter); 
Int_t doZ(TChain* chain, Int_t iter); 
Int_t doT0(TChain* chain); 
Int_t doStart(TChain* chain); 
Int_t doDelay(TChain* chain,Int_t iter); 

Int_t binning(TChain* chain, Int_t choice); 
Int_t initPVPD(); 
Int_t initTOFr(); 
Int_t getH(Int_t choice); 
Int_t getV(TH1D* histo); 

Int_t outputPVPD(TChain* chain); 
Int_t outputparPVPD(TChain* chain); 
Int_t outputDelay(TChain* chain, Int_t iter); 
Int_t outputTOT(TChain* chain, Int_t iter); 
Int_t outputZ(TChain* chain, Int_t iter); 
Int_t checkReso(TChain* chain, Int_t iter); 
float linearInter(float x1, float x2, float y1, float y2, float x);   
Bool_t xpisnan(Double_t x);

static bool FIRST_RUN = true;
