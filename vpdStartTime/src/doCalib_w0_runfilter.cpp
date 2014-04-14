/********************************************************
 *  MRPC TOF Calibration by using the TSplineFit class  *
 ********************************************************/

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

using namespace std;
using std::cout;   
using std::endl;  

#endif

#include "TSplineFit.h"
//#include "PVPDPicoDst.h"
#include "TOFrPicoDst.h"
#include "StartPicoDst.h"
#include "TOTPicoDst.h"
#include "ZPicoDst.h"
#include "DelayPicoDst.h"

#include "cuts.h"

/*
#include "T0PicoDst.h"
*/

//-----------------------------------------------------------------//
// Define common parameters/histograms and functions
const Int_t maxTOFrHits = 8000;
const Int_t minpVPDHits = 3;
const Double_t piMass = 0.13957;
const Double_t pMass = 0.938272;
const Double_t kMass = 0.493677;
const Double_t c_light = 29.9792458; // light speed in vacuum, unit: cm/ns
const Float_t phaseDiff = 0;//no phase diff in run 9; run 8 -33796;//47.8; // ns

const Int_t nPVPDChannel = 38; // 19 (east) + 19 (west)
const Int_t nTray = 1; // number of trays installed
const Int_t nModule = 4; // number of MRPC modules per tray
//const Int_t nTDIG =1; //number of TDIGs 
const Int_t nCell = 6; // number of cells per TDIG //was 6 per module 
const Int_t nBoard = 1; //number of TDIG board per tray
const Int_t ntCell = 24;//number of cells per TDIG...

const Double_t minTOT_pVPD = 10;//5;//3; // ns
const Double_t maxTOT_pVPD = 40;//31; // ns
//const Double_t minTOT_TOFr = 12.; // ns   zebo
//const Double_t maxTOT_TOFr = 18.; // ns   zebo
const Double_t minTOT_TOFr = 5; // ns
const Double_t maxTOT_TOFr = 40;//31;// 22.; // ns
const Double_t minZ_TOFr = -3.05; // cm
const Double_t maxZ_TOFr = 3.05; // cm
const Double_t rCut = 1.; //cm
const Double_t vzCut = 6.; //cm
const Int_t nHitsFitCut = 25;
const Double_t massCut = 0.4; //GeV
const Double_t tdiffCut = 0.8;//ns, reject hits from other collisions

const Int_t nIter = 5; // iteration times for TOFr correction
const Int_t maxpvpdNBin = 50; // most binning number in spline fit
const Int_t maxNBin = 40; // most binning number in spline fit
const Int_t maxNZBin = 40; // most binning number in TZ spline fit
const Bool_t fixBinWidth = kFALSE; // binning method with fixed bin width
const Bool_t fixNBin = !(fixBinWidth); // binning method with fixed bin number
const Double_t  mrmscut = 0.3;

const Int_t ntdcbin = 800;//400; 
const Double_t tlimit = 40.; //ns, TDC window
const Int_t ntoftdcbin = 100;
const Double_t toftlimit = 2.;
const Double_t smallvalue = 0.0005;//0.001; //ns
const Double_t bigtime = 200.0; // ns
const Int_t maxiteration = 10; // pvpd slewing repeat times
const Double_t pCut = 0.20; // % of the highest (tdc-cor) values to be rejected - rafael



const Char_t *pvpdOutput = "pvpdCali";
const Char_t *t0Output = "t0Cali";
const Char_t *totOutput = "totCali";
const Char_t *zOutput = "zCali";
const Char_t *delayOutput = "delay";
const Char_t *resoOutput = "reso";
Int_t trayStart=0;
Int_t TDIGStart=0;
Double_t dedxcut = 2.8;


TF1 *piBetaCut;//funtion for pion beta cut 

double vpdDelay[38] = {11.852, 11.928, 11.764, 11.899, 11.915, 23.927, 11.735, 11.925, -0.014, 0.038,
  11.920, 0.143, -0.075, 24.103, 24.046, 0.0619, 0.0657, 0.0325, -0.026,         // West vpd
  11.915, 10.695, 11.913, 11.915, 11.814, 24.026, 11.934, 11.800, 0.067, -0.500,
  11.876, 0.000, 0.0220, 24.016, 24.284, 0.046, 0.055, 0.019, -0.030             // east vpd
}; 

const Int_t tsfM = 5; // points used to fit one spline

//Int_t nBinTOFr[nTray][nBoard][nCell];
Int_t nBinPVPD[nPVPDChannel];
Double_t tofrTOTBins[nTray][nBoard][ntCell][maxNBin+1];
Double_t tofrTOTX[nTray][nBoard][ntCell][maxNBin];    //Tofr Tot center values
Double_t tofrZBins[nTray][nBoard][ntCell][maxNZBin+1];
Double_t pvpdTOTBins[nPVPDChannel][maxpvpdNBin+1];
Double_t pvpdTOTX[nPVPDChannel][maxpvpdNBin];   // vpd Tot center values 
Double_t pvpdTOTcorr[nPVPDChannel][maxpvpdNBin];  
Double_t tofrDelay[nTray][nBoard][nCell];
Double_t tofrRes[nTray][nBoard][nCell];//CONTINUE HERE

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

TH1D *hPVPDMeanvsIt[nPVPDChannel]; // rafael
TH1D *hPVPDFitMeanvsIt[nPVPDChannel]; // rafael
TH1D *hPVPDFitSigmavsIt[nPVPDChannel]; // rafael
TH1D *hPVPDResovsIt[nPVPDChannel]; // rafael

TH2D *hPVPD[nPVPDChannel];
TH2D *hPVPDCorr[nPVPDChannel]; // rafael
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
//TH2D *hVzCorrexp;
//TH2D *hVzCorr2exp;

TH2D *hVzCorrvzdiff;
TH2D *hVzCorr2vzdiff;
//TH2D *hVzCorrexpvzdiff;
//TH2D *hVzCorr2expvzdiff;
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

//check distortion
TH1D *hplus = new TH1D("hplus","dt positive", 200,-2.5,2.5);
TH1D *hminus = new TH1D("hminus","dt negtive", 200, -2.5, 2.5);

TH1D *hdeltatwest = new TH1D("hdeltatwest", "delta time between two channels: west", 20000, -100, 100);
TH1D *hdeltateast = new TH1D("hdeltateast", "delta time between two channels: east", 20000, -100, 100);


TFile *forwrite;

Int_t doTOT(TChain* chain, Int_t iter); // TOT slewing correction
Int_t doZ(TChain* chain, Int_t iter); // Correction on local Z hit position on MRPC cell
Int_t doT0(TChain* chain); // event start time without start detector
Int_t doStart(TChain* chain); //event start time from start detector
Int_t doDelay(TChain* chain,Int_t iter); // signal delay from TOFr pads to electronic

Int_t binning(TChain* chain, Int_t choice); // histogramming with proper bin
Int_t initPVPD(); // initialize necessary parameters for pVPD
Int_t initTOFr(); // initialize necessary parameters for TOFr
Int_t getH(Int_t choice); // fill histogram for calibration according to choice
Int_t getV(TH1D* histo); // get vectors for spline-fitting, from histogram h

Int_t outputPVPD(TChain* chain); // output pVPD calibration results with TSpline fit
Int_t outputparPVPD(TChain* chain); // output pVPD calibration results with linear interpolation
Int_t outputDelay(TChain* chain, Int_t iter); // output TOFr signal delay results
Int_t outputTOT(TChain* chain, Int_t iter); // output TOFr TOT calibration results
Int_t outputZ(TChain* chain, Int_t iter); // output TOFr Z calibration results
Int_t checkReso(TChain* chain, Int_t iter); // check pVPD and TOFr time resolution
float linearInter(float x1, float x2, float y1, float y2, float x);   //linearInter function
Bool_t xpisnan(Double_t x);


//-------------------------------------------------------------------------//

Int_t main(int argc, char **argv)
{
  if(argc!=4 && argc!=1) return 0;
  
  const char *FileInput = "test.list";
//  Char_t *FileOutput = "test.root";
 TString *string;
 TString *string2;
 
  if(argc==1){
    string = new TString("0");
  }
 
  if(argc==4){
    FileInput = argv[1];
    //FileOutput = argv[2];
    string = new TString(argv[2]);
    string2 = new TString(argv[3]);    
  }
  trayStart = string->Atoi();
  TDIGStart = string2->Atoi();

  if(argc!=1 && argc!=4) return -1;

  //----------------------------------
  // Open files and add to the chain
  //---------------------------------- 
  Int_t fileNumber = 0;
  Char_t FileList[512];

//  TChain *pvpdchain = new TChain("pvpd");
  TChain *tofrchain = new TChain("tof");
 
//        TFile *ftmp = new TFile(FileInput);
//        printf(" read in file %s\n",FileInput);
//        tofrchain->Add(FileInput);
//        delete ftmp;

  ifstream* inputStream = new ifstream;
  inputStream->open(FileInput);
  if (!(inputStream)) {
    printf("can not open list file\n");
    return 0;
  }
  for(;inputStream->good();) {
    inputStream->getline(FileList,512);
    if( inputStream->good() ) {
      TFile *ftmp = new TFile(FileList);
      if(!ftmp||!(ftmp->IsOpen())||!(ftmp->GetNkeys())) {
        printf(" file %s error in opening!!!\n",FileList);
      } else {
        printf(" read in file %s\n",FileList);
        tofrchain->Add(FileList);
        fileNumber++;
      }
      delete ftmp;
    }
  }
  
  printf(" %d files read in\n",fileNumber);
  cout<<"do calibration for tray "<<trayStart*nTray+1<<"-"<<(trayStart+1)*nTray<<" TDIG "<<TDIGStart*nBoard<<endl;
  //---------------------------------

  //---------------------------------------
  // main work

  if(fixNBin) {
    binning(tofrchain, 0); // also for for outputparpvpd  //for outputparpvpd
  }
  initPVPD();

  doStart(tofrchain);  // get event start time (T0) from pVPD
  outputPVPD(tofrchain); // write out start tree for TOFr
//   outputparPVPD(tofrchain); // write out start tree for TOFr //for outputparpvpd
//   return 1;  //for outputparpvpd
//   doDelay(tofrchain,0); 
   // doT0(tofrchain); // get event start time without pVPD information
   // return 0; 
  
  if(fixNBin) { 
//    binning(tofrchain, 1);
//    binning(tofrchain, 2);
  }
  
    initTOFr();
 
  for(Int_t iter=1; iter<=nIter; iter++) {
  	
  //  doTOT(tofrchain, iter);
  //  doZ(tofrchain, iter);
  //  doDelay(tofrchain, iter);
  //  checkReso(tofrchain, iter);
    
  	}
 
  cout<<"done!!!"<<endl;
  return 1;
}

Int_t binning(TChain* chain, Int_t choice)
{
  Int_t nevents = (int)chain->GetEntries();

  if(choice==0) { // for pVPD TOT bin
    TOFrPicoDst *t1 = new TOFrPicoDst(chain);
    cout<< " woking on binning! choice = "<<choice<<endl;
 
    for(Int_t i=0;i<nevents;i++) {
      t1->GetEntry(i);
     
      if(i%(nevents/100)==0) cout<< i <<" events processed!"<<endl;     
      Int_t ieast = t1->numberOfVpdEast;
      Int_t iwest = t1->numberOfVpdWest;
      
      if(ieast<minpVPDHits&&iwest<minpVPDHits) continue;
     
      if(iwest>=minpVPDHits){
        for(Int_t j=0; j<nPVPDChannel/2; j++) {
		
//		if(j==7) continue;//MASK		

	  Double_t tot = t1->vpdTotWest[j];
	  if(tot>minTOT_pVPD && tot<maxTOT_pVPD) pvpdtots[j].push_back(tot);
        }
      }
      if(ieast>=minpVPDHits){
        for(Int_t j=nPVPDChannel/2; j<nPVPDChannel; j++) {
          Double_t tot = t1->vpdTotEast[j-nPVPDChannel/2];
          if(tot>minTOT_pVPD && tot<maxTOT_pVPD) pvpdtots[j].push_back(tot);
        }
      }
    }
    
    for(Int_t i=0; i<nPVPDChannel; i++) {
      Int_t size = pvpdtots[i].size();
      cout<<"channel "<<i<<": "<< size << " hits"<<endl;
      if(size<maxpvpdNBin) {
	Double_t step = (maxTOT_pVPD-minTOT_pVPD)/maxpvpdNBin;
	for(Int_t j=0; j<=maxpvpdNBin; j++) {
          pvpdTOTBins[i][j] = step*j+minTOT_pVPD; // dummy bin bounds//added mintotpvpd so when the parameters are read-in, they do not accept bad hits in dead channels for TOT range.
//	  pvpdTOTBins[i][j] = step*j+10; // dummy bin bounds//added +10 so when the parameters are read-in, they do not accept bad hits in dead channels for TOT range.
//          pvpdTOTBins[i][j] = step*j; // dummy bin bounds

	}
	cout << " pVPD Channel " << i << " seems dead ! " << "( " << size << " hits)" <<endl;
      } else {
	Int_t step = size/maxpvpdNBin; // (int)size/maxNBin is less than (float)size/maxNBin -> more entries in last bin
	std::sort(pvpdtots[i].begin(), pvpdtots[i].end());
	pvpdTOTBins[i][0] = pvpdtots[i].at(0);
	pvpdTOTBins[i][maxpvpdNBin] = pvpdtots[i].at(size-1);
	Double_t tmp;
	for(Int_t j=1; j<maxpvpdNBin; j++) {
	  pvpdTOTBins[i][j] = (pvpdtots[i].at(step*j)+pvpdtots[i].at(step*j-1))/2.0;
	  tmp = 0.;
	  for(Int_t k=step*(j-1); k<step*j; k++) {
	    tmp += pvpdtots[i].at(k);
	  }
	  pvpdTOT[i][j-1] = tmp/step;
	}
	tmp = 0;
	for(Int_t k=step*(maxpvpdNBin-1); k<size; k++) {
	  tmp += pvpdtots[i].at(k);
	}
	pvpdTOT[i][maxpvpdNBin-1] = tmp/(size-step*(maxpvpdNBin-1));
      }
    } // pvpd bins filled over!
    
  } else { // for TOFr TOT and Z bin
    TOFrPicoDst *t2 = new TOFrPicoDst(chain);

    char rootfile[500];
   // sprintf(rootfile, "%s.root", pvpdOutput);
    /*TChain *startchain = new TChain("StartPicoDst");
    startchain->AddFile(rootfile);

    StartPicoDst *start = new StartPicoDst(startchain);
    
    if(nevents!=startchain->GetEntries()) {
      cout << " Inconsistent TOFr and Start entry ! " << endl;
      return -1;
    }
*/
    for(Int_t i=0;i<nevents;i++) {
      t2->GetEntry(i);
  //    start->GetEntry(i);

      Double_t t0 = t2->T0;
      Int_t hits = t2->nTofHits;
      if(hits<=0) continue;

      for(Int_t j=0; j<hits; j++) {
	Int_t tray = t2->tray[j]-1-nTray*trayStart;
	Int_t module = (t2->module[j]-1)%4;
//	Int_t tmodule = module%4;//-nModule*moduleStart;
	Int_t cell = t2->cell[j]-1;
	Int_t tcell = module*6+cell;//watch for zombies

	if(tray>=nTray || tray<0) continue;
	if(module>=nModule || module<0) continue;
	if(cell>=nCell || cell<0) continue;	
	
	Double_t p = t2->pt[j]*TMath::CosH(t2->eta[j]);
	Double_t dedx = t2->dedx[j];
	if(p<0.3 || p>0.6 || dedx > dedxcut) continue;
        if(sqrt(pow(t2->vertexX,2)+pow(t2->vertexY,2))>rCut) continue;    //for auau
//        if(sqrt(pow(t2->dcaX[j],2)+pow(t2->dcaY[j],2))>rCut) continue; //for pprun9
         if(TMath::Abs(t2->vertexZ-t2->vzVpd)>vzCut) continue; //for auau
//         if(TMath::Abs(t2->dcaZ[j]-t2->vzVpd)>vzCut) continue; //for pprun9

	Double_t tracklength = t2->length[j];
	Double_t betagamma = p/piMass;
	Double_t velocity = TMath::Sqrt(1.0/(1.0/betagamma/betagamma+1.0))*c_light;
	Double_t tofPi = tracklength/velocity;

	Double_t tdc = t2->leTime[j];
        tdc -= phaseDiff;
        while(tdc>51200) tdc-=51200;
	Double_t dt = tdc-t0-tofPi-tofrDelay[tray][0][tcell];//should be TDIG instead of nTDIG, but its 1 and lazy.
//	cout<<"binning dt, tdc, t0, tofPi, tofrDelay: "<<dt<<", "<<tdc<<", "<<t0<<", "<<tofPi<<", "<<tofrDelay[tray][module][cell]<<endl;
	if(TMath::Abs(dt)>toftlimit) continue;
 
	if(choice==1) { // TOFr TOT 
	  Double_t tot = t2->tot[j];
//	cout<<"choice1, tot: "<<tot<<endl;
	  if(tot>minTOT_TOFr && tot<maxTOT_TOFr) tofrtots[tray][0][tcell].push_back(tot);
	  //if(tot>totlolim[tray][module][cell] && tot<totuplim[tray][module][cell]) tofrtots[tray][module][cell].push_back(tot);
	}
	if(choice==2) { // TOFr Z 
	  Double_t z = t2->zLocal[j];
//	cout<<"choice2, z: "<<z<<endl;
	  if(z>minZ_TOFr && z<maxZ_TOFr) tofrzs[tray][0][tcell].push_back(z);
	}
      }
    }
    
    for(Int_t i=0; i<nTray; i++) {
      for(Int_t j=0; j<nBoard; j++) {
	for(Int_t k=0; k<ntCell; k++) {
	  if(choice==1) {
	    Int_t size = tofrtots[i][j][k].size();
            cout << "TOFr "<<trayStart*nTray<<"th Tray "<<TDIGStart*nBoard<<"th Board "<<k<<"th Cell ( "<<size<< " hits)"<<endl;
	    if(size<maxNBin) {
	      Double_t step = (maxTOT_TOFr-minTOT_TOFr)/maxNBin;
	      for(Int_t l=0; l<=maxNBin; l++) {
		tofrTOTBins[i][j][k][l] = step*l; // dummy bin bounds
	      }
	      cout << " TOFr Channel at Tray " << trayStart*nTray << " Board " << TDIGStart*nBoard <<" Cell "<<k<< " seems dead ! (" << size << " hits)"<<endl;
	    } else {
	      Int_t step = (int)(size*1./maxNBin);
	      Int_t binl = (int)(size*0./maxNBin);
	      Int_t binh = size-binl;
	      std::sort(tofrtots[i][j][k].begin(), tofrtots[i][j][k].end());
	      tofrTOTBins[i][j][k][0] = tofrtots[i][j][k].at(binl);
	      tofrTOTBins[i][j][k][maxNBin] = tofrtots[i][j][k].at(binh-1);
	      Double_t tmp;
	      for(Int_t l=1; l<maxNBin; l++) {
		tofrTOTBins[i][j][k][l] = (tofrtots[i][j][k].at(binl+step*l)+tofrtots[i][j][k].at(binl+step*l-1))/2.0;
		tmp = 0.;
		for(Int_t m=binl+step*(l-1); m<binl+step*l; m++) {
		  tmp += tofrtots[i][j][k].at(m);
		}
		tofrTOT[i][j][k][l-1] = tmp/step;
	      }
	      tmp = 0;
	      for(Int_t m=binl+step*(maxNBin-1); m<binh; m++) {
		tmp += tofrtots[i][j][k].at(m);
	      }
	      tofrTOT[i][j][k][maxNBin-1] = tmp/(size-2*binl-step*(maxNBin-1));
	    }
	  }
	  if(choice==2) {
	    Int_t size = tofrzs[i][j][k].size();
	    if(size<maxNZBin) {
	      Double_t step = (maxZ_TOFr-minZ_TOFr)/maxNZBin;
	      for(Int_t l=0; l<=maxNZBin; l++) {
		tofrZBins[i][j][k][l] = step*l; // dummy bin bounds
	      }
	      cout << " TOFr Channel at Tray " << i << " Board " << TDIGStart*nBoard << " Cell "<<k<<" seems dead ! " << endl;
	    } else {
	      Int_t step = size/maxNZBin;
	      std::sort(tofrzs[i][j][k].begin(), tofrzs[i][j][k].end());
	      tofrZBins[i][j][k][0] = tofrzs[i][j][k].at(0);
	      tofrZBins[i][j][k][maxNZBin] = tofrzs[i][j][k].at(size-1);
	      Double_t tmp;
	      for(Int_t l=1; l<maxNZBin; l++) {
		tofrZBins[i][j][k][l] = (tofrzs[i][j][k].at(step*l)+tofrzs[i][j][k].at(step*l-1))/2.0;
		tmp = 0.;
		for(Int_t m=step*(l-1); m<step*l; m++) {
		  tmp += tofrzs[i][j][k].at(m);
		}
		tofrZ[i][j][k][l-1] = tmp/step;
	      }
	      tmp = 0;
	      for(Int_t m=step*(maxNZBin-1); m<size; m++) {
		tmp += tofrzs[i][j][k].at(m);
	      }
	      tofrZ[i][j][k][maxNZBin-1] = tmp/(size-step*(maxNZBin-1));
	    }
	  }
	}//cell loop
      }//module loop
    } // end nTray loop, tofr bins filled over !
  } // end choice

  // clear vectors
  if(choice==0) {
    for(Int_t i=0; i<nPVPDChannel; i++) {
      pvpdtots[i].clear();
    }
  } else {
    for(Int_t i=0; i<nTray; i++) {      ///BUTTER LOOK HERE LOOPS OVER TDIG 
      for(Int_t j=0; j<nBoard; j++) {
	for(Int_t k=0; k<ntCell; k++) {
	  tofrtots[i][j][k].clear();
	  tofrzs[i][j][k].clear();	
	}//cell
      }//module
    }  //tray
  }
  cout<<"TOFr tot	"<<endl;
for(int tmpk=0; tmpk<maxNBin; tmpk++)
	cout<<setw(10)<<tofrTOT[0][0][0][tmpk];
  return 1;
}

Int_t initPVPD()
{
  // Initialize pVPD

  Char_t buf[100];
  for(Int_t i=0; i<nPVPDChannel; i++) {
    sprintf(buf, "pVPD_Channel_%d", i);
    if(fixNBin) {
      hPVPD[i] = new TH2D(buf, buf, maxpvpdNBin, pvpdTOTBins[i], ntdcbin, -tlimit, tlimit);
    } else if (fixBinWidth) {
      hPVPD[i] = new TH2D(buf, buf, maxpvpdNBin, minTOT_pVPD,  maxTOT_pVPD, ntdcbin, -tlimit, tlimit);
    } else {
      cout << " ? You have not decided which binning algrithm to be used !?" << endl;
      return -1;
    }  

	sprintf(buf, "upVPD_Channel_%d_TOT", i);
    htotPVPD[i] = new TH1D(buf, buf, maxpvpdNBin*10, minTOT_pVPD-2,  30); 	
    sprintf(buf, "pVPD1st_Channel_Time_%d", i);
    htPVPD1st[i] = new TH1D(buf, buf, ntdcbin, -1.5*tlimit, 1.5*tlimit); 
    sprintf(buf, "pVPD_Channel_Time_%d", i);
    htPVPD[i] = new TH1D(buf, buf, ntdcbin, -tlimit/20, tlimit/20);
  

    // begin - rafael
    sprintf(buf, "pVPDCorr_Channel_%d", i);
    if(fixNBin) {
      hPVPDCorr[i] = new TH2D(buf, buf, maxpvpdNBin, pvpdTOTBins[i], ntdcbin, -tlimit, tlimit);
    } else if (fixBinWidth) {
      hPVPDCorr[i] = new TH2D(buf, buf, maxpvpdNBin, minTOT_pVPD,  maxTOT_pVPD, ntdcbin, -tlimit, tlimit);
    } else {
      cout << " ? You have not decided which binning algorithm to be used !?" << endl;
      return -1;
    }
    sprintf(buf, "pVPDMeanvsIt_%d",i);
    hPVPDMeanvsIt[i] = new TH1D(buf,buf,maxiteration,0.5,maxiteration+0.5);
    sprintf(buf, "pVPDFitMeanvsIt_%d",i);
    hPVPDFitMeanvsIt[i] = new TH1D(buf,buf,maxiteration,0.5,maxiteration+0.5);
    sprintf(buf, "pVPDFitSigmavsIt_%d",i);
    hPVPDFitSigmavsIt[i] = new TH1D(buf,buf,maxiteration,0.5,maxiteration+0.5);
    sprintf(buf, "pVPDResovsIt_%d",i);
    hPVPDResovsIt[i] = new TH1D(buf,buf,maxiteration,0.5,maxiteration+0.5);
    // end - rafael
  }


  sprintf(buf, "pVPD_Res_11");
  hResoPVPD11 = new TH1D(buf, buf, 400, -20, 20);

  sprintf(buf, "pVPD_Res_22");
  hResoPVPD22 = new TH1D(buf, buf, 400, -20, 20);

  sprintf(buf, "pVPD_Res_33");
  hResoPVPD33 = new TH1D(buf, buf, 400, -20, 20);

  for(Int_t i=0; i<nTray; i++) {
    for(Int_t j=0; j<nBoard; j++) {
      for(Int_t k=0; k<nCell; k++) {
        sprintf(buf, "Delay_Tray_%d_Board_%d_Cell_%d", i+nTray*trayStart+1, j+(TDIGStart)*nBoard, k);
        hDelay[i][j][k] = new TH1D(buf, buf, 1600, -1.5*bigtime, 0.5*bigtime);
        sprintf(buf, "Delay2_Tray_%d_Board_%d_Cell_%d", i+nTray*trayStart+1, j+(TDIGStart)*nBoard, k);
        hDelay2[i][j][k] = new TH1D(buf, buf, ntoftdcbin, -2*toftlimit, 2*toftlimit);
      }
    }
  }

  sprintf(buf, "VzCorr");
  hVzCorr = new TH2D(buf, buf, 400, -200, 200, 400, -200, 200);

  sprintf(buf, "VzCorr2");
  hVzCorr2 = new TH2D(buf, buf, 400, -200, 200, 500, -50, 50);

//  sprintf(buf, "VzCorrexp");
//  hVzCorrexp = new TH2D(buf, buf, 1000, -1000, 1000, 1000, -1000,1000);

//  sprintf(buf, "VzCorr2exp");
//  hVzCorr2exp = new TH2D(buf, buf, 1000, -1000, 1000, 1000, -500, 500);
////////////////////////////////////
  sprintf(buf, "VzCorrvzdiff");
  hVzCorrvzdiff = new TH2D(buf, buf, 400, -200, 200, 400, -200, 200);

  sprintf(buf, "VzCorr2vzdiff");
  hVzCorr2vzdiff = new TH2D(buf, buf, 400, -200, 200, 500, -50, 50);

//  sprintf(buf, "VzCorrexpvzdiff");
//  hVzCorrexpvzdiff = new TH2D(buf, buf, 1000, -1000, 1000, 1000, -1000,1000);

//  sprintf(buf, "VzCorr2expvzdiff");
//  hVzCorr2expvzdiff = new TH2D(buf, buf, 1000, -1000, 1000, 1000, -500, 500);

  sprintf(buf, "vxvvy");
  hvxvvy = new TH2D(buf, buf, 500, -10, 10, 500, -10, 10);
  sprintf(buf, "vxvvyvzdiff");
  hvxvvyvzdiff = new TH2D(buf, buf, 500, -10, 10, 500, -10, 10);


/*
  cout<<"JOEY REPLACE THE FOLLOWING WITH THE ABOVE"<<endl;
  sprintf(buf, "VzCorr");
  hVzCorr = new TH2D(buf, buf, 2000, -15000, 15000, 2000, -15000, 15000);
  sprintf(buf, "VzCorr2");
  hVzCorr2 = new TH2D(buf, buf, 400, -1000, 1000, 500, -1000, 1000);
*/

  return 1;
}

Int_t initTOFr()
{
  // Initialize TOFr

  //1/beta cut for pion
  piBetaCut = new TF1("piBetaCut","sqrt([0]*[0]+x*x)/x",0,5);
  piBetaCut->SetParameter(0,massCut);


  Char_t buff[100];
  for(Int_t i=0; i<nTray; i++) {
    for(Int_t j=0; j<nBoard; j++) {    //BUTTER LOOK LOOPS OVER TDIG
      for(Int_t k=0; k<ntCell; k++) {
	if(fixNBin) {
	  sprintf(buff, "TOT Tray_%d_Board_%d_Cell_%d", i+nTray*trayStart+1, j+TDIGStart*nBoard,k);
	  hTOT[i][j][k] = new TH2D(buff, buff, maxNBin, tofrTOTBins[i][j][k], ntoftdcbin, -toftlimit, toftlimit);
	  sprintf(buff, "Z Tray_%d_Board_%d_Cell_%d", i+nTray*trayStart+1, j+TDIGStart*nBoard,k);
	  hZ[i][j][k] = new TH2D(buff, buff, maxNZBin, tofrZBins[i][j][k], ntoftdcbin, -toftlimit, toftlimit);
	} else if (fixBinWidth) {
	  sprintf(buff, "TOT Tray_%d_Board_%d_Cell_%d", i+nTray*trayStart+1, j+TDIGStart*nBoard,k);
	  hTOT[i][j][k] = new TH2D(buff, buff, maxNBin, minTOT_TOFr, maxTOT_TOFr, ntoftdcbin, -toftlimit, toftlimit);
	  sprintf(buff, "Z Tray_%d_Board_%d_Cell_%d", i+nTray*trayStart+1, j+TDIGStart*nBoard,k);
	  hZ[i][j][k] = new TH2D(buff, buff, maxNZBin, minZ_TOFr, maxZ_TOFr, ntoftdcbin, -toftlimit, toftlimit);
	} else {
	  cout << " ? You have not decided which binning algrithm to be used !?" << endl;
	  return -1;
	}
	sprintf(buff, "Tray_%d_Board_%d_Cell%d_Time", i+nTray*trayStart+1, j+TDIGStart*nBoard,k);
	htTOT[i][j][k] = new TH1D(buff, buff, ntoftdcbin, -toftlimit, toftlimit);
      }
    }
  }

  Char_t title[100];
  sprintf(title, "hResoTOT");
  hResoTOT = new TH2D(title, title, 500, 0, 5, ntoftdcbin, -toftlimit, toftlimit);
  sprintf(title, "hResoZ");
  hResoZ = new TH2D(title, title, 500, 0, 5, ntoftdcbin, -toftlimit, toftlimit);
  sprintf(title, "hResoD");
  hResoD = new TH2D(title, title, 500, 0, 5, ntoftdcbin, -toftlimit, toftlimit);

  sprintf(title, "TOFr_Time_Resolution_after_TOT_Calibration");
  hResoTOT1D = new TH1D(title, title,  ntoftdcbin, -toftlimit, toftlimit);
  sprintf(title, "TOFr_Time_Resolution_after_Z_Calibration");
  hResoZ1D = new TH1D(title, title,  ntoftdcbin, -toftlimit, toftlimit);
  sprintf(title, "TOFr_Time_Resolution_after_D_Calibration");
  hResoD1D = new TH1D(title, title,  ntoftdcbin, -toftlimit, toftlimit);

  sprintf(title, "hResoTOT_Pos");
  hResoTOT_Pos = new TH2D(title, title, 500, 0, 5, ntoftdcbin, -toftlimit, toftlimit);
  sprintf(title, "hResoZ_Pos");
  hResoZ_Pos = new TH2D(title, title, 500, 0, 5, ntoftdcbin, -toftlimit, toftlimit);
  sprintf(title, "hResoD_Pos");
  hResoD_Pos = new TH2D(title, title, 500, 0, 5, ntoftdcbin, -toftlimit, toftlimit);

  sprintf(title, "hResoTOT1D_Pos");
  hResoTOT1D_Pos = new TH1D(title, title,  ntoftdcbin, -toftlimit, toftlimit);
  sprintf(title, "hReso_Z_Pos");
  hResoZ1D_Pos = new TH1D(title, title,  ntoftdcbin, -toftlimit, toftlimit);
  sprintf(title, "hReso_D_Pos");
  hResoD1D_Pos = new TH1D(title, title,  ntoftdcbin, -toftlimit, toftlimit);

  sprintf(title, "hResoTOT_Neg");
  hResoTOT_Neg = new TH2D(title, title, 500, 0, 5, ntoftdcbin, -toftlimit, toftlimit);
  sprintf(title, "hResoZ_Neg");
  hResoZ_Neg = new TH2D(title, title, 500, 0, 5, ntoftdcbin, -toftlimit, toftlimit);
  sprintf(title, "hResoD_Neg");
  hResoD_Neg = new TH2D(title, title, 500, 0, 5, ntoftdcbin, -toftlimit, toftlimit);

  sprintf(title, "hResoTOT1D_Neg");
  hResoTOT1D_Neg = new TH1D(title, title,  ntoftdcbin, -toftlimit, toftlimit);
  sprintf(title, "hReso_Z_Neg");
  hResoZ1D_Neg = new TH1D(title, title,  ntoftdcbin, -toftlimit, toftlimit);
  sprintf(title, "hReso_D_Neg");
  hResoD1D_Neg = new TH1D(title, title,  ntoftdcbin, -toftlimit, toftlimit);



  hnull = new TH1D("null", "null", 1, 0, 1);
  return 1;
}



Int_t doStart(TChain* chain)
{
  cout << "Come to doStart from pVPD "<< endl;
  Double_t par[nPVPDChannel][5];
  Double_t cor[nPVPDChannel];
  Double_t res[nPVPDChannel];
  Double_t tmpres[nPVPDChannel], tmpmean[nPVPDChannel];
  Double_t mean, sigma;
  Double_t nfitsigma = 2.0;//2.5; commented out .. rafael
  Double_t tdcwidth = 1.0; // ns

  TF1 *f1 = new TF1("f1","gaus");
  TF1 *polfit = new TF1("polfit","pol4");

  for(Int_t j=0; j<nPVPDChannel; j++) {
    cor[j] = 0.;
  }

  //--------------------
  // loop the chain
  //--------------------
  Int_t nevents = (int)chain->GetEntries();
  //  Int_t nevents = 20000; // test sample
  cout << "== total entries : " << nevents << endl;

  TOFrPicoDst *t = new TOFrPicoDst(chain);

  Char_t title[500];
  sprintf(title, "%s.ps", pvpdOutput);
  // TPostScript *ps = new TPostScript(title,112);
  TCanvas *c0 = new TCanvas("c0","c0",0,0,780,600);
  c0->Divide(3,2);
  c0->SetFillColor(10);
  TCanvas *c1 = new TCanvas("c1","c1",0,0,800,600);
  c1->SetFillColor(10);
  c1->Divide(3,2);
  TCanvas *c2 = new TCanvas("c2","c2",0,0,800,600);
  c2->SetFillColor(10);
  c2->Divide(3,2);
  TPostScript *ps = new TPostScript(title,112);
  ps->Range(26,20);

  gStyle->SetOptFit(1111);

  Bool_t done = kFALSE;
  Int_t ir =0;
  while(!done && ir<maxiteration) {
    cout << ir+1 << "-th iteration ... " << endl;

    if(ir==0||ir>=4)
      {
        for(Int_t jj=0; jj<nPVPDChannel; jj++) {
          htotPVPD[jj]->Reset();
        }
      }

    ps->NewPage();

    for(Int_t j=0; j<nPVPDChannel; j++) {
      hPVPD[j]->Reset();
      hPVPDCorr[j]->Reset(); // rafael
      htPVPD[j]->Reset();
      htPVPD1st[j]->Reset();
    }

    for(Int_t i=0; i<nevents; i++) {
      chain->GetEntry(i);

      Int_t flag[nPVPDChannel]; // 0: rejected, 1: accepted - rafael
      for(int i=0;i<nPVPDChannel;i++){
		 flag[i]=1; // rafael
//		if(i==7) flag[i]=0;//MASK
		}
      /*
      if(i%10000==0) {
        cout << "Working on event #" << i << endl;
      }
      */

      Int_t Ieast = t->numberOfVpdEast;
      Int_t Iwest = t->numberOfVpdWest;
      Double_t tot[nPVPDChannel];
      Double_t tdc[nPVPDChannel], tdccorr[nPVPDChannel];

      if(Iwest>=minpVPDHits) {

        Double_t average=0;
        Double_t tdcsum=0;
        Int_t largermswest = 0;
        Double_t countwest = 0;

        for(Int_t j=0; j<nPVPDChannel/2; j++) {
//		if(j==7) continue;  //MASK
          tot[j] = t->vpdTotWest[j];
          if(tot[j]<=minTOT_pVPD || tot[j]>maxTOT_pVPD) continue;
          tdc[j] = t->vpdLeWest[j];//-vpdDelay[j];//remove cable delay
          if(ir==0) { }
          else cor[j] = tsfPVPD[j]->V(tot[j]);

          if(ir>=4){
            tdccorr[j] = tdc[j]-cor[j];
            if(tdccorr[j]>51200) tdccorr[j] -= 51200;
            if(tdccorr[j]<0) tdccorr[j] += 51200;
            tdcsum = tdcsum + tdccorr[j];
            countwest = countwest + 1.0;
            //  cout<<"westhit  "<<countwest<<endl;
          }

        }//j 0-nPVPDChannel/2 loop
        //rms calculation ir>=4

        if(ir==maxiteration-1) {
          for(Int_t j=0; j<nPVPDChannel/2; j++) {
//		if(j==7) continue;//MASK

            //tot[j] = t->vpdTotWest[j];
            if(tot[j]<=minTOT_pVPD || tot[j]>maxTOT_pVPD) continue;
            for(Int_t k=0; k<nPVPDChannel/2; k++) {

//		if(k==7) continue;//MASK

              if(k==j) continue;
              if(tot[k]<=minTOT_pVPD || tot[k]>maxTOT_pVPD) continue;
              double deltatime = tdccorr[k] - tdccorr[j];
              //  if (deltatime < 0) deltatime += 51200;
              if (deltatime >51200) deltatime -= 51200;
              hdeltatwest->Fill(deltatime);
            }//k
          }//j
        } //if ir==maxiteration-1

        if(ir>=4) {
          average = tdcsum/countwest;
          Double_t rmssum=0;
          Double_t vpdHit[nPVPDChannel/2] = {0}; // outlier rejection - rafael
          Int_t *hitIndex = new Int_t[nPVPDChannel/2]; // outlier rejection - rafael
          for(Int_t j=0; j<nPVPDChannel/2; j++) {
		
//		if(j==7) continue; //MASK

            if(tot[j]<=minTOT_pVPD || tot[j]>maxTOT_pVPD) continue;
            if(fabs(tdccorr[j]-average)>tdiffCut) flag[j] = 0; // outlier rejection - rafael
            //  cout<<"rmssum west      "<<rmssum<<endl;
            rmssum = (tdccorr[j]-average)*(tdccorr[j]-average);
            vpdHit[j] = tdc[j]-cor[j]; // outlier rejection - rafael
          }
          Double_t rmswest = sqrt(rmssum/countwest);
          //cout<<"rms west     "<<rmswest<<endl;
          if (rmswest>mrmscut) largermswest = 1;
          // begin - outlier rejection - rafael
          int countHits = 0;
          TMath::Sort(nPVPDChannel/2,vpdHit,hitIndex);
          for(int i=0;i<nPVPDChannel/2;i++) if(vpdHit[i]!=0) countHits++;
          int nrejected = (int)(pCut*countHits+0.5);
          for(int i=0;i<nrejected;i++) flag[hitIndex[i]] = 0;
          delete hitIndex;
          //      cout << endl;
          //      cout << "number of hits: " << countHits << "\t hits rejected: " << nrejected << endl;
          //      for(int i=0;i<nPVPDChannel/2;i++) {
          //        cout << vpdHit[i];
          //        if(flag[i]==0) cout << " --> rejected!" << endl;
          //        else cout << endl;
          //      }
          // end - outlier rejection - rafael
        }

        if(!largermswest) {
          for(Int_t j=0; j<nPVPDChannel/2; j++) {
//		if(j==7) continue; //MASK

            if(!flag[j]) continue; // outlier rejection - rafael
            if(tot[j]<=minTOT_pVPD || tot[j]>maxTOT_pVPD) continue;
            if(ir==0||ir>=4)htotPVPD[j]->Fill(tot[j]);
            Double_t dt = 0.0;
            Int_t dn = 0;
            for(Int_t k=0; k<nPVPDChannel/2; k++) {
//		if(j==7) continue;//MASK
              if(k==j) continue;
              if(ir>2 && tmpres[j]>5) continue;
              if(tot[k]<=minTOT_pVPD || tot[k]>maxTOT_pVPD) continue;
              if(!flag[k]) continue; // outlier rejection - rafael
              //if(tot[k]<1e-6);
              dt += (tdc[k]-cor[k]);
              dn++;
              if(dn==(minpVPDHits-1)&&(ir==maxiteration-1)) break;//select the first minpVPDHits as reference at the last iteration
            }
            //      if(ir>0) cout<<"dn = "<<dn<<"  tmpmean/res = "<<tmpmean[j]<<"/"<<tmpres[j]<< endl;
            if(dn<minpVPDHits-1) continue;
            //cout<<" ++ "<<tdc[j]-cor[j]-dt/dn<<"  "<<tmpmean[j]/1000<<"  "<<tmpres[j]/1000<<"  "<<endl;
            if(ir>=0) {
              hPVPD[j]->Fill(tot[j], tdc[j]-dt/dn);
              hPVPDCorr[j]->Fill(tot[j],tdc[j]-cor[j]-dt/dn); // rafael
            }
            else if((tdc[j]-cor[j]-dt/dn)>(tmpmean[j]-5*tmpres[j]) && (tdc[j]-cor[j]-dt/dn)<(tmpmean[j]+5*tmpres[j])) {
              hPVPD[j]->Fill(tot[j], tdc[j]-dt/dn);
              hPVPDCorr[j]->Fill(tot[j],tdc[j]-cor[j]-dt/dn);// rafael
            }
            if(ir<=1) htPVPD1st[j]->Fill(tdc[j]-cor[j]-dt/dn);
            else htPVPD[j]->Fill(tdc[j]-cor[j]-dt/dn);
          }
        }//rms condition
      }//Iwest>minhits

      if(Ieast>=minpVPDHits) {

        Double_t average=0;
        Double_t tdcsum=0;
        Int_t largermseast = 0;
        Double_t counteast = 0;

        for(Int_t j=nPVPDChannel/2; j<nPVPDChannel; j++) {
          tot[j] = t->vpdTotEast[j-nPVPDChannel/2];
          if(tot[j]<=minTOT_pVPD || tot[j]>maxTOT_pVPD) continue;
          tdc[j] = t->vpdLeEast[j-nPVPDChannel/2];//-vpdDelay[j];//remove cable delay

          if(ir==0) { }
          else cor[j] = tsfPVPD[j]->V(tot[j]);

          if(ir>=4){
            tdccorr[j] = tdc[j]-cor[j];
            if(tdccorr[j]>51200) tdccorr[j] -= 51200;
            if(tdccorr[j]<0) tdccorr[j] += 51200;
            tdcsum = tdcsum + tdccorr[j];
            counteast = counteast + 1.0;
            //  cout<<"easthit  "<<counteast<<endl;
          }

        } // j=nPVPDChannel/2; j<nPVPDChannel loop

        if(ir==maxiteration-1) {
          for(Int_t j=nPVPDChannel/2; j<nPVPDChannel; j++) {
            if(tot[j]<=minTOT_pVPD || tot[j]>maxTOT_pVPD) continue;
            for(Int_t k=nPVPDChannel/2; k<nPVPDChannel; k++) {
              if(k==j) continue;
              if(tot[k]<=minTOT_pVPD || tot[k]>maxTOT_pVPD) continue;
              double deltatime = tdccorr[k] - tdccorr[j];
              //  if (deltatime < 0) deltatime += 51200;
              if (deltatime >51200) deltatime -= 51200;
              hdeltateast->Fill(deltatime);
            }//k
          }//j
        } //if ir==maxiteration-1

        //rms calculation ir>=4
        if(ir>=4) {
          average = tdcsum/counteast;
          Double_t rmssum=0;
          Double_t vpdHit[nPVPDChannel/2] = {0}; // outlier rejection - rafael
          Int_t *hitIndex = new Int_t[nPVPDChannel/2]; // outlier rejection - rafael
          for(Int_t j=nPVPDChannel/2; j<nPVPDChannel; j++) {
            if(tot[j]<=minTOT_pVPD || tot[j]>maxTOT_pVPD) continue;
            if(fabs(tdccorr[j]-average)>tdiffCut) flag[j] = 0; // outlier rejection - rafael
            //  cout<<"rmssum west      "<<rmssum<<endl;
            rmssum = (tdccorr[j]-average)*(tdccorr[j]-average);
            vpdHit[j-nPVPDChannel/2] = tdc[j] - cor[j]; // outlier rejection - rafael
          }
          Double_t rmseast = sqrt(rmssum/counteast);
          //cout<<"rms east     "<<rmseast<<endl;
          if (rmseast>mrmscut) largermseast = 1;
          // begin - outlier rejection - rafael
          int countHits = 0;
          TMath::Sort(nPVPDChannel/2,vpdHit,hitIndex);
          for(int i=0;i<nPVPDChannel/2;i++) if(vpdHit[i]!=0) countHits++;
          int nrejected = (int)(pCut*countHits+0.5);
          for(int i=0;i<nrejected;i++) flag[hitIndex[i]+nPVPDChannel/2] = 0;
          delete hitIndex;
          //      cout << endl;
          //      cout << "number of hits: " << countHits << "\t hits rejected: " << nrejected << endl;
          //      for(int i=0;i<nPVPDChannel/2;i++) {
          //        cout << vpdHit[i];
          //        if(flag[i+nPVPDChannel/2]==0) cout << " --> rejected!" << endl;
          //        else cout << endl;
          //      }
          // end - outlier rejection - rafael
        }

        if(!largermseast) {
          for(Int_t j=nPVPDChannel/2; j<nPVPDChannel; j++) {
            if(!flag[j]) continue; // outlier rejection - rafael
            if(tot[j]<=minTOT_pVPD || tot[j]>maxTOT_pVPD) continue;
            if(ir==0||ir>=4)htotPVPD[j]->Fill(tot[j]);
            Double_t dt = 0.0;
            Int_t dn = 0;
            for(Int_t k=nPVPDChannel/2; k<nPVPDChannel; k++) {
              if(k==j) continue;
              if(ir>2&&tmpres[j]>5) continue;
              if(tot[k]<=minTOT_pVPD || tot[k]>maxTOT_pVPD) continue;
              if(!flag[k]) continue; // outlier rejection - rafael
              dt += (tdc[k]-cor[k]);
              dn++;
              if((dn==minpVPDHits-1)&&(ir==maxiteration-1)) break;//select the first minpVPDHits as reference at the last iteration
            }

            if(dn<minpVPDHits-1) continue;
            if(ir>=0) {
              hPVPD[j]->Fill(tot[j], tdc[j]-dt/dn);
              hPVPDCorr[j]->Fill(tot[j],tdc[j]-cor[j]-dt/dn); // rafael
            }
            else if((tdc[j]-cor[j]-dt/dn)>(tmpmean[j]-5*tmpres[j]) && (tdc[j]-cor[j]-dt/dn)<(tmpmean[j]+5*tmpres[j])) {
              hPVPD[j]->Fill(tot[j], tdc[j]-dt/dn);
              hPVPDCorr[j]->Fill(tot[j],tdc[j]-cor[j]-dt/dn); // rafael
            }
            if(ir<=1) htPVPD1st[j]->Fill(tdc[j]-cor[j]-dt/dn);
            else htPVPD[j]->Fill(tdc[j]-cor[j]-dt/dn);
          }
        }//rms condition east
      }//east hit >= minhit

    } //end nevent loop



    if(ir==0||ir>=4){
      ofstream totmaxfile;
      totmaxfile.open("totmaxpos.dat");

      gStyle->SetOptStat(kFALSE);
      for(Int_t k=0; k<nPVPDChannel; k++) {
        c0->cd(k%(nPVPDChannel/6)+1);
        gPad->SetFillColor(10);

        htotPVPD[k]->Draw();
        int maxbin=htotPVPD[k]->GetMaximumBin();
        double maxpostot=htotPVPD[k]->GetBinCenter(maxbin);
        //cout<<k<<"    "<<maxbin<<"    "<<maxpostot<<" "<<htotPVPD[k]->GetMaximum()<<endl;
        totmaxfile<<k<<"        "<<maxbin<<"    "<<maxpostot<<" "<<htotPVPD[k]->GetMaximum()<<endl;
        if((k+1)%(nPVPDChannel/6)==0){
          c0->Update();
          ps->NewPage();
        }
      }
      c0->Update();
      ps->NewPage();
      gStyle->SetOptStat(1111);
      gStyle->SetOptFit(111);
      totmaxfile.close();
    }

    for(Int_t k=0; k<nPVPDChannel; k++) {
      if(k>0&&k%6==0) {
        c1->Update();
        ps->NewPage();
      }
      c1->cd(k%6+1);
      if(ir<=1){
/*	if(k==7){
		mean=0;
		sigma=0;
		tmpres[k]=0;
		tmpmean[k]=0;
		}//MASK
*/  
        mean = htPVPD1st[k]->GetMean();
        sigma = htPVPD1st[k]->GetRMS();
        f1->SetRange(mean-nfitsigma*sigma, mean+nfitsigma*sigma);
        htPVPD1st[k]->Fit("f1","QR");
        tmpres[k] = f1->GetParameter(2)*tdcwidth;//ns -> ps
        tmpmean[k] = f1->GetParameter(1)*tdcwidth;//ns -> ps
      } else {
        mean = htPVPD[k]->GetMean();
        sigma = htPVPD[k]->GetRMS();
        f1->SetRange(mean-nfitsigma*sigma, mean+nfitsigma*sigma);
        htPVPD[k]->Fit("f1","NQR");
        tmpres[k] = f1->GetParameter(2)*tdcwidth;
        tmpmean[k] = f1->GetParameter(1)*tdcwidth;//ns -> ps
        f1->SetRange((tmpmean[k]-nfitsigma*tmpres[k])/tdcwidth, (tmpmean[k]+nfitsigma*tmpres[k])/tdcwidth);
        htPVPD[k]->Fit("f1","NQR");
        tmpres[k] = f1->GetParameter(2)*tdcwidth;
        tmpmean[k] = f1->GetParameter(1)*tdcwidth;//ns -> ps
        f1->SetRange((tmpmean[k]-nfitsigma*tmpres[k])/tdcwidth, (tmpmean[k]+nfitsigma*tmpres[k])/tdcwidth);
        htPVPD[k]->Fit("f1","QR");
        tmpres[k] = f1->GetParameter(2)*tdcwidth;
        tmpmean[k] = f1->GetParameter(1)*tdcwidth;//ns -> ps

/*        if(k==7){
                mean=0;
                sigma=0;
                tmpres[k]=0;
                tmpmean[k]=0;
                }//MASK
*/
      }
    }
    c1->Update();
    ps->NewPage();

    cout << "========================================" << endl;
    cout << " Resolutions after last correction (ps) " << endl;
    cout << "========================================" << endl;
    cout << "  WEST pVPDs:   " << endl;
    for(Int_t j=0; j<nPVPDChannel/2; j++) {
      if(ir<=1) cout << "    " << tmpres[j] << "    " << htPVPD1st[j]->GetEntries() << endl;
      else cout << "    " << tmpres[j] << "    " << htPVPD[j]->GetEntries() << endl;
    }
    cout << "  EAST pVPDs:   " << endl;
    for(Int_t j=nPVPDChannel/2; j<nPVPDChannel; j++) {
      if(ir<=1) cout << "    " << tmpres[j] << "    " << htPVPD1st[j]->GetEntries() << endl;
      else cout << "    " << tmpres[j] << "    " << htPVPD[j]->GetEntries() << endl;
    }
    cout << "========================================" << endl;
    cout <<endl;
    cout << " Start Slewing Correction "<<endl;

    char testrootfile[500];
    sprintf(testrootfile, "%s_iter_%i.root", pvpdOutput,ir);
    //    TFile *testforwrite = new TFile(testrootfile,"RECREATE");
    //  testforwrite->cd();

    for(Int_t j=0; j<nPVPDChannel; j++) {
      if(j>0&&j%6==0){
        c2->Update();
        ps->NewPage();
      }
      c2->cd(j%6+1);

      Char_t buf[500];
      Char_t buff[500];

      TH1D *htemp1;
      TH1D *htemp2;
      TH1D *htemp3;
      TH1D *htemp4;

      sprintf(buf, "pVPD_Channel_%d_0", j);
      htemp1 = (TH1D *)gDirectory->FindObject(buf);
      if(htemp1) delete htemp1;
      sprintf(buf, "pVPD_Channel_%d_1", j);
      htemp2 = (TH1D *)gDirectory->FindObject(buf);
      if(htemp2) delete htemp2;
      sprintf(buf, "pVPD_Channel_%d_2", j);
      htemp3 = (TH1D *)gDirectory->FindObject(buf);
      if(htemp3) delete htemp3;
      sprintf(buf, "pVPD_Channel_%d_chi2", j);
      htemp4 = (TH1D *)gDirectory->FindObject(buf);
      if(htemp4) delete htemp4;

      //      hPVPD[j]->Write();
      sprintf(buf, "pVPD_Channel_%d_1", j);
      hPVPD[j]->FitSlicesY();
      hpPVPD[j] = (TH1D *)gDirectory->Get(buf); //Center is the bin center, be careful when bin width is large
      for(int k=0; k<maxpvpdNBin; k++){
        pvpdTOTX[j][k] = hpPVPD[j]->GetBinCenter(k+1);
        double content=hpPVPD[j]->GetBinContent(k+1);
        double meanyvalue = hPVPD[j]->GetMean(2);
        hPVPD[j]->GetXaxis()->SetRange(k+1,k+1);
        double meanvalueproj = hPVPD[j]->GetMean(2);
        double meanvalueprojerr = hPVPD[j]->GetMeanError(2);
        //       cout<<"channel "<<j+1<<"       iter    "<<ir+1<<endl;
        //       cout<<"fitmean "<<meanyvalue<<" x      "<<hPVPD[j]->GetMean(1)<<" projymean    "<<meanvalueproj<<"     dif     "<<meanyvalue-meanvalueproj<<"  fiterr  "<<hpPVPD[j]->GetBinError(k+1)<<"       proj errr    "<<meanvalueprojerr<<endl;
        double deltat=10;
        if(ir<4) deltat=5;
        else if(ir>=4) deltat=0.3;
        //        cout << " Channel " << j << " bin " << k << " Mean = " << meanvalueproj << " content = " << content << endl;

        if(fabs(content)>100||content>(meanvalueproj+deltat)||content<(meanvalueproj-deltat)) {
          //        cout << " Channel " << j << " bin " << k << " Mean = " << meanvalueproj << " content = " << content << endl;
          hpPVPD[j]->SetBinContent(k+1,meanvalueproj);
          hpPVPD[j]->SetBinError(k+1,meanvalueprojerr);
        }

      }
      //  hpPVPD[j]->Write();
      sprintf(buf, "pVPD Fit Channel_%d", j);
      sprintf(buff, "pVPD Fit | Channel_%d", j);
      if(tsfPVPD[j]) delete tsfPVPD[j];
      tsfPVPD[j] = new TSplineFit(buf, buff, maxpvpdNBin, tsfM, hpPVPD[j]);
      //cout<<j<<"  "<<hpPVPD[j]->GetMinimum()<<"  "<<hpPVPD[j]->GetMaximum()<<endl;
      tsfPVPD[j]->DrawHere(max((Int_t)hpPVPD[j]->GetMinimum()-2, (Int_t) -tlimit), min(int(hpPVPD[j]->GetMaximum())+2, (Int_t) tlimit));
      /*
        if(tsfPVPD[j]->Chi2()>1000)
        {
        //     polfit->SetParameters(5,0.2,0.005,0.0001,0.00001);
        //      hpPVPD[j]->Fit(polfit, "ERMB", "", pvpdTOTX[j][0], pvpdTOTX[j][maxNBin-1]);
        //        polfit->GetParameters(par[j]);
        //        cout<<"channel  "<<j+1<<"   iter        "<<ir+1<<"    polfit  "<<endl;

        //  cout<<"par[0]       "<<par[j][0]<<" par[1]  "<<par[j][1]<<" par[2]  "<<par[j][2]<<" par[3]  "<<par[j][3]<<" par[4]  "<<par[j][4]<<endl;
        }

        //cout<<"channel        "<<j+1<<"       iter    "<<ir+1<<endl;
        //cout<<"chi2   "<<tsfPVPD[j]->Chi2()<<endl;
        for(int k=0; k<maxNBin; k++)
        {double meandeltat = hpPVPD[j]->GetBinContent(k+1);
        double splinefitdeltat = tsfPVPD[j]->V(pvpdTOTX[j][k]);
        //       cout<<"meandeltat      "<<setw(10)<<meandeltat<<"  splinefit   "<<setw(10)<<splinefitdeltat<<" diff    "
        //              <<setw(10)<<meandeltat-splinefitdeltat<<endl;
        }
      */

    }

    // testforwrite->Close();

    c2->Update();

    // begin check file - rafael
    TFile *f = new TFile("pvpdCheck.root","UPDATE");
    char buf[100];
    for(Int_t k=0; k<nPVPDChannel; k++) {
      sprintf(buf,"pVPDCorr_Ch%d_it%d",k,ir+1);
      hPVPDCorr[k]->SetName(buf);
      hPVPDCorr[k]->Write();
    }
    f->Close();
    delete f;
    // end check file - rafael


    done = kTRUE;
    for(Int_t k=0; k<nPVPDChannel; k++) {
      //compare time resolution to last iteration
      if(TMath::Abs(res[k]-tmpres[k])>smallvalue) { done = kFALSE; }
      res[k] = tmpres[k];
      // begin - rafael
      if(ir<=1) hPVPDMeanvsIt[k]->SetBinContent(ir+1,htPVPD1st[k]->GetMean());
      else hPVPDMeanvsIt[k]->SetBinContent(ir+1,htPVPD[k]->GetMean());
      hPVPDFitMeanvsIt[k]->SetBinContent(ir+1,f1->GetParameter(1));
      hPVPDFitSigmavsIt[k]->SetBinContent(ir+1,f1->GetParameter(2));
      hPVPDResovsIt[k]->SetBinContent(ir+1,res[k]);
      // end - rafael
    }
    ir++;
  } // end while

  ps->Close();

  // output calibration results

  cout << "==============================" << endl;
  cout << endl;
  cout << "  Correction over !  Writing out ... " << endl;

  return 1;
}




Int_t doDelay(TChain* chain, Int_t iter) 
{
  cout << "Come to doDelay from TOFr "<< endl; 

  char rootfile[500];

   Int_t nevt = (int)(chain->GetEntries());

  TOFrPicoDst *tofr = new TOFrPicoDst(chain);
  if(iter==0) {
/*    sprintf(rootfile, "%s.root", pvpdOutput);
    TChain *startchain = new TChain("StartPicoDst");
    startchain->AddFile(rootfile);
    StartPicoDst *start = new StartPicoDst(startchain);

    if(nevt!=startchain->GetEntries()) {
      cout << " Inconsistent TOFr and Start entry ! " << endl;
      return -1;
    }
*/
    cout << "== total entries : " << nevt << endl; 

    // fill delay histogram
    for (Int_t i=0;i<nevt;i++) {
      tofr->GetEntry(i);
  //    start->GetEntry(i);
      if(i%(nevt/10)==0) cout << i << " events processed!" << endl;

      Int_t nhits = tofr->nTofHits;
      Double_t t0 = tofr->T0;
      Float_t vzVpd = tofr->vzVpd;
//	cout<<"Do delay hits: "<<nhits<<endl;
      for(Int_t j=0; j<nhits; j++) {
        Double_t dedx = tofr->dedx[j];
        Double_t p = tofr->pt[j]*TMath::CosH(tofr->eta[j]);
        if(p<0.3 || p>0.6 || dedx > dedxcut) continue; // is pion?
        if(sqrt(pow(tofr->vertexX,2)+pow(tofr->vertexY,2))>rCut) continue;//for AuAu
//        if(sqrt(pow(tofr->dcaX[j],2)+pow(tofr->dcaY[j],2))>rCut) continue;//for pprun9

        if(TMath::Abs(tofr->vertexZ-vzVpd)>vzCut) continue;//for AuAu
//        if(TMath::Abs(tofr->dcaZ[j]-vzVpd)>vzCut) continue;//for pprun9


        Double_t tracklength = tofr->length[j];
        Double_t betagamma = p/piMass;
        Double_t velocity = TMath::Sqrt(1.0/(1.0/betagamma/betagamma+1.0))*c_light;
        Double_t tofPi = tracklength/velocity;

        Int_t tray = tofr->tray[j]-1-nTray*trayStart;
        Int_t module = (tofr->module[j]-1)%4;//-moduleStart*nModule;
//	Int_t tmodule = module%4;
        Int_t cell = tofr->cell[j]-1;
        Int_t tcell = module*6+cell;
//	cout <<"tray, module, cell, tcell: "<<tray<<", "<<module<<", "<<cell<<", "<<tcell<<endl;
        Double_t tdc = tofr->leTime[j];
        tdc -= phaseDiff;
        //while(tdc>51200) tdc -= 51200;
        if(tray>=nTray || tray<0) continue;
        if(module>=nModule || module<0) continue;
        if(cell>=nCell || cell<0) continue;
//	cout<<"Made it through trays,modules,cellscut"<<endl;
        Double_t dt = tdc-t0-tofPi;
        if(TMath::Abs(dt)<10*bigtime) {
//		cout<<"dt: "<<dt<<endl;
          hDelay[tray][0][tcell]->Fill(dt);
        }
      }
    } // end nevt loop
  } else {
    sprintf(rootfile, "tray_%d_board_%d_%s_%d.root",trayStart, TDIGStart,zOutput, iter);
    TChain *zhitchain = new TChain("ZPicoDst");
    zhitchain->AddFile(rootfile);
    ZPicoDst *zhit = new ZPicoDst(zhitchain);

    if(nevt!=zhitchain->GetEntries()) {
      cout << " Inconsistent TOFr and Zhit entry ! " << endl;
      return -1;
    }
    for (Int_t i=0;i<nevt;i++) {
      tofr->GetEntry(i);
      zhit->GetEntry(i);
    if(i%(nevt/10)==0) cout << i << " events processed!" << endl;

      Int_t nhits = tofr->nTofHits;

      for(Int_t j=0; j<nhits; j++) {
        Double_t dedx = tofr->dedx[j];
        Double_t p = tofr->pt[j]*TMath::CosH(tofr->eta[j]);
        if(p<0.3||p>1.) continue;
        if(p<0.6&&dedx > dedxcut) continue; // dedx cut at p<0.6
        if(sqrt(pow(tofr->vertexX,2)+pow(tofr->vertexY,2))>rCut) continue;//for AuAu--remember not really necessary
        if(TMath::Abs(tofr->vertexZ-zhit->vzVpd)>vzCut) continue;//for AuAu
//        if(sqrt(pow(tofr->dcaX[j],2)+pow(tofr->dcaY[j],2))>rCut) continue;//for pprun9
//        if(TMath::Abs(tofr->dcaZ[j]-zhit->vzVpd)>vzCut) continue;//for pprun9


        if(tofr->nHitsFit[j]<nHitsFitCut) continue;

        Double_t tracklength = tofr->length[j];
        Double_t betagamma = p/piMass;
        Double_t velocity = TMath::Sqrt(1.0/(1.0/betagamma/betagamma+1.0))*c_light;
        Double_t tofPi = tracklength/velocity;

        Int_t tray = tofr->tray[j]-1-nTray*trayStart;
        Int_t module = (tofr->module[j]-1)%4;//-moduleStart*nModule;
//	Int_t tmodule =module%4;
        Int_t cell = tofr->cell[j]-1;
	Int_t tcell = module*6+cell;
        if(tray>=nTray || tray<0) continue;
        if(module>=nModule || module<0) continue;
        if(cell>=nCell || cell<0) continue;
        Double_t tof = zhit->tof[j];

        Double_t invBeta = tof/tracklength*c_light;
        if(invBeta>piBetaCut->Eval(p)) continue; //TOF cut at p<1.;

        Double_t dt = tof-tofPi;
        if(TMath::Abs(dt)<20*toftlimit) {
          hDelay2[tray][0][tcell]->Fill(dt);
        }
      }
    } // end nevt loop
  }

  Char_t title[500];
  sprintf(title, "tray_%d_board_%d_%s_%d.ps",trayStart, TDIGStart,delayOutput, iter);
/*  TCanvas *c1 = new TCanvas("c1","c1",0,0,800,600);
  c1->SetFillColor(10);
  c1->SetBorderSize(2);
  c1->SetBorderMode(0);
  c1->Divide(4,3);
*/
  TPostScript *ps = new TPostScript(title,112);
  ps->Range(26,20);

  TCanvas *c1 = new TCanvas("c1","c1",0,0,780,600);
  c1->SetFillColor(10);
  c1->Divide(3,2);

  gStyle->SetOptStat(111111);
  gStyle->SetOptFit(1111);
  TF1 *g = new TF1("g","gaus");
  g->SetLineWidth(1);
  g->SetLineColor(2);
  for(Int_t i=0; i<nTray; i++) {
    for(Int_t j=0; j<nBoard; j++) {
      c1->Update();
      ps->NewPage();

      for(Int_t k=0; k<nCell; k++) {

 
        c1->cd(k+1);
        float tmpmean, tmpres, histmean, histres, nsigma = 2.5;
	    float tmpresfit, tmpmeanfit;
	    float maxres=1.5;
		float minres=0.2;
		float tmpcount=30;
		float fitlow, fithi, meandifcut;
		
		if(iter>0&&iter<3) {maxres=0.8; minres = 0.15; meandifcut = 0.5;}
		if(iter>=3&&iter<4) {maxres=0.4; minres = 0.05; meandifcut = 0.2;}
		if(iter>=4) {maxres=0.2; minres = 0.05; meandifcut = 0.1;}

		
        if(iter==0) {
		 hDelay[i][j][k]->GetXaxis()->UnZoom();	
		// hDelay[i][j][k]->GetXaxis()->SetRangeUser(-1.5*bigtime, 0.5*bigtime);			
	     tmpmean = hDelay[i][j][k]->GetMean();
         tmpres = hDelay[i][j][k]->GetRMS();
		 if(tmpres >= maxres)
		 {
		 	 if(hDelay[i][j][k]->GetMaximum()>1){
		 	 int tmpmaxbin = hDelay[i][j][k]->GetMaximumBin();
			 tmpmean = hDelay[i][j][k]->GetBinCenter(tmpmaxbin);}
			 
			 if(hDelay[i][j][k]->GetMaximum()<=1)
			 	{for(int tmpct=0; tmpct<7; tmpct++)
			 	 {hDelay[i][j][k]->Rebin(2);
			      if(hDelay[i][j][k]->GetMaximum()>1){
				  int tmpmaxbin = hDelay[i][j][k]->GetMaximumBin();
			      tmpmean = hDelay[i][j][k]->GetBinCenter(tmpmaxbin);
				  break;
			      } //if
			 	}//for loop
			  } //if maxbincount ==1
		 }
		 
		 for(int ltmp=0; ltmp<5; ltmp++)
		 {
		 tmpres = fabs(tmpres);
	     tmpres = TMath::Min(tmpres,maxres);		 
	     tmpres = TMath::Max(minres, tmpres);
		 hDelay[i][j][k]->GetXaxis()->SetRangeUser(tmpmean-nsigma*tmpres, tmpmean+nsigma*tmpres);
		 histmean = hDelay[i][j][k]->GetMean();
         histres = hDelay[i][j][k]->GetRMS(); 
		 if(fabs(histmean-tmpmean)<1e-3){ break;}
		 tmpmean = histmean;
		 tmpres = histres;
		 // cout<<"ltmp	"<<ltmp<<"	deltadealy	"<<histmean-tmpmean<<endl;
		 }
		 tmpmean = histmean;
		 tmpres = histres;
		 /*
		 tmpres = fabs(tmpres);
	     tmpres = TMath::Min(tmpres,maxres);		 
	     tmpres = TMath::Max(minres, tmpres);
         */


		// hDelay[i][j][k]->GetXaxis()->SetRangeUser(-1.5*bigtime, 0.5*bigtime);
         hDelay[i][j][k]->GetXaxis()->UnZoom();
		/* 
	     g->SetParameters(30,tmpmean,tmpres);
         hDelay[i][j][k]->Fit("g","NQ","",tmpmean-nsigma*tmpres, tmpmean+nsigma*tmpres);
	     tmpcount = g->GetParameter(0);
         tmpmeanfit = g->GetParameter(1);
         tmpresfit = g->GetParameter(2);
		 tmpresfit = fabs(tmpresfit);
		  
      	 g->SetParameters(tmpcount,tmpmeanfit,tmpresfit);		  
         hDelay[i][j][k]->Fit("g","NQ","",tmpmeanfit-nsigma*tmpresfit, tmpmeanfit+nsigma*tmpresfit);
		 tmpcount = g->GetParameter(0);
         tmpmeanfit = g->GetParameter(1);
         tmpresfit = g->GetParameter(2);
		 tmpresfit = fabs(tmpresfit);		 
          
      	 g->SetParameters(tmpcount,tmpmeanfit,tmpresfit);
         hDelay[i][j][k]->Fit("g","Q","",tmpmeanfit-nsigma*tmpresfit, tmpmeanfit+nsigma*tmpresfit);
         tmpmeanfit = g->GetParameter(1);
         tmpresfit = g->GetParameter(2);

          if(fabs(tmpmeanfit-tmpmean)<meandifcut) {tmpmean=tmpmeanfit; tmpres=tmpresfit;}
		  else{
		  cout<<"tray	"<<i<<"	module	"<<j<<"	cell	"<<k<<"	check delay fit iter= "<<iter<<endl;
		  cout<<"refit	"<<endl;		  	
          g->SetParameters(30,tmpmean,tmpres);
		  g->SetParLimits(1,tmpmean-nsigma*tmpres, tmpmean+nsigma*tmpres);
		  g->SetParLimits(2,0,nsigma*tmpres);
          hDelay[i][j][k]->Fit("g","NQ","",tmpmean-nsigma*tmpres, tmpmean+nsigma*tmpres);
		  tmpcount = g->GetParameter(0);
          tmpmeanfit = g->GetParameter(1);
          tmpresfit = g->GetParameter(2); 
		  tmpresfit = fabs(tmpresfit);		  

      	  g->SetParameters(tmpcount,tmpmeanfit,tmpresfit);	
		  g->SetParLimits(1,tmpmeanfit-nsigma*tmpresfit, tmpmeanfit+nsigma*tmpresfit);
		  g->SetParLimits(2,0,nsigma*tmpresfit);		  
          hDelay[i][j][k]->Fit("g","NQ","",tmpmeanfit-nsigma*tmpresfit, tmpmeanfit+nsigma*tmpresfit);
		  tmpcount = g->GetParameter(0);
          tmpmeanfit = g->GetParameter(1);
          tmpresfit = g->GetParameter(2);
		  tmpresfit = fabs(tmpresfit);		  
          
      	  g->SetParameters(tmpcount,tmpmeanfit,tmpresfit);
		  g->SetParLimits(1,tmpmeanfit-nsigma*tmpresfit, tmpmeanfit+nsigma*tmpresfit);
		  g->SetParLimits(2,0,nsigma*tmpresfit);
          hDelay[i][j][k]->Fit("g","Q","",tmpmeanfit-nsigma*tmpresfit, tmpmeanfit+nsigma*tmpresfit);
          tmpmeanfit = g->GetParameter(1);
          tmpresfit = g->GetParameter(2);		  

		  if(fabs(tmpmeanfit-tmpmean)<meandifcut) {tmpmean=tmpmeanfit; tmpres=tmpresfit;}
		  }
		  */
		  hDelay[i][j][k]->Draw();
          hDelay[i][j][k]->GetXaxis()->SetRangeUser(tmpmean-20*tmpres, tmpmean+20*tmpres);
		  
		  if(hDelay[i][j][k]->GetEntries()<1) { tmpmean=0; tmpres=0;}
          
	  //cout<<" +++ delay "<<i<<"  "<<j<<"  "<<k<<"  "<<tmpmean<<"  "<<hDelay[i][j][k]->GetMean()<<endl;
	  //tofrDelay[i][j][k] = hDelay[i][j][k]->GetMean();
	  
        }else {
        //  cout<<"1 "<<endl; 
// 		  hDelay2[i][j][k]->GetXaxis()->UnZoom();       
  		  hDelay2[i][j][k]->GetXaxis()->SetRangeUser(-2*toftlimit, 2*toftlimit);  
		 if(hDelay2[i][j][k]->GetMaximum()>1){
		 int tmpmaxbin = hDelay2[i][j][k]->GetMaximumBin();
		 tmpmean = hDelay2[i][j][k]->GetBinCenter(tmpmaxbin);}
         else{tmpmean = hDelay2[i][j][k]->GetMean();}
          tmpres = hDelay2[i][j][k]->GetRMS();
		  for(int ltmp=0; ltmp<5; ltmp++)
		  {
		  tmpres = fabs(tmpres);
	      tmpres = TMath::Min(tmpres,maxres);		 
	      tmpres = TMath::Max(minres, tmpres);
		  hDelay2[i][j][k]->GetXaxis()->SetRangeUser(tmpmean-nsigma*tmpres, tmpmean+nsigma*tmpres);
		  histmean = hDelay2[i][j][k]->GetMean();
		  histres = hDelay2[i][j][k]->GetRMS(); 
		  if(fabs(histmean-tmpmean)<1e-3){ break;}
		  tmpmean = histmean;
		  tmpres = histres;
		  // cout<<"ltmp	"<<ltmp<<"	deltadealy	"<<histmean-tmpmean<<endl;
		  }
		  tmpmean = histmean;
		  tmpres = histres;

		  
		//  hDelay2[i][j][k]->GetXaxis()->SetRangeUser(-2*toftlimit, 2*toftlimit);	
		  hDelay2[i][j][k]->GetXaxis()->UnZoom();
        //  cout<<"2 "<<endl;
		  
	    //  tmpres = TMath::Min(tmpres,maxres);
	  	//  tmpres = TMath::Max(minres, tmpres);
          /*
		  fitlow = tmpmean - fabs(nsigma*tmpres);
		  fithi = tmpmean + fabs(nsigma*tmpres);
		  if(fitlow < -2*toftlimit+0.01) fitlow = -2*toftlimit+0.01;
		  if(fitlow > 2*toftlimit-0.02) fitlow = 2*toftlimit-0.02;
		  if(fithi < -2*toftlimit+0.02) fithi = -2*toftlimit+0.02;
		  if(fithi > 2*toftlimit-0.01)  fithi = 2*toftlimit-0.01;	
		//  cout<<"fitlow	"<<fitlow<<"	fithi	"<<fithi<<endl;
		  
	      g->SetParameters(30,tmpmean,tmpres);		 
          hDelay2[i][j][k]->Fit("g","NQ","",fitlow, fithi);
	      tmpcount = g->GetParameter(0);
          tmpmeanfit = g->GetParameter(1);
          tmpresfit = g->GetParameter(2);
		  
          g->SetParameters(tmpcount,tmpmeanfit,tmpresfit);
		  fitlow = tmpmeanfit - fabs(nsigma*tmpresfit);
		  fithi = tmpmeanfit + fabs(nsigma*tmpresfit);	
		  if(fitlow < -2*toftlimit+0.01) fitlow = -2*toftlimit+0.01;
		  if(fitlow > 2*toftlimit-0.02) fitlow = 2*toftlimit-0.02;
		  if(fithi < -2*toftlimit+0.02) fithi = -2*toftlimit+0.02;
		  if(fithi > 2*toftlimit-0.01)  fithi = 2*toftlimit-0.01;		  
		  
          hDelay2[i][j][k]->Fit("g","NQ","",fitlow, fithi);
	      tmpcount = g->GetParameter(0);		  
          tmpmeanfit = g->GetParameter(1);
          tmpresfit = g->GetParameter(2);
          
          g->SetParameters(tmpcount,tmpmeanfit,tmpresfit);
		  fitlow = tmpmeanfit - fabs(nsigma*tmpresfit);
		  fithi = tmpmeanfit + fabs(nsigma*tmpresfit);	
		  if(fitlow < -2*toftlimit+0.01) fitlow = -2*toftlimit+0.01;
		  if(fitlow > 2*toftlimit-0.02) fitlow = 2*toftlimit-0.02;
		  if(fithi < -2*toftlimit+0.02) fithi = -2*toftlimit+0.02;
		  if(fithi > 2*toftlimit-0.01)  fithi = 2*toftlimit-0.01;		  
          hDelay2[i][j][k]->Fit("g","Q","",fitlow, fithi);
	      tmpcount = g->GetParameter(0);		  
          tmpmeanfit = g->GetParameter(1);
          tmpresfit = g->GetParameter(2);
		  
		  if(fabs(tmpmeanfit-tmpmean)<meandifcut) {tmpmean=tmpmeanfit; tmpres=tmpresfit;}
		  else{
			  cout<<"tray "<<i<<" module	"<<j<<" cell	"<<k<<" check delay fit iter= "<<iter<<endl;
			  cout<<"refit	"<<endl;		  	
			  g->SetParameters(30,tmpmean,tmpres);
			  g->SetParLimits(1,tmpmean-fabs(nsigma*tmpres), tmpmean+fabs(nsigma*tmpres));
			  g->SetParLimits(2,0,nsigma*fabs(tmpres));
			  fitlow = tmpmean - fabs(nsigma*tmpres);
		      fithi = tmpmean + fabs(nsigma*tmpres);
		  	  if(fitlow < -2*toftlimit+0.01) fitlow = -2*toftlimit+0.01;
		  	  if(fitlow > 2*toftlimit-0.02) fitlow = 2*toftlimit-0.02;
		  	  if(fithi < -2*toftlimit+0.02) fithi = -2*toftlimit+0.02;
		  	  if(fithi > 2*toftlimit-0.01)  fithi = 2*toftlimit-0.01;	 
			  hDelay2[i][j][k]->Fit("g","NQ","",fitlow, fithi);
			  tmpcount = g->GetParameter(0);
			  tmpmeanfit = g->GetParameter(1);
			  tmpresfit = g->GetParameter(2); 
			  g->SetParameters(tmpcount,tmpmeanfit,tmpresfit);	
			  g->SetParLimits(1,tmpmeanfit-nsigma*fabs(tmpresfit), tmpmeanfit+nsigma*fabs(tmpresfit));
			  g->SetParLimits(2,0,nsigma*fabs(tmpresfit));
			  fitlow = tmpmeanfit - fabs(nsigma*tmpresfit);
		 	  fithi = tmpmeanfit + fabs(nsigma*tmpresfit);	
		      if(fitlow < -2*toftlimit+0.01) fitlow = -2*toftlimit+0.01;			  
		      if(fitlow > 2*toftlimit-0.02) fitlow = 2*toftlimit-0.02;
		      if(fithi < -2*toftlimit+0.02) fithi = -2*toftlimit+0.02;
		  	  if(fithi > 2*toftlimit-0.01)  fithi = 2*toftlimit-0.01;
			  hDelay2[i][j][k]->Fit("g","NQ","",fitlow, fithi);
			  tmpcount = g->GetParameter(0);		  
			  tmpmeanfit = g->GetParameter(1);
			  tmpresfit = g->GetParameter(2);
			  
			  g->SetParameters(tmpcount,tmpmeanfit,tmpresfit);
			  g->SetParLimits(1,tmpmeanfit-nsigma*fabs(tmpresfit), tmpmeanfit+nsigma*fabs(tmpresfit));
			  g->SetParLimits(2,0,nsigma*fabs(tmpresfit));	
			  fitlow = tmpmeanfit - fabs(nsigma*tmpresfit);
		      fithi = tmpmeanfit + fabs(nsigma*tmpresfit);	
		 	  if(fitlow < -2*toftlimit+0.01) fitlow = -2*toftlimit+0.01;
		      if(fitlow > 2*toftlimit-0.02) fitlow = 2*toftlimit-0.02;
		      if(fithi < -2*toftlimit+0.02) fithi = -2*toftlimit+0.02;
		      if(fithi > 2*toftlimit-0.01)  fithi = 2*toftlimit-0.01;		  
			  hDelay2[i][j][k]->Fit("g","Q","",fitlow, fithi);
			  tmpcount = g->GetParameter(0);		  
			  tmpmeanfit = g->GetParameter(1);
			  tmpresfit = g->GetParameter(2);			  

			  if(fabs(tmpmeanfit-tmpmean)<meandifcut) {tmpmean=tmpmeanfit; tmpres=tmpresfit;}
		  }
		  */
		hDelay2[i][j][k]->Draw();  
         // hDelay2[i][j][k]->GetXaxis()->SetRangeUser(tmpmean-5*tmpres, tmpmean+5*tmpres);
		//hDelay2[i][j][k]->GetXaxis()->SetRangeUser(-2*toftlimit, 2*toftlimit);  
	    if(hDelay2[i][j][k]->GetEntries()<1) {tmpmean=0; tmpres=0;}

	  	}

        tofrDelay[i][j][k] = tmpmean;
	    tofrRes[i][j][k] = tmpres;
//     cout<<"tofrDelay okay?: "<<tofrDelay[i][j][k]<<endl;
      	}//cell loop
    }  // module loop
  } // end nTray loop
  ps->Close();
  for(int i=0;i<nTray;i++){
    for(int j=0;j<nBoard;j++){
      for(int k=0;k<nCell;k++){
        if(iter==0) hDelay[i][j][k]->Reset();
        else hDelay2[i][j][k]->Reset();
      }
    }
  }
  outputDelay(chain,iter);

  return 1;
}

Int_t doTOT(TChain* chain, Int_t iter) 
{
  cout << "Come to doTOT from TOFr, " << iter << " iteration. " << endl; 
  if(iter<0) {
    cout << "Wrong iteration number! exit! " << endl;
    return -1;
  }
  
  Int_t N;
  Double_t X[maxNBin];
  Double_t Y[maxNBin];
  Double_t EY[maxNBin];
  
  TOFrPicoDst *tofr = new TOFrPicoDst(chain);
  Int_t nevt = (int)(chain->GetEntries());

  char rootfile[500];
  sprintf(rootfile, "%s.root", pvpdOutput);
/*  TChain *startchain = new TChain("StartPicoDst");
  startchain->AddFile(rootfile);
  StartPicoDst *start = new StartPicoDst(startchain);
*/
  sprintf(rootfile, "tray_%d_board_%d_%s_%d.root",trayStart,TDIGStart, delayOutput, iter-1);

  TChain *delaychain = new TChain("DelayPicoDst");
  delaychain->AddFile(rootfile);
  DelayPicoDst *delay = new DelayPicoDst(delaychain);

  if(nevt!=delaychain->GetEntries()) {
    cout << " Inconsistent TOFr and Delay entry ! " << endl;
    return -1;
  }

  // fill slewing histogram
  for (Int_t i=0;i<nevt;i++) {
    tofr->GetEntry(i);
  //  start->GetEntry(i);
    delay->GetEntry(i);

    if(i%(nevt/10)==0) cout << i << " events processed!" << endl;

    Int_t ieast = tofr->Ieast;
    Int_t iwest = tofr->Iwest;
    if(ieast==0||iwest==0) continue;
    Int_t nhits = tofr->nTofHits;
    
    for(Int_t j=0; j<nhits; j++) {
      Double_t dedx = tofr->dedx[j];
      Double_t p = tofr->pt[j]*TMath::CosH(tofr->eta[j]);
      if(p<0.3 || p>0.6 || dedx > dedxcut) continue; // is pion?
      if(sqrt(pow(tofr->vertexX,2)+pow(tofr->vertexY,2))>rCut) continue;//for AuAu
      if(TMath::Abs(tofr->vertexZ-delay->vzVpd)>vzCut) continue;//for AuAu
//      if(sqrt(pow(tofr->dcaX[j],2)+pow(tofr->dcaY[j],2))>rCut) continue;//for pprun9
//      if(TMath::Abs(tofr->dcaZ[j]-delay->vzVpd)>vzCut) continue;//for pprun9


      //      Double_t nSigPi = tofr->nSigPi[j];
      //      Double_t nSigK = tofr->nSigK[j];
      //      Double_t nSigP = tofr->nSigP[j];
      //      Double_t nSigE = tofr->nSigE[j];
      //      if(p<0.2 || p>1.0 || TMath::Abs(nSigPi)>2.0 || TMath::Abs(nSigK)<2.0 || TMath::Abs(nSigP)<2.0 || TMath::Abs(nSigE)<2.0) continue; // another way to choose pion
      
      Double_t tracklength = tofr->length[j];
      Double_t betagamma = p/piMass;
      Double_t velocity = TMath::Sqrt(1.0/(1.0/betagamma/betagamma+1.0))*c_light;
      Double_t tofPi = tracklength/velocity;
      
      Int_t tray = tofr->tray[j]-1-nTray*trayStart;
      Int_t module = (tofr->module[j]-1)%4;//-moduleStart*nModule;
      Int_t cell = tofr->cell[j]-1;
      Int_t tcell = module*6+cell;
      Double_t tot = tofr->tot[j];
      
      if(tray>=nTray || tray<0) continue;
      if(module>=nModule || module<0) continue;
      if(cell>=nCell || cell<0) continue;
      
      Double_t tof = delay->tof[j];
      Double_t dt = tof-tofPi;
      if(TMath::Abs(dt)<toftlimit) {
        if(tofr->charge[j]==1) hplus->Fill(dt);
        if(tofr->charge[j]==-1) hminus->Fill(dt);

        hTOT[tray][0][tcell]->Fill(tot, dt);
        htTOT[tray][0][tcell]->Fill(dt);
      }
    }
  } // end nevt loop

  // do tot calibration
  Char_t title[500];
  sprintf(title, "tray_%d_board_%d_%s_%d.ps", trayStart, TDIGStart,totOutput, iter);
  TCanvas *c1 = new TCanvas("c1","c1",0,0,800,600);
  c1->SetFillColor(10);
  c1->SetBorderSize(2);
  c1->SetBorderMode(0);
  c1->Divide(3,3);
  TPostScript *ps = new TPostScript(title,112);
  ps->Range(26,20);
  for(Int_t i=0; i<nTray; i++) {
    ps->NewPage();
    for(Int_t j=0; j<nBoard; j++) {  ///BUTTER LOOK HERE TDIG LOOP
      // if(j%4==0) ps->NewPage(); 
      for(Int_t k=0; k<ntCell; k++) {
        if(k%4==0) ps->NewPage();
	c1->cd(2*(k%4)+1);
    tofrTOTbad[i][j][k]=0;
	     hTOT[i][j][k]->Draw();
        c1->cd(2*(k%4)+2);

	Int_t entry = (int)htTOT[i][j][k]->GetEntries();
	Double_t tm = htTOT[i][j][k]->GetMean();
	Double_t ts = htTOT[i][j][k]->GetRMS();
	N=0;
	cout << setw(4) << i+trayStart*nTray <<  setw(4) << j+TDIGStart*nBoard  <<  setw(4) << k << endl;
	cout << entry << " events" << endl;

	hTOT[i][j][k]->FitSlicesY();
	Char_t buf[500];
	sprintf(buf, "TOT Tray_%d_Board_%d_Cell_%d_1", i+nTray*trayStart+1, j+TDIGStart*nBoard,k);
	hpTOT[i][j][k] = (TH1D *)gDirectory->Get(buf);

	// convert histogram to array, identify empty histograms or bins
	double minres, maxres, nsigma=2.5;
    double deltat=1.0;
	if(iter<2) deltat = 0.3;
	else if(iter<4) deltat = 0.15;
	else deltat = 0.1;
	if(iter>0&&iter<3) {maxres=0.5; minres = 0.1;}
	if(iter>=3&&iter<4) {maxres=0.3; minres = 0.05; }
	if(iter>=4) {maxres=0.2; minres = 0.05; }	
	Int_t bincounter = 0;
	for (Int_t l=1; l<=maxNBin; l++ ) {  // l=0 means underflow
	  Double_t ytmp, eytmp;
	  ytmp = hpTOT[i][j][k]->GetBinContent(l);
	  eytmp = hpTOT[i][j][k]->GetBinError(l);
      hTOT[i][j][k]->GetXaxis()->SetRange(l,l);
	  double meanvalueproj = hTOT[i][j][k]->GetMean(2);
	  double projrms = fabs(hTOT[i][j][k]->GetRMS(2));
	  double histmeandt;
	  for(int ltmp = 0; ltmp <5; ltmp++)
	  {projrms = fabs(projrms);
	  projrms = TMath::Min(projrms,maxres);		 
	  projrms = TMath::Max(minres, projrms);	  	
	  hTOT[i][j][k]->GetYaxis()->SetRangeUser(meanvalueproj-nsigma*projrms, meanvalueproj+nsigma*projrms);	
      histmeandt = hTOT[i][j][k]->GetMean(2);
	  //cout<<"ltmp	"<<ltmp<<" dt	"<<histmeandt - meanvalueproj<<endl;
	  if(fabs(histmeandt - meanvalueproj)<1e-3) {break;}
	  meanvalueproj = histmeandt;
	  projrms = fabs(hTOT[i][j][k]->GetRMS(2));	  
	  	}
	  meanvalueproj = hTOT[i][j][k]->GetMean(2);
	  projrms = fabs(hTOT[i][j][k]->GetRMS(2));
	  
	  double meanvalueprojerr = hTOT[i][j][k]->GetMeanError(2);
	//  cout<<"tray	"<<i+1<<"   board	"<<j<<"	iter	"<<iter<<endl;
	//  cout<<" x	"<<hTOT[i][j]->GetMean(1)<<" yfit	"<<ytmp<<" projymean	"<<meanvalueproj<<"	dif	"<<ytmp-meanvalueproj<<"	fiterr	"<<eytmp<<"	projerrr	"<<meanvalueprojerr<<endl;

	  ytmp = meanvalueproj;
	  eytmp = meanvalueprojerr;
	  if(ytmp>=toftlimit) ytmp = toftlimit-0.05;
	  if(ytmp<=-toftlimit) ytmp = -toftlimit+0.05;
	  if(fabs(ytmp)>0&&fabs(eytmp)<1e-8) eytmp = 5.0;
/*
	  if(fabs(ytmp)>10||fabs(ytmp-meanvalueproj)>deltat) 
	  {ytmp = meanvalueproj;
	   eytmp = meanvalueprojerr;
	  }
*/	  
		  
	  if(!fixNBin) {
	    if(TMath::Abs(ytmp)<smallvalue*1.e-5 && TMath::Abs(eytmp)<smallvalue*1.e-5) continue;
	  	}
	  X[bincounter] =  hpTOT[i][j][k]->GetBinCenter(l);
	  Y[bincounter] = ytmp;
	  EY[bincounter] = eytmp;
	  bincounter++;
	}
	N = bincounter;
	
	cout << N << " data points" << endl;
	 hTOT[i][j][k]->GetXaxis()->UnZoom();
	 hTOT[i][j][k]->GetYaxis()->UnZoom();	
	
	if(fixNBin) {
	  if(N==maxNBin) {
	    for(Int_t l=0; l<maxNBin; l++) {
	      X[l] = tofrTOT[i][j][k][l];
	    }
	  }
	}

	Char_t title1[500], title2[500];
	sprintf(title1, "TOT Tray_%d Board_%d Cell_%d Fit", i+nTray*trayStart+1, j+TDIGStart*nBoard,k);
	sprintf(title2, "TOT Slewing Fit | Tray_%d Board_%d Cell_%d", i+nTray*trayStart+1, j+TDIGStart*nBoard,k);
	if(N<tsfM || entry<maxNBin) { 
	  tsfTOT[i][j][k] = new TSplineFit();
	} else { 
	  tsfTOT[i][j][k] = new TSplineFit(title1, title2, 51+iter*2, N, tsfM, X, Y, EY, kFALSE, tm-5.0*ts, kFALSE, tm+5.0*ts, tofrTOTBins[i][j][k][0], tofrTOTBins[i][j][k][maxNBin]); 
	}
	Bool_t ismin,ismax;
	Double_t xmintmp,ymintmp=9e9;
	Double_t xmaxtmp,ymaxtmp=-9e9;
	Double_t ymi,yma;
	tsfTOT[i][j][k]->MinMax(ismin,xmintmp,ymintmp,ismax,xmaxtmp,ymaxtmp);
	if(tsfTOT[i][j][k]->GetCategory()==-1) hnull->Draw();
	else if(ismin && ismax) {
	  ymi = ymintmp - 0.5*(ymaxtmp-ymintmp);
	  yma = ymaxtmp + 0.5*(ymaxtmp-ymintmp);
	  if(ymi<-toftlimit) ymi = -toftlimit;
	  if(yma>toftlimit) yma = toftlimit;
          if(iter==1) {
            ymi = -0.6;
            yma = 0.7;
          } else {
            ymi = -0.2;
            yma = 0.2;
          } 
	  tsfTOT[i][j][k]->DrawHere(ymi,yma);
	} else {   
                    if(iter==0) {
            ymi = -0.6;
            yma = 0.7;
          } else {
            ymi = -0.2;
            yma = 0.2;
          }
	  tsfTOT[i][j][k]->DrawHere(-toftlimit, toftlimit);
	}
      
	 if(hTOT[i][j][k]->GetEntries()>10&&fabs(tsfTOT[i][j][k]->V(tofrTOT[i][j][k][0]))<1e-8&&fabs(tsfTOT[i][j][k]->V(tofrTOT[i][j][k][int(maxNBin/2)]))<1e-8
	 	&&fabs(tsfTOT[i][j][k]->V(tofrTOT[i][j][k][maxNBin-1]))<1e-8)
	 	{
	  	tofrTOTbad[i][j][k]=1;
		cout<<"tray	"<<i+1+nTray*trayStart<<"	Board "<<j+1+TDIGStart*nBoard<<"  Cell "<<k+1<<" iter	"<<iter<<"	tot splinefit fails"<<endl;
        for(int ltmp=0; ltmp<maxNBin; ltmp++)
      	{//hTOT[i][j]->GetXaxis()->SetRange(ltmp+1,ltmp+1);
		 //tofrTOTbadcorr[i][j][ltmp] = hTOT[i][j]->GetMean(2);
		 tofrTOTbadcorr[i][j][k][ltmp] = Y[ltmp];
        }
	 	}

	if(xpisnan(tsfTOT[i][j][k]->V(tofrTOT[i][j][k][0]))||xpisnan(tsfTOT[i][j][k]->V(tofrTOT[i][j][k][int(maxNBin/2)]))
	   ||xpisnan(tsfTOT[i][j][k]->V(tofrTOT[i][j][k][maxNBin-1])))
	   {
	   tofrTOTbad[i][j][k]=1;
	   cout<<"tray "<<i+1+nTray*trayStart<<"   Board "<<j+1+TDIGStart*nBoard<<"   Cell "<<k+1<<" iter   "<<iter<<"  tot splinefit is nan"<<endl;
	   for(int ltmp=0; ltmp<maxNBin; ltmp++)
	   {
		tofrTOTbadcorr[i][j][k][ltmp] = Y[ltmp];
	   }
	   }   

    	 
    
	hTOT[i][j][k]->Reset(); // reset histogram for further use, if any
	TH1D *htemp1;
	TH1D *htemp2;
	TH1D *htemp3;
	TH1D *htemp4;
	
	sprintf(buf, "TOT Tray_%d_Board_%d_Cell_%d_0", i+nTray*trayStart+1, j+TDIGStart*nBoard,k);
	htemp1 = (TH1D *)gDirectory->FindObject(buf);
	if(htemp1) delete htemp1;
	sprintf(buf, "TOT Tray_%d_Board_%d_Cell_%d_1", i+nTray*trayStart+1, j+TDIGStart*nBoard,k);
	htemp2 = (TH1D *)gDirectory->FindObject(buf);
	if(htemp2) delete htemp2;
	sprintf(buf, "TOT Tray_%d_Board_%d_Cell_%d_2", i+nTray*trayStart+1, j+TDIGStart*nBoard,k);
	htemp3 = (TH1D *)gDirectory->FindObject(buf);
	if(htemp3) delete htemp3;
	sprintf(buf, "TOT Tray_%d_Board_%d_Cell_%d_chi2", i+nTray*trayStart+1, j+TDIGStart*nBoard,k);
	htemp4 = (TH1D *)gDirectory->FindObject(buf);
	if(htemp4) delete htemp4;
      
      c1->Update();
	}//cell
    	} //end module loop
  } // end nTray loop

  ps->Close();

  outputTOT(chain, iter);

  return 1;
}

Int_t doZ(TChain* chain, Int_t iter)
{
  cout << "Come to doZ from TOFr, " << iter << " iteration. " << endl; 

  Int_t nevt = (int)(chain->GetEntries());

  Int_t N;
  Double_t X[maxNZBin];
  Double_t Y[maxNZBin];
  Double_t EY[maxNZBin];

  char rootfile[500];
/*  sprintf(rootfile, "%s.root", pvpdOutput);
  TChain *startchain = new TChain("StartPicoDst");
  startchain->AddFile(rootfile);
*/  sprintf(rootfile, "tray_%d_board_%d_%s_%d.root",trayStart,TDIGStart, totOutput, iter);

  TChain *slewchain = new TChain("TOTPicoDst");
  slewchain->AddFile(rootfile);

  TOFrPicoDst *tofr = new TOFrPicoDst(chain);
//  StartPicoDst *start = new StartPicoDst(startchain);
  TOTPicoDst *slew = new TOTPicoDst(slewchain);

  if(nevt!=slewchain->GetEntries()) {
    cout << " Inconsistent TOFr and Slew entry ! " << endl;
    return -1;
  }

  // fill slewing histogram
  for (Int_t i=0;i<nevt;i++) {
    tofr->GetEntry(i);
  //  start->GetEntry(i);
    slew->GetEntry(i);

    if(i%(nevt/10)==0) cout << i << " events processed!" << endl;

    Int_t ieast = tofr->Ieast;
    Int_t iwest = tofr->Iwest;
    if(ieast<=0||iwest<=0) continue;
    Int_t nhits = tofr->nTofHits;

    for(Int_t j=0; j<nhits; j++) {
      Double_t dedx = tofr->dedx[j];
      Double_t p = tofr->pt[j]*TMath::CosH(tofr->eta[j]);

      if(p<0.3||p>1.) continue;
      if( p<0.6 && dedx > dedxcut) continue; // is pion?
      if(sqrt(pow(tofr->vertexX,2)+pow(tofr->vertexY,2))>rCut) continue;//for AuAu
      if(TMath::Abs(tofr->vertexZ-slew->vzVpd)>vzCut) continue;//for Auau

//      if(sqrt(pow(tofr->dcaX[j],2)+pow(tofr->dcaY[j],2))>rCut) continue;//for pprun9
//      if(TMath::Abs(tofr->dcaZ[j]-slew->vzVpd)>vzCut) continue;//for pprun9

      if(tofr->nHitsFit[j]<nHitsFitCut) continue;
      //      Double_t nSigPi = tofr->nSigPi[j];
      //      Double_t nSigK = tofr->nSigK[j];
      //      Double_t nSigP = tofr->nSigP[j];
      //      Double_t nSigE = tofr->nSigE[j];
      //      if(p<0.2 || p>1.0 || 
      //	 TMath::Abs(nSigPi)>2.0 || TMath::Abs(nSigK)<2.0 ||
      //	 TMath::Abs(nSigP)<2.0 || TMath::Abs(nSigE)<2.0) continue; // another way to choose pion

      Double_t tracklength = tofr->length[j];
      Double_t betagamma = p/piMass;
      Double_t velocity = TMath::Sqrt(1.0/(1.0/betagamma/betagamma+1.0))*c_light;
      Double_t tofPi = tracklength/velocity;

      Int_t tray = tofr->tray[j]-1-nTray*trayStart;
      Int_t module = (tofr->module[j]-1)%4;//-moduleStart*nModule;
      Int_t cell = tofr->cell[j]-1;
      Int_t tcell = module*6+cell;
      Double_t z = tofr->zLocal[j];
      Double_t tof = slew->tof[j];
    
      if(tray>=nTray || tray<0) continue;
      if(module>=nModule || module<0) continue;
      if(cell>=nCell || cell<0) continue;

      Double_t dt = tof-tofPi;
 
      Double_t invBeta = tof/tracklength*c_light;
      if(invBeta>piBetaCut->Eval(p)) continue;

      if(TMath::Abs(dt)<toftlimit) {
	hZ[tray][0][tcell]->Fill(z, dt);
      }
    }
  } // end nevt loop

  // do z calibration 
  Char_t title[500];
  sprintf(title, "tray_%d_board_%d_%s_%d.ps",trayStart,TDIGStart, zOutput, iter);
  TCanvas *c1 = new TCanvas("c1","c1",0,0,800,600);
  c1->SetFillColor(10);
  c1->SetBorderMode(0);
  c1->SetBorderSize(2);
  c1->Divide(3,3);
  TPostScript *ps = new TPostScript(title,112);

  for(Int_t i=0; i<nTray; i++) {
    ps->NewPage();
    for(Int_t j=0; j<nBoard; j++) {   //BUTTER TDIG LOOP
//		tofrZbad[i][j][k]=0;		
//      if(j%4==0) ps->NewPage();
      for(Int_t k=0; k<ntCell; k++) {

                tofrZbad[i][j][k]=0;
      if(k%4==0) ps->NewPage();


	c1->cd(2*(k%4)+1);
        hZ[i][j][k]->Draw();
        c1->cd(2*(k%4)+2);

	Int_t entry = (int)hZ[i][j][k]->GetEntries();
	Double_t tm = hZ[i][j][k]->GetMean(2); // Mean on axis Y 
	Double_t ts = hZ[i][j][k]->GetRMS(2); // RMS on axis Y
	N=0;

	hZ[i][j][k]->FitSlicesY();
	Char_t buf[500];
	sprintf(buf, "Z Tray_%d_Board_%d_Cell_%d_1", i+nTray*trayStart+1, j+TDIGStart*nBoard,k);
	hpZ[i][j][k] = (TH1D *)gDirectory->Get(buf);

	// convert histogram to array, identify empty histograms or bins
	double minres, maxres, nsigma=2.5;
    double deltat=1.0;
	if(iter<2) deltat = 0.3;
	else if(iter<4) deltat = 0.15;
	else deltat = 0.1;

	if(iter>0&&iter<3) {maxres=0.5; minres = 0.1;}
	if(iter>=3&&iter<4) {maxres=0.3; minres = 0.05; }
	if(iter>=4) {maxres=0.2; minres = 0.05; }
	
	Int_t bincounter = 0;
	for (Int_t l=1; l<=maxNZBin; l++ ) {   
	  Double_t ytmp, eytmp;
	  ytmp = hpZ[i][j][k]->GetBinContent(l);
	  eytmp = hpZ[i][j][k]->GetBinError(l);
	  hZ[i][j][k]->GetXaxis()->SetRange(l,l);
	  double meanvalueproj = hZ[i][j][k]->GetMean(2);
	  double projrms = fabs(hZ[i][j][k]->GetRMS(2));
	  double histmeandt;
	  for(int ltmp = 0; ltmp <5; ltmp++)
	  	{
	  projrms = fabs(projrms);
	  projrms = TMath::Min(projrms,maxres);		 
	  projrms = TMath::Max(minres, projrms);	  	
	  hZ[i][j][k]->GetYaxis()->SetRangeUser(meanvalueproj-nsigma*projrms, meanvalueproj+nsigma*projrms);	
      histmeandt = hZ[i][j][k]->GetMean(2);
	  //cout<<"ltmp	"<<ltmp<<" dtforz	"<<histmeandt - meanvalueproj<<endl;
	  if(fabs(histmeandt - meanvalueproj)<1e-3) {break;}
	  meanvalueproj = histmeandt;
	  projrms = fabs(hZ[i][j][k]->GetRMS(2));	  
	  	}
	  meanvalueproj = hZ[i][j][k]->GetMean(2);
	  projrms = fabs(hZ[i][j][k]->GetRMS(2));	  
	  double meanvalueprojerr = hZ[i][j][k]->GetMeanError(2);
	  //cout<<"tray	"<<i+1<<"Moduleoard	"<<j<<"	iter	"<<iter<<endl;
	 // cout<<" x	"<<hZ[i][j]->GetMean(1)<<" yfit	"<<ytmp<<" projymean	"<<meanvalueproj<<"	dif	"<<ytmp-meanvalueproj<<"	fiterr	"<<eytmp<<"	projerrr	"<<meanvalueprojerr<<endl;

	  ytmp = meanvalueproj;
	  eytmp = meanvalueprojerr;
	  if(ytmp>=toftlimit) ytmp = toftlimit-0.05;
	  if(ytmp<=-toftlimit) ytmp = -toftlimit+0.05;
	  if(fabs(ytmp)>0&&fabs(eytmp)<1e-8) eytmp = 5.0;	
/*
	  if(fabs(ytmp)>10||fabs(ytmp-meanvalueproj)>deltat) 
	  {ytmp = meanvalueproj;
	   eytmp = meanvalueprojerr;
	  }
*/	  
	  if(!fixNBin) {
	    if(TMath::Abs(ytmp)<smallvalue && TMath::Abs(eytmp)<smallvalue) continue;
	  }else{
            if(TMath::Abs(ytmp)<smallvalue && TMath::Abs(eytmp)<smallvalue) continue;
          }
	  X[bincounter] =  hpZ[i][j][k]->GetBinCenter(l);
	  Y[bincounter] = ytmp;
	  EY[bincounter] = eytmp;
	  bincounter++;
	}
	N = bincounter;
	 hZ[i][j][k]->GetXaxis()->UnZoom();
	 hZ[i][j][k]->GetYaxis()->UnZoom();		
	
	if(fixNBin) {
	  if(N==maxNZBin) {
	    for(Int_t l=0; l<maxNZBin; l++) {
	      X[l] = tofrZ[i][j][k][l];
	    }
	  }
	}

	Char_t title1[500], title2[500];
	sprintf(title1, "Z Tray_%d_Board_%d_Cell_%d Fit", i+nTray*trayStart+1, j+TDIGStart*nBoard,k);
	sprintf(title2, "Z Position Fit | Tray_%d Board_%d Cell_%d", i+nTray*trayStart+1, j+TDIGStart*nBoard,k);
	if(N<tsfM || entry<maxNZBin) { 
	  tsfZ[i][j][k] = new TSplineFit();
	} else { 
	  tsfZ[i][j][k] = new TSplineFit(title1, title2, 52+iter*2, N, tsfM, X, Y, EY, kFALSE, tm-5.0*ts, kFALSE, tm+5.0*ts, tofrZBins[i][j][k][0], tofrZBins[i][j][k][maxNZBin]); 
	}
	Bool_t ismin,ismax;
	Double_t xmintmp,ymintmp;
	Double_t xmaxtmp,ymaxtmp;
	Double_t ymi,yma;
	tsfZ[i][j][k]->MinMax(ismin,xmintmp,ymintmp,ismax,xmaxtmp,ymaxtmp);
	if(tsfZ[i][j][k]->GetCategory()==-1) hnull->Draw();
	else if(ismin && ismax) {
	  ymi = ymintmp - 0.5*(ymaxtmp-ymintmp);
	  yma = ymaxtmp + 0.5*(ymaxtmp-ymintmp);
	  if(ymi<-toftlimit) ymi = -toftlimit;
	  if(yma>toftlimit) yma = toftlimit;
          ymi = -0.2;
          yma = 0.2;
	  tsfZ[i][j][k]->DrawHere(ymi,yma);
	} else {  
          ymi = -0.2;
          yma = 0.2; 
	  tsfZ[i][j][k]->DrawHere(ymi, yma);
	}


	if(hZ[i][j][k]->GetEntries()>10&&fabs(tsfZ[i][j][k]->V(tofrZ[i][j][k][0]))<1e-8&&fabs(tsfZ[i][j][k]->V(tofrZ[i][j][k][int(maxNZBin/2)]))<1e-8
	   &&fabs(tsfZ[i][j][k]->V(tofrZ[i][j][k][maxNZBin-1]))<1e-8)
	   {
	   tofrZbad[i][j][k]=1;
	   cout<<"tray "<<i+1+nTray*trayStart<<"   Board "<<j+1+TDIGStart*nBoard<<"   Cell "<<k+1<<" iter    "<<iter<<"  Z splinefit fails"<<endl;
	   for(int ltmp=0; ltmp<maxNZBin; ltmp++)
	   {//tofrZbadcorr[i][j][ltmp] = hpZ[i][j]->GetBinContent(ltmp+1);
	    tofrZbadcorr[i][j][k][ltmp] = Y[ltmp];
		   }
	  }

	if(xpisnan(tsfZ[i][j][k]->V(tofrZ[i][j][k][0]))||xpisnan(tsfZ[i][j][k]->V(tofrZ[i][j][k][int(maxNZBin/2)]))
	   ||xpisnan(tsfZ[i][j][k]->V(tofrZ[i][j][k][maxNZBin-1])))
	   {
	   tofrZbad[i][j][k]=1;
	   cout<<"tray "<<i+1+nTray*trayStart<<"   Board "<<j+1+TDIGStart*nBoard<<"   Cell "<<k+1<<" iter   "<<iter<<"  Z splinefit is nan"<<endl;
	   for(int ltmp=0; ltmp<maxNZBin; ltmp++)
	   {
		tofrZbadcorr[i][j][k][ltmp] = Y[ltmp];
	   }
	}

	hZ[i][j][k]->Reset(); // reset histogram for further use, if any
	TH1D *htemp1;
	TH1D *htemp2;
	TH1D *htemp3;
	TH1D *htemp4;
	
	sprintf(buf, "Z Tray_%d_Board_%d_Cell_%d_0", i+nTray*trayStart+1, j+TDIGStart*nBoard,k);
	htemp1 = (TH1D *)gDirectory->FindObject(buf);
	if(htemp1) delete htemp1;
	sprintf(buf, "Z Tray_%d_Board_%d_Cell_%d_1", i+nTray*trayStart+1, j+TDIGStart*nBoard,k);
	htemp2 = (TH1D *)gDirectory->FindObject(buf);
	if(htemp2) delete htemp2;
	sprintf(buf, "Z Tray_%d_Board_%d_Cell_%d_2", i+nTray*trayStart+1, j+TDIGStart*nBoard,k);
	htemp3 = (TH1D *)gDirectory->FindObject(buf);
	if(htemp3) delete htemp3;
	sprintf(buf, "Z Tray_%d_Board_%d_Cell_%d_chi2", i+nTray*trayStart+1, j+TDIGStart*nBoard,k);
	htemp4 = (TH1D *)gDirectory->FindObject(buf);
	if(htemp4) delete htemp4;
      
      c1->Update();
     }//cell
    }//module
  } // end nTray loop

  ps->Close();

  outputZ(chain, iter);

  return 1;
}
 
Int_t doT0(TChain* chain)
{
  const Double_t sigmacut = 2.0;

  Float_t T0;
  Int_t hits;
  Int_t nPi;
  Float_t tag[maxTOFrHits];
  Float_t tofPi[maxTOFrHits];

  char rootfile[500];
  sprintf(rootfile, "%s.root", t0Output);
  TFile *forwrite = new TFile(rootfile,"RECREATE");

  TTree *t0 = new TTree("T0PicoDst","T0PicoDst");
  t0->Branch("T0", &T0, "T0/F");
  t0->Branch("hits", &hits, "hits/I");
  t0->Branch("nPi", &nPi, "nPi/I");
  t0->Branch("tag", tag, "tag[hits]/F");
  t0->Branch("tofPi", tofPi, "tofPi[hits]/F");

  TOFrPicoDst *tofr = new TOFrPicoDst(chain);

  Int_t nevt = (int)(chain->GetEntries());

  for (Int_t i=0;i<nevt;i++) {
    tofr->GetEntry(i);

    hits = tofr->nTofHits; // should >0

    Float_t t0Pi = 0.;
    nPi = 0;
    Float_t w = 0.; // normalize factor in average

    for(Int_t j=0; j<hits; j++) {
      Float_t tdc = tofr->leTime[j]; 
      Float_t p = tofr->pt[j]*TMath::CosH(tofr->eta[j]);
      Float_t tracklength = tofr->length[j];
      Float_t nsigPi = tofr->nSigPi[j];
      Float_t nsigK = tofr->nSigK[j];
      Float_t nsigP = tofr->nSigP[j];
      Float_t nsigE = tofr->nSigE[j];
        
      tag[j] = 0.;
      Float_t betagamma = p/piMass;
      Float_t velocity = TMath::Sqrt(1.0/(1.0/betagamma/betagamma+1.0))*c_light;
      tofPi[j] = tracklength/velocity;
      
      if(TMath::Abs(nsigPi)<sigmacut && TMath::Abs(nsigK)>sigmacut &&
	 TMath::Abs(nsigP)>sigmacut && TMath::Abs(nsigE)>sigmacut &&
	 p>0.2 && p<1.0) { // choose pion sample
	Float_t corTOT = 0.; // initial correction parameters ? here set to 0.0
	Float_t corZ = 0.;
	
	tag[j] = 1.0; // arithmatic average
	t0Pi += tag[j]*(tdc - tofPi[j] - corTOT - corZ);
	nPi++;
	w += tag[j];
      }
    }

    if(nPi>0) {
      T0 = t0Pi/w;
    } else {
      T0 = -1.0e6;
    }

    t0->Fill();
  }

  forwrite->cd();
  t0->Write();
  forwrite->Close();

  return 1;
}


Int_t outputPVPD(TChain* chain)
{
  // data member of output tree
  Int_t Ieast;
  Int_t Iwest;
  Double_t sumEast;
  Double_t sumWest;
  Double_t T0;
  Double_t T0vzdiff;
  Double_t vzVpd=-1e6;
//  const Float_t vzOffset = 4542.37-2.02792-4.59159+9.59631e-02-2093.13;//run 09 with run11 information.
  const Float_t vzOffset = 4542.37-2.02792-4.59159+9.59631e-02;//preliminary from initial hlt run11 data+auau19 calibration set+from august
//  const Float_t vzOffset =2446+11.4+0.9948-4.36+0.7216-0.011-0.5282-2.38469e-4;//12.29;//12.4-0.11; //11.86+0.54; //cm run10

  Double_t tot[nPVPDChannel];
  Double_t tdc[nPVPDChannel];
  Double_t cor[nPVPDChannel];

  char rootfile[500];
  sprintf(rootfile, "%s.root", pvpdOutput);
  TFile *forwrite = new TFile(rootfile,"RECREATE");

  cout << " Writing to pVPD calibration tree ... " << endl;

  TTree *pvpd = new TTree("StartPicoDst","StartPicoDst");
  pvpd->Branch("Ieast", &Ieast, "Ieast/I");
  pvpd->Branch("Iwest", &Iwest, "Iwest/I");
  pvpd->Branch("sumEast", &sumEast, "sumEast/D");
  pvpd->Branch("sumWest", &sumWest, "sumWest/D");
  pvpd->Branch("T0", &T0, "T0/D");
  pvpd->Branch("T0vzdiff", &T0vzdiff, "T0vzdiff/D");
  pvpd->Branch("vzVpd", &vzVpd, "vzVpd/D");


  TOFrPicoDst *tofr = new TOFrPicoDst(chain);
  Int_t nentries = (int)(chain->GetEntries());
  for(int i=0; i<nentries; i++) {
    tofr->GetEntry(i);
    Double_t zvtx = tofr->vertexZ;

    Ieast = 0;
    Iwest = 0;
    sumEast = 0.0;
    sumWest = 0.0;

    for(Int_t j=0; j<nPVPDChannel; j++) {
      tdc[j] = 0.0;
      tot[j] = 0.0;
      cor[j] = 0.0;
    }

    for(Int_t j=0; j<nPVPDChannel; j++) {
//	if(j==7) continue; //MASK
      tot[j] = j<nPVPDChannel/2? tofr->vpdTotWest[j]:tofr->vpdTotEast[j-nPVPDChannel/2];
      tdc[j] =  j<nPVPDChannel/2? tofr->vpdLeWest[j]:tofr->vpdLeEast[j-nPVPDChannel/2];
      //tot[j] /= 1000.;
      //tdc[j] /= 1000.;
      if(tot[j]>minTOT_pVPD && tot[j]<maxTOT_pVPD) {
        cor[j] = tsfPVPD[j]->V(tot[j]);
        if(j<nPVPDChannel/2) {
          /*      Int_t itmp = 0;
            for(Int_t k=0; k<nPVPDChannel/2; k++) {
            if(k==j) continue;
            if(TMath::Abs(tofr->pvpdLeadingEdgeTimeWest[k]-tdc[j])>10) itmp++;
            }
            if(itmp>1) continue;
          */
          sumWest += tdc[j] - cor[j];
          Iwest++;
        } else {
          /*      Int_t itmp = 0;
            for(Int_t k=nPVPDChannel/2; k<nPVPDChannel; k++) {
            if(k==j) continue;
            if(TMath::Abs(tofr->pvpdLeadingEdgeTimeEast[k]-tdc[j])>10) itmp++;
            }
            if(itmp>1) continue;
          */
          sumEast += tdc[j] - cor[j];
          Ieast++;
        }
      }
    }


    // begin outlier rejection - rafael
    int flag[nPVPDChannel];
    for(int i=0;i<nPVPDChannel;i++){
		 flag[i] = 1; // all hits are accepted at the beginning

//		if(i==7) flag[i]=0;//MASK

		}
    Double_t vpdWHit[nPVPDChannel/2] = {0};
    Double_t vpdEHit[nPVPDChannel/2] = {0};
    for(Int_t j=0; j<nPVPDChannel; j++) {
//	if(j==7) continue; //MASK
      if(tot[j]<=minTOT_pVPD || tot[j]>=maxTOT_pVPD) continue;
      if(j<nPVPDChannel/2) vpdWHit[j] = tdc[j]-cor[j];
      else vpdEHit[j-nPVPDChannel/2] = tdc[j]-cor[j];
    }
    // west
    int countWHits = 0;
    Int_t *hitIndexW = new Int_t[nPVPDChannel/2];
    TMath::Sort(nPVPDChannel/2,vpdWHit,hitIndexW);
    for(int i=0;i<nPVPDChannel/2;i++) if(vpdWHit[i]!=0) countWHits++;
    int nWRejected = (int)(pCut*countWHits+0.5);
    for(int i=0;i<nWRejected;i++) flag[hitIndexW[i]] = 0; // hit rejected
    delete hitIndexW;
    // east
    int countEHits = 0;
    Int_t *hitIndexE = new Int_t[nPVPDChannel/2];
    TMath::Sort(nPVPDChannel/2,vpdEHit,hitIndexE);
    for(int i=0;i<nPVPDChannel/2;i++) if(vpdEHit[i]!=0) countEHits++;
    int nERejected = (int)(pCut*countEHits+0.5);
    for(int i=0;i<nERejected;i++) flag[hitIndexE[i]+nPVPDChannel/2] = 0; // hit rejected
    delete hitIndexE;
    //     cout << endl;
    //     cout << "number of hits: " << countWHits << "\t hits rejected: " << nWRejected << endl;
    //     for(int i=0;i<nPVPDChannel/2;i++) {
    //       cout << vpdWHit[i];
    //       if(flag[i]==0) cout << " --> rejected!" << endl;
    //       else cout << endl;
    //     }
    // end outlier rejection - rafael


    if(i<10) cout<<"Iwest/Ieast = "<<Iwest<<"/"<<Ieast<<endl;
    Double_t vpdtime;
    //---- reject hits from other collisions
    for(int j=0;j<nPVPDChannel;++j){

//	if(j==7) continue;//MASK

      if(tot[j]<=minTOT_pVPD ||tot[j]>=maxTOT_pVPD) continue;
      if(j<nPVPDChannel/2&&Iwest>1) {
        vpdtime = ((tdc[j]-cor[j])*Iwest-sumWest)/(Iwest-1);
        //if(i<10) cout<<tdc[j]<<"  "<<cor[j]<<" vpdtime = "<<vpdtime<<endl;
        //      if(fabs(vpdtime)>tdiffCut) { // commented .. rafael
        if(fabs(vpdtime)>tdiffCut || flag[j]==0) { // added .. rafael
          sumWest -= tdc[j]-cor[j];
          Iwest--;
        }
      }
      if(j>=nPVPDChannel/2&&Ieast>1) {
        vpdtime = ((tdc[j]-cor[j])*Ieast-sumEast)/(Ieast-1);
        //if(i<10) cout<<"vpdtime="<<vpdtime<<endl;
        //      if(fabs(vpdtime)>tdiffCut) { // commented .. rafael
        if(fabs(vpdtime)>tdiffCut || flag[j]==0) { // added .. rafael
          sumEast -= tdc[j]-cor[j];
          Ieast--;
        }
      }
    }//rejected
    //if(i<10) cout<<"Iwest/Ieast2 = "<<Iwest<<"/"<<Ieast<<endl;

    sumEast += 2*vzOffset*Ieast/c_light;//correct global offset between two sides

    if(Ieast>0||Iwest>0) {
      //T0 = (sumEast+sumWest+(Iwest-Ieast)*zvtx/c_light)/(Ieast+Iwest);
      //float teast = sumEast;
      //while(teast>51200) teast -= 51200;
      //teast += 2*vzOffset/c_light;
      float tdiff;
      if(Ieast>0&&Iwest>0) {
//	cout<<"sE: "<<sumEast<<" iE: "<<Ieast<<" sW: "<<sumWest<<" iW: "<<Iwest<<endl;
        tdiff = sumEast/Ieast - sumWest/Iwest;
        vzVpd = tdiff/2*c_light;
      } else {
        vzVpd = -1e6;
      }
      /* float pttmp=0.;
         int itmp=-1;
         for(int k=0;k<tofr->nTofHits;k++) {
         if(tofr->pt[k]>pttmp) {//select track with highest pt
         pttmp = (tofr->pt[k]);
         itmp = k;
         }
         if(sqrt(pow(tofr->dcaX[itmp],2)+pow(tofr->dcaY[itmp], 2))<rCut) {//see correlation
         hVzCorr->Fill(vzVpd, tofr->dcaZ[itmp]);
         hVzCorr2->Fill(vzVpd, tofr->dcaZ[itmp]- vzVpd);
         }
         }
      */
/*      float dcaRmin = 9999;
      Int_t itmp = -1;
      for(int k=0;k<tofr->nTofHits;++k) {
      //    cout<<"k = "<<k<<"  "<<tofr->dcaX[k]<<"  "<<tofr->dcaY[k]<<"  nhits: "<< tofr->nHitsFit[k]<<"   dcaz: "<<tofr->dcaZ[k]<<"  vzvpd: "<< vzVpd<<endl;

        if(dcaRmin>sqrt(pow(tofr->dcaX[k],2)+pow(tofr->dcaY[k],2))&&tofr->nHitsFit[k]>=nHitsFitCut){
//          if(dcaRmin>sqrt(pow(tofr->dcaX[k],2)+pow(tofr->dcaY[k],2))&&TMath::Abs(tofr->dcaZ[k]-vzVpd)<=vzCut&&tofr->nHitsFit[k]>=nHitsFitCut){
          dcaRmin = sqrt(pow(tofr->dcaX[k],2)+pow(tofr->dcaY[k],2));
          itmp = k;
        //  cout<<"after cut k = "<<k<<"  "<<tofr->dcaX[k]<<"  "<<tofr->dcaY[k]<<"  "<<endl;
        }
          //      cout<<"dcaRmin = "<<dcaRmin<<endl;

      } //k loop
*/
      /*if(itmp>-1){
        if(sqrt(pow(tofr->dcaX[itmp],2)+pow(tofr->dcaY[itmp], 2))<rCut) {//see correlation
        ////cout<<"fill correlation"<<endl;
        hVzCorr->Fill(vzVpd, tofr->dcaZ[itmp]);
        hVzCorr2->Fill(vzVpd, vzVpd - tofr->dcaZ[itmp]);
        }
        }*/ //used for dAu and pp

//	if(itmp>-1){//added
//                  if(sqrt(pow(tofr->dcaX[itmp],2)+pow(tofr->dcaY[itmp], 2))<rCut) {//see correlation

      hVzCorr->Fill(vzVpd, zvtx);
      hVzCorr2->Fill(vzVpd, vzVpd-zvtx);//for AuAu
//      hVzCorrexp->Fill(vzVpd, zvtx);
//      hVzCorr2exp->Fill(vzVpd, vzVpd-zvtx);//for AuAu

	if(fabs(vzVpd-zvtx)<=vzCut){
      hVzCorrvzdiff->Fill(vzVpd, zvtx);
      hVzCorr2vzdiff->Fill(vzVpd, vzVpd-zvtx);//for AuAu
//      hVzCorrexpvzdiff->Fill(vzVpd, zvtx);
//      hVzCorr2expvzdiff->Fill(vzVpd, vzVpd-zvtx);//for AuAu
	}

//                  }
//        }
	Double_t vtxx = tofr->vertexX;
	Double_t vtxy = tofr->vertexY;
	hvxvvy->Fill(vtxx,vtxy);
	if(fabs(vzVpd-zvtx)<=vzCut){
		hvxvvyvzdiff->Fill(vtxx,vtxy);
	}

      //      if(tofr->nTofHits>0 && itmp>-1) T0 = (sumEast+sumWest+(Iwest-Ieast)*tofr->dcaZ[itmp]/c_light)/(Ieast+Iwest);
      if(tofr->nTofHits>0) T0 = (sumEast+sumWest+(Iwest-Ieast)*tofr->vertexZ/c_light)/(Ieast+Iwest); // added .
      else T0 = -1e6;
	if(tofr->nTofHits>0&&(fabs(vzVpd-zvtx)<=vzCut)) T0vzdiff = (sumEast+sumWest+(Iwest-Ieast)*tofr->vertexZ/c_light)/(Ieast+Iwest); // added .
      else T0vzdiff = -1e6;
      //cout<<" dcaZ[" <<itmp<<"] = "<<tofr->dcaZ[itmp]<<"  "<<T0<<"  "<<sumEast<<"  "<<sumWest<<"  "<<Iwest<<"  "<<Ieast<

      if(Ieast==1&&Iwest==1) {
        //        hResoPVPD11->Fill((sumEast-sumWest)/2.-tofr->dcaZ[itmp]/c_light); // commented .. rafael
        hResoPVPD11->Fill((sumEast-sumWest)/2.-tofr->vertexZ/c_light); // added .. rafael
      }
      if(Ieast==2&&Iwest==2) {
        //        hResoPVPD22->Fill((sumEast-sumWest)/4.-tofr->dcaZ[itmp]/c_light); // commented .. rafael
        hResoPVPD22->Fill((sumEast-sumWest)/4.-tofr->vertexZ/c_light); // added .. rafael
      }
      if(Ieast==3&&Iwest==3) {
        //        hResoPVPD33->Fill((sumEast-sumWest)/6.-tofr->dcaZ[itmp]/c_light); // commented .. rafael
        hResoPVPD33->Fill((sumEast-sumWest)/6.-tofr->vertexZ/c_light); // added .. rafael
      }
    } else {
      T0 = -1.0e6;
      T0vzdiff = -1.0e6;
      vzVpd = -1.0e6;
    }
    if(i<10) cout<<"T0/vzVpd/Iwest/Ieast/sumWest/sumEast = "<<T0<<"/"<<vzVpd<<"/"<<Iwest<<"/"<<Ieast<<"/"<<sumWest<<"/"<<sumEast<<endl;

    pvpd->Fill();
  } // end nevent loop

  forwrite->cd();
  pvpd->Write();
  hResoPVPD11->Write();
  hResoPVPD22->Write();
  hResoPVPD33->Write();
  hVzCorr->Write();
  hVzCorr2->Write();
//  hVzCorrexp->Write();
//  hVzCorr2exp->Write();
  hVzCorrvzdiff->Write();
  hVzCorr2vzdiff->Write();
//  hVzCorrexpvzdiff->Write();
//  hVzCorr2expvzdiff->Write();
  hvxvvy->Write();
  hvxvvyvzdiff->Write();
  hdeltatwest->Write();
  hdeltateast->Write();
  forwrite->Close();

  // begin check file - rafael
  TFile *f = new TFile("pvpdCheck.root","UPDATE");
  for(int i=0;i<nPVPDChannel;i++) {
    hPVPDMeanvsIt[i]->Write();
    hPVPDFitMeanvsIt[i]->Write();
    hPVPDFitSigmavsIt[i]->Write();
    hPVPDResovsIt[i]->Write();
    //    hPVPD[i]->Write();
    //    hPVPDCorr[i]->Write();
  }
  f->Close();
  delete f;
  // end check file - rafael

  cout << " Writing out pVPD calibration parameters ... " << endl;

  char datafile[500];
  sprintf(datafile, "%s.dat", pvpdOutput);
  ofstream outData;
  outData.open(datafile);

  for(Int_t i=0; i<nPVPDChannel; i++) {
    outData << setw(4) << i << endl;
    if(fixNBin) {
      for(Int_t j=0; j<=maxpvpdNBin; j++) {
        outData << setw(16) << pvpdTOTBins[i][j];
      }
      outData << endl;
      for(Int_t j=0; j<=maxpvpdNBin; j++) {
        outData << setw(16) << tsfPVPD[i]->V(pvpdTOTBins[i][j]);
      }
      outData << endl;

      for(Int_t j=0; j<maxpvpdNBin; j++) {
        outData << setw(16) << pvpdTOTX[i][j];
      }
      outData << endl;

      for(Int_t j=0; j<maxpvpdNBin; j++) {
        outData << setw(16) << tsfPVPD[i]->V(pvpdTOTX[i][j]);
      }
      outData << endl;
    }

    if(fixBinWidth) {
      Double_t step = (maxTOT_pVPD-minTOT_pVPD)/maxpvpdNBin;
      for(Int_t j=0; j<=maxpvpdNBin; j++) {
        outData << setw(16) << minTOT_pVPD+step*j;
      }
      outData << endl;
      for(Int_t j=0; j<=maxpvpdNBin; j++) {
        outData << setw(16) << tsfPVPD[i]->V(minTOT_pVPD+step*j);
      }
      outData << endl;
    }
  }
  outData.close();

  return 1;
}



Int_t outputparPVPD(TChain * chain)
{
  // data member of output tree
  Float_t pvpdTotEdge[nPVPDChannel][maxpvpdNBin+1];
  Float_t pvpdTotCorr[nPVPDChannel][maxpvpdNBin+1];
  ifstream indata;
  indata.open("./pvpdCali_4DB.dat");
  int nchl, nbin;
  for(int i=0;i<nPVPDChannel;i++){
    indata>>nchl;
    indata>>nbin;
        cout<<nchl<<"   "<<nbin<<endl;
    for(int j=0;j<=maxpvpdNBin;j++){
      indata>>pvpdTotEdge[i][j];
        //  cout<<
    }
    for(int j=0;j<=maxpvpdNBin;j++){
      indata>>pvpdTotCorr[i][j];
    }
  }
  indata.close();

    for(int i=0;i<nPVPDChannel;i++){
        cout<<pvpdTotEdge[i][0]<<"      "<<pvpdTotEdge[i][maxpvpdNBin]<<endl;
        cout<<pvpdTotCorr[i][0]<<"      "<<pvpdTotCorr[i][maxpvpdNBin]<<endl;
    }

  cout<<"read pvpd parameters"<<endl;

  Int_t Ieast;
  Int_t Iwest;
  Double_t sumEast;
  Double_t sumWest;
  Double_t T0;
  Double_t vzVpd=-1e6;
//  const Float_t vzOffset = 2446+11.4+0.9948-4.36+0.7216-0.011-0.5282-2.38469e-4;//12.29;//12.4-0.11; //11.86+0.54; //cm
  const Float_t vzOffset = 0;

  Double_t tot[nPVPDChannel];
  Double_t tdc[nPVPDChannel];
  Double_t cor[nPVPDChannel];
  Double_t mVPDLeTime[nPVPDChannel];//added
  char rootfile[500];
  sprintf(rootfile, "%s.root", pvpdOutput);
  TFile *forwrite = new TFile(rootfile,"RECREATE");

  cout << " Writing to pVPD calibration tree ... " << endl;

  TTree *pvpd = new TTree("StartPicoDst","StartPicoDst");
  pvpd->Branch("Ieast", &Ieast, "Ieast/I");
  pvpd->Branch("Iwest", &Iwest, "Iwest/I");
  pvpd->Branch("sumEast", &sumEast, "sumEast/D");
  pvpd->Branch("sumWest", &sumWest, "sumWest/D");
  pvpd->Branch("T0", &T0, "T0/D");
  pvpd->Branch("vzVpd", &vzVpd, "vzVpd/D");

   Int_t nentries = (int)(chain->GetEntries());
  cout<<" entries       "<<nentries<<endl;
  TOFrPicoDst *tofr = new TOFrPicoDst(chain);

  for(int i=0; i<nentries; i++) {
    tofr->GetEntry(i);

    //Double_t zvtx = tofr->vertexZ;

    Ieast = 0;
    Iwest = 0;
    sumEast = 0.0;
    sumWest = 0.0;

    for(Int_t j=0; j<nPVPDChannel; j++) {
      tdc[j] = 0.0;
      tot[j] = 0.0;
      cor[j] = 0.0;
      mVPDLeTime[j] = 0.0;
    }

    for(Int_t j=0; j<nPVPDChannel; j++) {

//	if(j==7) continue;//MASK

      tot[j] = j<nPVPDChannel/2? tofr->vpdTotWest[j]:tofr->vpdTotEast[j-nPVPDChannel/2];
      tdc[j] =  j<nPVPDChannel/2? tofr->vpdLeWest[j]:tofr->vpdLeEast[j-nPVPDChannel/2];

         //cout<<"tot   "<<tot[j]<<"    tdc     "<<tdc[j]<<endl;
      //tot[j] /= 1000.;
      //tdc[j] /= 1000.;
    //  if(tot[j]>minTOT_pVPD && tot[j]<maxTOT_pVPD) {
          if(tot[j]>=pvpdTotEdge[j][0] && tot[j]<pvpdTotEdge[j][maxpvpdNBin]) {
        //cor[j] = tsfPVPD[j]->V(tot[j]);
       int ibin = -1;
       for(int k=0;k<maxpvpdNBin;k++){
       if(tot[j]>=pvpdTotEdge[j][k]&&tot[j]<pvpdTotEdge[j][k+1]){
       ibin = k;
       break;}
       }
           if(ibin<0) cout<<"error in getting vpd calibration parameters!  tot =        "<<tot[j]<<"    channel "<<j+1<<endl;
       cor[j] = linearInter(pvpdTotEdge[j][ibin],pvpdTotEdge[j][ibin+1],
                             pvpdTotCorr[j][ibin],pvpdTotCorr[j][ibin+1], tot[j]);
        mVPDLeTime[j]=tdc[j]-cor[j];//added

        if(j<nPVPDChannel/2) {
/*        Int_t itmp = 0;
          for(Int_t k=0; k<nPVPDChannel/2; k++) {
            if(k==j) continue;
            if(TMath::Abs(tofr->pvpdLeadingEdgeTimeWest[k]-tdc[j])>10) itmp++;
          }
          if(itmp>1) continue;
*/
          sumWest += tdc[j] - cor[j];
          Iwest++;
        } else {
/*        Int_t itmp = 0;
          for(Int_t k=nPVPDChannel/2; k<nPVPDChannel; k++) {
            if(k==j) continue;
            if(TMath::Abs(tofr->pvpdLeadingEdgeTimeEast[k]-tdc[j])>10) itmp++;
          }
          if(itmp>1) continue;
*/
          sumEast += tdc[j] - cor[j];
          Ieast++;
        }
                } // tot check
    } //channel loop

    if(i<10) cout<<"Iwest/Ieast = "<<Iwest<<"/"<<Ieast<<endl;
    Double_t vpdtime;
    //---- reject hits from other collisions
    for(int j=0;j<nPVPDChannel;++j){
//	if(j==7) continue;//MASK
      if(tot[j]<pvpdTotEdge[j][0] ||tot[j]>=pvpdTotEdge[j][maxpvpdNBin]) continue;
      if(j<nPVPDChannel/2&&Iwest>1) {
        vpdtime = (( mVPDLeTime[j])*Iwest-sumWest)/(Iwest-1);
        //if(i<10) cout<<tdc[j]<<"  "<<cor[j]<<" vpdtime = "<<vpdtime<<endl;
        if(fabs(vpdtime)>tdiffCut) {//addedRafael outlierflag
          sumWest -=  mVPDLeTime[j];
           mVPDLeTime[j]=0.;
          Iwest--;
        }
      }
      if(j>=nPVPDChannel/2&&Ieast>1) {
        vpdtime = ( (mVPDLeTime[j])*Ieast-sumEast)/(Ieast-1);
        //if(i<10) cout<<"vpdtime="<<vpdtime<<endl;
        if(fabs(vpdtime)>tdiffCut) {//added Rafael outlierflag
          sumEast -=  mVPDLeTime[j];
           mVPDLeTime[j] =0.;
          Ieast--;
        }
      }
    }//rejected
    Int_t hitIndex[nPVPDChannel];
    Int_t nTube = nPVPDChannel/2;
    TMath::Sort(nTube, &mVPDLeTime[0], &hitIndex[0]);
    int nRejectedWest = (int)(pCut*Iwest+0.5);
//    cout << " NWest before = " << Iwest << " rejected = " << nRejectedWest << endl;
    for(int ii=0;ii<nRejectedWest;ii++) {
      int index = hitIndex[ii];
/*	if(index==7){
			cout<<"JOEY ERROR?"<<endl;
			continue;
			}
	if(index==7)cout<<"JOEY ERROR?lolcontinue"<<endl;
*/
      sumWest -= mVPDLeTime[index];
      mVPDLeTime[index] = 0.;
      Iwest--;
//      mVPDHitPatternWest &= ( 0x7ffff - (1<<index) );
//      mFlag[index] = 0;
    }

    TMath::Sort(nTube, &mVPDLeTime[nPVPDChannel/2], &hitIndex[nPVPDChannel/2]);
    int nRejectedEast = (int)(pCut*Ieast+0.5);
//   cout << " NEast before = " << Ieast << " rejected = " << nRejectedEast << endl;
    for(int ii=0;ii<nRejectedEast;ii++) {
      int index = hitIndex[ii+nPVPDChannel/2] + nPVPDChannel/2;
      sumEast -= mVPDLeTime[index];
      mVPDLeTime[index] = 0.;
      Ieast--;
//      mVPDHitPatternEast &= ( 0x7ffff - (1<<(index-NVPD)) );
//      mFlag[index] = 0;
    }




/*    // begin outlier rejection - rafael
    int flag[nPVPDChannel];
    for(int i=0;i<nPVPDChannel;i++) flag[i] = 1; // all hits are accepted at the beginning
    Double_t vpdWHit[nPVPDChannel/2] = {0};
    Double_t vpdEHit[nPVPDChannel/2] = {0};
    for(Int_t j=0; j<nPVPDChannel; j++) {
      if(tot[j]<=minTOT_pVPD || tot[j]>=maxTOT_pVPD) continue;
      if(j<nPVPDChannel/2) vpdWHit[j] = tdc[j]-cor[j];
      else vpdEHit[j-nPVPDChannel/2] = tdc[j]-cor[j];
    }
    // west
    int countWHits = 0;
    Int_t *hitIndexW = new Int_t[nPVPDChannel/2];
    TMath::Sort(nPVPDChannel/2,vpdWHit,hitIndexW);
    for(int i=0;i<nPVPDChannel/2;i++) if(vpdWHit[i]!=0) countWHits++;
    int nWRejected = (int)(pCut*countWHits+0.5);
    for(int ii=0;ii<nWRejected;ii++){
         flag[hitIndexW[ii]] = 0; // hit rejected
         sumWest -= vpdWHit[ii];
         Iwest--;
    }
    delete hitIndexW;
    // east
    int countEHits = 0;
    Int_t *hitIndexE = new Int_t[nPVPDChannel/2];
    TMath::Sort(nPVPDChannel/2,vpdEHit,hitIndexE);
    for(int i=0;i<nPVPDChannel/2;i++) if(vpdEHit[i]!=0) countEHits++;
    int nERejected = (int)(pCut*countEHits+0.5);
    for(int ii=0;ii<nERejected;ii++){
         flag[hitIndexE[ii]+nPVPDChannel/2] = 0; // hit rejected
         sumEast -= vpdEHit[ii];
         Ieast--;
    }
    delete hitIndexE;
    //     cout << endl;
    //     cout << "number of hits: " << countWHits << "\t hits rejected: " << nWRejected << endl;
    //     for(int i=0;i<nPVPDChannel/2;i++) {
    //       cout << vpdWHit[i];
    //       if(flag[i]==0) cout << " --> rejected!" << endl;
    //       else cout << endl;
    //     }
    // end outlier rejection - rafael

*/

/*    if(i<10) cout<<"Iwest/Ieast = "<<Iwest<<"/"<<Ieast<<endl;
    Double_t vpdtime;
    //---- reject hits from other collisions
    for(int j=0;j<nPVPDChannel;++j){
      if(tot[j]<pvpdTotEdge[j][0] ||tot[j]>=pvpdTotEdge[j][maxpvpdNBin]) continue;
      if(j<nPVPDChannel/2&&Iwest>1) {
        vpdtime = ((tdc[j]-cor[j])*Iwest-sumWest)/(Iwest-1);
        //if(i<10) cout<<tdc[j]<<"  "<<cor[j]<<" vpdtime = "<<vpdtime<<endl;
        if(fabs(vpdtime)>tdiffCut||flag[j]==0) {//addedRafael outlierflag
          sumWest -= tdc[j]-cor[j];
          Iwest--;
        }
      }
      if(j>=nPVPDChannel/2&&Ieast>1) {
        vpdtime = ((tdc[j]-cor[j])*Ieast-sumEast)/(Ieast-1);
        //if(i<10) cout<<"vpdtime="<<vpdtime<<endl;
        if(fabs(vpdtime)>tdiffCut||flag[j]==0) {//added Rafael outlierflag
          sumEast -= tdc[j]-cor[j];
          Ieast--;
        }
      }
    }//rejected
    //if(i<10) cout<<"Iwest/Ieast2 = "<<Iwest<<"/"<<Ieast<<endl;
*/
    sumEast += 2*vzOffset*Ieast/c_light;//correct global offset between two sides

    if(Ieast>0||Iwest>0) {
      //T0 = (sumEast+sumWest+(Iwest-Ieast)*zvtx/c_light)/(Ieast+Iwest);
      //float teast = sumEast;
      //while(teast>51200) teast -= 51200;
      //teast += 2*vzOffset/c_light;
      float tdiff;
     if(Ieast>0&&Iwest>0) {
      tdiff = sumEast/Ieast - sumWest/Iwest;
      vzVpd = tdiff/2*c_light;
     } else {
       vzVpd = -1e6;
     }
     /* float pttmp=0.;
      int itmp=-1;
      for(int k=0;k<tofr->nTofHits;k++) {
        if(tofr->pt[k]>pttmp) {//select track with highest pt
           pttmp = (tofr->pt[k]);
           itmp = k;
        }
       if(sqrt(pow(tofr->dcaX[itmp],2)+pow(tofr->dcaY[itmp], 2))<rCut) {//see correlation
          hVzCorr->Fill(vzVpd, tofr->vertexZ);
          hVzCorr2->Fill(vzVpd, tofr->vertexZ- vzVpd);
        }
      }
      */
      float dcaRmin = 9999;
      Int_t itmp = -1;
      for(int k=0;k<tofr->nTofHits;++k) {
        if(dcaRmin>sqrt(pow(tofr->dcaX[k],2)+pow(tofr->dcaY[k],2))&&TMath::Abs(tofr->dcaZ[k]-vzVpd)<=vzCut&&tofr->nHitsFit[k]>=nHitsFitCut){
          dcaRmin = sqrt(pow(tofr->dcaX[k],2)+pow(tofr->dcaY[k],2));
          itmp = k;
          //cout<<"k = "<<k<<"  "<<tofr->dcaX[k]<<"  "<<tofr->dcaY[k]<<"  "<<endl;
        }
                //cout<<"dcaRmin = "<<dcaRmin<<endl;

        } //k loop
//       if(itmp>-1){
//                  if(sqrt(pow(tofr->dcaX[itmp],2)+pow(tofr->dcaY[itmp], 2))<rCut) {//see correlation
          ////cout<<"fill correlation"<<endl;
          hVzCorr->Fill(vzVpd, tofr->vertexZ);
          hVzCorr2->Fill(vzVpd, vzVpd - tofr->vertexZ);
	 if(fabs(vzVpd-tofr->vertexZ)<=vzCut){
	hVzCorrvzdiff->Fill(vzVpd, tofr->vertexZ);
          hVzCorr2vzdiff->Fill(vzVpd, vzVpd - tofr->vertexZ);
	}
//                  }
//        }

//      if(tofr->nTofHits>0 && itmp>-1) T0 = (sumEast+sumWest+(Iwest-Ieast)*tofr->vertexZ/c_light)/(Ieast+Iwest);
	if(tofr->nTofHits>0 ) T0 = (sumEast+sumWest+(Iwest-Ieast)*tofr->vertexZ/c_light)/(Ieast+Iwest);
      else T0 = -1e6;
      //cout<<" dcaZ[" <<itmp<<"] = "<<tofr->dcaZ[itmp]<<"  "<<T0<<"  "<<sumEast<<"  "<<sumWest<<"  "<<Iwest<<"  "<<Ieast<<endl;

      if(Ieast==1&&Iwest==1) {
        hResoPVPD11->Fill((sumEast-sumWest)/2.-tofr->vertexZ/c_light);
      }
      if(Ieast==2&&Iwest==2) {
        hResoPVPD22->Fill((sumEast-sumWest)/4.-tofr->vertexZ/c_light);
      }
      if(Ieast==3&&Iwest==3) {
        hResoPVPD33->Fill((sumEast-sumWest)/6.-tofr->vertexZ/c_light);
      }
    } else {
      T0 = -1.0e6;
      vzVpd = -1.0e6;
    }
    if(i<10) cout<<"T0/vzVpd/Iwest/Ieast/sumWest/sumEast = "<<T0<<"/"<<vzVpd<<"/"<<Iwest<<"/"<<Ieast<<"/"<<sumWest<<"/"<<sumEast<<endl;

    pvpd->Fill();
  } // end nevent loop

  forwrite->cd();
  pvpd->Write();
  hResoPVPD11->Write();
  hResoPVPD22->Write();
  hResoPVPD33->Write();
  hVzCorr->Write();
  hVzCorr2->Write();
  hVzCorrvzdiff->Write();
  hVzCorr2vzdiff->Write();
  hdeltatwest->Write();
  hdeltateast->Write();
  forwrite->Close();

/*
  cout << " Writing out pVPD calibration parameters ... " << endl;

  char datafile[500];
  sprintf(datafile, "%s.dat", pvpdOutput);
  ofstream outData;
  outData.open(datafile);

  for(Int_t i=0; i<nPVPDChannel; i++) {
    outData << setw(4) << i << endl;
    if(fixNBin) {
      for(Int_t j=0; j<=maxNBin; j++) {
        outData << setw(16) << pvpdTOTBins[i][j];
      }
      outData << endl;
      for(Int_t j=0; j<=maxNBin; j++) {
        outData << setw(16) << tsfPVPD[i]->V(pvpdTOTBins[i][j]);
      }
      outData << endl;

      for(Int_t j=0; j<maxNBin; j++) {
        outData << setw(16) << pvpdTOTX[i][j];
      }
      outData << endl;

        for(Int_t j=0; j<maxNBin; j++) {
          outData << setw(16) << tsfPVPD[i]->V(pvpdTOTX[i][j]);
        }
        outData << endl;
        }

    if(fixBinWidth) {
      Double_t step = (maxTOT_pVPD-minTOT_pVPD)/maxNBin;
      for(Int_t j=0; j<=maxNBin; j++) {
        outData << setw(16) << minTOT_pVPD+step*j;
      }
      outData << endl;
      for(Int_t j=0; j<=maxNBin; j++) {
        outData << setw(16) << tsfPVPD[i]->V(minTOT_pVPD+step*j);
      }
      outData << endl;
    }
  }
  outData.close();
*/
  return 1;
}



Int_t outputTOT(TChain* chain, Int_t iter)
{
  Int_t hits;
  Float_t tof[maxTOFrHits];
  Float_t tag[maxTOFrHits];
  Float_t vzVpd;

  TOFrPicoDst *tofr = new TOFrPicoDst(chain);
  Int_t nentries = (int)(chain->GetEntries());

  char rootfile[500];
  sprintf(rootfile, "tray_%d_board_%d_%s_%d.root",trayStart,TDIGStart, totOutput, iter);
  TFile *forwrite = new TFile(rootfile,"RECREATE");
  
  cout << " Writing to TOT slewing correction tree ... " << endl;
  
  TTree *slew = new TTree("TOTPicoDst","TOTPicoDst");
  slew->Branch("hits", &hits, "hits/I");
  slew->Branch("tof", tof, "tof[hits]/F");
  slew->Branch("vzVpd", &vzVpd, "vzVpd/F");
  slew->Branch("tag", tag, "tag[hits]/F");

  sprintf(rootfile, "tray_%d_board_%d_%s_%d.root",trayStart, TDIGStart,delayOutput, iter-1);
  TChain *delaychain = new TChain("DelayPicoDst");
  delaychain->AddFile(rootfile);
  DelayPicoDst *delay = new DelayPicoDst(delaychain);
  
  for(int i=0; i<nentries; i++) {
    tofr->GetEntry(i);
    delay->GetEntry(i);

    hits = tofr->nTofHits;
    vzVpd = delay->vzVpd;
 
    for(Int_t j=0; j<hits; j++) {
      tag[j] = 0;
      tof[j] = -1.0e6;
      
      Int_t tray = tofr->tray[j]-1-nTray*trayStart;
      Int_t module = (tofr->module[j]-1)%4;//-moduleStart*nModule;
      Int_t cell = tofr->cell[j]-1;
      Double_t tot = tofr->tot[j];
      Int_t tcell = module*6+cell;
      
      if(tray>=nTray || tray<0) continue;
      if(module>=nModule || module<0) continue;
      if(cell>=nCell || cell<0)continue;
      
      Double_t oldtof = delay->tof[j] ;
	  Double_t tmptot=0;
	  Int_t kindex = 0;
      
      tag[j] = 1;
      if(tsfTOT[tray][0][tcell]->GetCategory()==-1) {
      } else {
         if(tot>=tofrTOT[tray][0][tcell][0]&&tot<=tofrTOT[tray][0][tcell][maxNBin-1])
         {
            if(tofrTOTbad[tray][0][tcell]<1)
            {tof[j] = oldtof-tsfTOT[tray][0][tcell]->V(tot);}
		    else{
                  int ibin = -1;
                  for(int k=0;k<=maxNBin;k++){
                    if(tot>=tofrTOT[tray][0][tcell][k]&&tot<tofrTOT[tray][0][tcell][k+1]){
                	ibin = k;
                	break;
                	} //if 
            	  } //for loop
              float tmpcor = linearInter(tofrTOT[tray][0][tcell][ibin],tofrTOT[tray][0][tcell][ibin+1],
              tofrTOTbadcorr[tray][0][tcell][ibin],tofrTOTbadcorr[tray][0][tcell][ibin+1], tot);
              tof[j] = oldtof - tmpcor;
		    } //else fit not good
         } //if tot in boundary
		 else{
		 	   if(tot<tofrTOT[tray][0][tcell][0]) 
		 		{tmptot = tofrTOT[tray][0][tcell][0];
			     kindex = 0;}
			    else
				{tmptot = tofrTOT[tray][0][tcell][maxNBin-1];
			     kindex = maxNBin-1;}
		
		 		if(tofrTOTbad[tray][0][tcell]<1)
            	{tof[j] = oldtof - tsfTOT[tray][0][tcell]->V(tmptot);}
				else {
				tof[j] = oldtof - tofrTOTbadcorr[tray][0][tcell][kindex];
				}
		 	} // deal with outliers, if tot goes out the boundary, choose the mean value of this bin			
      	
      } //check category
    }  // hit loop
    slew->Fill();
    
  } // end nevent loop

  forwrite->cd();
  hplus->Write();
  hminus->Write();
  slew->Write();
  forwrite->Close();

  cout << " Writing out TOT calibration parameters ... " << endl;

  char datafile[500];
  sprintf(datafile, "tray_%d_board_%d_%s_%d.dat",trayStart,TDIGStart, totOutput, iter);
  ofstream outData;
  outData.open(datafile); 

  for(Int_t i=0; i<nTray; i++) {
    for(Int_t j=0; j<nBoard; j++) {  //BUTTER LOOK HERE TDIG LOOP
      for(Int_t k=0; k<ntCell; k++) {
	outData << setw(4) << i+nTray*trayStart << setw(4) << j+TDIGStart*nBoard << setw(4) << k <<endl;
	if(tsfTOT[i][j][k]->GetCategory()==-1) {
	  outData << setw(4) << 0 << endl;;
	  continue;
	}
	outData << setw(4) << maxNBin << endl;
	if(fixNBin) {
	  for(Int_t l=0; l<=maxNBin; l++) {
	    outData << setw(16) << tofrTOTBins[i][j][k][l];
	  }
	  outData << endl;
	  for(Int_t l=0; l<=maxNBin; l++) {
	    outData << setw(16) << tsfTOT[i][j][k]->V(tofrTOTBins[i][j][k][l]);
	  }
	  outData << endl;	  
	  for(Int_t l=0; l<maxNBin; l++) {
	    outData << setw(16) << tofrTOT[i][j][k][l];
	  }
	  outData << endl;
	  if(tofrTOTbad[i][j][k]<1)
	  {
	  	for(Int_t l=0; l<maxNBin; l++) {
	    outData << setw(16) << tsfTOT[i][j][k]->V(tofrTOT[i][j][k][l]);
	  	}
	  	outData << endl;
	  }else{
	    for(Int_t l=0; l<maxNBin; l++) {
	    outData << setw(16) << tofrTOTbadcorr[i][j][k][l];
		}
	  	outData << endl;
	  }
	  
	}
	if(fixBinWidth) {
	  Double_t step = (maxTOT_TOFr-minTOT_TOFr)/maxNBin;
	  for(Int_t l=0; l<=maxNBin; l++) {
	    outData << setw(16) << minTOT_TOFr+step*l;
	  }
	  outData << endl;
	  for(Int_t l=0; l<=maxNBin; l++) {
	    outData << setw(16) << tsfTOT[i][j][k]->V(minTOT_TOFr+step*l);
	  }
	  outData << endl;
	}

      }//cell
    }
  } // end nTray loop
  outData.close();

  return 1;
}

Int_t outputZ(TChain* chain, Int_t iter)
{
  char rootfile[500];
  sprintf(rootfile, "tray_%d_board_%d_%s_%d.root",trayStart,TDIGStart, totOutput, iter);
  TChain *slewchain = new TChain("TOTPicoDst");
  slewchain->AddFile(rootfile);

  TOFrPicoDst *tofr = new TOFrPicoDst(chain);
  TOTPicoDst *slew = new TOTPicoDst(slewchain);

  Int_t hits;
  Float_t tof[maxTOFrHits];
  Float_t tag[maxTOFrHits];
  Float_t vzVpd;

  sprintf(rootfile, "tray_%d_board_%d_%s_%d.root",trayStart,TDIGStart, zOutput, iter);
  TFile *forwrite = new TFile(rootfile,"RECREATE");

  cout << " Writing to Z position correction tree ... " << endl;

  TTree *zhit = new TTree("ZPicoDst","ZPicoDst");
  zhit->Branch("hits", &hits, "hits/I");
  zhit->Branch("tof", tof, "tof[hits]/F");
  zhit->Branch("vzVpd",&vzVpd,"vzVpd/F");
  zhit->Branch("tag", tag, "tag[hits]/F");

  Int_t nentries = (int)(chain->GetEntries());
  for(int i=0; i<nentries; i++) {
    tofr->GetEntry(i);
    slew->GetEntry(i);

    hits = tofr->nTofHits;
    vzVpd = slew->vzVpd;

    for(Int_t j=0; j<hits; j++) {
      tag[j] = 0;
      tof[j] = -1.0e6;

      Int_t tray = tofr->tray[j]-1-nTray*trayStart;
      Int_t module = (tofr->module[j]-1)%4;//-moduleStart*nModule;
      Int_t cell = tofr->cell[j]-1;
      Double_t z = tofr->zLocal[j];//there was a bug here,  z has been subjected 1
      Int_t tcell = module*6+cell;

      if(tray>=nTray || tray<0) continue;
      if(module>=nModule || module<0) continue;
      if(cell>=nCell || cell<0) continue;

      Double_t oldtof = slew->tof[j];
	  Double_t tmpz=0;
	  Int_t kindex = 0;
      tag[j] = 1;
	  /*
      if(tsfZ[tray][module/(nModule/nBoard)]->GetCategory()==-1) {
      } else {
	  tof[j] = oldtof-tsfZ[tray][module/(nModule/nBoard)]->V(z);
      }
     */
     if(tsfZ[tray][0][tcell]->GetCategory()==-1) {
      } else {
         if(z>=tofrZ[tray][0][tcell][0]&&z<=tofrZ[tray][0][tcell][maxNZBin-1])
         {
            if(tofrZbad[tray][0][tcell]<1)
            {tof[j] = oldtof-tsfZ[tray][0][tcell]->V(z);}
		    else{
                  int ibin = -1;
                  for(int k=0;k<=maxNZBin;k++){
                    if(z>=tofrZ[tray][0][tcell][k]&&z<tofrZ[tray][0][tcell][k+1]){
                	ibin = k;
                	break;
                	} //if 
            	  } //for loop
              float tmpcor = linearInter(tofrZ[tray][0][tcell][ibin],tofrZ[tray][0][tcell][ibin+1],
              tofrZbadcorr[tray][0][tcell][ibin],tofrZbadcorr[tray][0][tcell][ibin+1], z);
              tof[j] = oldtof - tmpcor;
		    } //else fit not good
         } //if z in boundary
		 else{
		 	   if(z<tofrZ[tray][0][tcell][0]) 
		 		{tmpz = tofrZ[tray][0][tcell][0];
			     kindex = 0;}
			    else
				{tmpz = tofrZ[tray][0][tcell][maxNZBin-1];
			     kindex = maxNZBin-1;}
		
		 		if(tofrZbad[tray][0][tcell]<1)
            	{tof[j] = oldtof-tsfZ[tray][0][tcell]->V(tmpz);}
				else {
				tof[j] = oldtof - tofrZbadcorr[tray][0][tcell][kindex];
				}
		 	} // deal with outliers, if z goes out the boundary, choose the mean value of this bin			
      	
      	}  // else Category()
	  
    } //hit loop
    zhit->Fill();

  } // end nevent loop

  forwrite->cd();
  zhit->Write();
  forwrite->Close();

  cout << " Writing out Z calibration parameters ... " << endl;

  char datafile[500];
  sprintf(datafile, "tray_%d_board_%d_%s_%d.dat",trayStart,TDIGStart, zOutput, iter);
  ofstream outData;
  outData.open(datafile); 

  for(Int_t i=0; i<nTray; i++) {
    for(Int_t j=0; j<nBoard; j++) {   //BUTTER TDIG LOOP HERE
      for(Int_t k=0; k<ntCell; k++) {
	outData << setw(4) << i+nTray*trayStart << setw(4) << j+TDIGStart*nBoard << setw(4) <<k<< endl;
	if(tsfZ[i][j][k]->GetCategory()==-1) {
	  outData << setw(4) << 0 << endl;;
	  continue;
	}
	outData << setw(4) << maxNZBin << endl;
	if(fixNBin) {
	  for(Int_t l=0; l<=maxNZBin; l++) {
	    outData << setw(16) << tofrZBins[i][j][k][l];
	  }
	  outData << endl;
	  for(Int_t l=0; l<=maxNZBin; l++) {
	    outData << setw(16) << tsfZ[i][j][k]->V(tofrZBins[i][j][k][l]);
	  }
	  outData << endl;
	  
	  for(Int_t l=0; l<maxNZBin; l++) {
	    outData << setw(16) << tofrZ[i][j][k][l];
	  }
	  outData << endl;
	  /*
	  for(Int_t l=0; l<maxNZBin; l++) {
	    outData << setw(16) << tsfZ[i][j]->V(tofrZ[i][j][l]);
	  }
	  outData << endl;	
	  */
	  if(tofrZbad[i][j][k]<1)
	  {
	  	for(Int_t l=0; l<maxNZBin; l++) {
	    outData << setw(16) << tsfZ[i][j][k]->V(tofrZ[i][j][k][l]);
	  	}
	  	outData << endl;
	  }else
	  { for(Int_t l=0; l<maxNZBin; l++) {
	    outData << setw(16) << tofrZbadcorr[i][j][k][l];
	  	}
	  	outData << endl;
	  }

	  
	}
	if(fixBinWidth) {
	  Double_t step = (maxZ_TOFr-minZ_TOFr)/maxNZBin;
	  for(Int_t l=0; l<=maxNZBin; l++) {
	    outData << setw(16) << minZ_TOFr+step*l;
	  }
	  outData << endl;
	  for(Int_t l=0; l<=maxNZBin; l++) {
	    outData << setw(16) << tsfZ[i][j][k]->V(minZ_TOFr+step*l);
	  }
	  outData << endl;
	}
      }//cell
    }
  } // end nTray loop
  outData.close();

  return 1;
}

Int_t outputDelay(TChain* chain, Int_t iter)
{
  char rootfile[500];

  TOFrPicoDst *tofr = new TOFrPicoDst(chain);

  Int_t hits;
  Float_t tof[maxTOFrHits];
  Float_t tag[maxTOFrHits];
  Float_t vzVpd;

  sprintf(rootfile, "tray_%d_board_%d_%s_%d.root",trayStart,TDIGStart, delayOutput, iter);
  TFile *forwrite = new TFile(rootfile,"RECREATE");

  //cout << " "<<dt = "<<Writing to signal Delay correction tree ... " << endl;

  TTree *delay = new TTree("DelayPicoDst","DelayPicoDst");
  delay->Branch("hits", &hits, "hits/I");
  delay->Branch("tof", tof, "tof[hits]/F");
  delay->Branch("vzVpd", &vzVpd, "vzVpd/F");
  delay->Branch("tag", tag, "tag[hits]/F");

  Int_t nevt = (int)(chain->GetEntries());

  if(iter==0) {
    char rootfile[500];
/*    sprintf(rootfile, "%s.root", pvpdOutput);
    TChain *startchain = new TChain("StartPicoDst");
    startchain->AddFile(rootfile);
    StartPicoDst *start = new StartPicoDst(startchain);

    if(nevt!=startchain->GetEntries()) {
      cout << " Inconsistent TOFr and Start entry ! " << endl;
      return -1;
    }
*/
    for (Int_t i=0;i<nevt;i++) {
      tofr->GetEntry(i);
  //    start->GetEntry(i);

      hits = tofr->nTofHits;
      Double_t t0 = tofr->T0;
      vzVpd = tofr->vzVpd;

      for(Int_t j=0; j<hits; j++) {
        tof[j] = -1.0e6;
        tag[j] = 0.0;

        Int_t tray = tofr->tray[j]-1-nTray*trayStart;
        Int_t module = (tofr->module[j]-1)%4;//-moduleStart*nModule;
//	Int_t tmodule = module%4;
        Int_t cell = tofr->cell[j]-1;
        Double_t tdc = tofr->leTime[j];
	Int_t tcell = module*6+cell;
        tdc -= phaseDiff;
        while(tdc>51200) tdc -= 51200;
 
        if(tray>=nTray || tray<0) continue;
        if(module>=nModule || module<0) continue;
        if(cell>=nCell || cell<0) continue;

        tof[j] = tdc-t0-tofrDelay[tray][module][cell];
      }
      delay->Fill();
    } // end nevt loop
  } else {
    sprintf(rootfile, "tray_%d_board_%d_%s_%d.root",trayStart,TDIGStart, zOutput, iter);
    TChain *zhitchain = new TChain("ZPicoDst");
    zhitchain->AddFile(rootfile);
    ZPicoDst *zhit = new ZPicoDst(zhitchain);

    for (Int_t i=0;i<nevt;i++) {
      tofr->GetEntry(i);
      zhit->GetEntry(i);

      hits = tofr->nTofHits;
      vzVpd = zhit->vzVpd;

      for(Int_t j=0; j<hits; j++) {
        tof[j] = -1.0e6;
        tag[j] = 0.0;

        Int_t tray = tofr->tray[j]-1-nTray*trayStart;
        Int_t module = (tofr->module[j]-1)%4;//-moduleStart*nModule;
//	Int_t tmodule = module%4;
        Int_t cell = tofr->cell[j]-1;
        Double_t oldtof = zhit->tof[j];
	Int_t tcell=module*6+cell;

        if(tray>=nTray || tray<0) continue;
        if(module>=nModule || module<0) continue;
        if(cell>=nCell || cell<0) continue;

        tof[j] = oldtof-tofrDelay[tray][0][tcell];
      }
      delay->Fill();
    } // end nevt loop

  }
  forwrite->cd();
  delay->Write();
  forwrite->Close();  

  cout << " Writing out TOFr signal delay parameters ... " << endl;

  char datafile[500];
  sprintf(datafile, "tray_%d_board_%d_%s_%d.dat",trayStart,TDIGStart, delayOutput, iter);
  ofstream outData;
  outData.open(datafile); 

  for(Int_t i=0; i<nTray; i++) {
    for(Int_t j=0; j<nBoard; j++) {
      for(Int_t k=0; k<nCell; k++) {
	outData << setw(4) << i+nTray*trayStart << setw(4) << j+(TDIGStart%8)*4 << setw(4) << k << endl;
	outData << setw(16) << tofrDelay[i][j][k] << endl;
      }
    }
  } // end nTray loop
  outData.close();

  //char datafile[500];
  sprintf(datafile, "tray_%d_board_%d_%sRes_%d.dat",trayStart,TDIGStart, delayOutput, iter);
  ofstream outData2;
  outData2.open(datafile);

  for(Int_t i=0; i<nTray; i++) {
    for(Int_t j=0; j<nBoard; j++) {
      for(Int_t k=0; k<nCell; k++) {
	outData2 << setw(4) << i+nTray*trayStart << setw(4) << j+(TDIGStart%8)*4 << setw(4) << k << endl;
	outData2 << setw(16) << tofrRes[i][j][k] << endl;
      }
    }
  } // end nTray loop
  outData2.close();


  return 1;
}


Int_t checkReso(TChain* chain, Int_t iter)
{

  char rootfile[500];
  sprintf(rootfile, "%s.root", pvpdOutput);
//  TChain *startchain = new TChain("StartPicoDst");
//  startchain->AddFile(rootfile);
  sprintf(rootfile, "tray_%d_board_%d_%s_%d.root",trayStart,TDIGStart, totOutput, iter);
  TChain *slewchain = new TChain("TOTPicoDst");
  slewchain->AddFile(rootfile);
  sprintf(rootfile, "tray_%d_board_%d_%s_%d.root",trayStart,TDIGStart, zOutput, iter);
  TChain *zhitchain = new TChain("ZPicoDst");
  zhitchain->AddFile(rootfile);
  sprintf(rootfile, "tray_%d_board_%d_%s_%d.root",trayStart,TDIGStart, delayOutput, iter);
  TChain *delaychain = new TChain("DelayPicoDst");
  delaychain->AddFile(rootfile);

  TOFrPicoDst *tofr = new TOFrPicoDst(chain);
//  StartPicoDst *start = new StartPicoDst(startchain);
  TOTPicoDst *slew = new TOTPicoDst(slewchain);
  ZPicoDst *zhit = new ZPicoDst(zhitchain);
  DelayPicoDst *delay = new DelayPicoDst(delaychain);

  Int_t nevt = (int)(chain->GetEntries());
//  if(nevt!=startchain->GetEntries()) {
//    cout << " Inconsistent TOFr and Start entry ! " << endl;
//    return -1;
//  }
  if(nevt!=slewchain->GetEntries()) {
    cout << " Inconsistent TOFr and Slew entry ! " << endl;
    return -1;
  }
  if(nevt!=zhitchain->GetEntries()) {
    cout << " Inconsistent TOFr and Zhit entry ! " << endl;
    return -1;
  }

  // fill time resolution histogram
  for (Int_t i=0;i<nevt;i++) {
    tofr->GetEntry(i);
   // start->GetEntry(i);
    slew->GetEntry(i);
    zhit->GetEntry(i);
    delay->GetEntry(i);

    Int_t nhits = tofr->nTofHits;
    Int_t ieast = tofr->Ieast;
    Int_t iwest = tofr->Iwest;
    if(ieast<=0&&iwest<=0) continue;

    for(Int_t j=0; j<nhits; j++) {
      Double_t dedx = tofr->dedx[j];
      Double_t p = tofr->pt[j]*TMath::CosH(tofr->eta[j]);
//      if(sqrt(pow(tofr->vertexX,2)+pow(tofr->vertexY,2))>rCut) continue;
      if(TMath::Abs(tofr->vertexZ-zhit->vzVpd)>vzCut) continue;

      //Double_t dca = tofr->dcaX[j];
      Int_t fitpts = tofr->nHitsFit[j];
      Double_t dy = tofr->yLocal[j];

//      if(dca<0 || dca>1.5) continue;
      if(fitpts<25) continue;
      if(TMath::Abs(dy)>1.5) continue;

      //      Double_t nSigPi = tofr->nSigPi[j];
      //      Double_t nSigK = tofr->nSigK[j];
      //      Double_t nSigP = tofr->nSigP[j];
      //      Double_t nSigE = tofr->nSigE[j];
      //      if(p<0.2 || p>1.0 || 
      //	 TMath::Abs(nSigPi)>2.0 || TMath::Abs(nSigK)<2.0 ||
      //	 TMath::Abs(nSigP)<2.0 || TMath::Abs(nSigE)<2.0) continue; // another way to choose pion

      Double_t tracklength = tofr->length[j];
      Double_t betagamma = p/piMass;
      Double_t velocity = TMath::Sqrt(1.0/(1.0/betagamma/betagamma+1.0))*c_light;
      Double_t tofPi = tracklength/velocity;

      Int_t tray = tofr->tray[j]-1-nTray*trayStart;
      Int_t module = (tofr->module[j]-1)%4;//-moduleStart*nModule;
      Int_t cell = tofr->cell[j]-1;
      Int_t tcell = module*6+cell;
      Double_t tot = tofr->tot[j];
      Double_t z  = tofr->zLocal[j];	  

      Double_t tof_s = slew->tof[j];
      Double_t tof_z = zhit->tof[j];
      Double_t tof_d = delay->tof[j];
    
      if(tray>=nTray || tray<0) continue;
      if(module>=nModule || module<0) continue;
      if(cell>=nCell || cell<0) continue;
      if(tot<tofrTOTBins[tray][0][tcell][0] ||
	 tot>tofrTOTBins[tray][0][tcell][maxNBin]) continue;
	  if(z<tofrZBins[tray][0][tcell][0] ||
	 z>tofrZBins[tray][0][tcell][maxNZBin]) continue; 
	  
      //see 1/beta vs p distribution
      //after slewing correction
      Double_t invBeta_s = tof_s/tracklength*c_light;
      Double_t ms_s = pow(p,2)*(pow(invBeta_s,2)-1);
      hBetaP_TOT->Fill(p, invBeta_s);
      hMassP_TOT->Fill(p, ms_s);
      if(tofr->charge[j]==1) {
         hBetaP_TOT_Pos->Fill(p, invBeta_s);
         hMassP_TOT_Pos->Fill(p, ms_s);
      } else if(tofr->charge[j]==-1) {
         hBetaP_TOT_Neg->Fill(p, invBeta_s);
         hMassP_TOT_Neg->Fill(p, ms_s);
      }

      //after zlocal correction  
      Double_t invBeta_z = tof_z/tracklength*c_light;
      Double_t ms_z = pow(p,2)*(pow(invBeta_z,2)-1);
      hBetaP_Z->Fill(p, invBeta_z);
      hMassP_Z->Fill(p, ms_z);
      if(tofr->charge[j]==1) {
        hBetaP_Z_Pos->Fill(p, invBeta_z);
        hMassP_Z_Pos->Fill(p, ms_z);
      } else if(tofr->charge[j]==-1) {
        hBetaP_Z_Neg->Fill(p, invBeta_z);
        hMassP_Z_Neg->Fill(p, ms_z);
      }        

     //after T0 correction  
      Double_t invBeta_d = tof_d/tracklength*c_light;
      Double_t ms_d = pow(p,2)*(pow(invBeta_z,2)-1);
      hBetaP_D->Fill(p, invBeta_d);
      hMassP_D->Fill(p, ms_d);
      if(tofr->charge[j]==1) {
        hBetaP_D_Pos->Fill(p, invBeta_d);
        hMassP_D_Pos->Fill(p, ms_d);
      } else if(tofr->charge[j]==-1) {
        hBetaP_D_Neg->Fill(p, invBeta_d);
        hMassP_D_Neg->Fill(p, ms_d);
      }

      if(dedx > dedxcut) continue; // select pion to see resolution
      Double_t dt;
      dt = tof_s-tofPi;
      hResoTOT->Fill(p, dt);
      if(p>0.3&&p<0.6) hResoTOT1D->Fill(dt);
      dt = tof_z-tofPi;
      hResoZ->Fill(p, dt);
      if(p>0.3&&p<0.6)  hResoZ1D->Fill(dt);
      dt = tof_d-tofPi;
      hResoD->Fill(p,dt);
      if(p>0.3&&p<0.6) hResoD1D->Fill(dt);

      if(tofr->charge[j]==1) {//positive tracks
        dt = tof_s-tofPi;
        hResoTOT_Pos->Fill(p, dt);
        if(p>0.3&&p<0.6) hResoTOT1D_Pos->Fill(dt);

        dt = tof_z-tofPi;
        hResoZ_Pos->Fill(p, dt);
        if(p>0.3&&p<0.6) hResoZ1D_Pos->Fill(dt); 

        dt = tof_d-tofPi;
        hResoD_Pos->Fill(p,dt);
        if(p>0.3&&p<0.6) hResoD1D_Pos->Fill(dt);
      } else if(tofr->charge[j]==-1) {//negative tracks
        dt = tof_s-tofPi;
        hResoTOT_Neg->Fill(p, dt);
        if(p>0.3&&p<0.6) hResoTOT1D_Neg->Fill(dt);

        dt = tof_z-tofPi;
        hResoZ_Neg->Fill(p, dt);
        if(p>0.3&&p<0.6) hResoZ1D_Neg->Fill(dt);    

        dt = tof_d-tofPi;
        hResoD_Neg->Fill(p,dt);
        if(p>0.3&&p<0.6) hResoD1D_Neg->Fill(dt);
      }
    }
  } // end nevt loop

  TF1 *g = new TF1("g","gaus");
  float nsigma = 2.5;
  float mean_s, mean_z, mean_d, res_s, res_z, res_d;
  mean_s = hResoTOT1D->GetMean();
  res_s = hResoTOT1D->GetRMS();
  g->SetRange(mean_s-nsigma*res_s, mean_s+nsigma*res_s);
  hResoTOT1D->Fit("g","NRQ");
  mean_s = g->GetParameter(1);
  res_s = g->GetParameter(2);
  g->SetRange(mean_s-nsigma*res_s, mean_s+nsigma*res_s);
  hResoTOT1D->Fit("g","NRQ");
  mean_s = g->GetParameter(1);
  res_s = g->GetParameter(2);
  
  mean_z = hResoZ1D->GetMean();
  res_z = hResoZ1D->GetRMS();
  g->SetRange(mean_z-nsigma*res_z, mean_z+nsigma*res_z);
  hResoZ1D->Fit("g","NRQ");
  mean_z = g->GetParameter(1);
  res_z = g->GetParameter(2);
  g->SetRange(mean_z-nsigma*res_z, mean_z+nsigma*res_z);
  hResoZ1D->Fit("g","NRQ");
  mean_z = g->GetParameter(1);
  res_z = g->GetParameter(2);

  mean_d = hResoD1D->GetMean();
  res_d = hResoD1D->GetRMS();
  g->SetRange(mean_d-nsigma*res_d, mean_d+nsigma*res_d);
  hResoD1D->Fit("g","NRQ");
  mean_d = g->GetParameter(1);
  res_d = g->GetParameter(2);
  g->SetRange(mean_d-nsigma*res_d, mean_d+nsigma*res_d);
  hResoD1D->Fit("g","NRQ");
  mean_d = g->GetParameter(1);
  res_d = g->GetParameter(2);

  cout<<"***********************************************************"<<endl;
  cout<<"* "<<iter<<"th iteration"<<endl;
  cout<<"* after slewing correction: "<<mean_s<<"  "<<res_s<<" ns"<<endl;
  cout<<"* after TZ      correction: "<<mean_z<<"  "<<res_z<<" ns"<<endl;
  cout<<"* after T0      correction: "<<mean_d<<"  "<<res_d<<" ns"<<endl;
  cout<<"***********************************************************"<<endl;


  sprintf(rootfile, "tray_%d_board_%d_%s_%d.root",trayStart,TDIGStart, resoOutput, iter);
  TFile *forwrite = new TFile(rootfile,"RECREATE");

  forwrite->cd();
  hBetaP_TOT->Write();
  hBetaP_Z->Write();
  hBetaP_D->Write();
  hMassP_TOT->Write();
  hMassP_Z->Write();
  hMassP_D->Write();

  hBetaP_TOT_Pos->Write();
  hBetaP_Z_Pos->Write();
  hBetaP_D_Pos->Write();
  hMassP_TOT_Pos->Write();
  hMassP_Z_Pos->Write();
  hMassP_D_Pos->Write();

  hBetaP_TOT_Neg->Write();
  hBetaP_Z_Neg->Write();
  hBetaP_D_Neg->Write();
  hMassP_TOT_Neg->Write();
  hMassP_Z_Neg->Write();
  hMassP_D_Neg->Write();

  hResoTOT->Write();
  hResoZ->Write();
  hResoD->Write();
  hResoTOT1D->Write();
  hResoZ1D->Write();
  hResoD1D->Write();
  hResoTOT_Pos->Write();
  hResoZ_Pos->Write();
  hResoD_Pos->Write();
  hResoTOT1D_Pos->Write();
  hResoZ1D_Pos->Write();
  hResoD1D_Pos->Write();
  hResoTOT_Neg->Write();
  hResoZ_Neg->Write();
  hResoD_Neg->Write();
  hResoTOT1D_Neg->Write();
  hResoZ1D_Neg->Write();
  hResoD1D_Neg->Write();

  forwrite->Close();

  hBetaP_TOT->Reset();
  hBetaP_Z->Reset();
  hBetaP_D->Reset();
  hMassP_TOT->Reset();
  hMassP_Z->Reset();
  hMassP_D->Reset();
  hBetaP_TOT_Pos->Reset();
  hBetaP_Z_Pos->Reset();
  hBetaP_D_Pos->Reset();
  hMassP_TOT_Pos->Reset();
  hMassP_Z_Pos->Reset();
  hMassP_D_Pos->Reset();
  hBetaP_TOT_Neg->Reset();
  hBetaP_Z_Neg->Reset();
  hBetaP_D_Neg->Reset();
  hMassP_TOT_Neg->Reset();
  hMassP_Z_Neg->Reset();
  hMassP_D_Neg->Reset();

  hResoTOT->Reset();
  hResoZ->Reset();
  hResoD->Reset();
  hResoTOT1D->Reset();
  hResoZ1D->Reset();
  hResoD1D->Reset();
  hResoTOT_Pos->Reset();
  hResoZ_Pos->Reset();
  hResoD_Pos->Reset();
  hResoTOT1D_Pos->Reset();
  hResoZ1D_Pos->Reset();
  hResoD1D_Pos->Reset();
  hResoTOT_Neg->Reset();
  hResoZ_Neg->Reset();
  hResoD_Neg->Reset();
  hResoTOT1D_Neg->Reset();
  hResoZ1D_Neg->Reset();
  hResoD1D_Neg->Reset();

  return 1;
}

float linearInter(float x1, float x2, float y1, float y2, float x){
 // float y;
  if(x2==x1) 
  {cout<<"error in linear interation, x2==x1"<<endl;
   return 0;
  }
  return ((x-x1)*y2+(x2-x)*y1)/(x2-x1);
}

Bool_t xpisnan(Double_t x){
return (x!=x);
}

