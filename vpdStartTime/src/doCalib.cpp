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
 #include "TLeaf.h"

using namespace std;
using std::cout;   
using std::endl;  

#endif

#include "TSplineFit.h"

#include "TOFrPicoDst.h"
#include "StartPicoDst.h"
#include "TOTPicoDst.h"
#include "ZPicoDst.h"
#include "DelayPicoDst.h"

#include "cuts.h"

/*
#include "T0PicoDst.h"
*/



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
const Int_t maxiteration = 20; 
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


bool skipBadRun( TOFrPicoDst *t1 ) {
    if ( t1->run == 14128018){
      return true;
    }
    return false;
}
bool skipBadRun( TChain *chain ) {
    int run = chain->GetLeaf("run")->GetValue();
    if ( run == 14128018){
      return true;
    }
    return false;
}





Int_t main(int argc, char **argv)
{
  if(argc!=4 && argc!=1) return 0;
  
  const char *FileInput = "test.list";

 TString *string;
 TString *string2;
 
  if(argc==1){
    string = new TString("0");
  }
 
  if(argc==4){
    FileInput = argv[1];
    string = new TString(argv[2]);
    string2 = new TString(argv[3]);    
  }
  trayStart = string->Atoi();
  TDIGStart = string2->Atoi();

  if(argc!=1 && argc!=4) return -1;

  
  
  
  Int_t fileNumber = 0;
  Char_t FileList[512];

  TChain *tofrchain = new TChain("tof");

  ifstream* inputStream = new ifstream;
  inputStream->open(FileInput);

  if (!(inputStream)) {
    printf("can not open list file\n");
    return 0;
  } else {
    printf("Opened list file\n");
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
  

  
  

  if(fixNBin) {
    binning(tofrchain, 0); 
  }
  initPVPD();

  if ( FIRST_RUN ){
    doStart(tofrchain);  
    outputPVPD(tofrchain); 
    return 1;
  } else {
    outputparPVPD(tofrchain); 
    return 1;  
 }
 
}

Int_t binning(TChain* chain, Int_t choice)
{
  Int_t nevents = (int)chain->GetEntries();

  if(choice==0) { 
    TOFrPicoDst *t1 = new TOFrPicoDst(chain);
    cout<< " woking on binning! choice = "<<choice<<endl;
 
    for(Int_t i=0;i<nevents;i++) {
      t1->GetEntry(i);
      if ( skipBadRun( t1 ) )
        continue;

      if(i%(nevents/10)==0) cout<< i <<" events processed!"<<endl;     
      Int_t ieast = t1->numberOfVpdEast;
      Int_t iwest = t1->numberOfVpdWest;
      
      if(ieast<minpVPDHits&&iwest<minpVPDHits) continue;
     
      if(iwest>=minpVPDHits){
        for(Int_t j=0; j<nPVPDChannel/2; j++) {
		
	

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
          pvpdTOTBins[i][j] = step*j+minTOT_pVPD; 


	}
	cout << " pVPD Channel " << i << " seems dead ! " << "( " << size << " hits)" <<endl;
      } else {
	Int_t step = size/maxpvpdNBin; 
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
    } 
    
  } else { 
    TOFrPicoDst *t2 = new TOFrPicoDst(chain);

    char rootfile[500];

    for(Int_t i=0;i<nevents;i++) {
      t2->GetEntry(i);


      Double_t t0 = t2->T0;
      Int_t hits = t2->nTofHits;
      if(hits<=0) continue;

      for(Int_t j=0; j<hits; j++) {
	Int_t tray = t2->tray[j]-1-nTray*trayStart;
	Int_t module = (t2->module[j]-1)%4;
	Int_t cell = t2->cell[j]-1;
	Int_t tcell = module*6+cell;

	if(tray>=nTray || tray<0) continue;
	if(module>=nModule || module<0) continue;
	if(cell>=nCell || cell<0) continue;	
	
	Double_t p = t2->pt[j]*TMath::CosH(t2->eta[j]);
	Double_t dedx = t2->dedx[j];
	if(p<0.3 || p>0.6 || dedx > dedxcut) continue;
        if(sqrt(pow(t2->vertexX,2)+pow(t2->vertexY,2))>rCut) continue;    

         if(TMath::Abs(t2->vertexZ-t2->vzVpd)>vzCut) continue; 


	Double_t tracklength = t2->length[j];
	Double_t betagamma = p/piMass;
	Double_t velocity = TMath::Sqrt(1.0/(1.0/betagamma/betagamma+1.0))*c_light;
	Double_t tofPi = tracklength/velocity;

	Double_t tdc = t2->leTime[j];
        tdc -= phaseDiff;
        while(tdc>51200) tdc-=51200;
	Double_t dt = tdc-t0-tofPi-tofrDelay[tray][0][tcell];

	if(TMath::Abs(dt)>toftlimit) continue;
 
	if(choice==1) { 
	  Double_t tot = t2->tot[j];

	  if(tot>minTOT_TOFr && tot<maxTOT_TOFr) tofrtots[tray][0][tcell].push_back(tot);
	  
	}
	if(choice==2) { 
	  Double_t z = t2->zLocal[j];
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
		tofrTOTBins[i][j][k][l] = step*l; 
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
		tofrZBins[i][j][k][l] = step*l; 
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
	}
      }
    } 
  } 

  
  if(choice==0) {
    for(Int_t i=0; i<nPVPDChannel; i++) {
      pvpdtots[i].clear();
    }
  } else {
    for(Int_t i=0; i<nTray; i++) {      
      for(Int_t j=0; j<nBoard; j++) {
	for(Int_t k=0; k<ntCell; k++) {
	  tofrtots[i][j][k].clear();
	  tofrzs[i][j][k].clear();	
	}
      }
    }  
  }
  cout<<"TOFr tot	"<<endl;
  for(int tmpk=0; tmpk<maxNBin; tmpk++)
	cout<<setw(10)<<tofrTOT[0][0][0][tmpk];
  return 1;
}

Int_t initPVPD()
{
  

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

  sprintf(buf, "VzCorrvzdiff");
  hVzCorrvzdiff = new TH2D(buf, buf, 400, -200, 200, 400, -200, 200);

  sprintf(buf, "VzCorr2vzdiff");
  hVzCorr2vzdiff = new TH2D(buf, buf, 400, -200, 200, 500, -50, 50);

  sprintf(buf, "vxvvy");
  hvxvvy = new TH2D(buf, buf, 500, -10, 10, 500, -10, 10);
  sprintf(buf, "vxvvyvzdiff");
  hvxvvyvzdiff = new TH2D(buf, buf, 500, -10, 10, 500, -10, 10);

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
  Double_t nfitsigma = 2.0;
  Double_t tdcwidth = 1.0; 

  TF1 *f1 = new TF1("f1","gaus");
  TF1 *polfit = new TF1("polfit","pol4");

  for(Int_t j=0; j<nPVPDChannel; j++) {
    cor[j] = 0.;
  }

  
  Int_t nevents = (int)chain->GetEntries();
  
  cout << "== total entries : " << nevents << endl;

  TOFrPicoDst *t = new TOFrPicoDst(chain);

  Char_t title[500];
  sprintf(title, "%s.ps", pvpdOutput);
  
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
      hPVPDCorr[j]->Reset(); 
      htPVPD[j]->Reset();
      htPVPD1st[j]->Reset();
    }

    for(Int_t i=0; i<nevents; i++) {
      chain->GetEntry(i);
      
      if ( skipBadRun( chain ) )
        continue;

      Int_t flag[nPVPDChannel]; 
      for(int i=0;i<nPVPDChannel;i++){
		  flag[i]=1; 

		}


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

          tot[j] = t->vpdTotWest[j];
          if(tot[j]<=minTOT_pVPD || tot[j]>maxTOT_pVPD) continue;
          tdc[j] = t->vpdLeWest[j];
          if(ir==0) { }
          else cor[j] = tsfPVPD[j]->V(tot[j]);

          if(ir>=4){
            tdccorr[j] = tdc[j]-cor[j];
            if(tdccorr[j]>51200) tdccorr[j] -= 51200;
            if(tdccorr[j]<0) tdccorr[j] += 51200;
            tdcsum = tdcsum + tdccorr[j];
            countwest = countwest + 1.0;
            
          }

        }
        

        if(ir==maxiteration-1) {
          for(Int_t j=0; j<nPVPDChannel/2; j++) {
            if(tot[j]<=minTOT_pVPD || tot[j]>maxTOT_pVPD) continue;
            for(Int_t k=0; k<nPVPDChannel/2; k++) {


              if(k==j) continue;
              if(tot[k]<=minTOT_pVPD || tot[k]>maxTOT_pVPD) continue;
              double deltatime = tdccorr[k] - tdccorr[j];
              
              if (deltatime >51200) deltatime -= 51200;
              hdeltatwest->Fill(deltatime);
            }
          }
        } 

        if(ir>=4) {
          average = tdcsum/countwest;
          Double_t rmssum=0;
          Double_t vpdHit[nPVPDChannel/2] = {0}; 
          Int_t *hitIndex = new Int_t[nPVPDChannel/2]; 
          for(Int_t j=0; j<nPVPDChannel/2; j++) {
	

            if(tot[j]<=minTOT_pVPD || tot[j]>maxTOT_pVPD) continue;
            if(fabs(tdccorr[j]-average)>tdiffCut) flag[j] = 0; 
            
            rmssum = (tdccorr[j]-average)*(tdccorr[j]-average);
            vpdHit[j] = tdc[j]-cor[j]; 
          }
          Double_t rmswest = sqrt(rmssum/countwest);
          
          if (rmswest>mrmscut) largermswest = 1;
          
          int countHits = 0;
          TMath::Sort(nPVPDChannel/2,vpdHit,hitIndex);
          for(int i=0;i<nPVPDChannel/2;i++) if(vpdHit[i]!=0) countHits++;
          int nrejected = (int)(pCut*countHits+0.5);
          for(int i=0;i<nrejected;i++) flag[hitIndex[i]] = 0;
          delete hitIndex;

        }

        if(!largermswest) {
          for(Int_t j=0; j<nPVPDChannel/2; j++) {


            if(!flag[j]) continue; 
            if(tot[j]<=minTOT_pVPD || tot[j]>maxTOT_pVPD) continue;
            if(ir==0||ir>=4)htotPVPD[j]->Fill(tot[j]);
            Double_t dt = 0.0;
            Int_t dn = 0;
            for(Int_t k=0; k<nPVPDChannel/2; k++) {

              if(k==j) continue;
              if(ir>2 && tmpres[j]>5) continue;
              if(tot[k]<=minTOT_pVPD || tot[k]>maxTOT_pVPD) continue;
              if(!flag[k]) continue; 
              
              dt += (tdc[k]-cor[k]);
              dn++;
              if(dn==(minpVPDHits-1)&&(ir==maxiteration-1)) break;
            }
            
            if(dn<minpVPDHits-1) continue;
            
            if(ir>=0) {
              hPVPD[j]->Fill(tot[j], tdc[j]-dt/dn);
              hPVPDCorr[j]->Fill(tot[j],tdc[j]-cor[j]-dt/dn); 
            }
            else if((tdc[j]-cor[j]-dt/dn)>(tmpmean[j]-5*tmpres[j]) && (tdc[j]-cor[j]-dt/dn)<(tmpmean[j]+5*tmpres[j])) {
              hPVPD[j]->Fill(tot[j], tdc[j]-dt/dn);
              hPVPDCorr[j]->Fill(tot[j],tdc[j]-cor[j]-dt/dn);
            }
            if(ir<=1) htPVPD1st[j]->Fill(tdc[j]-cor[j]-dt/dn);
            else htPVPD[j]->Fill(tdc[j]-cor[j]-dt/dn);
          }
        }
      }

      if(Ieast>=minpVPDHits) {

        Double_t average=0;
        Double_t tdcsum=0;
        Int_t largermseast = 0;
        Double_t counteast = 0;

        for(Int_t j=nPVPDChannel/2; j<nPVPDChannel; j++) {
          tot[j] = t->vpdTotEast[j-nPVPDChannel/2];
          if(tot[j]<=minTOT_pVPD || tot[j]>maxTOT_pVPD) continue;
          tdc[j] = t->vpdLeEast[j-nPVPDChannel/2];

          if(ir==0) { }
          else cor[j] = tsfPVPD[j]->V(tot[j]);

          if(ir>=4){
            tdccorr[j] = tdc[j]-cor[j];
            if(tdccorr[j]>51200) tdccorr[j] -= 51200;
            if(tdccorr[j]<0) tdccorr[j] += 51200;
            tdcsum = tdcsum + tdccorr[j];
            counteast = counteast + 1.0;
            
          }

        } 

        if(ir==maxiteration-1) {
          for(Int_t j=nPVPDChannel/2; j<nPVPDChannel; j++) {
            if(tot[j]<=minTOT_pVPD || tot[j]>maxTOT_pVPD) continue;
            for(Int_t k=nPVPDChannel/2; k<nPVPDChannel; k++) {
              if(k==j) continue;
              if(tot[k]<=minTOT_pVPD || tot[k]>maxTOT_pVPD) continue;
              double deltatime = tdccorr[k] - tdccorr[j];
              
              if (deltatime >51200) deltatime -= 51200;
              hdeltateast->Fill(deltatime);
            }
          }
        } 

        
        if(ir>=4) {
          average = tdcsum/counteast;
          Double_t rmssum=0;
          Double_t vpdHit[nPVPDChannel/2] = {0}; 
          Int_t *hitIndex = new Int_t[nPVPDChannel/2]; 
          for(Int_t j=nPVPDChannel/2; j<nPVPDChannel; j++) {
            if(tot[j]<=minTOT_pVPD || tot[j]>maxTOT_pVPD) continue;
            if(fabs(tdccorr[j]-average)>tdiffCut) flag[j] = 0; 
            
            rmssum = (tdccorr[j]-average)*(tdccorr[j]-average);
            vpdHit[j-nPVPDChannel/2] = tdc[j] - cor[j]; 
          }
          Double_t rmseast = sqrt(rmssum/counteast);
          
          if (rmseast>mrmscut) largermseast = 1;
          
          int countHits = 0;
          TMath::Sort(nPVPDChannel/2,vpdHit,hitIndex);
          for(int i=0;i<nPVPDChannel/2;i++) if(vpdHit[i]!=0) countHits++;
          int nrejected = (int)(pCut*countHits+0.5);
          for(int i=0;i<nrejected;i++) flag[hitIndex[i]+nPVPDChannel/2] = 0;
          delete hitIndex;
          
          
          
          
          
          
          
          
        }

        if(!largermseast) {
          for(Int_t j=nPVPDChannel/2; j<nPVPDChannel; j++) {
            if(!flag[j]) continue; 
            if(tot[j]<=minTOT_pVPD || tot[j]>maxTOT_pVPD) continue;
            if(ir==0||ir>=4)htotPVPD[j]->Fill(tot[j]);
            Double_t dt = 0.0;
            Int_t dn = 0;
            for(Int_t k=nPVPDChannel/2; k<nPVPDChannel; k++) {
              if(k==j) continue;
              if(ir>2&&tmpres[j]>5) continue;
              if(tot[k]<=minTOT_pVPD || tot[k]>maxTOT_pVPD) continue;
              if(!flag[k]) continue; 
              dt += (tdc[k]-cor[k]);
              dn++;
              if((dn==minpVPDHits-1)&&(ir==maxiteration-1)) break;
            }

            if(dn<minpVPDHits-1) continue;
            if(ir>=0) {
              hPVPD[j]->Fill(tot[j], tdc[j]-dt/dn);
              hPVPDCorr[j]->Fill(tot[j],tdc[j]-cor[j]-dt/dn); 
            }
            else if((tdc[j]-cor[j]-dt/dn)>(tmpmean[j]-5*tmpres[j]) && (tdc[j]-cor[j]-dt/dn)<(tmpmean[j]+5*tmpres[j])) {
              hPVPD[j]->Fill(tot[j], tdc[j]-dt/dn);
              hPVPDCorr[j]->Fill(tot[j],tdc[j]-cor[j]-dt/dn); 
            }
            if(ir<=1) htPVPD1st[j]->Fill(tdc[j]-cor[j]-dt/dn);
            else htPVPD[j]->Fill(tdc[j]-cor[j]-dt/dn);
          }
        }
      }

    } 



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
 
        mean = htPVPD1st[k]->GetMean();
        sigma = htPVPD1st[k]->GetRMS();
        f1->SetRange(mean-nfitsigma*sigma, mean+nfitsigma*sigma);
        htPVPD1st[k]->Fit("f1","QR");
        tmpres[k] = f1->GetParameter(2)*tdcwidth;
        tmpmean[k] = f1->GetParameter(1)*tdcwidth;
      } else {
        mean = htPVPD[k]->GetMean();
        sigma = htPVPD[k]->GetRMS();
        f1->SetRange(mean-nfitsigma*sigma, mean+nfitsigma*sigma);
        htPVPD[k]->Fit("f1","NQR");
        tmpres[k] = f1->GetParameter(2)*tdcwidth;
        tmpmean[k] = f1->GetParameter(1)*tdcwidth;
        f1->SetRange((tmpmean[k]-nfitsigma*tmpres[k])/tdcwidth, (tmpmean[k]+nfitsigma*tmpres[k])/tdcwidth);
        htPVPD[k]->Fit("f1","NQR");
        tmpres[k] = f1->GetParameter(2)*tdcwidth;
        tmpmean[k] = f1->GetParameter(1)*tdcwidth;
        f1->SetRange((tmpmean[k]-nfitsigma*tmpres[k])/tdcwidth, (tmpmean[k]+nfitsigma*tmpres[k])/tdcwidth);
        htPVPD[k]->Fit("f1","QR");
        tmpres[k] = f1->GetParameter(2)*tdcwidth;
        tmpmean[k] = f1->GetParameter(1)*tdcwidth;


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

      
      sprintf(buf, "pVPD_Channel_%d_1", j);
      hPVPD[j]->FitSlicesY();
      hpPVPD[j] = (TH1D *)gDirectory->Get(buf); 
      for(int k=0; k<maxpvpdNBin; k++){
        pvpdTOTX[j][k] = hpPVPD[j]->GetBinCenter(k+1);
        double content=hpPVPD[j]->GetBinContent(k+1);
        double meanyvalue = hPVPD[j]->GetMean(2);
        hPVPD[j]->GetXaxis()->SetRange(k+1,k+1);
        double meanvalueproj = hPVPD[j]->GetMean(2);
        double meanvalueprojerr = hPVPD[j]->GetMeanError(2);
        
        
        double deltat=10;
        if(ir<4) deltat=5;
        else if(ir>=4) deltat=0.3;
        

        if(fabs(content)>100||content>(meanvalueproj+deltat)||content<(meanvalueproj-deltat)) {
          
          hpPVPD[j]->SetBinContent(k+1,meanvalueproj);
          hpPVPD[j]->SetBinError(k+1,meanvalueprojerr);
        }

      }
      
      sprintf(buf, "pVPD Fit Channel_%d", j);
      sprintf(buff, "pVPD Fit | Channel_%d", j);
      if(tsfPVPD[j]) delete tsfPVPD[j];
      tsfPVPD[j] = new TSplineFit(buf, buff, maxpvpdNBin, tsfM, hpPVPD[j]);
      
      tsfPVPD[j]->DrawHere(max((Int_t)hpPVPD[j]->GetMinimum()-2, (Int_t) -tlimit), min(int(hpPVPD[j]->GetMaximum())+2, (Int_t) tlimit));

    }

    

    c2->Update();

    
    TFile *f = new TFile("pvpdCheck.root","UPDATE");
    char buf[100];
    for(Int_t k=0; k<nPVPDChannel; k++) {
      sprintf(buf,"pVPDCorr_Ch%d_it%d",k,ir+1);
      hPVPDCorr[k]->SetName(buf);
      hPVPDCorr[k]->Write();
    }
    f->Close();
    delete f;
    


    done = kTRUE;
    for(Int_t k=0; k<nPVPDChannel; k++) {
      
      if(TMath::Abs(res[k]-tmpres[k])>smallvalue) { done = kFALSE; }

      res[k] = tmpres[k];
      
      if(ir<=1) hPVPDMeanvsIt[k]->SetBinContent(ir+1,htPVPD1st[k]->GetMean());
      else hPVPDMeanvsIt[k]->SetBinContent(ir+1,htPVPD[k]->GetMean());
      hPVPDFitMeanvsIt[k]->SetBinContent(ir+1,f1->GetParameter(1));
      hPVPDFitSigmavsIt[k]->SetBinContent(ir+1,f1->GetParameter(2));
      hPVPDResovsIt[k]->SetBinContent(ir+1,res[k]);
      
    }
    ir++;
  } 

  ps->Close();

  

  cout << "==============================" << endl;
  cout << endl;
  cout << "  Correction over !  Writing out ... " << endl;

  return 1;
}



Int_t outputPVPD(TChain* chain)
{
  
  Int_t Ieast;
  Int_t Iwest;
  Double_t sumEast;
  Double_t sumWest;
  Double_t T0;
  Double_t T0vzdiff;
  Double_t vzVpd=-1e6;

  const Float_t vzOffset = 9.67581530999990633e+01;


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
    if ( skipBadRun( tofr ))
      continue;

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

      tot[j] = j<nPVPDChannel/2? tofr->vpdTotWest[j]:tofr->vpdTotEast[j-nPVPDChannel/2];
      tdc[j] =  j<nPVPDChannel/2? tofr->vpdLeWest[j]:tofr->vpdLeEast[j-nPVPDChannel/2];
      
      
      if(tot[j]>minTOT_pVPD && tot[j]<maxTOT_pVPD) {
        cor[j] = tsfPVPD[j]->V(tot[j]);
        if(j<nPVPDChannel/2) {
          
          sumWest += tdc[j] - cor[j];
          Iwest++;
        } else {
         
          sumEast += tdc[j] - cor[j];
          Ieast++;
        }
      }
    }


    
    int flag[nPVPDChannel];
    for(int i=0;i<nPVPDChannel;i++){
		 flag[i] = 1; 



		}
    Double_t vpdWHit[nPVPDChannel/2] = {0};
    Double_t vpdEHit[nPVPDChannel/2] = {0};
    for(Int_t j=0; j<nPVPDChannel; j++) {

      if(tot[j]<=minTOT_pVPD || tot[j]>=maxTOT_pVPD) continue;
      if(j<nPVPDChannel/2) vpdWHit[j] = tdc[j]-cor[j];
      else vpdEHit[j-nPVPDChannel/2] = tdc[j]-cor[j];
    }
    
    int countWHits = 0;
    Int_t *hitIndexW = new Int_t[nPVPDChannel/2];
    TMath::Sort(nPVPDChannel/2,vpdWHit,hitIndexW);
    for(int i=0;i<nPVPDChannel/2;i++) if(vpdWHit[i]!=0) countWHits++;
    int nWRejected = (int)(pCut*countWHits+0.5);
    for(int i=0;i<nWRejected;i++) flag[hitIndexW[i]] = 0; 
    delete hitIndexW;
    
    int countEHits = 0;
    Int_t *hitIndexE = new Int_t[nPVPDChannel/2];
    TMath::Sort(nPVPDChannel/2,vpdEHit,hitIndexE);
    for(int i=0;i<nPVPDChannel/2;i++) if(vpdEHit[i]!=0) countEHits++;
    int nERejected = (int)(pCut*countEHits+0.5);
    for(int i=0;i<nERejected;i++) flag[hitIndexE[i]+nPVPDChannel/2] = 0; 
    delete hitIndexE;
    
    
    
    
    
    
    
    


    if(i<10) cout<<"Iwest/Ieast = "<<Iwest<<"/"<<Ieast<<endl;
    Double_t vpdtime;
    
    for(int j=0;j<nPVPDChannel;++j){



      if(tot[j]<=minTOT_pVPD ||tot[j]>=maxTOT_pVPD) continue;
      if(j<nPVPDChannel/2&&Iwest>1) {
        vpdtime = ((tdc[j]-cor[j])*Iwest-sumWest)/(Iwest-1);
        
        
        if(fabs(vpdtime)>tdiffCut || flag[j]==0) { 
          sumWest -= tdc[j]-cor[j];
          Iwest--;
        }
      }
      if(j>=nPVPDChannel/2&&Ieast>1) {
        vpdtime = ((tdc[j]-cor[j])*Ieast-sumEast)/(Ieast-1);
        
        
        if(fabs(vpdtime)>tdiffCut || flag[j]==0) { 
          sumEast -= tdc[j]-cor[j];
          Ieast--;
        }
      }
    }
    

    sumEast += 2*vzOffset*Ieast/c_light;

    if(Ieast>0||Iwest>0) {
      
      
      
      
      float tdiff;
      if(Ieast>0&&Iwest>0) {

        tdiff = sumEast/Ieast - sumWest/Iwest;
        vzVpd = tdiff/2*c_light;
      } else {
        vzVpd = -1e6;
      }
      



    
      hVzCorr->Fill(vzVpd, zvtx);
      hVzCorr2->Fill(vzVpd, vzVpd-zvtx);



	if(fabs(vzVpd-zvtx)<=vzCut){
      hVzCorrvzdiff->Fill(vzVpd, zvtx);
      hVzCorr2vzdiff->Fill(vzVpd, vzVpd-zvtx);


	}



	Double_t vtxx = tofr->vertexX;
	Double_t vtxy = tofr->vertexY;
	hvxvvy->Fill(vtxx,vtxy);
	if(fabs(vzVpd-zvtx)<=vzCut){
		hvxvvyvzdiff->Fill(vtxx,vtxy);
	}

      
      if(tofr->nTofHits>0) T0 = (sumEast+sumWest+(Iwest-Ieast)*tofr->vertexZ/c_light)/(Ieast+Iwest); 
      else T0 = -1e6;
	if(tofr->nTofHits>0&&(fabs(vzVpd-zvtx)<=vzCut)) T0vzdiff = (sumEast+sumWest+(Iwest-Ieast)*tofr->vertexZ/c_light)/(Ieast+Iwest); 
      else T0vzdiff = -1e6;
      

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
      T0vzdiff = -1.0e6;
      vzVpd = -1.0e6;
    }
    if(i<10) cout<<"T0/vzVpd/Iwest/Ieast/sumWest/sumEast = "<<T0<<"/"<<vzVpd<<"/"<<Iwest<<"/"<<Ieast<<"/"<<sumWest<<"/"<<sumEast<<endl;

    pvpd->Fill();
  } 

  forwrite->cd();
  pvpd->Write();
  hResoPVPD11->Write();
  hResoPVPD22->Write();
  hResoPVPD33->Write();
  hVzCorr->Write();
  hVzCorr2->Write();


  hVzCorrvzdiff->Write();
  hVzCorr2vzdiff->Write();


  hvxvvy->Write();
  hvxvvyvzdiff->Write();
  hdeltatwest->Write();
  hdeltateast->Write();
  forwrite->Close();

  
  TFile *f = new TFile("pvpdCheck.root","UPDATE");
  for(int i=0;i<nPVPDChannel;i++) {
    hPVPDMeanvsIt[i]->Write();
    hPVPDFitMeanvsIt[i]->Write();
    hPVPDFitSigmavsIt[i]->Write();
    hPVPDResovsIt[i]->Write();
    
    
  }
  f->Close();
  delete f;
  

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

  const Float_t vzOffset = 0;

  Double_t tot[nPVPDChannel];
  Double_t tdc[nPVPDChannel];
  Double_t cor[nPVPDChannel];
  Double_t mVPDLeTime[nPVPDChannel];
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
    if ( skipBadRun( tofr ))
      continue;

    

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



      tot[j] = j<nPVPDChannel/2? tofr->vpdTotWest[j]:tofr->vpdTotEast[j-nPVPDChannel/2];
      tdc[j] =  j<nPVPDChannel/2? tofr->vpdLeWest[j]:tofr->vpdLeEast[j-nPVPDChannel/2];

         
      
      
    
          if(tot[j]>=pvpdTotEdge[j][0] && tot[j]<pvpdTotEdge[j][maxpvpdNBin]) {
        
       int ibin = -1;
       for(int k=0;k<maxpvpdNBin;k++){
       if(tot[j]>=pvpdTotEdge[j][k]&&tot[j]<pvpdTotEdge[j][k+1]){
       ibin = k;
       break;}
       }
           if(ibin<0) cout<<"error in getting vpd calibration parameters!  tot =        "<<tot[j]<<"    channel "<<j+1<<endl;
       cor[j] = linearInter(pvpdTotEdge[j][ibin],pvpdTotEdge[j][ibin+1],
                             pvpdTotCorr[j][ibin],pvpdTotCorr[j][ibin+1], tot[j]);
        mVPDLeTime[j]=tdc[j]-cor[j];

        if(j<nPVPDChannel/2) {

          sumWest += tdc[j] - cor[j];
          Iwest++;
        } else {

          sumEast += tdc[j] - cor[j];
          Ieast++;
        }
                } 
    } 

    if(i<10) cout<<"Iwest/Ieast = "<<Iwest<<"/"<<Ieast<<endl;
    Double_t vpdtime;
    
    for(int j=0;j<nPVPDChannel;++j){

      if(tot[j]<pvpdTotEdge[j][0] ||tot[j]>=pvpdTotEdge[j][maxpvpdNBin]) continue;
      if(j<nPVPDChannel/2&&Iwest>1) {
        vpdtime = (( mVPDLeTime[j])*Iwest-sumWest)/(Iwest-1);
        
        if(fabs(vpdtime)>tdiffCut) {
          sumWest -=  mVPDLeTime[j];
           mVPDLeTime[j]=0.;
          Iwest--;
        }
      }
      if(j>=nPVPDChannel/2&&Ieast>1) {
        vpdtime = ( (mVPDLeTime[j])*Ieast-sumEast)/(Ieast-1);
        
        if(fabs(vpdtime)>tdiffCut) {
          sumEast -=  mVPDLeTime[j];
           mVPDLeTime[j] =0.;
          Ieast--;
        }
      }
    }
    Int_t hitIndex[nPVPDChannel];
    Int_t nTube = nPVPDChannel/2;
    TMath::Sort(nTube, &mVPDLeTime[0], &hitIndex[0]);
    int nRejectedWest = (int)(pCut*Iwest+0.5);

    for(int ii=0;ii<nRejectedWest;ii++) {
      int index = hitIndex[ii];

      sumWest -= mVPDLeTime[index];
      mVPDLeTime[index] = 0.;
      Iwest--;


    }

    TMath::Sort(nTube, &mVPDLeTime[nPVPDChannel/2], &hitIndex[nPVPDChannel/2]);
    int nRejectedEast = (int)(pCut*Ieast+0.5);

    for(int ii=0;ii<nRejectedEast;ii++) {
      int index = hitIndex[ii+nPVPDChannel/2] + nPVPDChannel/2;
      sumEast -= mVPDLeTime[index];
      mVPDLeTime[index] = 0.;
      Ieast--;


    }

    sumEast += 2*vzOffset*Ieast/c_light;

    if(Ieast>0||Iwest>0) {
      
      
      
      
      float tdiff;
     if(Ieast>0&&Iwest>0) {
      tdiff = sumEast/Ieast - sumWest/Iwest;
      vzVpd = tdiff/2*c_light;
     } else {
       vzVpd = -1e6;
     }
 
      float dcaRmin = 9999;
      Int_t itmp = -1;
      for(int k=0;k<tofr->nTofHits;++k) {
        if(dcaRmin>sqrt(pow(tofr->dcaX[k],2)+pow(tofr->dcaY[k],2))&&TMath::Abs(tofr->dcaZ[k]-vzVpd)<=vzCut&&tofr->nHitsFit[k]>=nHitsFitCut){
          dcaRmin = sqrt(pow(tofr->dcaX[k],2)+pow(tofr->dcaY[k],2));
          itmp = k;
          
        }
                

        } 


          
          hVzCorr->Fill(vzVpd, tofr->vertexZ);
          hVzCorr2->Fill(vzVpd, vzVpd - tofr->vertexZ);
	 if(fabs(vzVpd-tofr->vertexZ)<=vzCut){
	hVzCorrvzdiff->Fill(vzVpd, tofr->vertexZ);
          hVzCorr2vzdiff->Fill(vzVpd, vzVpd - tofr->vertexZ);
	}




	if(tofr->nTofHits>0 ) T0 = (sumEast+sumWest+(Iwest-Ieast)*tofr->vertexZ/c_light)/(Ieast+Iwest);
      else T0 = -1e6;
      

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
  } 

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

  return 1;
}


Int_t checkReso(TChain* chain, Int_t iter)
{

  char rootfile[500];
  sprintf(rootfile, "%s.root", pvpdOutput);


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

  TOTPicoDst *slew = new TOTPicoDst(slewchain);
  ZPicoDst *zhit = new ZPicoDst(zhitchain);
  DelayPicoDst *delay = new DelayPicoDst(delaychain);

  Int_t nevt = (int)(chain->GetEntries());




  if(nevt!=slewchain->GetEntries()) {
    cout << " Inconsistent TOFr and Slew entry ! " << endl;
    return -1;
  }
  if(nevt!=zhitchain->GetEntries()) {
    cout << " Inconsistent TOFr and Zhit entry ! " << endl;
    return -1;
  }

  
  for (Int_t i=0;i<nevt;i++) {
    tofr->GetEntry(i);
    if ( skipBadRun( tofr ))
      continue;
   
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

      if(TMath::Abs(tofr->vertexZ-zhit->vzVpd)>vzCut) continue;

      
      Int_t fitpts = tofr->nHitsFit[j];
      Double_t dy = tofr->yLocal[j];


      if(fitpts<25) continue;
      if(TMath::Abs(dy)>1.5) continue;

      
      
      
      
      
      
      

      Double_t tracklength = tofr->length[j];
      Double_t betagamma = p/piMass;
      Double_t velocity = TMath::Sqrt(1.0/(1.0/betagamma/betagamma+1.0))*c_light;
      Double_t tofPi = tracklength/velocity;

      Int_t tray = tofr->tray[j]-1-nTray*trayStart;
      Int_t module = (tofr->module[j]-1)%4;
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

      if(dedx > dedxcut) continue; 
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

      if(tofr->charge[j]==1) {
        dt = tof_s-tofPi;
        hResoTOT_Pos->Fill(p, dt);
        if(p>0.3&&p<0.6) hResoTOT1D_Pos->Fill(dt);

        dt = tof_z-tofPi;
        hResoZ_Pos->Fill(p, dt);
        if(p>0.3&&p<0.6) hResoZ1D_Pos->Fill(dt); 

        dt = tof_d-tofPi;
        hResoD_Pos->Fill(p,dt);
        if(p>0.3&&p<0.6) hResoD1D_Pos->Fill(dt);
      } else if(tofr->charge[j]==-1) {
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
  } 

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
 
  if(x2==x1) {
    cout<<"error in linear interation, x2==x1"<<endl;
    return 0;
  }
  return ((x-x1)*y2+(x2-x)*y1)/(x2-x1);
}

Bool_t xpisnan(Double_t x){
  return (x!=x);
}

