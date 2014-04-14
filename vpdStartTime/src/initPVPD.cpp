

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


#include "initPVPD.h"



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