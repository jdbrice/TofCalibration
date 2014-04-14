

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

#include "initTOFr.h"


Int_t initTOFr()
{
  

  
  piBetaCut = new TF1("piBetaCut","sqrt([0]*[0]+x*x)/x",0,5);
  piBetaCut->SetParameter(0,massCut);


  Char_t buff[100];
  for(Int_t i=0; i<nTray; i++) {
    for(Int_t j=0; j<nBoard; j++) {    
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